"""Optional high-top grid construction and passive top-sponge diagnostics.

The default solver path does not import or instantiate the observer.  The grid
helper copies the control faces exactly and appends a geometric continuation;
the final cell is truncated so the requested top is reached exactly.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from meteorological_flow.grid import Grid


def continued_top_faces(control_zf, new_top_m: float = 20_000.0) -> np.ndarray:
    """Return exact control faces followed by a geometric continuation.

    All input face values are copied, rather than regenerated or globally
    rescaled.  The ratio of the final two control-cell thicknesses is continued
    until the next complete cell would cross ``new_top_m``; one positive final
    cell terminates exactly at the requested height.
    """
    base = np.asarray(control_zf, dtype=np.float64)
    if base.ndim != 1 or len(base) < 4 or base[0] != 0.0:
        raise ValueError("control_zf must be a one-dimensional face array starting at zero")
    if not np.all(np.isfinite(base)) or not np.all(np.diff(base) > 0):
        raise ValueError("control_zf must be finite and strictly increasing")
    if not np.isfinite(new_top_m) or new_top_m <= base[-1]:
        raise ValueError("new_top_m must be finite and above the control top")
    dz = np.diff(base)
    ratio = dz[-1] / dz[-2]
    if ratio <= 0 or not np.isfinite(ratio):
        raise ValueError("the terminal control-cell ratio must be positive and finite")
    out = list(base)
    next_dz = dz[-1] * ratio
    while out[-1] + next_dz < new_top_m:
        out.append(out[-1] + next_dz)
        next_dz *= ratio
    out.append(float(new_top_m))
    return np.asarray(out, dtype=np.float64)


def configure_high_top(scfg, new_top_m: float = 20_000.0):
    """Configure an opt-in high top while preserving the control lower grid.

    The returned object is a deep copy.  The legacy scalar ``grid.dz`` used by
    LES/filter-width and microphysics code remains equal to the control value,
    avoiding a hidden lower-column treatment change.  The control's number of
    damping faces is retained explicitly.
    """
    import copy

    out = copy.deepcopy(scfg)
    sim = out.sim
    from meteorological_flow.backend import get_backend
    backend = get_backend(getattr(sim.performance, "device", "cpu"))
    control = Grid(
        nx=sim.grid.nx, ny=sim.grid.ny, nz=sim.grid.nz,
        Lx=sim.domain.Lx, Ly=sim.domain.Ly, Lz=sim.domain.Lz,
        z_stretch=sim.grid.z_stretch, backend=backend,
    )
    faces = continued_top_faces(backend.to_cpu(control.zf), new_top_m)
    reference_dz = float(control.dz)
    damping_faces = (sim.boundaries.damping_faces if sim.boundaries.damping_faces is not None
                     else max(2, sim.grid.nz // 10))
    sim.domain.Lz = float(new_top_m)
    sim.grid.nz = len(faces) - 1
    sim.grid.z_faces_m = faces.tolist()
    sim.grid.vertical_reference_dz_m = reference_dz
    sim.boundaries.damping_faces = int(damping_faces)
    return out


def _horizontal_profile(xp, field):
    mean = xp.mean(field, axis=(0, 1))
    rms = xp.sqrt(xp.mean(field * field, axis=(0, 1)))
    return mean, rms


class TopBoundaryObserver:
    """Lossless-in-statistics passive observer for every top damping call.

    It stores exact operator energy removal plus vertical profiles sufficient
    for propagation and lag analysis.  It never writes to solver-owned arrays.
    """
    contexts = ("pre_predictor", "post_predictor", "post_projection", "post_transport")

    def __init__(self, sim, path, metadata=None, flush_every=128,
                 full_field_start_s=None):
        import h5py

        self.path = Path(path)
        self.file = h5py.File(self.path, "x")
        self.file.attrs["schema"] = "storm-top-boundary-observer-v1"
        self.file.attrs["status"] = "running"
        self.file.attrs["passive"] = True
        self.file.attrs["metadata"] = json.dumps(metadata or {}, default=str)
        self.file.attrs["energy_definition"] = (
            "0.5*rho0_wface*dx*dy*dual_face_dz*(w_before^2-w_after^2)"
        )
        self.file.attrs["lag_sign"] = "positive lag means top activity precedes lower-level response"
        self.to = sim.grid.backend.to_cpu
        self.xp = sim.grid.xp
        self.grid = sim.grid
        self.rho0_wface = self.xp.asarray(sim.rho0_wface)
        dz = np.asarray(self.to(sim.grid.dz_c), dtype=float)
        dual = np.empty(sim.grid.nz + 1, dtype=float)
        dual[0] = 0.5 * dz[0]
        dual[-1] = 0.5 * dz[-1]
        dual[1:-1] = 0.5 * (dz[:-1] + dz[1:])
        self.dual_dz = self.xp.asarray(dual)
        zf = np.asarray(self.to(sim.grid.zf), dtype=np.float64)
        self.file.create_dataset("zf_m", data=zf)
        self.file.attrs["zf_sha256"] = hashlib.sha256(zf.tobytes()).hexdigest()
        self.file.attrs["vertical_reference_dz_m"] = float(sim.grid.dz)
        self.flush_every = int(flush_every)
        self.full_field_start_s = (None if full_field_start_s is None
                                   else float(full_field_start_s))
        self.file.attrs["full_field_start_s"] = (-1.0 if self.full_field_start_s is None
                                                  else self.full_field_start_s)
        self.buffer = []
        self.datasets = None
        self.full_buffer = []
        self.full_datasets = None
        self.calls = 0

    def _profiles(self, state):
        xp = self.xp
        w = 0.5 * (state.w[:, :, :-1] + state.w[:, :, 1:])
        w_mean, w_rms = _horizontal_profile(xp, w)
        p = getattr(state, "p_dyn", None)
        if p is None:
            p = xp.zeros_like(w)
        p_mean, p_rms = _horizontal_profile(xp, p)
        pw_mean, pw_rms = _horizontal_profile(xp, p * w)
        u = 0.5 * (state.u[:-1] + state.u[1:])
        v = 0.5 * (state.v[:, :-1] + state.v[:, 1:])
        if self.grid.periodic:
            dudx = (xp.roll(u, -1, axis=0) - xp.roll(u, 1, axis=0)) / (2 * self.grid.dx)
            dvdy = (xp.roll(v, -1, axis=1) - xp.roll(v, 1, axis=1)) / (2 * self.grid.dy)
        else:
            dudx = xp.gradient(u, self.grid.dx, axis=0)
            dvdy = xp.gradient(v, self.grid.dy, axis=1)
        conv_mean, conv_rms = _horizontal_profile(xp, -(dudx + dvdy))
        return tuple(np.asarray(self.to(a), dtype=float) for a in (
            w_mean, w_rms, p_mean, p_rms, pw_mean, pw_rms, conv_mean, conv_rms
        ))

    def record(self, *, state, grid, indices, multipliers, before,
               context, step, time_s, dt):
        if grid is not self.grid:
            raise RuntimeError("top-boundary observer attached to a different grid")
        xp = self.xp
        idx_desc = np.asarray(indices, dtype=int)
        order = np.argsort(idx_desc)
        idx = idx_desc[order]
        before_desc = before
        after_desc = state.w[:, :, list(indices)]
        delta_desc = after_desc - before_desc
        before_a = before_desc[:, :, order]
        after_a = after_desc[:, :, order]
        delta_a = delta_desc[:, :, order]
        multiplier_a = xp.asarray(multipliers)[order]

        mass = (self.rho0_wface[idx] * float(grid.dx) * float(grid.dy)
                * self.dual_dz[idx])
        energy = 0.5 * xp.sum(
            (before_a * before_a - after_a * after_a) * mass[None, None, :], axis=(0, 1)
        )
        delta_mean, delta_rms = _horizontal_profile(xp, delta_a)
        delta_absmax = xp.max(xp.abs(delta_a), axis=(0, 1))

        # curl(0,0,delta_w): xi=d(delta_w)/dy, eta=-d(delta_w)/dx, zeta=0.
        lower = xp.zeros(delta_a.shape[:2] + (1,), dtype=delta_a.dtype)
        face_delta = xp.concatenate((lower, delta_a), axis=2)
        centered_delta = 0.5 * (face_delta[:, :, :-1] + face_delta[:, :, 1:])
        if grid.periodic:
            dxi = (xp.roll(centered_delta, -1, axis=1)
                   - xp.roll(centered_delta, 1, axis=1)) / (2 * grid.dy)
            deta = -(xp.roll(centered_delta, -1, axis=0)
                     - xp.roll(centered_delta, 1, axis=0)) / (2 * grid.dx)
        else:
            dxi = xp.gradient(centered_delta, grid.dy, axis=1)
            deta = -xp.gradient(centered_delta, grid.dx, axis=0)
        xi_rms = xp.sqrt(xp.mean(dxi * dxi, axis=(0, 1)))
        eta_rms = xp.sqrt(xp.mean(deta * deta, axis=(0, 1)))
        xi_absmax = xp.max(xp.abs(dxi), axis=(0, 1))
        eta_absmax = xp.max(xp.abs(deta), axis=(0, 1))

        p = getattr(state, "p_dyn", None)
        if p is None:
            pface = xp.zeros_like(after_a)
        else:
            pieces = []
            for k in idx:
                pieces.append(p[:, :, -1] if k == grid.nz else
                              (p[:, :, 0] if k == 0 else 0.5 * (p[:, :, k - 1] + p[:, :, k])))
            pface = xp.stack(pieces, axis=2)
        area = float(grid.dx) * float(grid.dy)
        pressure_flux_before = xp.sum(pface * before_a, axis=(0, 1)) * area
        pressure_flux_after = xp.sum(pface * after_a, axis=(0, 1)) * area
        profiles = self._profiles(state)
        self.buffer.append(dict(
            time_s=float(time_s), dt=float(dt), step=int(step),
            ordinal=self.contexts.index(context), context=str(context),
            face_indices=idx, multipliers=np.asarray(self.to(multiplier_a), dtype=float),
            delta_w_mean=np.asarray(self.to(delta_mean), dtype=float),
            delta_w_rms=np.asarray(self.to(delta_rms), dtype=float),
            delta_w_absmax=np.asarray(self.to(delta_absmax), dtype=float),
            energy_removed_J=np.asarray(self.to(energy), dtype=float),
            delta_xi_rms=np.asarray(self.to(xi_rms), dtype=float),
            delta_xi_absmax=np.asarray(self.to(xi_absmax), dtype=float),
            delta_eta_rms=np.asarray(self.to(eta_rms), dtype=float),
            delta_eta_absmax=np.asarray(self.to(eta_absmax), dtype=float),
            pressure_flux_before_W=np.asarray(self.to(pressure_flux_before), dtype=float),
            pressure_flux_after_W=np.asarray(self.to(pressure_flux_after), dtype=float),
            w_mean=profiles[0], w_rms=profiles[1], p_dyn_mean=profiles[2],
            p_dyn_rms=profiles[3], p_dyn_w_mean=profiles[4], p_dyn_w_rms=profiles[5],
            convergence_mean=profiles[6], convergence_rms=profiles[7],
        ))
        if self.full_field_start_s is not None and float(time_s) >= self.full_field_start_s:
            # Lossless spatial w_before plus delta_w; w_after is reconstructed
            # exactly as their sum.  Restricting this to the analysis window
            # avoids tens of GB of redundant spin-up storage.
            self.full_buffer.append(dict(
                time_s=float(time_s), step=int(step), ordinal=self.contexts.index(context),
                face_indices=idx,
                w_before=np.asarray(self.to(before_a)),
                delta_w=np.asarray(self.to(delta_a)),
            ))
        self.calls += 1
        if len(self.buffer) >= self.flush_every:
            self.flush()

    def _create_datasets(self, row, group=None, chunk_rows=None):
        import h5py

        target = self.file if group is None else self.file.require_group(group)
        datasets = {}
        chunk_rows = min(self.flush_every, 128) if chunk_rows is None else int(chunk_rows)
        for name, value in row.items():
            if name == "context":
                dtype = h5py.string_dtype("ascii", length=24)
                shape = (); data_dtype = dtype
            else:
                array = np.asarray(value)
                shape = array.shape
                data_dtype = array.dtype
            datasets[name] = target.create_dataset(
                name, shape=(0,) + shape, maxshape=(None,) + shape,
                dtype=data_dtype, chunks=(chunk_rows,) + shape,
                compression="gzip" if shape else None, compression_opts=1 if shape else None,
                shuffle=bool(shape),
            )
        if group is None:
            self.file.attrs["profile_columns"] = (
                "w,p_dyn,p_dyn*w,convergence each stored as horizontal mean and RMS at cell centers"
            )
            self.file.attrs["curl_cell_indices"] = json.dumps(
                list(range(int(row["face_indices"][0]) - 1, self.grid.nz))
            )
        return datasets

    def flush(self):
        if not self.buffer:
            pass
        else:
            if self.datasets is None:
                self.datasets = self._create_datasets(self.buffer[0])
            start = self.datasets["step"].shape[0]
            stop = start + len(self.buffer)
            for name, ds in self.datasets.items():
                ds.resize((stop,) + ds.shape[1:])
                ds[start:stop] = [row[name] for row in self.buffer]
            self.buffer.clear()
        if self.full_buffer:
            if self.full_datasets is None:
                self.full_datasets = self._create_datasets(
                    self.full_buffer[0], group="full_fields", chunk_rows=1
                )
                self.file["full_fields"].attrs["w_after_reconstruction"] = "w_before + delta_w"
            start = self.full_datasets["step"].shape[0]
            stop = start + len(self.full_buffer)
            for name, ds in self.full_datasets.items():
                ds.resize((stop,) + ds.shape[1:])
                ds[start:stop] = [row[name] for row in self.full_buffer]
            self.full_buffer.clear()
        self.file.flush()

    def close(self, completed=True):
        self.flush()
        self.file.attrs["calls"] = self.calls
        self.file.attrs["status"] = "complete" if completed else "interrupted"
        self.file.flush()
        self.file.close()


__all__ = ["TopBoundaryObserver", "configure_high_top", "continued_top_faces"]
