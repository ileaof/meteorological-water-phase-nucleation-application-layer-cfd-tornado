"""Passive, opt-in provenance tracers for resolved velocity and vorticity.

The production tracer decomposes the native staggered velocity by source.  At
advection, limiter decisions are frozen from the total velocity and the same
linear MUSCL momentum-flux map is applied to every label.  Native split-stage
velocity increments are injected into the corresponding label.  Taking the
curl of each labelled velocity then gives an additive vorticity provenance.
Roundoff left by the frozen flux split is retained explicitly under
``advection_remainder``.

Nothing in this module is consulted by the prognostic solver.  The observer only
reads/copies solver arrays through the existing ``diagnostic_observer`` hooks.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


SOURCE_NAMES = (
    "initial",
    "buoyancy",
    "les",
    "surface_drag",
    "coriolis",
    "projection",
    "boundary",
    "other",
    "advection_remainder",
)

_STAGE_SOURCE = {
    "initial_bcs": "boundary",
    "predictor_bcs": "boundary",
    "projection_bcs": "boundary",
    "transport_microphysics_bcs": "boundary",
    "buoyancy": "buoyancy",
    "les": "les",
    "surface_drag": "surface_drag",
    "coriolis": "coriolis",
    "projection": "projection",
    "external_forcing": "other",
    "guard": "other",
}


def _centered(velocity):
    u, v, w = velocity
    return (
        0.5 * (u[:-1, :, :] + u[1:, :, :]),
        0.5 * (v[:, :-1, :] + v[:, 1:, :]),
        0.5 * (w[:, :, :-1] + w[:, :, 1:]),
    )


def curl_native(velocity, grid):
    """Curl of a C-grid velocity tuple at cell centres, using solver diagnostics."""
    u, v, w = _centered(velocity)
    return grid.xp.stack(
        (
            grid._central_y(w) - grid._central_z(v),
            grid._central_z(u) - grid._central_x(w),
            grid._central_x(v) - grid._central_y(u),
        )
    )


def velocity_kinematics(velocity, grid):
    """Return centred velocity, its gradient tensor, divergence and grad_h(w)."""
    centered = _centered(velocity)
    gradient = tuple(
        (grid._central_x(q), grid._central_y(q), grid._central_z(q))
        for q in centered
    )
    divergence = gradient[0][0] + gradient[1][1] + gradient[2][2]
    return centered, gradient, divergence, (gradient[2][0], gradient[2][1])


def homogeneous_vorticity_tendency(omega, centered_velocity, velocity_gradient,
                                   divergence, grid):
    """Resolved homogeneous vorticity equation for one source-labelled vector.

    Computes ``-u.grad(omega) + omega.grad(u) - omega div(u)``.  The operator is
    linear in ``omega`` once the resolved velocity is fixed, which is the key
    property required for a source-provenance decomposition.
    """
    xp = grid.xp
    u, v, w = centered_velocity
    out = xp.empty_like(omega)
    for component in range(3):
        ox, oy, oz = (
            grid._central_x(omega[component]),
            grid._central_y(omega[component]),
            grid._central_z(omega[component]),
        )
        transport = -(u * ox + v * oy + w * oz)
        stretching = sum(omega[j] * velocity_gradient[component][j] for j in range(3))
        out[component] = transport + stretching - omega[component] * divergence
    return out


def _frozen_slope(label, total, axis, grid, periodic):
    """Linear limited slope using minmod branch decisions from the total field."""
    xp = grid.xp
    if periodic:
        total_fwd = xp.roll(total, -1, axis=axis) - total
        total_bwd = total - xp.roll(total, 1, axis=axis)
        label_fwd = xp.roll(label, -1, axis=axis) - label
        label_bwd = label - xp.roll(label, 1, axis=axis)
        return xp.where(
            total_bwd * total_fwd <= 0.0,
            0.0,
            xp.where(
                xp.abs(total_bwd) < xp.abs(total_fwd),
                label_bwd,
                label_fwd,
            ),
        )

    slope = xp.zeros_like(label)
    interior = [slice(None)] * label.ndim
    minus = [slice(None)] * label.ndim
    plus = [slice(None)] * label.ndim
    interior[axis] = slice(1, -1)
    minus[axis] = slice(0, -2)
    plus[axis] = slice(2, None)
    total_here = total[tuple(interior)]
    total_fwd = total[tuple(plus)] - total_here
    total_bwd = total_here - total[tuple(minus)]
    label_here = label[tuple(interior)]
    label_fwd = label[tuple(plus)] - label_here
    label_bwd = label_here - label[tuple(minus)]
    slope[tuple(interior)] = xp.where(
        total_bwd * total_fwd <= 0.0,
        0.0,
        xp.where(
            xp.abs(total_bwd) < xp.abs(total_fwd),
            label_bwd,
            label_fwd,
        ),
    )
    return slope


def frozen_muscl_advection_tendency(label, total, velocity, grid, order=2):
    """Return stable -div(u label) with reconstruction frozen from total.

    Face velocities are the solver's projected staggered velocities. Upwind
    choices and MUSCL/minmod branches are diagnosed from the total field and
    reused for every source. The operator is consequently linear in label.
    """
    xp = grid.xp
    u, v, w = velocity
    periodic = bool(getattr(grid, "periodic", False))
    sx = xp.empty(grid.u_shape, dtype=label.dtype)
    sy = xp.empty(grid.v_shape, dtype=label.dtype)
    sz = xp.empty(grid.w_shape, dtype=label.dtype)

    if order >= 2:
        dx = _frozen_slope(label, total, 0, grid, periodic)
        dy = _frozen_slope(label, total, 1, grid, periodic)
        dz = _frozen_slope(label, total, 2, grid, False)
        sx[1:-1] = xp.where(
            u[1:-1] > 0.0,
            label[:-1] + 0.5 * dx[:-1],
            label[1:] - 0.5 * dx[1:],
        )
        sy[:, 1:-1] = xp.where(
            v[:, 1:-1] > 0.0,
            label[:, :-1] + 0.5 * dy[:, :-1],
            label[:, 1:] - 0.5 * dy[:, 1:],
        )
        sz[:, :, 1:-1] = xp.where(
            w[:, :, 1:-1] > 0.0,
            label[:, :, :-1] + 0.5 * dz[:, :, :-1],
            label[:, :, 1:] - 0.5 * dz[:, :, 1:],
        )
    else:
        sx[1:-1] = xp.where(u[1:-1] > 0.0, label[:-1], label[1:])
        sy[:, 1:-1] = xp.where(v[:, 1:-1] > 0.0, label[:, :-1], label[:, 1:])
        sz[:, :, 1:-1] = xp.where(
            w[:, :, 1:-1] > 0.0, label[:, :, :-1], label[:, :, 1:]
        )

    if periodic:
        if order >= 2:
            sx0 = xp.where(
                u[0] > 0.0,
                label[-1] + 0.5 * dx[-1],
                label[0] - 0.5 * dx[0],
            )
            sy0 = xp.where(
                v[:, 0] > 0.0,
                label[:, -1] + 0.5 * dy[:, -1],
                label[:, 0] - 0.5 * dy[:, 0],
            )
        else:
            sx0 = xp.where(u[0] > 0.0, label[-1], label[0])
            sy0 = xp.where(v[:, 0] > 0.0, label[:, -1], label[:, 0])
        sx[0] = sx0
        sx[-1] = sx0
        sy[:, 0] = sy0
        sy[:, -1] = sy0
    else:
        sx[0] = label[0]
        sx[-1] = label[-1]
        sy[:, 0] = label[:, 0]
        sy[:, -1] = label[:, -1]
    sz[:, :, 0] = label[:, :, 0]
    sz[:, :, -1] = label[:, :, -1]

    fx = u * sx
    fy = v * sy
    fz = w * sz
    spacing_z = (
        grid.dz
        if not getattr(grid, "stretched", False)
        else grid.dz_c[None, None, :]
    )
    return -(
        (fx[1:] - fx[:-1]) / grid.dx
        + (fy[:, 1:] - fy[:, :-1]) / grid.dy
        + (fz[:, :, 1:] - fz[:, :, :-1]) / spacing_z
    )


def discrete_homogeneous_vorticity_tendency(
    omega, total_omega, velocity, velocity_gradient, grid, order=2
):
    """Stable frozen-MUSCL form of -div(u omega) + (omega.grad)u."""
    xp = grid.xp
    out = xp.empty_like(omega)
    for component in range(3):
        transport = frozen_muscl_advection_tendency(
            omega[component],
            total_omega[component],
            velocity,
            grid,
            order=order,
        )
        stretching = sum(
            omega[j] * velocity_gradient[component][j] for j in range(3)
        )
        out[component] = transport + stretching
    return out


def _frozen_face(qm, qp, velocity, slope_m, slope_p, xp, order):
    if order >= 2:
        left = qm + 0.5 * slope_m
        right = qp - 0.5 * slope_p
        return xp.where(velocity > 0.0, left, right)
    return xp.where(velocity > 0.0, qm, qp)


def _frozen_u_tendency(label, total, grid, order, periodic):
    xp = grid.xp
    lu, _, _ = label
    u, v, w = total
    spacing_z = (
        grid.dz
        if not getattr(grid, "stretched", False)
        else grid.dz_c[None, None, :]
    )
    tendency = xp.zeros_like(lu)

    uc = 0.5 * (u[:-1] + u[1:])
    slope = _frozen_slope(lu, u, 0, grid, periodic)
    flux_x = uc * _frozen_face(
        lu[:-1], lu[1:], uc, slope[:-1], slope[1:], xp, order
    )
    tendency[1:-1] += -(flux_x[1:] - flux_x[:-1]) / grid.dx
    if periodic:
        wrap = -(flux_x[0] - flux_x[-1]) / grid.dx
        tendency[0] += wrap
        tendency[-1] += wrap

    v_xm = xp.roll(v, 1, axis=0) if periodic else xp.zeros_like(v)
    if not periodic:
        v_xm[1:] = v[:-1]
    v_on_u = 0.5 * (v_xm + v)
    v_full = xp.zeros((grid.nx + 1, grid.ny + 1, grid.nz))
    v_full[:-1] = v_on_u
    v_full[-1] = v_on_u[0] if periodic else 0.0
    slope_y = _frozen_slope(lu, u, 1, grid, periodic)
    if periodic:
        q_minus = xp.roll(lu, 1, axis=1)
        s_minus = xp.roll(slope_y, 1, axis=1)
        vc = v_full[:, :-1]
        flux_y = vc * _frozen_face(
            q_minus, lu, vc, s_minus, slope_y, xp, order
        )
        tendency += -(xp.roll(flux_y, -1, axis=1) - flux_y) / grid.dy
    else:
        q_minus = xp.zeros_like(lu)
        q_minus[:, 1:] = lu[:, :-1]
        s_minus = xp.zeros_like(slope_y)
        s_minus[:, 1:] = slope_y[:, :-1]
        vc = v_full[:, :-1]
        flux_full = xp.zeros((grid.nx + 1, grid.ny + 1, grid.nz))
        flux_full[:, :-1] = vc * _frozen_face(
            q_minus, lu, vc, s_minus, slope_y, xp, order
        )
        tendency[:, 1:-1] += -(
            flux_full[:, 2:-1] - flux_full[:, 1:-2]
        ) / grid.dy

    w_xm = xp.roll(w, 1, axis=0) if periodic else xp.zeros_like(w)
    if not periodic:
        w_xm[1:] = w[:-1]
    w_on_u = 0.5 * (w_xm + w)
    w_full = xp.zeros((grid.nx + 1, grid.ny, grid.nz + 1))
    w_full[:-1] = w_on_u
    w_full[-1] = w_on_u[0] if periodic else 0.0
    slope_z = _frozen_slope(lu, u, 2, grid, False)
    q_minus = xp.zeros_like(lu)
    q_minus[:, :, 1:] = lu[:, :, :-1]
    s_minus = xp.zeros_like(slope_z)
    s_minus[:, :, 1:] = slope_z[:, :, :-1]
    wc = w_full[:, :, 1:-1]
    flux_z = xp.zeros((grid.nx + 1, grid.ny, grid.nz + 1))
    flux_z[:, :, 1:-1] = wc * _frozen_face(
        q_minus[:, :, 1:],
        lu[:, :, 1:],
        wc,
        s_minus[:, :, 1:],
        slope_z[:, :, 1:],
        xp,
        order,
    )
    tendency += -(flux_z[:, :, 1:] - flux_z[:, :, :-1]) / spacing_z
    return tendency


def _frozen_v_tendency(label, total, grid, order, periodic):
    xp = grid.xp
    _, lv, _ = label
    u, v, w = total
    spacing_z = (
        grid.dz
        if not getattr(grid, "stretched", False)
        else grid.dz_c[None, None, :]
    )
    tendency = xp.zeros_like(lv)

    vc = 0.5 * (v[:, :-1] + v[:, 1:])
    slope = _frozen_slope(lv, v, 1, grid, periodic)
    flux_y = vc * _frozen_face(
        lv[:, :-1], lv[:, 1:], vc, slope[:, :-1], slope[:, 1:], xp, order
    )
    tendency[:, 1:-1] += -(flux_y[:, 1:] - flux_y[:, :-1]) / grid.dy
    if periodic:
        wrap = -(flux_y[:, 0] - flux_y[:, -1]) / grid.dy
        tendency[:, 0] += wrap
        tendency[:, -1] += wrap

    u_ym = xp.roll(u, 1, axis=1) if periodic else xp.zeros_like(u)
    if not periodic:
        u_ym[:, 1:] = u[:, :-1]
    u_on_v = 0.5 * (u_ym + u)
    u_full = xp.zeros((grid.nx + 1, grid.ny + 1, grid.nz))
    u_full[:, :-1] = u_on_v
    u_full[:, -1] = u_on_v[:, 0] if periodic else 0.0
    slope_x = _frozen_slope(lv, v, 0, grid, periodic)
    if periodic:
        q_minus = xp.roll(lv, 1, axis=0)
        s_minus = xp.roll(slope_x, 1, axis=0)
        uc = u_full[:-1]
        flux_x = uc * _frozen_face(
            q_minus, lv, uc, s_minus, slope_x, xp, order
        )
        tendency += -(xp.roll(flux_x, -1, axis=0) - flux_x) / grid.dx
    else:
        q_minus = xp.zeros_like(lv)
        q_minus[1:] = lv[:-1]
        s_minus = xp.zeros_like(slope_x)
        s_minus[1:] = slope_x[:-1]
        uc = u_full[:-1]
        flux_full = xp.zeros((grid.nx + 1, grid.ny + 1, grid.nz))
        flux_full[:-1] = uc * _frozen_face(
            q_minus, lv, uc, s_minus, slope_x, xp, order
        )
        tendency[1:-1] += -(
            flux_full[2:-1] - flux_full[1:-2]
        ) / grid.dx

    w_ym = xp.roll(w, 1, axis=1) if periodic else xp.zeros_like(w)
    if not periodic:
        w_ym[:, 1:] = w[:, :-1]
    w_on_v = 0.5 * (w_ym + w)
    w_full = xp.zeros((grid.nx, grid.ny + 1, grid.nz + 1))
    w_full[:, :-1] = w_on_v
    w_full[:, -1] = w_on_v[:, 0] if periodic else 0.0
    slope_z = _frozen_slope(lv, v, 2, grid, False)
    q_minus = xp.zeros_like(lv)
    q_minus[:, :, 1:] = lv[:, :, :-1]
    s_minus = xp.zeros_like(slope_z)
    s_minus[:, :, 1:] = slope_z[:, :, :-1]
    wc = w_full[:, :, 1:-1]
    flux_z = xp.zeros((grid.nx, grid.ny + 1, grid.nz + 1))
    flux_z[:, :, 1:-1] = wc * _frozen_face(
        q_minus[:, :, 1:],
        lv[:, :, 1:],
        wc,
        s_minus[:, :, 1:],
        slope_z[:, :, 1:],
        xp,
        order,
    )
    tendency += -(flux_z[:, :, 1:] - flux_z[:, :, :-1]) / spacing_z
    return tendency


def _frozen_w_tendency(label, total, grid, order, periodic):
    xp = grid.xp
    _, _, lw = label
    u, v, w = total
    centre_spacing = (
        grid.dz
        if not getattr(grid, "stretched", False)
        else grid.dzc_f[None, None, 1:-1]
    )
    tendency = xp.zeros_like(lw)

    wc = 0.5 * (w[:, :, :-1] + w[:, :, 1:])
    slope = _frozen_slope(lw, w, 2, grid, False)
    flux_z = wc * _frozen_face(
        lw[:, :, :-1],
        lw[:, :, 1:],
        wc,
        slope[:, :, :-1],
        slope[:, :, 1:],
        xp,
        order,
    )
    tendency[:, :, 1:-1] += -(
        flux_z[:, :, 1:] - flux_z[:, :, :-1]
    ) / centre_spacing

    u_on_w = xp.zeros((grid.nx + 1, grid.ny, grid.nz + 1))
    u_on_w[:, :, 1:-1] = 0.5 * (u[:, :, :-1] + u[:, :, 1:])
    u_on_w[:, :, 0] = u[:, :, 0]
    u_on_w[:, :, -1] = u[:, :, -1]
    slope_x = _frozen_slope(lw, w, 0, grid, periodic)
    if periodic:
        q_minus = xp.roll(lw, 1, axis=0)
        s_minus = xp.roll(slope_x, 1, axis=0)
        uc = u_on_w[:-1]
        flux_x = uc * _frozen_face(
            q_minus, lw, uc, s_minus, slope_x, xp, order
        )
        tendency += -(xp.roll(flux_x, -1, axis=0) - flux_x) / grid.dx
    else:
        q_minus = xp.zeros_like(lw)
        q_minus[1:] = lw[:-1]
        s_minus = xp.zeros_like(slope_x)
        s_minus[1:] = slope_x[:-1]
        uc = u_on_w[:-1]
        flux_full = xp.zeros((grid.nx + 1, grid.ny, grid.nz + 1))
        flux_full[:-1] = uc * _frozen_face(
            q_minus, lw, uc, s_minus, slope_x, xp, order
        )
        tendency[1:-1] += -(
            flux_full[2:-1] - flux_full[1:-2]
        ) / grid.dx

    v_on_w = xp.zeros((grid.nx, grid.ny + 1, grid.nz + 1))
    v_on_w[:, :, 1:-1] = 0.5 * (v[:, :, :-1] + v[:, :, 1:])
    v_on_w[:, :, 0] = v[:, :, 0]
    v_on_w[:, :, -1] = v[:, :, -1]
    slope_y = _frozen_slope(lw, w, 1, grid, periodic)
    if periodic:
        q_minus = xp.roll(lw, 1, axis=1)
        s_minus = xp.roll(slope_y, 1, axis=1)
        vc = v_on_w[:, :-1]
        flux_y = vc * _frozen_face(
            q_minus, lw, vc, s_minus, slope_y, xp, order
        )
        tendency += -(xp.roll(flux_y, -1, axis=1) - flux_y) / grid.dy
    else:
        q_minus = xp.zeros_like(lw)
        q_minus[:, 1:] = lw[:, :-1]
        s_minus = xp.zeros_like(slope_y)
        s_minus[:, 1:] = slope_y[:, :-1]
        vc = v_on_w[:, :-1]
        flux_full = xp.zeros((grid.nx, grid.ny + 1, grid.nz + 1))
        flux_full[:, :-1] = vc * _frozen_face(
            q_minus, lw, vc, s_minus, slope_y, xp, order
        )
        tendency[:, 1:-1] += -(
            flux_full[:, 2:-1] - flux_full[:, 1:-2]
        ) / grid.dy

    tendency[:, :, 0] = 0.0
    return tendency


def frozen_momentum_advection_tendency(
    label_velocity, total_velocity, grid, order=2, periodic=None
):
    """Linear source transport using the native staggered momentum flux geometry."""
    if periodic is None:
        periodic = bool(getattr(grid, "periodic", False))
    return (
        _frozen_u_tendency(label_velocity, total_velocity, grid, order, periodic),
        _frozen_v_tendency(label_velocity, total_velocity, grid, order, periodic),
        _frozen_w_tendency(label_velocity, total_velocity, grid, order, periodic),
    )


class DiagnosticObserverChain:
    """Fan out the solver's read-only diagnostic hooks to several observers."""

    def __init__(self, *observers):
        self.observers = tuple(observers)

    def mark(self, sim, name, dt):
        for observer in self.observers:
            observer.mark(sim, name, dt)

    def close(self, sim, completed=True):
        for observer in self.observers:
            close = getattr(observer, "close", None)
            if close is not None:
                close(sim, completed=completed)


class VorticityProvenanceTracer:
    """Evolve passive source-labelled vorticity alongside ``StormSimulation``.

    The fields are carried over the full column so vertical import into the
    0--2 km analysis layer retains its source label.  Only the requested lower
    layer plus output halos is written to HDF5.
    """

    def __init__(self, sim, path=None, metadata=None, interval=30.0, top=2000.0,
                 output_halo=3, storage_dtype="f8", record_stage_metrics=True,
                 source_checkpoint=None, write_full_vorticity=False):
        self.grid = sim.grid
        self.xp = sim.grid.xp
        self.to_cpu = sim.grid.backend.to_cpu
        self.interval = float(interval)
        self.next_output = float(sim.t) + self.interval
        self.record_stage_metrics = bool(record_stage_metrics)
        self.momentum_order = int(getattr(sim.dyn, "momentum_order", 2))
        self.source_names = SOURCE_NAMES
        self.source_index = {name: i for i, name in enumerate(self.source_names)}
        self.storage_dtype = np.dtype(storage_dtype)
        self.write_full_vorticity = bool(write_full_vorticity)

        z = np.asarray(self.to_cpu(sim.grid.zc), dtype=float)
        self.analysis_nk = max(1, min(sim.grid.nz, int(np.searchsorted(z, top, side="right"))))
        self.output_nk = min(sim.grid.nz, self.analysis_nk + int(output_halo))

        self.previous_velocity = self._native(sim)
        self.total_omega = curl_native(self.previous_velocity, self.grid)
        self.velocity_sources = tuple(
            self.xp.zeros(
                (len(self.source_names),) + tuple(component.shape),
                dtype=component.dtype,
            )
            for component in self.previous_velocity
        )
        initial = self.source_index["initial"]
        for labelled, component in zip(self.velocity_sources, self.previous_velocity):
            labelled[initial] = component
        self.history = []
        self.last_saved_time = float(sim.t)

        self.path = Path(path) if path is not None else None
        self.file = None
        self.snapshots = None
        self.count = 0
        if self.path is not None:
            import h5py
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.file = h5py.File(self.path, "x")
            self.file.attrs["schema"] = "storm-vorticity-provenance-v4"
            self.file.attrs["status"] = "running"
            self.file.attrs["passive"] = True
            self.file.attrs["source_names"] = json.dumps(self.source_names)
            self.file.attrs["analysis_top_m"] = float(top)
            self.file.attrs["analysis_nz"] = self.analysis_nk
            self.file.attrs["output_nz_with_halo"] = self.output_nk
            self.file.attrs["internal_full_column"] = True
            self.file.attrs["storage_dtype"] = self.storage_dtype.str
            self.file.attrs["write_full_vorticity"] = self.write_full_vorticity
            self.file.attrs["formulation"] = (
                "source-labelled staggered velocities advected with the native "
                "MUSCL momentum-flux geometry and limiter choices frozen from "
                "the total velocity; curls are diagnosed after every stage"
            )
            self.file.attrs["metadata"] = json.dumps(metadata or {}, default=str)
            grid_group = self.file.create_group("grid")
            for name in ("xc", "yc"):
                grid_group.create_dataset(name, data=self.to_cpu(getattr(sim.grid, name)))
            grid_group.create_dataset("zc", data=z[:self.output_nk])
            self.snapshots = self.file.create_group("snapshots")
            if source_checkpoint is not None:
                self.restore_full_checkpoint(source_checkpoint, sim)
            self.save(sim, initial=True)

    def _native(self, sim):
        return tuple(getattr(sim.state, name).copy() for name in ("u", "v", "w"))

    @property
    def omega_sources(self):
        return self.xp.stack(
            [
                curl_native(
                    tuple(component[index] for component in self.velocity_sources),
                    self.grid,
                )
                for index in range(len(self.source_names))
            ]
        )

    def _sum_sources(self, omega_sources=None):
        if omega_sources is None:
            omega_sources = self.omega_sources
        return self.xp.sum(omega_sources, axis=0)

    def _closure_metrics(self, velocity=None):
        omega_sources = self.omega_sources
        residual = self.total_omega - self._sum_sources(omega_sources)
        physical = residual[:, :, :, :self.analysis_nk]
        total = self.total_omega[:, :, :, :self.analysis_nk]
        rms = self.xp.sqrt(self.xp.mean(physical * physical))
        total_rms = self.xp.sqrt(self.xp.mean(total * total))
        max_abs = self.xp.max(self.xp.abs(physical))
        rel = rms / self.xp.maximum(total_rms, 1e-30)

        velocity = self.previous_velocity if velocity is None else velocity
        _, _, _, (wx, wy) = velocity_kinematics(velocity, self.grid)
        tilting_total = self.total_omega[0] * wx + self.total_omega[1] * wy
        tilting_sum = self.xp.sum(
            omega_sources[:, 0] * wx[None, ...]
            + omega_sources[:, 1] * wy[None, ...], axis=0
        )
        tr = (tilting_total - tilting_sum)[:, :, :self.analysis_nk]
        trms = self.xp.sqrt(self.xp.mean(tr * tr))
        tscale = self.xp.sqrt(self.xp.mean(tilting_total[:, :, :self.analysis_nk] ** 2))
        source_h = omega_sources[:, :2, :, :, :self.analysis_nk]
        source_rms = self.xp.sqrt(self.xp.mean(source_h * source_h, axis=(1, 2, 3, 4)))
        tilting_sources = (
            omega_sources[:, 0, :, :, :self.analysis_nk]
            * wx[None, :, :, :self.analysis_nk]
            + omega_sources[:, 1, :, :, :self.analysis_nk]
            * wy[None, :, :, :self.analysis_nk]
        )
        tilting_source_rms = self.xp.sqrt(
            self.xp.mean(tilting_sources * tilting_sources, axis=(1, 2, 3))
        )
        horizontal_scale = self.xp.sqrt(self.xp.mean(total[:2] * total[:2]))
        remainder = self.source_index["advection_remainder"]
        return {
            "omega_max_abs": float(max_abs),
            "omega_rms": float(rms),
            "omega_relative_rms": float(rel),
            "tilting_max_abs": float(self.xp.max(self.xp.abs(tr))),
            "tilting_rms": float(trms),
            "tilting_relative_rms": float(trms / self.xp.maximum(tscale, 1e-30)),
            "omega_condition_index": float(
                self.xp.sum(source_rms) / self.xp.maximum(horizontal_scale, 1e-30)
            ),
            "tilting_condition_index": float(
                self.xp.sum(tilting_source_rms) / self.xp.maximum(tscale, 1e-30)
            ),
            "omega_advection_remainder_fraction": float(
                source_rms[remainder] / self.xp.maximum(horizontal_scale, 1e-30)
            ),
            "omega_advection_remainder_source_fraction": float(
                source_rms[remainder] / self.xp.maximum(self.xp.sum(source_rms), 1e-30)
            ),
            "tilting_advection_remainder_fraction": float(
                tilting_source_rms[remainder] / self.xp.maximum(tscale, 1e-30)
            ),
            "tilting_advection_remainder_source_fraction": float(
                tilting_source_rms[remainder]
                / self.xp.maximum(self.xp.sum(tilting_source_rms), 1e-30)
            ),
        }

    def _record_metrics(self, sim, stage):
        metrics = self._closure_metrics()
        metrics.update(time_s=float(sim.state.t), step=int(sim.step), stage=str(stage))
        self.history.append(metrics)

    def mark(self, sim, name, dt):
        current = self._native(sim)
        if name == "begin":
            if any(not bool(self.xp.array_equal(a, b))
                   for a, b in zip(current, self.previous_velocity)):
                raise RuntimeError("Unobserved velocity mutation in provenance tracer")
            if self.record_stage_metrics:
                self._record_metrics(sim, name)
            return

        omega_new = curl_native(current, self.grid)
        actual_delta = tuple(
            after - before for after, before in zip(current, self.previous_velocity)
        )

        if name == "advection":
            predicted_sum = tuple(
                self.xp.zeros_like(component) for component in self.previous_velocity
            )
            for source in range(len(self.source_names)):
                labelled = tuple(
                    component[source] for component in self.velocity_sources
                )
                tendencies = frozen_momentum_advection_tendency(
                    labelled,
                    self.previous_velocity,
                    self.grid,
                    order=self.momentum_order,
                )
                for component, tendency, accumulated in zip(
                    self.velocity_sources, tendencies, predicted_sum
                ):
                    increment = float(dt) * tendency
                    component[source] += increment
                    accumulated += increment
            remainder = self.source_index["advection_remainder"]
            for component, native, predicted in zip(
                self.velocity_sources, actual_delta, predicted_sum
            ):
                component[remainder] += native - predicted
        else:
            source_name = _STAGE_SOURCE.get(name, "other")
            source = self.source_index[source_name]
            for component, increment in zip(self.velocity_sources, actual_delta):
                component[source] += increment

        self.previous_velocity = current
        self.total_omega = omega_new
        if self.record_stage_metrics:
            self._record_metrics(sim, name)

        if name == "transport_microphysics_bcs" and float(sim.state.t) >= self.next_output:
            self.save(sim)
            while self.next_output <= float(sim.state.t):
                self.next_output += self.interval

    def _dataset(self, group, name, value):
        array = np.asarray(self.to_cpu(value), dtype=self.storage_dtype)
        return group.create_dataset(
            name, data=array, compression="gzip", compression_opts=1, shuffle=True
        )

    def save(self, sim, initial=False):
        if self.file is None:
            return
        group = self.snapshots.create_group(f"{self.count:05d}")
        group.attrs["time_s"] = float(sim.state.t)
        group.attrs["step"] = int(sim.step + (not initial))
        group.attrs["complete"] = True
        metrics = self._closure_metrics()
        for key, value in metrics.items():
            group.attrs[key] = value

        nk = self.output_nk
        self._dataset(group, "omega_total_h", self.total_omega[:2, :, :, :nk])
        _, _, _, (wx, wy) = velocity_kinematics(self.previous_velocity, self.grid)
        self._dataset(group, "grad_w_h", self.xp.stack((wx[:, :, :nk], wy[:, :, :nk])))
        sources = group.create_group("omega_source_h")
        omega_sources = self.omega_sources
        for name, index in self.source_index.items():
            self._dataset(sources, name, omega_sources[index, :2, :, :, :nk])
        if self.write_full_vorticity:
            self._dataset(group, "omega_total", self.total_omega[:, :, :, :nk])
            full_sources = group.create_group("omega_source")
            for name, index in self.source_index.items():
                self._dataset(
                    full_sources, name, omega_sources[index, :, :, :, :nk]
                )
        self.file.flush()
        self.last_saved_time = float(sim.state.t)
        self.count += 1

    def memory_estimate(self):
        """Return measured persistent bytes and raw bytes per written snapshot."""
        persistent = int(sum(value.nbytes for value in self.velocity_sources)
                         + self.total_omega.nbytes
                         + sum(a.nbytes for a in self.previous_velocity))
        cells_out = self.grid.nx * self.grid.ny * self.output_nk
        saved_scalars = 2 * len(self.source_names) + 2 + 2
        if self.write_full_vorticity:
            saved_scalars += 3 * len(self.source_names) + 3
        return {
            "persistent_bytes": persistent,
            "raw_snapshot_bytes": int(cells_out * saved_scalars * self.storage_dtype.itemsize),
            "source_count": len(self.source_names),
            "internal_cells": int(self.grid.nx * self.grid.ny * self.grid.nz),
            "output_cells": int(cells_out),
        }

    def write_full_checkpoint(self, path, sim):
        """Persist all source vectors needed to fork/restart provenance exactly.

        Periodic output stays restricted to the low-level horizontal components;
        a restart checkpoint is intentionally different because vertical import
        and subsequent component rotation require all three components aloft.
        """
        import h5py
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(path, "x") as handle:
            handle.attrs["schema"] = "storm-vorticity-provenance-checkpoint-v3"
            handle.attrs["time_s"] = float(sim.state.t)
            handle.attrs["step"] = int(sim.step)
            handle.attrs["source_names"] = json.dumps(self.source_names)
            handle.attrs["full_column"] = True
            grid_group = handle.create_group("grid")
            for name in ("xc", "yc", "zc"):
                grid_group.create_dataset(name, data=self.to_cpu(getattr(sim.grid, name)))
            sources = handle.create_group("velocity_source")
            for name, index in self.source_index.items():
                source = sources.create_group(name)
                for component_name, values in zip(
                    ("u", "v", "w"), self.velocity_sources
                ):
                    value = np.asarray(
                        self.to_cpu(values[index]), dtype=np.float64
                    )
                    source.create_dataset(
                        component_name,
                        data=value,
                        compression="gzip",
                        compression_opts=1,
                        shuffle=True,
                    )
            total = handle.create_group("velocity_total")
            for component_name, value in zip(("u", "v", "w"), self.previous_velocity):
                total.create_dataset(
                    component_name,
                    data=np.asarray(self.to_cpu(value), dtype=np.float64),
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                )

    def restore_full_checkpoint(self, path, sim, atol=1e-13, closure_atol=1e-10):
        """Restore source fields at an already-restored prognostic state."""
        import h5py
        path = Path(path)
        with h5py.File(path, "r") as handle:
            if handle.attrs.get("schema") != "storm-vorticity-provenance-checkpoint-v3":
                raise ValueError("unsupported provenance checkpoint schema")
            names = tuple(json.loads(handle.attrs["source_names"]))
            if names != self.source_names:
                raise ValueError("provenance source list mismatch")
            if abs(float(handle.attrs["time_s"]) - float(sim.state.t)) > 1e-9:
                raise ValueError("provenance/prognostic checkpoint time mismatch")
            for name in ("xc", "yc", "zc"):
                expected = np.asarray(self.to_cpu(getattr(sim.grid, name)))
                if not np.array_equal(expected, handle[f"grid/{name}"][:]):
                    raise ValueError(f"provenance checkpoint grid mismatch: {name}")
            restored = tuple(
                self.xp.stack(
                    [
                        self.xp.asarray(
                            handle[f"velocity_source/{name}/{component}"][:]
                        )
                        for name in names
                    ]
                )
                for component in ("u", "v", "w")
            )
            archived_total = tuple(
                self.xp.asarray(handle[f"velocity_total/{component}"][:])
                for component in ("u", "v", "w")
            )
        current_velocity = self._native(sim)
        mismatch = max(
            float(self.xp.max(self.xp.abs(current - archived)))
            for current, archived in zip(current_velocity, archived_total)
        )
        if mismatch > float(atol):
            raise ValueError(
                f"provenance/prognostic velocity mismatch {mismatch:.3e} > {atol:.3e}"
            )
        self.velocity_sources = restored
        self.previous_velocity = current_velocity
        self.total_omega = curl_native(current_velocity, self.grid)
        closure = self._closure_metrics()
        if closure["omega_max_abs"] > float(closure_atol):
            raise ValueError(
                "restored provenance does not close: "
                f"{closure['omega_max_abs']:.3e} > {closure_atol:.3e}"
            )
        return closure

    def close(self, sim, completed=True):
        if self.file is None:
            return
        if float(sim.state.t) > self.last_saved_time + 1e-12:
            # ``close`` is called after the outer driver has incremented
            # ``sim.step``.  Periodic saves occur inside the final solver stage
            # and need ``+1``; a close-time save already has the completed step.
            self.save(sim, initial=True)
        if self.history:
            dtype = np.dtype([
                ("time_s", "f8"), ("step", "i8"), ("stage", "S32"),
                ("omega_max_abs", "f8"), ("omega_rms", "f8"),
                ("omega_relative_rms", "f8"), ("tilting_max_abs", "f8"),
                ("tilting_rms", "f8"), ("tilting_relative_rms", "f8"),
                ("omega_condition_index", "f8"), ("tilting_condition_index", "f8"),
                ("omega_advection_remainder_fraction", "f8"),
                ("omega_advection_remainder_source_fraction", "f8"),
                ("tilting_advection_remainder_fraction", "f8"),
                ("tilting_advection_remainder_source_fraction", "f8"),
            ])
            rows = np.empty(len(self.history), dtype=dtype)
            for i, row in enumerate(self.history):
                rows[i] = tuple(row[name] for name in dtype.names)
            self.file.create_dataset("closure_history", data=rows, compression="gzip")
        self.file.attrs["status"] = "complete" if completed else "interrupted"
        self.file.flush()
        self.file.close()
        self.file = None


__all__ = [
    "SOURCE_NAMES", "DiagnosticObserverChain", "VorticityProvenanceTracer",
    "curl_native", "homogeneous_vorticity_tendency",
    "frozen_muscl_advection_tendency",
    "frozen_momentum_advection_tendency",
    "discrete_homogeneous_vorticity_tendency", "velocity_kinematics",
]
