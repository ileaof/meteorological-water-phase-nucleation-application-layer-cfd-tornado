"""Interpolation of the unified state onto the model mesh (ROADMAP §3a).

Horizontal: bilinear on the (rectilinear, projected) source grid.  Vertical: linear or a
mass-conservative (integral-preserving) remap for scalars.  No **silent** extrapolation --
points outside the source are clamped to the edge and the clamped fraction is recorded
(task requirements 5, 6, 9).  The state's ``z`` coordinate is taken as HEIGHT in metres
(the pressure-level readers convert via geopotential before this step).
"""
from __future__ import annotations

import numpy as np

_INTERP_LOG = []


def interpolation_log():
    """Return (and one may clear) the accumulated per-variable interpolation diagnostics."""
    return list(_INTERP_LOG)


def _rgi(points, values):
    from scipy.interpolate import RegularGridInterpolator
    return RegularGridInterpolator(points, values, method="linear",
                                   bounds_error=False, fill_value=None)  # None -> extrapolate; we clamp first


def horizontal_interp(field2d, x_src, y_src, x_tgt, y_tgt):
    """Bilinear interp of a (y,x) field from source axes to target 2-D meshes (clamped)."""
    xq = np.clip(x_tgt, x_src.min(), x_src.max())
    yq = np.clip(y_tgt, y_src.min(), y_src.max())
    f = _rgi((y_src, x_src), field2d)
    return f(np.stack([yq.ravel(), xq.ravel()], -1)).reshape(x_tgt.shape)


def vertical_remap(col, z_src, z_tgt, conservative=False):
    """Interpolate a 1-D column from ``z_src`` to ``z_tgt`` (both ascending, m).  Linear, or a
    layer-integral-preserving remap when ``conservative`` (for positive scalars like qv)."""
    z_src = np.asarray(z_src, float); col = np.asarray(col, float)
    zt = np.clip(z_tgt, z_src.min(), z_src.max())
    if not conservative:
        return np.interp(zt, z_src, col)
    # conservative: match the cumulative integral (mass) between source and target layers
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (col[1:] + col[:-1]) * np.diff(z_src))])
    cum_t = np.interp(zt, z_src, cum)
    out = np.empty_like(zt)
    dz = np.gradient(zt)
    out[:] = np.gradient(cum_t, zt)                       # d(integral)/dz recovers the density
    return np.where(np.isfinite(out), out, np.interp(zt, z_src, col))


def regrid_surface(state, name, x_model, y_model):
    """Bilinearly regrid a 2-D ``(y,x)`` surface field (e.g. terrain) onto the model axes."""
    ds = state.ds
    x_src = np.asarray(ds["x"].values, float); y_src = np.asarray(ds["y"].values, float)
    XT, YT = np.meshgrid(x_model, y_model, indexing="xy")
    return horizontal_interp(np.asarray(ds[name].values, float), x_src, y_src, XT, YT)


def regrid_to_model(state, x_model, y_model, z_model, variables=None, conservative=True):
    """Regrid an :class:`AtmosphericState` (time,z,y,x native) onto the model mesh.

    ``x_model, y_model`` are 1-D model axes in the SAME projected frame as the state (centred
    at the domain centre); ``z_model`` model heights [m].  Returns ``{var: (time,nz,ny,nx)}``
    plus an entry ``'_log'`` of interpolation diagnostics.  Extrapolated points are clamped and
    the fraction recorded -- never silent."""
    ds = state.ds
    x_src = np.asarray(ds["x"].values, float); y_src = np.asarray(ds["y"].values, float)
    z_src = np.asarray(ds["z"].values, float)
    XT, YT = np.meshgrid(x_model, y_model, indexing="xy")   # (ny,nx)
    nt = ds.sizes["time"]; nzt = z_model.size; nyt, nxt = XT.shape
    names = variables or [v for v in ds.data_vars if ds[v].ndim == 4]
    out = {}
    clip_frac = float(np.mean((x_model < x_src.min()) | (x_model > x_src.max())) +
                      np.mean((y_model < y_src.min()) | (y_model > y_src.max()))) / 2.0
    for nm in names:
        src = np.asarray(ds[nm].values, float)              # (t,z,y,x)
        res = np.empty((nt, nzt, nyt, nxt))
        col_int_err = 0.0
        for it in range(nt):
            # horizontal at each source level, then vertical remap per column
            horiz = np.empty((z_src.size, nyt, nxt))
            for k in range(z_src.size):
                horiz[k] = horizontal_interp(src[it, k], x_src, y_src, XT, YT)
            cons = conservative and nm in ("qv", "qc", "qr", "qi", "qs", "qg")
            for j in range(nyt):
                for i in range(nxt):
                    res[it, :, j, i] = vertical_remap(horiz[:, j, i], z_src, z_model, conservative=cons)
            if cons:                                        # track column-integral (mass) drift
                _trapz = getattr(np, "trapezoid", np.trapz)
                src_int = _trapz(horiz, z_src, axis=0)
                tgt_int = _trapz(res[it], z_model, axis=0)
                col_int_err = max(col_int_err, float(np.nanmax(np.abs(tgt_int - src_int)
                                                               / (np.abs(src_int) + 1e-12))))
        out[nm] = res
        _INTERP_LOG.append({"variable": nm, "clamped_fraction": round(clip_frac, 4),
                            "column_integral_rel_err": round(col_int_err, 5),
                            "method": "bilinear+%s" % ("conservative-z" if nm.startswith("q") else "linear-z")})
    out["_log"] = interpolation_log()
    return out
