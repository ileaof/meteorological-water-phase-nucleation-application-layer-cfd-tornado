"""Rotation-focused visualisation for storm_dynamics.

These plots are the *visual* counterpart of :mod:`storm_dynamics.rotation`'s
numbers -- the pictures that show a storm has rotated:

* :func:`plot_rotation_slices` -- horizontal maps of the mid-level vertical
  vorticity (the mesocyclone / storm-split couplet), the vertical velocity, and
  the near-surface vertical vorticity, with the horizontal perturbation wind
  overlaid;
* :func:`plot_hodograph` -- the environmental hodograph with the 0-1 / 0-3 km SRH
  and bulk shear annotated;
* :func:`plot_rotation_timeseries` -- peak near-surface zeta, mid-level
  mesocyclone, updraft helicity and w_max vs. time.

Matplotlib is used with the non-interactive Agg backend (no display needed); each
function writes a PNG and returns its path.  These are demonstration diagnostics
at coarse resolution -- see docs/storm_dynamics_guide.md.
"""
from __future__ import annotations

import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

from . import rotation as rot     # noqa: E402


def _to_cpu(grid, a):
    return grid.backend.to_cpu(a)


def _level_index(grid, z_target):
    z = _to_cpu(grid, grid.zc)
    return int(np.argmin(np.abs(z - z_target)))


def _perturbation_wind_centers(sim):
    """Cell-centre horizontal perturbation wind (u-u0, v-v0)."""
    g = sim.grid
    up = sim.state.u - sim._u0_face
    vp = sim.state.v - sim._v0_face
    uc = 0.5 * (up[:-1, :, :] + up[1:, :, :])
    vc = 0.5 * (vp[:, :-1, :] + vp[:, 1:, :])
    return _to_cpu(g, uc), _to_cpu(g, vc)


def plot_rotation_slices(sim, outdir, z_mid=4000.0, z_near=500.0, tag="") -> str:
    """Horizontal maps of mid-level zeta, w, and near-surface zeta (PNG)."""
    os.makedirs(outdir, exist_ok=True)
    g = sim.grid
    _, _, zeta = rot.vorticity_3d(sim.state, g)
    _, _, wc = rot._centered_velocity(sim.state, g)
    zeta = _to_cpu(g, zeta); wc = _to_cpu(g, wc)
    uc, vc = _perturbation_wind_centers(sim)
    x = _to_cpu(g, g.xc) / 1000.0    # km
    y = _to_cpu(g, g.yc) / 1000.0
    kmid = _level_index(g, z_mid)
    knear = _level_index(g, z_near)
    zmid_km = _to_cpu(g, g.zc)[kmid] / 1000.0
    znear_km = _to_cpu(g, g.zc)[knear] / 1000.0

    def _mesh(field_k):
        return field_k.T   # (ny,nx) for pcolormesh(x,y)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), constrained_layout=True)
    # step for quiver so arrows are readable
    sx = max(1, g.nx // 20); sy = max(1, g.ny // 20)

    # panel 1: mid-level zeta (the mesocyclone / split couplet)
    ax = axes[0]
    zmax = max(1e-6, float(np.abs(zeta[:, :, kmid]).max()))
    pc = ax.pcolormesh(x, y, _mesh(zeta[:, :, kmid]), cmap="RdBu_r",
                       vmin=-zmax, vmax=zmax, shading="auto")
    ax.contour(x, y, _mesh(wc[:, :, kmid]), levels=[2, 5, 10], colors="k",
               linewidths=[0.6, 0.9, 1.2])
    ax.quiver(x[::sx], y[::sy], uc[::sx, ::sy, kmid].T, vc[::sx, ::sy, kmid].T,
              scale=250, width=0.003, color="0.25")
    fig.colorbar(pc, ax=ax, label=r"$\zeta$ [s$^{-1}$]")
    ax.set_title("mid-level $\\zeta$ @ %.1f km\n(red=cyclonic, blue=anticyclonic; "
                 "black=updraft w)" % zmid_km)
    ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]")

    # panel 2: vertical velocity at mid-level
    ax = axes[1]
    wmax = max(1e-3, float(np.abs(wc[:, :, kmid]).max()))
    pc = ax.pcolormesh(x, y, _mesh(wc[:, :, kmid]), cmap="RdBu_r",
                       vmin=-wmax, vmax=wmax, shading="auto")
    fig.colorbar(pc, ax=ax, label="w [m/s]")
    ax.set_title("vertical velocity w @ %.1f km" % zmid_km)
    ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]")

    # panel 3: near-surface zeta (tornadogenesis proxy)
    ax = axes[2]
    zmax = max(1e-6, float(np.abs(zeta[:, :, knear]).max()))
    pc = ax.pcolormesh(x, y, _mesh(zeta[:, :, knear]), cmap="RdBu_r",
                       vmin=-zmax, vmax=zmax, shading="auto")
    ax.quiver(x[::sx], y[::sy], uc[::sx, ::sy, knear].T, vc[::sx, ::sy, knear].T,
              scale=200, width=0.003, color="0.25")
    fig.colorbar(pc, ax=ax, label=r"$\zeta$ [s$^{-1}$]")
    ax.set_title("near-surface $\\zeta$ @ %.2f km\n(tornadogenesis proxy)" % znear_km)
    ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]")

    fig.suptitle("storm_dynamics rotation slices  (t = %.0f s)  -- IDEALISED, not a forecast"
                 % sim.t, fontsize=13)
    path = os.path.join(outdir, "rotation_slices%s.png" % (("_" + tag) if tag else ""))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def plot_hodograph(base, outdir, tag="") -> str:
    """Environmental hodograph (u vs v) with SRH / shear annotated (PNG)."""
    from . import soundings as snd
    os.makedirs(outdir, exist_ok=True)
    z = np.asarray(base.zc); u = np.asarray(base.u0); v = np.asarray(base.v0)
    sel = z <= 12000.0
    fig, ax = plt.subplots(figsize=(6.2, 6.0), constrained_layout=True)
    sc = ax.scatter(u[sel], v[sel], c=z[sel] / 1000.0, cmap="viridis", s=18, zorder=3)
    ax.plot(u[sel], v[sel], "-", color="0.5", lw=1, zorder=2)
    cx, cy = snd.bunkers_storm_motion(base)
    ax.plot(cx, cy, "r*", ms=16, label="storm motion (right-mover)", zorder=4)
    for h in (1000.0, 3000.0, 6000.0):
        k = int(np.argmin(np.abs(z - h)))
        ax.annotate("%.0f km" % (h / 1000.0), (u[k], v[k]),
                    textcoords="offset points", xytext=(6, 4), fontsize=9)
    fig.colorbar(sc, ax=ax, label="height [km]")
    ax.axhline(0, color="0.8", lw=0.8); ax.axvline(0, color="0.8", lw=0.8)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("u [m/s]"); ax.set_ylabel("v [m/s]")
    ax.set_title("environmental hodograph\nSRH$_{0-1}$=%.0f  SRH$_{0-3}$=%.0f m$^2$/s$^2$  |  "
                 "shear$_{0-6}$=%.0f m/s"
                 % (snd.storm_relative_helicity(base, 1000.0),
                    snd.storm_relative_helicity(base, 3000.0),
                    snd.bulk_shear(base, 0, 6000)))
    ax.legend(loc="best", fontsize=9)
    path = os.path.join(outdir, "hodograph%s.png" % (("_" + tag) if tag else ""))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def plot_rotation_timeseries(sim, outdir, tag="") -> str:
    """Peak near-surface zeta, mesocyclone, updraft helicity and w vs time (PNG)."""
    os.makedirs(outdir, exist_ok=True)
    tr = sim.tracker
    t = np.asarray(tr.t) / 60.0    # minutes
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    axes[0, 0].plot(t, tr.near_surface_zeta_max, "C3")
    axes[0, 0].set_title("near-surface $\\zeta$ max (tornadogenesis proxy)")
    axes[0, 0].set_ylabel(r"$\zeta$ [s$^{-1}$]")
    axes[0, 1].plot(t, tr.midlevel_meso, "C0")
    axes[0, 1].set_title("mid-level mesocyclone $|\\zeta|$")
    axes[0, 1].set_ylabel(r"$\zeta$ [s$^{-1}$]")
    axes[1, 0].plot(t, tr.updraft_helicity_max, "C2")
    axes[1, 0].set_title("updraft helicity max")
    axes[1, 0].set_ylabel(r"UH [m$^2$/s$^2$]")
    axes[1, 1].plot(t, tr.w_max, "C1")
    axes[1, 1].set_title("w max")
    axes[1, 1].set_ylabel("w [m/s]")
    for ax in axes.flat:
        ax.set_xlabel("time [min]"); ax.grid(alpha=0.3)
    fig.suptitle("storm_dynamics rotation time series -- IDEALISED, not a forecast",
                 fontsize=13)
    path = os.path.join(outdir, "rotation_timeseries%s.png" % (("_" + tag) if tag else ""))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def plot_all(sim, outdir, tag="") -> dict:
    """Write the full rotation figure set; return {name: path}."""
    return {
        "hodograph": plot_hodograph(sim.base, outdir, tag=tag),
        "slices": plot_rotation_slices(sim, outdir, tag=tag),
        "timeseries": plot_rotation_timeseries(sim, outdir, tag=tag),
    }


__all__ = [
    "plot_rotation_slices", "plot_hodograph", "plot_rotation_timeseries", "plot_all",
]
