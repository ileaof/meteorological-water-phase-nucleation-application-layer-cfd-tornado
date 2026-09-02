"""Diagnostic figure suite for tornadogenesis -- renders the new instruments from a live run:
the vertical-vorticity budget, the low-level vorticity, the cold pool, and the resolved (macro) vs
sub-grid (micro) temperature gradient with the w-vs-gradient relations.

    python examples/tornado_diagnostic_figures.py           # short sim -> docs/media/storm/tornado_diagnostics.png
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from storm_dynamics import vorticity_budget as vb, coldpool as cp, micro_gradient as mg
from storm_dynamics.rotation import _centered_velocity, vertical_vorticity

_TERMS = ("advection", "stretching", "tilting", "baroclinic", "divergence", "diffusion")


def make_diagnostic_figure(sim, out_path, z_low_m=500.0):
    g = sim.grid; to = g.backend.to_cpu
    sim.state.diagnose(sim.cfg)
    k = int(np.argmin(np.abs(np.asarray(to(g.zc)) - z_low_m)))
    terms = vb.zeta_budget(sim.state, g, Km=getattr(sim, "_Km", None))
    summ = vb.budget_layer_summary(terms, g, 0.0, 1000.0)
    _, _, wc = _centered_velocity(sim.state, g); wc = np.asarray(to(wc))
    zeta = np.asarray(to(vertical_vorticity(sim.state, g)))
    _, tvp = cp.coldpool_buoyancy(sim.state, g); tvp = np.asarray(to(tvp))
    macro = np.asarray(to(mg.macro_temperature_gradient(sim.state, g)))
    micro, _ = mg.micro_temperature_gradient(sim.state, g); micro = np.asarray(to(micro))
    xc = np.asarray(to(g.xc)) / 1e3; yc = np.asarray(to(g.yc)) / 1e3

    fig, ax = plt.subplots(2, 3, figsize=(15, 9), facecolor="white")

    # (a) vorticity-budget layer means (0-1 km)
    vals = [summ[t + "_absmean"] for t in _TERMS]
    ax[0, 0].bar(range(len(_TERMS)), vals, color="#2f6690")
    ax[0, 0].set_xticks(range(len(_TERMS))); ax[0, 0].set_xticklabels(_TERMS, rotation=35, ha="right", fontsize=8)
    ax[0, 0].set_ylabel("|term| mean  [s$^{-2}$]"); ax[0, 0].set_title("(a) 0-1 km vorticity budget", loc="left", fontweight="bold")

    # (b) low-level vertical vorticity slice
    zm = float(np.abs(zeta[:, :, k]).max()) or 1e-6
    im = ax[0, 1].pcolormesh(xc, yc, zeta[:, :, k].T, cmap="RdBu_r", vmin=-zm, vmax=zm, shading="auto")
    ax[0, 1].set_aspect("equal"); ax[0, 1].set_title("(b) low-level $\\zeta$ [s$^{-1}$]", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax[0, 1], fraction=0.046)

    # (c) cold-pool theta_v' slice
    tvm = float(np.abs(tvp[:, :, k]).max()) or 1e-6
    im = ax[0, 2].pcolormesh(xc, yc, tvp[:, :, k].T, cmap="Blues_r", vmin=-tvm, vmax=0.0, shading="auto")
    ax[0, 2].set_aspect("equal"); ax[0, 2].set_title("(c) cold pool $\\theta_v'$ [K]", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax[0, 2], fraction=0.046)

    # (d) w vs macro gradient, (e) w vs micro gradient
    wf = wc.ravel(); maf = macro.ravel(); mif = micro.ravel()
    sel = np.random.default_rng(0).choice(wf.size, size=min(4000, wf.size), replace=False)
    ax[1, 0].scatter(maf[sel], wf[sel], s=3, alpha=0.3, color="#2f6690")
    ax[1, 0].set_xlabel("|$\\nabla T$|$_{macro}$ [K/m]"); ax[1, 0].set_ylabel("w [m/s]")
    ax[1, 0].set_title("(d) w vs macro gradient", loc="left", fontweight="bold")
    ax[1, 1].scatter(mif[sel], wf[sel], s=3, alpha=0.3, color="#b8143c")
    ax[1, 1].set_xlabel("|$\\nabla T$|$_{micro}$ [K/m]"); ax[1, 1].set_ylabel("w [m/s]")
    ax[1, 1].set_title("(e) w vs micro gradient", loc="left", fontweight="bold")
    ax[1, 1].set_xscale("log")

    # (f) macro vs micro gradient
    ax[1, 2].scatter(maf[sel], mif[sel], s=3, alpha=0.3, color="#5b8c5a")
    ax[1, 2].set_xlabel("|$\\nabla T$|$_{macro}$"); ax[1, 2].set_ylabel("|$\\nabla T$|$_{micro}$")
    ax[1, 2].set_yscale("log"); ax[1, 2].set_title("(f) micro vs macro gradient", loc="left", fontweight="bold")

    fig.suptitle("Tornadogenesis diagnostics -- vorticity budget, cold pool, and the two-scale "
                 "temperature gradient (all measured, none imposed)", fontsize=12, y=1.0)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def main():
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm", nx=48, ny=48, nz=40, Lx=48000.0, Ly=48000.0,
                              Lz=15000.0, duration=1.0, dt_max=4.0, z_stretch=1.05, device="cpu")
    scfg.sim.physics.bubble_dtheta = 5.0
    sim = StormSimulation(scfg)
    for _ in range(200):
        sim._step(sim._dt()); sim.step += 1; sim.t = float(sim.state.t)
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "media", "storm", "tornado_diagnostics.png")
    print("saved", make_diagnostic_figure(sim, out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
