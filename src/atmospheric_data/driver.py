"""real_case driver: wire ingested real data into the storm CFD (ROADMAP §3a).

Ties the pipeline together WITHOUT changing the idealized mode:
    preprocess  -> load atmosphere, regrid to the model mesh, base state, IC/BC/surface, QC
    build_simulation -> a StormSimulation initialised from the real base + initial fields
    run_case    -> step it (CPU/GPU auto, with the Davies lateral relaxation toward the
                   real environment via storm_dynamics.limited_area)
    compare_radar -> project the CFD (u,v,w) to synthetic radial velocity vs NEXRAD

Backend: ``execution_backend: auto`` uses the GPU when present and **falls back to CPU**
automatically (never requires CUDA).
"""
from __future__ import annotations

import os

import numpy as np

from . import interpolate, basestate, ic_bc, qc, radial, validation
from .sources import load_atmosphere, load_radar


def _grid_dims(cfg, max_n=64):
    """Parent-grid (nx,ny,nz, Lx,Ly,Lz) from the domain + parent_dx (capped at ``max_n`` per
    axis so a default 400 km / 1.3 km case stays runnable in tests; real runs raise the cap)."""
    Lx = cfg.domain.width_km * 1000.0; Lz = cfg.domain.height_km * 1000.0
    nx = int(round(Lx / cfg.model.parent_dx_m))
    nz = max(16, int(round(Lz / (cfg.model.parent_dx_m))))          # coarse z by default
    nx = min(nx, max_n); nz = min(nz, max_n)
    return nx, nx, nz, Lx, Lx, Lz


def _make_grid(cfg, max_n=64):
    from meteorological_flow.grid import Grid
    from meteorological_flow.backend import get_backend
    nx, ny, nz, Lx, Ly, Lz = _grid_dims(cfg, max_n)
    dev = {"auto": "auto", "cpu": "cpu", "gpu": "gpu"}[cfg.model.execution_backend]
    backend = get_backend(dev)                                      # 'auto'/'gpu' -> CPU fallback inside
    return Grid(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz, z_stretch=1.03,
                periodic=True, backend=backend)


def preprocess(cfg, cache, outdir, logger=print, max_n=64):
    """Load -> regrid -> base state -> IC/BC/surface -> QC.  Returns a dict of artefacts."""
    os.makedirs(outdir, exist_ok=True)
    state = load_atmosphere(cfg, cache, logger=logger)
    grid = _make_grid(cfg, max_n)
    xc = np.asarray(grid.backend.to_cpu(grid.xc), float) - 0.5 * grid.Lx    # centre on the domain
    yc = np.asarray(grid.backend.to_cpu(grid.yc), float) - 0.5 * grid.Ly
    zc = np.asarray(grid.backend.to_cpu(grid.zc), float)
    logger("[preprocess] regridding to %dx%dx%d (parent dx=%.0f m) ..." % (grid.nx, grid.ny, grid.nz, grid.dx))
    fields = interpolate.regrid_to_model(state, xc, yc, zc,
                                         conservative=(cfg.processing.interpolation == "conservative"))
    base = basestate.base_state_from_fields(fields, zc)
    ic = ic_bc.write_initial_conditions(fields, xc, yc, zc, os.path.join(outdir, "initial_conditions.nc"),
                                        meta={"source": state.ds.attrs.get("source", "")})
    bcs = ic_bc.write_boundaries(fields, xc, yc, zc, state.ds["time"].values, outdir)
    terr = interpolate.regrid_surface(state, "terrain", xc, yc) if state.has("terrain") else None
    sf = ic_bc.write_surface_forcing(terr, xc, yc, os.path.join(outdir, "surface_forcing.nc"))
    report = qc.quality_control(state)
    qcj, qcm = qc.write_reports(report, os.path.join(outdir, "qc_report.json"),
                                os.path.join(outdir, "qc_report.md"))
    logger("[preprocess] QC %d/%d passed (%s); interp: %s"
           % (report["summary"]["passed"], report["summary"]["total"],
              "OK" if report["summary"]["ok"] else "REVIEW", fields["_log"][:1]))
    return {"state": state, "grid": grid, "fields": fields, "base": base,
            "initial_conditions": ic, "boundaries": bcs, "surface_forcing": sf,
            "qc": report, "qc_json": qcj, "qc_md": qcm, "interp_log": fields["_log"],
            "coords": (xc, yc, zc)}


def build_simulation(cfg, pre, logger=print):
    """Build a StormSimulation initialised from the real base + initial fields (staggered)."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    grid = pre["grid"]; base = pre["base"]; fields = pre["fields"]
    dev = {"auto": "auto", "cpu": "cpu", "gpu": "gpu"}[cfg.model.execution_backend]
    scfg = build_storm_config(preset="storm", nx=grid.nx, ny=grid.ny, nz=grid.nz,
                              Lx=grid.Lx, Ly=grid.Ly, Lz=grid.Lz, duration=1.0, dt_max=3.0,
                              drag=True, z_stretch=1.03, device=dev)
    scfg.sim.physics.bubble_dtheta = 0.0                            # real IC, not a warm bubble
    sim = StormSimulation(scfg, base=base)
    xp = sim.grid.xp
    tr = lambda a: np.transpose(np.asarray(a, float), (2, 1, 0))    # (z,y,x)->(x,y,z)
    th = tr(fields["theta"][0]); qv = tr(fields["qv"][0]) if "qv" in fields else None
    uc = tr(fields["u"][0]); vc = tr(fields["v"][0]); wc = tr(fields["w"][0])
    sim.state.theta = xp.asarray(th)
    if qv is not None:
        sim.state.qv = xp.asarray(np.clip(qv, 0.0, None))
    # cell-centred -> staggered faces (average neighbours; walls/edges clamped)
    u = np.zeros(sim.grid.u_shape); u[1:-1] = 0.5 * (uc[:-1] + uc[1:]); u[0] = uc[0]; u[-1] = uc[-1]
    v = np.zeros(sim.grid.v_shape); v[:, 1:-1] = 0.5 * (vc[:, :-1] + vc[:, 1:]); v[:, 0] = vc[:, 0]; v[:, -1] = vc[:, -1]
    w = np.zeros(sim.grid.w_shape); w[:, :, 1:-1] = 0.5 * (wc[:, :, :-1] + wc[:, :, 1:])
    sim.state.u = xp.asarray(u); sim.state.v = xp.asarray(v); sim.state.w = xp.asarray(w)
    sim.state.diagnose(sim.cfg)
    # lateral relaxation target = the environment (Davies zone; time-dependent driving optional)
    from storm_dynamics.limited_area import environment_target
    sim._lbc_target = environment_target(sim.grid, base)
    logger("[build] StormSimulation %dx%dx%d on %s, initialised from real fields"
           % (grid.nx, grid.ny, grid.nz, sim.grid.backend.name))
    return sim


def run_case(cfg, pre, sim=None, steps=None, logger=print, lbc_width=8, lbc_rate=1.0 / 300.0):
    """Step the real-case simulation with the Davies lateral relaxation each step."""
    from storm_dynamics.limited_area import apply_lateral_relaxation, lateral_relaxation_weight
    sim = sim or build_simulation(cfg, pre, logger=logger)
    n = steps or 5
    w = lateral_relaxation_weight(sim.grid, lbc_width, lbc_rate)
    for s in range(n):
        dt = float(sim._dt())
        sim._step(dt)
        apply_lateral_relaxation(sim.state, sim.grid, sim._lbc_target, dt, weight=w)
        sim.step += 1; sim.t = float(getattr(sim.state, "t", sim.t + dt))
    logger("[run] stepped %d steps; final t=%.1f s" % (n, sim.t))
    return sim


def run_multilevel_real_case(cfg, pre, window=None, mature_steps=0, half_frac=0.28,
                             sim=None, logger=print):
    """Drive the AMR **multi-level nest cascade** from the real initial conditions
    (parent Δx → nest → fine Δx), reusing ``storm_dynamics.run_multilevel_nest``.

    The nests are centred on the domain centre (where the real_case domain places the storm)
    with refinement factors from the config (``parent_dx_m → nest_dx_m → fine_dx_m``).  Ground
    frame; the storm-motion vector (``pre['fields']`` / the alert) can drive a Galilean frame
    later.  Returns ``(sims, rep)`` with ``sims[-1]`` the finest level."""
    from storm_dynamics import nesting as nst
    parent = sim or build_simulation(cfg, pre, logger=logger)
    for _ in range(int(mature_steps)):                             # optional short spin-up
        parent._step(float(parent._dt()))
    r1 = max(2, int(round(cfg.model.parent_dx_m / cfg.model.nest_dx_m)))
    r2 = max(2, int(round(cfg.model.nest_dx_m / cfg.model.fine_dx_m)))
    mkspec = lambda refine: (lambda g: nst.NestSpec.around(
        g, 0.5 * g.Lx, 0.5 * g.Ly, half=g.Lx * half_frac, refine=refine, nz=g.nz,
        z_stretch=getattr(g, "z_stretch", 1.0)))
    window = window or 20.0 * float(parent._dt())
    logger("[multilevel] cascade parent dx=%.0f -> nest /%d -> fine /%d over %.0f s"
           % (parent.grid.dx, r1, r2, window))
    sims, rep = nst.run_multilevel_nest(parent, [mkspec(r1), mkspec(r2)], window=window,
                                        les_boost=1.4, cfl=cfg_cfl(cfg))
    finest = sims[-1]
    rep.setdefault("nest", {})
    logger("[multilevel] finest dx=%.0f m, zeta_abs_max=%.3e"
           % (finest.grid.dx, rep["rotation"]["zeta_abs_max"]))
    return sims, rep


def cfg_cfl(cfg):
    return 0.2


def compare_radar(cfg, pre, cache, sim=None, logger=print):
    """Project the CFD (u,v,w) to synthetic radial velocity and score against NEXRAD."""
    radar = load_radar(cfg, cache, logger=logger)
    xc, yc, zc = pre["coords"]
    fields = pre["fields"]
    if sim is not None:                                            # use the simulated velocity
        xp = sim.grid.xp; to = sim.grid.backend.to_cpu
        uc = np.asarray(to(sim.state.u)); vc = np.asarray(to(sim.state.v)); wc = np.asarray(to(sim.state.w))
        ucc = 0.5 * (uc[:-1] + uc[1:]); vcc = 0.5 * (vc[:, :-1] + vc[:, 1:]); wcc = 0.5 * (wc[:, :, :-1] + wc[:, :, 1:])
        tr = lambda a: np.transpose(a, (2, 1, 0))[None]           # (x,y,z)->(1,z,y,x)
        fields = {"u": tr(ucc), "v": tr(vcc), "w": tr(wcc)}
    vr = radial.cfd_radial_velocity(fields, xc, yc, zc, radar)
    metrics = validation.radar_metrics(radar, vr)
    logger("[radar] Vr RMSE=%.2f corr=%.2f | mesocyclone displ=%s gates"
           % (metrics["radial_velocity"]["rmse"], metrics["radial_velocity"]["correlation"],
              metrics.get("mesocyclone_displacement_gates")))
    return {"radar": radar, "vr_sim": vr, "metrics": metrics}
