"""Attempt F -- ALTERNATIVE 1: resolve the UPDRAFT during maturation.

Attempts C-E kept the storm-generating parent coarse (750 m-1.5 km) and only refined the nest at
the end, so the parent updraft stayed weak (w ~6-14).  The old idealised model reached w ~25-56 m/s
by *resolving the updraft* (convergence study M8: convection is grid-dependent, needs dx <~ 250 m,
Bryan 2003).  So here the REAL KOUN storm is matured on a FINE 250 m parent from the start, then
cascaded to 28 m.  Question: does a resolved, strong updraft (more vertical stretching) break the
~6 m/s low-level-rotation ceiling that resolution (D) and SRH (E) could not?"""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "moore_fineparent"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress.log"); open(LOG, "w").close()
_t0 = time.time()
def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m)); print(m, flush=True)

import atmospheric_data as ad
from atmospheric_data.sources import iem_raob
from meteorological_flow.grid import Grid
from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config, MesoForcingConfig
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, nesting as nst, soundings as snd

def fair_low_vrot(u, v, zeta, zc, dx, zmax_m=1500.0):
    """Max low-level rotational velocity over z<zmax (same window everywhere)."""
    nb = zeta.shape[0] // 6; best = (0.0, 0.0, 0.0)
    for k in range(zeta.shape[2]):
        if zc[k] > zmax_m:
            break
        Z = np.abs(zeta[nb:-nb, nb:-nb, k]); i, j = np.unravel_index(np.argmax(Z), Z.shape); i += nb; j += nb
        R = max(4, int(1000 / dx))
        us = u[max(0, i - R):i + R, max(0, j - R):j + R, k]; vs = v[max(0, i - R):i + R, max(0, j - R):j + R, k]
        vr = float(np.sqrt((us - us.mean()) ** 2 + (vs - vs.mean()) ** 2).max())
        if vr > best[0]:
            best = (vr, float(zc[k]), float(np.abs(zeta[nb:-nb, nb:-nb, k]).max()))
    return best

DEV = os.environ.get("DEV", "gpu")
log("=== Attempt F: FINE 250 m parent (resolve the updraft), REAL KOUN, storm-relative, %s ===" % DEV)
cache = ad.Cache(os.path.join(REPO, "data", "cache"))
prof = iem_raob.download_sounding("KOUN", "2013-05-21T00:00:00Z", cache=cache)
nx = 120; nz = 48; Lx = 30000.0; Lz = 15000.0
g = Grid(nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, z_stretch=1.05, periodic=True); zc = np.asarray(g.zc)
base_ground = snd.from_observed_sounding(
    g, pressure_hPa=prof["pressure_Pa"] / 100.0, height_m=prof["height_m"],
    temperature_C=prof["temperature_K"] - 273.15, dewpoint_C=prof["dewpoint_K"] - 273.15,
    u_ms=prof["u_ms"], v_ms=prof["v_ms"])
d = sounding_diagnostics(base_ground); cx, cy = snd.bunkers_storm_motion(base_ground)
log("REAL KOUN env: CAPE=%.0f CIN=%.0f shear06=%.1f SRH03=%.0f | fine parent dx=%.0f m (%dx%dx%d)" %
    (d["CAPE_J_kg"], d["CIN_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(base_ground),
     g.dx, nx, nx, nz))
base = BaseState(zc=zc, theta0=base_ground.theta0, qv0=base_ground.qv0, p0=base_ground.p0,
                 T0=base_ground.T0, rho0=base_ground.rho0, u0=base_ground.u0 - cx, v0=base_ground.v0 - cy)

FORCE_S = 1500.0
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, duration=1.0,
                          dt_max=2.0, drag=True, z_stretch=1.05, C_s=0.18, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
scfg.dyn.forcing = MesoForcingConfig(enabled=True, heat_rate_K_s=0.008, moist_rate_kgkg_s=4e-6,
                                     radius_m=6000.0, z_top_m=2800.0, duration_s=FORCE_S)
sim = StormSimulation(scfg, base=base)
log("maturing on the FINE parent (forced first %.0fs) -- watching if w_max reaches ~25-40 ..." % FORCE_S)
T_MAT = 1650.0
step = 0; wmax_peak = 0.0
while sim.t < T_MAT:
    dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t); step += 1
    if step % 60 == 0:
        to = sim.grid.backend.to_cpu; w = np.asarray(to(sim.state.w))
        if not np.isfinite(w).all():
            log("  instability at t=%.0f -> stop" % sim.t); break
        wmax_peak = max(wmax_peak, float(w.max()))
        log("  parent t=%4.0f w_max=%5.1f (peak %5.1f)" % (sim.t, w.max(), wmax_peak))

to = sim.grid.backend.to_cpu
uc, vc, wc = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wc = np.asarray(to(wc)); zcp = np.asarray(to(sim.grid.zc))
zeta_p = np.gradient(vc, sim.grid.dx, axis=0) - np.gradient(uc, sim.grid.dy, axis=1)
pvr, pz, pzt = fair_low_vrot(uc, vc, zeta_p, zcp, sim.grid.dx)
kmeso = int(np.argmin(np.abs(zcp - 700))); nbp = nx // 6
zp7 = np.abs(zeta_p[nbp:-nbp, nbp:-nbp, kmeso])
im, jm = np.unravel_index(np.argmax(zp7), zp7.shape); im += nbp; jm += nbp
mx, my = float(sim.grid.xc[im]), float(sim.grid.yc[jm])
log("PARENT(250m) matured: w_max_peak=%.1f  low-lvl Vrot(fair)=%.1f at z=%.0f  meso@(%.1f,%.1f)km" %
    (wmax_peak, pvr, pz, mx / 1e3, my / 1e3))

# cascade 250 -> 83 -> 28 m (2x refine-3) centred on the low-level meso
fr = [0.10, 0.13]
def mkspec(i):
    def build(gg):
        c = (mx, my) if i == 0 else (0.5 * gg.Lx, 0.5 * gg.Ly)
        return nst.NestSpec.around(gg, c[0], c[1], half=gg.Lx * fr[i], refine=3,
                                   nz=gg.nz, z_stretch=getattr(gg, "z_stretch", 1.05))
    return build
WIN = 180.0
log("cascade 250 -> 83 -> ~28 m (refine 3 x2), window %.0f s ..." % WIN)
sims, rep = nst.run_multilevel_nest(sim, [mkspec(0), mkspec(1)], window=WIN, les_boost=1.5, cfl=0.2,
                                    restrict_momentum=True,
                                    progress=lambda t, w, s: log("  finest t=%5.1f/%.0f step=%d" % (t, w, s))
                                    if s % 80 == 0 else None)
finest = sims[-1]; ng = finest.grid; to = ng.backend.to_cpu
uc, vc, wcn = rot._centered_velocity(finest.state, ng)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wcn)); zcn = np.asarray(to(ng.zc))
zeta = np.gradient(vc, ng.dx, axis=0) - np.gradient(uc, ng.dy, axis=1)
fvr, fz, fzt = fair_low_vrot(uc, vc, zeta, zcn, ng.dx)
np.savez_compressed(os.path.join(OUT, "fields_cascade.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    xc=np.asarray(to(ng.xc)), yc=np.asarray(to(ng.yc)), zc=zcn, dx=float(ng.dx))
summary = {"attempt": "F fine-250m-parent", "source": "REAL KOUN 2013-05-21 00Z",
           "env": {k: float(d[k]) for k in ("CAPE_J_kg", "CIN_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(base_ground)), "parent_dx_m": float(sim.grid.dx),
           "parent_wmax_peak": wmax_peak, "parent_low_vrot": pvr,
           "levels_dx_m": [float(s.grid.dx) for s in sims], "finest_dx_m": float(ng.dx), "finest_nx": int(ng.nx),
           "finest_low_vrot_fair": fvr, "finest_low_vrot_z": fz, "finest_low_peak_zeta": fzt,
           "compare": {"D_ERA5_vrot": 6.7, "E_KOUN750_vrot": 6.0, "obs_vrot": 26.0},
           "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
log("=== DONE F: parent_dx=%.0fm w_max_peak=%.1f | finest=%.0fm low-lvl Vrot(fair)=%.1f (D 6.7, E 6.0, obs 26) (%.0fmin) ===" %
    (sim.grid.dx, wmax_peak, ng.dx, fvr, (time.time() - _t0) / 60))
