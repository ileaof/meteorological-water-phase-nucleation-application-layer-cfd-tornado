"""Moore Attempt E: the REAL KOUN radiosonde (SRH ~247, CIN -233) instead of ERA5 (SRH 152).

The resolution study (Attempt D) showed the ceiling is the SOURCE circulation, not the mesh, and
ranked *real low-level SRH* as the top lever.  KOUN 2013-05-21 00Z is the actual Norman, OK
radiosonde nearest the Moore tornado: SRH ~247 vs ERA5's 152, shear06 ~32 vs 27 -- much more
streamwise vorticity to tilt into the vertical.  Its strong cap (CIN -233) is no longer a blocker
now that the sustained-ascent forcing breaks caps.  Everything else is Attempt D: storm-relative,
run_multilevel_nest 750 -> 251 -> 84 -> 28 m, measure low-level V_rot vs observed KTLX (26 m/s)."""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "moore_koun_cascade"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress.log"); open(LOG, "w").close()
_t0 = time.time()
def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m)); print(m, flush=True)

import atmospheric_data as ad
from atmospheric_data.sources import iem_raob
from meteorological_flow.grid import Grid
from meteorological_flow.base_state import sounding_diagnostics
from storm_dynamics.config import build_storm_config, MesoForcingConfig
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, nesting as nst, soundings as snd

DEV = os.environ.get("DEV", "gpu")
log("=== Moore CASCADE from the REAL KOUN sounding (SRH ~247), storm-relative, device=%s ===" % DEV)
cache = ad.Cache(os.path.join(REPO, "data", "cache"))
prof = iem_raob.download_sounding("KOUN", "2013-05-21T00:00:00Z", cache=cache)
log("KOUN %s valid %s: %d levels" % (prof.get("station"), prof.get("valid"), prof["height_m"].size))

nx = 64; nz = 44; Lx = 48000.0; Lz = 15000.0
g = Grid(nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, z_stretch=1.05, periodic=True)
base_ground = snd.from_observed_sounding(
    g, pressure_hPa=prof["pressure_Pa"] / 100.0, height_m=prof["height_m"],
    temperature_C=prof["temperature_K"] - 273.15, dewpoint_C=prof["dewpoint_K"] - 273.15,
    u_ms=prof["u_ms"], v_ms=prof["v_ms"])
d = sounding_diagnostics(base_ground)
cx, cy = snd.bunkers_storm_motion(base_ground)
log("REAL KOUN env: CAPE=%.0f CIN=%.0f shear06=%.1f SRH03=%.0f | Bunkers C=(%.1f,%.1f)" %
    (d["CAPE_J_kg"], d["CIN_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(base_ground), cx, cy))

# storm-relative base (Galilean shift; tilting/stretching invariant; storm stays put)
from meteorological_flow.base_state import BaseState
zc = np.asarray(g.zc)
base = BaseState(zc=zc, theta0=base_ground.theta0, qv0=base_ground.qv0, p0=base_ground.p0,
                 T0=base_ground.T0, rho0=base_ground.rho0, u0=base_ground.u0 - cx, v0=base_ground.v0 - cy)

FORCE_S = 1800.0                                   # KOUN cap is stronger (-233) -> force longer/stronger
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, duration=1.0,
                          dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
scfg.dyn.forcing = MesoForcingConfig(enabled=True, heat_rate_K_s=0.009, moist_rate_kgkg_s=4e-6,
                                     radius_m=7000.0, z_top_m=2800.0, duration_s=FORCE_S)
sim = StormSimulation(scfg, base=base)
log("parent %dx%dx%d dx=%.0f m on %s; maturing (forced first %.0fs, breaking the -233 cap) ..." %
    (nx, nx, nz, sim.grid.dx, sim.grid.backend.name, FORCE_S))
T_MAT = 2000.0
step = 0
while sim.t < T_MAT:
    dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t); step += 1
    if step % 80 == 0:
        to = sim.grid.backend.to_cpu; w = np.asarray(to(sim.state.w))
        if not np.isfinite(w).all():
            log("  instability at t=%.0f -> stop" % sim.t); break
        log("  parent t=%4.0f w_max=%5.1f" % (sim.t, w.max()))

to = sim.grid.backend.to_cpu
uc, vc, wc = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wc = np.asarray(to(wc)); zcp = np.asarray(to(sim.grid.zc))
kmeso = int(np.argmin(np.abs(zcp - 700)))
zpar = np.gradient(vc[:, :, kmeso], sim.grid.dx, axis=0) - np.gradient(uc[:, :, kmeso], sim.grid.dy, axis=1)
b = nx // 6
im, jm = np.unravel_index(np.argmax(np.abs(zpar[b:-b, b:-b])), zpar[b:-b, b:-b].shape); im += b; jm += b
mx, my = float(sim.grid.xc[im]), float(sim.grid.yc[jm])
log("PARENT sustained: w_max=%.1f low-lvl meso zeta=%.3e at (%.1f,%.1f) km" %
    (float(np.abs(wc).max()), float(np.abs(zpar[b:-b, b:-b]).max()), mx / 1e3, my / 1e3))

fr = [0.17, 0.17, 0.18]
def mkspec(i):
    def build(gg):
        c = (mx, my) if i == 0 else (0.5 * gg.Lx, 0.5 * gg.Ly)
        return nst.NestSpec.around(gg, c[0], c[1], half=gg.Lx * fr[i], refine=3,
                                   nz=gg.nz, z_stretch=getattr(gg, "z_stretch", 1.05))
    return build
WIN = 180.0
log("cascade 750 -> 250 -> 83 -> ~28 m (refine 3 x3), window %.0f s, storm-relative ..." % WIN)
sims, rep = nst.run_multilevel_nest(sim, [mkspec(0), mkspec(1), mkspec(2)], window=WIN,
                                    les_boost=1.5, cfl=0.2, restrict_momentum=True,
                                    progress=lambda t, w, s: log("  finest t=%5.1f/%.0f step=%d" % (t, w, s))
                                    if s % 80 == 0 else None)
finest = sims[-1]; ng = finest.grid; to = ng.backend.to_cpu
uc, vc, wcn = rot._centered_velocity(finest.state, ng)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wcn)); zcn = np.asarray(to(ng.zc))
kl = int(np.argmin(np.abs(zcn - 500)))
zeta = np.gradient(vc, ng.dx, axis=0) - np.gradient(uc, ng.dy, axis=1)
bb = ng.nx // 6
zlow = float(np.abs(zeta[bb:-bb, bb:-bb, kl]).max())
i2, j2 = np.unravel_index(np.argmax(np.abs(zeta[bb:-bb, bb:-bb, kl])), zeta[bb:-bb, bb:-bb, kl].shape); i2 += bb; j2 += bb
R = max(4, int(1200 / ng.dx))
us = uc[max(0, i2 - R):i2 + R, max(0, j2 - R):j2 + R, kl]; vs = vc[max(0, i2 - R):i2 + R, max(0, j2 - R):j2 + R, kl]
vrot = float(np.sqrt((us - us.mean()) ** 2 + (vs - vs.mean()) ** 2).max())
prof_z = [(round(float(zcn[k])), float(np.abs(zeta[bb:-bb, bb:-bb, k]).max())) for k in range(0, ng.nz, 4)]
np.savez_compressed(os.path.join(OUT, "fields_cascade.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    xc=np.asarray(to(ng.xc)), yc=np.asarray(to(ng.yc)), zc=zcn, dx=float(ng.dx))
summary = {"source": "REAL KOUN 2013-05-21 00Z", "env": {k: float(d[k]) for k in ("CAPE_J_kg", "CIN_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(base_ground)), "bunkers_C": [float(cx), float(cy)],
           "frame": "storm-relative", "levels_dx_m": [float(s.grid.dx) for s in sims], "finest_dx_m": float(ng.dx),
           "finest_nx": int(ng.nx), "window_s": WIN, "forcing": {"heat_K_s": 0.009, "duration_s": FORCE_S},
           "low_level_zeta_interior_s": zlow, "low_level_Vrot_ms": vrot, "zeta_profile_z_val": prof_z,
           "observed": {"Vrot_ms": 26.0, "vorticity_s": 0.205},
           "compare_ERA5_D": {"Vrot_ms": 6.0, "SRH03": 152}, "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
log("=== DONE (KOUN): finest=%.0fm nx=%d LOW-LEVEL zeta=%.3e Vrot=%.1f (ERA5-D was 6.0; obs 26) (%.0fs) ===" %
    (ng.dx, ng.nx, zlow, vrot, time.time() - _t0))
