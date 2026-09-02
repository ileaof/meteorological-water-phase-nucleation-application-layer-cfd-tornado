"""Moore: SUSTAINED (forced) supercell -> deep AMR cascade to ~28 m at the low-level meso.

We finally have a LIVING parent (the sustained-forcing supercell).  Run it storm-relative
(subtract the Bunkers motion from the base wind -- a Galilean shift; tilting/stretching are
frame-invariant) so the storm stays quasi-stationary, then refine 750 -> 250 -> 83 -> ~28 m with
the concurrent multi-level driver (run_multilevel_nest: time-sub-cycled coarse->fine boundaries +
conservative restriction back up) over a long-enough window for vortex stretching to concentrate
the low-level rotation.  Measure the finest low-level V_rot vs the observed KTLX couplet (26 m/s)."""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "moore_cascade"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress.log"); open(LOG, "w").close()
_t0 = time.time()
def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m)); print(m, flush=True)

import atmospheric_data as ad
from atmospheric_data import thermo
from atmospheric_data.sources import era5
from meteorological_flow.grid import Grid
from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config, MesoForcingConfig
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, nesting as nst, soundings as snd

DEV = os.environ.get("DEV", "gpu")
log("=== Moore CASCADE to ~28 m (sustained forced supercell, storm-relative, device=%s) ===" % DEV)
cfg = ad.CaseConfig.from_yaml(os.path.join(REPO, "config", "moore_2013_real.yaml"))
cache = ad.Cache(cfg.data.cache_directory)
st = era5.load(cfg, cache); z_e = np.asarray(st.ds["z"].values)
mean = lambda n: np.asarray(st.ds[n].values)[0].mean(axis=(1, 2))
nx = 64; nz = 44; Lx = 48000.0; Lz = 15000.0
g = Grid(nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, z_stretch=1.05, periodic=True); zc = np.asarray(g.zc)
th0 = np.interp(zc, z_e, mean("theta")); qv0 = np.clip(np.interp(zc, z_e, mean("qv")), 0, None)
u0 = np.interp(zc, z_e, mean("u")); v0 = np.interp(zc, z_e, mean("v"))
psfc = float(np.interp(0, z_e, np.asarray(st.ds["p"].values)[0].mean(axis=(1, 2))))
p0, T0, r0 = thermo.hydrostatic_base_pressure(zc, th0, qv0, psfc)
base_ground = BaseState(zc=zc, theta0=th0, qv0=qv0, p0=p0, T0=T0, rho0=r0, u0=u0, v0=v0)
cx, cy = snd.bunkers_storm_motion(base_ground)
log("Bunkers storm motion C=(%.1f, %.1f) m/s -> shift into storm-relative frame" % (cx, cy))
# storm-relative base: subtract C from the mean wind (shear/CAPE unchanged; storm ~stationary)
base = BaseState(zc=zc, theta0=th0, qv0=qv0, p0=p0, T0=T0, rho0=r0, u0=u0 - cx, v0=v0 - cy)
d = sounding_diagnostics(base_ground)
log("REAL ERA5 env: CAPE=%.0f CIN=%.0f shear06=%.1f SRH03=%.0f" %
    (d["CAPE_J_kg"], d["CIN_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(base_ground)))

FORCE_S = 1500.0
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, duration=1.0,
                          dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20, device=DEV)
scfg.sim.physics.bubble_dtheta = 4.0
scfg.dyn.forcing = MesoForcingConfig(enabled=True, heat_rate_K_s=0.006, moist_rate_kgkg_s=3e-6,
                                     radius_m=7000.0, z_top_m=2500.0, duration_s=FORCE_S)
sim = StormSimulation(scfg, base=base)
log("parent %dx%dx%d dx=%.0f m on %s; maturing (forced first %.0fs) to a sustained supercell ..." %
    (nx, nx, nz, sim.grid.dx, sim.grid.backend.name, FORCE_S))
T_MAT = 1700.0
step = 0
while sim.t < T_MAT:
    dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t); step += 1
    if step % 80 == 0:
        to = sim.grid.backend.to_cpu; w = np.asarray(to(sim.state.w))
        if not np.isfinite(w).all():
            log("  instability at t=%.0f -> stop" % sim.t); break
        log("  parent t=%4.0f w_max=%5.1f" % (sim.t, w.max()))

to = sim.grid.backend.to_cpu
# low-level mesocyclone centre (interior), z ~ 700 m -- where we want to refine
uc, vc, wc = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wc = np.asarray(to(wc)); zcp = np.asarray(to(sim.grid.zc))
kmeso = int(np.argmin(np.abs(zcp - 700)))
zpar = np.gradient(vc[:, :, kmeso], sim.grid.dx, axis=0) - np.gradient(uc[:, :, kmeso], sim.grid.dy, axis=1)
b = nx // 6
im, jm = np.unravel_index(np.argmax(np.abs(zpar[b:-b, b:-b])), zpar[b:-b, b:-b].shape)
im += b; jm += b
mx, my = float(sim.grid.xc[im]), float(sim.grid.yc[jm])
log("PARENT sustained: w_max=%.1f low-lvl meso zeta=%.3e at (%.1f,%.1f) km" %
    (float(np.abs(wc).max()), float(np.abs(zpar[b:-b, b:-b]).max()), mx / 1e3, my / 1e3))

# deep cascade centred on the low-level meso.  L1 on the meso; L2/L3 on their parent centre
# (storm-relative -> stays put).  fracs ~0.17 keep each nest ~65^3 (tractable on the GPU).
fr = [0.17, 0.17, 0.18]
def mkspec(i):
    def build(gg):
        cxc = (mx, my) if i == 0 else (0.5 * gg.Lx, 0.5 * gg.Ly)
        return nst.NestSpec.around(gg, cxc[0], cxc[1], half=gg.Lx * fr[i], refine=3,
                                   nz=gg.nz, z_stretch=getattr(gg, "z_stretch", 1.05))
    return build
WIN = 180.0
log("cascade 750 -> 250 -> 83 -> ~28 m (refine 3 x3), window %.0f s, storm-relative ..." % WIN)
sims, rep = nst.run_multilevel_nest(sim, [mkspec(0), mkspec(1), mkspec(2)], window=WIN,
                                    les_boost=1.5, cfl=0.2, restrict_momentum=True,
                                    progress=lambda t, w, s: log("  finest t=%5.1f/%.0f step=%d" % (t, w, s))
                                    if s % 60 == 0 else None)
finest = sims[-1]; ng = finest.grid; to = ng.backend.to_cpu
uc, vc, wcn = rot._centered_velocity(finest.state, ng)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wcn)); zcn = np.asarray(to(ng.zc))
kl = int(np.argmin(np.abs(zcn - 500)))
zeta = np.gradient(vc, ng.dx, axis=0) - np.gradient(uc, ng.dy, axis=1)
bb = ng.nx // 6
zlow = float(np.abs(zeta[bb:-bb, bb:-bb, kl]).max())
i2, j2 = np.unravel_index(np.argmax(np.abs(zeta[bb:-bb, bb:-bb, kl])), zeta[bb:-bb, bb:-bb, kl].shape)
i2 += bb; j2 += bb; R = max(4, int(1200 / ng.dx))
us = uc[max(0, i2 - R):i2 + R, max(0, j2 - R):j2 + R, kl]; vs = vc[max(0, i2 - R):i2 + R, max(0, j2 - R):j2 + R, kl]
vrot = float(np.sqrt((us - us.mean()) ** 2 + (vs - vs.mean()) ** 2).max())
# vertical profile of peak |zeta| (is the rotation reaching the surface?)
prof = [(round(float(zcn[k])), float(np.abs(zeta[bb:-bb, bb:-bb, k]).max())) for k in range(0, ng.nz, 4)]
np.savez_compressed(os.path.join(OUT, "fields_cascade.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    xc=np.asarray(to(ng.xc)), yc=np.asarray(to(ng.yc)), zc=zcn, dx=float(ng.dx))
summary = {"env": {k: float(d[k]) for k in ("CAPE_J_kg", "CIN_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(base_ground)),
           "bunkers_C": [float(cx), float(cy)], "frame": "storm-relative",
           "levels_dx_m": [float(s.grid.dx) for s in sims], "finest_dx_m": float(ng.dx),
           "finest_nx": int(ng.nx), "window_s": WIN,
           "low_level_zeta_interior_s": zlow, "low_level_Vrot_ms": vrot,
           "zeta_profile_z_val": prof,
           "observed": {"Vrot_ms": 26.0, "vorticity_s": 0.205},
           "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
log("=== DONE: levels dx=%s finest=%.0fm nx=%d LOW-LEVEL zeta=%.3e Vrot=%.1f (obs 0.205/26) (%.0fs) ===" %
    ([round(x) for x in summary["levels_dx_m"]], ng.dx, ng.nx, zlow, vrot, time.time() - _t0))
