"""Moore 2013 REAL-ERA5-environment funnel attempt on the GPU: 3-level AMR cascade
1.25 km -> 417 -> 139 -> 46 m, low-memory pressure solver.  Logs progress + saves L3 fields."""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "moore_real_funnel"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress.log"); open(LOG, "w").close()
_t0 = time.time()
def log(m):
    line = "[%6.1fs] %s" % (time.time() - _t0, m)
    open(LOG, "a", encoding="utf-8").write(line + "\n"); print(line, flush=True)

import atmospheric_data as ad
from atmospheric_data import thermo
from atmospheric_data.sources import era5
from meteorological_flow.grid import Grid
from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, nesting as nst

DEV = os.environ.get("FUNNEL_DEVICE", "gpu")
log("=== Moore REAL-ERA5 funnel (device=%s) ===" % DEV)

# 1) real ERA5 environmental profile -> base state on a fine grid
cfg = ad.CaseConfig.from_yaml(os.path.join(REPO, "config", "moore_2013_real.yaml"))
cache = ad.Cache(cfg.data.cache_directory)
st = era5.load(cfg, cache); z_e = np.asarray(st.ds["z"].values)
mean = lambda n: np.asarray(st.ds[n].values)[0].mean(axis=(1, 2))
nx = 48; nz = 40; Lx = 60000.0; Lz = 15000.0
g = Grid(nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, z_stretch=1.05, periodic=True)
zc = np.asarray(g.zc)
th0 = np.interp(zc, z_e, mean("theta")); qv0 = np.clip(np.interp(zc, z_e, mean("qv")), 0, None)
u0 = np.interp(zc, z_e, mean("u")); v0 = np.interp(zc, z_e, mean("v"))
psfc = float(np.interp(0, z_e, np.asarray(st.ds["p"].values)[0].mean(axis=(1, 2))))
p0, T0, r0 = thermo.hydrostatic_base_pressure(zc, th0, qv0, psfc)
base = BaseState(zc=zc, theta0=th0, qv0=qv0, p0=p0, T0=T0, rho0=r0, u0=u0, v0=v0)
d = sounding_diagnostics(base)
log("REAL ERA5 base: CAPE=%.0f CIN=%.0f LCL=%.0f shear06=%.1f" %
    (d["CAPE_J_kg"], d["CIN_J_kg"], d["LCL_m"] or -1, d["shear_0_6km_m_s"]))

# 2) supercell in the real environment (strong bubble breaks the cap), mature to ~peak
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, duration=1.0,
                          dt_max=4.0, drag=True, z_stretch=1.05, C_s=0.22, device=DEV)
scfg.sim.physics.bubble_dtheta = 6.0
sim = StormSimulation(scfg, base=base)
log("parent %dx%dx%d on %s, maturing to peak ..." % (nx, nx, nz, sim.grid.backend.name))
sim.cfg.time.duration = 260.0
sim.run(progress=lambda t, dd, s: log("  parent t=%3.0f/%.0f" % (t, dd)) if int(t) % 60 < 4 else None)
rr = rot.rotation_report(sim.state, sim.grid, base=sim.base)
log("PARENT@peak: w_max=%.1f zeta=%.2e midmeso=%.2e" %
    (float(np.abs(np.asarray(sim.grid.backend.to_cpu(sim.state.w))).max()),
     rr["zeta_abs_max"], rr.get("midlevel_mesocyclone", 0)))

# 3) 3-level AMR cascade toward ~46 m (centred on the domain = the storm).  Small half-fractions
#    keep each nest ~60 cells wide (tractable): L1~58^3 (417 m), L2~62^3 (139 m), L3~67^3 (46 m).
half_fracs = [0.20, 0.18, 0.18]
mkspec = lambda i: (lambda gg: nst.NestSpec.around(gg, 0.5 * gg.Lx, 0.5 * gg.Ly,
                    half=gg.Lx * half_fracs[i], refine=3, nz=gg.nz,
                    z_stretch=getattr(gg, "z_stretch", 1.05)))
log("cascade 1.25km -> /3 -> /3 -> /3 (~46 m), window 90 s ...")
sims, rep = nst.run_multilevel_nest(sim, [mkspec(0), mkspec(1), mkspec(2)], window=90.0,
                                    les_boost=1.5, cfl=0.2,
                                    progress=lambda t, w, s: log("  finest t=%4.1f/%.0f step=%d" % (t, w, s))
                                    if s % 40 == 0 else None)
finest = sims[-1]
to = finest.grid.backend.to_cpu
uc, vc, wc = rot._centered_velocity(finest.state, finest.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wc))
zeta = np.gradient(vc, finest.grid.dx, axis=0) - np.gradient(uc, finest.grid.dy, axis=1)
np.savez_compressed(os.path.join(OUT, "fields_L3.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    xc=np.asarray(to(finest.grid.xc)), yc=np.asarray(to(finest.grid.yc)),
                    zc=np.asarray(to(finest.grid.zc)), dx=float(finest.grid.dx))
summary = {"env": {k: float(d[k]) for k in ("CAPE_J_kg", "CIN_J_kg", "LCL_m", "shear_0_6km_m_s") if d[k] is not None},
           "levels_dx_m": [float(s.grid.dx) for s in sims],
           "finest_dx_m": float(finest.grid.dx),
           "finest_zeta_abs_max": float(np.abs(zeta).max()),
           "finest_w_max": float(np.abs(wcn).max()),
           "rotation_report_zeta": float(rep["rotation"]["zeta_abs_max"]),
           "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
log("=== DONE: levels dx=%s finest=%.0fm zeta=%.3e w_max=%.1f (%.0fs) ===" %
    ([round(x) for x in summary["levels_dx_m"]], summary["finest_dx_m"],
     summary["finest_zeta_abs_max"], summary["finest_w_max"], summary["wall_clock_s"]))
