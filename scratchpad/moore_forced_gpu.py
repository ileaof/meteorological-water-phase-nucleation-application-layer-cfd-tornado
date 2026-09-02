"""Moore REAL-ERA5 env + SUSTAINED mesoscale-ascent forcing (dryline proxy) on the GPU.
Finer parent (dx 750 m) so the updraft is resolved; a low-level heating(+moistening)
cylinder held 25 min lifts parcels through the real CIN cap so the supercell SUSTAINS
(the single bubble decayed w9->1).  Then a storm-following nest measures the low-level
rotation vs the observed KTLX couplet (Vrot 26, zeta 0.21)."""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "moore_forced"); os.makedirs(OUT, exist_ok=True)
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
log("=== Moore FORCED supercell (REAL ERA5 env + sustained ascent, device=%s) ===" % DEV)
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
base = BaseState(zc=zc, theta0=th0, qv0=qv0, p0=p0, T0=T0, rho0=r0, u0=u0, v0=v0)
d = sounding_diagnostics(base)
log("REAL ERA5 env: CAPE=%.0f CIN=%.0f shear06=%.1f SRH03=%.0f" %
    (d["CAPE_J_kg"], d["CIN_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(base)))

FORCE_S = 1500.0
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=Lz, duration=1.0,
                          dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20, device=DEV)
scfg.sim.physics.bubble_dtheta = 4.0
scfg.dyn.forcing = MesoForcingConfig(enabled=True, heat_rate_K_s=0.006, moist_rate_kgkg_s=3e-6,
                                     radius_m=7000.0, z_top_m=2500.0, duration_s=FORCE_S)
sim = StormSimulation(scfg, base=base)
log("parent %dx%dx%d dx=%.0f m on %s; sustained forcing first %.0f s then free ..." %
    (nx, nx, nz, sim.grid.dx, sim.grid.backend.name, FORCE_S))
zlow = int(np.argmin(np.abs(zc - 500)))
hist = []
T_END = 1900.0
step = 0
while sim.t < T_END:
    dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t); step += 1
    if step % 40 == 0:
        to = sim.grid.backend.to_cpu; w = np.asarray(to(sim.state.w))
        if not np.isfinite(w).all():
            log("  instability at t=%.0f -> stop" % sim.t); break
        uc, vc, _ = rot._centered_velocity(sim.state, sim.grid)
        uc = np.asarray(to(uc)); vc = np.asarray(to(vc))
        b = nx // 6
        zl = np.gradient(vc[:, :, zlow], sim.grid.dx, axis=0) - np.gradient(uc[:, :, zlow], sim.grid.dy, axis=1)
        zli = float(np.abs(zl[b:-b, b:-b]).max())
        flag = "FORCED" if sim.t < FORCE_S else "free  "
        log("  t=%4.0f [%s] w_max=%5.1f low-lvl zeta(int)=%.2e" % (sim.t, flag, w.max(), zli))
        hist.append((sim.t, float(w.max()), zli))

to = sim.grid.backend.to_cpu
rr = rot.rotation_report(sim.state, sim.grid, base=sim.base)
log("PARENT after run: w_max=%.1f midmeso=%.2e zeta=%.2e" %
    (float(np.abs(np.asarray(to(sim.state.w))).max()), rr.get("midlevel_mesocyclone", 0), rr["zeta_abs_max"]))

# storm-following nest on the updraft
_, _, wc = rot._centered_velocity(sim.state, sim.grid); wc = np.asarray(to(wc))
i, j = np.unravel_index(np.argmax(wc.max(axis=2)), wc.shape[:2])
half = nx // 4; i0 = int(np.clip(i - half // 2, 1, nx - half - 1)); j0 = int(np.clip(j - half // 2, 1, nx - half - 1))
spec = nst.NestSpec.aligned(sim.grid, i0=i0, j0=j0, ncx=half, ncy=half, refine=3)
log("storm-following nest dx=%.0f m over the updraft, window 600 s ..." % (sim.grid.dx / 3))
nest, rep = nst.run_concurrent_nest(sim, spec, window=600.0, follow=True, les_boost=1.5, cfl=0.2,
                                    progress=lambda t, w, s: log("  nest t=%4.0f/%.0f" % (t, w)) if int(t) % 150 < 5 else None)
ng = nest.grid; to = ng.backend.to_cpu; uc, vc, wcn = rot._centered_velocity(nest.state, ng)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wcn)); zcn = np.asarray(to(ng.zc))
kl = int(np.argmin(np.abs(zcn - 500)))
zeta = np.gradient(vc, ng.dx, axis=0) - np.gradient(uc, ng.dy, axis=1)
b = ng.nx // 6
zlow_int = float(np.abs(zeta[b:-b, b:-b, kl]).max())
i2, j2 = np.unravel_index(np.argmax(np.abs(zeta[b:-b, b:-b, kl])), zeta[b:-b, b:-b, kl].shape)
i2 += b; j2 += b; R = max(3, int(1500 / ng.dx))
us = uc[max(0, i2 - R):i2 + R, max(0, j2 - R):j2 + R, kl]; vs = vc[max(0, i2 - R):i2 + R, max(0, j2 - R):j2 + R, kl]
vrot = float(np.sqrt((us - us.mean()) ** 2 + (vs - vs.mean()) ** 2).max())
np.savez_compressed(os.path.join(OUT, "fields_nest.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    xc=np.asarray(to(ng.xc)), yc=np.asarray(to(ng.yc)), zc=zcn, dx=float(ng.dx))
summary = {"env": {k: float(d[k]) for k in ("CAPE_J_kg", "CIN_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(base)), "parent_dx_m": float(sim.grid.dx),
           "forcing": {"heat_K_s": 0.006, "moist_kgkg_s": 3e-6, "duration_s": FORCE_S},
           "nest_dx_m": float(ng.dx), "low_level_zeta_interior_s": zlow_int, "low_level_Vrot_ms": vrot,
           "observed": {"Vrot_ms": 26.0, "vorticity_s": 0.205},
           "wall_clock_s": time.time() - _t0, "history": hist}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
log("=== DONE: nest dx=%.0fm LOW-LEVEL zeta=%.3e Vrot=%.1f (obs 0.205/26) (%.0fs) ===" %
    (ng.dx, zlow_int, vrot, time.time() - _t0))
