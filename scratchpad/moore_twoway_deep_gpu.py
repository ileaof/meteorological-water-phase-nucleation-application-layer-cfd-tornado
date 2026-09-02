"""Attempt H -- the decisive run: real KOUN env + sustained forcing + storm-relative deep cascade
to ~28 m WITH two-way coupling (now wired into run_multilevel_nest), then READ the low-level
vorticity budget on the result.  Combines the validated lever (two-way), the instrument (budget),
and real data.  Answers: does two-way in the deep cascade push low-level V_rot past 11.6 toward the
observed 26, and WHICH budget term (baroclinic vs tilting vs stretching) still limits it?"""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "moore_twoway_deep"); os.makedirs(OUT, exist_ok=True)
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
from storm_dynamics import vorticity_budget as vb, vortex_diagnostics as vd
from storm_dynamics import coldpool as cp, classification as cl

def fair_low_vrot(u, v, zeta, zc, dx, zmax=1500.0):
    nb = zeta.shape[0] // 6; best = (0.0, 0.0, 0.0)
    for k in range(zeta.shape[2]):
        if zc[k] > zmax:
            break
        Z = np.abs(zeta[nb:-nb, nb:-nb, k]); i, j = np.unravel_index(np.argmax(Z), Z.shape); i += nb; j += nb
        R = max(4, int(1000 / dx))
        us = u[max(0, i - R):i + R, max(0, j - R):j + R, k]; vs = v[max(0, i - R):i + R, max(0, j - R):j + R, k]
        vr = float(np.sqrt((us - us.mean()) ** 2 + (vs - vs.mean()) ** 2).max())
        if vr > best[0]:
            best = (vr, float(zc[k]), float(Z.max()))
    return best

DEV = os.environ.get("DEV", "gpu")
log("=== Attempt H: KOUN + TWO-WAY deep cascade + vorticity budget (%s) ===" % DEV)
cache = ad.Cache(os.path.join(REPO, "data", "cache"))
prof = iem_raob.download_sounding("KOUN", "2013-05-21T00:00:00Z", cache=cache)
nx = 64; nz = 44; Lx = 48000.0
g = Grid(nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=15000.0, z_stretch=1.05, periodic=True); zc = np.asarray(g.zc)
bg = snd.from_observed_sounding(g, pressure_hPa=prof["pressure_Pa"] / 100.0, height_m=prof["height_m"],
                                temperature_C=prof["temperature_K"] - 273.15, dewpoint_C=prof["dewpoint_K"] - 273.15,
                                u_ms=prof["u_ms"], v_ms=prof["v_ms"])
d = sounding_diagnostics(bg); cx, cy = snd.bunkers_storm_motion(bg)
log("REAL KOUN: CAPE=%.0f CIN=%.0f shear06=%.1f SRH03=%.0f" %
    (d["CAPE_J_kg"], d["CIN_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(bg)))
base = BaseState(zc=zc, theta0=bg.theta0, qv0=bg.qv0, p0=bg.p0, T0=bg.T0, rho0=bg.rho0, u0=bg.u0 - cx, v0=bg.v0 - cy)

scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=15000.0, duration=1.0,
                          dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
scfg.dyn.forcing = MesoForcingConfig(enabled=True, heat_rate_K_s=0.009, moist_rate_kgkg_s=4e-6,
                                     radius_m=7000.0, z_top_m=2800.0, duration_s=1800.0)
sim = StormSimulation(scfg, base=base)
log("maturing (forced, breaking the -247 cap) ...")
while sim.t < 2000.0:
    dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t)
    if sim.step % 200 == 0:
        w = np.asarray(sim.grid.backend.to_cpu(sim.state.w)); log("  parent t=%4.0f w_max=%.1f" % (sim.t, w.max()))

to = sim.grid.backend.to_cpu; uc, vc, wc = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); zcp = np.asarray(to(sim.grid.zc))
kk = int(np.argmin(np.abs(zcp - 500))); nbp = nx // 6
zl = np.abs(np.gradient(vc[:, :, kk], sim.grid.dx, axis=0) - np.gradient(uc[:, :, kk], sim.grid.dy, axis=1))
im, jm = np.unravel_index(np.argmax(zl[nbp:-nbp, nbp:-nbp]), zl[nbp:-nbp, nbp:-nbp].shape); im += nbp; jm += nbp
mx, my = float(sim.grid.xc[im]), float(sim.grid.yc[jm])
log("PARENT matured: w_max=%.1f meso@(%.1f,%.1f)km" %
    (float(np.abs(np.asarray(to(sim.state.w))).max()), mx / 1e3, my / 1e3))

fr = [0.17, 0.17, 0.18]
def mkspec(i):
    def build(gg):
        ncx = max(6, int(2 * fr[i] * gg.nx))              # parent cells covered (footprint ~ around's)
        if i == 0:
            ic = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.xc)) - mx)))
            jc = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.yc)) - my)))
        else:
            ic, jc = gg.nx // 2, gg.ny // 2
        i0 = int(np.clip(ic - ncx // 2, 1, gg.nx - ncx - 1)); j0 = int(np.clip(jc - ncx // 2, 1, gg.ny - ncx - 1))
        return nst.NestSpec.aligned(gg, i0=i0, j0=j0, ncx=ncx, ncy=ncx, refine=3)
    return build
WIN = 210.0
log("TWO-WAY deep cascade 750->251->84->~28 m (window %.0f s, two_way=True) ..." % WIN)
sims, rep = nst.run_multilevel_nest(sim, [mkspec(0), mkspec(1), mkspec(2)], window=WIN,
                                    two_way=True, two_way_rate=0.5, storm_motion=(0.0, 0.0),
                                    les_boost=1.5, cfl=0.2,
                                    progress=lambda t, w, s: log("  finest t=%5.1f/%.0f step=%d" % (t, w, s))
                                    if s % 80 == 0 else None)
finest = sims[-1]; ng = finest.grid; to = ng.backend.to_cpu
finest.state.diagnose(finest.cfg)                         # rho, P_total for the baroclinic term
uc, vc, wcn = rot._centered_velocity(finest.state, ng)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wcn)); zcn = np.asarray(to(ng.zc))
zeta = np.gradient(vc, ng.dx, axis=0) - np.gradient(uc, ng.dy, axis=1)
fvr = fair_low_vrot(uc, vc, zeta, zcn, ng.dx)

# SAVE FIELDS FIRST (before any fragile diagnostic) so the 2 h run can never be lost.
np.savez_compressed(os.path.join(OUT, "fields_deep.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    rho=np.asarray(to(finest.state.rho)), p=np.asarray(to(finest.state.P_total)),
                    xc=np.asarray(to(ng.xc)), yc=np.asarray(to(ng.yc)), zc=zcn, dx=float(ng.dx))
log("fields saved; low-lvl Vrot(fair)=%.2f at z=%.0f" % (fvr[0], fvr[1]))

# --- the instrument: low-level vorticity budget on the two-way deep-cascade result ---
low, dom, prod, vrep, crep, category = {}, "n/a", {}, {}, {}, "n/a"
try:
    terms = vb.zeta_budget(finest.state, ng, Km=getattr(finest, "_Km", None))
    low = vb.budget_layer_summary(terms, ng, 0.0, 1000.0)
    dom, prod = vb.dominant_mechanism(terms, ng, 0.0, 1000.0)
    log("BUDGET low 0-1km: dominant=%s  baroclinic=%.2e tilting=%.2e stretching=%.2e (|mean|)" %
        (dom, low["baroclinic_absmean"], low["tilting_absmean"], low["stretching_absmean"]))
    vrep = vd.vortex_report(finest.state, ng, z_m=150.0)
    crep = cp.coldpool_report(finest.state, ng, z_m=150.0)
    category = cl.classify_simulation(finest, z_surface_m=150.0)["category"]
except Exception as e:
    log("diagnostics error (fields are saved, re-run offline): %r" % e)
summary = {"attempt": "H KOUN two-way deep cascade", "two_way": True,
           "env": {k: float(d[k]) for k in ("CAPE_J_kg", "CIN_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(bg)), "finest_dx_m": float(ng.dx), "window_s": WIN,
           "low_level_vrot_fair": fvr[0], "low_level_vrot_z": fvr[1],
           "budget_low_0_1km": low, "dominant_mechanism": dom, "budget_production": prod,
           "vortex": vrep, "cold_pool": crep, "classification": category,
           "reference": {"G_oneway_83m": 8.1, "G_twoway_83m": 11.6, "D_oneway_28m": 6.7, "obs": 26.0},
           "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
log("=== DONE H: finest=%.0fm low-lvl Vrot(fair)=%.1f | class=%s | dominant low-lvl term=%s "
    "(G two-way 83m was 11.6, obs 26) (%.0fmin) ===" %
    (ng.dx, fvr[0], category, dom, (time.time() - _t0) / 60))
