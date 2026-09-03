"""Attempt J -- the culmination: resolve the tornado-scale vortex in the FREELY-EVOLVED supercell.

Attempt I showed a freely-evolving idealised supercell (no held forcing) builds the correct low-level
geometry (alignment 0.02->0.20, streamwise 0.40->0.64, low-level |zeta| 3.4x) as it occludes -- the
one thing the forced runs A-H lacked.  Here we mature that storm to its low-level-meso peak (~t=2800
s), then drop a deep storm-following AMR cascade (600->200->67->~22 m) centred on the low-level
mesocyclone to RESOLVE the vortex -- with the vortex-line geometry already aligned.  Then read
V_rot, circulation, pressure deficit, the tilting alignment, and classification on the finest level.
Hardened: fields saved BEFORE diagnostics; diagnostics guarded."""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "tornado_occlusion"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress.log"); open(LOG, "w").close()
_t0 = time.time()
def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m)); print(m, flush=True)

from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, soundings as snd, nesting as nst
from storm_dynamics import vorticity_budget as vb, vortex_diagnostics as vd, coldpool as cp, classification as cl

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

def low_align(sim, z_top=800.0):
    g = sim.grid; to = g.backend.to_cpu
    uc, vc, wc = rot._centered_velocity(sim.state, g)
    e = vb.tilting_efficiency(uc, vc, wc, g); zc = np.asarray(to(g.zc)); nb = g.nx // 6; km = zc < z_top
    c = lambda a: np.asarray(to(a))[nb:-nb, nb:-nb][:, :, km]
    return float(c(e["tilting"]).sum() / (np.abs(c(e["omega_h"]) * c(e["grad_h_w"])).sum() + 1e-20))

DEV = os.environ.get("DEV", "gpu")
log("=== Attempt J: tornado-scale vortex in the freely-evolved supercell (%s) ===" % DEV)
nx = 120; nz = 48; Lx = 72000.0
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=15000.0,
                          duration=1.0, dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20,
                          hodograph_kind="quarter_circle", U_max=30.0, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
sim0 = StormSimulation(scfg); b = sim0.base
d = sounding_diagnostics(b); cx, cy = snd.bunkers_storm_motion(b)
log("env: CAPE=%.0f shear06=%.1f SRH03=%.0f | Bunkers=(%.1f,%.1f) dx=%.0fm" %
    (d["CAPE_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(b), cx, cy, sim0.grid.dx))
base_sr = BaseState(zc=b.zc, theta0=b.theta0, qv0=b.qv0, p0=b.p0, T0=b.T0, rho0=b.rho0, u0=b.u0 - cx, v0=b.v0 - cy)
sim = StormSimulation(scfg, base=base_sr); sim.state.diagnose(sim.cfg)

T_MAT = 2800.0                                    # the low-level-meso peak (Attempt I)
log("maturing the freely-evolving supercell to its low-level-meso peak (t=%.0f s) ..." % T_MAT)
while sim.t < T_MAT:
    dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t)
    if sim.step % 300 == 0:
        w = np.asarray(sim.grid.backend.to_cpu(sim.state.w))
        if not np.isfinite(w).all():
            log("  instability at t=%.0f -> stop" % sim.t); break
        sim.state.diagnose(sim.cfg)
        log("  parent t=%4.0f w_max=%5.1f align=%+.3f" % (sim.t, w.max(), low_align(sim)))

# locate the low-level mesocyclone (z~500 m interior)
to = sim.grid.backend.to_cpu; uc, vc, _ = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); zcp = np.asarray(to(sim.grid.zc))
kk = int(np.argmin(np.abs(zcp - 500))); nbp = nx // 6
zl = np.abs(np.gradient(vc[:, :, kk], sim.grid.dx, axis=0) - np.gradient(uc[:, :, kk], sim.grid.dy, axis=1))
im, jm = np.unravel_index(np.argmax(zl[nbp:-nbp, nbp:-nbp]), zl[nbp:-nbp, nbp:-nbp].shape); im += nbp; jm += nbp
mx, my = float(sim.grid.xc[im]), float(sim.grid.yc[jm])
log("PARENT matured: w_max=%.1f align=%+.3f meso@(%.1f,%.1f)km" %
    (float(np.abs(np.asarray(to(sim.state.w))).max()), low_align(sim), mx / 1e3, my / 1e3))

# deep one-way cascade 600 -> 200 -> 67 -> ~22 m centred on the low-level meso (geometry already
# right).  ncx=28 -> ~84^3 nests with LARGER footprints (16.8/5.6/1.9 km) so the vortex stays
# interior despite residual storm drift (the ncx=18 run drifted the vortex to the finest edge).
NCX = 28
def mkspec(i):
    def build(gg):
        ncx = NCX
        if i == 0:
            ic = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.xc)) - mx)))
            jc = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.yc)) - my)))
        else:
            ic, jc = gg.nx // 2, gg.ny // 2
        i0 = int(np.clip(ic - ncx // 2, 1, gg.nx - ncx - 1)); j0 = int(np.clip(jc - ncx // 2, 1, gg.ny - ncx - 1))
        return nst.NestSpec.aligned(gg, i0=i0, j0=j0, ncx=ncx, ncy=ncx, refine=3)
    return build
WIN = 200.0
log("deep cascade 600 -> 200 -> 67 -> ~22 m (window %.0f s, one-way) ..." % WIN)
sims, rep = nst.run_multilevel_nest(sim, [mkspec(0), mkspec(1), mkspec(2)], window=WIN,
                                    restrict_up=True, restrict_momentum=True, les_boost=1.5, cfl=0.2,
                                    progress=lambda t, w, s: log("  finest t=%5.1f/%.0f step=%d" % (t, w, s))
                                    if s % 100 == 0 else None)
finest = sims[-1]; ng = finest.grid; to = ng.backend.to_cpu
finest.state.diagnose(finest.cfg)
uc, vc, wcn = rot._centered_velocity(finest.state, ng)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wcn)); zcn = np.asarray(to(ng.zc))
zeta = np.gradient(vc, ng.dx, axis=0) - np.gradient(uc, ng.dy, axis=1)
fvr = fair_low_vrot(uc, vc, zeta, zcn, ng.dx)
# vertical profile of the INTERIOR V_rot (surface-connected vs elevated?), + edge check
nbf = ng.nx // 5
prof = []
for zt in (50., 200., 500., 1000., 1500.):
    k = int(np.argmin(np.abs(zcn - zt)))
    Zi = np.abs(zeta[nbf:-nbf, nbf:-nbf, k])
    ii, jj = np.unravel_index(int(np.argmax(Zi)), Zi.shape); ii += nbf; jj += nbf
    R = max(4, int(1000 / ng.dx))
    us = uc[max(0, ii - R):ii + R, max(0, jj - R):jj + R, k]; vs = vc[max(0, ii - R):ii + R, max(0, jj - R):jj + R, k]
    vr = float(np.sqrt((us - us.mean()) ** 2 + (vs - vs.mean()) ** 2).max())
    edge = int(min(ii, jj, ng.nx - 1 - ii, ng.ny - 1 - jj))
    prof.append({"z": round(float(zcn[k])), "vrot": round(vr, 2), "zeta": float(Zi.max()), "edge_cells": edge})
    log("  Vrot profile z=%5.0fm  Vrot=%.2f  |zeta|=%.2e  [%d cells from edge]" % (zcn[k], vr, Zi.max(), edge))
np.savez_compressed(os.path.join(OUT, "fields_tornado.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    rho=np.asarray(to(finest.state.rho)), p=np.asarray(to(finest.state.P_total)),
                    xc=np.asarray(to(ng.xc)), yc=np.asarray(to(ng.yc)), zc=zcn, dx=float(ng.dx))
log("fields saved; finest dx=%.0fm  low-lvl Vrot(fair)=%.2f at z=%.0f" % (ng.dx, fvr[0], fvr[1]))

vrep, crep, category, align_fine = {}, {}, "n/a", 0.0
try:
    vrep = vd.vortex_report(finest.state, ng, z_m=100.0)
    crep = cp.coldpool_report(finest.state, ng, z_m=100.0)
    category = cl.classify_simulation(finest, z_surface_m=100.0)["category"]
    align_fine = low_align(finest)
    log("VORTEX: Vtheta_max=%.1f circ=%.2e Delta_p=%.1f core=%.0fm | class=%s align=%+.3f" %
        (vrep["v_theta_max_m_s"], vrep["circulation_m2_s"], vrep["pressure_deficit_Pa"],
         vrep["core_radius_m"], category, align_fine))
except Exception as e:
    log("diagnostics error (fields saved): %r" % e)

summary = {"attempt": "J tornado in freely-evolved supercell", "env": {k: float(d[k]) for k in ("CAPE_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(b)), "finest_dx_m": float(ng.dx), "window_s": WIN,
           "low_level_vrot_fair": fvr[0], "low_level_vrot_z": fvr[1], "fine_align": align_fine,
           "vrot_profile": prof,
           "vortex": vrep, "cold_pool": crep, "classification": category,
           "reference": {"obs_Vrot": 26.0, "D_oneway_28m": 6.7, "G_twoway_83m": 11.6},
           "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
log("=== DONE J: finest=%.0fm low-lvl Vrot=%.1f class=%s align=%+.3f (obs 26, D 6.7, G 11.6) (%.0fmin) ===" %
    (ng.dx, fvr[0], category, align_fine, (time.time() - _t0) / 60))
