"""Attempt K -- the intensity run: everything the A-J arc established, combined.

  * freely-evolving supercell (I) -> correct streamwise/aligned low-level geometry
  * surface-layer STRESS-DIVERGENCE drag + height-consistent LOG-LAW C_d (the closure fix that
    restored surface connection instead of stripping the tangential wind)
  * a nest that resolves the SURFACE LAYER vertically (dz1 ~ 5 m) as well as horizontally
  * a MOVING nest (follow_interval) so the vortex cannot drift out of the fine domain (Attempt J)

The vertical CFL forbids a fine near-surface mesh on the parent for a 45-min maturation
(dz1 ~ 4 m => dt ~ 0.05 s), so the surface layer is refined only inside the nest, which runs a short
window.  Question: with the geometry right AND the surface closure right AND the corner-flow layer
resolved, how intense does the surface-connected vortex get?  Scored by surface_connection_report
(sfc/aloft ratio), vortex_report (V_theta, circulation, Delta p) and the objective classifier."""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "tornado_intensity"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress.log"); open(LOG, "w").close()
_t0 = time.time()
def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m)); print(m, flush=True)

from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, soundings as snd, nesting as nst
from storm_dynamics import vorticity_budget as vb, vortex_diagnostics as vd
from storm_dynamics import coldpool as cp, classification as cl

def low_align(sim, z_top=800.0):
    g = sim.grid; to = g.backend.to_cpu
    uc, vc, wc = rot._centered_velocity(sim.state, g)
    e = vb.tilting_efficiency(uc, vc, wc, g); zc = np.asarray(to(g.zc)); nb = g.nx // 6; km = zc < z_top
    c = lambda a: np.asarray(to(a))[nb:-nb, nb:-nb][:, :, km]
    return float(c(e["tilting"]).sum() / (np.abs(c(e["omega_h"]) * c(e["grad_h_w"])).sum() + 1e-20))

DEV = os.environ.get("DEV", "gpu")
log("=== Attempt K: intensity -- free evolution + fixed surface closure + resolved surface layer ===")
nx = 120; nz = 48; Lx = 72000.0
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=15000.0,
                          duration=1.0, dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20,
                          hodograph_kind="quarter_circle", U_max=30.0, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
# THE SURFACE CLOSURE FIX (applies on the parent and is inherited by the nests)
scfg.dyn.drag.stress_divergence = True
scfg.dyn.drag.surface_layer_depth_m = 150.0
scfg.dyn.drag.use_log_law = True
scfg.dyn.drag.roughness_length_m = 0.1

sim0 = StormSimulation(scfg); b = sim0.base
d = sounding_diagnostics(b); cx, cy = snd.bunkers_storm_motion(b)
log("env: CAPE=%.0f shear06=%.1f SRH03=%.0f | parent dx=%.0fm dz1=%.1fm | drag=stress-div+log-law" %
    (d["CAPE_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(b), sim0.grid.dx,
     float(np.asarray(sim0.grid.backend.to_cpu(sim0.grid.zc))[0])))
base_sr = BaseState(zc=b.zc, theta0=b.theta0, qv0=b.qv0, p0=b.p0, T0=b.T0, rho0=b.rho0,
                    u0=b.u0 - cx, v0=b.v0 - cy)
sim = StormSimulation(scfg, base=base_sr); sim.state.diagnose(sim.cfg)

T_MAT = 2800.0
log("maturing the freely-evolving supercell to its low-level-meso peak ...")
while sim.t < T_MAT:
    dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t)
    if sim.step % 400 == 0:
        w = np.asarray(sim.grid.backend.to_cpu(sim.state.w))
        if not np.isfinite(w).all():
            log("  instability at t=%.0f -> stop" % sim.t); break
        sim.state.diagnose(sim.cfg)
        log("  parent t=%4.0f w_max=%5.1f align=%+.3f" % (sim.t, w.max(), low_align(sim)))

to = sim.grid.backend.to_cpu; uc, vc, _ = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); zcp = np.asarray(to(sim.grid.zc))
kk = int(np.argmin(np.abs(zcp - 500))); nbp = nx // 6
zl = np.abs(np.gradient(vc[:, :, kk], sim.grid.dx, axis=0) - np.gradient(uc[:, :, kk], sim.grid.dy, axis=1))
im, jm = np.unravel_index(np.argmax(zl[nbp:-nbp, nbp:-nbp]), zl[nbp:-nbp, nbp:-nbp].shape); im += nbp; jm += nbp
mx, my = float(sim.grid.xc[im]), float(sim.grid.yc[jm])
log("PARENT matured: w_max=%.1f align=%+.3f meso@(%.1f,%.1f)km" %
    (float(np.abs(np.asarray(to(sim.state.w))).max()), low_align(sim), mx / 1e3, my / 1e3))

# MOVING cascade 600 -> 200 -> 67 m, nests resolving the SURFACE LAYER (nz=64, z_stretch -> dz1~5 m)
NCX = 22; NEST_NZ = 64; NEST_ZS = 1.077
def mkspec(i):
    def build(gg):
        if i == 0:
            ic = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.xc)) - mx)))
            jc = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.yc)) - my)))
        else:
            ic, jc = gg.nx // 2, gg.ny // 2
        i0 = int(np.clip(ic - NCX // 2, 1, gg.nx - NCX - 1)); j0 = int(np.clip(jc - NCX // 2, 1, gg.ny - NCX - 1))
        return nst.NestSpec.aligned(gg, i0=i0, j0=j0, ncx=NCX, ncy=NCX, refine=3,
                                    nz=NEST_NZ, z_stretch=NEST_ZS)
    return build
WIN = 200.0
log("MOVING cascade 600 -> 200 -> 67 m, nests nz=%d z_stretch=%.3f (surface layer resolved), window %.0fs"
    % (NEST_NZ, NEST_ZS, WIN))
sims, rep = nst.run_multilevel_nest(sim, [mkspec(0), mkspec(1)], window=WIN,
                                    restrict_up=True, restrict_momentum=True,
                                    follow_interval=8, follow_field="zeta", follow_frac=0.4,
                                    follow_filter=0.5, les_boost=1.5, cfl=0.2,
                                    progress=lambda t, w, s: log("  finest t=%5.1f/%.0f step=%d" % (t, w, s))
                                    if s % 200 == 0 else None)
finest = sims[-1]; ng = finest.grid; to = ng.backend.to_cpu
finest.state.diagnose(finest.cfg)
uc, vc, wcn = rot._centered_velocity(finest.state, ng)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wcn)); zcn = np.asarray(to(ng.zc))
zeta = np.gradient(vc, ng.dx, axis=0) - np.gradient(uc, ng.dy, axis=1)
np.savez_compressed(os.path.join(OUT, "fields_intensity.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    rho=np.asarray(to(finest.state.rho)), p=np.asarray(to(finest.state.P_total)),
                    xc=np.asarray(to(ng.xc)), yc=np.asarray(to(ng.yc)), zc=zcn, dx=float(ng.dx))
log("fields saved; finest dx=%.0fm dz1=%.2fm nest_moves=%d" %
    (ng.dx, zcn[0], rep["nest"].get("nest_moves", 0)))

scon, vrep, crep, category, align_f = {}, {}, {}, "n/a", 0.0
try:
    scon = vd.surface_connection_report(finest.state, ng)
    for p in scon["profile"]:
        log("  z=%7.2fm  V_rot=%5.2f  |zeta|=%.2e  [%d cells from edge]" %
            (p["z_m"], p["v_rot_m_s"], p["zeta_max_s"], p["edge_cells"]))
    vrep = vd.vortex_report(finest.state, ng, z_m=max(20.0, float(zcn[0])))
    crep = cp.coldpool_report(finest.state, ng, z_m=max(20.0, float(zcn[0])))
    category = cl.classify_simulation(finest, z_surface_m=max(20.0, float(zcn[0])))["category"]
    align_f = low_align(finest)
    log("SURFACE: dz1=%.2fm ratio=%.2f connected=%s | V_theta=%.1f circ=%.2e dP=%.1f | class=%s align=%+.3f"
        % (scon["first_cell_height_m"], scon["surface_aloft_ratio"], scon["surface_connected"],
           vrep.get("v_theta_max_m_s", 0), vrep.get("circulation_m2_s", 0),
           vrep.get("pressure_deficit_Pa", 0), category, align_f))
except Exception as e:
    log("diagnostics error (fields saved): %r" % e)

summary = {"attempt": "K intensity", "env": {k: float(d[k]) for k in ("CAPE_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(b)), "finest_dx_m": float(ng.dx),
           "finest_dz1_m": float(zcn[0]), "nest_moves": rep["nest"].get("nest_moves", 0),
           "surface": scon, "vortex": vrep, "cold_pool": crep, "classification": category,
           "fine_align": align_f, "window_s": WIN,
           "reference": {"obs_Vrot": 26.0, "J_elevated_ratio": 0.0, "sens_fine_sdiv_loglaw": 0.82},
           "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
log("=== DONE K: dx=%.0fm dz1=%.2fm | sfc/aloft=%.2f connected=%s | V_theta=%.1f | class=%s (%.0fmin) ===" %
    (ng.dx, zcn[0], scon.get("surface_aloft_ratio", 0), scon.get("surface_connected", False),
     vrep.get("v_theta_max_m_s", 0), category, (time.time() - _t0) / 60))
