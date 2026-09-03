"""Matched-domain resolution test -- the controlled version of Attempt L's cross-level comparison.

Attempt L diagnosed every cascade level at the same instants, but a nest must sit INSIDE its
parent, so each finer level was also a smaller box (13.2 -> 4.4 -> 2.0 km).  dx and domain size
were confounded, and "finer mesh, weaker surface rotation" could not be concluded from it: with a
0.2 border margin the 22 m level had only ~1.2 km of trusted interior, so the vortex's feeding
inflow was boundary-controlled.

Here the FINEST domain is the SAME 4.0 km in both runs and only dx differs:

    MODE=coarse : 600 -> 200 -> 67 m , finest =  60 cells over 4.0 km
    MODE=fine   : 600 -> 200 -> 67 -> 22 m , finest = 180 cells over 4.0 km

A first attempt matched the domains at 2.0 km and the coarse branch NaN-ed on its first step: that
made the 67 m level only 30 cells wide, so with a 0.2 relaxation zone on each side the boundary
regions effectively meet and there is no interior left to evolve.  Every nest that has run stably
in this study is 60+ cells across.  4.0 km is the smallest domain that keeps BOTH branches in that
range, which is what makes the comparison possible at all.

Both branch from the same cached matured parent at the same instant, are centred on the same
low-level circulation, and are scored with the same fixed 400 m comparison radius and the same 0.2
interior margin.  The question is narrow and answerable: at a fixed physical domain, does resolving
the vortex at 22 m instead of 67 m intensify the surface rotation, or merely sharpen a shear band?

The window is short (intensity, not persistence -- Attempt L measures persistence).
"""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

MODE = os.environ.get("MODE", "fine").lower()
OUT = os.path.join(REPO, "outputs", "matched_domain"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress_%s.log" % MODE); open(LOG, "w").close()
_t0 = time.time()
def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m)); print(m, flush=True)

from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, soundings as snd, nesting as nst
from storm_dynamics import vortex_diagnostics as vd, classification as cl

BORDER = 0.2
CMP_RADIUS_M = 400.0          # identical measurement radius in both runs
TARGET_KM = 4.0               # identical finest-domain width in both runs

nx = 120; nz = 48; Lx = 72000.0; T_MAT = 2800.0
DEV = os.environ.get("DEV", "gpu")
log("=== matched-domain test, MODE=%s (finest domain %.1f km in BOTH modes) ===" % (MODE, TARGET_KM))

scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=15000.0,
                          duration=1.0, dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20,
                          hodograph_kind="quarter_circle", U_max=30.0, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
scfg.dyn.drag.stress_divergence = True
scfg.dyn.drag.surface_layer_depth_m = 150.0
scfg.dyn.drag.use_log_law = True
scfg.dyn.drag.roughness_length_m = 0.1

sim0 = StormSimulation(scfg); b = sim0.base
d = sounding_diagnostics(b); cx, cy = snd.bunkers_storm_motion(b)
base_sr = BaseState(zc=b.zc, theta0=b.theta0, qv0=b.qv0, p0=b.p0, T0=b.T0, rho0=b.rho0,
                    u0=b.u0 - cx, v0=b.v0 - cy)
sim = StormSimulation(scfg, base=base_sr); sim.state.diagnose(sim.cfg)

CACHE = os.path.join(REPO, "outputs", "parent_matured_%d_%d_%d.npz" % (nx, nz, int(T_MAT)))
if not os.path.exists(CACHE):
    log("FATAL: cached matured parent missing (%s) -- run tornado_intensity_L_gpu.py once" % CACHE)
    sys.exit(1)
z = np.load(CACHE); xpn = sim.grid.xp
for k in z.files:
    if k != "_t":
        setattr(sim.state, k, xpn.asarray(z[k]))
sim.t = float(z["_t"][0]); sim.state.t = sim.t; sim.state.diagnose(sim.cfg)
log("loaded cached matured parent t=%.0f" % sim.t)

to = sim.grid.backend.to_cpu; uc, vc, _ = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); zcp = np.asarray(to(sim.grid.zc)); nbp = nx // 6
k100 = int(np.argmin(np.abs(zcp - 100.0)))
zl = np.abs(np.gradient(vc[:, :, k100], sim.grid.dx, axis=0)
            - np.gradient(uc[:, :, k100], sim.grid.dy, axis=1))
a, bq = np.unravel_index(np.argmax(zl[nbp:-nbp, nbp:-nbp]), zl[nbp:-nbp, nbp:-nbp].shape)
a += nbp; bq += nbp
sx, sy = float(sim.grid.xc[a]), float(sim.grid.yc[bq])
log("low-level circulation @(%.1f,%.1f)km zeta=%.4f" % (sx / 1e3, sy / 1e3, float(zl[a, bq])))

# NCX per level chosen so the FINEST level is TARGET_KM wide in both modes.
#   coarse: 200 m level 18 km, then 67 m level over 20 x 200 m = 4.0 km   (60 cells)
#   fine  : 200 m level 18 km, 67 m level over 30 x 200 m = 6.0 km, then 22 m over
#           60 x 67 m = 4.0 km                                            (180 cells)
NCX = [30, 20] if MODE == "coarse" else [30, 30, 60]
NEST_NZ = 64; NEST_ZS = 1.077
def mkspec(i):
    def build(gg):
        ic = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.xc)) - sx)))
        jc = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.yc)) - sy)))
        n = NCX[i]
        i0 = int(np.clip(ic - n // 2, 1, gg.nx - n - 1)); j0 = int(np.clip(jc - n // 2, 1, gg.ny - n - 1))
        return nst.NestSpec.aligned(gg, i0=i0, j0=j0, ncx=n, ncy=n, refine=3,
                                    nz=NEST_NZ, z_stretch=NEST_ZS)
    return build

WIN = float(os.environ.get("WIN", 60.0))
SAMPLE_EVERY = 3
series = []; _n = {"i": 0}

def sample(sims, t):
    _n["i"] += 1
    if _n["i"] % SAMPLE_EVERY:
        return
    g = sims[-1].grid
    try:
        z1 = float(np.asarray(g.backend.to_cpu(g.zc))[0])
        sc = vd.surface_connection_report(sims[-1].state, g, border_frac=BORDER,
                                          radius_m=CMP_RADIUS_M)
        vr = vd.vortex_report(sims[-1].state, g, z_m=max(20.0, z1), border_frac=BORDER,
                              radius_m=CMP_RADIUS_M)
        vt = vr.get("v_theta_max_m_s") or 0.0
        dp = vr.get("pressure_deficit_Pa")
        rec = {"t": float(t), "dx_m": float(g.dx), "width_km": float(g.dx) * g.nx / 1e3,
               "z1_m": z1, "v_rot_sfc": sc["profile"][0]["v_rot_m_s"],
               "zeta_sfc": sc["profile"][0]["zeta_max_s"], "ratio": sc["surface_aloft_ratio"],
               "connected": bool(sc["surface_connected"]), "v_theta": vt, "dP": dp,
               "dP_cyclostrophic": -1.2 * vt * vt,      # the physical scale to judge dP against
               "circ": vr.get("circulation_m2_s"), "core_m": vr.get("core_radius_m")}
        series.append(rec)
        log("  t=%5.1f dx=%5.1fm dom=%.2fkm V_sfc=%5.2f zeta=%.3f r=%.2f Vth=%5.2f "
            "dP=%8.1f (cyclo %7.1f) core=%4.0fm"
            % (t, rec["dx_m"], rec["width_km"], rec["v_rot_sfc"], rec["zeta_sfc"], rec["ratio"],
               vt, dp if dp == dp else float("nan"), rec["dP_cyclostrophic"], rec["core_m"] or 0))
        json.dump(series, open(os.path.join(OUT, "series_%s.json" % MODE), "w"), indent=1)
    except Exception as e:
        log("  sample error t=%.1f: %r" % (t, e))

log("cascade %s: NCX=%s -> finest dx target %s, window %.0fs"
    % (MODE, NCX, "67 m" if MODE == "coarse" else "22 m", WIN))
sims, rep = nst.run_multilevel_nest(sim, [mkspec(i) for i in range(len(NCX))], window=WIN,
                                    restrict_up=True, restrict_momentum=True,
                                    follow_interval=8, follow_field="zeta", follow_frac=0.4,
                                    follow_filter=0.5, follow_z_lo=0.0, follow_z_hi=1500.0,
                                    les_boost=1.5, cfl=0.2, sample=sample)
finest = sims[-1]; ng = finest.grid; tog = ng.backend.to_cpu
finest.state.diagnose(finest.cfg)
ucn, vcn, wcn = rot._centered_velocity(finest.state, ng)
ucn = np.asarray(tog(ucn)); vcn = np.asarray(tog(vcn)); wcn = np.asarray(tog(wcn))
zcn = np.asarray(tog(ng.zc))
zeta = np.gradient(vcn, ng.dx, axis=0) - np.gradient(ucn, ng.dy, axis=1)
np.savez_compressed(os.path.join(OUT, "fields_%s.npz" % MODE), u=ucn, v=vcn, w=wcn, zeta=zeta,
                    p=np.asarray(tog(getattr(finest.state, "p_dyn", finest.state.p))),
                    xc=np.asarray(tog(ng.xc)), yc=np.asarray(tog(ng.yc)), zc=zcn, dx=float(ng.dx))

z1 = max(20.0, float(zcn[0]))
scon = vd.surface_connection_report(finest.state, ng, border_frac=BORDER, radius_m=CMP_RADIUS_M)
vrep = vd.vortex_report(finest.state, ng, z_m=z1, border_frac=BORDER, radius_m=CMP_RADIUS_M)
category = cl.classify_simulation(finest, z_surface_m=z1)["category"]
agg = {}
if series:
    agg = {"samples": len(series),
           "peak_v_rot_sfc": max(r["v_rot_sfc"] for r in series),
           "mean_v_rot_sfc": float(np.mean([r["v_rot_sfc"] for r in series])),
           "peak_zeta_sfc": max(r["zeta_sfc"] for r in series),
           "mean_zeta_sfc": float(np.mean([r["zeta_sfc"] for r in series])),
           "peak_v_theta": max(r["v_theta"] for r in series),
           "connected_fraction": sum(1 for r in series if r["connected"]) / len(series)}
    log("AGG %s dx=%.1fm dom=%.2fkm: peak V_sfc=%.2f mean=%.2f | peak zeta=%.3f mean=%.3f | "
        "peak Vth=%.2f | connected %.0f%%"
        % (MODE, ng.dx, float(ng.dx) * ng.nx / 1e3, agg["peak_v_rot_sfc"], agg["mean_v_rot_sfc"],
           agg["peak_zeta_sfc"], agg["mean_zeta_sfc"], agg["peak_v_theta"],
           100 * agg["connected_fraction"]))

json.dump({"mode": MODE, "finest_dx_m": float(ng.dx), "finest_width_km": float(ng.dx) * ng.nx / 1e3,
           "finest_nx": ng.nx, "first_cell_m": float(zcn[0]), "window_s": WIN,
           "comparison_radius_m": CMP_RADIUS_M, "border_frac": BORDER,
           "aggregate": agg, "surface": scon, "vortex": vrep, "classification": category,
           "series": series, "nest_moves": rep["nest"].get("nest_moves", 0),
           "wall_clock_s": time.time() - _t0},
          open(os.path.join(OUT, "summary_%s.json" % MODE), "w"), indent=2)
log("=== DONE %s: dx=%.1fm over %.2fkm | peak V_sfc=%.2f | ratio=%.2f | class=%s (%.0f min) ==="
    % (MODE, ng.dx, float(ng.dx) * ng.nx / 1e3, agg.get("peak_v_rot_sfc", float("nan")),
       scon["surface_aloft_ratio"], category, (time.time() - _t0) / 60))
