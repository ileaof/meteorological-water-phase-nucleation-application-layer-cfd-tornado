"""Attempt L -- INTENSITY.  Attempt K solved the STRUCTURE (a surface-connected, surface-intensified
vortex: V_rot 7.76 at z=5.1 m decreasing upward, ratio 0.82, connected=True) but not the INTENSITY:
7.76 m/s against the observed KTLX 26 m/s, and the surface feature was still a linear shear BAND
(strongest anticyclonic) rather than a compact axisymmetric core.

Three changes, each aimed at that gap:

  1. A THIRD cascade level: 600 -> 200 -> 67 -> 22 m.  A 0.3 km observed TVS core spans ~4 cells at
     67 m -- far too few to hold an axisymmetric vortex against numerical diffusion.  At 22 m it
     spans ~14.  This is the resolution axis Attempt K could not test.
  2. A LONGER window with PERSISTENCE tracking (`sample` hook, new): the vortex is measured every
     few parent steps, not only at the end, so a lasting vortex is distinguishable from a
     transient -- persistence is one of the completion criteria and was never actually measured.
  3. PRESSURE DEFICIT, now measurable at all.  The low-memory projection discarded the pressure
     potential, so `state.p` was stale on every nest and every deficit read exactly 0 -- which
     silently disabled the classifier tiers that require one.  Fixed; Delta p is a real number here
     for the first time on a nest.

Scored against the observed Moore 2013 KTLX couplet: V_rot 26 m/s, deltaV 52, vorticity 0.21 /s,
core ~0.3 km.  Honest reporting: a shear band that fails to close into a core is reported as such.
"""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "tornado_intensity_L"); os.makedirs(OUT, exist_ok=True)
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

BORDER = 0.2          # interior margin: MUST exclude the nest relaxation zone (Attempt K artifact)
OBS = {"v_rot": 26.0, "delta_v": 52.0, "zeta": 0.21, "core_km": 0.3}

def low_align(sim, z_top=800.0):
    g = sim.grid; to = g.backend.to_cpu
    uc, vc, wc = rot._centered_velocity(sim.state, g)
    e = vb.tilting_efficiency(uc, vc, wc, g); zc = np.asarray(to(g.zc)); nb = g.nx // 6; km = zc < z_top
    c = lambda a: np.asarray(to(a))[nb:-nb, nb:-nb][:, :, km]
    return float(c(e["tilting"]).sum() / (np.abs(c(e["omega_h"]) * c(e["grad_h_w"])).sum() + 1e-20))

DEV = os.environ.get("DEV", "gpu")
log("=== Attempt L: INTENSITY -- 3-level cascade to 22 m, persistence, real pressure deficit ===")
nx = 120; nz = 48; Lx = 72000.0
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=15000.0,
                          duration=1.0, dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20,
                          hodograph_kind="quarter_circle", U_max=30.0, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
scfg.dyn.drag.stress_divergence = True          # the closure that restored surface connection
scfg.dyn.drag.surface_layer_depth_m = 150.0
scfg.dyn.drag.use_log_law = True
scfg.dyn.drag.roughness_length_m = 0.1

sim0 = StormSimulation(scfg); b = sim0.base
d = sounding_diagnostics(b); cx, cy = snd.bunkers_storm_motion(b)
log("env: CAPE=%.0f shear06=%.1f SRH03=%.0f | parent dx=%.0fm | drag=stress-div+log-law" %
    (d["CAPE_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(b), sim0.grid.dx))
base_sr = BaseState(zc=b.zc, theta0=b.theta0, qv0=b.qv0, p0=b.p0, T0=b.T0, rho0=b.rho0,
                    u0=b.u0 - cx, v0=b.v0 - cy)
sim = StormSimulation(scfg, base=base_sr); sim.state.diagnose(sim.cfg)

T_MAT = 2800.0
# The maturation is 13 minutes of identical work on every attempt, and the cascade is what is
# actually being varied -- cache the matured parent so a failed cascade can be retried cheaply.
CACHE = os.path.join(REPO, "outputs", "parent_matured_%d_%d_%d.npz" % (nx, nz, int(T_MAT)))

def _state_arrays(st):
    return {k: v for k, v in vars(st).items() if isinstance(v, (np.ndarray,))
            or type(v).__module__.split(".", 1)[0] == "cupy"}

if os.path.exists(CACHE):
    z = np.load(CACHE)
    xpn = sim.grid.xp
    for k in z.files:
        if k == "_t":
            continue
        setattr(sim.state, k, xpn.asarray(z[k]))
    sim.t = float(z["_t"][0]); sim.state.t = sim.t
    sim.state.diagnose(sim.cfg)
    log("loaded cached matured parent t=%.0f from %s" % (sim.t, os.path.basename(CACHE)))
else:
    log("maturing the freely-evolving supercell to its low-level-meso peak ...")
    while sim.t < T_MAT:
        dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t)
        if sim.step % 400 == 0:
            w = np.asarray(sim.grid.backend.to_cpu(sim.state.w))
            if not np.isfinite(w).all():
                log("  instability at t=%.0f -> stop" % sim.t); break
            sim.state.diagnose(sim.cfg)
            log("  parent t=%4.0f w_max=%5.1f align=%+.3f" % (sim.t, w.max(), low_align(sim)))
    tc = sim.grid.backend.to_cpu
    np.savez_compressed(CACHE, _t=np.array([sim.t]),
                        **{k: np.asarray(tc(v)) for k, v in _state_arrays(sim.state).items()})
    log("cached matured parent -> %s" % os.path.basename(CACHE))

to = sim.grid.backend.to_cpu; uc, vc, _ = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); zcp = np.asarray(to(sim.grid.zc))
nbp = nx // 6
def _zeta_peak_xy(z_target):
    """(x, y) of the interior |zeta| maximum at the level nearest z_target, on the parent."""
    k = int(np.argmin(np.abs(zcp - z_target)))
    zl = np.abs(np.gradient(vc[:, :, k], sim.grid.dx, axis=0)
                - np.gradient(uc[:, :, k], sim.grid.dy, axis=1))
    a, bq = np.unravel_index(np.argmax(zl[nbp:-nbp, nbp:-nbp]), zl[nbp:-nbp, nbp:-nbp].shape)
    a += nbp; bq += nbp
    return float(sim.grid.xc[a]), float(sim.grid.yc[bq]), float(zl[a, bq])

mx, my, zmid = _zeta_peak_xy(500.0)          # mid-level mesocyclone (what J and K centred on)
sx, sy, zsfc = _zeta_peak_xy(100.0)          # the LOW-LEVEL circulation -- the target here
sep_km = float(np.hypot(sx - mx, sy - my)) / 1e3
log("PARENT matured: w_max=%.1f align=%+.3f | meso(500m)@(%.1f,%.1f)km zeta=%.3f | "
    "low-level(100m)@(%.1f,%.1f)km zeta=%.3f | SEPARATION %.2f km" %
    (float(np.abs(np.asarray(to(sim.state.w))).max()), low_align(sim), mx / 1e3, my / 1e3, zmid,
     sx / 1e3, sy / 1e3, zsfc, sep_km))
# Centre EVERY level on the low-level circulation.  Centring the finer levels on their parent
# grid's geometric centre (i.e. on the mid-level meso) left the 2 km / 22 m nest looking at quiet
# air -- |zeta| 0.002 against 0.238 on its own 67 m parent.  The surface vortex is the subject of
# this run, so it is what the mesh must be placed on.

# ---- 3-level MOVING cascade 600 -> 200 -> 67 -> 22 m, surface layer resolved on every nest ----
LEVELS = int(os.environ.get("LEVELS", 3))   # 3 -> 22 m finest; 2 reproduces Attempt K's 67 m
NCX = [22, 22, 30][:LEVELS]     # parent cells covered per level; the finest is WIDER (2.0 km) so a
NEST_NZ = 64; NEST_ZS = 1.077   # 0.3 km core plus its inflow fits inside the 0.2 interior margin
def mkspec(i):
    def build(gg):
        ic = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.xc)) - sx)))
        jc = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.yc)) - sy)))
        n = NCX[i]
        i0 = int(np.clip(ic - n // 2, 1, gg.nx - n - 1)); j0 = int(np.clip(jc - n // 2, 1, gg.ny - n - 1))
        return nst.NestSpec.aligned(gg, i0=i0, j0=j0, ncx=n, ncy=n, refine=3,
                                    nz=NEST_NZ, z_stretch=NEST_ZS)
    return build

WIN = float(os.environ.get("WIN", 240.0))
SAMPLE_EVERY = 4                              # parent steps between diagnostic samples
# ONE physical radius for every level: the finest domain is 2.0 km and the 0.2 interior margin
# leaves ~1.2 km, so 400 m fits inside every level's trusted interior.  Keeping it fixed is what
# makes "same storm, same instant, only dx differs" true.
CMP_RADIUS_M = 400.0
series = []
_state = {"n": 0}

def _diag_level(sim_l, t):
    """Diagnose ONE nest level.  The measurement radius is the SAME PHYSICAL LENGTH on every
    level -- that is what makes the cross-level comparison a controlled experiment.

    (First version scaled the radius to each level's own domain, which silently broke the
    comparison: V_rot is the peak deviation within the radius, so a larger radius sweeps in more
    of the broad flow and reports a larger "rotation" for purely geometric reasons.)"""
    g = sim_l.grid
    z1 = float(np.asarray(g.backend.to_cpu(g.zc))[0])
    width = float(g.dx) * g.nx
    rad = CMP_RADIUS_M
    sc = vd.surface_connection_report(sim_l.state, g, border_frac=BORDER, radius_m=rad)
    vr = vd.vortex_report(sim_l.state, g, z_m=max(20.0, z1), border_frac=BORDER, radius_m=rad)
    return {"t": float(t), "dx_m": float(g.dx), "z1_m": z1, "width_km": width / 1e3,
            "radius_m": rad,
            "v_rot_sfc": sc["profile"][0]["v_rot_m_s"],
            "zeta_sfc": sc["profile"][0]["zeta_max_s"],
            "ratio": sc["surface_aloft_ratio"], "connected": bool(sc["surface_connected"]),
            "conv": sc["near_surface_convergence_s"],
            "v_theta": vr.get("v_theta_max_m_s"), "dP": vr.get("pressure_deficit_Pa"),
            "circ": vr.get("circulation_m2_s"), "core_m": vr.get("core_radius_m")}


def sample(sims, t):
    """Persistence tracking on EVERY nest level (new `sample` hook).

    Diagnosing all levels at the same instant turns the run into a controlled resolution
    experiment: identical storm, identical time, only dx differs.  That is the direct test of
    whether a finer horizontal mesh intensifies the surface vortex -- the open question left by
    Attempt K, which could only compare across separate runs."""
    _state["n"] += 1
    if _state["n"] % SAMPLE_EVERY:
        return
    try:
        recs = [_diag_level(s, t) for s in sims[1:]]
        series.append(recs)
        log("  t=%5.1f | " % t + " | ".join(
            "dx=%4.0fm V_sfc=%5.2f z=%.3f r=%.2f dP=%7.1f core=%4.0fm" %
            (r["dx_m"], r["v_rot_sfc"], r["zeta_sfc"], r["ratio"],
             r["dP"] if r["dP"] == r["dP"] else float("nan"), r["core_m"] or 0) for r in recs))
        json.dump(series, open(os.path.join(OUT, "series.json"), "w"), indent=1)
    except Exception as e:
        log("  sample error t=%.1f: %r" % (t, e))

log("MOVING cascade %d levels (NCX=%s), nests nz=%d z_stretch=%.3f, window %.0fs"
    % (LEVELS, NCX, NEST_NZ, NEST_ZS, WIN))
sims, rep = nst.run_multilevel_nest(sim, [mkspec(i) for i in range(LEVELS)], window=WIN,
                                    restrict_up=True, restrict_momentum=True,
                                    follow_interval=8, follow_field="zeta", follow_frac=0.4,
                                    follow_filter=0.5, follow_z_lo=0.0, follow_z_hi=1500.0,
                                    les_boost=1.5, cfl=0.2, sample=sample,
                                    progress=lambda t, w, s: log("  finest t=%5.1f/%.0f step=%d" % (t, w, s))
                                    if s % 400 == 0 else None)
finest = sims[-1]; ng = finest.grid; to = ng.backend.to_cpu
finest.state.diagnose(finest.cfg)
uc, vc, wcn = rot._centered_velocity(finest.state, ng)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); wcn = np.asarray(to(wcn)); zcn = np.asarray(to(ng.zc))
zeta = np.gradient(vc, ng.dx, axis=0) - np.gradient(uc, ng.dy, axis=1)
np.savez_compressed(os.path.join(OUT, "fields_L.npz"), u=uc, v=vc, w=wcn, zeta=zeta,
                    rho=np.asarray(to(finest.state.rho)), p=np.asarray(to(finest.state.p)),
                    xc=np.asarray(to(ng.xc)), yc=np.asarray(to(ng.yc)), zc=zcn, dx=float(ng.dx))
log("fields saved; finest dx=%.1fm dz1=%.2fm nest_moves=%d" %
    (ng.dx, zcn[0], rep["nest"].get("nest_moves", 0)))

scon, vrep, crep, category, align_f = {}, {}, {}, "n/a", 0.0
try:
    z1 = max(20.0, float(zcn[0]))
    _rad = CMP_RADIUS_M
    scon = vd.surface_connection_report(finest.state, ng, border_frac=BORDER, radius_m=_rad)
    for p in scon["profile"]:
        log("  z=%7.2fm  V_rot=%5.2f  |zeta|=%.2e  [%d cells from edge]" %
            (p["z_m"], p["v_rot_m_s"], p["zeta_max_s"], p["edge_cells"]))
    vrep = vd.vortex_report(finest.state, ng, z_m=z1, border_frac=BORDER, radius_m=_rad)
    crep = cp.coldpool_report(finest.state, ng, z_m=z1)
    category = cl.classify_simulation(finest, z_surface_m=z1)["category"]
    align_f = low_align(finest)
    log("SURFACE: dz1=%.2fm ratio=%.2f connected=%s | V_theta=%.1f circ=%.2e dP=%.1f Pa | class=%s"
        % (scon["first_cell_height_m"], scon["surface_aloft_ratio"], scon["surface_connected"],
           vrep.get("v_theta_max_m_s", 0), vrep.get("circulation_m2_s", 0),
           vrep.get("pressure_deficit_Pa", float("nan")), category))
except Exception as e:
    log("diagnostics error (fields saved): %r" % e)

# persistence + the resolution comparison, per level (series entries are one list per sample)
pers = {}
if series:
    nlev = len(series[0])
    for li in range(nlev):
        col = [s[li] for s in series if len(s) > li]
        con = [c for c in col if c["connected"]]
        dxl = col[0]["dx_m"]
        # longest CONSECUTIVE connected run -> the honest persistence time
        best = cur = 0.0
        for a, c in zip(col, col[1:]):
            cur = cur + (c["t"] - a["t"]) if c["connected"] else 0.0
            best = max(best, cur)
        pers["dx_%.0fm" % dxl] = {
            "dx_m": dxl, "width_km": col[0]["width_km"], "samples": len(col),
            "connected_fraction": len(con) / max(1, len(col)),
            "longest_connected_s": best,
            "peak_v_rot_sfc": max(c["v_rot_sfc"] for c in col),
            "mean_v_rot_sfc": float(np.mean([c["v_rot_sfc"] for c in col])),
            "peak_zeta_sfc": max(c["zeta_sfc"] for c in col),
            "peak_v_theta": max((c["v_theta"] or 0) for c in col),
            "min_dP_Pa": min((c["dP"] for c in col if c["dP"] == c["dP"]), default=float("nan")),
            "span_s": col[-1]["t"] - col[0]["t"]}
        p = pers["dx_%.0fm" % dxl]
        log("PERSISTENCE dx=%4.0fm (%.1fkm dom): %d samples/%.0fs | connected %.0f%% "
            "(longest run %.0fs) | peak V_sfc=%.2f Vth=%.2f zeta=%.3f | min dP=%.1f Pa"
            % (dxl, p["width_km"], p["samples"], p["span_s"], 100 * p["connected_fraction"],
               p["longest_connected_s"], p["peak_v_rot_sfc"], p["peak_v_theta"],
               p["peak_zeta_sfc"], p["min_dP_Pa"]))
    fine = pers.get("dx_%.0fm" % series[0][-1]["dx_m"], {})
    coarse = pers.get("dx_%.0fm" % series[0][0]["dx_m"], {})
    if fine and coarse and coarse.get("peak_v_rot_sfc"):
        log("RESOLUTION EFFECT (same storm, same instants): peak V_sfc %.2f at dx=%.0fm vs "
            "%.2f at dx=%.0fm -> x%.2f"
            % (fine["peak_v_rot_sfc"], fine["dx_m"], coarse["peak_v_rot_sfc"], coarse["dx_m"],
               fine["peak_v_rot_sfc"] / coarse["peak_v_rot_sfc"]))

summary = {"attempt": "L intensity (3-level, 22 m, persistence, real dP)",
           "env": {k: float(d[k]) for k in ("CAPE_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(b)), "finest_dx_m": float(ng.dx),
           "finest_dz1_m": float(zcn[0]), "nest_moves": rep["nest"].get("nest_moves", 0),
           "levels": rep["nest"].get("levels"), "window_s": WIN,
           "surface": scon, "vortex": vrep, "cold_pool": crep, "classification": category,
           "fine_align": align_f, "persistence": pers, "series": series,
           "observed_KTLX": OBS,
           "reference": {"K_dx_m": 67.0, "K_v_rot_sfc": 7.76, "K_ratio": 0.82},
           "note": "series entries are one list per sample, one record per nest level -- the "
                   "levels are diagnosed at the SAME instant, so dx is the only variable",
           "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
_fin = pers.get("dx_%.0fm" % float(ng.dx), {})
log("=== DONE L: dx=%.1fm | sfc/aloft=%.2f connected=%s | V_sfc peak=%.2f (K:7.76, obs:26) | "
    "dP=%.1f Pa | class=%s (%.0f min) ===" %
    (ng.dx, scon.get("surface_aloft_ratio", 0), scon.get("surface_connected", False),
     _fin.get("peak_v_rot_sfc", float("nan")), vrep.get("pressure_deficit_Pa", float("nan")),
     category, (time.time() - _t0) / 60))
