"""Does refinement intensify the surface vortex?  Asked WITHOUT nests.

Every attempt to answer this through the AMR cascade has produced an unusable measurement:
v1 contaminated (the 22 m vortex sat 111 m from its wall), v2 contaminated differently
(178 m, deeper inside the widened sponge), and the L2 fix -- which DID repair the intermediate's
placement (peak-in-sponge 7/14 -> 0/14) -- left the finest level worse, edge/interior 4.41 ->
21.41.  Measured separately: a nest's spurious edge vorticity does NOT track its boundary
divergence (direct-solver nest residual 1.3e-04 with zeta edge/int 323, against lowmem 2.8e-01
with 778 and 2.2e-01 with 217), so it is not a projection bug to fix -- it looks intrinsic to
Davies-zone nesting, which is exactly what `border_frac` exists to exclude.

The PARENT, by contrast, is clean and verified: edge/interior 0.18, |div(rho u)| 3.2e-05, and a
direct check showed its inflow is unmodified environmental air (1.8% of the 500 m level >1 K
cold; dtheta -0.23..+0.14 K from 5-60 km upstream).  So ask the question where the numerics are
trusted: run the SAME storm at two uniform resolutions with no nesting at all.

    RES=120 -> dx = 600 m   (the resolution every attempt A-L used)
    RES=240 -> dx = 300 m   (2x refinement, ~35 min to t=2800 on this GPU)

Identical in every other respect: same analytic environment, same periodic laterals (vindicated
by measurement), same storm-relative Bunkers frame, same stress-divergence + log-law surface
closure, same fixed 400 m comparison radius and 0.2 interior margin.

HONEST SCOPE: 300 m does not resolve a tornado, and this cannot replace the 22 m question.  What
it can do is establish the SIGN and MAGNITUDE of the resolution trend on numerics we trust,
which no nested run in this study has managed.  Chaos means the two storms differ in detail, so
the comparison is over the MATURE WINDOW (statistics), never a single instant.
"""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
os.chdir(REPO)
import numpy as np, warnings; warnings.filterwarnings("ignore")

RES = int(os.environ.get("RES", 120))
OUT = os.path.join(REPO, "outputs", "uniform_resolution"); os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("TAG", "%d" % RES)
LOG = os.path.join(OUT, "progress_%s.log" % TAG); open(LOG, "w").close()
_t0 = time.time()


def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m))
    print(m, flush=True)


from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, soundings as snd
from storm_dynamics import vortex_diagnostics as vd, classification as cl

# VERTICAL near-surface resolution -- the lever this project MEASURED to dominate (surface/aloft
# ratio 0.08 at dz1=49.6 m vs 0.53 at dz1=6.2 m, a 6.6x effect) while horizontal refinement gave
# only +22% for 9x (Attempt D).  The default parent grid puts ONE cell below 100 m (dz1 79.8 m),
# so the corner-flow layer that drives the surface vortex sits inside a single cell -- refining
# dx while that holds cannot change surface behaviour.  Measured cost at nx=120:
#     nz=48 zs=1.05 -> dz1 79.8 m, dt 2.12 s,  1 cell <100 m,    9 min to t=2800
#     nz=64 zs=1.09 -> dz1  5.5 m, dt 0.14 s, 11 cells <100 m, 195 min   <- affordable
#     nz=80 zs=1.09 -> dz1  1.4 m, dt 0.01 s, 23 cells <100 m,  65 h     <- not
# So dz1 5.5 m is runnable with a FULL maturation, which the earlier surface-sensitivity study
# never had (it refined inside a nest over a short window).  CAVEAT: nz and z_stretch both change
# between the two settings, so this is a near-surface-resolution comparison, not a single
# isolated variable -- z_stretch redistributes every level, not just the lowest.
NZ = int(os.environ.get("NZ", 48))
Z_STRETCH = float(os.environ.get("Z_STRETCH", 1.05))
LX = 72000.0
T_MAT = float(os.environ.get("T_MAT", 2800.0))
SAMPLE_FROM = float(os.environ.get("SAMPLE_FROM", 1800.0))   # mature phase only
SAMPLE_EVERY_S = float(os.environ.get("SAMPLE_EVERY_S", 50.0))
# COMPARISON RADIUS.  surface_connection_report converts this to CELLS with a floor,
#     R = max(3, int(radius_m / grid.dx))
# so a radius smaller than 3*dx silently becomes 3 cells -- a DIFFERENT physical disk at each
# resolution.  With 400 m it would have been 1800 m at dx=600 and 900 m at dx=300: a 2x mismatch
# in exactly the quantity under test, biasing the coarse run upward.  Pick a radius the floor
# cannot bind at the COARSEST mesh: 3000 m = 5 cells at 600 m and 10 cells at 300 m, the same
# physical disk in both.  (This is the "one fixed physical radius" rule, which the floor
# defeats unless radius_m >= 3*dx_coarsest.)
CMP_RADIUS_M = float(os.environ.get("CMP_RADIUS_M", 3000.0))
BORDER = 0.2
DEV = os.environ.get("DEV", "gpu")

log("=== uniform-resolution test, RES=%d (dx=%.0f m) nz=%d zs=%.3f, NO NESTS ==="
    % (RES, LX / RES, NZ, Z_STRETCH))

scfg = build_storm_config(preset="storm", nx=RES, ny=RES, nz=NZ, Lx=LX, Ly=LX, Lz=15000.0,
                          duration=1.0, dt_max=3.0, drag=True, z_stretch=Z_STRETCH, C_s=0.20,
                          hodograph_kind="quarter_circle", U_max=30.0, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
scfg.dyn.drag.stress_divergence = True
scfg.dyn.drag.surface_layer_depth_m = 150.0
scfg.dyn.drag.use_log_law = True
scfg.dyn.drag.roughness_length_m = 0.1

sim0 = StormSimulation(scfg); b = sim0.base
d = sounding_diagnostics(b)
cx, cy = snd.bunkers_storm_motion(b)
base_sr = BaseState(zc=b.zc, theta0=b.theta0, qv0=b.qv0, p0=b.p0, T0=b.T0, rho0=b.rho0,
                    u0=b.u0 - cx, v0=b.v0 - cy)
sim = StormSimulation(scfg, base=base_sr)

# ENSEMBLE SEED.  A supercell is chaotic, so a single realisation per configuration cannot
# distinguish an effect from run-to-run spread -- and the A/B pair (600 vs 300 m) currently rests
# on exactly one run each.  SEED>0 adds a physically negligible theta perturbation (amplitude
# SEED_K, default 0.01 K, ~0.003% of theta) whose only role is to send the run down a different
# chaotic trajectory.  SEED=0 reproduces the original member bit-for-bit.
SEED = int(os.environ.get("SEED", 0))
SEED_K = float(os.environ.get("SEED_K", 0.01))
if SEED:
    _rng = np.random.default_rng(SEED)
    _pert = _rng.normal(0.0, SEED_K, size=tuple(sim.state.theta.shape))
    sim.state.theta = sim.state.theta + sim.grid.xp.asarray(_pert)
    log("ENSEMBLE MEMBER seed=%d, theta perturbation amplitude %.3f K (chaotic divergence only)"
        % (SEED, SEED_K))
sim.state.diagnose(sim.cfg)
g = sim.grid
log("grid %dx%dx%d dx=%.1f m dz1=%.1f m periodic=%r | CAPE=%.0f shear06=%.1f | Bunkers=(%.1f,%.1f)"
    % (g.nx, g.ny, g.nz, g.dx, float(np.asarray(g.backend.to_cpu(g.zc))[0]) * 2,
       getattr(g, "periodic", None), d.get("CAPE_J_kg", float("nan")),
       d.get("shear_0_6km_m_s", float("nan")), cx, cy))

series = []
_next = {"t": SAMPLE_FROM}


def sample(t):
    """Low-level rotation, scored identically at both resolutions."""
    try:
        zc = np.asarray(g.backend.to_cpu(g.zc))
        sc_ = vd.surface_connection_report(sim.state, g, border_frac=BORDER,
                                           radius_m=CMP_RADIUS_M)
        vr = vd.vortex_report(sim.state, g, z_m=max(20.0, float(zc[0])), border_frac=BORDER,
                              radius_m=CMP_RADIUS_M)
        prof = sc_["profile"]
        rec = {"t": float(t), "dx_m": float(g.dx),
               "v_rot_sfc": prof[0]["v_rot_m_s"], "zeta_sfc": prof[0]["zeta_max_s"],
               "ratio": sc_["surface_aloft_ratio"], "connected": bool(sc_["surface_connected"]),
               "v_theta": vr.get("v_theta_max_m_s") or 0.0,
               "circ": vr.get("circulation_m2_s"), "core_m": vr.get("core_radius_m"),
               "profile": [{"z": p["z_m"], "v": p["v_rot_m_s"], "zeta": p["zeta_max_s"]}
                           for p in prof]}
        series.append(rec)
        log("  t=%6.1f dx=%5.1fm V_sfc=%6.2f zeta=%.4f r=%.2f Vth=%5.2f circ=%.2e core=%sm"
            % (t, g.dx, rec["v_rot_sfc"], rec["zeta_sfc"], rec["ratio"], rec["v_theta"],
               rec["circ"] or 0.0, rec["core_m"]))
        json.dump(series, open(os.path.join(OUT, "series_%s.json" % TAG), "w"), indent=1)
    except Exception as e:
        log("  sample error t=%.1f: %r" % (t, e))


n = 0
while sim.t < T_MAT:
    dt = float(sim._dt())
    if sim.t + dt > T_MAT:
        dt = T_MAT - sim.t
    sim._step(dt)
    sim.step += 1; sim.t = float(getattr(sim.state, "t", sim.t + dt)); n += 1
    if sim.t >= _next["t"]:
        sample(sim.t)
        _next["t"] += SAMPLE_EVERY_S
    elif n % 300 == 0:
        u = np.abs(np.asarray(g.backend.to_cpu(sim.state.u)))
        w = np.abs(np.asarray(g.backend.to_cpu(sim.state.w)))
        log("  t=%6.1f (spin-up) max|u|=%5.2f max|w|=%5.2f" % (sim.t, u.max(), w.max()))

if series:
    pk = max(r["v_rot_sfc"] for r in series); mn = float(np.mean([r["v_rot_sfc"] for r in series]))
    pz = max(r["zeta_sfc"] for r in series); mz = float(np.mean([r["zeta_sfc"] for r in series]))
    pc = max(abs(r["circ"] or 0.0) for r in series)
    log("AGG dx=%.0fm over t=%.0f..%.0f (n=%d): peak V_sfc=%.2f mean=%.2f | peak |zeta|=%.4f "
        "mean=%.4f | peak |circ|=%.2e | connected %.0f%%"
        % (g.dx, series[0]["t"], series[-1]["t"], len(series), pk, mn, pz, mz, pc,
           100.0 * sum(1 for r in series if r["connected"]) / len(series)))
    json.dump({"res": RES, "dx_m": float(g.dx), "n_samples": len(series),
               "peak_v_rot_sfc": pk, "mean_v_rot_sfc": mn, "peak_zeta_sfc": pz,
               "mean_zeta_sfc": mz, "peak_abs_circ": pc,
               "connected_fraction": sum(1 for r in series if r["connected"]) / len(series),
               "classification": cl.classify_simulation(sim, z_surface_m=20.0)["category"],
               "series": series},
              open(os.path.join(OUT, "summary_%s.json" % TAG), "w"), indent=1)
log("DONE RES=%d dx=%.0fm in %.0f min" % (RES, g.dx, (time.time() - _t0) / 60))
