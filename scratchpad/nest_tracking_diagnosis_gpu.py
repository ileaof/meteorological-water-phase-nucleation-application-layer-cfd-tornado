"""WHICH LEVEL loses the vortex? -- a cheap diagnosis before paying for another full window.

The 4.0 km matched-domain pair could not answer the resolution question because the 22 m branch
never held the low-level vortex in trusted interior.  Fixing the sponge width and the moving-nest
tracker (4e4a1ce) did NOT fix it: measured on the v2 final field, edge/interior mean |zeta| got
WORSE (4.98 -> 5.61) and the one coherent vortex moved from 111 m to 178 m from the boundary --
i.e. deeper INSIDE the now-267 m sponge.  The 67 m branch stays healthy (edge/interior 0.77,
vortex 933 m in).

So the dominant defect is neither the sponge width nor the tracker: the 3-level fine cascade
cannot keep the vortex interior, and we do not know WHICH level drops it -- the 22 m box, or the
6 km / 67 m intermediate whose placement it inherits.  A full window costs 317 min; this run is
SHORT and instruments EVERY level each sample:

  * where the low-level |zeta| peak sits in that level's own box (metres from the nearest edge,
    and as a fraction of the half-width),
  * whether the peak is inside that level's sponge band,
  * the level's edge/interior mean |zeta| ratio (the contamination measure),
  * the box centre, so a drifting box shows up as a trajectory.

Read the output as: the FIRST level whose peak crosses into its sponge is the one to fix.
Diagnosis only -- no physics changed.
"""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "tracking_diagnosis_%s" % os.environ.get("TAG","orig")); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress.log"); open(LOG, "w").close()
_t0 = time.time()


def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m))
    print(m, flush=True)


from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, soundings as snd, nesting as nst

WIN = float(os.environ.get("WIN", 12.0))          # SHORT -- diagnosis, not a result
SAMPLE_EVERY = int(os.environ.get("SAMPLE_EVERY", 2))
SPONGE_M = float(os.environ.get("SPONGE_M", 267.0))
FOLLOW_BORDER = os.environ.get("FOLLOW_BORDER", "auto")
FOLLOW_FILTER = float(os.environ.get("FOLLOW_FILTER", 0.7))
# NCX_ENV lets the L2 fix be tested without editing the script: [30,30,60] is the ORIGINAL
# ladder (67 m intermediate = 6 km, the one that loses the vortex), [30,50,60] widens that
# intermediate to 10 km so the 4 km child sits inside its trusted interior.
NCX = [int(x) for x in os.environ.get("NCX", "30,30,60").split(",")]
NEST_NZ = 64
NEST_ZS = 1.077

nx = 120; nz = 48; Lx = 72000.0; T_MAT = 2800.0
DEV = os.environ.get("DEV", "gpu")
log("=== nest-tracking diagnosis: which level loses the vortex? ===")

scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=15000.0,
                          duration=1.0, dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20,
                          hodograph_kind="quarter_circle", U_max=30.0, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
scfg.dyn.drag.stress_divergence = True
scfg.dyn.drag.surface_layer_depth_m = 150.0
scfg.dyn.drag.use_log_law = True
scfg.dyn.drag.roughness_length_m = 0.1

sim0 = StormSimulation(scfg); b = sim0.base
sounding_diagnostics(b)
cx, cy = snd.bunkers_storm_motion(b)
base_sr = BaseState(zc=b.zc, theta0=b.theta0, qv0=b.qv0, p0=b.p0, T0=b.T0, rho0=b.rho0,
                    u0=b.u0 - cx, v0=b.v0 - cy)
sim = StormSimulation(scfg, base=base_sr); sim.state.diagnose(sim.cfg)

CACHE = os.path.join(REPO, "outputs", "parent_matured_%d_%d_%d.npz" % (nx, nz, int(T_MAT)))
if not os.path.exists(CACHE):
    log("FATAL: cached matured parent missing (%s)" % CACHE)
    sys.exit(1)
z = np.load(CACHE); xpn = sim.grid.xp
for k in z.files:
    if k != "_t":
        setattr(sim.state, k, xpn.asarray(z[k]))
sim.t = float(z["_t"][0]); sim.state.t = sim.t; sim.state.diagnose(sim.cfg)
log("loaded cached matured parent t=%.0f" % sim.t)

to = sim.grid.backend.to_cpu
uc, vc, _ = rot._centered_velocity(sim.state, sim.grid)
uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); zcp = np.asarray(to(sim.grid.zc))
nbp = nx // 6
k100 = int(np.argmin(np.abs(zcp - 100.0)))
zl = np.abs(np.gradient(vc[:, :, k100], sim.grid.dx, axis=0)
            - np.gradient(uc[:, :, k100], sim.grid.dy, axis=1))
a, bq = np.unravel_index(np.argmax(zl[nbp:-nbp, nbp:-nbp]), zl[nbp:-nbp, nbp:-nbp].shape)
a += nbp; bq += nbp
sx, sy = float(sim.grid.xc[a]), float(sim.grid.yc[bq])
log("low-level circulation @(%.1f,%.1f)km zeta=%.4f" % (sx / 1e3, sy / 1e3, float(zl[a, bq])))


def mkspec(i):
    def build(gg):
        ic = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.xc)) - sx)))
        jc = int(np.argmin(np.abs(np.asarray(gg.backend.to_cpu(gg.yc)) - sy)))
        n = NCX[i]
        i0 = int(np.clip(ic - n // 2, 1, gg.nx - n - 1))
        j0 = int(np.clip(jc - n // 2, 1, gg.ny - n - 1))
        return nst.NestSpec.aligned(gg, i0=i0, j0=j0, ncx=n, ncy=n, refine=3,
                                    nz=NEST_NZ, z_stretch=NEST_ZS, relax_width_m=SPONGE_M)
    return build


def level_report(s, lvl):
    """Where does the low-level vortex sit inside THIS level's box, and how contaminated is it?"""
    g = s.grid
    zeta = np.asarray(g.backend.to_cpu(rot.vertical_vorticity(s.state, g)))
    zc = np.asarray(g.backend.to_cpu(g.zc))
    k = int(np.argmin(np.abs(zc - 100.0)))
    aa = np.abs(zeta[:, :, k])
    n = aa.shape[0]
    idx = np.arange(n)
    db = np.minimum(np.minimum(idx[:, None], n - 1 - idx[:, None]),
                    np.minimum(idx[None, :], n - 1 - idx[None, :]))
    spec = getattr(s, "spec", None)
    band = nst.effective_relax_width(spec, g) if spec is not None else 0
    trust = np.where(db >= band, aa, 0.0)                 # peak OUTSIDE the sponge = trustworthy
    i, j = np.unravel_index(int(np.argmax(trust)), aa.shape)
    gi, gj = np.unravel_index(int(np.argmax(aa)), aa.shape)   # peak anywhere, sponge included
    half = (n - 1) / 2.0
    interior = float(aa[db >= 12].mean()) if (db >= 12).any() else float("nan")
    edge = float(aa[db < 4].mean()) if (db < 4).any() else float("nan")
    return {"level": lvl, "dx_m": float(g.dx), "nx": int(n), "width_km": float(g.dx) * n / 1e3,
            "sponge_cells": int(band), "sponge_m": float(band * g.dx),
            "peak_zeta_trusted": float(trust[i, j]), "peak_zeta_anywhere": float(aa[gi, gj]),
            "peak_edge_dist_m": float(int(db[i, j]) * g.dx),
            "peak_frac_of_halfwidth": float(max(abs(i - half), abs(j - half)) / half),
            "global_peak_in_sponge": bool(int(db[gi, gj]) < band),
            "global_peak_edge_dist_m": float(int(db[gi, gj]) * g.dx),
            "edge_over_interior": (edge / interior) if interior and interior == interior else float("nan")}


series = []
_n = {"i": 0}


def sample(sims, t):
    _n["i"] += 1
    if _n["i"] % SAMPLE_EVERY:
        return
    try:
        recs = [level_report(s, lvl) for lvl, s in enumerate(sims)]
        series.append({"t": float(t), "levels": recs})
        for r in recs:
            log("  t=%5.1f L%d dx=%6.1fm dom=%5.2fkm sponge=%4.0fm | trusted peak|z|=%.4f at "
                "%6.0fm from edge (%.2f of half-width) | edge/interior=%5.2f%s"
                % (t, r["level"], r["dx_m"], r["width_km"], r["sponge_m"],
                   r["peak_zeta_trusted"], r["peak_edge_dist_m"], r["peak_frac_of_halfwidth"],
                   r["edge_over_interior"],
                   "  <-- GLOBAL PEAK IS IN THE SPONGE" if r["global_peak_in_sponge"] else ""))
        json.dump(series, open(os.path.join(OUT, "tracking_diagnosis.json"), "w"), indent=1)
    except Exception as e:
        log("  sample error t=%.1f: %r" % (t, e))


log("cascade DIAGNOSIS: NCX=%s window=%.0fs sponge=%.0fm follow_border=%s follow_filter=%.2f"
    % (NCX, WIN, SPONGE_M, FOLLOW_BORDER, FOLLOW_FILTER))
sims, rep = nst.run_multilevel_nest(sim, [mkspec(i) for i in range(len(NCX))], window=WIN,
                                    restrict_up=True, restrict_momentum=True,
                                    follow_interval=8, follow_field="zeta", follow_frac=0.4,
                                    follow_filter=FOLLOW_FILTER, follow_border=FOLLOW_BORDER,
                                    follow_z_lo=0.0, follow_z_hi=1500.0,
                                    les_boost=1.5, cfl=0.2, sample=sample)

log("")
log("=== VERDICT (the FIRST level whose vortex enters its sponge is the one to fix) ===")
nlev = max(len(s["levels"]) for s in series) if series else 0
for lvl in range(nlev):
    rs = [s["levels"][lvl] for s in series if lvl < len(s["levels"])]
    if not rs:
        continue
    inside = sum(1 for r in rs if r["global_peak_in_sponge"])
    frac = [r["peak_frac_of_halfwidth"] for r in rs]
    eoi = [r["edge_over_interior"] for r in rs if r["edge_over_interior"] == r["edge_over_interior"]]
    log("L%d dx=%6.1fm dom=%5.2fkm sponge=%4.0fm: global peak in sponge %2d/%2d samples | "
        "mean trusted-peak position %.2f of half-width (1.00 = at the wall) | mean edge/interior %.2f"
        % (lvl, rs[0]["dx_m"], rs[0]["width_km"], rs[0]["sponge_m"], inside, len(rs),
           sum(frac) / len(frac), (sum(eoi) / len(eoi)) if eoi else float("nan")))

json.dump({"window_s": WIN, "sponge_m": SPONGE_M, "follow_border": FOLLOW_BORDER,
           "follow_filter": FOLLOW_FILTER, "nest_moves": rep.get("nest", {}).get("nest_moves"),
           "series": series},
          open(os.path.join(OUT, "tracking_diagnosis_summary.json"), "w"), indent=1)
log("wrote %s" % os.path.join(OUT, "tracking_diagnosis.json"))
