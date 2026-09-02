"""Attempt I -- the falsifiable test of the misalignment diagnosis.

A FREELY-EVOLVING idealised supercell (strong CAPE + curved-hodograph shear, warm-bubble trigger,
NO held forcing), integrated long enough (~90 min) for a rear-flank downdraft / occlusion to
develop, tracking over time the two numbers the tilting-efficiency diagnosis identified:

  * low-level tilting ALIGNMENT  cos(theta) = sum(omega_h . grad_h w) / sum(|omega_h||grad_h w|)
  * low-level STREAMWISE FRACTION of the horizontal vorticity.

Hypothesis: if this model can make a tornado, these must RISE as the storm occludes -- before V_rot
does.  If they stay ~0.1, the model's RFD/occlusion is not reorienting the vortex lines (a deeper
structural gap).  Resolves the storm-scale RFD/occlusion at dx~600 m (no nest needed for the
alignment question)."""
import os, sys, time, json
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
import numpy as np, warnings; warnings.filterwarnings("ignore")

OUT = os.path.join(REPO, "outputs", "supercell_alignment"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "progress.log"); open(LOG, "w").close()
_t0 = time.time()
def log(m):
    open(LOG, "a", encoding="utf-8").write("[%6.1fs] %s\n" % (time.time() - _t0, m)); print(m, flush=True)

from meteorological_flow.base_state import BaseState, sounding_diagnostics
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, soundings as snd, vorticity_budget as vb, coldpool as cp

def low_level_metrics(sim, z_top=800.0):
    g = sim.grid; to = g.backend.to_cpu
    uc, vc, wc = rot._centered_velocity(sim.state, g)              # device arrays
    e = vb.tilting_efficiency(uc, vc, wc, g)                        # computed on-device
    xi, eta = vb.horizontal_vorticity(uc, vc, wc, g)
    sw, cw = vb.streamwise_crosswise(uc, vc, xi, eta, g)
    zc = np.asarray(to(g.zc)); nb = g.nx // 6; kmask = zc < z_top
    cpu = lambda a: np.asarray(to(a))[nb:-nb, nb:-nb][:, :, kmask]  # host + interior + low layer
    omh = cpu(e["omega_h"]); gw = cpu(e["grad_h_w"]); tilt = cpu(e["tilting"])
    align = float(tilt.sum() / (np.abs(omh * gw).sum() + 1e-20))    # activity-weighted cos(theta)
    s = np.abs(cpu(sw)).sum(); c = np.abs(cpu(cw)).sum()
    swf = float(s / (s + c + 1e-20))
    ucn = np.asarray(to(uc)); vcn = np.asarray(to(vc))
    zeta = np.gradient(vcn, g.dx, axis=0) - np.gradient(ucn, g.dy, axis=1)
    zlow = float(np.abs(zeta[nb:-nb, nb:-nb, zc < 500.0]).max())
    return {"align": align, "streamwise_frac": swf, "w_max": float(np.asarray(to(wc)).max()),
            "low_zeta": zlow, "midmeso": rot.rotation_report(sim.state, g)["midlevel_mesocyclone"]}

DEV = os.environ.get("DEV", "gpu")
log("=== Attempt I: freely-evolving supercell, tracking tilting alignment + streamwise (%s) ===" % DEV)
# idealised strong-supercell environment (curved hodograph -> SRH), freely evolving
nx = 120; nz = 48; Lx = 72000.0
scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=Lx, Ly=Lx, Lz=15000.0,
                          duration=1.0, dt_max=3.0, drag=True, z_stretch=1.05, C_s=0.20,
                          hodograph_kind="quarter_circle", U_max=30.0, device=DEV)
scfg.sim.physics.bubble_dtheta = 5.0
sim0 = StormSimulation(scfg)                       # build once to get the WK base + Bunkers motion
b = sim0.base
d = sounding_diagnostics(b)
cx, cy = snd.bunkers_storm_motion(b)
log("env: CAPE=%.0f CIN=%.0f shear06=%.1f SRH03=%.0f | Bunkers C=(%.1f,%.1f) dx=%.0fm" %
    (d["CAPE_J_kg"], d["CIN_J_kg"], d["shear_0_6km_m_s"], snd.storm_relative_helicity(b), cx, cy, sim0.grid.dx))
# storm-relative base (Galilean shift) so the freely-evolving storm stays framed for the long run
base_sr = BaseState(zc=b.zc, theta0=b.theta0, qv0=b.qv0, p0=b.p0, T0=b.T0, rho0=b.rho0,
                    u0=b.u0 - cx, v0=b.v0 - cy)
sim = StormSimulation(scfg, base=base_sr); sim.state.diagnose(sim.cfg)

T_END = 5400.0; SAMPLE = 300.0; next_s = 300.0
hist = []
while sim.t < T_END:
    dt = float(sim._dt()); sim._step(dt); sim.step += 1; sim.t = float(sim.state.t)
    w = np.asarray(sim.grid.backend.to_cpu(sim.state.w))
    if not np.isfinite(w).all():
        log("  instability at t=%.0f -> stop" % sim.t); break
    if sim.t >= next_s:
        sim.state.diagnose(sim.cfg)
        m = low_level_metrics(sim); m["t"] = round(sim.t)
        hist.append(m)
        log("  t=%4.0f w_max=%5.1f midmeso=%.2e low_zeta=%.2e | ALIGN=%+.3f streamwise=%.2f" %
            (sim.t, m["w_max"], m["midmeso"], m["low_zeta"], m["align"], m["streamwise_frac"]))
        next_s += SAMPLE

# trend of the two decisive numbers over the mature phase (second half)
mature = [h for h in hist if h["t"] >= T_END * 0.4]
al = [h["align"] for h in mature]; sf = [h["streamwise_frac"] for h in mature]
trend_align = (al[-1] - al[0]) if len(al) > 1 else 0.0
trend_sf = (sf[-1] - sf[0]) if len(sf) > 1 else 0.0
cprep = cp.coldpool_report(sim.state, sim.grid)
summary = {"attempt": "I freely-evolving supercell", "env": {k: float(d[k]) for k in ("CAPE_J_kg", "CIN_J_kg", "shear_0_6km_m_s") if d[k] is not None},
           "SRH03": float(snd.storm_relative_helicity(sim.base)), "dx_m": float(sim.grid.dx),
           "history": hist, "align_trend_mature": trend_align, "streamwise_trend_mature": trend_sf,
           "final_align": al[-1] if al else 0.0, "final_streamwise": sf[-1] if sf else 0.0,
           "cold_pool": cprep, "wall_clock_s": time.time() - _t0}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
verdict = ("ALIGNMENT RISES as it occludes -> model CAN build the geometry" if trend_align > 0.05
           else "alignment stays low -> RFD/occlusion does NOT reorient the vortex lines (structural gap)")
log("=== DONE I: final align=%+.3f (trend %+.3f) streamwise=%.2f (trend %+.3f) | %s (%.0fmin) ===" %
    (summary["final_align"], trend_align, summary["final_streamwise"], trend_sf, verdict, (time.time() - _t0) / 60))
