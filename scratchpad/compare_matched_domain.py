"""Verdict for the matched-domain resolution test: 67 m vs 22 m over the SAME 2.0 km domain.

Prints the controlled comparison and writes a two-panel figure.  The claim under test is narrow:
holding the finest domain, the storm, the instant, the centring, the measurement radius and the
interior margin fixed, does dx alone intensify the surface vortex?

Reported honestly either way -- including the case where refinement changes nothing, or makes the
measured rotation worse.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
OUT = os.path.join(REPO, "outputs", "matched_domain")
MEDIA = os.path.join(REPO, "docs", "media", "storm"); os.makedirs(MEDIA, exist_ok=True)

S = {}
for mode in ("coarse", "fine"):
    fp = os.path.join(OUT, "summary_%s.json" % mode)
    if not os.path.exists(fp):
        raise SystemExit("missing %s -- run tornado_matched_domain_gpu.py MODE=%s" % (fp, mode))
    S[mode] = json.load(open(fp))

print("=" * 96)
print("MATCHED-DOMAIN RESOLUTION TEST -- same 2.0 km finest domain, same storm/instant, "
      "only dx differs")
print("=" * 96)
hdr = ("%-8s %8s %7s %8s %9s %9s %8s %8s %8s %9s"
       % ("mode", "dx(m)", "cells", "dom(km)", "V_sfc_pk", "V_sfc_mn", "zeta_pk", "Vth_pk",
          "ratio", "conn%"))
print(hdr); print("-" * len(hdr))
for mode in ("coarse", "fine"):
    s = S[mode]; a = s.get("aggregate", {})
    print("%-8s %8.1f %7d %8.2f %9.2f %9.2f %8.3f %8.2f %8.2f %9.0f"
          % (mode, s["finest_dx_m"], s["finest_nx"], s["finest_width_km"],
             a.get("peak_v_rot_sfc", float("nan")), a.get("mean_v_rot_sfc", float("nan")),
             a.get("peak_zeta_sfc", float("nan")), a.get("peak_v_theta", float("nan")),
             s["surface"]["surface_aloft_ratio"], 100 * a.get("connected_fraction", 0)))

c, f = S["coarse"]["aggregate"], S["fine"]["aggregate"]
print()
def ratio(k):
    a, b = f.get(k), c.get(k)
    return (a / b) if (a and b) else float("nan")
print("dx %.1f m -> %.1f m (x%.1f refinement) at FIXED 2.0 km domain:"
      % (S["coarse"]["finest_dx_m"], S["fine"]["finest_dx_m"],
         S["coarse"]["finest_dx_m"] / S["fine"]["finest_dx_m"]))
for k, lbl in (("peak_v_rot_sfc", "peak surface V_rot"), ("mean_v_rot_sfc", "mean surface V_rot"),
               ("peak_zeta_sfc", "peak surface |zeta|"), ("peak_v_theta", "peak v_theta")):
    print("   %-22s x%.2f   (%.3f -> %.3f)" % (lbl, ratio(k), c.get(k, float('nan')),
                                               f.get(k, float('nan'))))

# The pressure sanity check: a nest's projection pressure carries a boundary-imposed component
# that grows like 1/dt, so a deficit far exceeding -rho*v_theta^2 is an artifact, not a tornado.
print("\npressure-deficit sanity (|dP| should be within ~a few x the cyclostrophic scale):")
for mode in ("coarse", "fine"):
    ser = S[mode]["series"]
    vals = [(r["dP"], r["dP_cyclostrophic"]) for r in ser
            if r["dP"] is not None and r["dP"] == r["dP"] and r["dP_cyclostrophic"]]
    if vals:
        rr = [abs(a) / max(abs(b), 1e-9) for a, b in vals]
        print("   %-8s median |dP|/|rho v^2| = %6.2f   (n=%d, min %.2f max %.2f)"
              % (mode, float(np.median(rr)), len(rr), min(rr), max(rr)))
    else:
        print("   %-8s no finite dP samples" % mode)

OBS_V, OBS_Z = 26.0, 0.205
print("\nagainst the observed Moore-2013 KTLX couplet (V_rot 26 m/s, zeta 0.205 /s):")
for mode in ("coarse", "fine"):
    a = S[mode]["aggregate"]
    print("   %-8s V_sfc %.2f = %3.0f%% of observed | zeta %.3f = %3.0f%% of observed | class %s"
          % (mode, a.get("peak_v_rot_sfc", float("nan")),
             100 * a.get("peak_v_rot_sfc", 0) / OBS_V, a.get("peak_zeta_sfc", float("nan")),
             100 * a.get("peak_zeta_sfc", 0) / OBS_Z, S[mode]["classification"]))

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), facecolor="white")
fig.suptitle("Matched-domain resolution test — 67 m vs 22 m over the same 2.0 km domain",
             fontsize=13, fontweight="bold")
cols = {"coarse": "#6a994e", "fine": "#c1121f"}
for mode in ("coarse", "fine"):
    ser = S[mode]["series"]; t = [r["t"] for r in ser]
    lbl = "dx = %.1f m (%d cells)" % (S[mode]["finest_dx_m"], S[mode]["finest_nx"])
    axes[0].plot(t, [r["v_rot_sfc"] for r in ser], "-", color=cols[mode], lw=1.8, label=lbl)
    axes[1].plot(t, [r["zeta_sfc"] for r in ser], "-", color=cols[mode], lw=1.8, label=lbl)
    axes[2].plot(t, [r["ratio"] for r in ser], "-", color=cols[mode], lw=1.8, label=lbl)
axes[0].axhline(OBS_V, ls="--", color="k", lw=1.1)
axes[0].text(0, OBS_V * 1.02, "observed 26 m/s", fontsize=8)
axes[0].set_ylabel("V$_{rot}$ surface (m s$^{-1}$)"); axes[0].set_title("surface rotation", fontsize=10)
axes[1].axhline(OBS_Z, ls="--", color="k", lw=1.1)
axes[1].text(0, OBS_Z * 1.02, r"observed 0.205 s$^{-1}$", fontsize=8)
axes[1].set_ylabel(r"$|\zeta|$ surface (s$^{-1}$)"); axes[1].set_title("surface vorticity", fontsize=10)
axes[2].axhline(0.8, ls="--", color="k", lw=1.1)
axes[2].text(0, 0.82, "connected threshold", fontsize=8)
axes[2].set_ylabel("surface / aloft ratio"); axes[2].set_title("surface connection", fontsize=10)
for ax in axes:
    ax.set_xlabel("t within window (s)"); ax.grid(alpha=.3); ax.legend(fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.93])
dst = os.path.join(MEDIA, "matched_domain_resolution.png")
fig.savefig(dst, dpi=140, facecolor="white")
print("\nwrote", dst)
