"""Figures for Attempt L (intensity): the persistence traces, the surface vorticity structure, the
surface pressure field (measurable on a nest for the first time), and the V_rot profile against
Attempts J/K and the observed KTLX couplet.

The question these must answer honestly: does the 22 m mesh turn Attempt K's linear shear BAND into
a compact axisymmetric CORE, and does it hold?
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
OUT = os.path.join(REPO, "outputs", "tornado_intensity_L")
MEDIA = os.path.join(REPO, "docs", "media", "storm"); os.makedirs(MEDIA, exist_ok=True)

f = np.load(os.path.join(OUT, "fields_L.npz"))
summary = json.load(open(os.path.join(OUT, "summary.json")))
series = summary.get("series") or json.load(open(os.path.join(OUT, "series.json")))

xc, yc, zc = f["xc"], f["yc"], f["zc"]
zeta, u, v, w, p = f["zeta"], f["u"], f["v"], f["w"], f["p"]
dx = float(f["dx"])
OBS_V = 26.0

fig = plt.figure(figsize=(16.5, 9.2), facecolor="white")
fig.suptitle("Attempt L — intensity: 3-level cascade to %.0f m, %.0f s window "
             "(Moore-2013-like environment)" % (dx, summary.get("window_s", 0)),
             fontsize=14, fontweight="bold")

# ---- (a) persistence traces -------------------------------------------------------------
ax = fig.add_subplot(2, 3, 1)
t = [s["t"] for s in series]
ax.plot(t, [s["v_rot_sfc"] for s in series], "-", color="#c1121f", lw=1.8, label="V$_{rot}$ surface")
ax.plot(t, [s["v_theta"] or np.nan for s in series], "-", color="#264653", lw=1.2, alpha=.8,
        label=r"V$_\theta$ (vortex report)")
ax.axhline(OBS_V, ls="--", color="k", lw=1.2)
ax.text(t[0] if t else 0, OBS_V * 1.02, "observed KTLX 26 m/s", fontsize=8)
ax.axhline(7.76, ls=":", color="#6a994e", lw=1.4)
ax.text(t[0] if t else 0, 8.1, "Attempt K (67 m)", fontsize=8, color="#6a994e")
ax.set_xlabel("t within window (s)"); ax.set_ylabel("m s$^{-1}$")
ax.set_title("(a) Persistence — is the vortex lasting or transient?", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=.3)

# ---- (b) pressure deficit trace (NEW: impossible before the nest-pressure fix) -----------
ax = fig.add_subplot(2, 3, 2)
dP = [s["dP"] if s["dP"] is not None and s["dP"] == s["dP"] else np.nan for s in series]
ax.plot(t, dP, "-", color="#5a189a", lw=1.8)
ax.axhline(0, color="k", lw=.8)
ax.set_xlabel("t within window (s)"); ax.set_ylabel(r"$\Delta p$ (Pa)")
ax.set_title(r"(b) Core pressure deficit — first measurable on a nest", fontsize=10)
ax.grid(alpha=.3)

# ---- (c) surface zeta with storm-relative streamlines ------------------------------------
ax = fig.add_subplot(2, 3, 3)
k0 = 0
Z = zeta[:, :, k0]
X, Y = np.meshgrid((xc - xc.mean()) / 1e3, (yc - yc.mean()) / 1e3, indexing="ij")
lim = float(np.abs(Z).max())
m = ax.pcolormesh(X, Y, Z, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto")
plt.colorbar(m, ax=ax, label=r"$\zeta$ (s$^{-1}$)")
us = u[:, :, k0] - u[:, :, k0].mean(); vs = v[:, :, k0] - v[:, :, k0].mean()
ax.streamplot(X.T, Y.T, us.T, vs.T, color="k", density=1.1, linewidth=.5, arrowsize=.6)
ax.set_title(r"(c) Surface $\zeta$ + storm-relative flow, z=%.1f m" % zc[k0], fontsize=10)
ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)"); ax.set_aspect("equal")

# ---- (d) surface perturbation pressure ---------------------------------------------------
ax = fig.add_subplot(2, 3, 4)
P = p[:, :, k0] - np.median(p[:, :, k0])
pl = float(np.percentile(np.abs(P), 99))
m = ax.pcolormesh(X, Y, P, cmap="viridis", vmin=-pl, vmax=pl, shading="auto")
plt.colorbar(m, ax=ax, label="p' (Pa)")
ax.set_title("(d) Surface perturbation pressure", fontsize=10)
ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)"); ax.set_aspect("equal")

# ---- (e) vertical cross-section of zeta through the peak ---------------------------------
ax = fig.add_subplot(2, 3, 5)
nb = max(2, int(0.2 * Z.shape[0]))
sub = np.abs(Z[nb:-nb, nb:-nb])
ii, jj = np.unravel_index(np.argmax(sub), sub.shape); ii += nb; jj += nb
CS = zeta[:, jj, :]
lim2 = float(np.abs(CS).max())
m = ax.pcolormesh((xc - xc.mean()) / 1e3, zc, CS.T, cmap="RdBu_r", vmin=-lim2, vmax=lim2,
                  shading="auto")
plt.colorbar(m, ax=ax, label=r"$\zeta$ (s$^{-1}$)")
ax.set_ylim(0, min(2000, float(zc.max())))
ax.set_xlabel("x (km)"); ax.set_ylabel("z (m)")
ax.set_title(r"(e) Vertical section through peak $|\zeta|$", fontsize=10)

# ---- (f) V_rot profile: J vs K vs L vs observed -------------------------------------------
ax = fig.add_subplot(2, 3, 6)
J_z = [40, 208, 522, 968, 1488]; J_v = [2.72, 4.05, 5.30, 6.71, 7.67]
K_z = [5.11, 52.8, 203, 522, 968, 1462]; K_v = [7.76, 7.54, 6.94, 6.40, 5.56, 6.71]
ax.plot(J_v, J_z, "o-", color="#8d99ae", label="J — elevated (22 m, fixed nest)")
ax.plot(K_v, K_z, "s-", color="#6a994e", label="K — surface-intensified (67 m)")
prof = (summary.get("surface") or {}).get("profile") or []
if prof:
    ax.plot([q["v_rot_m_s"] for q in prof], [q["z_m"] for q in prof], "^-", color="#c1121f",
            lw=2.2, label="L — %.0f m, 3 levels" % dx)
ax.axvline(OBS_V, ls="--", color="k", lw=1.2)
ax.text(OBS_V * 0.63, 1200, "observed\nKTLX 26 m/s", fontsize=8)
ax.set_xlabel("V$_{rot}$ (m s$^{-1}$)"); ax.set_ylabel("z (m)")
ax.set_title("(f) Rotation profile vs the observed couplet", fontsize=10)
ax.legend(fontsize=7.5); ax.grid(alpha=.3)

pers = summary.get("persistence") or {}
fig.text(0.5, 0.008,
         "peak surface V$_{rot}$ = %.2f m/s (K: 7.76, observed: 26) · connected %.0f%% of %d samples "
         "· min Δp = %.1f Pa · class: %s"
         % (pers.get("peak_v_rot_sfc", float("nan")),
            100 * pers.get("connected_fraction", 0), pers.get("samples", 0),
            pers.get("min_dP_Pa", float("nan")), summary.get("classification", "n/a")),
         ha="center", fontsize=9.5)
fig.tight_layout(rect=[0, 0.022, 1, 0.96])
dst = os.path.join(MEDIA, "tornado_intensity_L.png")
fig.savefig(dst, dpi=140, facecolor="white")
print("wrote", dst)
