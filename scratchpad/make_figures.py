"""Figures for the README + x.com: the honest tornadogenesis journey, and a hero of the
resolved low-level mesocyclone at 28 m from real data."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
MED = os.path.join(REPO, "docs", "media", "storm"); os.makedirs(MED, exist_ok=True)

f = np.load(os.path.join(REPO, "outputs", "moore_cascade", "fields_cascade.npz"))
u, v, w, zeta = f["u"], f["v"], f["w"], f["zeta"]
xc, yc, zc, dx = f["xc"], f["yc"], f["zc"], float(f["dx"])
summ = json.load(open(os.path.join(REPO, "outputs", "moore_cascade", "summary.json")))

# ---- low-level slice at ~500 m: the full nest INTERIOR (crop the sponge border) ------
kl = int(np.argmin(np.abs(zc - 500)))
nb = zeta.shape[0] // 6
Z = zeta[:, :, kl]
sub = np.abs(Z[nb:-nb, nb:-nb]); i0, j0 = np.unravel_index(np.argmax(sub), sub.shape)
i0 += nb; j0 += nb                                  # interior vorticity max (for centring)
sl = (slice(nb, -nb), slice(nb, -nb))               # show the whole interior, honestly
Zc = Z[sl].T
uu = (u[:, :, kl] - u[:, :, kl].mean())[sl].T      # storm-relative -> reveals rotation
vv = (v[:, :, kl] - v[:, :, kl].mean())[sl].T
xx = (xc[sl[0]] - xc[i0]) / 1e3; yy = (yc[sl[1]] - yc[j0]) / 1e3
zmax = float(np.percentile(np.abs(Zc), 99.5))       # robust colour scale
attempts = ["A\nbubble\n46 m", "B\nbubble\n521 m", "C\nforced\n250 m", "D\ncascade\n28 m"]
vrots = [3.2, 2.3, 4.9, 6.0]; obs = 26.0

# =====================================================================================
# FIGURE 1 -- README: the honest journey (3 panels)
# =====================================================================================
plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#888"})
fig = plt.figure(figsize=(15, 5.0), facecolor="white")
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.8, 1.15], wspace=0.32)

# (a) V_rot progression toward the observed tornado
axa = fig.add_subplot(gs[0])
cols = ["#9db4c0", "#c98a8a", "#5b8c5a", "#2f6690"]
bars = axa.bar(range(4), vrots, color=cols, width=0.7, zorder=3)
axa.axhline(obs, ls="--", lw=2, color="#d1495b", zorder=2)
axa.text(3.45, obs - 1.4, "observed KTLX\ntornado  V_rot 26 m/s", color="#d1495b",
         ha="right", va="top", fontsize=10, fontweight="bold")
for b, val in zip(bars, vrots):
    axa.text(b.get_x() + b.get_width() / 2, val + 0.4, "%.1f" % val, ha="center", fontsize=10, fontweight="bold")
axa.set_xticks(range(4)); axa.set_xticklabels(attempts, fontsize=9)
axa.set_ylabel("low-level rotational velocity  V_rot  [m/s]")
axa.set_ylim(0, 28); axa.set_title("(a)  from real data toward the tornado", fontsize=11, loc="left", fontweight="bold")
axa.grid(axis="y", alpha=0.3, zorder=0)

# (b) vertical structure -- surface-connected low-level mesocyclone
axb = fig.add_subplot(gs[1])
zp = np.array(summ["zeta_profile_z_val"], float)
axb.plot(zp[:, 1] * 1e3, zp[:, 0] / 1e3, "-o", color="#2f6690", lw=2, ms=4)
axb.fill_betweenx(zp[:, 0] / 1e3, 0, zp[:, 1] * 1e3, color="#2f6690", alpha=0.15)
axb.set_ylim(0, 6); axb.set_xlim(0, None)
axb.set_xlabel("peak |vorticity|  [10$^{-3}$ s$^{-1}$]"); axb.set_ylabel("height  [km]")
axb.set_title("(b)  a low-level meso,\n      reaching the ground", fontsize=11, loc="left", fontweight="bold")
axb.annotate("peak 0.5-1 km", xy=(11.7, 0.5), xytext=(6, 3.2), fontsize=9,
             arrowprops=dict(arrowstyle="->", color="#444"))
axb.grid(alpha=0.3)

# (c) the resolved vortex at 28 m
axc = fig.add_subplot(gs[2])
norm = TwoSlopeNorm(vmin=-zmax, vcenter=0, vmax=zmax)
pcm = axc.pcolormesh(xx, yy, Zc, cmap="RdBu_r", norm=norm, shading="auto")
sp = axc.streamplot(xx, yy, uu, vv, color="k", density=1.1, linewidth=0.6, arrowsize=0.7)
axc.set_aspect("equal"); axc.set_xlabel("x  [km]"); axc.set_ylabel("y  [km]")
axc.set_title("(c)  resolved low-level vortex, 28 m", fontsize=11, loc="left", fontweight="bold")
cb = fig.colorbar(pcm, ax=axc, fraction=0.046, pad=0.03); cb.set_label("vertical vorticity  [s$^{-1}$]")

fig.suptitle("Tornadogenesis from real data (Moore, OK, 20 May 2013) — the honest edge: a sustained "
             "supercell + AMR cascade builds a surface-connected low-level mesocyclone at ~23% of the "
             "observed tornado intensity",
             fontsize=11.5, y=1.02, x=0.5)
fig.savefig(os.path.join(MED, "tornadogenesis_journey.png"), dpi=140, bbox_inches="tight", facecolor="white")
print("saved docs/media/storm/tornadogenesis_journey.png")

# =====================================================================================
# FIGURE 2 -- x.com hero: the resolved low-level mesocyclone (dark, striking)
# =====================================================================================
plt.rcParams.update({"font.size": 13})
figh = plt.figure(figsize=(9, 9), facecolor="#0d1117")
axh = figh.add_axes([0.02, 0.065, 0.96, 0.845]); axh.set_facecolor("#0d1117")
pcm = axh.pcolormesh(xx, yy, Zc, cmap="RdBu_r", norm=norm, shading="gouraud")
axh.streamplot(xx, yy, uu, vv, color="white", density=1.5, linewidth=0.8, arrowsize=0.9)
axh.set_aspect("equal"); axh.set_xticks([]); axh.set_yticks([])
for s in axh.spines.values():
    s.set_visible(False)
axh.text(0.5, 1.075, "A low-level mesocyclone, resolved at 28 m",
         transform=axh.transAxes, ha="center", va="bottom",
         color="white", fontsize=18, fontweight="bold")
axh.text(0.5, 1.028, "grown from the REAL Moore 2013 supercell environment — anelastic LES + storm-relative AMR cascade, GPU",
         transform=axh.transAxes, ha="center", va="bottom", color="#9db4c0", fontsize=10)
axh.text(0.5, -0.035, "vertical vorticity (red/blue = opposite spin) · storm-relative streamlines · the honest edge of tornadogenesis, not yet a tornado",
         transform=axh.transAxes, ha="center", va="top", color="#8b98a5", fontsize=9)
figh.savefig(os.path.join(MED, "meso_hero_28m.png"), dpi=150, bbox_inches="tight", facecolor="#0d1117")
print("saved docs/media/storm/meso_hero_28m.png")
