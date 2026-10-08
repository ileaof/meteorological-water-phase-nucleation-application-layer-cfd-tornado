"""Figures for a real non-US convective environment (Brazil, Centro-Oeste).

These are **diagnostics of a real GFS profile, not a simulation** -- no CFD is integrated
here.  The question they answer is whether the measured environment supports the
supercell -> mesocyclone -> occlusion chain the A-L study models, and the honest answer at
these points is no, for one dominant reason: the deep-layer shear is an order of magnitude
short of the supercell threshold.

Two design choices carry the argument, and both are deliberate:

* **Small multiples with a SHARED axis range**, set by the project's own idealised tornadic
  hodograph (CAPE ~2225, shear ~41 m/s, SRH ~648).  Auto-scaling each panel would draw a
  1.4 m/s hodograph the same size as a 20 m/s one -- the figure would then show structure
  and hide the magnitude, which is the whole finding.  The reference panel is what makes
  the other four legible as *small*.
* **The ill-conditioning flag is written as text**, never as colour alone: below ~2 m/s of
  layer shear the critical angle is noise, and a panel must say so rather than imply
  precision (see ``soundings.WELL_POSED_SHEAR_MS``).

Read the hodograph dots carefully: they are **model-grid cells (dz = 100 m), not GFS
levels**.  GFS reports every 25 hPa near the surface (~200-500 m apart), so the straight
runs between kinks are interpolation, not observed structure -- the kinks are where the
data actually is.

Height is encoded sequentially with **viridis**, matching ``plotting.plot_hodograph`` so
this figure reads alongside the project's existing ones.  That is a deliberate deviation
from a single-hue default: viridis is monotonic in lightness and CVD-safe (it is not the
rainbow/jet the rule targets) and it ships with a colourbar.  Series colours are the
validated categorical slots 1-2 and the ``critical`` status step for threshold rules --
hand-picked hexes failed the lightness-band and chroma-floor checks and were replaced.

Run::

    python examples/real_case_centro_oeste.py                     # live GFS
    python examples/real_case_centro_oeste.py --when 2026-10-08T18:00:00Z
"""
from __future__ import annotations

import argparse
import os
import sys
import warnings

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from meteorological_flow.grid import Grid            # noqa: E402
from storm_dynamics import soundings as snd          # noqa: E402

#: Points sampled, with latitude so the hemisphere is inferred rather than assumed.
POINTS = [
    ("Dourados, MS", -22.22, -54.81),
    ("Campo Grande, MS", -20.47, -54.67),
    ("Cuiaba, MT", -15.65, -56.12),
    ("Goiania, GO", -16.63, -49.22),
]

#: Literature thresholds the measurement is scored against, with what they gate.
SUPERCELL_SHEAR_MS = 18.0        # 0-6 km bulk shear below this does not organise supercells
TORNADIC_LCL_M = 1200.0          # tornadic environments want a LOW cloud base

#: CAPE is grid-sensitive by a factor of ~7; dz <= 100 m is the converged regime.
GRID = Grid(nx=4, ny=4, nz=160, Lx=16000, Ly=16000, Lz=16000)


def collect(when):
    """Download each profile and compute the thermodynamics + low-level geometry."""
    from atmospheric_data.sources import gfs, sounding
    rows = []
    for name, lat, lon in POINTS:
        prof = gfs.download_profile(lat, lon, when=when)
        base = sounding.profile_to_basestate(GRID, prof)
        diag = sounding.profile_diagnostics(GRID, prof)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")          # the flag is read from the report
            geom = snd.low_level_geometry_report(base, z_top=500.0, latitude_deg=lat)
        rows.append({"name": name, "lat": lat, "lon": lon, "base": base,
                     "diag": diag, "geom": geom, "gfs": prof["gfs_diagnostics"]})
    return rows


def _reference():
    """The project's own idealised tornadic hodograph -- the scale bar for the others."""
    base = snd.build_sounding(GRID)
    geom = snd.low_level_geometry_report(base, z_top=500.0)      # northern by construction
    return {"name": "idealised tornadic\n(this project's reference)", "lat": 35.2,
            "base": base, "geom": geom, "diag": None, "gfs": None}


def _panel(ax, row, lim, reference=False, annotate=True, center=(0.0, 0.0)):
    base, geom = row["base"], row["geom"]
    z = np.asarray(base.zc); u = np.asarray(base.u0); v = np.asarray(base.v0)
    sel = z <= 6000.0
    ax.plot(u[sel], v[sel], "-", color="0.65", lw=1.2, zorder=2)
    sc = ax.scatter(u[sel], v[sel], c=z[sel] / 1000.0, cmap="viridis",
                    s=14, zorder=3, vmin=0, vmax=6)
    cx, cy = geom["storm_motion_ms"]
    mover = "left" if geom["hemisphere"] == "south" else "right"
    ax.plot(cx, cy, "*", color="#d03b3b", ms=13, zorder=5,
            label="storm motion (%s-mover)" % mover)

    # the 0-500 m shear vector: the quantity that decides the verdict
    u_lo, u_hi = np.interp([0.0, 500.0], z, u)
    v_lo, v_hi = np.interp([0.0, 500.0], z, v)
    ax.annotate("", xy=(u_hi, v_hi), xytext=(u_lo, v_lo), zorder=4,
                arrowprops=dict(arrowstyle="-|>", lw=1.8, color="#2a78d6",
                                shrinkA=0, shrinkB=0))

    ax.axhline(0, color="0.9", lw=0.8, zorder=1)
    ax.axvline(0, color="0.9", lw=0.8, zorder=1)
    ax.set_xlim(center[0] - lim, center[0] + lim)
    ax.set_ylim(center[1] - lim, center[1] + lim)
    ax.set_aspect("equal", "box")
    if annotate:
        ax.set_title(row["name"], fontsize=9.5, pad=6)

    txt = ("shear$_{0-6}$ %.1f m/s\nshear$_{0-500}$ %.2f m/s"
           % (geom["shear_0_6km_m_s"], geom["shear_layer_m_s"]))
    if geom["geometry_is_ill_conditioned"]:
        txt += "\ncrit. angle: ILL-COND."      # stated, not implied by colour
    else:
        txt += "\ncrit. angle %.0f$^\\circ$ (f=%.2f)" % (
            geom["critical_angle_deg"], geom["streamwise_fraction"])
    if annotate:
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left",
                fontsize=8, color="0.15",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.8", lw=0.6, alpha=0.9))
    if reference:
        for spine in ax.spines.values():
            spine.set_edgecolor("#2a78d6")
            spine.set_linewidth(1.6)
    return sc


def _native_box(row):
    """Tight square box around this hodograph AND its storm motion -- ``(xc, yc, half)``.

    Deliberately **not** centred on the origin: these hodographs sit well off (0, 0), so a
    symmetric range spends most of the axis on empty space and the zoom stops at ~x2.
    Boxing the data itself is what makes the shape -- and the 0-500 m shear arrow --
    actually visible.  The storm motion stays inside because the critical angle is the
    angle to it; drop it and the panel would show a shape with nothing to measure against.
    """
    base = row["base"]
    z = np.asarray(base.zc); sel = z <= 6000.0
    u = np.asarray(base.u0)[sel]; v = np.asarray(base.v0)[sel]
    cx, cy = row["geom"]["storm_motion_ms"]
    xs = np.r_[u, cx]; ys = np.r_[v, cy]
    xc = 0.5 * (np.nanmax(xs) + np.nanmin(xs))
    yc = 0.5 * (np.nanmax(ys) + np.nanmin(ys))
    half = 0.5 * max(np.nanmax(xs) - np.nanmin(xs), np.nanmax(ys) - np.nanmin(ys))
    return float(xc), float(yc), float(max(1.0, half * 1.18))


def figure_hodographs(rows, outdir):
    """Two rows, because one cannot carry both facts honestly.

    The shared-range row answers "how much rotation is on offer": the four points are
    specks beside the tornadic reference, which IS the finding.  But at that range their
    structure collapses into a blob and the 0-500 m shear arrow is sub-pixel, so a
    single-row figure would advertise a vector it never actually shows.  The second row
    redraws each panel at its OWN range with the zoom factor printed, so shape becomes
    readable without letting the reader mistake it for magnitude.
    """
    panels = rows + [_reference()]
    lim = max(25.0, max(np.nanmax(np.abs(np.asarray(r["base"].u0))) for r in panels) * 1.05)
    n = len(panels)
    fig, axes = plt.subplots(2, n, figsize=(3.05 * n, 7.3), constrained_layout=True)
    sc = None
    for j, row in enumerate(panels):
        is_ref = row["diag"] is None
        sc = _panel(axes[0, j], row, lim, reference=is_ref, annotate=True)
        xc, yc, half = _native_box(row)
        _panel(axes[1, j], row, half, reference=is_ref, annotate=False, center=(xc, yc))
        axes[1, j].set_title("zoom ×%.0f   (box %.0f m/s wide)" % (lim / half, 2 * half),
                             fontsize=8.5, pad=4, color="0.35")
        axes[1, j].set_xlabel("u [m/s]")
    axes[0, 0].set_ylabel("v [m/s] — shared range")
    axes[1, 0].set_ylabel("v [m/s] — native range")
    fig.colorbar(sc, ax=axes, label="height [km]", shrink=0.6, pad=0.01)
    fig.suptitle("Centro-Oeste hodographs vs a tornadic reference\n"
                 "TOP one shared range — magnitude.   BOTTOM each at its own range "
                 "— shape.   arrow = 0–500 m shear · star = favoured mover",
                 fontsize=11)
    path = os.path.join(outdir, "centro_oeste_hodographs.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def figure_verdict(rows, outdir):
    """Measured shear and cloud base against the thresholds each one gates."""
    names = [r["name"] for r in rows]
    y = np.arange(len(rows))
    shear = [r["geom"]["shear_0_6km_m_s"] for r in rows]
    lcl = [r["diag"]["LCL_m"] for r in rows]
    cape = [r["diag"]["CAPE_J_kg"] for r in rows]
    cape_gfs = [r["gfs"]["CAPE_J_kg"] for r in rows]

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.5), constrained_layout=True)

    ax = axes[0]
    ax.barh(y, shear, height=0.55, color="#2a78d6")
    ax.axvline(SUPERCELL_SHEAR_MS, color="#d03b3b", lw=1.8)
    # pinned to the axes fraction, not a data row: in data coords this label fell below
    # the last bar and was clipped by the axes
    ax.text(SUPERCELL_SHEAR_MS, 0.03, "  supercell threshold ~18 m/s",
            transform=ax.get_xaxis_transform(), color="#d03b3b", fontsize=8.5,
            va="bottom", ha="left", zorder=5,
            bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85))
    ax.set_xlim(0, max(SUPERCELL_SHEAR_MS * 1.35, max(shear) * 1.2))
    ax.set_title("0–6 km bulk shear [m/s]\nmeasured far below the threshold", fontsize=9.5)
    for i, val in enumerate(shear):
        ax.text(val, i, " %.1f" % val, va="center", fontsize=8.5, color="0.15")

    ax = axes[1]
    ax.barh(y, lcl, height=0.55, color="#2a78d6")
    ax.axvline(TORNADIC_LCL_M, color="#d03b3b", lw=1.8)
    ax.text(TORNADIC_LCL_M, 0.03, "  tornadic LCL ~1200 m",
            transform=ax.get_xaxis_transform(), color="#d03b3b", fontsize=8.5,
            va="bottom", ha="left", zorder=5,
            bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85))
    ax.set_xlim(0, max(lcl) * 1.25)
    ax.set_title("LCL [m] — cloud base\nfar HIGHER than tornadic environments want",
                 fontsize=9.5)
    for i, val in enumerate(lcl):
        ax.text(val, i, " %.0f" % val, va="center", fontsize=8.5, color="0.15")

    # CAPE is the one quantity with an independent value, so show both and let the
    # disagreement be visible rather than averaged away.
    ax = axes[2]
    ax.barh(y - 0.17, cape, height=0.32, color="#2a78d6", label="engine (dz=100 m)")
    ax.barh(y + 0.17, cape_gfs, height=0.32, color="#eb6834", label="GFS own field")
    ax.set_xlim(0, max(max(cape), max(cape_gfs)) * 1.3)
    for i, (a, b) in enumerate(zip(cape, cape_gfs)):
        ax.text(a, i - 0.17, " %.0f" % a, va="center", fontsize=8, color="0.15")
        ax.text(b, i + 0.17, " %.0f" % b, va="center", fontsize=8, color="0.15")
    ax.set_title("CAPE [J/kg] — engine vs GFS\nindependent cross-check, not averaged",
                 fontsize=9.5)
    ax.legend(fontsize=8, loc="lower right", frameon=False)

    for ax in axes:
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=9)
        ax.invert_yaxis()
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="0.9", lw=0.7)
        ax.set_axisbelow(True)
    path = os.path.join(outdir, "centro_oeste_verdict.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--when", default="present", help="NCSS time, e.g. 2026-10-08T18:00:00Z")
    ap.add_argument("--outdir", default=os.path.join("outputs", "real_case_centro_oeste"))
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    rows = collect(args.when)
    print("%-20s %7s %7s %7s %7s %8s %s" % ("point", "CAPE", "gfsCAPE", "LCL", "shr6",
                                            "shr500", "geometry"))
    for r in rows:
        print("%-20s %7.0f %7.0f %7.0f %7.1f %8.2f %s" % (
            r["name"], r["diag"]["CAPE_J_kg"], r["gfs"]["CAPE_J_kg"], r["diag"]["LCL_m"],
            r["geom"]["shear_0_6km_m_s"], r["geom"]["shear_layer_m_s"],
            "ILL-CONDITIONED" if r["geom"]["geometry_is_ill_conditioned"] else "usable"))
    for p in (figure_hodographs(rows, args.outdir), figure_verdict(rows, args.outdir)):
        print("wrote", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
