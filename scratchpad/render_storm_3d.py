"""3-D rendering of the simulated storm -- cloud, precipitation, updraft, vortex, and the FUNNEL.

The goal is a picture of the storm with a condensation funnel.  A funnel is PHYSICS, not
graphics: it is cloud forming BELOW the ambient cloud base because the vortex's pressure deficit
lowers the local saturation level.  So this script does two things and keeps them separate:

  1. renders what the model ACTUALLY produced -- cloud (ql+qi), rain shaft (qr), updraft (w) and
     the vortex core (|zeta|), which together are the supercell and its mesocyclone;
  2. evaluates an explicit FUNNEL CRITERION -- condensate present below the ambient LCL, inside
     the vortex -- and reports honestly whether one exists, rather than drawing something
     suggestive and letting the picture imply it.

The funnel test is written now so that the moment a run produces the pressure deficit for one, the
image appears by itself.  Current status is printed with the numbers behind it.

    python scratchpad/render_storm_3d.py [path/to/fields.npz]

Env: FIELD=<npz>  OUT=<png stem>  AZIM=<deg>  ZEXAG=<factor>
"""
import os, sys
REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
sys.path.insert(0, os.path.join(REPO, "src"))
os.chdir(REPO)
import numpy as np, warnings; warnings.filterwarnings("ignore")

FIELD = os.environ.get("FIELD", sys.argv[1] if len(sys.argv) > 1
                       else "outputs/parent_matured_120_48_2800.npz")
OUTSTEM = os.environ.get("OUT", "docs/media/storm/storm_3d")
ZEXAG = float(os.environ.get("ZEXAG", 2.5))       # vertical exaggeration (storms are wide + flat)
os.makedirs(os.path.dirname(OUTSTEM), exist_ok=True)


def load_fields(path):
    """Return (fields dict, x, y, z) with z the true stretched model levels."""
    z = np.load(path)
    f = {k: z[k] for k in z.files if k != "_t"}
    if "xc" in f and "zc" in f:                    # nest dumps carry their own coordinates
        return f, np.asarray(f["xc"]), np.asarray(f["yc"]), np.asarray(f["zc"])
    # parent maturation caches store only the state -> rebuild the grid it was written on
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    nx, ny, nz = f["theta"].shape
    cfg = build_storm_config(preset="storm", nx=nx, ny=ny, nz=nz, Lx=nx * 600.0, Ly=ny * 600.0,
                             Lz=15000.0, duration=1.0, dt_max=3.0, z_stretch=1.05, device="cpu")
    g = StormSimulation(cfg).grid
    to = g.backend.to_cpu
    return f, np.asarray(to(g.xc)), np.asarray(to(g.yc)), np.asarray(to(g.zc))


def centred(f):
    """Cell-centred u, v, w from the staggered state."""
    u, v, w = f["u"], f["v"], f["w"]
    uc = 0.5 * (u[:-1] + u[1:]) if u.shape[0] > f["theta"].shape[0] else u
    vc = 0.5 * (v[:, :-1] + v[:, 1:]) if v.shape[1] > f["theta"].shape[1] else v
    wc = 0.5 * (w[:, :, :-1] + w[:, :, 1:]) if w.shape[2] > f["theta"].shape[2] else w
    return uc, vc, wc


def funnel_report(f, x, y, z, lcl_m=1068.0, ql_thresh=1e-5):
    """Is there a CONDENSATION FUNNEL?

    A funnel is cloud hanging BELOW the ambient cloud base in a column only a few hundred metres
    wide, because the vortex's core pressure deficit lowers the local saturation level.  Two
    weaker criteria were tried first and BOTH gave false positives, so they are recorded here:

      * "condensate exists below the ambient LCL" -- fires trivially, because the model's cloud
        base (706 m) is already below the analytic LCL (1068 m) over much of the domain;
      * "cloud base over the vortex is lower than a ring 10-20 cells away" -- reported a 2610 m
        'depression', but that ring is 6-12 km out, i.e. OUTSIDE the convective tower where only
        anvil exists at 3-7 km.  It compared "under the storm" with "beside the storm".

    The honest test is whether the cloud base over the vortex is BELOW THE AMBIENT LCL, measured
    against the surrounding STORM BODY rather than against clear air.  And note the hard limit:
    a funnel is 100-300 m across, so at dx >= 300 m it is sub-grid and the model cannot represent
    one however strong the vortex.
    """
    ql = f.get("ql"); qi = f.get("qi")
    cond = (ql if ql is not None else 0.0) + (qi if qi is not None else 0.0)
    uc, vc, _ = centred(f)
    dx = float(x[1] - x[0]); dy = float(y[1] - y[0])
    zeta = np.gradient(vc, dx, axis=0) - np.gradient(uc, dy, axis=1)

    # cloud base per column
    base = np.full(cond.shape[:2], np.nan)
    lit_any = cond > ql_thresh
    for i in range(cond.shape[0]):
        for j in range(cond.shape[1]):
            lit = np.where(lit_any[i, j])[0]
            if lit.size:
                base[i, j] = z[lit[0]]

    k = int(np.argmin(np.abs(z - 100.0)))
    zl = np.abs(zeta[:, :, k])
    vi, vj = np.unravel_index(int(np.argmax(zl)), zl.shape)
    ii, jj = np.mgrid[0:base.shape[0], 0:base.shape[1]]
    d = np.hypot(ii - vi, jj - vj)
    # STORM BODY = nearby columns that have a LOW cloud base (inside the tower), not the anvil
    body = (d > 2) & (d < 10) & np.isfinite(base) & (base < lcl_m + 1500.0)
    body_base = float(np.nanmedian(base[body])) if body.any() else float("nan")
    vort_base = float(base[vi, vj])
    sub_grid = dx > 300.0
    below_lcl = np.isfinite(vort_base) and (vort_base < lcl_m - 100.0)
    localized = (np.isfinite(body_base) and np.isfinite(vort_base)
                 and (body_base - vort_base) > 200.0)
    return {"lcl_m": lcl_m, "dx_m": dx, "funnel_is_subgrid": bool(sub_grid),
            "vortex_ij": (int(vi), int(vj)), "zeta_at_vortex": float(zl[vi, vj]),
            "cloud_base_over_vortex_m": vort_base,
            "cloud_base_storm_body_m": body_base,
            "local_depression_m": (body_base - vort_base
                                   if np.isfinite(body_base) and np.isfinite(vort_base)
                                   else float("nan")),
            "below_ambient_lcl": bool(below_lcl),
            "funnel": bool(below_lcl and localized and not sub_grid)}


def render(f, x, y, z, stem):
    import pyvista as pv
    pv.OFF_SCREEN = True
    uc, vc, wc = centred(f)
    dx = float(x[1] - x[0]); dy = float(y[1] - y[0])
    zeta = np.gradient(vc, dx, axis=0) - np.gradient(uc, dy, axis=1)
    cond = (f.get("ql", 0.0) + f.get("qi", 0.0))
    rain = f.get("qr", np.zeros_like(zeta))

    X, Y, Z = np.meshgrid(x / 1000.0, y / 1000.0, z / 1000.0 * ZEXAG, indexing="ij")
    grid = pv.StructuredGrid(X, Y, Z)
    grid["cloud"] = cond.ravel(order="F")
    grid["rain"] = np.asarray(rain).ravel(order="F")
    grid["w"] = np.asarray(wc).ravel(order="F")
    grid["zeta"] = np.abs(zeta).ravel(order="F")

    p = pv.Plotter(off_screen=True, window_size=(1600, 1200))
    p.set_background("white")
    added = []
    def iso(name, val, color, opacity, label):
        if val is None or not np.isfinite(val) or val <= 0:
            return
        try:
            s = grid.contour([val], scalars=name)
            if s.n_points:
                p.add_mesh(s, color=color, opacity=opacity, smooth_shading=True)
                added.append("%s (%s)" % (label, color))
        except Exception:
            pass

    iso("cloud", max(1e-5, 0.25 * float(np.max(cond))), "white", 0.28, "cloud ql+qi")
    iso("rain", max(1e-5, 0.35 * float(np.max(rain))), "royalblue", 0.35, "rain qr")
    iso("w", max(2.0, 0.45 * float(np.max(wc))), "orangered", 0.55, "updraft w")
    iso("zeta", max(2e-3, 0.5 * float(np.max(np.abs(zeta)))), "gold", 0.85, "vortex |zeta|")

    p.add_axes(xlabel="x [km]", ylabel="y [km]", zlabel="z (x%.1f)" % ZEXAG)
    p.camera_position = "xz"
    p.camera.azimuth = float(os.environ.get("AZIM", 35.0))
    p.camera.elevation = 18.0
    p.reset_camera()
    out = stem + ".png"
    p.screenshot(out)
    p.close()
    return out, added


if __name__ == "__main__":
    print("field: %s" % FIELD)
    f, x, y, z = load_fields(FIELD)
    print("grid %s, dx=%.0f m, z %.0f..%.0f m" % (f["theta"].shape, x[1] - x[0], z[0], z[-1]))
    r = funnel_report(f, x, y, z)
    print()
    print("FUNNEL CRITERION -- cloud below the ambient LCL, localised on the vortex")
    print("   ambient LCL                        %.0f m" % r["lcl_m"])
    print("   model dx                           %.0f m %s" % (
        r["dx_m"], "(a funnel is 100-300 m across => SUB-GRID)" if r["funnel_is_subgrid"] else ""))
    print("   low-level vortex at                (%d,%d), |zeta|=%.4f 1/s"
          % (r["vortex_ij"][0], r["vortex_ij"][1], r["zeta_at_vortex"]))
    print("   cloud base OVER the vortex         %.0f m" % r["cloud_base_over_vortex_m"])
    print("   cloud base in the storm body       %.0f m" % r["cloud_base_storm_body_m"])
    print("   local depression                   %+.0f m" % r["local_depression_m"])
    print("   below the ambient LCL?             %s" % r["below_ambient_lcl"])
    print("   => FUNNEL PRESENT: %s" % r["funnel"])
    if not r["funnel"]:
        why = []
        if r["funnel_is_subgrid"]:
            why.append("dx=%.0f m cannot represent a 100-300 m funnel" % r["dx_m"])
        if not r["below_ambient_lcl"]:
            why.append("cloud base over the vortex (%.0f m) is NOT below the ambient LCL (%.0f m)"
                       % (r["cloud_base_over_vortex_m"], r["lcl_m"]))
        print("      reason: " + "; ".join(why))
    print()
    out, added = render(f, x, y, z, OUTSTEM)
    print("rendered %s" % out)
    for a in added:
        print("   isosurface: %s" % a)
