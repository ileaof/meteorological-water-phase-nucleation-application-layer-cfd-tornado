"""First 3-D look at the refined nest: the rotating vortex, as a point cloud OR as
3-D **streamlines with arrows** spiralling inflow up into the updraft.

Runs a storm-following nest (as in examples/tornado_nest.py), saves the 3-D velocity /
vorticity fields to fields.npz, and renders (re-tune instantly with ``--load``):
  --mode streamlines : storm-relative streamlines (RK4 on the direction field) seeded in
                       a ring around the vortex axis, coloured by vertical velocity, with
                       arrowheads along each line -- the spiralling inflow + updraft.
  --mode points      : the cyclonic-vorticity column + updraft tower as a point cloud.

IDEALISED, under-resolved (dx ~ 444 m): this shows the storm's 3-D rotational structure,
NOT a resolved funnel-to-ground (that needs O(10-100 m) / adaptive AMR). Not a forecast.

    python examples/render_tornado_3d.py --device gpu --mode streamlines
    python examples/render_tornado_3d.py --load outputs/tornado_3d/fields.npz --mode streamlines
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np                                             # noqa: E402
import matplotlib                                             # noqa: E402
matplotlib.use("Agg")                                        # noqa: E402
import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib import cm                                     # noqa: E402
from PIL import Image                                         # noqa: E402


def _run_fields(args):
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    from storm_dynamics import nesting as nst, rotation as rot
    scfg = build_storm_config(
        preset="storm", nx=args.parent_nx, ny=args.parent_nx, nz=args.parent_nz,
        Lx=32000.0, Ly=32000.0, Lz=15000.0, duration=args.parent_duration, dt_max=3.0,
        hodograph_kind="quarter_circle", drag=True, z_stretch=1.05, U_max=args.u_max,
        z_turn=2000.0, C_s=0.22, device=args.device)
    parent = StormSimulation(scfg)
    print("maturing parent %dx%dx%d (dx=%.0f m) for %.0f s ..."
          % (parent.grid.nx, parent.grid.ny, parent.grid.nz, parent.grid.dx, args.parent_duration))
    parent.run(progress=lambda t, d, s: print("  parent t=%6.0f/%.0f" % (t, d), end="\r"))
    print(" " * 60, end="\r")
    _, _, wc = rot._centered_velocity(parent.state, parent.grid)
    wc = np.asarray(parent.grid.backend.to_cpu(wc))
    i, j = np.unravel_index(np.argmax(wc.max(axis=2)), wc.shape[:2])
    nc = parent.grid.nx; half = nc // 3
    i0 = int(np.clip(i - half // 2, 1, nc - half - 1))
    j0 = int(np.clip(j - half // 2, 1, nc - half - 1))
    spec = nst.NestSpec.aligned(parent.grid, i0=i0, j0=j0, ncx=half, ncy=half, refine=args.refine)
    print("nest (aligned) over the updraft; refine %d, storm-following%s ..."
          % (args.refine, " + COMPOSITE projection" if args.composite else ""))
    nest, rep = nst.run_concurrent_nest(parent, spec, window=args.window, follow=True,
                                        composite_projection=args.composite,
                                        les_boost=args.les_boost, cfl=args.cfl,
                                        progress=lambda t, d, s: print("  nest   t=%6.0f/%.0f" % (t, d), end="\r"))
    if args.composite:
        print("composite |div(rho0 u)| at interface: %.2e (max over sub-steps)"
              % rep["nest"]["composite_div_interface"])
    print(" " * 60, end="\r")
    g = nest.grid; to = g.backend.to_cpu
    uc, vc, wcn = rot._centered_velocity(nest.state, g)
    uc = np.asarray(to(uc)); vc = np.asarray(to(vc)); w = np.asarray(to(wcn))
    zeta = np.gradient(vc, g.dx, axis=0) - np.gradient(uc, g.dy, axis=1)
    cond = np.zeros_like(w)
    for nm in ("ql", "qi", "qr", "qs", "qg", "qh"):
        q = getattr(nest.state, nm, None)
        if q is not None:
            cond = cond + np.asarray(to(q))
    return dict(u=uc, v=vc, w=w, zeta=zeta, cond=cond,
                xc=np.asarray(to(g.xc)) / 1000.0, yc=np.asarray(to(g.yc)) / 1000.0,
                zc=np.asarray(to(g.zc)) / 1000.0, dx=float(g.dx))


def _streamlines(F, args):
    """RK4 on the storm-relative direction field, seeded in a ring around the vortex axis."""
    from scipy.interpolate import RegularGridInterpolator as RGI
    xc, yc, zc = F["xc"], F["yc"], F["zc"]
    grid = (xc, yc, zc)
    # reveal the mesocyclone: subtract the per-height horizontal MEAN wind (the ambient
    # storm-relative flow) so the streamlines follow the rotational + updraft perturbation
    up = F["u"] - F["u"].mean(axis=(0, 1), keepdims=True)
    vp = F["v"] - F["v"].mean(axis=(0, 1), keepdims=True)
    fu = RGI(grid, up, bounds_error=False, fill_value=0.0)
    fv = RGI(grid, vp, bounds_error=False, fill_value=0.0)
    fw = RGI(grid, F["w"], bounds_error=False, fill_value=0.0)
    gain = args.w_gain                                     # amplify vertical pitch (viz only)

    def vel_dir(p):                                        # direction field the line follows
        q = p[None, :]
        return np.array([float(fu(q)), float(fv(q)), gain * float(fw(q))])

    def vel(p):                                            # TRUE velocity (for w-colouring)
        q = p[None, :]
        return np.array([float(fu(q)), float(fv(q)), float(fw(q))])

    # vortex axis (horizontal centroid of the strong per-height vorticity)
    zeta = F["zeta"]; perlev = np.abs(zeta).max(axis=(0, 1), keepdims=True) + 1e-30
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    m = zeta / perlev > 0.4
    ax0, ay0 = float(X[m].mean()), float(Y[m].mean())

    rng = np.random.default_rng(0)
    seeds = []
    for zs in np.linspace(0.2, 6.0, 9):                    # heights
        for rr in (0.4, 0.8, 1.3, 2.0):                    # radii around the axis [km]
            for th in np.linspace(0, 2 * np.pi, 12, endpoint=False):
                seeds.append((ax0 + rr * np.cos(th), ay0 + rr * np.sin(th), zs))
    seeds = np.array(seeds) + rng.normal(0, 0.05, (len(seeds), 3))

    ds = 0.5 * F["dx"] / 1000.0                            # arc-length step [km]
    xmin, xmax = xc[0], xc[-1]; ymin, ymax = yc[0], yc[-1]; zmax = zc[-1]

    def dirf(q):
        V = vel_dir(q); n = np.linalg.norm(V)
        return V / n if n > 1e-6 else None

    lines = []
    for s in seeds:
        p = s.astype(float).copy(); path = [p.copy()]
        for _ in range(args.stream_steps):
            k1 = dirf(p)
            if k1 is None:
                break
            k2 = dirf(p + 0.5 * ds * k1); k3 = dirf(p + 0.5 * ds * (k2 if k2 is not None else k1))
            k4 = dirf(p + ds * (k3 if k3 is not None else k1))
            k2 = k2 if k2 is not None else k1; k3 = k3 if k3 is not None else k1
            k4 = k4 if k4 is not None else k1
            p = p + (ds / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            if not (xmin <= p[0] <= xmax and ymin <= p[1] <= ymax and 0 <= p[2] <= zmax):
                break
            path.append(p.copy())
        if len(path) > 8:
            lines.append((np.array(path), vel))
    return lines, (ax0, ay0), (fu, fv, fw)


def _flow_animation(F, args, outdir):
    """Animate the FLOW rotating (camera FIXED): particles advected by the
    mesocyclone's storm-relative velocity field, trails coloured by vertical
    velocity.  The field is a steady snapshot, so this shows the circulation
    pattern of that snapshot (the true spiralling inflow + updraft), not the
    storm's time evolution; the vertical pitch gain is a viz exaggeration."""
    from scipy.interpolate import RegularGridInterpolator as RGI
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    xc, yc, zc = F["xc"], F["yc"], F["zc"]
    grid = (xc, yc, zc)
    # same trick as the streamlines: remove the per-height MEAN wind so particles
    # follow the rotational + updraft perturbation, not the ambient flow
    up = F["u"] - F["u"].mean(axis=(0, 1), keepdims=True)
    vp = F["v"] - F["v"].mean(axis=(0, 1), keepdims=True)
    fu = RGI(grid, up, bounds_error=False, fill_value=0.0)
    fv = RGI(grid, vp, bounds_error=False, fill_value=0.0)
    fw = RGI(grid, F["w"], bounds_error=False, fill_value=0.0)
    gain = args.flow_w_gain                              # vertical pitch (viz only)
    wmax = max(float(np.nanmax(F["w"])), 1e-3)

    # vortex axis (horizontal centroid of the strong per-height vorticity)
    zeta = F["zeta"]; perlev = np.abs(zeta).max(axis=(0, 1), keepdims=True) + 1e-30
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    m = zeta / perlev > 0.4
    ax0, ay0 = float(X[m].mean()), float(Y[m].mean())
    span = 4.0
    x0, x1 = ax0 - span, ax0 + span; y0, y1 = ay0 - span, ay0 + span
    ztop = max(6.0, min(10.0, float(zc[-1])))            # keep the view on the vortex
    zmax = zc[-1]

    rng = np.random.default_rng(0)
    N, TRAIL = args.flow_particles, args.flow_trail
    ds = 0.12 * F["dx"] / 1000.0                          # RK4 arc step [km]

    def new_particles(n):                                 # rings around the axis
        zs = rng.uniform(0.15, 2.2, n)
        rr = rng.uniform(0.25, 1.5, n)
        th = rng.uniform(0.0, 2.0 * np.pi, n)
        return np.stack([ax0 + rr * np.cos(th), ay0 + rr * np.sin(th), zs], axis=1)

    def vel(P):                                           # (N,3): pitch-scaled
        return np.stack([fu(P), fv(P), gain * fw(P)], axis=1)

    def vdir(P):                                          # UNIT direction (arc-length RK4)
        V = vel(P)
        n = np.linalg.norm(V, axis=1, keepdims=True)
        return np.where(n > 1e-9, V / np.maximum(n, 1e-9), 0.0)

    def rk4(P):
        k1 = vdir(P); k2 = vdir(P + 0.5 * ds * k1)
        k3 = vdir(P + 0.5 * ds * k2); k4 = vdir(P + ds * k3)
        return P + (ds / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    def inside(P):
        return ((P[:, 0] > x0) & (P[:, 0] < x1) & (P[:, 1] > y0) & (P[:, 1] < y1)
                & (P[:, 2] > 0.0) & (P[:, 2] < ztop))

    cmap = plt.get_cmap("coolwarm")
    P = new_particles(N)
    trails = np.empty((N, TRAIL, 3))
    trails[:, 0, :] = P
    for t in range(1, TRAIL):                        # pre-advance a real trail history
        trails[:, t, :] = rk4(trails[:, t - 1, :])

    def draw(ax, trails):
        ax.clear()
        Tw = fw(trails.reshape(-1, 3)).reshape(N, TRAIL)
        cols = cmap(np.clip(0.5 + Tw[:, :-1] / (2 * wmax), 0, 1)).reshape(-1, 4)
        segs = np.stack([trails[:, :-1, :].reshape(-1, 3),
                         trails[:, 1:, :].reshape(-1, 3)], axis=1)
        ax.add_collection3d(Line3DCollection(segs, colors=cols, lw=1.1, alpha=0.85))
        ax.scatter(trails[:, -1, 0], trails[:, -1, 1], trails[:, -1, 2],
                   c=cols[-N:], s=4, depthshade=False)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_zlim(0, ztop)
        ax.set_box_aspect((1, 1, args.zexag * ztop / zc[-1])); ax.view_init(elev=45, azim=-60)
        ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]"); ax.set_zlabel("z [km]")
        ax.set_title("storm_dynamics 3-D: the FLOW rotating -- particles advected by the\n"
                     "mesocyclone circulation (colour = w; vertical pitch x%.1f)\n"
                     "IDEALISED (dx=%.0f m); steady-snapshot advection, camera fixed"
                     % (gain, F["dx"]), fontsize=9)

    fig = plt.figure(figsize=(8, 8.5)); ax = fig.add_subplot(111, projection="3d")
    draw(ax, trails); fig.tight_layout()
    still = os.path.join(outdir, "tornado_3d_flow_still.png")
    fig.savefig(still, dpi=130); print("still     ->", still)
    frames = []
    for _ in range(args.flow_frames):
        for _ in range(args.flow_substeps):
            P = rk4(P)
        bad = ~inside(P)
        if bad.any():
            nb = new_particles(int(bad.sum())); P[bad] = nb; trails[bad] = nb[:, None, :]
        trails = np.concatenate([trails[:, 1:, :], P[:, None, :]], axis=1)
        draw(ax, trails); fig.canvas.draw()
        buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
        w_, h_ = fig.canvas.get_width_height()
        frames.append(Image.fromarray(buf.reshape(h_, w_, 4)[:, :, :3].copy()))
    gif = os.path.join(outdir, "tornado_3d_flow_rotating.gif")
    frames[0].save(gif, save_all=True, append_images=frames[1:], duration=80, loop=0)
    plt.close(fig); print("animation ->", gif)
    print("IDEALISED, NOT a forecast; particle advection of a steady snapshot "
          "(vertical pitch x%.1f), not the storm's time evolution." % args.flow_w_gain)
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parent-nx", type=int, default=24)
    p.add_argument("--parent-nz", type=int, default=40)
    p.add_argument("--parent-duration", type=float, default=1200.0)
    p.add_argument("--u-max", type=float, default=18.0)
    p.add_argument("--refine", type=int, default=3)
    p.add_argument("--window", type=float, default=300.0)
    p.add_argument("--les-boost", type=float, default=1.4)
    p.add_argument("--cfl", type=float, default=0.20)
    p.add_argument("--composite", action="store_true",
                   help="docs/ROADMAP.md section 1: composite parent+nest mass-flux "
                        "projection every nest sub-step (div(rho0 u)=0 across the "
                        "coarse-fine interface; footprint already cell-aligned)")
    p.add_argument("--device", choices=["cpu", "gpu", "auto"], default="cpu")
    p.add_argument("--load", default=None, help="load fields.npz and just re-render")
    p.add_argument("--mode", choices=["streamlines", "points", "flow"], default="streamlines")
    p.add_argument("--flow-frames", type=int, default=40,
                   help="flow mode: animation frames (particle advection)")
    p.add_argument("--flow-particles", type=int, default=500,
                   help="flow mode: number of advected particles")
    p.add_argument("--flow-substeps", type=int, default=2,
                   help="flow mode: RK4 advection sub-steps per frame")
    p.add_argument("--flow-trail", type=int, default=30,
                   help="flow mode: trail length in advection steps")
    p.add_argument("--flow-w-gain", type=float, default=0.6,
                   help="flow mode: vertical pitch exaggeration (viz only)")
    p.add_argument("--stream-steps", type=int, default=180)
    p.add_argument("--w-gain", type=float, default=4.0,
                   help="amplify the vertical pitch of streamlines (viz only; rotation stays real)")
    p.add_argument("--zexag", type=float, default=1.6)
    p.add_argument("--outdir", default="outputs/tornado_3d")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    print("=" * 72)
    print("storm_dynamics -- 3-D %s of the refined storm-following nest" % args.mode)
    print("=" * 72)
    if args.load:
        d = np.load(args.load); F = {k: d[k] for k in d.files}
        print("loaded fields from", args.load)
    else:
        F = _run_fields(args)
        np.savez(os.path.join(args.outdir, "fields.npz"), **F)
        print("fields    ->", os.path.join(args.outdir, "fields.npz"))
    xc, yc, zc, dx = F["xc"], F["yc"], F["zc"], float(F["dx"])
    wmax = float(np.nanmax(F["w"]))
    print("nest dx=%.0f m | zeta_max=%.2e s^-1 | w_max=%.1f m/s" % (dx, float(np.nanmax(F["zeta"])), wmax))

    if args.mode == "flow":
        return _flow_animation(F, args, args.outdir)

    if args.mode == "streamlines":
        lines, (ax0, ay0), _ = _streamlines(F, args)
        print("streamlines integrated: %d" % len(lines))
        cmap = plt.get_cmap("coolwarm")
        span = 4.0
        x0, x1 = ax0 - span, ax0 + span; y0, y1 = ay0 - span, ay0 + span

        def draw(ax, azim):
            ax.clear()
            for path, vel in lines:
                wv = np.array([vel(q)[2] for q in path])
                col = cmap(np.clip(0.5 + wv[:-1] / (2 * max(wmax, 1e-3)), 0, 1))
                for a in range(len(path) - 1):
                    ax.plot(path[a:a + 2, 0], path[a:a + 2, 1], path[a:a + 2, 2],
                            color=col[a], lw=1.3, alpha=0.8)
                for a in range(10, len(path) - 1, 45):    # a couple of arrowheads per line
                    d = path[a + 1] - path[a]; n = np.linalg.norm(d)
                    if n > 1e-6:
                        d = d / n * 0.55
                        ax.quiver(path[a, 0], path[a, 1], path[a, 2], d[0], d[1], d[2],
                                  color=col[a], arrow_length_ratio=0.6, lw=1.1)
            ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_zlim(0, zc[-1])
            ax.set_box_aspect((1, 1, args.zexag)); ax.view_init(elev=12, azim=azim)
            ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]"); ax.set_zlabel("z [km]")
            ax.set_title("storm_dynamics 3-D mesocyclone streamlines: rotation + updraft "
                         "(blue=down, red=up; vertical pitch x%.0f)\nIDEALISED, under-resolved "
                         "(dx=%.0f m) -- rotational inflow, NOT a resolved funnel" % (args.w_gain, dx),
                         fontsize=9)
    else:
        zeta = F["zeta"]; w = F["w"]
        perlev = np.abs(zeta).max(axis=(0, 1), keepdims=True) + 1e-30
        znorm = zeta / perlev
        X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
        mv = znorm > 0.32; mu = w > 0.22 * wmax
        cx0 = float(X[mv | mu].mean()); cy0 = float(Y[mv | mu].mean()); span = 4.0

        def draw(ax, azim):
            ax.clear()
            if mu.any():
                ax.scatter(X[mu], Y[mu], Z[mu], c="royalblue", s=55, alpha=0.16, edgecolors="none")
            if mv.any():
                zz = znorm[mv]
                ax.scatter(X[mv], Y[mv], Z[mv], c=zz, cmap="autumn_r", vmin=0.32, vmax=1.0,
                           s=26 + 60 * (zz - 0.32) / 0.68, alpha=0.72, edgecolors="none")
            ax.set_xlim(cx0 - span, cx0 + span); ax.set_ylim(cy0 - span, cy0 + span); ax.set_zlim(0, zc[-1])
            ax.set_box_aspect((1, 1, args.zexag)); ax.view_init(elev=12, azim=azim)
            ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]"); ax.set_zlabel("z [km]")
            ax.set_title("storm_dynamics 3-D vorticity column (warm) + updraft (blue)\n"
                         "IDEALISED, under-resolved (dx=%.0f m)" % dx, fontsize=9)

    fig = plt.figure(figsize=(8, 8.5)); ax = fig.add_subplot(111, projection="3d")
    draw(ax, azim=-60); fig.tight_layout()
    still = os.path.join(args.outdir, "tornado_3d_%s_still.png" % args.mode)
    fig.savefig(still, dpi=130); print("still     ->", still)
    frames = []
    for az in range(-90, 270, 12):
        draw(ax, azim=az); fig.canvas.draw()
        buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
        w_, h_ = fig.canvas.get_width_height()
        frames.append(Image.fromarray(buf.reshape(h_, w_, 4)[:, :, :3].copy()))
    gif = os.path.join(args.outdir, "tornado_3d_%s_rotating.gif" % args.mode)
    frames[0].save(gif, save_all=True, append_images=frames[1:], duration=90, loop=0)
    plt.close(fig); print("animation ->", gif)
    print("IDEALISED, NOT a forecast; under-resolved -- rotational structure, not a funnel-to-ground.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
