"""3-D render of the resolved low-level mesocyclone from the real-data cascade (28 m):
the vertical-vorticity tube (isosurface) threaded by storm-relative streamlines coloured by
vertical velocity (red updraft / blue downdraft).  Off-screen PyVista -> PNG."""
import os
import numpy as np
import pyvista as pv

REPO = r"c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
MED = os.path.join(REPO, "docs", "media", "storm")
f = np.load(os.path.join(REPO, "outputs", "moore_cascade", "fields_cascade.npz"))
u, v, w, zeta = f["u"], f["v"], f["w"], f["zeta"]
xc, yc, zc, dx = f["xc"], f["yc"], f["zc"], float(f["dx"])

# crop: interior (drop the sponge border) + the low/mid levels where the meso lives (z < 3.5 km)
nb = zeta.shape[0] // 6
kz = int(np.argmin(np.abs(zc - 3500)))
sx = slice(nb, zeta.shape[0] - nb); sy = slice(nb, zeta.shape[1] - nb); sz = slice(0, kz + 1)
u = u[sx, sy, sz]; v = v[sx, sy, sz]; w = w[sx, sy, sz]; zeta = zeta[sx, sy, sz]
X = (xc[sx] - xc[sx][len(xc[sx]) // 2]) / 1e3          # km, centred
Y = (yc[sy] - yc[sy][len(yc[sy]) // 2]) / 1e3
Z = zc[sz] / 1e3                                        # km
# storm-relative horizontal wind (reveals the rotation)
u = u - u.mean(axis=(0, 1), keepdims=True); v = v - v.mean(axis=(0, 1), keepdims=True)

Xg, Yg, Zg = np.meshgrid(X, Y, Z, indexing="ij")
grid = pv.StructuredGrid(Xg, Yg, Zg)
grid["zeta"] = zeta.flatten(order="F")
grid["w"] = w.flatten(order="F")
vec = np.column_stack([u.flatten(order="F"), v.flatten(order="F"), w.flatten(order="F")])
grid["vel"] = vec

zpk = float(np.percentile(np.abs(zeta), 99.7))
pv.set_plot_theme("dark")
pl = pv.Plotter(off_screen=True, window_size=(1200, 1500))
pl.set_background("#0b0e14", top="#1a2030")

# (1) the vorticity tube -- one clean isosurface of cyclonic (positive) vorticity, largest
#     connected body only (drops floating fragments aloft)
iso = grid.contour(isosurfaces=[0.60 * zpk], scalars="zeta")
try:
    iso = iso.connectivity("largest")
except Exception:
    pass
if iso.n_points:
    pl.add_mesh(iso, color="#e8663c", opacity=0.30, show_scalar_bar=False, smooth_shading=True)

# (2) storm-relative streamlines through the vortex, coloured by vertical velocity
seed = pv.Disc(center=(0, 0, 0.25), inner=0.0, outer=0.55, normal=(0, 0, 1), r_res=4, c_res=28)
strm = grid.streamlines_from_source(seed, vectors="vel", integration_direction="both",
                                    max_length=30.0)
tube = strm.tube(radius=0.008)
wm = max(float(np.nanmax(np.abs(w))), 1e-3)
pl.add_mesh(tube, scalars="w", cmap="RdBu_r", clim=[-wm, wm],
            scalar_bar_args=dict(title="w  [m/s]", color="white", title_font_size=18,
                                 label_font_size=14, position_x=0.83, position_y=0.06, width=0.12))

# a faint ground plane for depth
ground = pv.Plane(center=(0, 0, 0.0), direction=(0, 0, 1),
                  i_size=X[-1] - X[0], j_size=Y[-1] - Y[0])
pl.add_mesh(ground, color="#20242e", opacity=0.5, show_scalar_bar=False)

pl.set_scale(zscale=0.85)                              # gentle vertical exaggeration control
pl.camera.azimuth = 25; pl.camera.elevation = 18
pl.camera.zoom(1.5)
pl.add_text("Low-level mesocyclone in 3-D  |  real Moore 2013 env  |  AMR cascade to 28 m",
            position="upper_left", font_size=11, color="white")
pl.add_text("vorticity tube (red) + storm-relative streamlines coloured by updraft/downdraft",
            position="lower_left", font_size=9, color="#9db4c0")
out = os.path.join(MED, "meso_3d_28m.png")
pl.screenshot(out)
print("saved", out, "| iso pts", iso.n_points, "| streamline pts", strm.n_points, "| zeta_pk", zpk)
