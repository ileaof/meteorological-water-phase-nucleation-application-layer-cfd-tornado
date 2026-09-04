import numpy as np
d = np.load('outputs/nexrad_moore/ktlx_velocity_dealiased.npz', allow_pickle=True)
x, y, z = d['x_m'], d['y_m'], d['z_m']
vr, vd = d['vr_raw'].astype(float), d['vr_dealiased'].astype(float)
nyq = float(d['nyquist_m_s']); elev = float(d['elevation_deg'])
mx, my = d['moore_xy']
r = np.hypot(x, y)
az = np.degrees(np.arctan2(x, y)) % 360.0
print("N gates", x.size, " nyq", nyq, " 2*nyq", 2*nyq, " elev", elev)
print("range  %.1f .. %.1f m" % (r.min(), r.max()))
print("z      %.1f .. %.1f m" % (z.min(), z.max()))
print("raw    %+.2f .. %+.2f" % (vr.min(), vr.max()))
print("dealia %+.2f .. %+.2f" % (vd.min(), vd.max()))
print("Moore at (%.0f,%.0f) range %.2f km az %.2f deg" % (mx, my, np.hypot(mx,my)/1e3,
      np.degrees(np.arctan2(mx,my))%360))
# gate/azimuth spacing
ru = np.unique(np.round(r, 1))
print("unique-ish range count", ru.size, "median drange", np.median(np.diff(ru)))
au = np.unique(np.round(az, 3))
print("unique az count", au.size)
