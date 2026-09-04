import numpy as np
d = np.load('outputs/nexrad_moore/ktlx_velocity_dealiased.npz', allow_pickle=True)
x, y = d['x_m'], d['y_m']
r = np.hypot(x, y); az = np.degrees(np.arctan2(x, y)) % 360.0
# azimuth index
au = np.sort(np.unique(np.round(az, 4)))
print("first 8 az:", au[:8]); print("diff stats:", np.median(np.diff(au)), np.diff(au).min(), np.diff(au).max())
r0 = r.min()
ir_f = (r - r0)/249.99
print("range index residual max", np.abs(ir_f-np.round(ir_f)).max())
