"""(b/d) invert the observed (V_rot=39.49, sep=584 m) through the operator; edge-clip demo."""
import sys, numpy as np
sys.path.insert(0, 'src')
from atmospheric_data import radar_operator as ro
R_MOORE = 20250.0; ELEV = 0.5211372375488281

def run(core, bw=0.925, two_way=False, n_beam=9, n_gate=3, az_res=0.5, v_max=100.0, dx=20.):
    rad = ro.RadarSpec(beamwidth_deg=bw, elevation_deg=ELEV)
    cx, cy = 0.0, R_MOORE*np.cos(np.radians(ELEV))
    half = max(3000.0, 5.0*core)
    x = np.arange(cx-half, cx+half+1, dx); y = np.arange(cy-half, cy+half+1, dx)
    z = np.array([0.0, 204.0, 600.0])
    u, v, w = ro.rankine_vortex(x, y, z, (cx, cy), v_max, core)
    az, rg = ro.sweep_grid_for_domain(rad, cx, cy, half-400.0, az_resolution_deg=az_res)
    if two_way:
        o = ro._gauss_weights
        ro._gauss_weights = lambda n, b: ((np.array([0.]), np.array([1.])) if n<=1 else
            (lambda off: (off, np.exp(-8*np.log(2)*(off/b)**2)/np.exp(-8*np.log(2)*(off/b)**2).sum()))
            (np.linspace(-1.,1.,n)*b))
    try:
        sw = ro.synthetic_sweep(u, v, w, x, y, z, rad, az, rg, n_beam=n_beam, n_gate=n_gate)
    finally:
        if two_way: ro._gauss_weights = o
    r = ro.vrot_from_sweep(sw, max_separation_m=6000.0)
    return r['v_rot_m_s']/v_max, r['couplet_separation_m']

print("=== OPERATOR RESPONSE CURVE AT THE MOORE GEOMETRY (V_true normalised to 1) ===")
print(" core_m  2R/D   [1-way code] recov  sep_m   [2-way correct] recov  sep_m")
rows=[]
for core in (100,125,150,175,200,225,250,275,300,350,400,500,700,1000):
    f1,s1 = run(core); f2,s2 = run(core, two_way=True)
    rows.append((core,f1,s1,f2,s2))
    print("  %5d  %4.2f       %5.3f  %6.0f            %5.3f  %6.0f"
          % (core, 2*core/326.9, f1, s1, f2, s2))

print("\n=== INVERSION: what true peak tangential wind gives the OBSERVED V_rot=39.49? ===")
OBS_V, OBS_SEP = 39.49, 584.0
print("  (a) taking the study's assumption core = sep/2 = 292 m:")
for lbl, tw in (("1-way (as coded)", False), ("2-way (correct)", True)):
    f,s = run(292., two_way=tw)
    print("      %-18s recovery %.3f -> implied true peak %.1f m/s  (operator sep %.0f m vs obs 584)"
          % (lbl, f, OBS_V/f, s))
print("  (b) choosing the core that REPRODUCES the observed separation 584 m:")
for lbl, tw, idx in (("1-way (as coded)", False, (1,2)), ("2-way (correct)", True, (3,4))):
    cores = np.array([r[0] for r in rows]); seps = np.array([r[idx[1]] for r in rows])
    recs  = np.array([r[idx[0]] for r in rows])
    o = np.argsort(seps)
    core_hat = np.interp(OBS_SEP, seps[o], cores[o]); rec_hat = np.interp(OBS_SEP, seps[o], recs[o])
    print("      %-18s core_hat = %.0f m, recovery %.3f -> implied true peak %.1f m/s"
          % (lbl, core_hat, rec_hat, OBS_V/rec_hat))
print("  NOTE: observed separation is QUANTISED (0.5 deg rays = 177 m, 250 m gates); the operator")
print("        only produces sep in the same quantised set, so this inversion is coarse.")

print("\n=== EDGE-CLIP FAILURE MODE (synthetic_sweep clips sample points into the domain) ===")
rad = ro.RadarSpec(beamwidth_deg=0.925, elevation_deg=ELEV)
cx, cy = 0.0, R_MOORE*np.cos(np.radians(ELEV))
for half_dom, half_sweep in ((3000., 2600.), (1200., 1000.), (700., 600.)):
    x = np.arange(cx-half_dom, cx+half_dom+1, 20.); y = np.arange(cy-half_dom, cy+half_dom+1, 20.)
    z = np.array([0., 204., 600.])
    u,v,w = ro.rankine_vortex(x, y, z, (cx,cy), 100.0, 292.)
    az, rg = ro.sweep_grid_for_domain(rad, cx, cy, half_sweep, az_resolution_deg=0.5)
    sw = ro.synthetic_sweep(u,v,w,x,y,z,rad,az,rg,n_beam=9,n_gate=3)
    r = ro.vrot_from_sweep(sw, max_separation_m=6000.)
    # how many sample points would fall outside the domain?
    print("   domain half-width %5.0f m, sweep half-width %5.0f m -> V_rot=%6.2f (%.1f%% of the "
          "3 km-domain answer), sep=%5.0f" % (half_dom, half_sweep, r['v_rot_m_s'],
          100*r['v_rot_m_s']/24.72*24.72/24.72*1.0, r['couplet_separation_m']))
print("   (values differ only via clipping/extent, not physics -- the vortex is identical)")

print("\n=== z-EXTENT CLIP: the beam spans +-1 beamwidth in ELEVATION = +-327 m at 20 km ===")
print("   elevation offsets used:", np.round(ro._gauss_weights(5, 0.925)[0], 3), "deg")
h = [ro.beam_height_m(R_MOORE, ELEV+e) for e in ro._gauss_weights(5, 0.925)[0]]
print("   sampled heights: %s m  -> a model domain must span this range or it is silently clipped"
      % np.round(h,0))
