"""(b) forward-operator audit: geometry checks, weighting-form bug, quadrature convergence."""
import sys, os, numpy as np
sys.path.insert(0, 'src')
from atmospheric_data import radar_operator as ro

R_MOORE = 20250.0; ELEV = 0.5211372375488281
print("=== GEOMETRY CROSS-CHECK AGAINST THE PY-ART GATE ALTITUDES ===")
for r in (19300., 20120., 20250.):
    print("  slant %.0f m, elev %.4f -> beam_height_m = %.1f m (Py-ART gate z at the couplet: 208-211 m)"
          % (r, ELEV, ro.beam_height_m(r, ELEV)))
rad = ro.RadarSpec(beamwidth_deg=0.925, elevation_deg=ELEV)
print("  beam diameter at 20.25 km, 0.925 deg = %.1f m" % rad.beam_diameter_m(R_MOORE))
print("  beam diameter at 20.25 km, 1.02  deg = %.1f m  (super-res EFFECTIVE beamwidth)"
      % ro.RadarSpec(beamwidth_deg=1.02).beam_diameter_m(R_MOORE))
print("  NOTE beam_height_m returns height above the RADAR ANTENNA, not AGL at the target.")

print("\n=== ONE-WAY vs TWO-WAY BEAM WEIGHTING ===")
off, w = ro._gauss_weights(101, 0.925)
half = off[np.argmin(np.abs(w/w.max() - 0.5))]
print("  code weight falls to 0.5 at phi = %.4f deg -> implied full width %.4f deg" % (half, 2*abs(half)))
print("  that is the ONE-WAY power pattern exp(-4 ln2 (phi/theta)^2).")
w2 = np.exp(-8*np.log(2)*(off/0.925)**2)
half2 = off[np.argmin(np.abs(w2/w2.max() - 0.5))]
print("  correct TWO-WAY exp(-8 ln2 (phi/theta)^2) halves at %.4f deg -> width %.4f deg (= theta/sqrt2)"
      % (half2, 2*abs(half2)))
sd1 = np.sqrt((w*off**2).sum()/w.sum()); sd2 = np.sqrt((w2*off**2).sum()/w2.sum())
print("  effective sigma: code %.4f deg vs two-way %.4f deg  -> code smooths %.2fx too wide"
      % (sd1, sd2, sd1/sd2))

# ---- quantify the effect on the recovery factor at the ACTUAL Moore geometry -------------
def sweep_vrot(core_m, v_max, beamwidth, two_way=False, n_beam=5, n_gate=3, dx=25.0,
               range_m=R_MOORE, elev=ELEV, az_res=0.5):
    rad = ro.RadarSpec(beamwidth_deg=beamwidth, elevation_deg=elev)
    cx, cy = 0.0, range_m*np.cos(np.radians(elev))
    half = max(3000.0, 4.0*core_m)
    x = np.arange(cx-half, cx+half+1, dx); y = np.arange(cy-half, cy+half+1, dx)
    z = np.array([0.0, 204.0, 500.0, 1200.0])
    u, v, w = ro.rankine_vortex(x, y, z, (cx, cy), v_max, core_m)
    az, rg = ro.sweep_grid_for_domain(rad, cx, cy, half-400.0, az_resolution_deg=az_res)
    if two_way:
        orig = ro._gauss_weights
        def tw(n, bw):
            if n <= 1: return np.array([0.0]), np.array([1.0])
            o = np.linspace(-1., 1., n)*float(bw)
            ww = np.exp(-8.0*np.log(2.0)*(o/float(bw))**2)
            return o, ww/ww.sum()
        ro._gauss_weights = tw
    try:
        sw = ro.synthetic_sweep(u, v, w, x, y, z, rad, az, rg, n_beam=n_beam, n_gate=n_gate)
    finally:
        if two_way: ro._gauss_weights = orig
    return ro.vrot_from_sweep(sw, max_separation_m=4000.0)

print("\n=== RECOVERY FACTOR AT THE MOORE GEOMETRY (range 20.25 km, elev 0.5211) ===")
print("  core   2R/D    V_rot(code,1-way)  V_rot(true 2-way)   ratio   sep_code")
for core in (125., 164., 250., 292., 500., 1000., 2000.):
    a = sweep_vrot(core, 39.49, 0.925, two_way=False)
    b = sweep_vrot(core, 39.49, 0.925, two_way=True)
    D = ro.RadarSpec(beamwidth_deg=0.925).beam_diameter_m(R_MOORE)
    print("  %5.0f  %5.2f      %6.2f (%.0f%%)      %6.2f (%.0f%%)     %.3f   %6.0f"
          % (core, 2*core/D, a['v_rot_m_s'], 100*a['v_rot_m_s']/39.49,
             b['v_rot_m_s'], 100*b['v_rot_m_s']/39.49,
             a['v_rot_m_s']/b['v_rot_m_s'], a['couplet_separation_m']))

print("\n=== QUADRATURE CONVERGENCE (core 292 m, the value implied by sep/2=292) ===")
for nb_ in (1,3,5,7,9,13,17):
    r_ = sweep_vrot(292., 39.49, 0.925, n_beam=nb_, n_gate=3)
    print("   n_beam=%2d  V_rot=%6.3f  sep=%6.0f" % (nb_, r_['v_rot_m_s'], r_['couplet_separation_m']))
for ng in (1,3,5,9):
    r_ = sweep_vrot(292., 39.49, 0.925, n_beam=9, n_gate=ng)
    print("   n_gate=%2d  V_rot=%6.3f  sep=%6.0f" % (ng, r_['v_rot_m_s'], r_['couplet_separation_m']))

print("\n=== SENSITIVITY TO THE ASSUMED BEAMWIDTH (super-res effective width) ===")
for bw in (0.89, 0.925, 1.02, 1.25, 1.39):
    r_ = sweep_vrot(292., 39.49, bw, n_beam=9)
    print("   beamwidth %.3f deg (D=%4.0f m): V_rot=%6.2f (%.0f%% of truth)"
          % (bw, ro.RadarSpec(beamwidth_deg=bw).beam_diameter_m(R_MOORE),
             r_['v_rot_m_s'], 100*r_['v_rot_m_s']/39.49))

print("\n=== AZIMUTHAL SAMPLING (the obs is 0.5 deg super-res) ===")
for azr in (0.25, 0.5, 1.0):
    r_ = sweep_vrot(292., 39.49, 0.925, n_beam=9, az_res=azr)
    print("   az_resolution %.2f deg: V_rot=%6.2f sep=%6.0f" % (azr, r_['v_rot_m_s'], r_['couplet_separation_m']))
