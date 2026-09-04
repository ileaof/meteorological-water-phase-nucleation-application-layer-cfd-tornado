"""(d) uncertainty budget inputs for the corrected 39.49 m/s target."""
import numpy as np
P = np.load('scratchpad/_polar.npz')
Gr, Gd, GX, GY, GZ, Gk = P['Gr'], P['Gd'], P['GX'], P['GY'], P['GZ'], P['Gk']
NY=float(P['NY']); TWO=2*NY; mx,my=P['mxy']
D = np.hypot(GX-mx, GY-my); fin = np.isfinite(Gd)

print("=== VELOCITY QUANTISATION ===")
u = np.unique(Gr[np.isfinite(Gr)])
print("distinct raw values: %d ; min step %.4f ; median step %.4f ; max |v| %.2f"
      % (u.size, np.diff(u).min(), np.median(np.diff(u)), np.abs(u).max()))
print("Nyquist %.4f ; Nyquist/step = %.1f (88D 8-bit velocity: 0.5 m/s LSB)" % (NY, NY/np.median(np.diff(u))))

print("\n=== COUPLET DEFINITION SENSITIVITY (dealiased, 5 km search) ===")
sel = fin & (D<=5000); idx = np.argwhere(sel)
XY = np.stack([GX[sel], GY[sel]], 1); V = Gd[sel]
best = {}
for maxsep in (300., 500., 750., 1000., 1500., 2000., 3000., 4000., 1e9):
    # exhaustive best pair under the separation constraint
    order_hi = np.argsort(V)[::-1][:400]; order_lo = np.argsort(V)[:400]
    bd, bi = -1e9, None
    for i in order_hi:
        d2 = np.hypot(XY[order_lo,0]-XY[i,0], XY[order_lo,1]-XY[i,1])
        ok = d2 <= maxsep
        if not ok.any(): continue
        j = order_lo[ok][np.argmin(V[order_lo][ok])]
        dv = V[i]-V[j]
        if dv > bd: bd, bi = dv, (i, j, float(np.hypot(XY[i,0]-XY[j,0], XY[i,1]-XY[j,1])))
    if bi:
        print("  max_sep <= %7.0f m : dV=%6.1f  V_rot=%6.2f  sep=%6.0f  in/out %+.1f/%+.1f"
              % (maxsep, bd, 0.5*bd, bi[2], V[bi[1]], V[bi[0]]))

print("\n=== STRICT GATE-TO-GATE (adjacent pair) DELTA-V, within 3 km of Moore ===")
for axis, nm in ((0,'azimuthal (ray to ray, same gate)'), (1,'radial (gate to gate, same ray)')):
    dv = np.diff(Gd, axis=axis); Dm = D[:,:-1] if axis==1 else D[:-1,:]
    f = np.isfinite(dv) & (Dm<=3000)
    print("  %-34s max |dV| = %6.1f -> V_rot = %5.2f" % (nm, np.abs(dv[f]).max(), 0.5*np.abs(dv[f]).max()))

print("\n=== SENSITIVITY TO THE ASSUMED 'MOORE' CENTRE (5 km search) ===")
for dxk in (-2000., -1000., 0., 1000., 2000.):
    for dyk in (-2000., 0., 2000.):
        DD = np.hypot(GX-(mx+dxk), GY-(my+dyk)); s = fin & (DD<=5000)
        idx2 = np.argwhere(s); v2 = Gd[s]
        ih = tuple(idx2[np.argmax(v2)]); il = tuple(idx2[np.argmin(v2)])
        sep = np.hypot(GX[ih]-GX[il], GY[ih]-GY[il])
        print("   centre offset (%+5.0f,%+5.0f) m : V_rot=%6.2f sep=%6.0f" %
              (dxk, dyk, 0.5*(Gd[ih]-Gd[il]), sep))

print("\n=== COUPLET AXIS ORIENTATION (rotation vs convergence contamination) ===")
s = fin & (D<=5000); idx3 = np.argwhere(s); v3 = Gd[s]
ih = tuple(idx3[np.argmax(v3)]); il = tuple(idx3[np.argmin(v3)])
dx_, dy_ = GX[ih]-GX[il], GY[ih]-GY[il]
# beam direction at the couplet
bx, by = GX[ih]+GX[il], GY[ih]+GY[il]; bn = np.hypot(bx,by); bx, by = bx/bn, by/bn
sep = np.hypot(dx_, dy_)
radial_comp = abs(dx_*bx + dy_*by); azim_comp = abs(-dx_*by + dy_*bx)
print("  hi=(ray %d,gate %d)  lo=(ray %d,gate %d)  -> d_ray=%d  d_gate=%d"
      % (ih[0], ih[1], il[0], il[1], ih[0]-il[0], ih[1]-il[1]))
print("  separation %.0f m = %.0f m ACROSS-beam (azimuthal) + %.0f m ALONG-beam (radial)"
      % (sep, azim_comp, radial_comp))
print("  tilt of the couplet axis from beam-perpendicular = %.1f deg" %
      np.degrees(np.arctan2(radial_comp, azim_comp)))
print("  => %.0f%% of dV is a pure azimuthal (rotational) signature; the rest mixes in "
      "along-beam convergence" % (100*azim_comp/sep))

print("\n=== HOW WELL SUPPORTED ARE THE TWO EXTREME GATES? ===")
for lbl, ij in (('outbound', ih), ('inbound', il)):
    a,g = ij
    nbrs = [Gd[(a-1)%Gd.shape[0],g], Gd[(a+1)%Gd.shape[0],g], Gd[a,g-1], Gd[a,g+1]]
    print("  %s v=%+.2f ; 4-neighbour values %s ; |v - median(nbrs)| = %.1f"
          % (lbl, Gd[ij], np.round(nbrs,1), abs(Gd[ij]-np.nanmedian(nbrs))))
# count of gates within 10% of the extreme
for lbl, ij, sgn in (('outbound', ih, +1), ('inbound', il, -1)):
    thr = 0.9*Gd[ij]
    n = ((Gd*sgn >= abs(thr)) & (D<=3000)).sum()
    print("  gates within 3 km of Moore reaching >=90%% of the %s peak: %d" % (lbl, n))
