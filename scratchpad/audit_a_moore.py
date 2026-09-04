"""(a) the MOORE couplet itself: is it a coherent, adjacent, non-artefactual pair?"""
import numpy as np
P = np.load('scratchpad/_polar.npz')
Gr, Gd, GX, GY, GZ, Gk = P['Gr'], P['Gd'], P['GX'], P['GY'], P['GZ'], P['Gk']
NY = float(P['NY']); mx, my = P['mxy']
D = np.hypot(GX-mx, GY-my)
fin = np.isfinite(Gd)

def couplet(R, G=Gd, name=''):
    sel = fin & (D <= R)
    if sel.sum() < 10: return None
    idx = np.argwhere(sel)
    vals = G[sel]
    ih = idx[np.nanargmax(vals)]; il = idx[np.nanargmin(vals)]
    hi, lo = tuple(ih), tuple(il)
    sep = float(np.hypot(GX[hi]-GX[lo], GY[hi]-GY[lo]))
    dv = float(G[hi]-G[lo])
    rng = 0.5*(np.hypot(GX[hi],GY[hi])+np.hypot(GX[lo],GY[lo]))
    return dict(hi=hi, lo=lo, sep=sep, dv=dv, vrot=0.5*dv, rng=rng,
                vhi=float(G[hi]), vlo=float(G[lo]),
                khi=float(Gk[hi]), klo=float(Gk[lo]),
                zhi=float(GZ[hi]), zlo=float(GZ[lo]), n=int(sel.sum()))

print("=== SEARCH-RADIUS SENSITIVITY (dealiased field, Moore-centred) ===")
print(" R_km  ngates    dV   V_rot   sep_m  range_km  in/out (fold k)      z_hi/z_lo")
for R in (1e3,1.5e3,2e3,2.5e3,3e3,4e3,5e3,6e3,8e3,10e3,15e3):
    c = couplet(R)
    if c is None: print("  %.1f  --"%(R/1e3)); continue
    print(" %5.1f %6d %6.1f %6.2f %7.0f %8.2f   %+6.1f/%+6.1f (%+d/%+d) %5.0f/%5.0f"
          % (R/1e3, c['n'], c['dv'], c['vrot'], c['sep'], c['rng']/1e3,
             c['vlo'], c['vhi'], c['klo'], c['khi'], c['zlo'], c['zhi']))

print("\n=== SAME, on the RAW (aliased) field, for contrast ===")
for R in (3e3,5e3,8e3):
    c = couplet(R, Gr)
    print(" %5.1f %6d %6.1f %6.2f %7.0f %8.2f   %+6.1f/%+6.1f"
          % (R/1e3, c['n'], c['dv'], c['vrot'], c['sep'], c['rng']/1e3, c['vlo'], c['vhi']))

c = couplet(5e3)
hi, lo = c['hi'], c['lo']
print("\n=== LOCAL NEIGHBOURHOOD OF THE 5 km COUPLET MEMBERS ===")
print(" (ray, gate) indices: hi=%s lo=%s   -> ray gap %d, gate gap %d"
      % (hi, lo, hi[0]-lo[0], hi[1]-lo[1]))
def patch(c_, name, G):
    a, g = c_
    print("\n %s gate: v=%+.2f  k=%+d  x=%.0f y=%.0f z=%.0f  d_Moore=%.0f m" %
          (name, G[a,g], Gk[a,g], GX[a,g], GY[a,g], GZ[a,g], D[a,g]))
    print("   5x5 dealiased patch (rows=rays, cols=range gates):")
    for aa in range(a-2, a+3):
        print("    " + " ".join("%7.2f" % G[aa % G.shape[0], gg] for gg in range(g-2, g+3)))
    print("   5x5 RAW patch:")
    for aa in range(a-2, a+3):
        print("    " + " ".join("%7.2f" % Gr[aa % G.shape[0], gg] for gg in range(g-2, g+3)))
    print("   5x5 fold k:")
    for aa in range(a-2, a+3):
        print("    " + " ".join("%+7.0f" % Gk[aa % G.shape[0], gg] for gg in range(g-2, g+3)))
patch(hi, 'OUTBOUND(max)', Gd)
patch(lo, 'INBOUND (min)', Gd)

# robustness: 3x3 median-smoothed field
print("\n=== ROBUSTNESS: couplet from a 3x3 MEDIAN-smoothed field (kills single-gate spikes) ===")
from numpy.lib.stride_tricks import sliding_window_view
def med3(G):
    Gp = np.pad(G, 1, constant_values=np.nan)
    W = sliding_window_view(Gp, (3,3)).reshape(G.shape[0], G.shape[1], 9)
    return np.nanmedian(W, axis=2)
Gm = med3(Gd); Grm = med3(Gr)
finm = np.isfinite(Gm)
def couplet2(R, G, F):
    sel = F & (D<=R); idx = np.argwhere(sel); vals = G[sel]
    ih = tuple(idx[np.nanargmax(vals)]); il = tuple(idx[np.nanargmin(vals)])
    sep = float(np.hypot(GX[ih]-GX[il], GY[ih]-GY[il])); dv = float(G[ih]-G[il])
    return dv, 0.5*dv, sep, float(G[il]), float(G[ih]), ih, il
for R in (3e3,5e3,8e3):
    dv, vrot, sep, vlo, vhi, ih, il = couplet2(R, Gm, finm)
    print(" median3  R=%.0fkm  dV=%6.1f  V_rot=%6.2f  sep=%6.0f m  in/out %+.1f/%+.1f"
          % (R/1e3, dv, vrot, sep, vlo, vhi))

print("\n=== ROBUSTNESS: drop ALL unfolded (k!=0) gates, then re-extract ===")
Gd_nofold = np.where(Gk==0, Gd, np.nan)
for R in (3e3,5e3,8e3):
    dv, vrot, sep, vlo, vhi, ih, il = couplet2(R, Gd_nofold, np.isfinite(Gd_nofold))
    print(" no-fold R=%.0fkm  dV=%6.1f  V_rot=%6.2f  sep=%6.0f m  in/out %+.1f/%+.1f"
          % (R/1e3, dv, vrot, sep, vlo, vhi))

print("\n=== HOW MANY GATES NEAR MOORE WERE UNFOLDED? ===")
for R in (1e3,2e3,3e3,5e3,8e3):
    s = fin & (D<=R)
    print("  R=%.0f km : %5d gates, %4d unfolded (%.2f%%), k values %s"
          % (R/1e3, s.sum(), (Gk[s]!=0).sum(), 100*(Gk[s]!=0).mean(),
             np.unique(Gk[s][~np.isnan(Gk[s])]).astype(int)))
