import numpy as np
np.set_printoptions(linewidth=250)
P = np.load('scratchpad/_polar.npz')
Gr, Gd, GX, GY, GZ, Gk = P['Gr'], P['Gd'], P['GX'], P['GY'], P['GZ'], P['Gk']
NY=float(P['NY']); mx,my=P['mxy']
D = np.hypot(GX-mx, GY-my)
a0, g0 = 539, 72          # outbound couplet member
A = range(a0-9, a0+8); G = range(g0-8, g0+9)
def show(F, name, fmt="%6.1f"):
    print("\n--- %s  (rows = rays %d..%d, cols = gates %d..%d ; ray step 0.5deg=175m, gate 250m) ---"
          % (name, a0-9, a0+7, g0-8, g0+8))
    print("        " + " ".join("%6d" % g for g in G))
    for a in A:
        mark = " <hi" if a==539 else (" <lo" if a==536 else "")
        print("ray%4d " % a + " ".join((fmt % F[a,g]) if np.isfinite(F[a,g]) else "     ." for g in G) + mark)
show(Gd, "DEALIASED v_r [m/s]")
show(Gr, "RAW v_r [m/s]")
show(Gk, "fold multiple k", "%6.0f")
show(D/1000.0, "distance from assumed Moore point [km]", "%6.2f")

print("\n=== COUPLET SYMMETRY / AMBIENT CHECK ===")
sel = np.isfinite(Gd) & (D<=3000)
print("median dealiased v_r within 3 km of Moore     = %+.2f m/s (n=%d)" % (np.median(Gd[sel]), sel.sum()))
sel0 = sel & (Gk==0)
print("median over UNFOLDED-ONLY (k=0) gates         = %+.2f m/s (n=%d)" % (np.median(Gd[sel0]), sel0.sum()))
print("couplet midpoint (v_hi+v_lo)/2                = %+.2f m/s" % (0.5*(39.24-39.74)))
print("  -> a symmetric couplet about the local ambient supports the unfold being self-consistent")

print("\n=== RESIDUAL LARGE JUMPS IN THE DEALIASED FIELD: WHERE ARE THEY? ===")
for axis, nm, dd in ((1,'range',None),(0,'azimuth',None)):
    dv = np.diff(Gd, axis=axis); f = np.isfinite(dv)
    Dm = (D[:, :-1] if axis==1 else D[:-1, :])
    big = f & (np.abs(dv) > 30)
    print(" %s: %d jumps >30 m/s ; distance from Moore (km): %s"
          % (nm, big.sum(), np.round(np.sort(Dm[big])/1e3, 1)))
