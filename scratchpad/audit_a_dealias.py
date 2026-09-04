"""(a) VERIFY THE DEALIASING - failure-mode audit on the saved KTLX field."""
import numpy as np
d = np.load('outputs/nexrad_moore/ktlx_velocity_dealiased.npz', allow_pickle=True)
x, y, z = d['x_m'], d['y_m'], d['z_m']
vr = d['vr_raw'].astype(float); vd = d['vr_dealiased'].astype(float)
NY = float(d['nyquist_m_s']); TWO = 2*NY
mx, my = d['moore_xy']
r = np.hypot(x, y); az = np.degrees(np.arctan2(x, y)) % 360.0

# ---------- 1. fold multiples --------------------------------------------------------
k = (vd - vr)/TWO
kr = np.round(k)
print("=== 1. FOLD MULTIPLE k = (v_dealiased - v_raw)/(2*Nyquist) ===")
print("max |k - round(k)| = %.3e   (must be ~0 if only integer folds were added)" % np.abs(k-kr).max())
u, c = np.unique(kr.astype(int), return_counts=True)
for uu, cc in zip(u, c):
    print("   k=%+d : %7d gates (%6.3f%%)" % (uu, cc, 100*cc/kr.size))
print("gates changed: %d (%.3f%%)" % ((kr!=0).sum(), 100*(kr!=0).mean()))
# where in raw space are folded gates?
for uu in u:
    if uu == 0: continue
    m = kr.astype(int)==uu
    print("   k=%+d raw v range %+.1f..%+.1f (mean %+.1f), dealiased %+.1f..%+.1f, "
          "range from radar %.1f..%.1f km" % (uu, vr[m].min(), vr[m].max(), vr[m].mean(),
          vd[m].min(), vd[m].max(), r[m].min()/1e3, r[m].max()/1e3))

# ---------- 2. rebuild the polar (ray, gate) lattice ---------------------------------
uaz, iaz = np.unique(np.round(az, 6), return_inverse=True)
r0 = r.min(); ir = np.round((r - r0)/249.99).astype(int)
nA, nR = uaz.size, ir.max()+1
print("\n=== 2. POLAR LATTICE %d rays x %d gates (%.1f%% populated) ===" %
      (nA, nR, 100*x.size/(nA*nR)))
def grid(v):
    G = np.full((nA, nR), np.nan); G[iaz, ir] = v; return G
Gr, Gd, GX, GY, GZ = grid(vr), grid(vd), grid(x), grid(y), grid(z)

def jumps(G, axis, label):
    dv = np.diff(G, axis=axis)
    fin = np.isfinite(dv)
    a = np.abs(dv[fin])
    out = dict(n=fin.sum(), p50=np.percentile(a,50), p99=np.percentile(a,99),
               p999=np.percentile(a,99.9), mx=a.max(),
               n20=(a>20).sum(), n30=(a>30).sum(), n40=(a>40).sum(),
               nfold=((a>0.75*TWO)&(a<1.25*TWO)).sum(),
               n2fold=(a>1.75*TWO).sum())
    print("  %-28s n=%7d  |dv| p50=%5.2f p99=%6.2f p99.9=%6.2f max=%6.2f | "
          ">20:%5d >30:%4d >40:%4d  ~2Nyq(39-65):%4d  >91:%3d"
          % (label, out['n'], out['p50'], out['p99'], out['p999'], out['mx'],
             out['n20'], out['n30'], out['n40'], out['nfold'], out['n2fold']))
    return dv, fin

print("\n=== 3. ADJACENT-GATE JUMPS (spatial coherence) ===")
print(" 2*Nyquist = %.2f m/s ; a residual/added fold shows as a ~52 m/s neighbour jump" % TWO)
for G, nm in ((Gr,'RAW'),(Gd,'DEALIASED')):
    jumps(G, 1, nm+' along-range (gate to gate)')
    jumps(G, 0, nm+' along-azimuth (ray to ray)')

# ---------- 4. implausible added folds: fold seams that are NOT at a real discontinuity
print("\n=== 4. FOLD-SEAM AUDIT (did dealiasing create a jump where raw had none?) ===")
Gk = grid(kr)
for axis, nm in ((1,'along-range'),(0,'along-azimuth')):
    dk = np.diff(Gk, axis=axis); dr_ = np.diff(Gr, axis=axis); dd = np.diff(Gd, axis=axis)
    seam = np.isfinite(dk) & (np.abs(dk) >= 1)
    ns = seam.sum()
    if ns:
        # at a CORRECT seam the raw jump is ~ -k*2Nyq, so the dealiased jump becomes small
        print("  %-14s seams=%6d  |raw jump| at seam: p50=%5.1f  "
              "|dealiased jump| at seam: p50=%5.2f p95=%6.2f max=%6.2f  worse-than-raw: %d (%.1f%%)"
              % (nm, ns, np.median(np.abs(dr_[seam])), np.median(np.abs(dd[seam])),
                 np.percentile(np.abs(dd[seam]),95), np.abs(dd[seam]).max(),
                 (np.abs(dd[seam])>np.abs(dr_[seam])).sum(),
                 100*(np.abs(dd[seam])>np.abs(dr_[seam])).mean()))

# ---------- 5. global coherence score ------------------------------------------------
print("\n=== 5. TOTAL-VARIATION (smoothness) SCORE, whole sweep ===")
for G, nm in ((Gr,'RAW'),(Gd,'DEALIASED')):
    tv = 0.0; n = 0
    for axis in (0,1):
        dv = np.diff(G, axis=axis); f = np.isfinite(dv)
        tv += np.abs(dv[f]).sum(); n += f.sum()
    print("  %-10s mean |grad| over %d neighbour pairs = %.4f m/s" % (nm, n, tv/n))
np.savez(r'scratchpad/_polar.npz', Gr=Gr, Gd=Gd, GX=GX, GY=GY, GZ=GZ, Gk=Gk,
         uaz=uaz, r0=r0, NY=NY, mxy=np.array([mx,my]))
