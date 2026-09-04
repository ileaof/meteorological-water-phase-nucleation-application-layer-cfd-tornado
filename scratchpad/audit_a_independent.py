"""(a) INDEPENDENT dealiasing (my own BFS continuity unfold) vs Py-ART's, gate by gate."""
import numpy as np
from collections import deque
P = np.load('scratchpad/_polar.npz')
Gr, Gd, GX, GY, Gk = P['Gr'], P['Gd'], P['GX'], P['GY'], P['Gk']
NY = float(P['NY']); TWO = 2*NY; mx, my = P['mxy']
D = np.hypot(GX-mx, GY-my)
nA, nR = Gr.shape
fin = np.isfinite(Gr)

# --- seed: gates that are locally smooth AND well inside the Nyquist interval -------------
def nb(a, g):
    return ((a-1) % nA, g), ((a+1) % nA, g), (a, g-1), (a, g+1)
smooth = np.zeros_like(fin)
for a in range(nA):
    for g in range(1, nR-1):
        if not fin[a, g] or abs(Gr[a, g]) > 0.55*NY: continue
        vs = [Gr[p] for p in nb(a, g) if 0 <= p[1] < nR and fin[p]]
        if len(vs) == 4 and max(abs(np.array(vs) - Gr[a, g])) < 6.0:
            smooth[a, g] = True
print("seed (smooth, |v|<0.55*Nyq, all 4 nbrs within 6 m/s): %d gates (%.1f%% of data)"
      % (smooth.sum(), 100*smooth.sum()/fin.sum()))

K = np.full(Gr.shape, np.nan)
K[smooth] = 0.0
q = deque(np.argwhere(smooth).tolist())
V = np.where(smooth, Gr, np.nan)           # current unfolded estimate
KS = np.arange(-3, 4)
n_assigned = int(smooth.sum())
while q:
    a, g = q.popleft()
    for (pa, pg) in nb(a, g):
        if not (0 <= pg < nR) or not fin[pa, pg] or np.isfinite(K[pa, pg]):
            continue
        # reference = mean of already-assigned neighbours of the candidate
        refs = [V[p] for p in nb(pa, pg) if 0 <= p[1] < nR and np.isfinite(V[p])]
        ref = float(np.mean(refs))
        cand = Gr[pa, pg] + KS*TWO
        kbest = KS[int(np.argmin(np.abs(cand - ref)))]
        K[pa, pg] = kbest
        V[pa, pg] = Gr[pa, pg] + kbest*TWO
        q.append([pa, pg])
        n_assigned += 1
print("assigned by BFS: %d of %d finite gates (%.2f%% unreached)"
      % (n_assigned, fin.sum(), 100*(1 - n_assigned/fin.sum())))

both = fin & np.isfinite(K) & np.isfinite(Gk)
agree = (K[both] == Gk[both])
print("\n=== INDEPENDENT vs PY-ART fold multiple ===")
print("agreement over %d co-assigned gates: %.4f%%  (disagree %d)"
      % (both.sum(), 100*agree.mean(), (~agree).sum()))
dis = both & (K != Gk)
if dis.any():
    print("disagreement k-differences:", np.unique((K-Gk)[dis]).astype(int))
    print("their distance from Moore (km):", np.round(np.sort(D[dis])/1e3, 1)[:40])
    print("their distance from radar (km):", np.round(np.sort(np.hypot(GX[dis],GY[dis]))/1e3, 1)[:40])

# --- couplet from MY independently unfolded field -----------------------------------------
print("\n=== COUPLET FROM THE INDEPENDENT UNFOLD (Moore-constrained) ===")
Vi = np.where(np.isfinite(K), Gr + K*TWO, np.nan)
for R in (2e3, 3e3, 5e3, 8e3):
    s = np.isfinite(Vi) & (D <= R)
    idx = np.argwhere(s); vals = Vi[s]
    ih = tuple(idx[np.argmax(vals)]); il = tuple(idx[np.argmin(vals)])
    sep = np.hypot(GX[ih]-GX[il], GY[ih]-GY[il]); dv = Vi[ih]-Vi[il]
    print("  R=%.0f km : dV=%6.1f  V_rot=%6.2f  sep=%6.0f m  in/out %+.1f/%+.1f  (k %+d/%+d)"
          % (R/1e3, dv, 0.5*dv, sep, Vi[il], Vi[ih], K[il], K[ih]))

# --- residual-jump enrichment near Moore --------------------------------------------------
print("\n=== RESIDUAL >30 m/s JUMPS: ENRICHMENT NEAR MOORE (Py-ART field) ===")
tot_pairs = 0; tot_big = 0; near_pairs = 0; near_big = 0
for axis in (0, 1):
    dv = np.diff(Gd, axis=axis); f = np.isfinite(dv)
    Dm = D[:, :-1] if axis == 1 else D[:-1, :]
    tot_pairs += f.sum(); tot_big += (f & (np.abs(dv) > 30)).sum()
    n = f & (Dm <= 2000); near_pairs += n.sum(); near_big += (n & (np.abs(dv) > 30)).sum()
print("  whole sweep : %6d big jumps / %7d pairs = %.4f%%" % (tot_big, tot_pairs, 100*tot_big/tot_pairs))
print("  <=2 km Moore: %6d big jumps / %7d pairs = %.4f%%" % (near_big, near_pairs, 100*near_big/near_pairs))
print("  ENRICHMENT  : %.0fx" % ((near_big/near_pairs)/(tot_big/tot_pairs)))

# ---- focused: what did the independent unfold do AT the couplet gates? ------------------
print("\n=== INDEPENDENT vs PY-ART, GATE BY GATE IN THE COUPLET BLOCK ===")
print(" ray gate  raw     pyart_k  pyart_v   mine_k   mine_v")
for a in range(534, 543):
    for g in range(70, 76):
        if not fin[a, g]: continue
        if Gk[a, g] == 0 and (not np.isfinite(K[a, g]) or K[a, g] == 0): continue
        print("  %3d %4d %+6.1f  %+6.0f %+8.2f  %+6s %+8s"
              % (a, g, Gr[a, g], Gk[a, g], Gd[a, g],
                 ("%.0f" % K[a, g]) if np.isfinite(K[a, g]) else "unset",
                 ("%.2f" % (Gr[a, g] + K[a, g]*TWO)) if np.isfinite(K[a, g]) else "-"))
print("\n disagreements INSIDE 2 km of Moore: %d" % ((both & (K != Gk) & (D <= 2000)).sum()))
for a, g in np.argwhere(both & (K != Gk) & (D <= 2000)):
    print("   ray%4d gate%4d raw %+6.1f : pyart k=%+d -> %+7.2f | independent k=%+d -> %+7.2f"
          % (a, g, Gr[a,g], Gk[a,g], Gd[a,g], K[a,g], Gr[a,g]+K[a,g]*TWO))
