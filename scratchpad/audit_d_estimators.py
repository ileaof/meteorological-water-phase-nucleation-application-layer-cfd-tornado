import numpy as np
P=np.load('scratchpad/_polar.npz')
Gr,Gd,GX,GY,Gk=P['Gr'],P['Gd'],P['GX'],P['GY'],P['Gk']
NY=float(P['NY']);TWO=2*NY;mx,my=P['mxy']
D=np.hypot(GX-mx,GY-my);fin=np.isfinite(Gd)
s=fin&(D<=3000); V=Gd[s]; XY=np.stack([GX[s],GY[s]],1)
hi=np.sort(V)[::-1]; lo=np.sort(V)
print("=== ESTIMATOR SPREAD ON THE SAME (Py-ART) FIELD, 3 km search ===")
print(" top-5 outbound:", np.round(hi[:5],2), "  bottom-5 inbound:", np.round(lo[:5],2))
for n in (1,2,3,4,5):
    print("  mean of top-%d / bottom-%d  -> V_rot = %6.2f" % (n,n,0.5*(hi[:n].mean()-lo[:n].mean())))
print("  2nd-ranked pair only        -> V_rot = %6.2f" % (0.5*(hi[1]-lo[1])))
print("  median-3x3 smoothed (earlier) -> V_rot = 33.49")
print()
print("=== DEALIASING-BRANCH SCENARIOS FOR THE COUPLET ===")
print("  Py-ART as-is                       : in -39.74 out +39.24 -> V_rot 39.49  sep 584 m")
print("  + core gates unfolded (indep. BFS) : in -51.24 out +39.24 -> V_rot %5.2f  sep %.0f m"
      % (0.5*(39.24-(-51.24)), 2*0.5*np.pi/180*20250))
print("  core gates treated as bad/ignored  : in -39.74 out +39.24 -> V_rot 39.49  sep 584 m")
print("  independent BFS as-is              : in -51.24 out +31.24 -> V_rot 41.24  sep 753 m")
print()
print("=== IMPLIED VORTICITY / SHEAR (definition matters) ===")
for dv,sep,lbl in ((79.0,584.,'corrected'),(52.0,253.,'old (aliased, wrong feature)')):
    print("  %-30s dV/sep = %.3f s^-1 (= Omega) ; zeta = 2*dV/sep = %.3f s^-1"
          % (lbl, dv/sep, 2*dv/sep))
print("  docs/TORNADOGENESIS_FINDINGS.md calls dV/sep 'vorticity'; for solid-body rotation it is")
print("  Omega = zeta/2, so the observed-vs-model vorticity table compares unlike quantities.")
