import numpy as np
from collections import deque
P = np.load('scratchpad/_polar.npz')
Gr,Gd,GX,GY,Gk = P['Gr'],P['Gd'],P['GX'],P['GY'],P['Gk']
NY=float(P['NY']); TWO=2*NY; mx,my=P['mxy']
D=np.hypot(GX-mx,GY-my); nA,nR=Gr.shape
print("=== CONNECTED COMPONENTS OF UNFOLDED (k != 0) GATES ===")
lab=np.zeros(Gr.shape,int); cur=0; comps=[]
M=np.isfinite(Gk)&(Gk!=0)
for a in range(nA):
    for g in range(nR):
        if M[a,g] and lab[a,g]==0:
            cur+=1; q=deque([(a,g)]); lab[a,g]=cur; cells=[]
            while q:
                p=q.popleft(); cells.append(p)
                for pa,pg in (((p[0]-1)%nA,p[1]),((p[0]+1)%nA,p[1]),(p[0],p[1]-1),(p[0],p[1]+1)):
                    if 0<=pg<nR and M[pa,pg] and lab[pa,pg]==0 and Gk[pa,pg]==Gk[a,g]:
                        lab[pa,pg]=cur; q.append((pa,pg))
            comps.append((len(cells),int(Gk[a,g]),cells))
comps.sort(reverse=True,key=lambda c:c[0])
print(" %d components; sizes (top 15): %s" % (len(comps),[c[0] for c in comps[:15]]))
print(" singletons: %d (%.0f%% of components), gates in components of size 1-2: %d"
      % (sum(1 for c in comps if c[0]==1), 100*sum(1 for c in comps if c[0]==1)/len(comps),
         sum(c[0] for c in comps if c[0]<=2)))
print("\n top-8 components:")
for n,k,cells in comps[:8]:
    ca=np.array(cells); rr=np.hypot(GX[ca[:,0],ca[:,1]],GY[ca[:,0],ca[:,1]])/1e3
    dm=D[ca[:,0],ca[:,1]]/1e3
    print("  n=%4d k=%+d  radar-range %5.1f-%5.1f km  d_Moore %5.1f-%5.1f km  "
          "raw %+6.1f..%+6.1f -> dealiased %+6.1f..%+6.1f"
          % (n,k,rr.min(),rr.max(),dm.min(),dm.max(),
             Gr[ca[:,0],ca[:,1]].min(),Gr[ca[:,0],ca[:,1]].max(),
             Gd[ca[:,0],ca[:,1]].min(),Gd[ca[:,0],ca[:,1]].max()))

print("\n=== VAD CONSISTENCY CHECK (large-scale: no whole-region 2*Nyquist offset?) ===")
az=np.degrees(np.arctan2(GX,GY)); rr=np.hypot(GX,GY)
print(" ring_km  n   VAD fit u,v (m/s)  rms resid  max|resid|  n resid>0.7*2Nyq")
for r0 in (5,10,15,20,30,40,60,80):
    m=np.isfinite(Gd)&(np.abs(rr-r0*1000)<1000)
    if m.sum()<200: continue
    A=np.radians(az[m]); V=Gd[m]
    Mx=np.stack([np.sin(A),np.cos(A),np.ones_like(A)],1)
    coef,*_=np.linalg.lstsq(Mx,V,rcond=None)
    res=V-Mx@coef
    print("  %5d %6d   u=%+6.1f v=%+6.1f   %6.2f     %6.1f        %4d"
          % (r0,m.sum(),coef[0],coef[1],np.sqrt((res**2).mean()),np.abs(res).max(),
             (np.abs(res)>0.7*TWO).sum()))
