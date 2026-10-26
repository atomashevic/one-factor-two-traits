import numpy as np
exec(open("fold_defs.py").read())
mA,mB,lA,lB=2,3,.85,.75
prev=None
for phi in np.round(np.arange(.5,.96,.025),3):
    R=popg([lA]*mA,[lB]*mB,phi); res=locmins(R,30) if False else None
    b=full_best(R,25); p=mA+mB; l=np.abs(b.x[:p])
    M=np.array([[1,phi],[phi,1]]); w=np.linalg.solve(M,[l[:mA].mean()/lA,l[mA:].mean()/lB]); w=w/np.abs(w).sum()
    print(phi, "loadings A %.4f B %.4f"%(l[:mA].mean(),l[mA:].mean()), "share A %.4f"%w[0], "F %.6f"%b.fun)
