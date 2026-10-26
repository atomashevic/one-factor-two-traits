import numpy as np
exec(open("fold_defs.py").read())
for lA,lB in [(.68,.72),(.6,.8),(.69,.71)]:
    bad=[]
    for phi in np.arange(0.01,0.96,0.02):
        R=popg([lA]*4,[lB]*4,phi); b=full_best(R,12); l=np.abs(b.x[:8])
        M=np.array([[1,phi],[phi,1]]); w=np.linalg.solve(M,[l[:4].mean()/lA,l[4:].mean()/lB])
        if not (w[1]>w[0]): bad.append(round(phi,2))
    print(lA,lB,"phis where global min does not favour stronger block:",bad)
