import numpy as np, itertools
from scipy.stats import norm
from scipy.optimize import minimize
exec(open('binary4.py').read().split("b=np.array")[0])
b=np.array([0,0,0,0.])
print("4 binary items (2+2), loading .7, thresholds 0; exact KL fit of 1-trait normal ogive")
for phi in (0.0,0.2,0.4,0.5,0.6,0.8):
    px=true_dist(.7,phi,b)
    best=None
    for s in range(12):
        rs=np.random.RandomState(s); x0=np.r_[rs.uniform(.2,2,4),np.zeros(4)]
        obj=lambda th: -np.sum(px*np.log((onef(th)[0]*gh_w).sum(1)))
        r=minimize(obj,x0,method='L-BFGS-B',bounds=[(0.001,8)]*4+[(-4,4)]*4)
        if best is None or r.fun<best.fun: best=r
    a=best.x[:4]
    # symmetric constrained fit
    objs=lambda t: -np.sum(px*np.log((onef(np.r_[[t[0]]*4,[0]*4])[0]*gh_w).sum(1)))
    rs_=minimize(objs,[1.0],bounds=[(0.001,8)])
    print(f" phi={phi:.1f}: slopes {np.round(a,3)}  KLgap(sym-best)={(rs_.fun-best.fun)/np.log(2):.2e} bits")
