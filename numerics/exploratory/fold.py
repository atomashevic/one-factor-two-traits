import numpy as np
from scipy.optimize import minimize
exec(open('localmin.py').read().split("for lab,lA,lB in")[0])
def track(lA,lB,phis):
    # start on the branch favouring block A (the weaker one)
    R=popg([lA]*4,[lB]*4,phis[0])
    x=np.r_[[lA]*4,[0.05]*4,np.log([1-lA**2]*4),np.log([.99]*4)]
    last=None
    for phi in phis:
        R=popg([lA]*4,[lB]*4,phi)
        r=minimize(F,x,args=(R,),jac=grad,method='L-BFGS-B',options={'maxiter':20000,'ftol':1e-16,'gtol':1e-12})
        l=np.abs(r.x[:8])
        if l[:4].mean()<l[4:].mean(): return last
        x=r.x; last=phi
    return last
for lab,lA,lB in [("eps.0025",.6975,.7025),("eps.01",.69,.71),("mild .68/.72",.68,.72),("eps.04",.66,.74),(".65/.75",.65,.75),("strong .6/.8",.6,.8)]:
    phis=np.arange(0.0005,0.42,0.0005)
    f=track(lA,lB,phis)
    print(lab, "weaker-trait local minimum persists up to phi =", None if f is None else round(f,4))
