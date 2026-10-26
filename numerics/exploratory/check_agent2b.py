import numpy as np
from scipy.optimize import minimize
exec(open("fold_defs.py").read())
def polish(S,x0): return minimize(F,x0,args=(S,),jac=grad,method='L-BFGS-B',options={'maxiter':20000,'ftol':1e-16,'gtol':1e-13})
Rm=popg([.68]*4,[.72]*4,.05)
for x0 in (np.r_[[.68]*4,[.05]*4,np.log([1-.68**2]*4),np.log([.99]*4)], np.r_[[.05]*4,[.72]*4,np.log([.99]*4),np.log([1-.72**2]*4)]):
    r=polish(Rm,x0); l=np.abs(r.x[:8]); ps=np.exp(r.x[8:]); al,be=l[:4].mean(),l[4:].mean(); x,y=ps[:4].mean(),ps[4:].mean()
    print("solution F %.6f  v=%.6f"%(r.fun,(be/y)/(al/x)))
for lA in (.565,.575,.59):
    for phi in (.001,.005):
        res=locmins(popg([lA]*4,[.8]*4,phi),40)
        print(f"lamA={lA} phi={phi}: {len(res)} local minima",[(k,round(f,5)) for k,f,ev in res])
