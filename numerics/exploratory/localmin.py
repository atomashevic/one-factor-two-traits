import numpy as np
from scipy.optimize import minimize
exec(open('fullparam.py').read().split("# 1) reproduce")[0])
rng=np.random.default_rng(5)
def popg(lA,lB,phi):
    lA=np.asarray(lA,float); lB=np.asarray(lB,float); mA,mB=len(lA),len(lB); p=mA+mB
    L=np.zeros((p,2)); L[:mA,0]=lA; L[mA:,1]=lB
    S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1); return S
def locmins(R,n=60):
    out=[]
    for s in range(n):
        x0=np.r_[rng.uniform(-0.1,1,8), np.log(rng.uniform(.2,1,8))]
        r=minimize(F,x0,args=(R,),jac=grad,method='L-BFGS-B',options={'maxiter':20000,'ftol':1e-16,'gtol':1e-12})
        l=np.abs(r.x[:8]); key=(round(l[:4].mean(),3),round(l[4:].mean(),3))
        # check hessian PD
        h=1e-6; H=np.array([(grad(r.x+h*e,R)-grad(r.x-h*e,R))/(2*h) for e in np.eye(16)]); H=(H+H.T)/2
        ev=np.linalg.eigvalsh(H).min()
        if not any(abs(k[0]-key[0])<2e-3 and abs(k[1]-key[1])<2e-3 for k,_,_ in out): out.append((key,r.fun,ev))
    return out
for lab,lA,lB in [("mild",.68,.72),("strong",.6,.8),("eps.01",.69,.71)]:
    for phi in (.05,.2,.3,.38):
        res=locmins(popg([lA]*4,[lB]*4,phi))
        print(lab,phi,[(k,round(f,6),f"{ev:.1e}") for k,f,ev in sorted(res,key=lambda t:t[1])])
