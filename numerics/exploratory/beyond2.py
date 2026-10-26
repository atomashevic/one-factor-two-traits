import numpy as np
from scipy.optimize import minimize, brentq
exec(open('fullparam.py').read().split("# 1) reproduce")[0])
exec(open('beyond.py').read().split("# 1) extended")[0].split("rng=")[0])
def popg(lA,lB,phi):
    lA=np.asarray(lA,float); lB=np.asarray(lB,float); mA,mB=len(lA),len(lB); p=mA+mB
    L=np.zeros((p,2)); L[:mA,0]=lA; L[mA:,1]=lB
    S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1); return S
def sym_solution(lam,phi):
    k=len(lam); R=popg(lam,lam,phi)
    def full(t): return np.r_[t[:k],t[:k],t[k:],t[k:]]
    f=lambda t: F(full(t),R)
    g=lambda t: (lambda G: np.r_[G[:k]+G[k:2*k], G[2*k:3*k]+G[3*k:]])(grad(full(t),R))
    t0=np.r_[np.array(lam)*.8, np.log(1-np.array(lam)**2*.64)]
    r=minimize(f,t0,jac=g,method='BFGS',options={'gtol':1e-13})
    return full(r.x),R
def anti_min_eig(lam,phi,h=1e-6):
    th,R=sym_solution(lam,phi); k=len(lam); n=len(th)
    H=np.zeros((n,n))
    for i in range(n):
        e=np.zeros(n); e[i]=h; H[:,i]=(grad(th+e,R)-grad(th-e,R))/(2*h)
    H=(H+H.T)/2
    # antisymmetric directions: (v,-v) in loadings and log-psi
    B=[]
    for i in range(k):
        v=np.zeros(n); v[i]=1; v[k+i]=-1; B.append(v/np.sqrt(2))
        v=np.zeros(n); v[2*k+i]=1; v[3*k+i]=-1; B.append(v/np.sqrt(2))
    B=np.array(B).T
    return np.linalg.eigvalsh(B.T@H@B).min()
def phistar_exch(L,m): return (1+(m-1)*L)/(m*(2-L))
for lam in ([.7]*4,[.6,.65,.75,.8],[.5,.6,.7,.8]):
    lam=np.array(lam); m=len(lam)
    ps=brentq(lambda f: anti_min_eig(lam,f),.05,.9,xtol=1e-6)
    c=lam**2/(1-lam**2)
    print(lam, f"exact threshold {ps:.4f}",
          f"| mean lambda {phistar_exch(lam.mean()**2,m):.4f}",
          f"| mean lambda^2 {phistar_exch(np.mean(lam**2),m):.4f}",
          f"| mean c {phistar_exch(c.mean()/(1+c.mean()),m):.4f}")
