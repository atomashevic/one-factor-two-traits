import numpy as np, itertools
from scipy.stats import norm
from scipy.optimize import minimize
pats=np.array(list(itertools.product([0,1],repeat=4)))
def setup(n):
    x,w=np.polynomial.hermite_e.hermegauss(n); return x,w/w.sum()
def true_dist(lam,phi,b,x,w):
    a=lam/np.sqrt(1-lam**2); L=np.linalg.cholesky(np.array([[1,phi],[phi,1]])+1e-15*np.eye(2))
    T1,T2=np.meshgrid(x,x,indexing='ij'); W=np.outer(w,w); th1=L[0,0]*T1; th2=L[1,0]*T1+L[1,1]*T2
    P=np.stack([norm.cdf(a*th1-b[0]),norm.cdf(a*th1-b[1]),norm.cdf(a*th2-b[2]),norm.cdf(a*th2-b[3])],-1)
    return np.array([np.sum(W*np.prod(np.where(xx==1,P,1-P),-1)) for xx in pats])
def fit(px,x,w,starts=20):
    def q(th):
        P=norm.cdf(np.outer(x,th[:4])-th[4:]); lik=np.prod(np.where(pats[:,None,:]==1,P[None],1-P[None]),-1); return (lik*w).sum(1)
    obj=lambda th: -np.sum(px*np.log(q(th)))
    best=None
    for s in range(starts):
        rs=np.random.RandomState(s); r=minimize(obj,np.r_[rs.uniform(.05,2,4),rs.normal(0,.3,4)],method='L-BFGS-B',bounds=[(0,8)]*4+[(-4,4)]*4)
        if best is None or r.fun<best.fun: best=r
    return best, np.sum(px*np.log2(px/q(best.x)))
b=np.zeros(4)
for n in (41,81,121):
    x,w=setup(n); px=true_dist(.7,0.0,b,x,w)
    m12=np.array([[px[(pats[:,0]==i)&(pats[:,1]==j)].sum() for j in (0,1)] for i in (0,1)])
    mi=np.sum(m12*np.log2(m12/np.outer(m12.sum(1),m12.sum(0))))
    best,kl=fit(px,x,w)
    print(f"nodes={n}: phi=0 KL(one trait)={kl:.6f} bits; I(X1;X2)={mi:.6f} bits; slopes {np.round(best.x[:4],3)}")
for n in (61,121):
    x,w=setup(n)
    for phi in (.2,.3,.35,.4):
        px=true_dist(.7,phi,b,x,w); best,kl=fit(px,x,w); print(f"nodes={n} phi={phi}: slopes {np.round(best.x[:4],3)}, KL {kl:.6f}")
