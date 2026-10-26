import numpy as np, itertools
from scipy.stats import norm
from scipy.optimize import minimize
pats=np.array(list(itertools.product([0,1],repeat=4)))            # 16 patterns
gh_x,gh_w=np.polynomial.hermite_e.hermegauss(61); gh_w=gh_w/gh_w.sum()   # N(0,1) nodes
def true_dist(lam,phi,b):
    a=lam/np.sqrt(1-lam**2)                                       # normal-ogive slope matching loading lam
    L=np.linalg.cholesky(np.array([[1,phi],[phi,1]]))
    T1,T2=np.meshgrid(gh_x,gh_x,indexing='ij'); W=np.outer(gh_w,gh_w)
    th1=L[0,0]*T1; th2=L[1,0]*T1+L[1,1]*T2
    P=np.stack([norm.cdf(a*th1-b[0]),norm.cdf(a*th1-b[1]),norm.cdf(a*th2-b[2]),norm.cdf(a*th2-b[3])],-1)
    px=np.zeros(16)
    for k,x in enumerate(pats):
        px[k]=np.sum(W*np.prod(np.where(x==1,P,1-P),-1))
    return px
def onef(params):                                                # 1-trait normal ogive on the 1D grid
    a=params[:4]; b=params[4:]
    P=norm.cdf(np.outer(gh_x,a)-b)                                # grid x items
    lik=np.prod(np.where(pats[:,None,:]==1,P[None],1-P[None]),-1) # patterns x grid
    return lik, P
def H(p): p=p[p>0]; return -np.sum(p*np.log2(p))
def fit1(px):
    obj=lambda th: -np.sum(px*np.log((onef(th)[0]*gh_w).sum(1)))
    best=None
    for s in range(6):
        x0=np.r_[np.random.RandomState(s).uniform(.5,1.5,4), np.zeros(4)]
        r=minimize(obj,x0,method='L-BFGS-B',bounds=[(0.01,8)]*4+[(-4,4)]*4)
        if best is None or r.fun<best.fun: best=r
    return best.x
b=np.array([-.5,.5,-.3,.3])
print("4 binary items, 2 per trait (normal ogive, loading .7), 1-trait model fitted by exact KL minimisation (N -> inf)")
print(" phi |  H(X) bits  TC bits | KL_1F bits  infoCFI | decomposition: TC(X|T)  item  prior  sum-KL | N for power .8 (df=7)")
from scipy.stats import ncx2, chi2
for phi in (.3,.5,.7,.8,.9,.95):
    px=true_dist(.7,phi,b)
    marg=np.array([[px[pats[:,i]==v].sum() for v in (0,1)] for i in range(4)])
    TCx=sum(H(marg[i]) for i in range(4))-H(px)
    th=fit1(px); lik,P=onef(th); qx=(lik*gh_w).sum(1)
    KL=np.sum(px*np.log2(px/qx))
    # posterior coupling P(x,t)=p(x)q(t|x) on the grid
    post=lik*gh_w/qx[:,None]; J=px[:,None]*post                    # patterns x grid
    Pt=J.sum(0)                                                    # P(theta) on grid
    # conditional joint P(x|t), conditional marginals P(x_i|t)
    Px_t=J/np.where(Pt>0,Pt,1)[None]
    tc_cond=0; item=0
    for g in range(len(gh_x)):
        if Pt[g]<1e-300: continue
        pj=Px_t[:,g]
        mi=np.array([[pj[pats[:,i]==v].sum() for v in (0,1)] for i in range(4)])
        tc_cond+=Pt[g]*(sum(H(mi[i]) for i in range(4))-H(pj))
        qi=np.stack([1-P[g],P[g]],1)
        qc=np.clip(qi,1e-300,1); item+=Pt[g]*np.sum(np.where(mi>0, mi*np.log2(np.where(mi>0,mi,1)/qc), 0.0))
    prior=np.sum(Pt[Pt>0]*np.log2(Pt[Pt>0]/gh_w[Pt>0]))
    # sample size for 80% power of the LR (G2) test, noncentrality = 2N*KL(nats), df = 15-8 = 7
    KLn=KL*np.log(2); crit=chi2.ppf(.95,7)
    N=next(n for n in range(10,10**7,10) if ncx2.sf(crit,7,2*n*KLn)>=.8) if KLn>0 else np.inf
    print(f" {phi:.2f} | {H(px):9.4f} {TCx:8.4f} | {KL:10.6f} {1-KL/TCx:8.4f} | {tc_cond:14.6f} {item:6.6f} {prior:6.6f} {tc_cond+item+prior-KL:8.1e} | {N}")
