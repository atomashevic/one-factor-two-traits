import numpy as np
from scipy.optimize import minimize
rng=np.random.default_rng(1)
def pop(m,lam,phi):
    p=2*m; L=np.zeros((p,2)); L[:m,0]=lam; L[m:,1]=lam
    S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1); return S
def F(theta,S):
    p=S.shape[0]; l=theta[:p]; psi=np.exp(theta[p:])
    Sig=np.outer(l,l)+np.diag(psi); A=np.linalg.solve(Sig,S)
    return np.trace(A)-np.linalg.slogdet(A)[1]-p
def grad(theta,S):
    p=S.shape[0]; l=theta[:p]; psi=np.exp(theta[p:])
    Sig=np.outer(l,l)+np.diag(psi); Si=np.linalg.inv(Sig); D=Si-Si@S@Si
    return np.r_[2*D@l, np.diag(D)*psi]
def full_best(S,starts=40):
    p=S.shape[0]; best=None
    for s in range(starts):
        x0=np.r_[rng.uniform(-0.1,1,p), np.log(rng.uniform(.2,1,p))]
        r=minimize(F,x0,args=(S,),jac=grad,method='L-BFGS-B',options={'maxiter':5000,'ftol':1e-15,'gtol':1e-12})
        if best is None or r.fun<best.fun: best=r
    return best
# 1) reproduce reviewer's lambda .8, m=4
for phi in (0,.3,.6):
    S=pop(4,.8,phi); b=full_best(S); p=8; l=np.abs(b.x[:p])
    FB=-np.linalg.slogdet(S)[1]
    print(f"lam .8 m=4 phi={phi}: block loadings {l[:4].mean():.3f}/{l[4:].mean():.3f} (within-block sd {l[:4].std():.1e},{l[4:].std():.1e}), F={b.fun:.4f}, CFI={1-b.fun/FB:.3f}")
# 2) symmetric (PCA-type, equal psi) solution F at phi .3 vs best
m=4;lam=.8;phi=.3;a=1-lam**2;c=lam**2/a;km=m*c*(1-phi);pp=8
Fsym=(pp-1)*np.log1p(km/(pp-1))-np.log1p(km); print(f"  symmetric-solution F at phi .3 = {Fsym:.3f}  (phi* = {(1+(m-1)*lam**2)/(m*(2-lam**2)):.3f})")
# 3) grid: is the full-parameter optimum always block-uniform, and does it match symmetric formula above phi*?
worst_sd=0; worst_gap=0; n=0
for m in (2,3,4):
    for lam in (.4,.6,.8,.9):
        ps=(1+(m-1)*lam**2)/(m*(2-lam**2))
        for phi in np.linspace(0,.95,12):
            S=pop(m,lam,phi); b=full_best(S,25); p=2*m; l=np.abs(b.x[:p]); psi=np.exp(b.x[p:])
            sd=max(l[:m].std(),l[m:].std(),psi[:m].std(),psi[m:].std()); worst_sd=max(worst_sd,sd)
            if phi>=ps+1e-9:
                a=1-lam**2;c=lam**2/a;km=m*c*(1-phi)
                Fs=(p-1)*np.log1p(km/(p-1))-np.log1p(km); worst_gap=max(worst_gap,abs(Fs-b.fun))
            n+=1
print(f"grid of {n} populations (m 2..4, lam .4..9, phi 0..95; 25 random starts each, full 2p parameters):")
print(f"  max within-block spread of fitted loadings/uniquenesses: {worst_sd:.1e}")
print(f"  max |F_best - F_symmetric| above phi*: {worst_gap:.1e}")
