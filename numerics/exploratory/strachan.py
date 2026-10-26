import numpy as np
from scipy.optimize import minimize
exec(open('fullparam.py').read().split("# 1) reproduce")[0])
rng=np.random.default_rng(22)  # after exec: fullparam.py resets rng to seed 1
def design(J):
    h=J//2
    a1=np.r_[rng.uniform(1.2,1.6,h), rng.uniform(.2,.6,h)]
    a2=np.r_[rng.uniform(.2,.6,h), rng.uniform(1.2,1.6,h)]
    return np.c_[a1,a2]/1.702          # logistic -> probit slopes
def phi_eff(p,s,rho): return (2*p*s+rho*(p*p+s*s))/(p*p+s*s+2*p*s*rho)
p_,s_=1.4/1.702,.4/1.702
lam2=lambda rho:(p_**2+s_**2+2*p_*s_*rho)/(1+p_**2+s_**2+2*p_*s_*rho)
for J in (40,):
    A=design(J); m=J//2
    print(f"Strachan et al. slope design, J={J} items, Gaussian latent-response analogue, one-factor ML (25 starts)")
    for rho in (.8,.4,0,-.4,-.8):
        Phi=np.array([[1,rho],[rho,1]]); C=A@Phi@A.T+np.eye(J); d=np.sqrt(np.diag(C)); R=C/np.outer(d,d)
        b=full_best(R,25); l=np.abs(b.x[:J]); FB=-np.linalg.slogdet(R)[1]
        fe=phi_eff(p_,s_,rho); L2=lam2(rho); ps=(1+(m-1)*L2)/(m*(2-L2))
        print(f" rho={rho:+.1f}: effective phi={fe:+.2f}, effective loading {np.sqrt(L2):.2f}, ML phi*={ps:.2f} -> block mean loadings {l[:m].mean():.3f} / {l[m:].mean():.3f}, CFI {1-b.fun/FB:.3f}")
print("signed loadings check:")
rng=np.random.default_rng(22); A=design(40); m=20
for rho in (0,-.4,-.8):
    Phi=np.array([[1,rho],[rho,1]]); C=A@Phi@A.T+np.eye(40); d=np.sqrt(np.diag(C)); R=C/np.outer(d,d)
    b=full_best(R,25); l=b.x[:40]; s=np.sign(l[:m].mean()) or 1; l=l*s
    print(f" rho={rho:+.1f}: block A mean {l[:m].mean():+.3f}, block B mean {l[m:].mean():+.3f}; mean observed corr A-B {R[:m,m:].mean():+.3f}, within A {R[:m,:m][~np.eye(m,dtype=bool)].mean():.3f}")
