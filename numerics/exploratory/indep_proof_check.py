import numpy as np
exec(open('fullparam.py').read().split("# 1) reproduce")[0])   # my own full-parameter optimizer (written earlier, independent of theirs)
rng=np.random.default_rng(99)
def closed(m,lam,phi):
    L=lam*lam; a=1-L; n=m-1; N=2*m-1; b=1+n*L; q=m*L*phi/b; qs=L/(2-L); F0=-n*np.log(a)-np.log(b)
    ps=b/(m*(2-L))
    if phi<ps:
        S=L+(2-L)*q*q; al=(np.sqrt(S+2*q)+np.sqrt(max(S-2*q,0)))/2; be=(np.sqrt(S+2*q)-np.sqrt(max(S-2*q,0)))/2
        return F0+n*np.log1p(-q*q), al, be, ps
    km=m*L*(1-phi)/a; ls=np.sqrt(L*(n+m*phi)/N)
    return N*np.log1p(km/N)-np.log1p(km), ls, ls, ps
worst=0; worstl=0; rows=0
for _ in range(60):
    m=int(rng.integers(2,7)); lam=float(rng.uniform(.25,.95)); ps=(1+(m-1)*lam**2)/(m*(2-lam**2))
    phi=float(np.clip(ps+rng.normal(0,.15),0.005,.97)) if rng.random()<.6 else float(rng.uniform(.005,.97))
    S=pop(m,lam,phi); b=full_best(S,30); Fc,al,be,_=closed(m,lam,phi)
    l=np.abs(b.x[:2*m]); la,lb=sorted([l[:m].mean(),l[m:].mean()],reverse=True)
    worst=max(worst,abs(b.fun-Fc)); worstl=max(worstl,abs(la-al),abs(lb-be)); rows+=1
print(f"{rows} random populations (m 2..6, 30 starts, my optimizer): max |F_opt - F_closed| = {worst:.1e}, max loading error = {worstl:.1e}")
# below-threshold branch: optimizer must never beat the closed form
# C1-not-C2 check by finite differences of F_min(phi) around phi*
m,lam=3,.6; ps=(1+(m-1)*lam**2)/(m*(2-lam**2)); h=1e-3
f=lambda p: closed(m,lam,p)[0]
d1l=(f(ps-1e-7)-f(ps-1e-7-h))/h; d1r=(f(ps+1e-7+h)-f(ps+1e-7))/h
d2l=(f(ps-2*h)-2*f(ps-h)+f(ps))/h**2; d2r=(f(ps)-2*f(ps+h)+f(ps+2*h))/h**2
print(f"at phi*={ps:.4f}: one-sided slopes {d1l:.4f} / {d1r:.4f}; one-sided curvatures {d2l:.3f} / {d2r:.3f}")
# Heywood boundary derivative (4.3) vs finite difference, my own construction
for (m,lam,phi) in [(3,.7,.4),(2,.5,.2),(4,.8,.0)]:
    R=pop(m,lam,phi); p=2*m; j=0; l=R[:,j].copy(); psi=1-l*l; psi[j]=0.0
    Sig=np.outer(l,l)+np.diag(psi); Fn=lambda Sg: np.trace(np.linalg.solve(Sg,R))-np.linalg.slogdet(np.linalg.solve(Sg,R))[1]-p
    e=1e-7; E=np.zeros((p,p)); E[j,j]=e; fd=(Fn(Sig+E)-Fn(Sig))/e
    s=0
    for i in range(p):
        for k in range(p):
            if i!=k and i!=j and k!=j: s+=R[i,j]*R[k,j]*(R[i,k]-R[i,j]*R[k,j])/((1-R[i,j]**2)*(1-R[k,j]**2))
    print(f"Heywood inward derivative m={m} lam={lam} phi={phi}: finite diff {fd:.6f}, formula (4.3) {-s:.6f}")
