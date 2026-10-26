import numpy as np
from scipy.optimize import minimize
def F_block(params, m, lam, phi):
    lA,lB,lpA,lpB = params; psA,psB=np.exp(lpA),np.exp(lpB); a=1-lam**2
    S2=np.array([[a+m*lam**2, m*lam**2*phi],[m*lam**2*phi, a+m*lam**2]])
    T2=np.array([[psA+m*lA**2, m*lA*lB],[m*lA*lB, psB+m*lB**2]])
    A=np.linalg.solve(T2,S2); f2=np.trace(A)-np.linalg.slogdet(A)[1]-2
    g=lambda s: a/s-np.log(a/s)-1
    return (m-1)*(g(psA)+g(psB))+f2
def F_sym(m,lam,phi):
    a=1-lam**2; c=lam**2/a; km=m*c*(1-phi); p=2*m
    return (p-1)*np.log1p(km/(p-1))-np.log1p(km)
def F_global(m,lam,phi):
    best=np.inf
    for st in [(lam,lam*0.9),(lam,0.05),(0.05,lam),(lam,-lam*.5),(.3,.6),(.9,.1)]:
        for ps in [(1-lam**2,1-lam**2),(1-lam**2,1.0),(.5,.9)]:
            r=minimize(F_block,[st[0],st[1],np.log(ps[0]),np.log(ps[1])],args=(m,lam,phi),method='Nelder-Mead',options={'xatol':1e-12,'fatol':1e-14,'maxiter':20000,'maxfev':40000})
            r=minimize(F_block,r.x,args=(m,lam,phi),method='BFGS',options={'gtol':1e-12})
            best=min(best,r.fun)
    return best
def F0(m,lam,phi):
    a=1-lam**2; c=lam**2/a; p=2*m; return -(p*np.log(a)+np.log1p(m*c*(1+phi))+np.log1p(m*c*(1-phi)))
phistar=lambda m,lam:(1+(m-1)*lam**2)/(m*(2-lam**2))
lam=.7
for m in (2,4,8,10**6):
    ps=phistar(m,lam)
    print(f"m={m}: phi*={ps:.4f}")
    for phi in (0.0,0.2,ps-0.02,ps+0.02,0.6,0.8):
        fg=F_global(m,lam,phi); fs=F_sym(m,lam,phi)
        print(f"   phi={phi:.3f}  F_sym={fs:.6f}  F_global={fg:.6f}  gap={fs-fg:.2e}  CFI_sym={1-fs/F0(m,lam,phi):.4f} CFI_global={1-fg/F0(m,lam,phi):.4f}")
# exact limit formula above phi*_inf
cinf=lambda lam,phi:1-np.log1p(lam**2*(1-phi)/(2*(1-lam**2)))/(-np.log(1-lam**2))
print("limit check phi=.5:", cinf(.7,.5), " numeric m=1e6:", 1-F_global(10**6,.7,.5)/F0(10**6,.7,.5))
