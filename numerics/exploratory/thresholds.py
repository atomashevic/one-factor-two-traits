import numpy as np
from scipy.optimize import minimize, brentq
exec(open('global1f.py').read().split("lam=.7\nfor m in")[0])
# ULS threshold check: block-uniform MINRES objective
def uls(m,lam,phi):
    C=m*(m-1)/2
    f=lambda x: C*((lam**2-x[0]**2)**2+(lam**2-x[1]**2)**2)+m*m*(lam**2*phi-x[0]*x[1])**2
    best=min((minimize(f,s,method='BFGS') for s in ([lam,.01],[.4,.6],[.01,lam])), key=lambda r:r.fun)
    return np.abs(best.x)
print("ULS symmetry-breaking threshold predicted 1/p:")
for m in (2,3,5):
    for lam in (.5,.8):
        lo=1/(2*m)
        a=uls(m,lam,lo-.02); b=uls(m,lam,lo+.02)
        print(f" m={m} lam={lam}: phi=1/p-.02 -> {np.round(a,3)} ; phi=1/p+.02 -> {np.round(b,3)}")
# fit index thresholds (global ML)
def idx(m,lam,phi):
    p=2*m; F=F_global(m,lam,phi) if phi<phistar(m,lam) else F_sym(m,lam,phi); FB=F0(m,lam,phi)
    dfM=p*(p+1)//2-2*p; dfB=p*(p-1)//2
    cfi=1-F/FB; tli=1-(F/dfM)/(FB/dfB); rmsea=np.sqrt(F/dfM)
    srmr=lam**2*(1-phi)*np.sqrt(m*(m-1)/(4*m*m-1))
    return cfi,tli,rmsea,srmr
def crit(m,lam,k,cut,up):
    g=lambda f:(idx(m,lam,f)[k]-cut)*(1 if up else -1)
    return brentq(g,0.0,0.9999) if g(0.0)<0 else 0.0
print("\nphi above which a merged one-factor model passes each cutoff (population, Hu & Bentler cutoffs)")
for lam in (.5,.7,.85):
    for m in (2,3,4,6,8):
        r=[crit(m,lam,0,.95,True),crit(m,lam,1,.95,True),crit(m,lam,2,.06,False),crit(m,lam,3,.08,False)]
        print(f" lam={lam} m={m} (p={2*m:2d}): CFI .95 -> {r[0]:.2f}  TLI .95 -> {r[1]:.2f}  RMSEA .06 -> {r[2]:.2f}  SRMR .08 -> {r[3]:.2f}   [phi* ML {phistar(m,lam):.2f}]")
    cinf=lambda f:1-np.log1p(lam**2*(1-f)/(2*(1-lam**2)))/(-np.log(1-lam**2))
    print(f" lam={lam} m=inf: CFI=TLI -> {brentq(lambda f:cinf(f)-.95,0.3,.9999):.2f}  RMSEA -> 0 (always passes)  SRMR -> {max(0,1-.16/lam**2):.2f}")
