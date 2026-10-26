import numpy as np
from scipy.optimize import brentq
def fit(m,lam,phi):
    L=lam*lam; a=1-L; n=m-1; N=2*m-1; b=1+n*L; q=m*L*phi/b; qs=L/(2-L); F0=-n*np.log(a)-np.log(b); p=2*m
    ps=b/(m*(2-L))
    if phi<ps:
        F=F0+n*np.log1p(-q*q); S=L+(2-L)*q*q; al=(np.sqrt(S+2*q)+np.sqrt(max(S-2*q,0)))/2; be=q/al if al>0 else 0
    else:
        km=m*L*(1-phi)/a; F=N*np.log1p(km/N)-np.log1p(km); al=be=np.sqrt(L*(n+m*phi)/N)
    FB=2*F0-np.log1p(-q*q); dfM=p*(p-3)//2; dfB=p*(p-1)//2
    srmr=np.sqrt((m*(m-1)/2*((L-al**2)**2+(L-be**2)**2)+m*m*(L*phi-al*be)**2)/(m*(2*m+1)))
    return dict(cfi=1-F/FB, tli=1-(F/dfM)/(FB/dfB), rmsea=np.sqrt(F/dfM), srmr=srmr, ps=ps)
def crit(m,lam,key,cut,up):
    g=lambda f:(fit(m,lam,f)[key]-cut)*(1 if up else -1)
    return 0.0 if g(0.0)>=0 else brentq(g,0.0,0.99999)
for lam in (.5,.7,.85):
    L=lam*lam
    for m in (2,4,8):
        r=[crit(m,lam,'cfi',.95,True),crit(m,lam,'tli',.95,True),crit(m,lam,'rmsea',.06,False),crit(m,lam,'srmr',.08,False)]
        print(f"lam {lam} items {2*m:2d}: CFI {r[0]:.2f} TLI {r[1]:.2f} RMSEA {r[2]:.2f} SRMR {r[3]:.2f} phi* {fit(m,lam,0)['ps']:.2f}")
    cinf=lambda f:1-np.log1p(L*(1-f)/(2*(1-L)))/(-np.log(1-L))
    print(f"lam {lam} many: CFI=TLI {brentq(lambda f:cinf(f)-.95,L/(2-L),.99999):.2f}  RMSEA 0  SRMR {max(0,1-.16/L):.2f}  phi* {L/(2-L):.2f}")
# sanity: SRMR closed form vs lavaan values earlier (m=4, lam .7, phi .6 -> .08554; m=2 lam .5 phi=0 -> .079)
print(round(fit(4,.7,.6)['srmr'],5), round(fit(2,.5,0)['srmr'],4), round(fit(3,.5,0)['srmr'],4), round(fit(4,.5,0)['srmr'],4))
