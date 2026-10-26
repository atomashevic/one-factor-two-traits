import numpy as np
from scipy.optimize import minimize
exec(open('global1f.py').read().split("lam=.7\nfor m in")[0])
g=lambda x: x-np.log(x)-1
def Fconc(psA,psB,m,lam,phi):
    a=1-lam**2; M=np.array([[(a+m*lam**2)/psA, m*lam**2*phi/np.sqrt(psA*psB)],[m*lam**2*phi/np.sqrt(psA*psB),(a+m*lam**2)/psB]])
    th=np.r_[[a/psA]*(m-1),[a/psB]*(m-1),np.linalg.eigvalsh(M)]
    th=np.sort(th)[::-1]
    return np.sum(g(th[1:])) + (g(th[0]) if th[0]<1 else 0)
def F2(m,lam,phi):
    best=np.inf
    for s in ([.5,.5],[.5,.99],[.99,.5],[.3,.9],[.9,.3]):
        r=minimize(lambda x:Fconc(np.exp(x[0]),np.exp(x[1]),m,lam,phi),np.log(s),method='Nelder-Mead',options={'xatol':1e-12,'fatol':1e-15,'maxiter':20000})
        best=min(best,r.fun)
    return best
for (m,lam,phi) in [(2,.7,0),(2,.7,.2),(4,.8,.3),(4,.8,.6),(3,.5,.1),(5,.9,.65)]:
    print(f"m={m} lam={lam} phi={phi}: concentrated-2D {F2(m,lam,phi):.6f}  vs 4-param {F_global(m,lam,phi):.6f}")
rng=np.random.default_rng(0); bad=0
for _ in range(10000):
    m=rng.integers(2,30); lam=rng.uniform(.05,.99); a=1-lam**2; c=lam**2/a
    ps=(1+(m-1)*lam**2)/(m*(2-lam**2)); phi=rng.uniform(0,1); kp=m*c*(1+phi); km=m*c*(1-phi); p=2*m
    cond=(p-1)*(kp-km)>=km*(2+kp+km)
    if cond!=(phi>=ps) and abs(phi-ps)>1e-9: bad+=1
print("threshold-form equivalence mismatches in 10000 draws:", bad)
for m in (2,5): 
    for lam in (.4,.8):
        F0v=-(m-1)*np.log(1-lam**2)-np.log(1+(m-1)*lam**2); print(f"phi=0 m={m} lam={lam}: formula F={F0v:.6f}, numeric {F_global(m,lam,0):.6f}")
