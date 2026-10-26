import numpy as np
from scipy.optimize import minimize
exec(open("fold_defs.py").read())
# 1. counterexample (1.19)
R=np.array([[1,.1,.5,.5],[.1,1,.5,.5],[.5,.5,1,.1],[.5,.5,.1,1]])
def Fpsi(l,psi,S):
    Sig=np.outer(l,l)+np.diag(psi); A=np.linalg.solve(Sig,S); return np.trace(A)-np.linalg.slogdet(A)[1]-len(l)
l0=np.array([1,.1,.5,.5]); p0=np.array([0,.99,.75,.75])
F0=Fpsi(l0,p0,R); print("F star",F0,"log(275/84)",np.log(275/84))
h=1e-7; print("inward derivative",(Fpsi(l0,p0+np.r_[h,0,0,0],R)-F0)/h,"2/165=",2/165)
# local optimization with psi>=0 bounds from small perturbations of star point
fun=lambda th: Fpsi(th[:4],th[4:],R)
best=[]
rng=np.random.default_rng(3)
for _ in range(30):
    x0=np.r_[l0,p0]+rng.normal(0,.01,8); x0[4:]=np.maximum(x0[4:],0)
    r=minimize(fun,x0,method='L-BFGS-B',bounds=[(None,None)]*4+[(0,2)]*4,options={'ftol':1e-15,'gtol':1e-12,'maxiter':5000})
    best.append((round(r.fun,8),round(r.x[4],6)))
print("local runs near star (F, psi1):",sorted(set(best))[:5])
print("uniform F",np.log(6859/2187), "global from random starts:", min(minimize(fun,np.r_[rng.uniform(0,1,4),rng.uniform(.1,1,4)],method='L-BFGS-B',bounds=[(None,None)]*4+[(0,2)]*4).fun for _ in range(30)))
# 2. existence criterion: m=4, lamB=.8, threshold lamA*
m=4;LB=.64; LAstar=(m-1)*LB/(m+(m-1)*LB); print("lamA* =",np.sqrt(LAstar))
for lA in (.55,.565,.575,.59):
    f=None
    try: f=track(lA,.8,np.arange(0.0002,0.2,0.0002))
    except Exception as e: f='err'
    print(f"lamA={lA}: weaker-block branch survives continuation up to phi={f}")
# 3. quartic roots vs fitted solutions, mild design phi=.05
def quartic(n,k,lA,lB,phi):
    LA,LB=lA**2,lB**2; aA,aB=1-LA,1-LB; bA,bB=1+(n-1)*LA,1+(k-1)*LB; c=lA*lB*phi; q=n+k-1
    A=(n-1)*aA;B=(k-1)*aB;C=k*c;D=n*c;E=n*k*(1-c*c)-A*B;J=q*bA-E;K=q*bB-E;f=q-B;g=q-A
    P=[-bB*C*B, bB*K-D*C*B-C*C*f, bB*D*g+D*K-bA*C*f-C*J, D*D*g-bA*J+C*D*A, bA*D*A]
    return np.roots(P)
rts=quartic(4,4,.68,.72,.05); print("quartic roots",np.sort_complex(rts))
for x0 in (np.r_[[.68]*4,[.05]*4,np.log([1-.68**2]*4),np.log([.99]*4)], None):
    Rm=popg([.68]*4,[.72]*4,.05)
    r=polish(Rm,x0) if x0 is not None else full_best(Rm,20)
    l=np.abs(r.x[:8]); ps=np.exp(r.x[8:]); al,be=l[:4].mean(),l[4:].mean(); x,y=ps[:4].mean(),ps[4:].mean()
    print("solution F %.6f  v=(beta/y)/(alpha/x)=%.6f"%(r.fun,(be/y)/(al/x)))
