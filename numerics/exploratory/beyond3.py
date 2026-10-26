import numpy as np, json, time
from scipy.optimize import minimize
exec(open('fullparam.py').read().split("# 1) reproduce")[0])
rng=np.random.default_rng(7)
def popg(lA,lB,phi):
    lA=np.asarray(lA,float); lB=np.asarray(lB,float); mA,mB=len(lA),len(lB); p=mA+mB
    L=np.zeros((p,2)); L[:mA,0]=lA; L[mA:,1]=lB
    S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1); return S
def share(lAbar,lBbar,lamA,lamB,phi):
    M=np.array([[1,phi],[phi,1]]); w=np.abs(np.linalg.solve(M,[lAbar/lamA,lBbar/lamB])); return max(w)/w.sum()
def polish(R,x0):
    r=minimize(F,x0,args=(R,),jac=grad,method='L-BFGS-B',options={'maxiter':20000,'ftol':1e-16,'gtol':1e-13}); return r
t0=time.time()
# epsilon scaling at the exchangeable threshold
m=4; L=.49; ps=(1+(m-1)*L)/(m*(2-L))
eps=[.0005,.001,.002,.004,.008,.016]; tilts=[]
for e in eps:
    R=popg([.7-e]*4,[.7+e]*4,ps); b=full_best(R,30); b=polish(R,b.x); l=np.abs(b.x[:8])
    t=share(l[:4].mean(),l[4:].mean(),.7-e,.7+e,ps)-.5; tilts.append(t); print(f"eps={e}: tilt {t:.5f}",flush=True)
for i in range(len(eps)-2):
    sl=np.polyfit(np.log(eps[i:i+3]),np.log(tilts[i:i+3]),1)[0]; print(f"local slope over eps {eps[i]}..{eps[i+2]}: {sl:.3f}")
print("overall slope",np.polyfit(np.log(eps),np.log(tilts),1)[0], "time",time.time()-t0)
json.dump(dict(eps=eps,tilts=tilts),open('eps_scaling.json','w'))
