import numpy as np, itertools, sys, time
def set_partitions(n):
    # restricted growth strings
    a=[0]*n
    def rec(i,mx):
        if i==n:
            yield tuple(a); return
        for v in range(mx+2):
            a[i]=v
            yield from rec(i+1, max(mx,v))
    yield from rec(1,0) if n>0 else iter([()])
def pop2(m1,m2,lam,phi):
    p=m1+m2; L=np.zeros((p,2)); L[:m1,0]=lam; L[m1:,1]=lam
    Phi=np.array([[1,phi],[phi,1]]); S=L@Phi@L.T; np.fill_diagonal(S,1); return S
def H_ew(B):  # EGAnet: -trace(rho %*% log(rho)) elementwise, na.rm (rows with zeros dropped)
    rho=np.abs(B)/B.shape[0]
    with np.errstate(divide='ignore',invalid='ignore'):
        d=np.einsum('ij,ji->i',rho,np.log(rho))
    return -np.nansum(d)
def H_vn(B):
    ev=np.linalg.eigvalsh(np.abs(B)/B.shape[0]); ev=ev[ev>1e-15]; return -np.sum(ev*np.log(ev))
def tefi(R,s,H,cache):
    K=max(s)+1; S=cache['S_'+H.__name__]
    Sk=[]
    for k in range(K):
        idx=tuple(i for i,v in enumerate(s) if v==k)
        key=(H.__name__,idx)
        if key not in cache: cache[key]=H(R[np.ix_(idx,idx)])
        Sk.append(cache[key])
    Sk=np.array(Sk); return (Sk.mean()-S)+(S-Sk.sum())*np.sqrt(K)
def describe(s,m1):
    K=max(s)+1; blocks=[[i for i,v in enumerate(s) if v==k] for k in range(K)]
    out=[]
    for b in blocks:
        a=sum(1 for i in b if i<m1); c=len(b)-a; out.append(f"{a}A+{c}B" if a and c else (f"{a}A" if a else f"{c}B"))
    return "{"+", ".join(sorted(out))+"}"
m1=m2=int(sys.argv[1]) if len(sys.argv)>1 else 4
lam=float(sys.argv[2]) if len(sys.argv)>2 else .7
P=list(set_partitions(m1+m2)); print(f"p={m1+m2}, lambda={lam}, partitions={len(P)}")
for phi in [0.05,0.2,0.3,0.5,0.7,0.9,0.99]:
    R=pop2(m1,m2,lam,phi); cache={'S_H_ew':H_ew(R),'S_H_vn':H_vn(R)}
    res={}
    for H in (H_ew,H_vn):
        vals=np.array([tefi(R,s,H,cache) for s in P]); order=np.argsort(vals)
        true=tuple([0]*m1+[1]*m2); tv=tefi(R,true,H,cache)
        rank=int((vals<tv-1e-12).sum())+1
        best=P[order[0]]
        res[H.__name__]=(describe(best,m1),vals[order[0]],tv,rank)
    e=res['H_ew']; v=res['H_vn']
    print(f"phi={phi:4.2f} | EGAnet-TEFI argmin {e[0]:<22} ({e[1]:7.3f}); true={e[2]:7.3f}, rank of true {e[3]:>6} | VN-TEFI argmin {v[0]:<26} ({v[1]:7.3f}); true rank {v[3]}")
