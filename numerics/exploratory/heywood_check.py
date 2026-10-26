import numpy as np
rng=np.random.default_rng(4)
def randcorr(p):
    A=rng.normal(size=(p,p+2)); S=A@A.T; d=np.sqrt(np.diag(S)); return S/np.outer(d,d)
worst=dict(K=0,step=0,final=0,fd=0,var=0)
for trial in range(400):
    p=int(rng.integers(3,9)); R=randcorr(p); j=int(rng.integers(p)); o=[i for i in range(p) if i!=j]
    l=R[:,j].copy(); l[j]=1.0; psi=1-l**2; psi[j]=0.0
    Sig=np.outer(l,l)+np.diag(psi); K=np.linalg.inv(Sig)
    w=np.zeros(p); w[o]=R[o,j]/psi[o]; s=np.sum(w[o]*R[o,j])
    # (4.0): closed-form precision from the regression factorization
    Kc=np.zeros((p,p)); Kc[j,j]=1+s
    for i in o: Kc[i,j]=Kc[j,i]=-w[i]; Kc[i,i]=1/psi[i]
    worst['K']=max(worst['K'],np.abs(K-Kc).max())
    # intermediate identity: H_jj = s + s^2 - sum w^2 - sum_{i!=k} w_i w_k R_ik
    Hjj=(K-K@R@K)[j,j]
    off=sum(w[i]*w[k]*R[i,k] for i in o for k in o if i!=k)
    worst['step']=max(worst['step'],abs(Hjj-(s+s*s-np.sum(w[o]**2)-off)))
    # final formula (4.3), general R
    f43=-sum(w[i]*w[k]*(R[i,k]-R[i,j]*R[k,j]) for i in o for k in o if i!=k)
    worst['final']=max(worst['final'],abs(Hjj-f43))
    # variance form: sum w_i^2 psi_i - Var(sum w_i r_i), r_i = X_i - R_ij X_j
    C=np.array([[R[i,k]-R[i,j]*R[k,j] for k in o] for i in o]); wo=w[o]
    worst['var']=max(worst['var'],abs(Hjj-(np.sum(wo**2*psi[o])-wo@C@wo)))
    # one-sided finite difference in psi_j
    F=lambda S_: np.trace(np.linalg.solve(S_,R))+np.linalg.slogdet(S_)[1]-np.linalg.slogdet(R)[1]-p
    e=1e-6; E=np.zeros((p,p)); E[j,j]=e
    worst['fd']=max(worst['fd'],abs((F(Sig+E)-F(Sig))/e-f43)/max(1,abs(f43)))
print("400 random correlation matrices, p = 3..8:")
for k,v in worst.items(): print(f"  max error [{k}]: {v:.1e}")
