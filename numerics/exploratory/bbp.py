import numpy as np
from scipy.stats import ncx2, chi2
rng=np.random.default_rng(7)
m=10; p=2*m; lam=.5; a=1-lam**2; c=lam**2/a; N=300; g=p/N
print(f"p={p}, N={N}, gamma={g:.4f}, sqrt(gamma)={np.sqrt(g):.3f}; threshold kappa_- > sqrt(gamma) <=> phi < {1-np.sqrt(g)/(m*c):.3f}")
print(" phi  kappa- | mean l2/a  predicted  bulk-edge | overlap^2 sim  pred | LRT power(1F)")
for phi in (.70,.80,.85,.88,.90,.92,.94,.96):
    L=np.zeros((p,2)); L[:m,0]=lam; L[m:,1]=lam
    Phi=np.array([[1,phi],[phi,1]]); S=L@Phi@L.T; np.fill_diagonal(S,1)
    km=m*c*(1-phi); u=np.r_[np.ones(m),-np.ones(m)]/np.sqrt(p)
    C=np.linalg.cholesky(S); l2=[]; ov=[]
    for r in range(300):
        X=rng.standard_normal((N,p))@C.T; Sh=X.T@X/N
        w,V=np.linalg.eigh(Sh); l2.append(w[-2]); ov.append((V[:,-2]@u)**2)
    pred=(1+km)*(1+g/km) if km>np.sqrt(g) else (1+np.sqrt(g))**2
    ovp=(1-g/km**2)/(1+g/km) if km>np.sqrt(g) else 0.0
    F=(p-1)*np.log1p(km/(p-1))-np.log1p(km); df=p*(p+1)//2-2*p
    pw=ncx2.sf(chi2.ppf(.95,df),df,N*F)
    print(f" {phi:.2f} {km:6.3f} | {np.mean(l2)/a:9.3f} {pred:9.3f} {(1+np.sqrt(g))**2:9.3f} | {np.mean(ov):12.3f} {ovp:5.3f} | {pw:.3f}")
# TEFI element-wise identity and VN small-correlation expansion
def H_ew(B):
    rho=np.abs(B)/B.shape[0]; return -np.sum(rho*np.log(rho))
def H_vn(B):
    ev=np.linalg.eigvalsh(B/B.shape[0]); ev=ev[ev>1e-15]; return -np.sum(ev*np.log(ev))
A=rng.uniform(.05,.6,(9,9)); A=(A+A.T)/2; np.fill_diagonal(A,1); pp=9
off=~np.eye(pp,dtype=bool)
ident=np.log(pp)+np.sum(A[off]*np.log(pp/A[off]))/pp
print(f"\nElement-wise identity: H_ew={H_ew(A):.10f}  log p + (1/p) sum a log(p/a) = {ident:.10f}")
for s in (.05,.1,.2,.3):
    B=np.eye(pp)+s*(rng.uniform(-1,1,(pp,pp))); B=(B+B.T)/2; np.fill_diagonal(B,1)
    if np.linalg.eigvalsh(B).min()<=0: continue
    approx=np.log(pp)-np.sum((B-np.eye(pp))**2)/(2*pp)
    print(f" VN expansion, |r|<~{s}: exact {H_vn(B):.6f}  log p - ||R-I||^2/(2p) {approx:.6f}")
