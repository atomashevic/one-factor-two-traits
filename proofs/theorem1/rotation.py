import numpy as np
def pop(m,lam,phi):
    p=2*m; L=np.zeros((p,2)); L[:m,0]=lam; L[m:,1]=lam
    S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1); return S
def F(l,psi,R):
    S=np.outer(l,l)+np.diag(psi); A=np.linalg.solve(S,R); return np.trace(A)-np.linalg.slogdet(A)[1]-len(R)
for m,lam in [(2,.7),(3,.5),(5,.9)]:
    R=pop(m,lam,0.0); L2=lam**2; n=m-1; N=2*m-1; z=L2*n/N; psi=np.full(2*m,1-z)
    eA=np.r_[np.ones(m),np.zeros(m)]/np.sqrt(m); eB=np.r_[np.zeros(m),np.ones(m)]/np.sqrt(m); c=np.sqrt(z*2*m)
    Fs=F(np.sqrt(z)*np.ones(2*m),psi,R)
    vals=[F(c*(np.cos(t)*eA+np.sin(t)*eB),psi,R)-Fs for t in np.linspace(0,np.pi/2,7)]
    t=np.pi/4-1e-3; l=c*(np.cos(t)*eA+np.sin(t)*eB); Sg=np.outer(l,l)+np.diag(psi); K=np.linalg.inv(Sg); H=K-K@R@K
    gl=2*H@l; gpsi=np.diag(H); pred=(np.diag(Sg)-1)/psi**2
    print(f"m={m} lam={lam}: F on rotation circle minus F_sym, max |.| = {max(abs(v) for v in vals):.1e}; "
          f"at t=pi/4-1e-3: |dF/dl| = {np.abs(gl).max():.1e}, |dF/dpsi| = {np.abs(gpsi).max():.2e}, "
          f"matches (Sigma_ii-1)/psi^2 to {np.abs(gpsi-pred).max():.1e}")
