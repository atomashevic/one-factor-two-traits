import numpy as np
exec(open('fullparam.py').read().split("# 1) reproduce")[0])
def popg(lA,lB,phi):
    lA=np.asarray(lA); lB=np.asarray(lB); mA,mB=len(lA),len(lB); p=mA+mB
    L=np.zeros((p,2)); L[:mA,0]=lA; L[mA:,1]=lB
    S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1); return S
print("Probe 1: block-specific loadings and sizes (A: 3 items at .6, B: 5 items at .8)")
prev=None
for phi in np.round(np.arange(0,0.96,0.05),2):
    R=popg([.6]*3,[.8]*5,phi); b=full_best(R,30); l=np.abs(b.x[:8])
    sd=max(l[:3].std(),l[3:].std())
    print(f"  phi={phi:.2f}: A {l[:3].mean():.3f}  B {l[3:].mean():.3f}  within-block spread {sd:.1e}")
print("Probe 2: mirrored heterogeneous loadings (.5,.6,.7,.8) in both blocks")
lam=np.array([.5,.6,.7,.8])
for phi in np.round(np.arange(0.20,0.56,0.025),3):
    R=popg(lam,lam,phi); b=full_best(R,30); l=np.abs(b.x[:8])
    tilt=np.abs(l[:4]-l[4:]).max()
    print(f"  phi={phi:.3f}: max |A_i - B_i| = {tilt:.4f}")
c=lam**2/(1-lam**2); m=4
for name,L in [("mean lambda^2",np.mean(lam**2)),("mean of lambda, squared",np.mean(lam)**2),("from mean signal-to-noise c",np.mean(c)/(1+np.mean(c)))]:
    print(f"  exchangeable phi* using {name} (L={L:.3f}): {(1+(m-1)*L)/(m*(2-L)):.3f}")
