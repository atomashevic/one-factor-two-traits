import numpy as np
exec(open('fullparam.py').read().split("# 1) reproduce")[0])
exec(open('global1f.py').read().split("lam=.7\nfor m in")[0])
for m in (2,3,4):
    for lam in (.4,.6,.8,.9):
        for phi in np.linspace(0,.95,12):
            S=pop(m,lam,phi); b=full_best(S,25); p=2*m; l=np.abs(b.x[:p]); psi=np.exp(b.x[p:])
            sd=max(l[:m].std(),l[m:].std(),psi[:m].std(),psi[m:].std())
            if sd>1e-4:
                Fb=F_global(m,lam,phi)
                print(f"m={m} lam={lam} phi={phi:.3f}: spread {sd:.2e}, F_full={b.fun:.8f}, F_blockuniform={Fb:.8f}, diff={Fb-b.fun:.1e}; loadings {np.round(l,3)} psi {np.round(psi,3)}")
