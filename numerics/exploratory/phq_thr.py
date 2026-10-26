import numpy as np
from scipy.optimize import brentq
exec(open("fold_defs.py").read())
df=104
def fit(lam,phi):
    R=popg([lam]*9,[lam]*7,phi); b=full_best(R,10); FB=-np.linalg.slogdet(R)[1]
    return b.fun, np.sqrt(b.fun/df), 1-b.fun/FB, FB
for lam in (.7,.75,.8):
    r=brentq(lambda f: fit(lam,f)[1]-.06,.5,.99,xtol=1e-4)
    c=brentq(lambda f: fit(lam,f)[2]-.95,.5,.99,xtol=1e-4)
    print(f"lam {lam}: RMSEA<=.06 once phi>={r:.3f}; CFI>=.95 once phi>={c:.3f}; F_B at phi=.865: {fit(lam,.865)[3]:.2f}")
