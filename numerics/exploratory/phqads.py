import numpy as np
exec(open("fold_defs.py").read())
df=16*13//2
print("Merger-only population ML fit, PHQ-ADS layout (9 depression items, 7 anxiety items), equal loadings lam in both blocks")
print("phi    lam   F_min   RMSEA   CFI    favoured-share   | reported 1F RMSEA/CFI, 2F RMSEA/CFI")
rep={.912:(".091/.956",".085/.962"),.865:(".083/.954",".071/.967"),.653:(".130/.862",".085/.941")}
for phi in (.912,.865,.653):
    for lam in (.7,.75,.8):
        R=popg([lam]*9,[lam]*7,phi); b=full_best(R,12); FB=-np.linalg.slogdet(R)[1]
        l=np.abs(b.x[:16]); M=np.array([[1,phi],[phi,1]]); w=np.linalg.solve(M,[l[:9].mean()/lam,l[9:].mean()/lam]); w=np.abs(w)/np.abs(w).sum()
        print(f"{phi:.3f}  {lam:.2f}  {b.fun:.3f}   {np.sqrt(b.fun/df):.3f}   {1-b.fun/FB:.3f}   {w.max():.3f} ({'dep' if w[0]>w[1] else 'anx'})   | {rep[phi][0]}, {rep[phi][1]}")
