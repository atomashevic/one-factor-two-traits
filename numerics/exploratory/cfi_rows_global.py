import numpy as np
exec(open('global1f.py').read().split("lam=.7\nfor m in")[0])
print("cross-check vs independent full-parameter optimizer: m=2,lam=.7,phi=.3 ->", round(F_global(2,.7,.3),4), "(checker .2349); m=4,lam=.9,phi=.5 ->", round(F_global(4,.9,.5),3), "(checker 2.992)")
lam=.7; grid=list(np.round(np.arange(0.30,0.981,0.02),2))+[0.99]
cinf=lambda phi:1-np.log1p(lam**2*(1-phi)/(2*(1-lam**2)))/(-np.log(1-lam**2))
rows=[]
for slug,lab,m in (('items4','4 items (2 per trait)',2),('items8','8 items',4),('items16','16 items',8),('limit','Many items (limit)',None)):
    for phi in grid:
        if m is None:
            ps=lam**2/(2-lam**2)
            v=cinf(phi) if phi>=ps else 1-F_global(10**6,lam,phi)/F0(10**6,lam,phi)
        else:
            ps=phistar(m,lam)
            f=F_sym(m,lam,phi) if phi>=ps else F_global(m,lam,phi)
            v=1-f/F0(m,lam,phi)
        rows.append("{phi:%.2f, series:'%s', slug:'%s', cfi:%.4f}"%(phi,lab,slug,v))
s="["+",".join(rows)+"]"; open('cfi_rows_global.txt','w').write(s); print(len(rows)); print(s)
from scipy.optimize import brentq
print("limit threshold exact:", round(brentq(lambda f: cinf(f)-.95,.4,.999),4))
for m in (2,4,8,16): print("phi* m=%d:"%m, round(phistar(m,lam),3))
