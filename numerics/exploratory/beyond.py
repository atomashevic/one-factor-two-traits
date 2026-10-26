import numpy as np, json
from scipy.optimize import minimize
exec(open('fullparam.py').read().split("# 1) reproduce")[0])
rng=np.random.default_rng(2026)
def popg(lA,lB,phi):
    lA=np.asarray(lA,float); lB=np.asarray(lB,float); mA,mB=len(lA),len(lB); p=mA+mB
    L=np.zeros((p,2)); L[:mA,0]=lA; L[mA:,1]=lB
    S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1); return S
# 1) extended block-uniformity lemma: random block-specific designs, full parameter space
worst=0; n=0
for _ in range(60):
    mA,mB=rng.integers(2,6,2); la,lb=rng.uniform(.3,.9,2); phi=rng.uniform(.01,.95)
    R=popg([la]*mA,[lb]*mB,phi); b=full_best(R,25); l=np.abs(b.x[:mA+mB]); ps=np.exp(b.x[mA+mB:])
    worst=max(worst,l[:mA].std(),l[mA:].std(),ps[:mA].std(),ps[mA:].std()); n+=1
print(f"extended lemma check: {n} random block-specific designs (sizes 2 to 5, loadings .3 to .9), max within-block spread {worst:.1e}")
# 2) tilt measure: composite weights w solving l_A/lam_A = wA + phi wB, l_B/lam_B = phi wA + wB; share of favoured trait
def share(lAbar,lBbar,lamA,lamB,phi):
    M=np.array([[1,phi],[phi,1]]); w=np.linalg.solve(M,[lAbar/lamA,lBbar/lamB]); w=np.abs(w)
    return max(w)/(w.sum())
grid=np.round(np.arange(0,0.901,0.025),3)
designs=[("Equal blocks, loadings .70","eq",[.7]*4,[.7]*4),
         ("Mirrored loadings .60 to .80","mir",[.6,.65,.75,.8],[.6,.65,.75,.8]),
         ("Unequal blocks, .68 vs .72","mild",[.68]*4,[.72]*4),
         ("Unequal blocks, .60 vs .80","strong",[.6]*4,[.8]*4)]
rows=[]
for lab,slug,lA,lB in designs:
    for phi in grid:
        R=popg(lA,lB,phi); b=full_best(R,20); l=np.abs(b.x[:8])
        sh=share(l[:4].mean(),l[4:].mean(),np.mean(lA),np.mean(lB),phi)
        rows.append(dict(phi=float(phi),design=lab,slug=slug,share=round(float(sh),4)))
    print(lab, [r['share'] for r in rows if r['slug']==slug][::4])
json.dump(rows,open('beyond_rows.json','w'))
# 3) epsilon^(1/3) check at the exchangeable threshold
m=4; L=.49; ps=(1+(m-1)*L)/(m*(2-L))
print(f"tilt at phi* = {ps:.4f} for loadings .7-eps vs .7+eps:")
eps=[.0025,.005,.01,.02,.04]; tilts=[]
for e in eps:
    R=popg([.7-e]*4,[.7+e]*4,ps); b=full_best(R,25); l=np.abs(b.x[:8])
    t=share(l[:4].mean(),l[4:].mean(),.7-e,.7+e,ps)-.5; tilts.append(t); print(f"  eps={e}: tilt {t:.4f}")
sl=np.polyfit(np.log(eps),np.log(tilts),1)[0]; print(f"  log-log slope {sl:.3f} (pitchfork unfolding predicts 1/3)")
