import numpy as np
exec(open("fold_defs.py").read())
TC=lambda m,L: -(m-1)*np.log(1-L)-np.log(1+(m-1)*L)
kap=lambda m,L: m*L/(1-L)
# find disagreement cases
cands=[]
for mA in range(2,9):
  for mB in range(mA+1,10):
    for lA in np.arange(.4,.95,.05):
      for lB in np.arange(.3,.95,.05):
        LA,LB=lA**2,lB**2
        if np.sign(TC(mA,LA)-TC(mB,LB))!=np.sign(kap(mA,LA)-kap(mB,LB)) and abs(TC(mA,LA)-TC(mB,LB))>.15 and abs(kap(mA,LA)-kap(mB,LB))>.4:
          cands.append((mA,mB,round(lA,2),round(lB,2),round(TC(mA,LA)-TC(mB,LB),3),round(kap(mA,LA)-kap(mB,LB),3)))
print(len(cands)); print(cands[:8])
def favoured(mA,mB,lA,lB,phi):
    R=popg([lA]*mA,[lB]*mB,phi); b=full_best(R,20); p=mA+mB; l=np.abs(b.x[:p])
    M=np.array([[1,phi],[phi,1]]); w=np.linalg.solve(M,[l[:mA].mean()/lA,l[mA:].mean()/lB]); return 'A' if w[0]>w[1] else 'B', np.round(w/np.abs(w).sum(),3)
for cse in cands[:4]:
    mA,mB,lA,lB=cse[:4]
    print(cse, "TC favours", 'A' if cse[4]>0 else 'B', "kappa favours", 'A' if cse[5]>0 else 'B', [ (phi,)+favoured(mA,mB,lA,lB,phi) for phi in (.01,.2,.5,.8)])
