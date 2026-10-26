import numpy as np
exec(open("fold_defs.py").read())
g=lambda t: t-np.log(t)-1
def Fred(x,y,mA,mB,lA,lB,phi):
    LA,LB=lA**2,lB**2; aA,aB=1-LA,1-LB; bA,bB=1+(mA-1)*LA,1+(mB-1)*LB; h=lA*lB*phi*np.sqrt(mA*mB)
    M=np.array([[bA/x,h/np.sqrt(x*y)],[h/np.sqrt(x*y),bB/y]]); th=np.linalg.eigvalsh(M)
    return (mA-1)*g(aA/x)+(mB-1)*g(aB/y)+g(th[0]), th
def popg2(mA,mB,lA,lB,phi): return popg([lA]*mA,[lB]*mB,phi)
def stat_eqs(al,be,mA,mB,lA,lB,phi):
    LA,LB=lA**2,lB**2; bA,bB=1+(mA-1)*LA,1+(mB-1)*LB; c=lA*lB*phi
    x,y=1-al**2,1-be**2; gam=mA*al**2/x+mB*be**2/y
    return bA*al/x+c*mB*be/y-(1+gam)*al, bB*be/y+c*mA*al/x-(1+gam)*be
# check on solutions
for (mA,mB,lA,lB,phi) in [(4,4,.68,.72,.05),(4,4,.68,.72,.5),(3,5,.6,.8,.3),(2,6,.85,.55,.4)]:
    R=popg2(mA,mB,lA,lB,phi); b=full_best(R,20); p=mA+mB; l=np.abs(b.x[:p]); ps=np.exp(b.x[p:])
    al,be=l[:mA].mean(),l[mA:].mean(); x,y=ps[:mA].mean(),ps[mA:].mean()
    Fr,th=Fred(x,y,mA,mB,lA,lB,phi)
    print((mA,mB,lA,lB,phi),"full F",round(b.fun,10),"reduced F",round(Fr,10),"psi=1-l^2 err",abs(x-(1-al**2)),"eqs",np.round(stat_eqs(al,be,mA,mB,lA,lB,phi),9))
