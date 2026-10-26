import numpy as np, itertools
exec(open('partitions.py').read().split("m1=m2=int")[0])   # reuse functions
P_all=list(set_partitions(8))
bal=[s for s in P_all if max(s)==1 and sum(1 for v in s if v==0)==4]   # 35 balanced 4+4 splits
print("2021 Simulation-1 design (2 traits x 4 items) in the population")
print(" lam   phi | implemented TEFI: rank of true among 35 balanced / among 4140 all | genuine VN: rank among 35 / 4140")
for lam in (.40,.55,.70):
    for phi in (0.001,.5,.7):
        R=pop2(4,4,lam,phi); cache={'S_H_ew':H_ew(R),'S_H_vn':H_vn(R)}; true=tuple([0]*4+[1]*4)
        out=[]
        for H in (H_ew,H_vn):
            tv=tefi(R,true,H,cache)
            rb=1+sum(tefi(R,s,H,cache)<tv-1e-12 for s in bal)
            ra=1+sum(tefi(R,s,H,cache)<tv-1e-12 for s in P_all)
            out+= [rb,ra]
        print(f" {lam:.2f} {phi:5.3f} | {out[0]:>3} / {out[1]:>4} | {out[2]:>3} / {out[3]:>4}")
# O-information sign for Gaussian one-factor vs two-factor models
def hG(S): return 0.5*np.linalg.slogdet(2*np.pi*np.e*S)[1]
def Oinfo(S):
    n=S.shape[0]; HX=hG(S)
    return (n-2)*HX+sum(hG(S[[j]][:,[j]])-hG(np.delete(np.delete(S,j,0),j,1)) for j in range(n))
rng=np.random.default_rng(3); neg1=0; neg2=0; mins=[]
for r in range(2000):
    p=rng.integers(3,12); l=rng.uniform(.1,.95,p); S=np.outer(l,l); np.fill_diagonal(S,1)
    o=Oinfo(S); mins.append(o); neg1+= o<0
for r in range(2000):
    m1,m2=rng.integers(2,7,2); phi=rng.uniform(-.9,.9); l=rng.uniform(.1,.95,m1+m2)
    L=np.zeros((m1+m2,2)); L[:m1,0]=l[:m1]; L[m1:,1]=l[m1:]; S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1)
    neg2+= Oinfo(S)<0
print(f"\nO-information: one-factor Gaussian models with Omega<0: {neg1}/2000 (min {min(mins):.2e}); two-factor: {neg2}/2000")
