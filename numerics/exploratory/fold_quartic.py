import numpy as np
def quartic(n,k,lA,lB,phi):
    LA,LB=lA**2,lB**2; aA,aB=1-LA,1-LB; bA,bB=1+(n-1)*LA,1+(k-1)*LB; c=lA*lB*phi; q=n+k-1
    A=(n-1)*aA;B=(k-1)*aB;C=k*c;D=n*c;E=n*k*(1-c*c)-A*B;J=q*bA-E;K=q*bB-E;f=q-B;g=q-A
    return [-bB*C*B, bB*K-D*C*B-C*C*f, bB*D*g+D*K-bA*C*f-C*J, D*D*g-bA*J+C*D*A, bA*D*A]
def npos(lA,lB,phi,m=4):
    r=np.roots(quartic(m,m,lA,lB,phi)); return sum(1 for z in r if abs(z.imag)<1e-9*max(1,abs(z)) and z.real>0)
for lA in (.565,.575,.59,.6,.68):
    lB=.8 if lA<.68 else .72
    phis=np.logspace(-6,np.log10(.45),4000); cnt=[npos(lA,lB,p) for p in phis]
    last3=[p for p,c in zip(phis,cnt) if c==3]
    print(lA,lB,"three positive roots up to phi =", (max(last3) if last3 else None), "| criterion holds:", lA**2>3*lB**2/(4+3*lB**2))
