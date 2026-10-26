import sympy as sp, mpmath as mp
mp.mp.dps=60
def pieces(mv, Lv, phiv):
    m=mp.mpf(mv); L=mp.mpf(Lv); phi=mp.mpf(phiv); n=m-1; N=2*m-1; a=1-L; b=1+n*L; h=m*L*phi
    z=L*(n+m*phi)/N; t=mp.sqrt(z); r=1-z; U=b+h; V=b-h
    def F_mine(e):
        al,be=t+e,t-e; x,y=1-al**2,1-be**2
        g=lambda q:q-mp.log(q)-1
        S=mp.matrix([[x+m*al**2,m*al*be],[m*al*be,y+m*be**2]]); R=mp.matrix([[b,h],[h,b]])
        return n*g(a/x)+n*g(a/y)+(mp.inverse(S)*R)[0,0]+(mp.inverse(S)*R)[1,1]+mp.log(mp.det(S))-mp.log(mp.det(R))-2
    def path_theirs(e):   # copied structure from verify_symmetry_breaking.py critical_checks()
        x,y=1-(t+e)**2,1-(t-e)**2
        pp,mm,off=U-e*e,r+N*e*e,2*n*t*e
        det=pp*mm-off*off
        return n*(a*(1/x+1/y)+mp.log(x*y))+(U*mm+V*pp)/det+mp.log(det)
    e=mp.mpf('1e-15')
    Em=(F_mine(e)-F_mine(0))/e**2; Et=(path_theirs(e)-path_theirs(0))/e**2
    Ec=4*n*m*(1+z)*((2-L)*z-L)/((1-z)**3*(1+N*z))
    return Em,Et,Ec
for (mv,Lv,phiv) in [(2,'0.49','0'),(2,'0.49','0.3'),(4,'0.64','0.2'),(5,'0.3','0')]:
    Em,Et,Ec=pieces(mv,Lv,phiv)
    print(f"m={mv} L={Lv} phi={phiv}: E(mine)={mp.nstr(Em,12)}  E(their path)={mp.nstr(Et,12)}  E(claimed)={mp.nstr(Ec,12)}")
