import numpy as np, json
from scipy.optimize import minimize
exec(open('fullparam.py').read().split("# 1) reproduce")[0])
rng=np.random.default_rng(11)
def popg(lA,lB,phi):
    lA=np.asarray(lA,float); lB=np.asarray(lB,float); mA,mB=len(lA),len(lB); p=mA+mB
    L=np.zeros((p,2)); L[:mA,0]=lA; L[mA:,1]=lB
    S=L@np.array([[1,phi],[phi,1]])@L.T; np.fill_diagonal(S,1); return S
def share(lAbar,lBbar,lamA,lamB,phi):
    M=np.array([[1,phi],[phi,1]]); w=np.abs(np.linalg.solve(M,[lAbar/lamA,lBbar/lamB])); return max(w)/w.sum()
designs=[("Equal loadings .70","eq",[.7]*4,[.7]*4,.4089),
         ("Mirrored loadings .60 to .80","mir",[.6,.65,.75,.8],[.6,.65,.75,.8],.4411),
         ("Unequal blocks, .68 vs .72","mild",[.68]*4,[.72]*4,None),
         ("Unequal blocks, .60 vs .80","strong",[.6]*4,[.8]*4,None)]
grid=np.round(np.arange(0,0.9001,0.02),3)
rows=[]; worst=0
for lab,slug,lA,lB,thr in designs:
    g=sorted(set(list(grid)+([round(thr,4)] if thr else [])+([.39,.40,.405] if slug=='eq' else [])+([.42,.43,.435] if slug=='mir' else [])))
    for phi in g:
        R=popg(lA,lB,phi); b=full_best(R,16)
        b=minimize(F,b.x,args=(R,),jac=grad,method='L-BFGS-B',options={'maxiter':20000,'ftol':1e-16,'gtol':1e-13})
        l=np.abs(b.x[:8]); sh=share(l[:4].mean(),l[4:].mean(),np.mean(lA),np.mean(lB),phi)
        if slug=='eq':  # closed-form check
            L=.49; m=4; n=3; bb=1+n*L; q=m*L*phi/bb; qs=L/(2-L)
            if q<qs:
                s2=L+(2-L)*q*q; al=np.sqrt((s2+2*q)/1)/1; 
                apb=np.sqrt(s2+2*q); amb=np.sqrt(max(s2-2*q,0)); a_=(apb+amb)/2; b_=(apb-amb)/2
            else:
                a_=b_=np.sqrt(L*(n+m*phi)/(2*m-1))
            shc=share(a_,b_,.7,.7,phi); worst=max(worst,abs(shc-sh))
        rows.append(dict(phi=float(phi),slug=slug,share=round(float(sh),4)))
print("max |share - closed form| for equal design:",worst)
json.dump(rows,open('fig2_rows.json','w')); print(len(rows))
for s in ['eq','mir','mild','strong']:
    print(s,[(r['phi'],r['share']) for r in rows if r['slug']==s and r['phi'] in (0.0,0.2,0.4,0.4089,0.44,0.4411,0.5,0.7,0.9)])
