import numpy as np, json
exec(open("fold_defs.py").read())
def share(lAbar,lBbar,lamA,lamB,phi):
    M=np.array([[1,phi],[phi,1]]); w=np.abs(np.linalg.solve(M,[lAbar/lamA,lBbar/lamB])); return w[0]/w.sum()
lA,lB=.68,.72
x=np.r_[[lA]*4,[0.01]*4,np.log([1-lA**2]*4),np.log([.99]*4)]
out=[]; phis=np.round(np.arange(0.0005,0.2416,0.0005),4)
for phi in phis:
    R=popg([lA]*4,[lB]*4,phi)
    r=minimize(F,x,args=(R,),jac=grad,method='L-BFGS-B',options={'maxiter':20000,'ftol':1e-16,'gtol':1e-12})
    l=np.abs(r.x[:8]); x=r.x
    if l[:4].mean()<l[4:].mean(): print("jumped at",phi); break
    if abs(phi*50-round(phi*50))<1e-9 or phi==0.2415:
        out.append(dict(phi=float(phi),share=round(float(share(l[:4].mean(),l[4:].mean(),lA,lB,phi)),4),F=float(r.fun)))
out=[dict(phi=0.0,share=1.0)]+[dict(phi=o['phi'],share=o['share']) for o in out]
print(out)
json.dump(out,open('fig2_local.json','w'))
