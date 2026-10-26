import numpy as np
exec(open("fold_defs.py").read())
ps=(1+3*.49)/(4*1.51)
E=[.001,.0025,.005,.01,.02,.04]; D=[]
for e in E:
    f=track(.7-e,.7+e,np.arange(0.0005,0.42,0.0005)); D.append(ps-f); print(e, round(f,4), round(ps-f,4))
E=np.array(E); D=np.array(D)
for i in range(len(E)-2): print("local slope", E[i], E[i+2], round(np.polyfit(np.log(E[i:i+3]),np.log(D[i:i+3]),1)[0],3))
