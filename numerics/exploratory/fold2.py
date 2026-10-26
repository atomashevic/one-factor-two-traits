import numpy as np
exec(open("fold_defs.py").read())
def track2(lA,lB,phis):
    try: return track(lA,lB,phis)
    except Exception as e: return 'err:'+str(e)
for lab,lA,lB in [(".65/.75",.65,.75),("strong .6/.8",.6,.8)]:
    print(lab, track2(lA,lB,np.arange(0.002,0.42,0.001)))
for phi in (.23,.25):
    res=locmins(popg([.68]*4,[.72]*4,phi),40); print("mild",phi,len(res),[k for k,_,_ in res])
for phi in (.005,.02):
    res=locmins(popg([.6]*4,[.8]*4,phi),40); print("strong",phi,len(res),[k for k,_,_ in res])
    res=locmins(popg([.65]*4,[.75]*4,phi),40); print(".65/.75",phi,len(res),[k for k,_,_ in res])
