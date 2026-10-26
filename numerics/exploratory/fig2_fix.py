import numpy as np, json
from scipy.optimize import brentq
exec(open('beyond2.py').read().split("for lam in (")[0])
psm=brentq(lambda f: anti_min_eig(np.array([.6,.65,.75,.8]),f),.4,.5,xtol=1e-8)
pse=(1+3*.49)/(4*1.51)
print(f"eq {pse:.6f} mir {psm:.6f}")
rows=json.load(open('fig2_rows.json'))
rows=[r for r in rows if not (r['slug']=='eq' and r['phi']==0.4089) and not (r['slug']=='mir' and r['phi']==0.4411)]
rows+= [dict(phi=round(pse,4)+0.0,slug='eq',share=0.5), dict(phi=round(psm,4),slug='mir',share=0.5)]
rows.sort(key=lambda r:(['eq','mir','mild','strong'].index(r['slug']),r['phi']))
json.dump(rows,open('fig2_rows.json','w'))
print(len(rows), round(pse,4), round(psm,4))
