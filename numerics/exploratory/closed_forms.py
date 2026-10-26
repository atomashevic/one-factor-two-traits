import numpy as np
from scipy.optimize import brentq
def parts(m, lam, phi):
    a=1-lam**2; c=lam**2/a; kp=m*c*(1+phi); km=m*c*(1-phi); p=2*m
    return a,c,kp,km,p
def KL1F(m,lam,phi):
    a,c,kp,km,p=parts(m,lam,phi); return 0.5*((p-1)*np.log1p(km/(p-1))-np.log1p(km))
def TC(m,lam,phi):
    a,c,kp,km,p=parts(m,lam,phi); return 0.5*(-p*np.log(a)-np.log1p(kp)-np.log1p(km))
def CFI(m,lam,phi): return 1-KL1F(m,lam,phi)/TC(m,lam,phi)
def CFI_inf(lam,phi):
    a=1-lam**2; return 1-(lam**2*(1-phi)/a)/(2*(-np.log(a)))
def deff(m,lam,phi):
    a,c,kp,km,p=parts(m,lam,phi); return kp/(1+kp)+km/(1+km)
print("Population CFI of the 1-factor model, two traits, m items each")
for lam in (.5,.7):
    for m in (2,4,8,16):
        row=[CFI(m,lam,phi) for phi in (.3,.5,.7,.8,.9)]
        th=brentq(lambda f: CFI(m,lam,f)-.95, 0.01, .999)
        print(f" lam={lam} m={m:2d} p={2*m:2d}: CFI@phi .3/.5/.7/.8/.9 = "+" ".join(f"{v:.3f}" for v in row)+f" | phi needed for CFI=.95: {th:.3f}")
    thi=brentq(lambda f: CFI_inf(lam,f)-.95, 0.01,.999)
    print(f" lam={lam} m=inf: CFI@phi .3/.5/.7/.8/.9 = "+" ".join(f"{CFI_inf(lam,f):.3f}" for f in (.3,.5,.7,.8,.9))+f" | phi needed: {thi:.3f}")
print("\nNon-invariance: fix kappa_minus=1 and p=8 (m=4); vary lambda (phi solved so that kappa_minus=1)")
for lam in (.4,.5,.6,.7,.8,.9):
    a=1-lam**2; c=lam**2/a; phi=1-1/(4*c)
    if phi<=-1 or phi>=1: continue
    print(f" lam={lam}: phi={phi:.3f}, KL={KL1F(4,lam,phi):.4f} (same for all), CFI={CFI(4,lam,phi):.4f}, RMSEA={np.sqrt(2*KL1F(4,lam,phi)/20):.4f}")
print("\nEffective latent dimension d_eff = sum kappa_j/(1+kappa_j), lam=.7")
for phi in (.3,.5,.8,.95):
    print(f" phi={phi}: "+" ".join(f"m={m}:{deff(m,.7,phi):.2f}" for m in (1,2,3,4,6,8,12,16,32)))
# items per trait needed for the contrast to reach reliability .7 / .8  (kappa/(1+kappa) = r -> kappa = r/(1-r))
print("\nItems per trait needed for contrast reliability >= .80 (kappa_minus >= 4), lam=.7 (c=.96):")
for phi in (.3,.5,.7,.8,.9,.95):
    print(f" phi={phi}: m >= {4/(0.7**2/0.51*(1-phi)):.1f}")
# chart data
import json
chart1=[]
for phi in np.round(np.arange(0.30,0.991,0.01),2):
    for m,lab in ((2,'2 per trait (4 items)'),(4,'4 per trait (8 items)'),(8,'8 per trait (16 items)')):
        chart1.append({'phi':float(phi),'series':lab,'cfi':round(float(CFI(m,.7,phi)),4)})
    chart1.append({'phi':float(phi),'series':'limit, many items','cfi':round(float(CFI_inf(.7,phi)),4)})
json.dump(chart1,open('chart_cfi.json','w'))
chart3=[]
for m in (1,2,3,4,5,6,8,10,12,16,20,24,32):
    for phi in (.5,.8,.95):
        chart3.append({'m':m,'series':f'phi {phi}','deff':round(float(deff(m,.7,phi)),3)})
json.dump(chart3,open('chart_deff.json','w'))
print(len(chart1),len(chart3))
