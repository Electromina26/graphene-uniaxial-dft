import numpy as np, json, os, sys
from analyze import analyze, tb_gap, d0
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'font.family':'serif','axes.linewidth':0.8,'figure.dpi':200})
TH=[0,10,20,30]; DATA={th:analyze(th) for th in TH}; DATA={k:v for k,v in DATA.items() if v is not None}
cols={0:'#1f4e79',10:'#2e8b57',20:'#d2691e',30:'#8b1a1a'}
out={}
# ---- ideal strain/stress: quadratic through the 3 points around max (use spline over all)
from scipy.interpolate import CubicSpline
for th,D in DATA.items():
    e,P=D['eps'],D['Pn']; cs=CubicSpline(e,P); x=np.linspace(0,e.max(),3001); y=cs(x); i=np.argmax(y)
    out.setdefault(th,{}).update(eps_m=float(x[i]),P_m=float(y[i]),peak_inside=bool(i<len(x)-5))
# ---- global symmetric fit: P = Y e + (D0+D6 c6) e^2 + (E0+E6 c6 + E12 c12) e^3 + (F0) e^4 , eps<=EMAX
def gfit(EMAX,quart=True):
    rows=[];rhs=[]
    for th,D in DATA.items():
        c6=np.cos(np.radians(6*th)); c12=np.cos(np.radians(12*th))
        for e,P in zip(D['eps'],D['Pn']):
            if 0<e<=EMAX+1e-9:
                r=[e,e**2,c6*e**2,e**3,c6*e**3,c12*e**3]+([e**4,c6*e**4] if quart else []); rows.append(r); rhs.append(P)
    c,res,*_=np.linalg.lstsq(np.array(rows),np.array(rhs),rcond=None)
    pred=np.array(rows)@c; rms=np.sqrt(np.mean((pred-rhs)**2)); return c,rms
fits={}
for EMAX,q in [(0.12,False),(0.15,True),(0.18,True),(0.21,True)]:
    c,rms=gfit(EMAX,q); fits[f'{EMAX}_{q}']=(c,rms); print('global fit',EMAX,q,'Y=%.1f D0=%.0f D6=%.0f E0=%.0f E6=%.0f E12=%.0f'%tuple(c[:6]),'rms=%.3f'%rms)
# ---- Poisson: nu(e) = nu1 + nu2(th) e + nu3 e^2 global
rows=[];rhs=[]
for th,D in DATA.items():
    c6=np.cos(np.radians(6*th))
    for e,et in zip(D['eps'],D['et']):
        if 0<e<=0.12+1e-9: rows.append([1,e,c6*e,e**2,c6*e**2]); rhs.append(-et/e)
cn,*_=np.linalg.lstsq(np.array(rows),np.array(rhs),rcond=None); print('nu fit nu1=%.4f nu20=%.3f nu26=%.3f nu30=%.2f nu36=%.2f'%tuple(cn))
# ---- shear gamma/e^2 vs sin6th
for th,D in DATA.items():
    m=(D['eps']>0)&(D['eps']<=0.12); print('theta',th,'gamma/e^2 mean',np.round(np.mean(D['g'][m]/D['eps'][m]**2),4) if m.any() else None)
json.dump(dict(ideal=out,fits={k:(list(map(float,v[0])),float(v[1])) for k,v in fits.items()},nu=list(map(float,cn))),open('summary.json','w'),indent=1)
print(json.dumps(out,indent=0))
