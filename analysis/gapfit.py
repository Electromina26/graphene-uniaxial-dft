import os
import json,numpy as np
from analyze import analyze
from scipy.optimize import least_squares
TH=[0,10,20,30]; rows=[]
for th in TH:
    D=analyze(th); G=json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'dft','gap',f'trim_th{th}.json')))
    for k,e in enumerate(D['eps']):
        key=f'{e:.3f}'
        if key in G:
            g=min(x[1] for x in G[key]['pigap']); s=D['bs'][k]
            sig=[v for n,l in G[key]['sigma_near_EF'] for v in l]
            rows.append((th,e,g,s,min(sig) if sig else None))
# signed gap for zigzag beyond crossing: inversion -> negative
def signed(th,e,g): return -g if (th==0 and e>0.256) else g
def model(p,s): t0,b=p; t=t0*np.exp(-b*np.asarray(s)); return 2*min(abs(t[i]-t[j]-t[l]) for i,j,l in [(0,1,2),(1,0,2),(2,0,1)])
def smodel(p,s):
    t0,b=p; t=np.sort(t0*np.exp(-b*np.asarray(s))); return 2*(t[0]+t[1]-t[2])
res=least_squares(lambda p:[smodel(p,s)-signed(th,e,g) for th,e,g,s,_ in rows],[2.7,3.37])
print('fit t0=%.3f beta=%.3f'%tuple(res.x),'rms=%.3f'%np.sqrt(np.mean(res.fun**2)))
crit={}
for th in TH:
    R=[(e,signed(th,e,g),sg) for t,e,g,s,sg in rows if t==th]; R.sort()
    e=np.array([r[0] for r in R]); g=np.array([r[1] for r in R])
    # zero crossing: linear interp if crosses else quadratic extrapolation of last 4
    i=np.where(np.diff(np.sign(g))<0)[0]
    if len(i): j=i[0]; ec=e[j]-g[j]*(e[j+1]-e[j])/(g[j+1]-g[j]); how='interp'
    else:
        c=np.polyfit(e[-4:],g[-4:],2); r=np.roots(c); r=r[np.isreal(r)].real; r=r[r>e[-1]]; ec=float(r.min()) if len(r) else None; how='extrap'
    print(th,'eps_gap',ec,how,[(round(a,3),round(b,3),c) for a,b,c in R])
    crit[th]=dict(eps_gap=ec,how=how,data=[(float(a),float(b),c) for a,b,c in R])
json.dump(dict(fit=list(res.x),crit=crit),open('gapsummary.json','w'),indent=1)
