import json, numpy as np, glob, os
from numpy.polynomial import polynomial as P
a0=2.4669; d0=a0/np.sqrt(3); beta=3.37; t0=2.7
A0=a0*a0*np.sqrt(3)/2  # area of primitive cell (Å^2)
ref_bonds=np.array([[a0/2, a0*np.sqrt(3)/6],[a0/2,-a0*np.sqrt(3)/6],[0, a0/np.sqrt(3)]])  # from A: B, B-a2+a1.., etc (directions 30,-30,90 deg)
def load(th):
    fn=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'dft','run',f'th{th}','results.json')
    if not os.path.exists(fn): return None
    R=sorted(json.load(open(fn)),key=lambda r:r['eps'])
    return R
def bonds(r):
    cell=np.array(r['cell'])[:2,:2]; pos=np.array(r['pos'])[:,:2]
    A,B=pos[0],pos[1]; L=[]
    for i in range(-2,3):
        for j in range(-2,3):
            v=B+i*cell[0]+j*cell[1]-A; L.append(v)
    L=np.array(L); n=np.linalg.norm(L,axis=1); idx=np.argsort(n)[:3]
    return L[idx]
def analyze(th):
    R=load(th)
    if R is None: return None
    thr=np.radians(th); n=np.array([np.cos(thr),np.sin(thr)])
    eps=np.array([0]+[r['eps'] for r in R]); et=np.array([0]+[r['et'] for r in R]); g=np.array([0]+[r['g'] for r in R])
    snn=np.array([0]+[r['snn'] for r in R]); E=np.array([r['E'] for r in R])
    Pn=snn*(1+et)  # nominal stress, N/m (force per reference width)
    # bond strains (match to reference bonds by orientation)
    bs=[np.zeros(3)]; ba=[np.zeros(3)]
    for r in R:
        b=bonds(r); U=np.eye(2)+r['eps']*np.outer(n,n)+r['et']*np.outer([-n[1],n[0]],[-n[1],n[0]])+r['g']*(np.outer(n,[-n[1],n[0]])+np.outer([-n[1],n[0]],n))
        aff=(U@ref_bonds.T).T
        # match each relaxed bond to affine bond with max |cos|
        out=np.zeros(3); used=set()
        for k,v in enumerate(aff):
            c=[abs(np.dot(v,w))/np.linalg.norm(v)/np.linalg.norm(w) if j not in used else -1 for j,w in enumerate(b)]
            j=int(np.argmax(c)); used.add(j); out[k]=np.linalg.norm(b[j])/d0-1
        bs.append(out); ba.append(np.linalg.norm(aff,axis=1)/d0-1)
    return dict(th=th,eps=eps,et=et,g=g,snn=snn,Pn=Pn,E=E,bs=np.array(bs),ba=np.array(ba))
def fitPoly(eps,y,emax,deg=3):
    m=(eps<=emax+1e-9)&(eps>0)
    X=np.vstack([eps[m]**k for k in range(1,deg+1)]).T
    c,*_=np.linalg.lstsq(X,y[m],rcond=None); return c
def tb_gap(bstr):
    t=t0*np.exp(-beta*np.asarray(bstr)); t=np.sort(t)
    return 2*max(0,t[2]-t[0]-t[1])  # gap (eV) when triangle inequality violated
if __name__=='__main__':
    for th in [0,10,20,30]:
        D=analyze(th)
        if D is None: continue
        print(f"\n== theta={th}  points={len(D['eps'])-1}")
        for k in range(len(D['eps'])):
            print(f"{D['eps'][k]:.3f} et={D['et'][k]:+.5f} g={D['g'][k]:+.5f} P={D['Pn'][k]:7.3f} bonds={np.round(D['bs'][k]*100,2)} aff={np.round(D['ba'][k]*100,2)} gapTB={tb_gap(D['bs'][k]):.3f}")
        for emax in [0.06,0.10,0.15]:
            if D['eps'].max()>=emax:
                c=fitPoly(D['eps'],D['Pn'],emax,3); print('fit P<=',emax,'Y,D,E=',np.round(c,1))
