import os, numpy as np, json
os.makedirs('out',exist_ok=True)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.interpolate import CubicSpline
from analyze import analyze
plt.rcParams.update({'font.size':9,'font.family':'serif','mathtext.fontset':'dejavuserif','axes.linewidth':0.8,
                     'xtick.direction':'in','ytick.direction':'in'})
S=json.load(open('summary.json')); G=json.load(open('gapsummary.json'))
TH=[0,10,20,30]; D={t:analyze(t) for t in TH}
col={0:'#1f4e79',10:'#2e8b57',20:'#d2691e',30:'#8b1a1a'}
NAVY='#16324f'; RED='#c0392b'; BLUE='#1f6fb2'

fig=plt.figure(figsize=(8.0,3.55))
bg=fig.add_axes([0,0,1,1],zorder=-10); bg.axis('off'); bg.set_xlim(0,1); bg.set_ylim(0,1)
# ---------- header band
fig.text(0.5,0.955,'Graphene under true uniaxial stress: does strain open a band gap before the lattice fails?',
         ha='center',va='center',fontsize=10.5,weight='bold',color=NAVY)
# panel boxes
for x0 in [0.010,0.341,0.672]:
    bg.add_patch(FancyBboxPatch((x0,0.135),0.318,0.775,boxstyle='round,pad=0.004,rounding_size=0.012',fc='#f7f9fc',ec='#c9d3df',lw=0.9,zorder=0))
def ptitle(x,s): fig.text(x,0.875,s,ha='center',va='center',fontsize=9.2,weight='bold',color=NAVY)
ptitle(0.169,'1  Stress-controlled DFT'); ptitle(0.500,'2  Direction-dependent strength'); ptitle(0.831,'3  Band gap at the M point')

# ---------- panel 1: lattice
a=fig.add_axes([0.02,0.20,0.30,0.64]); a.set_aspect('equal'); a.axis('off')
a1=np.array([1.0,0]); a2=np.array([0.5,np.sqrt(3)/2]); B=(a1+a2)/3
for i in range(-4,5):
    for j in range(-4,5):
        A=i*a1+j*a2
        if np.linalg.norm(A)>2.35: continue
        for nb in [B,B-a1,B-a2]:
            P=A+nb; a.plot([A[0],P[0]],[A[1],P[1]],'-',color='0.55',lw=1.3,zorder=1)
        a.plot(*A,'o',ms=4.2,color='#1f4e79',zorder=3); a.plot(*(A+B),'o',ms=4.2,color='#8b1a1a',zorder=3)
# highlighted bonds around origin
A0=np.zeros(2)
for nb,c,l,pos in [(B,'C2',r'$\delta_1$',(0.22,0.36)),(B-a1,'C2',r'$\delta_2$',(-0.62,0.36)),(B-a2,'C1',r'$\delta_3$',(0.10,-0.45))]:
    a.plot([0,nb[0]],[0,nb[1]],'-',color=c,lw=3.0,zorder=2); a.text(pos[0],pos[1],l,fontsize=9.5,color=c,weight='bold',zorder=6,bbox=dict(fc='w',ec='none',pad=0.5,alpha=0.85))
th=np.radians(17); n=np.array([np.cos(th),np.sin(th)])
a.plot([-2.25*n[0],2.25*n[0]],[-2.25*n[1],2.25*n[1]],'--',color='#e67e22',lw=1)
for sg in [1,-1]: a.annotate('',xy=sg*2.95*n,xytext=sg*2.25*n,arrowprops=dict(arrowstyle='-|>,head_width=0.35,head_length=0.6',lw=2.4,color='#e67e22'))
a.plot([0,2.3],[0,0],'-',color='k',lw=1.1,zorder=4); a.text(2.35,-0.05,'ZZ',fontsize=8,va='center',zorder=6,bbox=dict(fc='w',ec='none',pad=0.4)); t=np.linspace(0,th,60); R=1.75; a.plot(R*np.cos(t),R*np.sin(t),'-',color='k',lw=1.3,zorder=5)
a.text(2.02*np.cos(th/2),2.02*np.sin(th/2),r'$\theta$',fontsize=12,ha='center',va='center',zorder=7,bbox=dict(fc='w',ec='none',pad=0.3,alpha=0.9)); a.text(2.88*n[0]-0.15,2.88*n[1]+0.22,r'$\sigma_{nn}$',fontsize=9,color='#e67e22')
a.text(0,-2.75,r'$\sigma_{tt}=\sigma_{nt}=0$, atoms relaxed'+'\n'+r'$\theta=0^\circ$ (zigzag) $\rightarrow$ $30^\circ$ (armchair)',
       ha='center',va='top',fontsize=8,color=NAVY)
a.set_xlim(-3.1,3.1); a.set_ylim(-3.5,2.5)

# ---------- panel 2: stress-strain
b=fig.add_axes([0.400,0.30,0.248,0.52])
for t,l in zip(TH,['0° ZZ','10°','20°','30° AC']):
    e=D[t]['eps']; m=e<=0.28; e=e[m]; P=D[t]['Pn'][m]; x=np.linspace(0,e.max(),200)
    b.plot(x*100,CubicSpline(e,P)(x),'-',lw=1.6,color=col[t],label=l)
    em=S['ideal'][str(t)]; b.plot(em['eps_m']*100,em['P_m'],'*',ms=10,color=col[t],mec='k',mew=0.5,zorder=5)
b.set_xlim(0,28); b.set_ylim(0,40); b.set_xlabel(r'strain $\varepsilon$ (%)',labelpad=1); b.set_ylabel(r'$P$ (N m$^{-1}$)',labelpad=1)
b.legend(loc='lower right',fontsize=7,frameon=False,handlelength=1.4,labelspacing=0.25)
b.tick_params(labelsize=8)
fig.text(0.500,0.175,r'$\varepsilon_m$: 23.9% (ZZ) $\rightarrow$ 18.5% (AC)'+'\n'+r'$D(\theta)=D_0+D_6\cos6\theta$,  $\gamma\propto\varepsilon^2\sin6\theta$',
         ha='center',va='center',fontsize=7.8,color=NAVY)

# ---------- panel 3: gap at M (zigzag)
c=fig.add_axes([0.727,0.30,0.255,0.52])
Z=G['crit']['0']['data']; e=np.array([r[0] for r in Z])*100; g=np.array([r[1] for r in Z])
emz=S['ideal']['0']['eps_m']*100; egz=G['crit']['0']['eps_gap']*100
c.axvspan(emz,31,color='#f6dcdc',lw=0); c.axhline(0,color='0.5',lw=0.7)
c.plot(e,g,'o-',color=BLUE,ms=3.6,lw=1.6,mfc='w',mew=1.1)
c.axvline(emz,color=RED,lw=1.3); c.plot(egz,0,'*',ms=12,color='#f1c40f',mec='k',mew=0.6,zorder=6)
c.text(emz-0.4,1.85,'lattice\nfails',color=RED,fontsize=7.8,ha='right',va='center')
c.text(28.0,1.85,'gap\nopens',color='#9a7d0a',fontsize=7.8,ha='center',va='center')
c.annotate('',xy=(egz,-0.55),xytext=(emz,-0.55),arrowprops=dict(arrowstyle='<->',lw=1.0,color='k'))
c.text((emz+egz)/2,-0.95,r'$+1.7\%$',ha='center',fontsize=8,weight='bold')
c.set_xlim(12,30.5); c.set_ylim(-1.25,2.6); c.set_xlabel(r'strain $\varepsilon$ (%), zigzag',labelpad=1); c.set_ylabel(r'$\Delta_\pi(M)$ (eV)',labelpad=1)
c.tick_params(labelsize=8)
fig.text(0.831,0.175,'Other directions: no gap up to 27%\n'+r'TB recalibrated to DFT: $\beta=2.89$',ha='center',va='center',fontsize=7.8,color=NAVY)

# arrows between panels
for x in [0.3255,0.6565]:
    bg.add_patch(FancyArrowPatch((x,0.53),(x+0.018,0.53),arrowstyle='-|>,head_width=4,head_length=5',color=NAVY,lw=1.5,zorder=20))
# bottom message
bg.add_patch(FancyBboxPatch((0.17,0.022),0.66,0.075,boxstyle='round,pad=0.004,rounding_size=0.015',fc=NAVY,ec='none'))
fig.text(0.5,0.060,'Lattice instability precedes gap opening for every loading direction',ha='center',va='center',fontsize=9.8,weight='bold',color='w')
fig.savefig('out/Graphical_abstract.pdf'); fig.savefig('out/Graphical_abstract.png',dpi=400); fig.savefig('out/Graphical_abstract.tif',dpi=400)
