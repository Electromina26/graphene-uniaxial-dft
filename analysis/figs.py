import numpy as np, json
from analyze import analyze
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
plt.rcParams.update({'font.size':10,'font.family':'serif','mathtext.fontset':'dejavuserif','axes.linewidth':0.8,
    'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':True,
    'legend.framealpha':0.92,'legend.edgecolor':'0.7','legend.fancybox':False,'legend.fontsize':8.5,'lines.markersize':4.5})
TH=[0,10,20,30]; D={t:analyze(t) for t in TH}
for t in TH:
    m=D[t]['eps']<=0.28
    for k in ['eps','et','g','snn','Pn','bs','ba']: D[t][k]=D[t][k][m]
col={0:'#1f4e79',10:'#2e8b57',20:'#d2691e',30:'#8b1a1a'}
lab={0:r'$\theta=0^\circ$ (zigzag)',10:r'$\theta=10^\circ$',20:r'$\theta=20^\circ$',30:r'$\theta=30^\circ$ (armchair)'}
mk={0:'o',10:'s',20:'^',30:'D'}
S=json.load(open('summary.json')); G=json.load(open('gapsummary.json'))
def plab(a,s,x=0.03,y=0.95): a.text(x,y,s,transform=a.transAxes,fontweight='bold',va='top')
XL=r'applied strain $\varepsilon$ (%)'

# ---------- Fig 1
fig,ax=plt.subplots(1,2,figsize=(7.4,3.3))
for t in TH:
    e=D[t]['eps']; x=np.linspace(0,e.max(),300)
    ax[0].plot(x*100,CubicSpline(e,D[t]['Pn'])(x),'-',lw=1,color=col[t])
    ax[0].plot(e*100,D[t]['Pn'],mk[t],color=col[t],mfc='w',mew=1.1,label=lab[t])
    em=S['ideal'][str(t)]; ax[0].plot(em['eps_m']*100,em['P_m'],'*',ms=12,color=col[t],mec='k',mew=0.5,zorder=5)
    m=e>0; ax[1].plot(e[m]*100,-D[t]['et'][m]/e[m],mk[t]+'-',color=col[t],mfc='w',mew=1.1,lw=1,label=lab[t])
x=np.linspace(0,0.115,50); ax[0].plot(x*100,336.5*x,':',color='0.4',lw=1.2,label=r'harmonic, $Y\varepsilon$')
ax[0].plot([],[],'*',ms=10,color='0.5',mec='k',mew=0.5,ls='none',label=r'maximum $(\varepsilon_m,P_m)$')
ax[0].set_xlabel(XL); ax[0].set_ylabel(r'nominal stress $P$ (N m$^{-1}$)'); ax[0].set_ylim(0,42)
ax[0].legend(loc='lower right',fontsize=7.8)
ax[1].set_xlabel(XL); ax[1].set_ylabel(r'contraction ratio $-\varepsilon_t/\varepsilon$'); ax[1].legend(loc='lower left',fontsize=7.8)
plab(ax[0],'(a)'); plab(ax[1],'(b)',x=0.88)
fig.tight_layout(); fig.savefig('fig/fig1.pdf'); fig.savefig('fig/fig1.png',dpi=200)

# ---------- Fig 2
fig,ax=plt.subplots(1,3,figsize=(7.6,3.6))
for t in [10,20]:
    e=D[t]['eps']; m=(e>=0.02)
    ax[0].plot(e[m]*100,D[t]['g'][m]/e[m]**2/np.sin(np.radians(6*t)),mk[t]+'-',color=col[t],mfc='w',mew=1.1,lw=1,label=lab[t])
ax[0].axhline(0,color='0.6',lw=0.6)
ax[0].set_xlabel(XL); ax[0].set_ylabel(r'$\gamma\,/\,(\varepsilon^2\sin6\theta)$'); ax[0].legend(loc='upper center',bbox_to_anchor=(0.5,-0.24),fontsize=7.3,frameon=False)
th=np.linspace(0,60,300)
for w,ls,l in [('0.15_True','--',r'fit window $\varepsilon\leq15\%$'),('0.18_True','-',r'fit window $\varepsilon\leq18\%$'),('0.21_True',':',r'fit window $\varepsilon\leq21\%$')]:
    c=S['fits'][w][0]; ax[1].plot(th,c[1]+c[2]*np.cos(np.radians(6*th)),ls,color='k',lw=1.1,label=l)
ax[1].set_xlim(0,60); ax[1].set_xticks([0,15,30,45,60]); ax[1].set_ylim(-870,-470)
ax[1].set_xlabel(r'loading angle $\theta$ (deg)'); ax[1].set_ylabel(r'$D(\theta)$ (N m$^{-1}$)'); ax[1].legend(loc='upper center',bbox_to_anchor=(0.5,-0.24),fontsize=7.3,frameon=False)
ax[1].text(1,-860,'ZZ',fontsize=8); ax[1].text(27,-860,'AC',fontsize=8); ax[1].text(57,-860,'ZZ',fontsize=8)
em=[S['ideal'][str(t)]['eps_m']*100 for t in TH]; pm=[S['ideal'][str(t)]['P_m'] for t in TH]
l1,=ax[2].plot(TH,em,'o-',color='C0',lw=1.1,label=r'ideal strain $\varepsilon_m$')
a2=ax[2].twinx(); l2,=a2.plot(TH,pm,'s--',color='C1',lw=1.1,label=r'ideal strength $P_m$')
ax[2].set_xlabel(r'loading angle $\theta$ (deg)'); ax[2].set_ylabel(r'$\varepsilon_m$ (%)',color='C0'); a2.set_ylabel(r'$P_m$ (N m$^{-1}$)',color='C1')
ax[2].set_xticks([0,10,20,30]); ax[2].legend(handles=[l1,l2],loc='upper center',bbox_to_anchor=(0.5,-0.24),fontsize=7.3,frameon=False)
plab(ax[0],'(a)'); plab(ax[1],'(b)',x=0.05); plab(ax[2],'(c)',x=0.84)
fig.tight_layout(); fig.savefig('fig/fig2.pdf'); fig.savefig('fig/fig2.png',dpi=200)

# ---------- Fig 3
fig,ax=plt.subplots(1,2,figsize=(7.4,3.9))
spec={0:[(0,'C0','o','bonds 1, 2 (inclined, $\\pm30^\\circ$)'),(2,'C2','s','bond 3 (normal to load)')],
      30:[(0,'C0','o','bond 1 (parallel to load)'),(1,'C1','^','bonds 2, 3 (at $60^\\circ$ to load)')]}
for a,t in zip(ax,[0,30]):
    e=D[t]['eps']*100
    for k,c,m_,l in spec[t]:
        a.plot(e,D[t]['bs'][:,k]*100,m_+'-',color=c,mfc='w',mew=1.1,lw=1.1,label=l+', relaxed')
        a.plot(e,D[t]['ba'][:,k]*100,'--',color=c,lw=1.1,label=l+', affine')
    a.axhline(0,color='0.6',lw=0.6); a.set_xlabel(XL); a.set_ylabel(r'bond strain $\delta_i$ (%)')
    a.legend(loc='upper center',bbox_to_anchor=(0.5,-0.2),ncol=1,fontsize=7.6,frameon=False)
plab(ax[0],'(a) zigzag'); plab(ax[1],'(b) armchair')
ax[0].set_ylim(-6,23)
ins=ax[0].inset_axes([0.10,0.53,0.30,0.22]); e=D[0]['eps']; m=e>0
ins.plot(e[m]*100,(D[0]['bs'][m,0]-D[0]['ba'][m,0])*100,'o-',ms=2.5,lw=0.9,color='C0')
ins.axhline(0,color='0.5',lw=0.5); ins.tick_params(labelsize=7,direction='in')
ins.set_title(r'$\delta_1-\delta_1^{\rm aff}$ (%)',fontsize=7.5,pad=2)
fig.tight_layout(); fig.savefig('fig/fig3.pdf'); fig.savefig('fig/fig3.png',dpi=200)

# ---------- Fig 4
fig,ax=plt.subplots(figsize=(5.0,3.7)); Z=G['crit']['0']['data']
e=np.array([r[0] for r in Z])*100; g=np.array([r[1] for r in Z]); sg=np.array([np.nan if r[2] is None else r[2] for r in Z])
ax.plot(e,g,'o-',color='C0',mfc='w',mew=1.2,lw=1.2,label=r'DFT: signed $\pi$ splitting $\Delta_\pi(M)$')
ax.plot(e,sg,'s--',color='C2',mfc='w',mew=1.2,lw=1.1,label=r'DFT: lowest $\sigma^*$ state at $M$ ($-E_F$)')
t0,b=G['fit']; ee=D[0]['eps']; m=ee>=0.12
for bb,tt,ls,l in [(3.37,2.7,':',r'TB, $\beta=3.37$, $t_0=2.7$ eV'),(b,t0,'-.',r'TB, fitted $\beta=%.2f$, $t_0=%.2f$ eV'%(b,t0))]:
    h=np.array([2*tt*(lambda q:q[0]+q[1]-q[2])(np.sort(np.exp(-bb*s))) for s in D[0]['bs']]); ax.plot(ee[m]*100,h[m],ls,color='k',lw=1.2,label=l)
emz=S['ideal']['0']['eps_m']*100; ax.axvline(emz,color='C3',lw=1.1); ax.text(emz-0.3,-1.35,r'ideal strain $\varepsilon_m$',color='C3',ha='right',fontsize=8.5)
egz=G['crit']['0']['eps_gap']*100; ax.plot(egz,0,'*',ms=12,color='C3',mec='k',mew=0.5,zorder=5); ax.annotate(r'$\varepsilon_{\rm gap}=%.1f\%%$'%egz,xy=(egz,0),xytext=(egz+0.4,0.75),fontsize=8.5,arrowprops=dict(arrowstyle='->',lw=0.8))
ax.axhline(0,color='0.6',lw=0.6); ax.set_xlabel(XL); ax.set_ylabel('energy (eV)'); ax.set_xlim(11,31); ax.set_ylim(-1.6,3.0)
ax.legend(loc='upper right',fontsize=7.6)
fig.tight_layout(); fig.savefig('fig/fig4.pdf'); fig.savefig('fig/fig4.png',dpi=200)

# ---------- Fig 5
fig,ax=plt.subplots(figsize=(5.0,3.9))
emv=np.array([S['ideal'][str(t)]['eps_m'] for t in TH])*100
cs=CubicSpline(TH,emv,bc_type='clamped'); x=np.linspace(0,30,200)
ax.fill_between(x,cs(x),33,color='#f2dede',lw=0)
ax.plot(x,cs(x),'-',color='C3',lw=1.4,label=r'ideal strain $\varepsilon_m(\theta)$ (DFT)'); ax.plot(TH,emv,'o',color='C3')
ax.plot([0],[G['crit']['0']['eps_gap']*100],'D',color='C0',ms=7,label=r'gap opening, DFT ($M$ point)')
for t in [10,20,30]: ax.annotate('',xy=(t,30.6),xytext=(t,27.1),arrowprops=dict(arrowstyle='->',color='C0',lw=1.1))
ax.plot([10,20,30],[27,27,27],'_',color='C0',ms=12,mew=2,ls='none',label=r'DFT: no gap up to $\varepsilon=27\%$')
ax.plot([0],[24.26],'v',color='0.4',mfc='w',mew=1.2,ms=7,label=r'TB ($\beta=3.37$), affine bonds')
ax.plot([0,10],[22.62,25.75],'^',color='k',mfc='w',mew=1.2,ms=7,label=r'TB ($\beta=3.37$), relaxed bonds')
ax.set_xlabel(r'loading angle $\theta$ (deg)'); ax.set_ylabel(r'strain $\varepsilon$ (%)'); ax.set_xlim(-2,32); ax.set_ylim(9,33)
ax.set_xticks([0,10,20,30]); ax.set_xticklabels(['0\n(ZZ)','10','20','30\n(AC)'])
ax.text(16,31.7,'beyond the stress maximum',fontsize=8.5,ha='center',color='#8b1a1a')
ax.text(22,12.0,'mechanically stable, gapless',fontsize=8.5,ha='center',color='0.3')
ax.legend(loc='lower left',fontsize=7.3)
fig.tight_layout(); fig.savefig('fig/fig5.pdf'); fig.savefig('fig/fig5.png',dpi=200)
