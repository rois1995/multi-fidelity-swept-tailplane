import sys; sys.path.insert(0, sys.argv[1])
exec(open(sys.argv[1]+'/table3.py').read().split("rngs={")[0])
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','mathtext.fontset':'cm','font.size':9})
out=sys.argv[2]
ints=[(0,10,0,10),(10,20,10,20),(20,30,20,30),(30,34,30,35)]   # (WT d0,d1, num d0,d1)
bands={0:(-10,-4),1:(-10,-4)}
fig,ax=plt.subplots(2,4,figsize=(7.0,3.9),sharex=True,sharey='row')
sty={'RANS':dict(color='#1f5fa8',ls='-',lw=1.3),'LL':dict(color='#d0711c',ls='--',lw=1.3)}
for j,(w0,w1,n0,n1) in enumerate(ints):
    for i,key in enumerate(['CL','CM']):
        A=ax[i,j]
        if j in bands: A.axvspan(*bands[j],color='0.88',zorder=0)
        D=deriv('WT',w0,w1,key); sk='sCL' if key=='CL' else 'sCM'
        s0=C('WT',w0,sk); s1=C('WT',w1,sk)
        a=np.array(sorted(D)); v=np.array([D[k] for k in a]); e=np.array([np.hypot(s0[k],s1[k]) for k in a])/np.deg2rad(w1-w0)
        A.errorbar(a,v,e,color='k',marker='o',ms=2.5,lw=1.0,capsize=1.5,elinewidth=0.6,label='Wind tunnel',zorder=3)
        for n in ['RANS','LL']:
            D=deriv(n,n0,n1,key); a=np.array(sorted(D)); A.plot(a,[D[k] for k in a],label='RANS' if n=='RANS' else 'DUST lifting line',**sty[n])
        A.axhline(0,color='0.5',lw=0.5); A.grid(True,lw=0.3,color='0.8')
        A.set_xlim(-26,26); A.set_xticks([-20,-10,0,10,20])
        if i==0:
            lab=f'$\\delta\\in[{w0}^\\circ,{w1}^\\circ]$' if j<3 else '$\\delta\\in[30^\\circ,34^\\circ]$ (WT)\n$[30^\\circ,35^\\circ]$ (num.)'
            A.set_title(lab,fontsize=9)
        if i==1: A.set_xlabel(r'$\alpha$ [deg]')
ax[0,0].set_ylabel(r'$C_{L,\delta}$ [rad$^{-1}$]'); ax[1,0].set_ylabel(r'$C_{M_y,\delta}$ [rad$^{-1}$]')
ax[0,0].set_ylim(-2.0,3.5); ax[1,0].set_ylim(-1.6,0.6)
h,l=ax[0,0].get_legend_handles_labels()
import matplotlib.patches as mp
h.append(mp.Patch(color='0.88')); l.append(r'Averaging range $\mathcal{A}$')
fig.legend(h,l,loc='lower center',ncol=4,frameon=False,bbox_to_anchor=(0.5,-0.01))
fig.tight_layout(rect=(0,0.06,1,1))
fig.savefig(out); print('saved',out)
