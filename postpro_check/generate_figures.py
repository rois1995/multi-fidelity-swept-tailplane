"""Rebuild manuscript coefficient plots from the existing CSVs; run from any directory."""
from pathlib import Path
import csv
import hashlib
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/tailplane-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from audit import C, series, load, np, ROOT

OUT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Serif', 'mathtext.fontset':'stix',
    'font.size':9, 'pdf.fonttype':42, 'axes.spines.top':False,
    'axes.spines.right':False})
STYLE = {'WT':dict(color='black', marker='o', ms=2.5, lw=1),
         'RANS':dict(color='#1f5fa8', ls='-', lw=1.3),
         'LL':dict(color='#d0711c', ls='--', lw=1.3)}
LABEL = {'WT':'Wind tunnel', 'RANS':'RANS', 'LL':'DUST lifting line'}
KEYS = ['CL', 'CD', 'CM']
YLAB = [r'$C_L$', r'$C_D$', r'$C_{M_y}$']
ledger = []
sources = ['WT/WTCoeffs_DA_0_LE_1.csv', 'RANS/Coeffs_DA_0.csv',
           'DUST/LIFTINGLINE/Coeffs_DA_0.csv']

def draw(ax, data, name, figure, panel, delta, uncertainty=None, scale=1):
    a = np.array(sorted(data)); v = np.array([data[k] for k in a])*scale
    if uncertainty is not None:
        e = np.array([uncertainty[k] for k in a])
        assert np.all(e >= 0)
        line = ax.errorbar(a, v, e, capsize=1.4, elinewidth=.6, **STYLE[name]).lines[0]
    else:
        e = np.full(len(a), np.nan)
        line, = ax.plot(a, v, **STYLE[name])
    assert np.array_equal(line.get_xdata(), a)
    assert np.array_equal(line.get_ydata(), v, equal_nan=True)
    for aa, vv, ee in zip(a, v, e):
        ledger.append(dict(figure=figure, panel=panel, method=name, delta=delta,
                           alpha=aa, value=vv, uncertainty=ee))

def frame(ax):
    ax.set_xlim(-26,26); ax.set_xticks([-20,-10,0,10,20])
    ax.grid(True, lw=.35, color='0.8'); ax.axhline(0, color='0.5', lw=.5)

def finish(fig, name, handles=None):
    if handles is None:
        handles = [Line2D([], [], label=LABEL[n], **STYLE[n]) for n in STYLE]
    fig.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False,
               bbox_to_anchor=(.5,.005))
    fig.tight_layout(rect=(0,.075,1,1))
    fig.savefig(OUT/(name+'.pdf'))
    fig.savefig(OUT/'postpro_check'/(name+'.png'), dpi=140)
    plt.close(fig)

def title(d):
    return rf'$\delta={d}^\circ$' if d != 34 else '$\\delta=34^\\circ$ (WT)\n$35^\\circ$ (num.)'

for name, deltas in [('Coefficients_baseline',[0]),
                     ('Coefficients_selected',[0,10,20,34]),
                     ('Coefficients_complete',[0,10,20,30,34])]:
    fig, axes = plt.subplots(3,len(deltas),figsize=(7,5.5),sharex=True,sharey='row',squeeze=False)
    for j,d in enumerate(deltas):
        for i,key in enumerate(KEYS):
            ax=axes[i,j]; frame(ax)
            for n in STYLE:
                dn=35 if d==34 and n!='WT' else d
                err=C('WT',d,{'CL':'sCL','CD':'sigmaCD','CM':'sCM'}[key]) if n=='WT' else None
                draw(ax,C(n,dn,key),n,name,key,dn,err)
            if j==0: ax.set_ylabel(YLAB[i])
            if i==0: ax.set_title(title(d),fontsize=9)
            if i==2: ax.set_xlabel(r'$\alpha$ [deg]')
    finish(fig,name)

# Same finite differences and uncertainty propagation as the existing Fig. 7 renderer.
fig, axes = plt.subplots(2,4,figsize=(7,3.9),sharex=True,sharey='row')
intervals=[(0,10),(10,20),(20,30),(30,34)]
for j,(d0,d1) in enumerate(intervals):
    for i,key in enumerate(['CL','CM']):
        ax=axes[i,j];frame(ax)
        if j<2: ax.axvspan(*[(-6,0),(-14,-4)][j],color='0.88',zorder=0)
        for n in STYLE:
            dn=35 if d1==34 and n!='WT' else d1
            c0,c1=C(n,d0,key),C(n,dn,key)
            diff={a:(c1[a]-c0[a])/np.deg2rad(dn-d0) for a in c0 if a in c1}
            err=None
            if n=='WT':
                sk='sCL' if key=='CL' else 'sCM';s0,s1=C(n,d0,sk),C(n,d1,sk)
                err={a:np.hypot(s0[a],s1[a])/np.deg2rad(d1-d0) for a in diff}
            draw(ax,diff,n,'ControlEffectiveness_revised',key,f'{d0}:{dn}',err)
        if i==0:
            ax.set_title(rf'$\delta\in[{d0}^\circ,{d1}^\circ]$'+(' (WT)\n$[30^\circ,35^\circ]$ (num.)' if j==3 else ''),fontsize=9)
        if i==1:ax.set_xlabel(r'$\alpha$ [deg]')
axes[0,0].set_ylabel(r'$C_{L,\delta}$ [rad$^{-1}$]');axes[0,0].set_ylim(-2,3.5)
axes[1,0].set_ylabel(r'$C_{M_y,\delta}$ [rad$^{-1}$]');axes[1,0].set_ylim(-1.6,.6)
finish(fig,'ControlEffectiveness_revised',
       [Line2D([],[],label=LABEL[n],**STYLE[n]) for n in STYLE]+[Patch(color='0.88',label=r'Averaging range $\mathcal{A}$')])

# Preserve the patch accounting used in comparison_DistinctSurfaces.m.
patches={p:load(f'RANS/Coeffs_DA_0_{p}.csv') for p in
         ['Total','UpWing','DownWing','UpFlap','DownFlap','FlapJunction','BluntTE','Tip']}
sources += [f'RANS/Coeffs_DA_0_{p}.csv' for p in patches]
fig, axes=plt.subplots(2,2,figsize=(7,4.4),sharex=True,sharey='row')
parts=[('Total',['Total']),('Fixed surface',['UpWing','DownWing']),
       ('Control surface',['UpFlap','DownFlap','FlapJunction','BluntTE'])]
partstyles=[dict(color='black',lw=1.3),dict(color='#1f5fa8',ls='--',lw=1.3),dict(color='#d0711c',ls=':',lw=1.3)]
for j,d in enumerate([20,35]):
    for i,key in enumerate(['CL','CD']):
        ax=axes[i,j];frame(ax)
        for (label,ps),style in zip(parts,partstyles):
            ds=[series(patches[p],d,key) for p in ps]
            a=np.array(sorted(set.intersection(*(set(z) for z in ds))))
            v=np.array([sum(z[aa] for z in ds) for aa in a])
            if label=='Control surface' and key=='CL':v*=4
            line,=ax.plot(a,v,**style)
            assert np.array_equal(line.get_ydata(),v,equal_nan=True)
            for aa,vv in zip(a,v):ledger.append(dict(figure='Load_decomposition',panel=key,method=label,delta=d,alpha=aa,value=vv,uncertainty=np.nan))
        if i==0:ax.set_title(rf'$\delta={d}^\circ$')
        if j==0:ax.set_ylabel(YLAB[i])
        if i==1:ax.set_xlabel(r'$\alpha$ [deg]')
finish(fig,'Load_decomposition',[Line2D([],[],label=label+(' ($4C_L$)' if label=='Control surface' else ''),**style) for (label,_),style in zip(parts,partstyles)])

# Numerical regression checks against the manuscript and an independently written difference.
assert np.isclose(C('WT',0,'CL')[12],.7195,atol=5e-5)
assert np.isclose(np.mean([(C('WT',10,'CL')[a]-C('WT',0,'CL')[a])/np.deg2rad(10) for a in [-6,-4,-2,0]]),2.124,atol=.0005)
assert np.isclose(C('RANS',0,'CL')[12],.7875,atol=5e-5)
assert all(any(r['figure']==f and r['method']==n for r in ledger) for f in ['Coefficients_baseline','Coefficients_selected','Coefficients_complete','ControlEffectiveness_revised'] for n in STYLE)
with (OUT/'postpro_check/figure_data.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(ledger[0]));writer.writeheader();writer.writerows(ledger)
(OUT/'postpro_check/figure_sources.json').write_text(json.dumps({p:hashlib.sha256(Path(ROOT,p).read_bytes()).hexdigest() for p in sources},indent=2)+'\n')
print('Generated five figures; data identity and coefficient regression checks passed.')
