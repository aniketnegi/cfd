"""Side-by-side size and worst-junction views, using geometry study data only."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

out=Path(sys.argv[1])
names=('diameters.png','worst-faces.png','quality.png','operator.png','index.html','plot.json')
if any((out/n).exists() for n in names):
    raise SystemExit('refusing to overwrite diameter pictures')
rows=list(csv.DictReader((out/'quality.tsv').open(),delimiter='\t'))
diameters=sorted({int(r['diameter_mm']) for r in rows})

def row(d,l):
    return next(r for r in rows if int(r['diameter_mm'])==d and int(r['level'])==l)

def patch(d,l):
    a=np.genfromtxt(out/f'patch-{d}-{l}.tsv',names=True,delimiter='\t')
    xy=np.stack([np.stack([a[f'x{i}_mm'],a[f'z{i}_mm']],axis=-1) for i in range(1,5)],axis=1)
    return a,xy

def panel(ax,d,l,field,limits):
    a,xy=patch(d,l)
    cx,cy=xy.mean(axis=1).T
    keep=(cx>limits[0][0]-3)&(cx<limits[0][1]+3)&(cy>limits[1][0]-3)&(cy<limits[1][1]+3)
    pc=PolyCollection(xy[keep],array=a[field][keep],cmap='viridis',clim=(.8,1),edgecolors='#29383d',
                      linewidths=.35,antialiased=False,rasterized=True)
    ax.add_collection(pc)
    t=np.linspace(0,2*np.pi,256)
    ax.plot(d/2*np.cos(t),d/2*np.sin(t),color='#dc267f',lw=1.2)
    ax.set(xlim=limits[0],ylim=limits[1],aspect='equal',xlabel='X [mm]',ylabel='z − midpoint [mm]')
    return pc

fig,axes=plt.subplots(2,4,figsize=(18,10),layout='constrained')
for ax,d in zip(axes.flat,diameters):
    pc=panel(ax,d,1,'corner',((-45,45),(-45,45)))
    r=row(d,4)
    accepted=all(row(d,l)['quality_pass']=='true' for l in (1,2,4))
    ax.set_title(f'D = {d} mm · '+('all levels pass' if accepted else 'quality gate fails')+
                 f"\nfinest corner {float(r['minimum_corner']):.3f}, alignment {float(r['minimum_alignment']):.3f}",fontsize=11)
axes.flat[-1].axis('off')
axes.flat[-1].text(.02,.9,'Same bearing and cell budgets\n\nCoarse overview: 8,416 land cells\nRim fits the actual drilling circle\n+Y projection only\n\nColours: surface corner score\nNOT 3D Hex8 OQ\n\nNo load/pressure tuning',va='top',fontsize=12)
fig.colorbar(pc,ax=list(axes.flat[:-1]),shrink=.7,label='Surface corner score')
fig.suptitle('One diameter-relative construction rule — physical axes shared across all panels',fontsize=16)
fig.savefig(out/'diameters.png',dpi=150);plt.close(fig)

fig,axes=plt.subplots(2,4,figsize=(18,10),layout='constrained')
for ax,d in zip(axes.flat,diameters):
    r=row(d,4);x=float(r['worst_face_x_mm']);z=float(r['worst_face_z_from_mid_mm'])
    pc=panel(ax,d,4,'alignment',((x-3.5,x+3.5),(z-3.5,z+3.5)))
    ax.plot(x,z,'o',mfc='none',mec='#dc267f',ms=9,mew=1.5)
    ax.set_title(f"D = {d} mm · worst face {float(r['minimum_alignment']):.3f}\nmax neighbour-area ratio {float(r['max_area_ratio']):.3f}",fontsize=11)
axes.flat[-1].axis('off')
axes.flat[-1].text(.02,.9,'Finest level: 134,656 land cells\n\nEach panel follows its own\nworst-aligned face (magenta).\n\nSame 7 × 7 mm window size.\nColours: worst face alignment\namong each cell’s faces.\n\nNo claim of 3D OQ.',va='top',fontsize=12)
fig.colorbar(pc,ax=list(axes.flat[:-1]),shrink=.7,label='Cell minimum face alignment')
fig.suptitle('Close-ups of the actual worst faces — not selectively chosen attractive patches',fontsize=16)
fig.savefig(out/'worst-faces.png',dpi=150);plt.close(fig)

fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
for l,c in zip((1,2,4),('#92b9bd','#258d9b','#183e52')):
    for ax,key,title,bound in zip(axes,('minimum_corner','minimum_alignment','max_area_ratio'),
                                  ('Minimum corner ↑','Minimum face alignment ↑','Maximum neighbour area ratio ↓'),(.9,.9,2)):
        ax.plot(diameters,[float(row(d,l)[key]) for d in diameters],'o-',color=c,label=f'level {l}')
        ax.set(xlabel='Drill diameter [mm]',title=title,xticks=diameters);ax.grid(alpha=.2)
for ax,bound in zip(axes,(.9,.9,2)):
    ax.axhline(bound,color='#bd403a',ls='--',label='unchanged gate')
axes[-1].legend(fontsize=9)
fig.suptitle('Same quality gates for every diameter; positive/conforming cells checked separately')
fig.savefig(out/'quality.png',dpi=170);plt.close(fig)

if (out/'operator.tsv').exists():
    fields=list(csv.DictReader((out/'operator.tsv').open(),delimiter='\t'))
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for d in diameters:
        a=[r for r in fields if int(r['diameter_mm'])==d and r['kind']=='pressure']
        axes[0].loglog([int(r['level']) for r in a],[float(r['pressure_error']) for r in a],'o-',label=f'{d} mm')
        a=[r for r in fields if int(r['diameter_mm'])==d and r['kind']=='JFO' and float(r['rpm'])==2000]
        axes[1].loglog([int(r['level']) for r in a],[float(r['pressure_error']) for r in a],'o-')
        axes[2].loglog([int(r['level']) for r in a],[float(r['fill_error']) for r in a],'o-')
    for ax,title in zip(axes,('Smooth pressure: relative error','Curved JFO pressure error (+rotation)','Curved JFO content error (+rotation)')):
        ax.set(title=title,xlabel='Refinement level',xticks=(1,2,4));ax.grid(alpha=.2,which='both')
    axes[0].legend(fontsize=8,ncol=2)
    fig.suptitle('Manufactured fields, not physical diameter/load predictions; both rotation signs tested')
    fig.savefig(out/'operator.png',dpi=170);plt.close(fig)

html='<!doctype html><meta charset="utf-8"><title>Inlet diameter study</title><style>body{max-width:1600px;margin:2em auto;font:17px sans-serif}img{width:100%}</style><h1>4–10 mm inlet diameter study</h1><p>One construction rule, fixed bearing and cell budgets. Surface metrics, NOT Hex8 OQ. No physical operating-point or equilibrium predictions.</p>'
for name in names[:4]:
    if (out/name).exists():html+=f'<h2>{name[:-4]}</h2><a href="{name}"><img src="{name}"></a>'
(out/'index.html').write_text(html)
inputs=[out/'quality.tsv']+list(out.glob('patch-*.tsv'))
if (out/'operator.tsv').exists():inputs.append(out/'operator.tsv')
(out/'plot.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}},indent=2)+'\n')
