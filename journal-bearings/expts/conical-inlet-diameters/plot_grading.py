"""Visualize neighbour-area jumps directly from saved coarse chart connectivity."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

out=Path(sys.argv[1])
outputs=('grading-map.png','area-jump-closeup.png','grading-plot.json')
if any((out/n).exists() for n in outputs):
    raise SystemExit('refusing to overwrite grading pictures')
meta=tomllib.loads((out/'geometry.toml').read_text())
rows=list(csv.DictReader((out/'quality.tsv').open(),delimiter='\t'))
R=.05005; ca=np.cos(np.deg2rad(10)); period=2*np.pi*R
meshes={};inputs={}
for d in range(4,11):
    paths=[out/f'coarse-{d}-{suffix}.tsv' for suffix in ('nodes','quads')]
    for path in paths:
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest==meta['coarse_sha256'][path.name]
        inputs[path.name]=digest
    nodes=np.loadtxt(paths[0]);quads=np.loadtxt(paths[1],dtype=int)-1
    xy=nodes[quads].copy()
    seam=np.ptp(xy[:,:,0],axis=1)>period/2
    xy[seam,:,0]=np.where(xy[seam,:,0]<period/2,xy[seam,:,0]+period,xy[seam,:,0])
    local=xy-xy[:,0:1,:]
    nxt=np.roll(local,-1,axis=1)
    areas=.5*(local[:,:,0]*nxt[:,:,1]-local[:,:,1]*nxt[:,:,0]).sum(axis=1)
    assert (areas>0).all()
    owners={}
    for k,q in enumerate(quads):
        for a,b in zip(q,np.roll(q,-1)):
            owners.setdefault(tuple(sorted((a,b))),[]).append(k)
    ratios=np.ones(len(quads));rear_axial=[]
    for edge,cells in owners.items():
        assert len(cells) in (1,2)
        if len(cells)!=2:continue
        i,j=cells;r=max(areas[i],areas[j])/min(areas[i],areas[j])
        ratios[i]=max(ratios[i],r);ratios[j]=max(ratios[j],r)
        p=nodes[list(edge)]
        if abs(p[0,1]-p[1,1])<1e-12 and np.pi*R<p[:,0].mean()<period:
            rear_axial.append(r)
    reported=float(next(r for r in rows if int(r['diameter_mm'])==d and r['level']=='1')['max_area_ratio'])
    assert np.isclose(ratios.max(),reported,rtol=1e-8)
    predicted=9*d/(50-3*d);predicted=max(predicted,1/predicted)
    assert np.isclose(max(rear_axial),predicted,rtol=1e-8)
    meshes[d]=(xy,ratios,predicted)

fig,axes=plt.subplots(4,2,figsize=(15,13),layout='constrained')
for ax,d in zip(axes.flat,range(4,11)):
    xy,r,_=meshes[d];poly=xy.copy();poly[:,:,0]*=180/(np.pi*R);poly[:,:,1]*=ca*1000
    pc=PolyCollection(poly,array=r,cmap='magma',clim=(1,4.5),edgecolors='none',antialiased=False)
    ax.add_collection(pc)
    ax.set(xlim=(0,360),ylim=(0,100),xlabel='Circumferential angle [degrees]',ylabel='z [mm]',
           title=f'D = {d} mm · maximum ratio {r.max():.3f}')
    ax.axvline(90,color='#24b4b4',lw=.8,ls=':')
axes.flat[-1].axis('off')
axes.flat[-1].text(.03,.9,'Entire bearing, unwrapped · coarse level\n\nColour = largest area ratio to any neighbour.\nBrighter bands show abrupt size changes.\n\nDotted line: inlet angular position (+Y).\nThe rear half is 180–360°.\n\nThese are chart-area ratios, NOT Hex8 OQ.',va='top',fontsize=12)
fig.colorbar(pc,ax=list(axes.flat[:-1]),label='Neighbouring chart-area ratio',shrink=.8)
fig.suptitle('Where the size jumps occur — including the remote background, not just the inlet',fontsize=16)
fig.savefig(out/'grading-map.png',dpi=150,bbox_inches='tight',pad_inches=.15);plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(12,6),layout='constrained')
for ax,d in zip(axes,(5,10)):
    xy,r,ratio=meshes[d];poly=xy.copy()
    poly[:,:,0]=(poly[:,:,0]-R*3*np.pi/2)*1000
    step=50-3*d;poly[:,:,1]=poly[:,:,1]*ca*1000-step
    c=poly.mean(axis=1);keep=(abs(c[:,0])<10)&(abs(c[:,1])<10)
    pc=PolyCollection(poly[keep],array=r[keep],cmap='magma',clim=(1,4.5),edgecolors='#6f8a99',
                      linewidths=.8,antialiased=False)
    ax.add_collection(pc)
    ax.axhline(0,color='#2bbec5',ls='--',lw=1.4)
    ax.set(xlim=(-5,5),ylim=(-5,5),aspect='equal',xlabel='Reference arc distance from rear centre [mm]',
           ylabel=f'z − {step} mm',title=f'D = {d} mm · rear axial step ratio {ratio:.3f}\n'
           f'Inside spacing {3*d/8:.4f} mm / outside {(50-3*d)/24:.4f} mm')
    ax.text(.03,.96,'towards inlet mid-plane ↑',transform=ax.transAxes,va='top',color='white',
            bbox={'facecolor':'#222222','alpha':.8,'edgecolor':'none'})
fig.colorbar(pc,ax=axes,label='Neighbouring chart-area ratio',shrink=.75)
fig.suptitle('Same 10 × 10 mm windows on the rear side — why 10 mm reaches a 4.5:1 jump\n'
             'Coarse cells shown for clarity; uniform subdivision preserves this spacing ratio',fontsize=14)
fig.savefig(out/'area-jump-closeup.png',dpi=170,bbox_inches='tight',pad_inches=.15);plt.close(fig)
inputs['quality.tsv']=hashlib.sha256((out/'quality.tsv').read_bytes()).hexdigest()
(out/'grading-plot.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_sha256':inputs,'checks':'Recomputed coarse chart areas and adjacency; all seven maxima match quality.tsv; rear axial ratios match the spacing formula',
    'rear_axial_ratios':{str(d):m[2] for d,m in meshes.items()}},indent=2)+'\n')
