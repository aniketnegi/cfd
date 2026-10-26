"""Plot only converged radial states; do not label them free-axial equilibria."""
from pathlib import Path
import csv,hashlib,json,sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
out=Path(sys.argv[1])
dest=Path(sys.argv[2]) if len(sys.argv)>2 else out
if dest!=out:dest.mkdir(parents=True,exist_ok=False)
vtk_base=os.path.relpath(out,dest)
outputs=('load-balance.png','operating-branch.png','fixed-seat-fields.png','gallery.html','plot.json')
if any((dest/f).exists() for f in outputs):raise SystemExit('refusing to overwrite fitted load plots')
rows=list(csv.DictReader((out/'states.tsv').open(),delimiter='\t'))
fig,axes=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
for l,col in zip((1,2,4),('#78aaa6','#24738d','#182e4e')):
    a=sorted([r for r in rows if r['level']==str(l) and r['kind']=='branch' and r['radial_converged']=='true'],key=lambda r:float(r['clearance_um']))
    x=[float(r['clearance_um']) for r in a]
    for ax,key in zip(axes,('Fz_N','G_N','eccentricity_ratio')):
        ax.plot(x,[float(r[key]) for r in a],'o' if l==4 else 'o-',color=col,label=f'level {l}: {8416*l*l:,} cells')
axes[0].axhline(200,color='#b23842',ls='--',label='axial target 200 N')
axes[1].axhline(0,color='#b23842',ls='--')
for ax,label in zip(axes,('Axial film force [N]','Unbalanced axial force Fz − 200 [N]','Eccentricity / operating clearance')):
    ax.set(xlabel='Operating radial clearance [µm]',ylabel=label);ax.grid(alpha=.2)
axes[0].legend(fontsize=8)
fig.suptitle('5 mm fitted inlet: radial load balanced at each point; axial balance requires a zero crossing\nFinest mesh shown as markers only; connecting lines are guides, not resolved continuous branches')
fig.savefig(dest/'load-balance.png',dpi=170,bbox_inches='tight');plt.close(fig)

fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
keys=('eccentricity_um','min_normal_gap_um','net_feed_L_min','land_peak_MPa','radial_residual_N','minimum_fill')
labels=('Transverse offset [µm]','Minimum normal gap [µm]','Net supply [L/min]','Peak land-cell pressure [MPa]','Radial residual magnitude [N]','Minimum liquid content')
for l in (1,2,4):
    a=sorted([r for r in rows if r['level']==str(l) and r['kind']=='branch' and r['radial_converged']=='true'],key=lambda r:float(r['clearance_um']))
    for ax,key,label in zip(axes.flat,keys,labels):
        ax.plot([float(r['clearance_um']) for r in a],[float(r[key]) for r in a],'o' if l==4 else 'o-',label=f'level {l}')
        ax.set(xlabel='Operating radial clearance [µm]',ylabel=label);ax.grid(alpha=.2)
axes.flat[0].legend();fig.suptitle('Rigid, isoviscous JFO; fixed 0.5 MPa mouth pressure, 2000 rpm, μ = 0.0277 Pa·s')
fig.savefig(dest/'operating-branch.png',dpi=170,bbox_inches='tight');plt.close(fig)

field=out/'fixed-seat-fields.tsv'
field_vtk=None
if field.exists():
    a=np.genfromtxt(field,names=True,delimiter='\t')
    state=next(r for r in rows if r['level']=='4' and r['kind']=='branch' and np.isclose(float(r['clearance_um']),50))
    field_vtk=out/(state['label']+'.vtk');lines=field_vtk.read_text().splitlines()
    i=next(i for i,s in enumerate(lines) if s.startswith('POINTS '));nn=int(lines[i].split()[1])
    nodes=np.array([np.fromstring(s,sep=' ') for s in lines[i+1:i+1+nn]])
    i=next(i for i,s in enumerate(lines) if s.startswith('CELLS '));nc=int(lines[i].split()[1])
    quads=np.array([np.fromstring(s,sep=' ',dtype=int)[1:] for s in lines[i+1:i+1+nc]])
    theta=np.mod(np.rad2deg(np.arctan2(nodes[:,1],nodes[:,0])),360)[quads]
    seam=np.ptp(theta,axis=1)>180
    theta[seam]=np.where(theta[seam]<180,theta[seam]+360,theta[seam])
    polygons=np.stack((theta,nodes[quads,2]),axis=-1)
    assert nc==len(a)
    fig,axes=plt.subplots(3,1,figsize=(12,10),layout='constrained')
    for ax,key,title,cmap in zip(axes,('pressure_MPa','fill','normal_gap_um'),('Land-cell pressure [MPa]','Liquid content','Normal gap [µm]'),('viridis','cividis','magma')):
        sc=PolyCollection(polygons,array=a[key],edgecolors='none',antialiased=False,cmap=cmap,rasterized=True)
        ax.add_collection(sc)
        ax.set(xlim=(0,360),ylim=(0,100),xlabel='Circumferential angle [degrees]',ylabel='z [mm]')
        fig.colorbar(sc,ax=ax,label=title)
    fig.suptitle('Finest fixed-seat radial equilibrium: c = 50 µm; NOT simultaneous axial equilibrium\nActual exported land-cell polygons; reservoir interior omitted')
    fig.savefig(dest/'fixed-seat-fields.png',dpi=150,bbox_inches='tight');plt.close(fig)
html='<!doctype html><meta charset="utf-8"><title>Fitted inlet load balance</title><style>body{max-width:1500px;margin:2em auto;font:18px sans-serif;line-height:1.4}img{width:100%}</style><h1>5 mm fitted-inlet load continuation</h1><p>Radial equilibrium at prescribed axial positions. A free-axial equilibrium requires all three residuals to close. A missing sign crossing at sampled points is not a global nonexistence or stability result. Supply, viscosity, loads and mouth diameter are unchanged.</p>'
for name in outputs[:3]:
    if (dest/name).exists():html+=f'<h2>{name[:-4]}</h2><a href="{name}"><img src="{name}"></a>'
html+='<h2>Actual pressure/content/gap exports</h2><ul>'
for r in rows:
    if (out/(r['label']+'.vtk')).exists():html+=f'<li><a href="{vtk_base}/{r["label"]}.vtk">{r["label"]}</a>: Fz={float(r["Fz_N"]):.3f} N, axial residual={float(r["G_N"]):.3f} N</li>'
html+='</ul><p>Gross rim exchange is not pump-line injection/return. VTK coordinates are housing-local millimetres; pressure is Pa, gap m.</p>'
(dest/'gallery.html').write_text(html)
inputs=[out/'states.tsv']+([field,field_vtk] if field.exists() else [])
(dest/'plot.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    'output_sha256':{n:hashlib.sha256((dest/n).read_bytes()).hexdigest() for n in outputs[:-1] if (dest/n).exists()}},indent=2)+'\n')
