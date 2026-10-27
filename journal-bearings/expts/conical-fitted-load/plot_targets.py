from pathlib import Path
import csv,hashlib,json,sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1])
dest=Path(sys.argv[2]) if len(sys.argv)>2 else out
if dest!=out:dest.mkdir(parents=True,exist_ok=False)
vtk_base=os.path.relpath(out,dest)
names=('axial-targets.png','branch-diagnostics.png','gallery.html','plot.json')
if any((dest/n).exists() for n in names):raise SystemExit('refusing to overwrite target plots')
read=lambda p:list(csv.DictReader(p.open(),delimiter='\t'))
rows=read(out/'states.tsv');prior=read(out/'prior-states.tsv');targets=read(out/'targets.tsv')
fig,axes=plt.subplots(1,2,figsize=(14,5),layout='constrained')
for l,col in ((1,'#70aaa7'),(2,'#193e60')):
    a={round(float(r['clearance_um']),8):r for r in prior+rows if int(r['level'])==l and r['radial_converged']=='true'}
    a=[a[k] for k in sorted(a)];x=[float(r['clearance_um']) for r in a]
    axes[0].plot(x,[float(r['Fz_N']) for r in a],'o' if l==2 else 'o-',color=col,label=f'level {l}, {8416*l*l:,} cells',ms=4)
    axes[1].plot(x,[float(r['eccentricity_ratio']) for r in a],'o' if l==2 else 'o-',color=col,ms=4)
for target,col in ((500,'#bb533b'),(750,'#9651a2'),(1000,'#9e7900')):
    axes[0].axhline(target,color=col,ls='--',label=f'{target} N axial target')
for r in targets:
    if r['converged']=='true':axes[0].plot(float(r['clearance_um']),float(r['Fz_N']),'*',ms=13,color='#d52f48')
axes[0].set(ylabel='Axial film force [N]');axes[1].set(ylabel='Eccentricity / operating clearance')
for ax in axes:ax.set(xlabel='Operating radial clearance [µm]');ax.grid(alpha=.2)
axes[0].legend(fontsize=8)
fig.suptitle('1000 N radial load; 500/750/1000 N axial targets; 0.5 MPa mouth pressure\nFormal model only: physical validity/stability unverified; level2 markers show only computed samples')
fig.savefig(dest/names[0],dpi=170,bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for l,col in ((1,'#70aaa7'),(2,'#193e60')):
    a=sorted([r for r in rows if int(r['level'])==l and r['radial_converged']=='true'],key=lambda r:float(r['clearance_um']))
    for ax,key,title in zip(axes.flat,('min_normal_gap_um','net_feed_L_min','max_normal_gap_over_mouth_diameter','radial_residual_N'),
                           ('Minimum normal gap [µm]','Net feed [L/min]','Maximum normal gap / 5 mm mouth diameter','Radial residual magnitude [N]')):
        ax.plot([float(r['clearance_um']) for r in a],[float(r[key]) for r in a],'o' if l==2 else 'o-',color=col,label=f'level {l}',ms=4)
        ax.set(xlabel='Operating radial clearance [µm]',ylabel=title);ax.grid(alpha=.2)
axes.flat[0].legend();fig.suptitle('Converged radial states only; thin-film validity is not certified by force balance')
fig.savefig(dest/names[1],dpi=170,bbox_inches='tight');plt.close(fig)
html='<!doctype html><meta charset="utf-8"><title>Axial load targets</title><style>body{max-width:1500px;margin:2em auto;font:18px sans-serif;line-height:1.5}img{width:100%}td,th{padding:.4em;border:1px solid #ccc}table{border-collapse:collapse}</style><h1>Axial loads 500, 750 and 1000 N</h1><p>Radial load remains 1000 N and mouth pressure 0.5 MPa. The earlier 12.5–100 µm branch is reused because the axial target does not enter the fixed-seat radial equations. A formal steady root is not proof of thin-film validity or dynamic stability.</p><table><tr><th>Level</th><th>Target N</th><th>Result</th><th>Clearance µm</th><th>Axial residual N</th></tr>'
for r in targets:
    html+=f'<tr><td>{r["level"]}</td><td>{r["target_axial_N"]}</td><td>{r["reason"]}</td><td>{r["clearance_um"]}</td><td>{r["residual_N"]}</td></tr>'
html+='</table>'
for name in names[:2]:html+=f'<h2>{name[:-4]}</h2><a href="{name}"><img src="{name}"></a>'
html+='<h2>Computed fields</h2><ul>'
for r in rows:
    if (out/(r['label']+'.vtk')).exists():html+=f'<li><a href="{vtk_base}/{r["label"]}.vtk">{r["label"]}</a>: Fz={float(r["Fz_N"]):.5f} N</li>'
html+='</ul><p>VTKs are housing-local mm; pressure Pa, gap m. Gross rim exchanges are not pump-line injection/return.</p>'
(dest/'gallery.html').write_text(html)
inputs=['states.tsv','prior-states.tsv','targets.tsv']
(dest/'plot.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_sha256':{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in inputs},
    'output_sha256':{n:hashlib.sha256((dest/n).read_bytes()).hexdigest() for n in names[:-1]}},indent=2)+'\n')
