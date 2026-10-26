"""Independent force reconstruction from exported Cartesian nodes/connectivity.
Uses polygon moments for Fz and Green's theorem edge integrals for Fx/Fy,
not the solver's 2x2 volume quadrature. Pressure is constant in each land cell.
"""
from pathlib import Path
import csv,hashlib,json,sys,tomllib
import numpy as np
out=Path(sys.argv[1]);target=out/'export-verification.json'
if target.exists():raise SystemExit('refusing to overwrite export verification')
rows=list(csv.DictReader((out/'states.tsv').open(),delimiter='\t'))
meta=tomllib.loads((out/'run.toml').read_text())
assert meta['status']=='completed'
for name,digest in meta['source_sha256'].items():
    assert hashlib.sha256((out/'source-snapshot'/name).read_bytes()).hexdigest()==digest,name
u,w=np.polynomial.legendre.leggauss(8);u=(u+1)/2;w=w/2
slope=np.tan(np.deg2rad(10.));results=[]
for row in rows:
    if row['radial_converged']!='true':continue # Failed candidates are not operating states.
    assert max(abs(float(row[k])) for k in ('Rx_N','Ry_N'))<=.02
    assert float(row['land_min_Pa'])>=-1e-5
    assert -1e-10<=float(row['minimum_fill'])<=float(row['maximum_fill'])<=1+1e-10
    assert max(float(row[k]) for k in ('local_flux_error','global_imbalance'))<1e-9
    q=float(row['net_feed_L_min'])
    assert np.isclose(float(row['rim_out_L_min'])-float(row['rim_in_L_min']),q,rtol=1e-9,atol=1e-13)
    assert np.isclose(float(row['large_exit_L_min'])+float(row['small_exit_L_min']),q,rtol=1e-9,atol=1e-13)
    path=out/(row['label']+'.vtk');lines=path.read_text().splitlines()
    i=next(i for i,s in enumerate(lines) if s.startswith('POINTS '));nn=int(lines[i].split()[1])
    nodes=np.array([np.fromstring(s,sep=' ') for s in lines[i+1:i+1+nn]])*.001
    i=next(i for i,s in enumerate(lines) if s.startswith('CELLS '));nc=int(lines[i].split()[1])
    quads=np.array([np.fromstring(s,sep=' ',dtype=int)[1:] for s in lines[i+1:i+1+nc]])
    i=next(i for i,s in enumerate(lines) if s.startswith('SCALARS pressure_pa '))
    p=np.array(lines[i+2:i+2+nc],dtype=float)
    theta=np.mod(np.arctan2(nodes[:,1],nodes[:,0]),2*np.pi)[quads];z=nodes[quads,2]
    seam=np.ptp(theta,axis=1)>np.pi
    theta[seam]=np.where(theta[seam]<np.pi,theta[seam]+2*np.pi,theta[seam])
    x=theta-theta[:,0:1];y=z-z[:,0:1]
    cross=x*np.roll(y,-1,axis=1)-np.roll(x,-1,axis=1)*y
    area=cross.sum(axis=1)/2
    assert (area>0).all()
    centroid_z=z[:,0]+((y+np.roll(y,-1,axis=1))*cross).sum(axis=1)/(6*area)
    fz=np.sum(p*slope*(.05005+(.05-centroid_z)*slope)*area)
    fx=0.;fy=0.
    for e in range(4):
        de=(e+1)%4;dt=theta[:,de]-theta[:,e];dz=z[:,de]-z[:,e]
        t=theta[:,e,None]+dt[:,None]*u
        zz=z[:,e,None]+dz[:,None]*u
        rh=.05005+(.05-zz)*slope
        fx+=np.sum(p*dz*np.sum(rh*np.sin(t)*w,axis=1))
        fy-=np.sum(p*dz*np.sum(rh*np.cos(t)*w,axis=1))
    force=np.array([fx,fy,fz]);reported=np.array([float(row[k]) for k in ('Fx_N','Fy_N','Fz_N')])
    error=force-reported
    assert np.max(np.abs(error))<.001,(row['label'],error)
    assert fz>=slope*np.hypot(fx,fy)-1e-7
    results.append({'label':row['label'],'force_N':force.tolist(),'difference_N':error.tolist(),
        'vtk_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
return_difference=None
if any(r['kind']=='return-check' for r in rows):
    base=next(r for r in rows if r['level']=='2' and r['kind']=='branch' and np.isclose(float(r['clearance_um']),50,rtol=0,atol=1e-10))
    check=next(r for r in rows if r['kind']=='return-check')
    return_difference={k:float(check[k])-float(base[k]) for k in ('bx_um','by_um','Fz_N','net_feed_L_min')}
    assert abs(return_difference['Fz_N'])<.02
    assert np.hypot(return_difference['bx_um'],return_difference['by_um'])<.001
record={'status':'passed','force_method':'Polygon moments (Fz), eight-point edge Green integrals (Fx/Fy); actual exported Cartesian nodes and pressure',
    'checked_states':len(results),'maximum_force_disagreement_N':max(max(abs(x) for x in r['difference_N']) for r in results),
    'return_check_difference':return_difference,'states':results,
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'states_sha256':hashlib.sha256((out/'states.tsv').read_bytes()).hexdigest()}
target.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k!='states'},indent=2))
