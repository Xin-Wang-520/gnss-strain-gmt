"""Verify Xu source commands, tensor convention and analytic derivatives."""
from pathlib import Path
import sys, tempfile, subprocess, json, os, shutil
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from run4_calculate_strain_eu import ROOT,RATE_KEYS,xu_strain,read_grid_binary
from GNSS_analysis.gnss_processing import locate_gmt,gmt_environment
import numpy as np


def main():
    output=ROOT/'GNSS_analysis/output/run4_strain_eu/bm30m_c25_i20m_pixel'
    with np.load(output/'strain_eu.npz') as saved:
        fields={k:saved[k] for k in saved.files}
    manifest=json.loads((output/'run4_summary.json').read_text())
    assert not manifest['spherical_curvature_terms']
    gmt,env=locate_gmt(None),gmt_environment()
    tensor=np.empty(fields['ee'].shape+(2,2))
    tensor[...,0,0]=fields['ee'];tensor[...,1,1]=fields['nn']
    tensor[...,0,1]=tensor[...,1,0]=fields['en']
    eig=np.linalg.eigvalsh(tensor)
    np.testing.assert_allclose(fields['max_shear'],(eig[...,1]-eig[...,0])/2,atol=5e-5,rtol=1e-6)
    np.testing.assert_allclose(fields['dilatation'],eig.sum(axis=-1),atol=5e-5,rtol=1e-5)
    for key in RATE_KEYS:
        lon,lat,values=read_grid_binary(gmt,output/f'{key}_eu.nc',env)
        np.testing.assert_array_equal(values,fields[key])
        np.testing.assert_array_equal(lon,fields['longitude'])
        np.testing.assert_array_equal(lat,fields['latitude'])
    if not shutil.which('csh'):
        raise RuntimeError('csh required to verify original Xu scripts')
    source=Path(manifest['input_dir'])
    with tempfile.TemporaryDirectory(prefix='check_xu_original_') as tmp:
        base=Path(tmp)
        original_env=env.copy();original_env['PATH']=str(Path(gmt).parent)+os.pathsep+original_env['PATH']
        for script,products in [('plot_tensor.csh',{'exx':'ee','eyy':'nn','e12':'en','e_2nd':'second_invariant'}),
                                ('plot_shear.csh',{'shear':'max_shear'}),
                                ('plot_dilatation.csh',{'dilatation':'dilatation'})]:
            folder=base/Path(script).stem;folder.mkdir()
            def call(*args):
                return subprocess.run([gmt,*map(str,args)],cwd=folder,env=original_env,check=True,capture_output=True)
            call('set','PROJ_ELLIPSOID','6378000','PROJ_MEAN_RADIUS','authalic')
            for component in ('u','v'):
                call('grdmath',source/f'velocity_eu_{component}.nc','1000','DIV','=',f'{component}.grd')
            call('grdmath','u.grd','0','MUL','1','ADD','=','dem.grd')
            text=(ROOT/'xu'/script).read_text()
            start=text.index('# compute strain\n');stop=text.index('gmt grdgradient dem_tmp.grd',start)
            driver=folder/'original_commands.csh'
            driver.write_text('set uu = u.grd\nset vv = v.grd\n'+text[start:stop])
            subprocess.run([shutil.which('csh'),'-ef',str(driver)],cwd=folder,env=original_env,check=True,capture_output=True)
            for product,key in products.items():
                _,_,actual=read_grid_binary(gmt,folder/f'{product}.grd',env)
                np.testing.assert_array_equal(actual,fields[key])
    print('PASS: six saved fields match independently executed original Xu calculation blocks; tensor convention and NC/NPZ agree.')
    # Construct a known smooth field with nonzero normal and cross gradients.
    # Its planar derivatives have no spherical curvature terms.
    radius=6378000.;step=.01
    lon=90+np.arange(-4,5)*step;lat=32+np.arange(-4,5)*step
    lam,phi=np.meshgrid(np.deg2rad(lon),np.deg2rad(lat))
    x=radius*np.cos(np.deg2rad(32))*(lam-np.deg2rad(90))
    y=radius*(phi-np.deg2rad(32))
    a,b,c,d=80e-9,40e-9,20e-9,-50e-9
    ve=(a*x+b*y)*1000;vn=(c*x+d*y)*1000
    gridlon,gridlat=np.meshgrid(lon,lat)
    cache=ROOT/'GNSS_analysis/.cache';cache.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='check_run4_xu_',dir=cache) as tmp:
        folder=Path(tmp)
        for name,array in [('u.nc',ve),('v.nc',vn)]:
            xyz=np.column_stack((gridlon.ravel(),gridlat.ravel(),array.ravel())).astype(np.float64)
            subprocess.run([gmt,'xyz2grd','-bi3d','-R89.955/90.045/31.955/32.045',
                '-I0.01','-r','-fg','-G'+str(folder/name)],input=xyz.tobytes(),env=env,
                check=True,capture_output=True)
        actual_lon,actual_lat,result,_,_=xu_strain(gmt,folder/'u.nc',folder/'v.nc',folder/'work',6378,env)
        scale=np.cos(np.deg2rad(32))/np.cos(np.deg2rad(actual_lat))[:,None]
        expected_ee=np.broadcast_to(a*scale*1e9,result['ee'].shape)
        expected_nn=np.full_like(expected_ee,d*1e9)
        expected_en=np.broadcast_to((b+c*scale)/2*1e9,expected_ee.shape)
        core=np.s_[1:-1,1:-1]
        for key,expected in [('ee',expected_ee),('nn',expected_nn),('en',expected_en)]:
            np.testing.assert_allclose(result[key][core],expected[core],atol=0.001,rtol=0)
    print('PASS: analytic planar field verifies east/north signs, mm-to-m units and tensor EN factor 1/2.')


if __name__=='__main__':main()
