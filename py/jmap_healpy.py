# coding: utf-8

get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')


#--------general
import sys, os
from os.path import join

#--------plotting
#magics

import matplotlib
import plt_config
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#matplotlib.use('Qt4Agg')

#--------array
import numpy as np
from scipy.io import FortranFile

#--------yt
import yt
import yt.visualization.volume_rendering.api as vr
from yt.mods import *
from yt.frontends.stream.api import load_uniform_grid
#--------
import aplpy
from aplpy.image_util import percentile_function
#--------astro
from astropy.cosmology import Planck13 as cosmo
from astropy.io import fits
import healpy as hp

#------my routines
import dotunf as unf
from r3n_plotting import *

#---------------------------

#------ dirs --
home=os.getenv('HOME')
scratch='/global/cscratch1/sd/fantaye' #os.getenv('SCRATCH')
srcdir=join(scratch,'seedlets/fits1024/SurfDensMap/')
#os.listdir(srcdir)
print('home:',home)
print('srcdir',srcdir)


#---- original file info
#columns: ##1:snap 2:redshift 3:Nparts 4:shot_noise
id_sim, z_sim, nobj_sim,snoise_sim = np.loadtxt('../data/3Dneedlets/std_total_parts_found_from_049_to_049.dat',unpack=True)
r_sim=cosmo.comoving_distance(z_sim)
rho_file=lambda x:   '../data/3Dneedlets/HealPixMap_nside1024_NoRnd_Radius0_std/SurfDensMap_snap_049.{}.fits'.format(x)
kappa_file=lambda x: '../data/3Dneedlets/LightCones_nside1024_NoRnd_Radius0_std/KappaMap_snap_049.{}.fits'.format(x)
orig_file=rho_file

#---- maps file path
b2b_file=lambda x: join(srcdir, 'beta2b//maps/map.unf_{}'.format(x))
a2b_file=lambda x: join(srcdir, 'a2b//maps//map.unf_{}'.format(x))
beta_file=lambda x, y:join(srcdir, 'a2beta//maps/map.unf_gln1_j{}_r{}.unf'.format(x,y))

#---- power spectra file path
nlist_cln=['b2a/cln.unf','beta2b/cln.unf']
nlist_cl=[ 'b2a/meancl.unf','beta2b/meancl.unf']
clx_file=lambda x: join(srcdir, '{}'.format(x))

#----- figure output dir
rootdir='../'
figdir=rootdir+'figures/'
if not os.path.exists(figdir):
    os.makedirs(figdir)


#=======================


#print os.path.abspath(pview.__file__)

nj=10
nmax=128
lmax=3000
nshell=2*nmax
nside=1024
npix=hp.nside2npix(nside)

#3d projection of healpix maps
(x,y,z),(th,ph),grid_pix = sphere_grid(nside)

#two ways to get mollwide image
vec2pix_func = lambda x,y,z: hp.pixelfunc.vec2pix(nside,x,y,z,nest=False)
mw2image=hp.projector.MollweideProj() # 1) for any map mw2image(map,vec2pix_func) 
mw2img=MollweideImg()  #2) for any map - mw2img.projmap(map)


#---------------------

#------------ plot cln ---
def load_clx():
    cln = unf_plot(flambda=clx_file, nlist=nlist_cln,names=['b2a','beta2b'])
    cls = unf_plot(flambda=clx_file, nlist=nlist_cl,names=['b2a','beta2b'],shape=(lmax+1,2))

    return cln, cls
#-------------


#============= Needlet window =============
def gln_data(j=4,view=True):
    fname=join(srcdir,'a2beta/gln.unf')
    cnt=(lmax+1)*nshell*nj
    gln=unf.runf(fname,count=cnt, shape=(lmax+1,nshell,nj))

    x = np.linspace(0, lmax, lmax+1)
    y = np.linspace(0, nmax, nmax+1)
    x,y=np.meshgrid(y,x)
    z = gln[:,0:nmax+1,j]

    xyz=(x,y,z)

    if view:
        gln_view(xyz)

    return xyz

def gln_view(xyz):
    fig = plt.figure(figsize=(13,10))
    ax_gln = fig.gca(projection='3d')
    ax_gln.plot_surface(x, y, z,  rstride=4, cstride=4,cmap='jet',
                        linewidth=0, antialiased=False)




#============= Orignal ball =============
def band_limit(map, nlmax=None,nsmax=None,iter=3,fwhm=20):
    if nlmax is None: nlmax=lmax
    if nsmax is None: nsmax=nside

    alm=hp.map2alm(map, lmax=nlmax,iter=iter,use_weights=True,pol=False)
    return hp.alm2map(alm, nsmax,pol=False,pixwin=False,fwhm=np.deg2rad(fwhm/60.0))

def orig_vs_recon(nshell=3,irfirst=2,view=True):

    #morig=np.ndarray((1000,1000))
    #mb2b=np.ndarray((1000,1000))
    morig=0.0
    mb2b=0.0
    for ir in range(irfirst,nshell):
        mr=band_limit(hp.read_map(orig_file(ir)))
        
        #mr=unf.runf(a2b_file(ir),verbose=1)
        #r dip
        mr= hp.remove_dipole(mr,fitval=False,verbose=True)        

        mr_b2b=unf.runf(a2b_file(ir),verbose=1)
        #mr_b2b=unf.runf(b2b_file(ir))
        #r dip
        mr_b2b= hp.remove_dipole(mr_b2b,fitval=False,verbose=True)        

        morig = morig + mr[grid_pix] #- mr.mean()
        mb2b =  mb2b  + mr_b2b[grid_pix]

        if ir%10 ==0: print('ir={}'.format(ir))
    #- end ir loop
    if view: 
        view_balls(morig, mb2b,nshell=nshell)

    return morig, mb2b, (mr, mr_b2b)

def view_balls(morig, mb2b,nshell=3):

    fig,(ax1,ax2,ax3) = plt.subplots(ncols=3,figsize=(18,6))
    h1=ax1.pcolormesh(x, y, morig,cmap='RdBu_r') #, alpha=0.3
    ax1.set_title('Orignal')
    ax1.set_aspect("equal")
    show_colorbar(h1, ax=ax1)

    h2=ax2.pcolormesh(x, y, mb2b,cmap='RdBu_r') #, alpha=0.3
    ax2.set_title('Recovered')
    ax2.set_aspect("equal")
    show_colorbar(h2, ax=ax2)

    h3=ax3.pcolormesh(x, y, morig-mb2b,cmap='RdBu_r') #, alpha=0.3
    ax3.set_title('Difference')
    ax3.set_aspect("equal")
    show_colorbar(h3, ax=ax3)

    plt.tight_layout()
    plt.savefig(figdir+'sum_nshells{}_orig_ball.png'.format(nshell))

#===========================================

#============== Needlet shells ===================
def jmaps_data(jvec=[1,2,3],nshell=3, nside=1024, 
               view=True):

    npix=hp.nside2npix(nside)
    grid_map=[]
    for j in jvec:
        mm=np.ndarray((1000,1000))
        for ir in range(nshell):
            fname=beta_file(ir, j)
            mr=unf.runf(fname,count=npix)
            mm = mm+mr[grid_pix]
        #- end ir loop

        print('j={}'.format(j))
        grid_map.append(mm)

    if view:
        jmaps_view(grid_map)
    
    return grid_map

#-----------------

def jmaps_view(grid_map):

    ncols=min([len(grid_map),3])

    fig,axes = plt.subplots(ncols=ncols,figsize=(18,6))
    ax=axes.reshape(-1)

    for i,axj in enumerate(ax):
        h1=axj.pcolormesh(grid_map[i],cmap='RdBu_r') #, alpha=0.3
        axj.set_title('j1')
        axj.set_aspect("equal")
        show_colorbar(h1, ax=axj)

    plt.tight_layout()
    plt.savefig(figdir+'sum_nj{}_jmap.png'.format(ncols))

    #mapj.append(m)
            
#mapj[0].shape

#---------------------


#yt vis
#ds, cam = cube_vis(mapj[0])


#get spherical grid. x,y,z draws a sphere with unit radius
#th and ph are the projection of the sphere onto mollwide plane



#fig = plt.figure(figsize=(15,12))
#ax = fig.gca(projection='3d')
#
# h=ax.plot_surface(
#     x, y, z,  rstride=1, cstride=1,cmap='RdBu_r', alpha=0.3, linewidth=0)
#
#ax.set_zlim([-1,1])


#ax.set_xlim([-1,1])
#ax.set_ylim([-1,1])
#ax_ball.plot_surface(th, ph, grid_map,  rstride=4, cstride=4,cmap='jet',
#    linewidth=0, antialiased=False)

#===========================================


#---------------------
# In[34]:

z='''
The following code is taken from 
http://zonca.github.io/2013/03/interactive-3d-plot-of-sky-map.html


fig, (ax1, ax2)=plt.subplots(ncols=2,figsize=(15,6))
# load a WMAP map
npix=hp.nside2npix(512)
m = np.random.randn(npix)
#m = hp.read_map("data/wmap_band_iqumap_r9_7yr_W_v4.fits", 0) * 1e3 # muK
hp.mollview(m,ax=ax1)

nside = hp.npix2nside(len(m))

#value limit
vmin = -1e3; vmax = 1e3
# size of the grid
xsize = ysize = 1000

theta = np.linspace(np.pi, 0, ysize)
phi   = np.linspace(-np.pi, np.pi, xsize)
longitude = np.radians(np.linspace(-180, 180, xsize))
latitude = np.radians(np.linspace(-90, 90, ysize))

# project the map to a rectangular matrix xsize x ysize
PHI, THETA = np.meshgrid(phi, theta)
grid_pix = hp.ang2pix(nside, THETA, PHI)
grid_map = m[grid_pix]

# Create a sphere
r = 0.3
x = r*np.sin(THETA)*np.cos(PHI)
y = r*np.sin(THETA)*np.sin(PHI)
z = r*np.cos(THETA)

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
mlab.clf()

mlab.mesh(x, y, z, scalars=grid_map, colormap="jet", vmin=vmin, vmax=vmax) 
'''


# In[ ]:



