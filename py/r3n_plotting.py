#general
import sys, os
from os.path import join

#plotting
import matplotlib
import plt_config
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#matplotlib.use('Qt4Agg')

#array
import numpy as np
import pandas as pd
import scipy as sp

#yt
import yt
import yt.visualization.volume_rendering.api as vr
from yt.mods import *
from yt.frontends.stream.api import load_uniform_grid
#
import aplpy
from aplpy.image_util import percentile_function
#astro
from astropy.cosmology import Planck13 as cosmo
from astropy.io import fits
import healpy as hp

#my routines
import dotunf as unf

#---------------------
figdir='./figures/'
if not os.path.exists(figdir):
    os.makedirs(figdir)

#------------------------

class unf_plot(object):

    def __init__(self,fname=None,flambda=None,shape=None,
                 nlist=None,index=None, names=None,verbose=1):


        dd={}

        if fname:            
            dd['data'] = unf.runf(fname, shape=shape)
            
        elif flambda and nlist:

            if names is None:
                names=nlist

            dd={}
            for x in nlist:
                fname=flambda(x)
                if verbose>0:
                    print('x, fname: ',x, fname)
                dd[x]=unf.runf(fname, shape=shape)

        #here make the DF
        self.df=pd.DataFrame(dd)

        if index is None:
            try:
                self.pd.index=index
            except:
                print('passed index can not be used')


    #--------------------

    def plot(self,ix=None,x=None,iy=None, ax=None, func=None,
             *kargs, **kwargs):
        
        if ix is None:
            if x is None:
                x=range(len(self.data[:,0]))
        else:
            x=self.data[:,ix]


        if iy is None:
            y=self.data
        else:
            y=self.data[:,iy]

        #define what tool to use to plot
        if func is None:
            if ax is None:
                func=plt.plot
        else:
            func=ax.plot

        #now plot
        h=func(x,y,*kargs, **kargs)

        return h

#============================

class MollweideImg(hp.projector.MollweideProj):
    def projmap(self,map,nest=False,**kwds):
        nside = hp.pixelfunc.npix2nside(hp.pixelfunc.get_map_size(map))
        f = lambda x,y,z: hp.pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(MollweideImg,self).projmap(map,f,**kwds)


def sphere_grid(nside,xsize=None,ysize=None):

    if xsize is None: xsize=1000
    if ysize is None: ysize=1000

    theta = np.linspace(np.pi, 0, ysize)
    phi   = np.linspace(-np.pi, np.pi, xsize)
    longitude = np.radians(np.linspace(-180, 180, xsize))
    latitude = np.radians(np.linspace(-90, 90, ysize))

    # project the map to a rectangular matrix xsize x ysize
    PHI, THETA = np.meshgrid(phi, theta)
    grid_pix = hp.ang2pix(nside, THETA, PHI)

    # Create a sphere
    r = 0.5
    x = r*np.sin(THETA)*np.cos(PHI)
    y = r*np.sin(THETA)*np.sin(PHI)
    z = r*np.cos(THETA)

    return (x,y,z),(PHI,THETA),grid_pix

def volume_opaque(ds,figdir=figdir,fkey=''):
    # We start by building a transfer function, and initializing a camera.

    #if figdir is None: figdir=figdir

    tf = yt.ColorTransferFunction((-30, -22))
    cam = ds.camera([0.5, 0.5, 0.5], [0.2, 0.3, 0.4], 0.10, 256, tf)

    # Now let's add some isocontours, and take a snapshot.
    
    tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5], colormap = 'RdBu_r')
    cam.snapshot(figdir+fkey+"v1.png", clip_ratio=6.0)

    # In this case, the default alphas used (np.logspace(-3,0,Nbins)) does not
    # accentuate the outer regions of the galaxy. Let's start by bringing up the
    # alpha values for each contour to go between 0.1 and 1.0

    tf.clear()
    tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
                  alpha=np.logspace(0,0,4), colormap = 'RdBu_r')
    cam.snapshot(figdir+fkey+"v2.png", clip_ratio=6.0)
    
    # Now let's set the grey_opacity to True.  This should make the inner portions
    # start to be obcured
    
    tf.grey_opacity = True
    cam.snapshot(figdir+fkey+"v3.png", clip_ratio=6.0)
    
    # That looks pretty good, but let's start bumping up the opacity.
    
    tf.clear()
    tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
        alpha=10.0*np.ones(4,dtype='float64'), colormap = 'RdBu_r')
    cam.snapshot(figdir+fkey+"v4.png", clip_ratio=6.0)
    
    # Let's bump up again to see if we can obscure the inner contour.    
    tf.clear()
    tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
                  alpha=30.0*np.ones(4,dtype='float64'), colormap = 'RdBu_r')
    cam.snapshot(figdir+fkey+"v5.png", clip_ratio=6.0)

    # Now we are losing sight of everything.  Let's see if we can obscure the next
    # layer
    
    tf.clear()
    tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
        alpha=100.0*np.ones(4,dtype='float64'), colormap = 'RdBu_r')
    cam.snapshot(figdir+fkey+"v6.png", clip_ratio=6.0)
    
    # That is very opaque!  Now lets go back and see what it would look like with
    # grey_opacity = False
    
    tf.grey_opacity=False
    cam.snapshot(figdir+fkey+"v7.png", clip_ratio=6.0)
    
    # That looks pretty different, but the main thing is that you can see that the
    # inner contours are somewhat visible again.  


def volume_plot(ds,figdir=figdir, fkey=''):
#from http://yt-project.org/doc/cookbook/complex_plots.html#opaque-volume-rendering

    # Create a data container (like a sphere or region) that
    # represents the entire domain.
    ad = ds.all_data()

    # Get the minimum and maximum densities.
    mi, ma = ad.quantities.extrema("density")


    # Create a transfer function to map field values to colors.
    # We bump up our minimum to cut out some of the background fluid
    tf = yt.ColorTransferFunction((np.log10(mi)+2.0, np.log10(ma)))
    
    # Add three Gaussians, evenly spaced between the min and
    # max specified above with widths of 0.02 and using the
    # gist_stern colormap.
    tf.add_layers(3, w=0.02, colormap="gist_stern")
    
    # Choose a center for the render.
    c = [0.5, 0.5, 0.5]
    
    # Choose a vector representing the viewing direction.
    L = [0.5, 0.2, 0.7]
    
    # Set the width of the image.
    # Decreasing or increasing this value
    # results in a zoom in or out.
    W = 1.0
    
    # The number of pixels along one side of the image.
    # The final image will have Npixel^2 pixels.
    Npixels = 512
    
    # Create a camera object.
    # This object creates the images and
    # can be moved and rotated.
    cam = ds.camera(c, L, W, Npixels, tf)
    
    # Create a snapshot.
    # The return value of this function could also be accepted, modified (or saved
    # for later manipulation) and then put written out using write_bitmap.
    # clip_ratio applies a maximum to the function, which is set to that value
    # times the .std() of the array.
    im = cam.snapshot(figdir+fkey+"%s_volume_rendered.png" % ds, clip_ratio=8.0)
    
    # Add the domain edges, with an alpha blending of 0.3:
    nim = cam.draw_domain(im, alpha=0.3)
    nim.write_png(figdir+fkey+'%s_vr_domain.png' % ds)
    
    # Add the grids, colored by the grid level with the algae colormap
    nim = cam.draw_grids(im, alpha=0.3, cmap='algae')
    nim.write_png(figdir+fkey+'%s_vr_grids.png' % ds)
    
    # Here we can draw the coordinate vectors on top of the image by processing
    # it through the camera. Then save it out.
    cam.draw_coordinate_vectors(nim)
    nim.write_png(figdir+fkey+"%s_vr_vectors.png" % ds)

    return cam 

def slice_vis(ds,L = [0.5, 0.2, 0.9], fname=None):
    # The vector normal to the slicing plane
    #L = [0.5, 0.2, 0.9]

    # The command to generate the off axis projection plot
    p = OffAxisProjectionPlot(ds, L, 'density')
    
    # Applying a grey colormap to all data
    p.set_cmap(field='all', cmap='Greys')

    # Saving the image
    if not fname is None: 
        p.save(fname)

    return p

def vol_vis(ds, figdir=figdir, fkey=''):
    ad = ds.all_data()

    # Custom min and max to show using function in APLpy.
    # Haven't used APLpy yet? Check it out at aplpy.github.com.
    #auto_v = percentile_function(ad['density'])
    #mi, ma = auto_v(0.15), auto_v(99.95)
    
    # Get the minimum and maximum densities.
    mi, ma = ad.quantities.extrema("density")

    # Setting up parameters for volume rendering. See the following
    # link for more details on the parameters:
    # http://yt-project.org/doc/cookbook/simple_plots.html#cookbook-simple-volume-rendering
    tf = ColorTransferFunction((mi, ma))
    tf.add_layers(4, w=0.3)
    c = [0.5, 0.5, 0.5]
    L = [0.5, 0.2, 0.9]
    W = 1.0
    Nvec = 512
    cam = ds.h.camera(c, L, W, Nvec, tf)
    image = cam.snapshot(figdir+fkey+"%s_volume_rendered.png" % ds, 8.0)

    return cam

def cube_vis(cube):
    # Loading data into yt structure
    R_unit="mpc" #r_sim.unit.to_string()
    R=r_sim.value.max()

    print('Box Radius and unit:',R,R_unit)
    data = {}
    data["density"] = (cube,R_unit)

    bbox = np.array([[-0.5,0.5],[-0.5,0.5],[0,1]]) # bbox of width 1 on a side with center (0,0,0)
    ds = yt.load_uniform_grid(data, cube.shape, length_unit=(2*R,R_unit), 
                              nprocs=1, bbox=bbox,geometry="spherical")

    #cam = volume_plot(ds)
    cam = vol_vis(ds)

    return ds,cam

def show_colorbar(h,ax=plt.gca(),**kargs):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    #Source:
    #
    #http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    #
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    if h is None: h=ax.get_data()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(h, cax=cax,**kargs)
#==============================================    
