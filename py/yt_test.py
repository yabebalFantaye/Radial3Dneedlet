import numpy as np
from yt.mods import *
from yt.frontends.stream.api import load_uniform_grid
from astropy.io import fits
from aplpy.image_util import percentile_function
 
# Loading fits data
cube = fits.getdata('PerA_12coFCRAO_F_vxy.fits')
 
 # Masking out nan elements
cube[np.isnan(d)] = np.nanmin(d)

# Loading data into yt structure
data = dict(Density = cube)
pf = load_uniform_grid(data, cube.shape, 9e16, nprocs=16)

# Custom min and max to show using function in APLpy.
# Haven't used APLpy yet? Check it out at aplpy.github.com.
auto_v = percentile_function(data['Density'])
mi, ma = auto_v(0.15), auto_v(99.95)

# Setting up parameters for volume rendering. See the following
# link for more details on the parameters:
# http://yt-project.org/doc/cookbook/simple_plots.html#cookbook-simple-volume-rendering
tf = ColorTransferFunction((mi, ma))
tf.add_layers(4, w=0.3)
c = [0.5, 0.5, 0.5]
L = [0.5, 0.2, 0.9]
W = 1.0
Nvec = 512
cam = pf.h.camera(c, L, W, Nvec, tf)
image = cam.snapshot("%s_volume_rendered.png" % pf, 8.0)
