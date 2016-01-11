from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


lmax=200
nmax=200
j0=0
bb=2d0
nj=10


def unf(inputfilename):
    shape = (lmax+1,nmax+1,nj)
    fd = open(fname, 'rb')
    data = np.fromfile(file=fd, dtype=np.double).reshape(shape)
    fd.close()
    return data





fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(0, lmax, lmax+1)
y = np.linspace(0, lmax, lmax+1)


gln=unf('gln_data.unf')

j=6
fname = '../figures/gln_3d_lmax'+str(lmax)+'_nmax'+str(nmax)+'_j'+str(j)+'.pdf'
z = gln[:,:,6]

ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


# Actually savethe figure
plt.savefig(gln_)

# Close it
plt.close()


