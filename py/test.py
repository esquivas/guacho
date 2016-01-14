#!/usr/bin/python
import numpy as np
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()
rmin=100.
rmax=1e6

path = "/datos_europa/esquivel/Guacho-Test/BIN/"
nout = 15

rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=True)
plt.figure(1)
plt.clf()
plt.imshow(rho[::,50,::], norm=LogNorm(), vmin=rmin, vmax=rmax, origin='lower', cmap='Blues' )
plt.colorbar()
print '****'

cut= get_2d_cut(2,50,nout=nout,neq=0,path=path,verbose=False)
plt.figure(2)
plt.clf()
plt.imshow(cut, norm=LogNorm(), vmin=rmin, vmax=rmax, origin='lower', cmap='Blues' )
plt.colorbar()


#x_ax, y_ax, z_ax = get_axis(nout=0, path=path,verbose=False)

