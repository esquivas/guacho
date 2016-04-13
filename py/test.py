#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()
rmin=100.
rmax=1e6

path = '/datos_diable/esquivel/Guacho-1.2/chemH-M2/BIN/'
nout = 0

rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False)
plt.figure(1)
plt.clf()
plt.imshow(rho[::,96,::], norm=LogNorm(), vmin=rmin, vmax=rmax, origin='lower', cmap='Blues' )
plt.colorbar()
print '****'

cut= get_2d_cut(2,96,nout=nout,neq=0,path=path,verbose=False)
plt.figure(2)
plt.clf()
plt.imshow(cut, norm=LogNorm(), vmin=rmin, vmax=rmax, origin='lower', cmap='Blues' )
plt.colorbar()

vz = readbin3d_all(nout=nout,neq=3,path=path,verbose=False)
Pgas = readbin3d_all(nout=nout,neq=4,path=path,verbose=False)

#x_ax, y_ax, z_ax = get_axis(nout=0, path=path,verbose=False)

