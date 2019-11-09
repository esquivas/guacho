#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()

path = '/datos_diable/esquivel/Guacho-1.2/chemH-M2/BIN/'
nout = 0

# this reads the 3D data for the first equation (1 in fortran 0 here)
rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False)
plt.figure(1)
plt.clf()
#  plots the 2D cut at midplane of a 256^3 sim (y=128)
plt.imshow(rho[::,128,::], norm=LogNorm(), origin='lower', cmap='Blues' )
plt.colorbar()
print '****'

#  alternatively we can extract the 2D cut with this, without storing the
#  entire 3D data
cut= get_2d_cut(2,128,nout=nout,neq=0,path=path,verbose=False)
plt.figure(2)
plt.clf()
plt.imshow(cut, norm=LogNorm(),  origin='lower', cmap='Blues' )
plt.colorbar()

#  for a component of the velocity and thermal pressure one can use
vz = readbin3d_all(nout=nout,neq=3,path=path,verbose=False)
Pgas = readbin3d_all(nout=nout,neq=4,path=path,verbose=False)
