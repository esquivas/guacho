#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()

path = '/datos_diable/esquivel/Guacho-1.2/chemH-M2/BIN/'
path = '/home/oden/programs/guachoDev/Jorge/BIN/'
nout = 50

# this reads the 3D data for the first equation (1 in fortran 0 here)
phi = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,base="ph-")
rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False)
nx = rho.shape[0]
ny = rho.shape[1]
# for H_Alpha
f   = open(path+'h_alpha-'+str(nout).zfill(3)+'.bin','rb')
halpha = np.fromfile(f, dtype = 'd', count=(nx*ny)).reshape(nx,ny,order='F').T
f.close()

# for xrays
f   = open(path+'xray_'+str(nout).zfill(3)+'.bin','rb')
xrays = np.fromfile(f, dtype = 'd', count=(2*nx*ny)).reshape(nx,ny,2,order='F').T
f.close()


'''
plt.figure(1)
plt.clf()
plt.imshow(np.sum(rho,axis=0), norm=LogNorm(), origin='lower', cmap='plasma')
plt.colorbar()
#  plots the 2D cut at midplane of a 256^3 sim (y=128)

plt.figure(2)
plt.clf()
plt.imshow(np.sum(phi,axis=0), norm=LogNorm(), origin='lower', cmap='plasma')
plt.colorbar()
'''

plt.figure(3)
plt.clf()
plt.imshow(halpha, norm=LogNorm(), origin='lower', cmap='gist_heat',vmin=1e-20)
plt.colorbar()

plt.figure(4)
plt.clf()
plt.imshow(xrays[1,:,:], norm=LogNorm(), origin='lower', cmap='plasma',vmin=1e-3)
plt.colorbar()


