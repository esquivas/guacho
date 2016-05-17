#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()
#rmin=100.
#rmax=1e6

#path = '/datos_diable/esquivel/Guacho-1.2/chemH-M2/BIN/'
path = '../runaway/HLLD-CD/BIN/'
nout = 2

rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False)
vx  = readbin3d_all(nout=nout,neq=1,path=path,verbose=False)
vy  = readbin3d_all(nout=nout,neq=2,path=path,verbose=False)
vz  = readbin3d_all(nout=nout,neq=3,path=path,verbose=False)
p   = readbin3d_all(nout=nout,neq=4,path=path,verbose=False)
bx  = readbin3d_all(nout=nout,neq=5,path=path,verbose=False)
by  = readbin3d_all(nout=nout,neq=6,path=path,verbose=False)
bz  = readbin3d_all(nout=nout,neq=7,path=path,verbose=False)

#calculates mins and max
rmin  = np.min(rho)
rmax  = np.max(rho)
vxmax = np.max(vx)
vxmin = np.min(vx)
vymax = np.max(vy)
vymin = np.min(vy)
vzmax = np.max(vz)
vzmin = np.min(vz)
bxmax = np.max(bx)
bxmin = np.min(bx)
bymax = np.max(by)
bymin = np.min(by)
bymax = np.max(bz)
bymin = np.min(bz)


plt.figure(1)
plt.clf()

plt.imshow(rho[1,::,::], norm=LogNorm(), vmin=rmin, vmax=rmax, origin='lower', cmap='Blues' )
plt.colorbar()

plt.figure(2)
plt.clf()

plt.imshow(by[1,::,::], vmin=bymin, vmax=bymax, origin='lower', cmap='Blues' )
plt.colorbar()



#print '****'
#
#cut= get_2d_cut(2,96,nout=nout,neq=0,path=path,verbose=False)
#plt.figure(2)
#plt.clf()
#plt.imshow(cut, norm=LogNorm(), vmin=rmin, vmax=rmax, origin='lower', cmap='Blues' )
#plt.colorbar()
#
#vz = readbin3d_all(nout=nout,neq=3,path=path,verbose=False)
#Pgas = readbin3d_all(nout=nout,neq=4,path=path,verbose=False)

#x_ax, y_ax, z_ax = get_axis(nout=0, path=path,verbose=False)

