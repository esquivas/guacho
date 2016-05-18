#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
import matplotlib.image as mpimg
from guacho_utils import *
from matplotlib.colors import LogNorm


plt.ion()
#rmin=100.
#rmax=1e6

#path = '/datos_diable/esquivel/Guacho-1.2/chemH-M2/BIN/'
path = '../runaway/output/BIN/'
nout = 20

rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False)
vx  = readbin3d_all(nout=nout,neq=1,path=path,verbose=False)
vy  = readbin3d_all(nout=nout,neq=2,path=path,verbose=False)
vz  = readbin3d_all(nout=nout,neq=3,path=path,verbose=False)
p   = readbin3d_all(nout=nout,neq=4,path=path,verbose=False)
bx  = readbin3d_all(nout=nout,neq=5,path=path,verbose=False)
by  = readbin3d_all(nout=nout,neq=6,path=path,verbose=False)
bz  = readbin3d_all(nout=nout,neq=7,path=path,verbose=False)

bt  = np.sqrt(bx*bx+by*by+bz*bz) 
vt  = np.sqrt(vx*vx+vy*vy+vz*vz) 


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
bzmax = np.max(bz)
bzmin = np.min(bz)
btmax = np.max(bt)
btmin = np.min(bt)
vtmax = np.max(vt)
vtmin = np.min(vt)

#plt.figure(1)
#plt.clf()
#
#plt.imshow(rho[1,::,::], norm=LogNorm(), vmin=rmin, vmax=rmax, origin='lower', cmap='Blues' )
#plt.colorbar()
#
#plt.figure(2)
#plt.clf()
#
#plt.imshow(by[1,::,::], vmin=bymin, vmax=bymax, origin='lower', cmap='Blues' )
#plt.colorbar()
#
#plt.figure(3)
#plt.clf()
#
#plt.imshow(bt[1,::,::], vmin=btmin, vmax=btmax, origin='lower', cmap='Blues' )
#plt.colorbar()
#
#plt.figure(4)
#plt.clf()
#
#plt.imshow(vx[1,::,::], vmin=vxmin, vmax=vxmax, origin='lower', cmap='Blues' )
#plt.colorbar()


fig = plt.figure(5)
a   = fig.add_subplot(2,3,1)
plt.imshow(vx[1,::,::], vmin=vxmin, vmax=vxmax, origin='lower', cmap='Blues' )
a.set_title("Vx")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,2)
plt.imshow(rho[1,::,::],norm= LogNorm(), vmin=rmin, vmax=rmax, origin='lower', cmap='Blues' )
a.set_title("Density")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,3)
plt.imshow(bt[1,::,::],norm= LogNorm(), vmin=btmin, vmax=btmax, origin='lower', cmap='Blues' )
a.set_title("B tot")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,4)
plt.imshow(by[1,::,::], vmin=bymin, vmax=bymax, origin='lower', cmap='Blues' )
a.set_title("By")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,5)
plt.imshow(bz[1,::,::], vmin=bzmin, vmax=bzmax, origin='lower', cmap='Blues' )
a.set_title("bz")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,6)
plt.imshow(vt[1,::,::], vmin=vtmin, vmax=vtmax, origin='lower', cmap='Blues' )
a.set_title("V tot")
plt.colorbar(orientation = "horizontal")

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

