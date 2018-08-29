#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
import matplotlib.image as mpimg
from guacho_utils import *
from matplotlib.colors import LogNorm

jp =1
kp =1
ip =1

mpi_x = 2
mpi_y = 1
mpi_z = 1

plt.ion()
#rmin=100.
#rmax=1e6

#path = '/datos_diable/esquivel/Guacho-1.2/chemH-M2/BIN/'
#path = '../runaway/output/BIN/'
#path = '../runaway/output_snr/BIN/'
#path = '../ran/output/M1/BIN/'
#path = '../TFE/EXO/BIN/'
path = '/datos_diable/esquivel/guacho-working/TFE/'
nout = 0

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

plt.clf()


nh = 50

fig = plt.figure(1)

a   = fig.add_subplot(2,3,1)
plt.imshow(vx[nh,::,::], origin='lower', cmap='Blues' )
a.set_title("Vx")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,2)
plt.imshow(rho[nh,::,::]/np.amax(rho),norm= LogNorm(),origin='lower', cmap='coolwarm' )
a.set_title("Density")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,3)
plt.imshow(bt[nh,::,::],norm= LogNorm(), origin='lower', cmap='Blues' )
a.set_title("B tot")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,4)
plt.imshow(by[nh,::,::],  origin='lower', cmap='Blues' )
a.set_title("By")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,5)
plt.imshow(bz[nh,::,::], origin='lower', cmap='Blues' )
a.set_title("bz")
plt.colorbar(orientation = "horizontal")

a   = fig.add_subplot(2,3,6)
plt.imshow(vt[nh,::,::], origin='lower', cmap='Blues' )
a.set_title("V tot")
plt.colorbar(orientation = "horizontal")
