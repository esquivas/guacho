#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()
rmin=.5
rmax=1.

path = '../pic/BIN/'
nout = 0

rho = get_2d_cut(3,1,nout=nout,neq=0,path=path,verbose=False)
vx  = get_2d_cut(3,1,nout=nout,neq=1,path=path,verbose=False)
vy  = get_2d_cut(3,1,nout=nout,neq=2,path=path,verbose=False)

X,Y = np.meshgrid( np.linspace(-5.,5.,num=64),np.linspace(-5.,5.,num=64) )

plt.figure(2)
plt.clf()
plt.imshow(rho, extent=[-5,5,-5,5], origin='lower', cmap='Spectral', vmin=rmin, vmax=rmax )
plt.colorbar()

Q= plt.quiver(X[::4],Y[::4],vx[::4],vy[::4], pivot='mid')

def readpic(file_in,out):
    f = open(file_in,'rb')
    npart = struct.unpack('1i',f.read(4))[0]
    x  = np.zeros(npart)
    y  = np.zeros(npart)
    z  = np.zeros(npart)
    vx = np.zeros(npart)
    vy = np.zeros(npart)
    vz = np.zeros(npart)
    for i in range(npart):
        x[i],y[i],z[i] = struct.unpack('3d',f.read(24))
        vx[i],vy[i],vz[i]=struct.unpack('3d',f.read(24))
    f.close()
    return npart, x, y, z, vx, vy, vz
part = np.zeros(shape=(50,512,2) )
for (nout) in range(50):
    filepic = path+'pic000.'+str(nout).zfill(3)+'.bin'
    npart, x, y, z, vx, vy, vz = readpic(filepic,nout)
    part[nout,:,0] = x[:]
    part[nout,:,1] = y[:]

for ip in range(npart):
    plt.plot(part[:,ip,0]-5,part[:,ip,1]-5,'o',markersize= 2)
