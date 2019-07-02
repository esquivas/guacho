#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()
rmin=.5
rmax=1.

path = '../picBlast/BIN/'
nout = 10

nproc = get_Nproc(nout,path=path)

rho  = get_2d_cut(3,1,nout=nout,neq=0,path=path,verbose=False)
vx   = get_2d_cut(3,1,nout=nout,neq=1,path=path,verbose=False)
vy   = get_2d_cut(3,1,nout=nout,neq=2,path=path,verbose=False)

divV = get_2d_cut(3,1,nout=nout,neq=0,path=path,base='divV-',verbose=False)


X,Y = np.meshgrid( np.linspace(-5.,5.,num=rho.shape[0]),np.linspace(-5.,5.,num=rho.shape[1]) )

plt.figure(1)
plt.clf()
plt.imshow(divV, extent=[0,1,0,1], origin='lower', cmap='Spectral' )
plt.colorbar()

#Q= plt.quiver(X[::4],Y[::4],vx[::4],vy[::4], pivot='mid',scale=20)


part = np.zeros(shape=(10,4096,2) )
for (nout) in range(nout,nout+1):
    id,xp, yp, zp, vxp, vyp, vzp = readpic(nout,path=path)
    part[nout,:,0] = xp[:]
    part[nout,:,1] = yp[:]
    plt.plot(xp,yp,'o',markersize= 1)

#for ip in range(1024):
#    plt.plot(part[:,ip,0],part[:,ip,1],'-o',markersize= 2)
#    plt.xlim([0,1])
#    plt.ylim([0,1])
