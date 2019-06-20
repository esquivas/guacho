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
nout = 50

rho = get_2d_cut(3,1,nout=nout,neq=0,path=path,verbose=False)
vx  = get_2d_cut(3,1,nout=nout,neq=1,path=path,verbose=False)
vy  = get_2d_cut(3,1,nout=nout,neq=2,path=path,verbose=False)

X,Y = np.meshgrid( np.linspace(-5.,5.,num=64),np.linspace(-5.,5.,num=64) )

plt.figure(1)
plt.clf()
plt.imshow(rho, extent=[-5,5,-5,5], origin='lower', cmap='Spectral', vmin=rmin, vmax=rmax )
plt.colorbar()

Q= plt.quiver(X[::4],Y[::4],vx[::4],vy[::4], pivot='mid')
