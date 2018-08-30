#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()

path = '../TwoFluid/BIN/'
nout = 6

rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False)
plt.figure(1) ;plt.clf()
plt.imshow(rho[::,rho.shape[1]/2,::], norm=LogNorm(), origin='lower', cmap='jet', vmin=1e-20, vmax=1e-18 )
plt.colorbar()
plt.title("Density of Ions")

rhoN = readbin3d_all(nout=nout,neq=0,path=path,base='pointsN',verbose=False)
plt.figure(2) ;plt.clf()
plt.imshow(rhoN[::,rho.shape[1]/2,::], norm=LogNorm(), origin='lower', cmap='jet', vmin=1e-20, vmax=1e-18 )
plt.colorbar()
plt.title("Density of Neutrals")

Pi = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=True)
plt.figure(3) ;plt.clf()
plt.imshow(Pi[::,rho.shape[1]/2,::], norm=LogNorm(), origin='lower', cmap='Spectral', vmin=1e-7, vmax=5e-6 )
plt.colorbar()
plt.title("Pressure of Ions")


Pn = readbin3d_all(nout=nout,neq=4,path=path,base='pointsN',verbose=False, mhd=True)
plt.figure(4) ;plt.clf()
plt.imshow(Pn[::,rho.shape[1]/2,::], norm=LogNorm(), origin='lower', cmap='Spectral', vmin=1e-7)
plt.colorbar()
plt.title("Pressure of Neutrals")

vxn = readbin3d_all(nout=nout,neq=1,path=path,base='pointsN',verbose=False)/1e5
plt.figure(5) ;plt.clf()
plt.imshow(vxn[::,rho.shape[1]/2,::], origin='lower', cmap='seismic')
plt.colorbar()
plt.title("Neutral Velocity in x")

vx = readbin3d_all(nout=nout,neq=1,path=path,base='points',verbose=False)/1e5
plt.figure(6) ;plt.clf()
plt.imshow(vx[::,rho.shape[1]/2,::], origin='lower', cmap='seismic')
plt.colorbar()
plt.title("Ion Velocity in x")
