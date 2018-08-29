#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()

path = '../TwoFluid/BIN/'
nout = 10

rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False)
plt.figure(1) ;plt.clf()
plt.imshow(rho[::,25,::], norm=LogNorm(), origin='lower', cmap='viridis', vmin=1e-20, vmax=1e-18 )
plt.colorbar()
plt.title("Density of Ions")

rhoN = readbin3d_all(nout=nout,neq=0,path=path,base='pointsN',verbose=False)
plt.figure(2) ;plt.clf()
plt.imshow(rhoN[::,25,::], norm=LogNorm(), origin='lower', cmap='viridis', vmin=1e-20, vmax=1e-18 )
plt.colorbar()
plt.title("Density of Neutrals")

Pi = readbin3d_all(nout=nout,neq=4,path=path,verbose=False)
plt.figure(3) ;plt.clf()
plt.imshow(Pi[::,25,::], norm=LogNorm(), origin='lower', cmap='Spectral', vmin=1e-7, vmax=5e-6 )
plt.colorbar()
plt.title("Pressure of Ions")


Pn = readbin3d_all(nout=nout,neq=4,path=path,base='pointsN',verbose=False)
plt.figure(4) ;plt.clf()
plt.imshow(Pn[::,25,::], norm=LogNorm(), origin='lower', cmap='Spectral', vmin=1e-7)
plt.colorbar()
plt.title("Pressure of Neutrals")
