#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()
rmin=100.
rmax=1e6

path = '../pic/BIN/'
nout = 0

rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False)
plt.figure(1)
plt.clf()
plt.imshow(rho[1,:,:], norm=LogNorm(), origin='lower', cmap='Blues' )
plt.colorbar()

