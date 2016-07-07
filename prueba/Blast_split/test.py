#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm  

plt.ion()

path1 = './HLLE-Athena/BIN/'
path2 = './HLLE-Split-Athena/BIN/'

dtprint=0.01

for nout in range(1,100,1):
  
  rho_sp = readbin3d_all(nout=nout,neq=0,path=path1,verbose=False, mhd=True)
  vx_sp  = readbin3d_all(nout=nout,neq=1,path=path1,verbose=False, mhd=True)
  vy_sp  = readbin3d_all(nout=nout,neq=2,path=path1,verbose=False, mhd=True)
  vz_sp  = readbin3d_all(nout=nout,neq=3,path=path1,verbose=False, mhd=True)
  P_sp   = readbin3d_all(nout=nout,neq=4,path=path1,verbose=False, mhd=True)
  bx_sp  = readbin3d_all(nout=nout,neq=5,path=path1,verbose=False, mhd=True)
  by_sp  = readbin3d_all(nout=nout,neq=6,path=path1,verbose=False, mhd=True)
  bz_sp  = readbin3d_all(nout=nout,neq=7,path=path1,verbose=False, mhd=True)
  
  rho = readbin3d_all(nout=nout,neq=0,path=path2,verbose=False, mhd=True)
  vx  = readbin3d_all(nout=nout,neq=1,path=path2,verbose=False, mhd=True)
  vy  = readbin3d_all(nout=nout,neq=2,path=path2,verbose=False, mhd=True)
  vz  = readbin3d_all(nout=nout,neq=3,path=path2,verbose=False, mhd=True)
  P   = readbin3d_all(nout=nout,neq=4,path=path2,verbose=False, mhd=True)
  bx  = readbin3d_all(nout=nout,neq=5,path=path2,verbose=False, mhd=True)
  by  = readbin3d_all(nout=nout,neq=6,path=path2,verbose=False, mhd=True)
  bz  = readbin3d_all(nout=nout,neq=7,path=path2,verbose=False, mhd=True)

  tiempo=dtprint*nout

  plt.figure(1)
  plt.clf()
  plt.imshow(rho_sp[1,::,::], norm=LogNorm(),origin='lower',vmin=0.08, vmax=6.5)
  plt.colorbar()
  #plt.plot(rho[1,::,50],label='Densidad final')
  #plt.legend()
  plt.title('Densidad MHD split all tiempo='+str(tiempo)+'seg')
  
  plt.figure(2)
  plt.clf()
  plt.imshow(rho[1,::,::], norm=LogNorm(),origin='lower',vmin=0.08, vmax=6.5)
  plt.colorbar()
  #plt.plot(rho[1,::,50],label='Densidad final')
  #plt.legend()
  plt.title('Densidad MHD tiempo='+str(tiempo)+'seg')
  
  plt.draw()
