#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()

nytot=100
nxtot=100

mu=0.6
Rg=8.314e7

dtprint=1.

path1 = '../prueba/grav/split_all/mhd/2temp/BIN/'
path2 = '../prueba/grav/non_split/mhd/2temp/BIN/'
#path11 = '../prueba/grav/split_all/mhd/10G/BIN/'
#path2 = '../prueba/grav/non_split/mhd/depulso/10G/BIN/'
path3 = '../prueba/grav/split/hd/2temp/BIN/'
#path4 = '../prueba/grav/non_split/hd/pulso/BIN/'

for nout in range(0,1):
  
  rho_sp = readbin3d_all(nout=nout,neq=0,path=path1,verbose=False, mhd=True)
  vx_sp  = readbin3d_all(nout=nout,neq=1,path=path1,verbose=False, mhd=True)
  vy_sp  = readbin3d_all(nout=nout,neq=2,path=path1,verbose=False, mhd=True)
  vz_sp  = readbin3d_all(nout=nout,neq=3,path=path1,verbose=False, mhd=True)
  P_sp   = readbin3d_all(nout=nout,neq=4,path=path1,verbose=False, mhd=True)
  bx_sp  = readbin3d_all(nout=nout,neq=5,path=path1,verbose=False, mhd=True)
  by_sp  = readbin3d_all(nout=nout,neq=6,path=path1,verbose=False, mhd=True)
  bz_sp  = readbin3d_all(nout=nout,neq=7,path=path1,verbose=False, mhd=True)
  
  #P_i   = readbin3d_all(nout=0,neq=4,path=path11,verbose=False, mhd=True)
  rho_spf = readbin3d_all(nout=100,neq=0,path=path1,verbose=False, mhd=True)
  P_spf   = readbin3d_all(nout=100,neq=4,path=path1,verbose=False, mhd=True)

  rho = readbin3d_all(nout=nout,neq=0,path=path2,verbose=False, mhd=True)
  vx  = readbin3d_all(nout=nout,neq=1,path=path2,verbose=False, mhd=True)
  vy  = readbin3d_all(nout=nout,neq=2,path=path2,verbose=False, mhd=True)
  vz  = readbin3d_all(nout=nout,neq=3,path=path2,verbose=False, mhd=True)
  P   = readbin3d_all(nout=nout,neq=4,path=path2,verbose=False, mhd=True)
  bx  = readbin3d_all(nout=nout,neq=5,path=path2,verbose=False, mhd=True)
  by  = readbin3d_all(nout=nout,neq=6,path=path2,verbose=False, mhd=True)
  bz  = readbin3d_all(nout=nout,neq=7,path=path2,verbose=False, mhd=True)
  
  #rho_f = readbin3d_all(nout=100,neq=0,path=path2,verbose=False, mhd=True)
  #P_f   = readbin3d_all(nout=100,neq=4,path=path2,verbose=False, mhd=True)
  
  rho_hdsp = readbin3d_all(nout=nout,neq=0,path=path3,verbose=False, mhd=False)
  vx_hdsp  = readbin3d_all(nout=nout,neq=1,path=path3,verbose=False, mhd=False)
  vy_hdsp  = readbin3d_all(nout=nout,neq=2,path=path3,verbose=False, mhd=False)
  vz_hdsp  = readbin3d_all(nout=nout,neq=3,path=path3,verbose=False, mhd=False)
  P_hdsp   = readbin3d_all(nout=nout,neq=4,path=path3,verbose=False, mhd=False)
  
  rho_f = readbin3d_all(nout=100,neq=0,path=path3,verbose=False, mhd=False)
  P_f   = readbin3d_all(nout=100,neq=4,path=path3,verbose=False, mhd=False)
  
  #rho_hd = readbin3d_all(nout=nout,neq=0,path=path4,verbose=False, mhd=False)
  #vx_hd  = readbin3d_all(nout=nout,neq=1,path=path4,verbose=False, mhd=False)
  #vy_hd  = readbin3d_all(nout=nout,neq=2,path=path4,verbose=False, mhd=False)
  #vz_hd  = readbin3d_all(nout=nout,neq=3,path=path4,verbose=False, mhd=False)
  #P_hd   = readbin3d_all(nout=nout,neq=4,path=path4,verbose=False, mhd=False)
 
  #V=np.sqrt(vx*vx+vy*vy+vz*vz)
  #V_sp=np.sqrt(vx_sp*vx_sp+vy_sp*vy_sp+vz_sp*vz_sp)
  B_sp=np.sqrt(bx_sp*bx_sp+by_sp*by_sp+bz_sp*bz_sp)  
  B=np.sqrt(bx*bx+by*by+bz*bz)  


  Temp=mu*P/rho/Rg
  Temp_sp=mu*P_sp/rho_sp/Rg
  Temp_hdsp=mu*P_hdsp/rho_hdsp/Rg
  Temp_spf=mu*P_spf/rho_spf/Rg
  Temp_f=mu*P_f/rho_f/Rg
  #Temp_f=mu*P_f/rho_f/Rg
  
  tiempo=dtprint*nout
  
  if (nout==0):
    Pmin=P.min()
    Pmax=P.max()
    rhomin=rho.min()
    rhomax=rho.max()

  #plt.figure(1)
  #plt.clf()
  #plt.imshow(rho[1,::,::],norm=LogNorm(),origin='lower', cmap='Blues' )
  #plt.colorbar()
  ##plt.plot(rho[1,::,50],label='Densidad final')
  ##plt.legend()
  #plt.title('Densidad')
  #print '****'
  #
  #plt.figure(2)
  #plt.clf()
  #plt.imshow(Temp[1,::,::], origin='lower', cmap='Blues' )
  #plt.colorbar()
  #plt.plot(Temp[1,::,50],label='Temperatura inicial')
  #plt.legend()
  #
  #plt.figure(1)
  #plt.clf()
  #plt.imshow(P_sp[1,::,::], vmin=Pmin, vmax=Pmax, origin='lower', norm=LogNorm(),cmap='Blues')
  #plt.colorbar()
  ##cbar=plt.colorbar()
  ##cbar.set_ticks([P_sp.min(), P_sp.max()])
  ##cbar.set_ticklabels([P_sp.min(), P_sp.max()])
  ##plt.legend()
  #plt.title('Presion MHD split all tiempo='+str(tiempo)+'seg')
  
  #plt.figure(2)
  #plt.clf()
  #plt.imshow(P[1,::,::], vmin=P.min(), vmax=P.max(), origin='lower', cmap='Blues')
  #plt.colorbar()
  ##cbar=plt.colorbar()
  ##cbar.set_ticks([P.min(), P.max()])
  ##cbar.set_ticklabels([P.min(), P.max()])
  ##plt.plot(P[1,::,50],label='Presion final')
  ##plt.legend()
  #plt.title('Presion MHD non_split tiempo='+str(tiempo)+'seg')
  
  #plt.figure(3)
  #plt.clf()
  #plt.imshow(P_hdsp[1,::,::], vmin=Pmin, vmax=Pmax, origin='lower', cmap='Blues')
  #plt.colorbar()
  ##cbar=plt.colorbar()
  ##cbar.set_ticks([P_hd.min(), P_hd.max()])
  ##cbar.set_ticklabels([P_hd.min(), P_hd.max()])
  ##plt.plot(P[1,::,50],label='Presion final')
  ##plt.legend()
  #plt.title('Presion HD split')
  
  #plt.figure(4)
  #plt.clf()
  #plt.imshow(P_hd[1,::,::], vmin=Pmin, vmax=Pmax, origin='lower', cmap='Blues')
  #plt.colorbar()
  ##cbar=plt.colorbar()
  ##cbar.set_ticks([P_hd.min(), P_hd.max()])
  ##cbar.set_ticklabels([P_hd.min(), P_hd.max()])
  ##plt.plot(P[1,::,50],label='Presion final')
  ##plt.legend()
  #plt.title('Presion HD non_split')
  
  plt.figure(1)
  plt.clf()
  plt.plot(P_sp[1,::,50]/Pmax,'o-',label='mhd split inicial')
  #plt.plot(P[1,::,50]/Pmax,label='sin split')
  plt.plot(P_hdsp[1,::,50]/Pmax,'.-',label='hd split inicial')
  plt.plot(P_spf[1,::,50]/Pmax,label='mhd split final')
  plt.plot(P_f[1,::,50]/Pmax,'--',label='hd split final')
  plt.yscale('log')
  plt.title('Presion')
  plt.legend()
  
  plt.figure(2)
  plt.clf()
  plt.plot(rho_sp[1,::,50]/rhomax,'o-',label='mhd split inicial')
  #plt.plot(rho[1,::,50]/rhomax,label='sin split')
  plt.plot(rho_hdsp[1,::,50]/rhomax,'.-',label='hd split inicial')
  plt.plot(rho_spf[1,::,50]/rhomax,label='mhd split final')
  plt.plot(rho_f[1,::,50]/rhomax,'--',label='hd split final')
  plt.yscale('log')
  #plt.plot(P_i[1,::,50],label='split sin pulso')
  #plt.plot(P_f[1,::,50],label='split sin pulso')
  #plt.plot(P_hdsp[1,::,50],label='hd split')
  #plt.plot(P[1,::,50],label='sin split')
  plt.title('Densidad')
  plt.legend()
  
  #plt.figure(2)
  #plt.clf()
  #plt.imshow(P_hdsp[1,::,::], norm=LogNorm(), origin='lower', cmap='Blues')
  #cbar=plt.colorbar()
  #cbar.set_ticks([P_hdsp.min()-0.1, P_hdsp.max()+0.1])
  #cbar.set_ticklabels([P_hdsp.min()-0.1, P_hdsp.max()+0.1])
  ##plt.legend()
  #plt.title('Presion HD split all')
  
  plt.figure(3)
  plt.clf()
  plt.plot(Temp_sp[1,::,50],'o-',label='mhd split inicial')
  #plt.plot(Temp[1,::,50],label='sin split')
  plt.plot(Temp_hdsp[1,::,50],'.-',label='hd split inicial')
  plt.plot(Temp_spf[1,::,50],label='mhd split final')
  plt.plot(Temp_f[1,::,50],'--',label='hd split final')
  plt.legend(loc='lower right')
  plt.title('Temperatura')
  
  #plt.figure(4)
  #plt.clf()
  #plt.imshow(Temp[1,::,::], vmin=Temp.min(),vmax=Temp.max(), norm=LogNorm(), origin='lower' )
  #plt.colorbar(ticks=(Temp.min(), Temp.max()))
  ##plt.plot(P[1,::,50],label='Presion final')
  ##plt.legend()
  #plt.title('Temp')

  X,Y = np.meshgrid( np.arange(0.,nxtot,4.),np.arange(0.,nytot,4.))
  X,Z = np.meshgrid( np.arange(0.,nxtot,4.),np.arange(0.,nytot,4.))

  ##plt.figure(3)
  ##plt.clf()
  ##U=vx[1,::4,::4]
  ##W=vy[1,::4,::4]

  ##plt.imshow(V[1,:,:], origin='lower')
  ##plt.colorbar()
  ##plt.quiver(X,Y,U,W, pivot='mid')
  ##plt.title('Velocidad')
  
  #plt.figure(3)
  #plt.clf()
  #U=bx_sp[1,::4,::4]
  #W=by_sp[1,::4,::4]

  #plt.imshow(B_sp[1,:,:], origin='lower')
  #plt.colorbar()
  #plt.quiver(X,Y,U,W, pivot='mid')
  #plt.title('Campo split')

  #plt.figure(4)
  #plt.clf()
  #U=bx[1,::4,::4]
  #W=by[1,::4,::4]

  #plt.imshow(B[1,:,:], origin='lower')
  #plt.colorbar()
  #plt.quiver(X,Y,U,W, pivot='mid')
  #plt.title('Campo')

  plt.draw()