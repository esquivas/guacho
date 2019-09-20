#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()

path = '../snr/BIN/'
nout = 10

nproc = get_Nproc(nout,path=path)

rho  = get_2d_cut(3,64,nout=nout,neq=0,path=path,verbose=False)
vx   = get_2d_cut(3,64,nout=nout,neq=1,path=path,verbose=False)
vy   = get_2d_cut(3,64,nout=nout,neq=2,path=path,verbose=False)
bx   = get_2d_cut(3,64,nout=nout,neq=5,path=path,verbose=False, mhd=True)
by   = get_2d_cut(3,64,nout=nout,neq=6,path=path,verbose=False, mhd=True)
shockF = get_2d_cut(3,64,nout=nout,neq=0,path=path,base='shock-',verbose=False)

extent =  [0,1,0,1]

plt.figure(2)
plt.clf()
plt.imshow(shockF, extent=extent, origin='lower', cmap='rainbow',interpolation = 'none' )
#plt.colorbar()
plt.title('Shock detector (1 if shocked)')

X,Y = np.meshgrid( np.linspace(0.,1.,num=rho.shape[0]),np.linspace(0.,1.,num=rho.shape[1]) )

'''
plt.figure(1)
plt.clf()
plt.imshow(rho, extent=extent, origin='lower', cmap='Spectral')
plt.colorbar()
plt.title('Density')
'''

plt.figure(1) ; plt.clf()
plt.figure(3) ; plt.clf()

npart = 266

GeV = 1.e9*1.60218E-12 # 1 GeV in ergs

for (noutp) in range(nout+1):
    picData, SED, P1, P2= readpic(noutp,path=path)
    SED[npart,:,0] = SED[npart,:,0]# / GeV  # convert from erg to GeV
    #picData= readpic(noutp,path=path)
    print(picData.shape)
    plt.figure(1)
    plt.plot(picData[:,1],picData[:,2],'o',markersize= 1)
    plt.plot(picData[npart,1],picData[npart,2],'*',markersize=15, color='green')
    plt.figure(3)
    plt.loglog(SED[npart,:,0] / GeV,SED[npart,:,1],'-o',label=str(noutp),markersize=2)
   # plt.loglog(SED[npart,:,0],SED[npart,:,1]*SED[npart,:,0]**2,label=str(noutp))
    plt.xlabel(r'E [Gev]')
plt.ylabel(r'N(E)')
picData, SED, P1, P2= readpic(nout,path=path)

plt.legend()
plt.figure(2)
sc = plt.scatter(picData[:,1],picData[:,2], c=picData[:,4],alpha=1.,cmap='inferno',vmin=1.)
cb = plt.colorbar(sc)
Q= plt.streamplot( X,Y,bx,by, density = 2., color='silver',linewidth=1., minlength=0.8)

#picDataS = readpic(nout,path=path, base='pic-shocked-')
#fig = plt.figure(2)
#plt.plot(picDataS[:,1], picDataS[:,2], '*', markersize=5, color='blue',alpha=0.5)
#Q= plt.quiver(X[::8],Y[::8],bx[::8],by[::8], pivot='mid', scale=40, color='w', alpha=0.5)

plt.figure(4)
plt.clf()
map = plt.imshow(rho, extent=extent, origin='lower', cmap='inferno',interpolation = 'none', norm=LogNorm() )
sc  = plt.scatter(picData[:,1],picData[:,2], c=(picData[:,5]),alpha=0.5,cmap='viridis')
plt.xlim(0.,1.)
plt.ylim(0.,1.)
cb = plt.colorbar(sc)

#  get scalings, working on units problem
rsc,vsc,rhosc,Bsc = get_scalings(nout,path=path,verbose=False)
Psc = rhosc*vsc**2
Esc = rsc**3*Psc
