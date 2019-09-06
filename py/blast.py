#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm

plt.ion()

path = '../picBlast/BIN/'
nout = 10

nproc = get_Nproc(nout,path=path)

rho  = get_2d_cut(3,1,nout=nout,neq=0,path=path,verbose=False)
vx   = get_2d_cut(3,1,nout=nout,neq=1,path=path,verbose=False)
vy   = get_2d_cut(3,1,nout=nout,neq=2,path=path,verbose=False)
bx   = get_2d_cut(3,1,nout=nout,neq=5,path=path,verbose=False, mhd=True)
by   = get_2d_cut(3,1,nout=nout,neq=6,path=path,verbose=False, mhd=True)
shockF = get_2d_cut(3,1,nout=nout,neq=0,path=path,base='shock-',verbose=False)

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

plt.figure(3) ; plt.clf()

npart = 475#250


for (noutp) in range(nout+1):
    picData, SED, P1, P2= readpic(noutp,path=path)
    #picData= readpic(noutp,path=path)
    print(picData.shape)
    plt.figure(1)
    plt.plot(picData[:,1],picData[:,2],'o',markersize= 1)
    plt.plot(picData[npart,1],picData[npart,2],'*',markersize=15, color='green')
    plt.figure(3)
    plt.loglog(SED[npart,:,0],SED[npart,:,1]*SED[npart,:,0]**2,label=str(noutp))

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
plt.imshow(shockF, extent=extent, origin='lower', cmap='viridis',interpolation = 'none' )
sc = plt.scatter(picData[:,1],picData[:,2], c=(picData[:,5]),alpha=0.5,cmap='viridis')
cb = plt.colorbar(sc)
