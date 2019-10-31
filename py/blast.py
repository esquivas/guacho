#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter

plt.ion()

path = '../snr/BIN/'
nout = 10

nproc = get_Nproc(nout,path=path)
boxsize= get_boxsize(nout,path=path)

cut = int(boxsize[2]/2)

rho    = get_2d_cut(3,cut,nout=nout,neq=0,path=path,verbose=False)
vx     = get_2d_cut(3,cut,nout=nout,neq=1,path=path,verbose=False)
vy     = get_2d_cut(3,cut,nout=nout,neq=2,path=path,verbose=False)
bx     = get_2d_cut(3,cut,nout=nout,neq=5,path=path,verbose=False, mhd=True)
by     = get_2d_cut(3,cut,nout=nout,neq=6,path=path,verbose=False, mhd=True)
shockF = get_2d_cut(3,cut,nout=nout,neq=0,path=path,base='shock-',verbose=False)

f = open(path + 'stokes-'+str(nout).zfill(3)+'.bin','rb')
stokes = np.fromfile(f, dtype = 'd', count=(boxsize[0]*boxsize[1]*3)).reshape(boxsize[0],boxsize[1],3,order='F').T
f.close()

extent =  [-12,12,-12,12]

plt.figure(2)
plt.clf()
plt.imshow(shockF, extent=extent, origin='lower', cmap='rainbow',interpolation = 'none' )
#plt.colorbar()
plt.title('Shock detector (red if shocked)')

X,Y = np.meshgrid( np.linspace(-12.,12.,num=rho.shape[0]),np.linspace(-12.,12.,num=rho.shape[1]) )

fig = plt.figure(1) ; plt.clf()
ax =  fig.add_subplot(111); ax.set_aspect(aspect=1.0)
plt.figure(3) ; plt.clf()

npart = 1169

GeV = 1.e9*1.60218E-12 # 1 GeV in ergs

for (noutp) in range(nout+1):
    lmpData, SED, P1, P2= read_lmp(noutp,path=path)
    lmpData[:,1:3] = (lmpData[:,1:3]-0.5) * 24.0
    SED[npart,:,0] = SED[npart,:,0]# / GeV  # convert from erg to GeV
    #lmpData= read_lmp(noutp,path=path)
    print(lmpData.shape)
    fig = plt.figure(1)
    ax.plot(lmpData[:,1],lmpData[:,2],'o',markersize= 1)
    ax.plot(lmpData[npart,1],lmpData[npart,2],'*',markersize=15, color='green')
    plt.figure(3)
    plt.loglog(SED[npart,:,0] / GeV,SED[npart,:,1],'-o',label=str(noutp),markersize=2)
   # plt.loglog(SED[npart,:,0],SED[npart,:,1]*SED[npart,:,0]**2,label=str(noutp))
plt.xlabel(r'E [Gev]')
plt.ylabel(r'N(E)')

lmpData, SED, P1, P2= read_lmp(nout,path=path)
lmpData[:,1:3] = (lmpData[:,1:3]-0.5) * 24.0

plt.legend()
plt.figure(2)
sc = plt.scatter(lmpData[:,1],lmpData[:,2], c=lmpData[:,5],alpha=0.3,cmap='inferno')
cb = plt.colorbar(sc)
cb.set_label(r'$\theta_{B1}$')
Q= plt.streamplot( X,Y,bx,by, density = 2., color='silver',linewidth=1., minlength=0.8)
plt.xlim(-12.,12.) ; plt.xlabel(r'x (pc)')
plt.ylim(-12.,12.) ; plt.ylabel(r'y (pc)')


plt.figure(4)
plt.clf()
map = plt.imshow(rho, extent=extent, origin='lower', cmap='binary_r',interpolation = 'none', norm=LogNorm() )
sc  = plt.scatter(lmpData[:,1],lmpData[:,2], c=lmpData[:,4],s =2.,alpha=0.5,cmap='inferno',vmin=1.,vmax=4.)
plt.xlim(-12.,12.) ; plt.xlabel(r'x (pc)')
plt.ylim(-12.,12.) ; plt.ylabel(r'y (pc)')
cb1 = plt.colorbar(map)
cb2 = plt.colorbar(sc)
cb1.set_label(r'$\rho$ [g cm$^{-3}$]')
cb2.set_label(r'Compression factor')

stokes_Is = gaussian_filter(stokes[0,:,:], sigma=2)
plt.figure(5); plt.clf()
plt.imshow(stokes_Is[:,:], origin='lower', extent=extent )
plt.colorbar()
plt.xlabel(r'x (pc)')
plt.ylabel(r'y (pc)')
plt.title(r'Stokes I')

stokes_Qs = gaussian_filter(stokes[1,:,:], sigma=2)
plt.figure(6); plt.clf()
plt.imshow(stokes_Is[:,:], origin='lower', extent=extent )
plt.colorbar()
plt.xlabel(r'x (pc)')
plt.ylabel(r'y (pc)')
plt.title(r'Stokes Q')

stokes_Us = gaussian_filter(stokes[2,:,:], sigma=2)
plt.figure(7); plt.clf()
plt.imshow(stokes_Us[:,:], origin='lower', extent=extent )
plt.colorbar()
plt.xlabel(r'x (pc)')
plt.ylabel(r'y (pc)')
plt.title(r'Stokes U')

PI = np.sqrt(stokes[1,:,:]**2+stokes[2,:,:]**2)/(stokes[0,:,:]+1e-20)
PIs = gaussian_filter(PI[:,:], sigma=2)
plt.figure(8); plt.clf()
plt.imshow(PIs[:,:], origin='lower',  extent=extent )
plt.colorbar()
plt.xlabel(r'x (pc)')
plt.ylabel(r'y (pc)')
plt.title(r'Polarization Degree')

#  plot Polarization vectors
X, Y = np.meshgrid(np.arange(-12., 12.,24./boxsize[0]), np.arange(-12.,12.,24./boxsize[1]) )
px = np.sin(np.arctan2(stokes_Us[:,:],stokes_Qs[:,:])/2. + np.pi/2.)
py = np.cos(np.arctan2(stokes_Us[:,:],stokes_Qs[:,:])/2. + np.pi/2.)
plt.figure(9); plt.clf()
map = plt.imshow(rho, extent=extent, origin='lower', cmap='binary_r',interpolation = 'none', norm=LogNorm() )
cb1 = plt.colorbar(map)
Xp  = X[::4,::4]
Yp  = Y[::4,::4]
pxp = px[::4,::4]
pyp = py[::4,::4]
Isp = stokes_Is[::4,::4]
mask = (Isp > 1e-20)
vec = plt.quiver(Xp[mask], Yp[mask], pxp[mask], pyp[mask], Isp[mask],
pivot='mid',scale=20., headlength=0, headwidth=1., width=0.006, cmap='inferno_r')
plt.title(r'Polarization vectors')
plt.xlabel(r'x (pc)')
plt.ylabel(r'y (pc)')

#plt.figure(1) ; plt.savefig('fig1.png', transparent=True, bbox_inches='tight',dpi=300)
#plt.figure(2) ; plt.savefig('fig2.png', transparent=True, bbox_inches='tight',dpi=300)
#plt.figure(3) ; plt.savefig('fig3.png', transparent=True, bbox_inches='tight',dpi=300)
#plt.figure(4) ; plt.savefig('fig4.png', transparent=True, bbox_inches='tight',dpi=300)
