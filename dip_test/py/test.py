#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm
import scipy.ndimage
import vectorplot as vp

def minmax(a):
  print ( a.min(), a.max() )

plt.ion()

path = '/Users/esquivel/Desktop/datos_diable/guacho-working/NC-HLLE/BIN/'
nout = 0
firstRead = True
video   = True

Hydro = False

for nout in range(0,101):

  if firstRead :
    rhosc = get_scalings(nout=0, path=path, verbose=True )[2]
    nxtot, nytot, nztot = get_boxsize(nout=0, path=path, verbose=True )
  
    if (Hydro):
      rho   = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=False)
      vx    = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=False)
      vy    = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=False)
      vz    = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=False)
      pgas  = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=False)
      rho_n = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=False)*rhosc
    else:
      rho   = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=True)
      vx    = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=True)
      vy    = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=True)
      vz    = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=True)
      pgas  = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=True)
      bx    = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=True)
      by    = readbin3d_all(nout=nout,neq=6,path=path,verbose=False,mhd=True)
      bz    = readbin3d_all(nout=nout,neq=7,path=path,verbose=False,mhd=True)
      rho_n = readbin3d_all(nout=nout,neq=8,path=path,verbose=False,mhd=True)*rhosc
      magb=np.sqrt(bx*bx+by*by+bz*bz)    #magb=np.sqrt(bx*bx+by*by+bz*bz)

    magv=np.sqrt(vx*vx+vy*vy+vz*vz)
    Temp=pgas*0.5/(8.3145e7*(rho))  # adiabatic EOS

  rL = 1
  texture = np.random.rand(nxtot/rL,nytot/rL).astype(np.float32)
  texture =scipy.ndimage.zoom(texture,rL, order=2)

  if video :
    kernellen=50
    nframes = 3
    frame = nout*nframes

    for t in np.linspace(0,1,nframes):
      kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)*(1+np.sin(2*np.pi*5*(np.arange(kernellen)/float(kernellen)+t)))
      kernel = kernel.astype(np.float32)

      image = vp.line_integral_convolution(bx[nztot/2,::,::].astype(np.float32), 
                                           by[nztot/2,::,::].astype(np.float32), texture, kernel)
      image = image / image.mean()

      plt.figure(3) ; plt.clf()
      plt.imshow(image*rho[nztot/2,::,::], norm = LogNorm(), origin='lower', 
        vmin=1e-18, vmax=2e-16, cmap = 'viridis' )
      plt.colorbar()
      plt.savefig("PNG/image-%04d.png"%frame,dpi=100, bbox_inches='tight' )

      frame += 1
      print ('generating frame:',frame)
  else:
    '''
    kernellen=50
    kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)
    kernel = kernel.astype(np.float32)

    image = vp.line_integral_convolution(bx[nztot/2,::,::].astype(np.float32), 
                                         by[nztot/2,::,::].astype(np.float32), 
                                         texture, kernel)
    image = image / image.mean()

    plt.figure(1) ; plt.clf()
    plt.imshow(Temp[nztot/2,::,::], norm = LogNorm(), origin = 'lower'
      , vmin=Temp.min(), vmax=Temp.max(), cmap = 'inferno', interpolation = 'none' )
    plt.title(r'Temperature [K]')
    plt.colorbar()

    plt.figure(2) ; plt.clf()
    plt.imshow(image, origin='lower' , cmap = 'seismic',vmin=0., vmax=2.)
    plt.colorbar()
    
    plt.figure(3) ; plt.clf()
    plt.imshow(image*rho[nztot/2,::,::], norm = LogNorm(), origin='lower', cmap = 'viridis'
     ,vmin=1e-20, vmax=1e-16 )
    plt.colorbar()
    '''

    plt.figure(4) ; plt.clf()
    plt.imshow(rho[nztot/2,::,::], norm = LogNorm(), origin='lower', cmap = 'viridis'
     ,vmin=1e-18, vmax=2e-16 )
    plt.colorbar()
    plt.title(r'Density')

    plt.figure(5) ; plt.clf()
    plt.imshow(rho_n[nztot/2,::,::],  norm = LogNorm(), origin='lower', cmap = 'viridis'
      , vmin = 1e-24, vmax=1e-20 )
    plt.colorbar()
    plt.title(r'Passive Scalar')

  #plt.savefig("image.png", dpi = 300 , bbox_inches='tight'  )