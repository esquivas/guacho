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

path = '/Users/esquivel/Desktop/datos_diable/guacho-working/HR-CH/BIN/'

firstRead = True
video   = False

MHD = True

for nout in range(38,39):

  if firstRead :
    rhosc = get_scalings(nout=0, path=path, verbose=True )[2]
    nxtot, nytot, nztot = get_boxsize(nout=0, path=path, verbose=True )
    extent_axis = get_extent(nout=0, path=path, verbose=True )
    #scale to solar radius & center
    xax = (extent_axis[0] - 0.5) * 0.05 * 1.496e13 / 6.955e10
    yax = (extent_axis[1] - 0.5) * 0.05 * 1.496e13 / 6.955e10
    zax = (extent_axis[2] - 0.5) * 0.05 * 1.496e13 / 6.955e10
    extent_axis = (xax, yax, zax)

    if (MHD):
      rho   = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=MHD)
      vx    = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=MHD)
      vy    = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=MHD)
      vz    = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=MHD)
      pgas  = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=MHD)
      bx    = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=MHD)
      by    = readbin3d_all(nout=nout,neq=6,path=path,verbose=False,mhd=MHD)
      bz    = readbin3d_all(nout=nout,neq=7,path=path,verbose=False,mhd=MHD)
      rho_n = readbin3d_all(nout=nout,neq=8,path=path,verbose=False,mhd=MHD)*rhosc
      magb=np.sqrt(bx*bx+by*by+bz*bz)
    else:
      rho   = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=MHD)
      vx    = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=MHD)
      vy    = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=MHD)
      vz    = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=MHD)
      pgas  = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=MHD)
      rho_n = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=MHD)*rhosc
 
    magv=np.sqrt(vx*vx+vy*vy+vz*vz)/1e5
    #Temp=pgas*0.5/(8.3145e7*(rho))  # adiabatic EOS
    Temp = pgas/( 8.3145e7*(2.*rho-rho_n) ) #  EOS_H_RATE

    # ionization fraction
    yp= 1.-rho_n/rho

    # density and B field cuts
    rho_c = rho[nztot/2,::,::]
    bx_c  = bx [nztot/2,::,::]
    by_c  = by [nztot/2,::,::]

    rho_c = scipy.ndimage.zoom(rho_c, 2 )
    bx_c  = scipy.ndimage.zoom(bx_c,  2 )
    by_c  = scipy.ndimage.zoom(by_c,  2 )

    if (nout == 0 or firstRead ):
      rL = 1
      texture = np.random.rand(rho_c.shape[0]/rL,rho_c.shape[1]/rL).astype(np.float32)
      texture =scipy.ndimage.zoom(texture,rL)
      XX = np.linspace(extent_axis[0][0],extent_axis[0][1],rho.shape[0])
      YY = np.linspace(extent_axis[1][0],extent_axis[1][1],rho.shape[0])
      ZZ = np.linspace(extent_axis[2][0],extent_axis[2][1],rho.shape[0])
      nskip = 4
      K,Y = np.meshgrid(ZZ[::nskip],XX[::nskip])
      #K = K + nskip/2
      #Y = Y + nskip/2
      #extent = [0.,Lx,0.,Lx]

  if video :
    kernellen=50
    nframes = 3
    frame = nout*nframes

    for t in np.linspace(0,1,nframes):
      kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)*(1+np.sin(2*np.pi*5*(np.arange(kernellen)/float(kernellen)+t)))
      kernel = kernel.astype(np.float32)

      image = vp.line_integral_convolution(bx_c[::,::].astype(np.float32), 
                                           by_c[::,::].astype(np.float32), texture, kernel)
      image = image / image.mean()

      plt.figure(3) ; plt.clf()
      plt.imshow(image*rho[nztot/2,::,::], norm = LogNorm(), origin='lower', 
        vmin=1e-18, vmax=2e-16, cmap = 'viridis' )
      plt.colorbar()
      plt.savefig("PNG/HLLE-ad-%04d.png"%frame,dpi=100, bbox_inches='tight' )

      frame += 1
      print ('generating frame:',frame)
  else:
    
    kernellen=50
    kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)
    kernel = kernel.astype(np.float32)

    image = vp.line_integral_convolution(bx_c[::,::].astype(np.float32), 
                                         by_c[::,::].astype(np.float32), 
                                         texture, kernel)
    image = image / image.mean()

    plt.figure(1) ; plt.clf()
    plt.imshow(Temp[nztot/2,::,::], origin = 'lower',
       cmap = 'gist_heat',
       extent= [ extent_axis[0][0],extent_axis[0][1],extent_axis[1][0],extent_axis[1][1] ] )
    plt.xlabel(r'$X~[\mathrm{R_\odot}]$')
    plt.ylabel(r'$Y~[\mathrm{R_\odot}]$')
    plt.title(r'Temperature')
    plt.colorbar(label='[K]')
    

    plt.figure(2) ; plt.clf()
    plt.imshow(image, origin='lower' , cmap = 'seismic',vmin=0., vmax=2., 
       extent= [ extent_axis[0][0],extent_axis[0][1],extent_axis[1][0],extent_axis[1][1] ] )
    plt.xlabel(r'$X~[\mathrm{R_\odot}]$')
    plt.ylabel(r'$Y~[\mathrm{R_\odot}]$')
    plt.colorbar()
    

    plt.figure(3) ; plt.clf()
    plt.imshow(image*rho_c[::,::], norm = LogNorm(), origin='lower', cmap = 'viridis'
     ,vmin=8e-20, vmax=1e-16,
      extent= [ extent_axis[0][0],extent_axis[0][1],extent_axis[1][0],extent_axis[1][1] ] )
    plt.xlabel(r'$X~[\mathrm{R_\odot}]$')
    plt.ylabel(r'$Y~[\mathrm{R_\odot}]$')
    plt.title(r'density/Bfield')
    plt.colorbar(label=(r'[$\mathrm{g~cm^{-3}}$]') )
    

    plt.figure(4) ; plt.clf()
    plt.imshow(rho[nztot/2,::,::], norm = LogNorm(), origin='lower', cmap = 'Blues'
     ,vmin=8e-20, vmax=2e-16,
        extent= [ extent_axis[0][0],extent_axis[0][1],extent_axis[1][0],extent_axis[1][1] ] )
    plt.xlabel(r'$X~[\mathrm{R_\odot}]$')
    plt.ylabel(r'$Y~[\mathrm{R_\odot}]$')
    plt.title(r'Density')
    plt.colorbar(label=(r'[$\mathrm{g~cm^{-3}}$]') )
    

    plt.figure(5) ; plt.clf()
    plt.imshow(rho_n[nztot/2,::,::],  norm = LogNorm(), origin='lower', cmap = 'Blues',
      extent= [ extent_axis[0][0],extent_axis[0][1],extent_axis[1][0],extent_axis[1][1] ] )
    plt.xlabel(r'$X~[\mathrm{R_\odot}]$')
    plt.ylabel(r'$Y~[\mathrm{R_\odot}]$')
    plt.title(r'Neutral density')
    plt.colorbar(label=(r'[$\mathrm{g~cm^{-3}}$]') )


    plt.figure(6) ; plt.clf()
    plt.imshow(magv[nztot/2,::,::], origin='lower', cmap = 'magma',
      extent= [ extent_axis[0][0],extent_axis[0][1],extent_axis[1][0],extent_axis[1][1] ] )
    plt.xlabel(r'$X~[\mathrm{R_\odot}]$')
    plt.ylabel(r'$Y~[\mathrm{R_\odot}]$')
    plt.title(r'Magnitude of the velocity')
    plt.colorbar(label=(r'[$\mathrm{km~s^{-1}}$]') )


    plt.figure(7) ; plt.clf()
    plt.imshow(yp[nztot/2,::,::], origin='lower', cmap = 'cubehelix', vmin =0.01 , vmax =1.,
      norm = LogNorm(),
      extent= [ extent_axis[0][0],extent_axis[0][1],extent_axis[1][0],extent_axis[1][1] ] )
    plt.xlabel(r'$X~[\mathrm{R_\odot}]$')
    plt.ylabel(r'$Y~[\mathrm{R_\odot}]$')
    plt.title(r'Ionization fraction')
    plt.colorbar(label=(r'') )


    plt.figure(8) ; plt.clf()
    plt.streamplot(K,Y,bx[nxtot/2,::nskip,::nskip],by[nztot/2,::nskip,::nskip],color='k', density=[2,2])
    plt.imshow(magb[nztot/2,::,::], origin='lower', cmap = 'viridis',
      norm = LogNorm(),
      extent= [ extent_axis[0][0],extent_axis[0][1],extent_axis[1][0],extent_axis[1][1] ] )
    plt.xlabel(r'$X~[\mathrm{R_\odot}]$')
    plt.ylabel(r'$Y~[\mathrm{R_\odot}]$')
    plt.title(r'Magnetic Field')
    plt.colorbar(label=(r'') )
    


    print('minmax v:') ; minmax(magv)
    print('minmax T:') ; minmax(Temp)

  #plt.savefig("image.png", dpi = 300 , bbox_inches='tight'  )