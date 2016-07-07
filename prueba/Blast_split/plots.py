#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm
#import yt

def minmax(u) :
  print u.min(), u.max()

plt.ion()
nout = 100 ; figN = 1

#models = ('HLLE-CD','HLLE-8W','HLLE-NC','HLLD-CD','HLLD-8W','HLLD-NC',)
models = ('HLLE-Split-Athena','HLLE-Athena')

for model in models:

  path = './'+model+'/BIN/'
  rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False, mhd=True)
  Pg  = readbin3d_all(nout=nout,neq=4,path=path,verbose=False, mhd=True)
  vx  = readbin3d_all(nout=nout,neq=1,path=path,verbose=False, mhd=True)
  vy  = readbin3d_all(nout=nout,neq=2,path=path,verbose=False, mhd=True)
  vz  = readbin3d_all(nout=nout,neq=3,path=path,verbose=False, mhd=True)
  bx  = readbin3d_all(nout=nout,neq=5,path=path,verbose=False, mhd=True)
  by  = readbin3d_all(nout=nout,neq=6,path=path,verbose=False, mhd=True)
  bz  = readbin3d_all(nout=nout,neq=7,path=path,verbose=False, mhd=True)
  V    = np.sqrt(vx*vx+vy*vy+vz*vz)
  B2   = bx*bx+by*by+bz*bz
  Pm   = B2/8./np.pi
  Temp = Pg/rho
  divB = readbin3d_all(nout=nout,neq=0,path=path,base='divB-',verbose=False, mhd=True)
  print 'divB:' ; minmax(divB)

  plt.figure(figN) ; plt.clf() ; figN += 1
  plt.imshow(rho[1,::,::], origin='lower', vmin=0.08, vmax=6.5 )
  plt.colorbar() ; plt.title(r'$\rho$ ('+model+')')
  plt.savefig('rho-'+model+'-'+str(nout).zfill(3)+'.png', bbox_inches='tight')

  #plt.figure(figN) ; plt.clf() ; figN += 1
  ##plt.imshow(Temp[1,::,::], origin='lower', vmin=0.15, vmax=1.24)
  #plt.imshow(Pg[1,::,::], origin='lower')
  #plt.colorbar() ; plt.title('Pg  ('+model+')')
  #plt.savefig('Pg-'+model+'-'+str(nout).zfill(3)+'.png', bbox_inches='tight')
  
  #plt.figure(figN) ; plt.clf() ; figN += 1
  ##plt.imshow(Temp[1,::,::], origin='lower', vmin=0.15, vmax=1.24)
  #plt.imshow(V[1,::,::], origin='lower')
  #plt.colorbar() ; plt.title('Velocidad  ('+model+')')
  #plt.savefig('V-'+model+'-'+str(nout).zfill(3)+'.png', bbox_inches='tight')
  
  #plt.figure(figN) ; plt.clf() ; figN += 1
  ##plt.imshow(Temp[1,::,::], origin='lower', vmin=0.15, vmax=1.24)
  #plt.imshow(Pm[1,::,::], origin='lower')
  #plt.colorbar() ; plt.title('Presion magnetica  ('+model+')')
  #plt.savefig('Pm-'+model+'-'+str(nout).zfill(3)+'.png', bbox_inches='tight')

  plt.figure(figN) ; plt.clf() ; figN += 1
  plt.imshow(divB[1,::,::], origin='lower')
  plt.colorbar() ; plt.title(r' $\nabla \cdot B$ '+'('+model+')')
  plt.savefig('divB-'+model+'-'+str(nout).zfill(3)+'.png', bbox_inches='tight')


