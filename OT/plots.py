#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm
import yt

def minmax(u) :
  print u.min(), u.max()

plt.ion()
nout = 5 ; figN = 1

models = ('HLLE-CD','HLLE-8W','HLLE-NC','HLLD-CD','HLLD-8W','HLLD-NC',)

for model in models:

  path = './'+model+'/BIN/'
  rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False, mhd=True)
  Pg  = readbin3d_all(nout=nout,neq=4,path=path,verbose=False, mhd=True)
  Temp = Pg/rho
  divB = readbin3d_all(nout=nout,neq=0,path=path,base='divB-',verbose=False, mhd=True)
  print 'divB:' ; minmax(divB)

  plt.figure(figN) ; plt.clf() ; figN += 1
  plt.imshow(rho[1,::,::], vmin=0.05, vmax=0.5, origin='lower')
  plt.colorbar() ; plt.title(r'$\rho$ ('+model+')')
  plt.savefig('rho-'+model+'-'+str(nout).zfill(3)+'.png', bbox_inches='tight')

  plt.figure(figN) ; plt.clf() ; figN += 1
  plt.imshow(Temp[1,::,::], origin='lower', vmin=0.15, vmax=1.24)
  plt.colorbar() ; plt.title('Temp  ('+model+')')
  plt.savefig('Temp-'+model+'-'+str(nout).zfill(3)+'.png', bbox_inches='tight')

  plt.figure(figN) ; plt.clf() ; figN += 1
  plt.imshow(divB[1,::,::], origin='lower')
  plt.colorbar() ; plt.title(r' $\nabla \cdot B$ '+'('+model+')')
  plt.savefig('divB-'+model+'-'+str(nout).zfill(3)+'.png', bbox_inches='tight')


