#!/usr/bin/python
import numpy as np
import pylab
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as mpl
mpl.rc('text',usetex=True)
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
#plt.ion()

data = np.fromfile('BIN/blast-009.bin', dtype=np.dtype('d'),count=-1).reshape(510*2,510*2,4,order='C')

fig= plt.imshow(data[::,::,0],cmap='Spectral',origin='lower',vmin=0.05
                ,extent=(0,1,0,1), vmax=2.,norm=LogNorm() )
plt.colorbar(shrink=0.99,extend='both')#.add_lines(con)

X,Y = np.meshgrid( np.arange(0.,1.,1./64.),np.arange(0.,1.,1./64) )

U=data[::16,::16,1]/data[::16,::16,0]
V=data[::16,::16,2]/data[::16,::16,0]
Q=plt.quiver(X,Y,U,V, pivot='mid',scale=10.)

plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
#plt.title(r'Pressure')
plt.title(r'density') 


#plt.savefig('./contour_ex.pdf',dpi=300,transparent=True,bbox_inches='tight')
plt.ioff()

plt.show()


