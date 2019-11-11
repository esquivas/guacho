#!/usr/bin/python
import fortranfile
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d

nxmap=300
nymap=300
nvmap=250
endian='<'      # little endian
kind='d'        # double precision


#time=np.zeros(39)
time=np.zeros(74)
#intensity=np.zeros(39)
intensity=np.zeros(74)

#run='0.3G-G1.01-RD-300'
run=''
#path='/datos_diable/esquivel/EXO-GUACHO/'+run+'/BIN/'
#path='/datos_europa/matias/mhd_guacho'+run+'/BIN/'
#path='/datos_europa/matias/mhd_guacho/P3/BIN/'
path='../P3/BIN/'
#plt.ion()
ff=open(path+'ly-abs.dat','w')
for nout in range(1,74):
	
	#name of the file to read
	filein=path+'LA_tau-'+str(nout).zfill(3)+'.bin'

	dx=0.15*1.49e13/real(nxmap) #  dx in cm ( 1AU = 1.49e13 cm)
	rstar= 1.1*6.955e10*2.         #  star radius in cm
	plt.ion()

	f=fortranfile.FortranFile(filein,endian=endian)
	data=f.readReals(kind).reshape(nxmap,nymap,nvmap,order='F').T

	emtaunu=np.exp(-data)
	emtau= np.sum(emtaunu,0)/real(nvmap)

	Emission=np.zeros(shape=(nxmap,nymap))
	Emission[::,::]=1e-3
	TotalEmission=0.
	Ref = 0.

	for i in range(nxmap):
	    for j in range(nxmap):
	        rad=np.sqrt( (real(i-nxmap/2)+0.5)**2+(real(j-nymap/2)+0.5)**2 )*dx 
	        if (rad <= rstar):
	            Emission[i,j]=emtau[i,j]
	            Ref = Ref + 1.
	            TotalEmission=TotalEmission + Emission[i,j]
	        else:
	            Emission[i,j]=emtau[i,j]*.3# * 1e-2
	print ''
	print 'Observed intensity: ', TotalEmission/Ref, ' from original'
	print '(the total absorptoin of the disk is', 1.-TotalEmission/Ref,')'

	time[nout]     = nout*0.05#0.2
	intensity[nout]= TotalEmission/Ref
        ff.write("intensity[nout]")
	plt.figure(3)
	plt.clf()
	#plt.imshow(Emission[120:170,120:190], origin='Lower', cmap='gist_heat', vmin=0.5, vmax=1.)
	plt.imshow(Emission[110:180,::], origin='Lower', cmap='gist_heat', vmin=0.0, vmax=1.)
	plt.colorbar(orientation='horizontal')

	plt.savefig('proj-'+str(nout).zfill(3)+'.png',dpi=300,transparent=False,bbox_inches='tight')

ff.close()
#plt.imshow(data, origin='lower', norm=LogNorm(), vmin=1.e9, vmax=1e13)
#plt.imshow(data, origin='lower', cmap='spectral', vmin=-1E9,vmax=1e9)
#plt.colorbar()
#plt.show()
