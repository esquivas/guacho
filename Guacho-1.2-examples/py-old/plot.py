import numpy as np
import pylab as pl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
#import sys

path = '/home/matis/scratch/matis/mhd_phot/kai/SS/BIN/'

input_file0 = path+'points016.bin'
#input_file1 = path+'points001.bin'
#input_file2 = path+'points002.bin'
#input_file3 = path+'points003.bin'
#input_file4 = path+'points004.bin'
#input_file5 = path+'points005.bin'

f0 = open(input_file0, 'rb')
#f1 = open(input_file1, 'rb')
#f2 = open(input_file2, 'rb')
#f3 = open(input_file3, 'rb')
#f4 = open(input_file4, 'rb')
#f5 = open(input_file5, 'rb')

data = np.fromfile(f0, count=7*64*64*64, dtype='float16')
d0 = data.reshape((7,64,64,64), order='Fortran')
f0.close()

###data = np.fromfile(f1, count=7*64*64*64, dtype='float64')
###d1 = data.reshape((7,64,64,64), order='Fortran')
###f1.close()
###
###data = np.fromfile(f2, count=7*64*64*64, dtype='float64')
###d2 = data.reshape((7,64,64,64), order='Fortran')
###f2.close()
###
###data = np.fromfile(f3, count=7*64*64*64, dtype='float64')
###d3 = data.reshape((7,64,64,64), order='Fortran')
###f3.close()
###
###data = np.fromfile(f4, count=7*64*64*64, dtype='float64')
###d4 = data.reshape((7,64,64,64), order='Fortran')
###f4.close()
###
###data = np.fromfile(f5, count=7*64*64*64, dtype='float64')
###d5 = data.reshape((7,64,64,64), order='Fortran')
###f5.close()
'''
pl.figure()
pl.imshow(d5[4,:,:,32])
pl.show()
'''
print np.amax(d0[4,:,:,:]),np.amin(d0[4,:,:,:])
print d0
#print np.amax(d5[1,:,:,:]),np.amin(d5[1,:,:,:])
#print np.amax(d5[2,:,:,:]),np.amin(d5[2,:,:,:])
#print np.amax(d5[3,:,:,:]),np.amin(d5[3,:,:,:])
#print np.amax(d5[4,:,:,:]),np.amin(d5[4,:,:,:])
#print np.amax(d5[5,:,:,:]),np.amin(d5[5,:,:,:])
'''
fig = pl.figure(figsize=(16.0,5.0))
grid = ImageGrid(fig, 111, # similar to subplot(111)
                nrows_ncols = (3, 8),
                axes_pad=0.1, # pad between axes in inch.
                cbar_mode='single',
                )

im=grid[0].imshow(d1[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid.cbar_axes[0].colorbar(im)
grid[1].imshow(d2[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[2].imshow(d3[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[3].imshow(d4[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[4].imshow(d5[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[5].imshow(d6[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[6].imshow(d7[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[7].imshow(d8[0,:,:],norm=LogNorm(vmin=100,vmax=100000))
grid[8].imshow(d1[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[9].imshow(d2[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[10].imshow(d3[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[11].imshow(d4[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[12].imshow(d5[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[13].imshow(d6[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[14].imshow(d7[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[15].imshow(d8[:,0,:],norm=LogNorm(vmin=100,vmax=100000))
grid[16].imshow(d1[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[17].imshow(d2[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[18].imshow(d3[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[19].imshow(d4[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[20].imshow(d5[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[21].imshow(d6[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[22].imshow(d7[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
grid[23].imshow(d8[:,:,0],norm=LogNorm(vmin=100,vmax=100000))
pl.show()
'''

'''
pl.figure(figsize=(16.0,5.0))
pl.subplot(3,8,1)
pl.imshow(d1[29,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,2)
pl.imshow(d2[29,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,3)
pl.imshow(d3[29,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,4)
pl.imshow(d4[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,5)
pl.imshow(d5[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,6)
pl.imshow(d6[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,7)
pl.imshow(d7[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,8)
pl.imshow(d8[14,:,:],norm=LogNorm(vmin=100,vmax=100000))

pl.colorbar()

pl.subplot(3,8,9)
pl.imshow(d1[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,10)
pl.imshow(d2[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,11)
pl.imshow(d3[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,12)
pl.imshow(d4[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,13)
pl.imshow(d5[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,14)
pl.imshow(d6[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,15)
pl.imshow(d7[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,16)
pl.imshow(d8[:,14,:],norm=LogNorm(vmin=100,vmax=100000))

pl.colorbar()

pl.subplot(3,8,17)
pl.imshow(d1[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,18)
pl.imshow(d2[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,19)
pl.imshow(d3[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,20)
pl.imshow(d4[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,21)
pl.imshow(d5[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,22)
pl.imshow(d6[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,23)
pl.imshow(d7[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.subplot(3,8,24)
pl.imshow(d8[:,:,14],norm=LogNorm(vmin=100,vmax=100000))

pl.colorbar()

pl.show()
'''

#pl.imshow(d8[:,:,14])
#pl.imshow(d8[:,14,:])
#pl.imshow(d8[14,:,:])



#pl.imshow(dh[:,:,14])
#pl.show()

#pl.semilogy(data1D[:,0],d8[:,14,14])

'''
pl.plot(data1D[:,0],d8[:,14,14])

pl.show()
'''

'''
pl.figure()

pl.plot(mesh_file1[:],d1[0,0,:])
pl.plot(mesh_file1[:],d2[0,0,:])
pl.plot(mesh_file1[:],d3[0,0,:])
pl.plot(mesh_file1[:],d4[0,0,:])
pl.plot(mesh_file1[:],d5[0,0,:])
pl.plot(mesh_file1[:],d6[0,0,:])
pl.plot(mesh_file1[:],d7[0,0,:])
pl.plot(mesh_file1[:],d8[0,0,:])
pl.show()
'''

'''
pl.plot(data1D[:,0],d1[14,14,:])
pl.plot(data1D[:,0],d2[14,14,:])
pl.plot(data1D[:,0],d3[14,14,:])
pl.plot(data1D[:,0],d4[14,14,:])
pl.plot(data1D[:,0],d5[14,14,:])
pl.plot(data1D[:,0],d6[14,14,:])
pl.plot(data1D[:,0],d7[14,14,:])
pl.plot(data1D[:,0],d8[14,14,:])


pl.plot(data1D[:,0],d1[46,46,:])
pl.plot(data1D[:,0],d2[46,46,:])
pl.plot(data1D[:,0],d3[46,46,:])
pl.plot(data1D[:,0],d4[46,46,:])
pl.plot(data1D[:,0],d5[46,46,:])
pl.plot(data1D[:,0],d6[46,46,:])
pl.plot(data1D[:,0],d7[46,46,:])
pl.plot(data1D[:,0],d8[46,46,:])
pl.title('hihi')
'''
#pl.plot(data1D[:,0],d7[14,14,:])
#pl.plot(data1D[:,0],d8[14,14,:])
#pl.plot(data1D[:,0],d1[14,:,14])
#pl.plot(data1D[:,0],d1[:,14,14])

#pl.show()
