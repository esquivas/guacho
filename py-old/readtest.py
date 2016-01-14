#!/usr/bin/python
import numpy as np
import struct

path = "/Users/esquivel/Desktop/Guacho-Testing/Temp/"
file_in = "test.bin"

f = open(path+file_in,'rb')
s = f.readline().rstrip()
print s

while s != "-end of header-" :
    s = f.readline().rstrip()
    print s

nx, ny, nz = struct.unpack('3i',f.read(12))
dx, dy, dz = struct.unpack('3f',f.read(12))
x0, y0, z0 = struct.unpack('3i',f.read(12))
mpi_x, mpi_y, mpi_z = struct.unpack('3i',f.read(12))
neqs   = struct.unpack('1i',f.read(4))[0]
nghost = struct.unpack('1i',f.read(4))[0]

nelem= (nx+2*nghost)*(ny+2*nghost)*(nz+2*nghost)*neqs

data = np.fromfile(f, dtype='f',count=nelem).reshape(nx+2*nghost,ny+2*nghost,nz+2*nghost,neqs,order='F')

f.close()
