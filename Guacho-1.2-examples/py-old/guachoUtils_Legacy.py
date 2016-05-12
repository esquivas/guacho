#!/usr/bin/python
import fortranfile
import numpy as np

'''
 Returns a 2D cut perpenticular to the x, y or z axes (cut ==1, 2, 3, respectively)
 pos denotes the position of the cut in cells, neq the equation to be retrieved
 By default it reads double precision numbers, if 4byte floats are employed set kind='f'
 It assumes 2 ghost cells, if other number is used, set it as the optional parameter
 By default uses little endian, if data is written in big endian set endian='>'
'''
def coplot3d(cut,pos,neq,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout, \
             path='',kind='d',nghost=2,endian='<',base='points'):

    nx=nxtot/mpiX
    ny=nytot/mpiY
    nz=nztot/mpiZ
    
    if cut == 1 :
        map=np.zeros(shape=(nytot,nztot))
        print 'YZ cut, equation: ', neq
    elif cut == 2 :
        map=np.zeros(shape=(nxtot,nztot))
        print 'XZ cut, equation: ', neq
    elif cut == 3 : 
        map=np.zeros(shape=(nxtot,nytot))
        print 'XY cut, equation: ', neq
        
    proc=0

    for ip in range(mpiX) :
        for jp in range(mpiY) :
            for kp in range(mpiZ) :
                
                if cut == 1 :
                    if ( (pos >= ip*nx) and (pos <= (ip+1)*nx-1)) :
                        filein = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
#                        filein = path+base+str(nout).zfill(3)+'.bin'
                        print filein
                        f=fortranfile.FortranFile(filein,endian=endian)
                        data=f.readReals(kind).reshape(neqs,nx+2*nghost,ny+2*nghost,nz+2*nghost,order='F')
                        offset=pos-ip*nx+nghost
                        map[jp*ny:(jp+1)*ny,kp*nz:(kp+1)*nz]=data[neq,offset,nghost:(ny+nghost),nghost:(nz+nghost)]
                        
                elif cut == 2:
                    if ( (pos >= jp*ny) and (pos <= (jp+1)*ny-1)) :
#                        filein = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
                        filein = path+base+str(nout).zfill(3)+'.bin'
                        print filein
                        f=fortranfile.FortranFile(filein,endian=endian)
                        data=f.readReals(kind).reshape(neqs,nx+2*nghost,ny+2*nghost,nz+2*nghost,order='F')
                        offset=pos-jp*ny+nghost
                        map[ip*nx:(ip+1)*nx,kp*nz:(kp+1)*nz]=data[neq,nghost:(nx+nghost),offset,nghost:(nz+nghost)]
                        
                elif cut == 3:
                    if ( (pos >= kp*nz) and (pos <= (kp+1)*nz-1)) :
#                        filein = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
                        filein = path+base+str(nout).zfill(3)+'.bin'
                        print filein
                        f=fortranfile.FortranFile(filein,endian=endian)
                        data=f.readReals(kind).reshape(neqs,nx+2*nghost,ny+2*nghost,nz+2*nghost,order='F')
                        offset=pos-kp*nz+nghost
                        map[ip*nx:(ip+1)*nx,jp*ny:(jp+1)*ny]=data[neq,nghost:(nx+nghost),nghost:(ny+nghost),offset]
                proc=proc+1
                
    return map.T

'''
   Returns a temperature cut
'''
def coplotTemp(cut,pos,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout, \
             path='',kind='d',nghost=2,endian='<', cv=1.5, Tempsc=1.,base='points',neqdyn=5):


    nx=nxtot/mpiX
    ny=nytot/mpiY
    nz=nztot/mpiZ
    
    if cut == 1 :
        map=np.zeros(shape=(nytot,nztot))
        print 'YZ Temperature cut'
    elif cut == 2 :
        map=np.zeros(shape=(nxtot,nztot))
        print 'XZ  Temperature cut'
    elif cut == 3 : 
        map=np.zeros(shape=(nxtot,nytot))
        print 'XY  Temperature cut'
        
    proc=0

    for ip in range(mpiX) :
        for jp in range(mpiY) :
            for kp in range(mpiZ) :
                
                if cut == 1 :
                    if ( (pos >= ip*nx) and (pos <= (ip+1)*nx-1)) :
                        filein = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
#                        filein = path+base+str(nout).zfill(3)+'.bin'
                        print filein
                        f=fortranfile.FortranFile(filein,endian=endian)
                        data=f.readReals(kind).reshape(neqs,nx+2*nghost,ny+2*nghost,nz+2*nghost,order='F')
                        offset=pos-ip*nx+nghost
                        data2d=data[::,offset,nghost:(ny+nghost),nghost:(nz+nghost)]
                                                
                        pgas= (data2d[4,::,::]-0.5*(np.square(data2d[1,::,::])+np.square(data2d[2,::,::])+np.square(data2d[3,::,::])) / data2d[0,::,::] )/cv
                        dentot=2.*data2d[0,::,::]-data2d[neqdyn,::,::]
                        temp=pgas/dentot

                        
                        map[jp*ny:(jp+1)*ny,kp*nz:(kp+1)*nz]=temp[::,::]*Tempsc
                        
                elif cut == 2:
                    if ( (pos >= jp*ny) and (pos <= (jp+1)*ny-1)) :
                        filein = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
#                        filein = path+base+str(nout).zfill(3)+'.bin'
                        print filein
                        f=fortranfile.FortranFile(filein,endian=endian)
                        data=f.readReals(kind).reshape(neqs,nx+2*nghost,ny+2*nghost,nz+2*nghost,order='F')
                        offset=pos-jp*ny+nghost
                        data2d=data[::,nghost:(nx+nghost),offset,nghost:(nz+nghost)]

                        pgas= (data2d[4,::,::]-0.5*(np.square(data2d[1,::,::])+np.square(data2d[2,::,::])+np.square(data2d[3,::,::])) / data2d[0,::,::] )/cv
                        dentot=2.*data2d[0,::,::]-data2d[neqdyn,::,::]
                        temp=pgas/dentot
                        
                        map[ip*nx:(ip+1)*nx,kp*nz:(kp+1)*nz]=temp[::,::]*Tempsc
                        
                elif cut == 3:
                    if ( (pos >= kp*nz) and (pos <= (kp+1)*nz-1)) :
                        filein = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
#                        filein = path+base+str(nout).zfill(3)+'.bin'
                        print filein
                        f=fortranfile.FortranFile(filein,endian=endian)
                        data=f.readReals(kind).reshape(neqs,nx+2*nghost,ny+2*nghost,nz+2*nghost,order='F')
                        offset=pos-kp*nz+nghost
                        data2d=data[::,nghost:(nx+nghost),nghost:(ny+nghost),offset]

                        pgas= (data2d[4,::,::]-0.5*(np.square(data2d[1,::,::])+np.square(data2d[2,::,::])+np.square(data2d[3,::,::])) / data2d[0,::,::] )/cv
                        dentot=2.*data2d[0,::,::]-data2d[neqdyn,::,::]
                        temp=pgas/dentot
                        
                        map[ip*nx:(ip+1)*nx,jp*ny:(jp+1)*ny]=temp[::,::]*Tempsc
                proc=proc+1
                
    return map.T

             


'''
    returns the 3D array for equation neq
'''
def readbin3d(neq,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,\
              path='',kind='d',nghost=2,endian='<',base='points'):

    nx=nxtot/mpiX
    ny=nytot/mpiY
    nz=nztot/mpiZ

    proc=0
    
    map=np.zeros(shape=(nxtot,nytot,nztot))
    print 'Retrieving 3D map for eqn',neq
    
    for ip in range(mpiX) :
        for jp in range(mpiY) :
            for kp in range(mpiZ) :    
#                filein = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
                filein = path+base+str(nout).zfill(3)+'.bin'
                print filein
                f=fortranfile.FortranFile(filein,endian=endian)
                data=f.readReals(kind).reshape(neqs,nx+2*nghost,ny+2*nghost,nz+2*nghost,order='F')

                map[ip*nx:(ip+1)*nx,jp*ny:(jp+1)*ny,kp*nz:(kp+1)*nz]=\
                data[neq,nghost:(nx+nghost),nghost:(ny+nghost),nghost:(nz+nghost)]

                proc=proc+1
    return map.T

'''
    ?
'''
