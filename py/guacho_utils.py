#!/usr/bin/python
import numpy as np
import struct

'''
  Reads headers ASCII and BINARY
'''
def read_header(file_in, verbose=True):
  f = open(file_in,'rb')
  while 1==1 :
      s = f.readline().rstrip()
      if s == b'\xff' :
        break
      if (verbose) :
        print (s)
  f_kind = struct.unpack('s',f.read(1))[0]
  nx, ny, nz = struct.unpack('3i',f.read(12))
  if (f_kind == b'd'):
    dx, dy, dz = struct.unpack('3d',f.read(24))
  if (f_kind == b'f'):
    dx, dy, dz = struct.unpack('3f',f.read(12))
  x0, y0, z0 = struct.unpack('3i',f.read(12))
  mpi_x, mpi_y, mpi_z = struct.unpack('3i',f.read(12))
  neqs   = struct.unpack('1i',f.read(4))[0]
  neqdyn = struct.unpack('1i',f.read(4))[0]
  nghost = struct.unpack('1i',f.read(4))[0]
  if (f_kind == b'd'):
    rsc, vsc, rhosc = struct.unpack('3d',f.read(24))
    cv              = struct.unpack( 'd',f.read(8))
  if (f_kind == b'f'):
    rsc, vsc, rhosc = struct.unpack('3f',f.read(12))
    cv              = struct.unpack( 'f',f.read(4))
  return (f,f_kind,(nx,ny,nz),(dx,dy,dz),(x0,y0,z0),\
         (mpi_x,mpi_y,mpi_z),neqs,neqdyn,nghost,(rsc,vsc,rhosc),cv)

  '''
    converts conserved variables to primitives (and de-scales to cgs)
  '''
def u2prim(file_in, ublock, equation, mhd = False, entropy=False, scale=True, conserved=False) :
  head_info = read_header(file_in,verbose=False)
  if scale:
      rsc, vsc, rhosc = head_info[9]
  else:
       rsc, vsc, rhosc = 1., 1., 1.
  cv              = head_info[10]
  pblock = ublock
  if equation == 0 :
    #  density
    pblock[0,::] = pblock[0,::] * rhosc
  elif (equation >=1 and equation <= 3) :
    #  velocity components
    pblock[equation,::] = pblock[equation,::] / pblock[0,::] * vsc
  elif equation == 4 :
    #  thermal pressure
    pblock[1:4,::] = ublock[1:4,::] / ublock[0,::]
    pblock[4,::]   = ( ublock[4,::] - 0.5 * pblock[0,::]*( pblock[1,::]**2
                                                       + pblock[2,::]**2
                                                       + pblock[3,::]**2 ) ) /cv
    if mhd :
      pblock[4,::] = pblock[4,::] - 0.5 *( pblock[5,::]**2
                                          +pblock[6,::]**2
                                          +pblock[7,::]**2 )  /cv
    pblock[4,::] = pblock[4,::] * rhosc * vsc**2
  if (mhd and equation >= 5 and equation <= 7) :
    Bsc = np.sqrt(4.*np.pi*rhosc*vsc**2)
    pblock[equation,::] = pblock[equation,::] * Bsc
  if (entropy and equation == 9):
    pblock[equation,::] = ublock[9,::] * (pblock[0,::]**(1/cv[0]))* rhosc * vsc**2
  return pblock[equation,::]

'''
  Returns a 3D array of a single block for the selected equation
  consider that the index in python is that used in fortran -1
  e.g. density corresponds to equation 0
'''
def readbin3d_block(file_in, equation, verbose=False, mhd = False, entropy=False, scale=True, conserved = False):
  head_info = read_header(file_in,verbose=verbose)
  f                   = head_info[0]
  f_kind              = head_info[1]
  nx, ny, nz          = head_info[2]
  neqs                = head_info[6]
  nghost              = head_info[8]
  data = np.fromfile(f, dtype=f_kind, \
    count=(nx+2*nghost)*(ny+2*nghost)*(nz+2*nghost)*neqs)\
    .reshape(neqs,nx+2*nghost,ny+2*nghost,nz+2*nghost,order='F')
  f.close()
  if (conserved):
    primit = data[equation,nghost:(nx+nghost),nghost:(ny+nghost),nghost:(nz+nghost)]
  else :
    primit = u2prim(file_in,data[::,nghost:(nx+nghost),nghost:(ny+nghost),nghost:(nz+nghost)],equation, mhd = mhd, entropy = entropy, scale=scale, conserved=conserved)
  return primit

'''
  Returns the 3D array for equation neq in all the domain
'''
def readbin3d_all(nout,neq,path='',base='points',verbose=False, mhd=False, entropy=False, scale=True, conserved = False):

  print ('Retrieving 3D map for eqn',neq)
  file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
  head_info = read_header(file_in,verbose=False)
  f                   = head_info[0]
  mpi_x, mpi_y, mpi_z = head_info[5]
  f.close()
  proc= 0
  for ip in range(mpi_x) :
    for jp in range(mpi_y) :
      for kp in range(mpi_z) :
        file_in = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
        print (file_in)
        head_info = read_header(file_in,verbose=False)
        f                   = head_info[0]
        nx, ny, nz          = head_info[2]
        x0, y0, z0          = head_info[4]
        f.close()
        if proc ==0 :
          map3d=np.zeros(shape=(nx*mpi_x,ny*mpi_y,nz*mpi_z))
        proc = proc+1
        map3d[x0:x0+nx,y0:y0+ny,z0:z0+nz] = readbin3d_block(file_in, neq, verbose=verbose, mhd= mhd, entropy = entropy, scale=scale, conserved=conserved)
  return map3d.T

'''
  Returns the axis of the simulation box
'''
def get_axis(nout,path='',base='points',verbose=False):

  file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
  head_info = read_header(file_in,verbose=False)
  f                   = head_info[0]
  f_kind              = head_info[1]
  nx, ny, nz          = head_info[2]
  dx, dy, dz          = head_info[3]
  x0, y0, z0          = head_info[4]
  mpi_x, mpi_y, mpi_z = head_info[5]
  f.close()
  x_axis = np.linspace(0.5, nx*mpi_x-0.5 , nx*mpi_x)*dx
  y_axis = np.linspace(0.5, ny*mpi_y-0.5 , ny*mpi_y)*dy
  z_axis = np.linspace(0.5, nz*mpi_z-0.5 , nz*mpi_z)*dz
  return (x_axis, y_axis, z_axis)

'''
  Returns the extent of the computational box
'''
def get_extent(nout,path='',base='points',verbose=False, centered=False):

  file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
  head_info = read_header(file_in,verbose=False)
  f                   = head_info[0]
  f_kind              = head_info[1]
  nx, ny, nz          = head_info[2]
  dx, dy, dz          = head_info[3]
  x0, y0, z0          = head_info[4]
  mpi_x, mpi_y, mpi_z = head_info[5]
  f.close()
  if centered :
      x_extent = np.asarray( [ -(nx*mpi_x-0.5)/2, (nx*mpi_x-0.5)/2 ] )*dx
      y_extent = np.asarray( [ -(ny*mpi_y-0.5)/2, (ny*mpi_y-0.5)/2 ] )*dy
      z_extent = np.asarray( [ -(nz*mpi_z-0.5)/2, (nz*mpi_z-0.5)/2 ] )*dz
  else :
      x_extent = np.asarray( [ 0.5, nx*mpi_x-0.5 ] )*dx
      y_extent = np.asarray( [ 0.5, ny*mpi_y-0.5 ] )*dy
      z_extent = np.asarray( [ 0.5, nz*mpi_z-0.5 ] )*dz
  return ( x_extent, y_extent, z_extent)

'''
  Returns the box size
'''
def get_boxsize(nout,path='',base='points',verbose=False):

  file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
  head_info = read_header(file_in,verbose=False)
  nx, ny, nz          = head_info[2]
  mpi_x, mpi_y, mpi_z = head_info[5]
  return (nx*mpi_x, ny*mpi_y, nz*mpi_z)

'''
  Returns gamma = (cv + 1) / cv
'''
def get_gamma_ad(nout,path='',base='points',verbose=False):

  file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
  head_info = read_header(file_in,verbose=False)
  cv = head_info[10][0]
  gamma = (cv +1.0)  / cv
  return gamma

'''
  Returns the basic scalings
'''
def get_scalings(nout,path='',base='points',verbose=False):

  file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
  head_info = read_header(file_in,verbose=False)
  rsc, vsc, rhosc = head_info[9]
  Bsc = np.sqrt(4.*np.pi*rhosc*vsc**2)
  return (rsc, vsc, rhosc, Bsc)

'''
  Returns the number of cores read from HD output
'''
def get_Nproc(nout,path='',base='points',verbose=False):

  file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
  head_info = read_header(file_in,verbose=False)
  mpi_x, mpi_y, mpi_z = head_info[5]
  n_proc = mpi_x*mpi_y*mpi_z
  return n_proc

'''
 Returns a 2D cut perpenticular to the x, y or z axes (cut ==1, 2, 3, respectively)
 pos denotes the position of the cut in cells, neq the equation to be retrieved
'''
def get_2d_cut(cut,pos,nout,neq,path='',base='points',verbose=False,mhd=False, entropy=False, scale=True, conserved = False):

  file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
  head_info = read_header(file_in,verbose=False)
  f                   = head_info[0]
  nx, ny, nz          = head_info[2]
  x0, y0, z0          = head_info[4]
  mpi_x, mpi_y, mpi_z = head_info[5]
  f.close()
  nxtot = nx * mpi_x
  nytot = ny * mpi_y
  nztot = nz * mpi_z
  if cut == 1 :
      map2d=np.zeros(shape=(nytot,nztot))
      print ('YZ cut, equation: ', neq)
  elif cut == 2 :
      map2d=np.zeros(shape=(nxtot,nztot))
      print ('XZ cut, equation: ', neq)
  elif cut == 3 :
      map2d=np.zeros(shape=(nxtot,nytot))
      print ('XY cut, equation: ', neq)
  proc=0
  for ip in range(mpi_x) :
    for jp in range(mpi_y) :
      for kp in range(mpi_z) :

        if cut == 1 :
          if ( (pos >= ip*nx) and (pos <= (ip+1)*nx-1)) :
            file_in = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
            print (file_in)
            head_info = read_header(file_in,verbose=False)
            f                   = head_info[0]
            nx, ny, nz          = head_info[2]
            x0, y0, z0          = head_info[4]
            f.close()
            offset=pos-ip*nx
            block = readbin3d_block(file_in, neq, verbose=verbose, mhd= mhd, entropy = entropy, scale = scale, conserved=conserved )
            map2d[y0:y0+ny,z0:z0+nz]= block[offset,::,::]

        elif cut == 2:
          if ( (pos >= jp*ny) and (pos <= (jp+1)*ny-1)) :
            file_in = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
            print (file_in)
            head_info = read_header(file_in,verbose=False)
            f                   = head_info[0]
            nx, ny, nz          = head_info[2]
            x0, y0, z0          = head_info[4]
            f.close()
            offset=pos-jp*ny
            block = readbin3d_block(file_in, neq, verbose=verbose, mhd= mhd, entropy = entropy, scale = scale, conserved=conserved )
            map2d[x0:x0+nx,z0:z0+nz]= block[::,offset,::]

        elif cut == 3:
          if ( (pos >= kp*nz) and (pos <= (kp+1)*nz-1)) :
            file_in = path+base+str(proc).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
            print (file_in)
            head_info = read_header(file_in,verbose=False)
            f                   = head_info[0]
            nx, ny, nz          = head_info[2]
            x0, y0, z0          = head_info[4]
            f.close()
            offset=pos-kp*nz
            block = readbin3d_block(file_in, neq, verbose=verbose, mhd= mhd, entropy = entropy, scale = scale, conserved=conserved )
            map2d[x0:x0+nx,y0:y0+ny]= block[::,::,offset]

        proc = proc +1

  return  map2d.T

'''
   prints minumum and maximum from numpy array
'''
def minmax(q):
    print('min=',q.min(),' max=', q.max())


'''
   reads tracer particles
'''
def read_lmp(nout,path='', base='lmp',trim=True, verbose=False):
    #  read first output to determine N_MP and NP
    file_in = path+base+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
    f = open(file_in,'rb')
    nproc, n_mp, n_act, n_bins = struct.unpack('4i',f.read(16))
    x  = np.zeros(nproc*n_mp)
    y  = np.zeros(nproc*n_mp)
    z  = np.zeros(nproc*n_mp)
    #vx = np.zeros(nproc*n_mp)
    #vy = np.zeros(nproc*n_mp)
    #vz = np.zeros(nproc*n_mp)
    id = np.zeros(nproc*n_mp, dtype = 'i4')-1
    if (n_bins > 0) :
      SED = np.zeros( shape=(nproc*n_mp,n_bins,2) )
      r  = np.zeros(nproc*n_mp)
      th = np.zeros(nproc*n_mp)
      P1  = np.zeros( shape=(nproc*n_mp,8) )
      P2  = np.zeros( shape=(nproc*n_mp,8) )
    f.close()
    #  now reopen files one by one and read data onto the arrays
    for ip in range(nproc):
        file_in = path+base+str(ip).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
        f = open(file_in,'rb')
        nproc,n_mp,n_activeMP, n_bins = struct.unpack('4i',f.read(16))
        if verbose: print(n_activeMP, ' active particles in processor ', ip)
        for i_mp in range(n_activeMP):
            ii = struct.unpack('1i',f.read(4))[0]
            ii -= 1
            id[ii] = ii
            x [ii], y [ii], z [ii] = struct.unpack('3d',f.read(24))
            #vx[ii], vy[ii], vz[ii] = struct.unpack('3d',f.read(24))
            if (n_bins  > 0):
              r[ii], th[ii] = struct.unpack('2d',f.read(16))
              for i_bin in range (n_bins):
                SED[ii,i_bin,0]=struct.unpack('1d',f.read(8))[0]
              for i_bin in range (n_bins):
                SED[ii,i_bin,1]=struct.unpack('1d',f.read(8))[0]
              for i_eq in range(8):
                P1[ii,i_eq] = struct.unpack('1d', f.read(8) )[0]
                P2[ii,i_eq] = struct.unpack('1d', f.read(8) )[0]
        f.close()
    if (n_bins > 0):
        indices = np.argsort(id)
        SED   = np.array(SED[indices,:,:])
        P1    = np.array(P1[indices,:])
        P2    = np.array(P2[indices,:])
        array = np.array ( sorted(zip(id,x,y,z,r,th)) )
        #array = np.array ( sorted(zip(id,x,y,z,vx,vy,vz)) )
    else:
        array = np.array ( sorted(zip(id,x,y,z)) )
        #array = np.array ( sorted(zip(id,x,y,z,vx,vy,vz)) )
    if (trim) :
        #  trim
        n_mp = np.size(np.where(array[:,0] < 0 ))
        array = array[n_mp::,:]
        if (n_bins > 0):
            SED   = SED  [n_mp::,:]
            P1    = P1[n_mp::,:]
            P2    = P2[n_mp::,:]
    if (n_bins > 0):
        return array, SED, P1, P2
    else:
        return array

'''
   Converts the conserved variables to primitives for the
   relativistic case with an EOS of ideal gas with constant
   gamma, on a single cell
'''
def u2prim_rel_ideal(D, mx, my, mz, E, gamma_ad):
  #  this are the expressions in Ryu et al. 2006)
  M  = np.sqrt(mx**2 + my**2 + mz**2 )

  if (M > 0.0):
    denom = (M**2 + D**2)*(gamma_ad - 1.0)**2

    b1 =-2.0*gamma_ad*(gamma_ad - 1.0)*M*E / denom
    b2 = ( (gamma_ad**2)*(E**2) + 2.0*(gamma_ad-1.0)*(M**2) - ((gamma_ad-1.0)**2)*(D**2) ) /denom
    b3 =-2.0*gamma_ad*M*E / denom
    b4 = M**2 / denom

    a1 = -b2
    a2 =  b1*b3 -4.0*b4
    a3 =  4.0*b2*b4 - b3**2 -b1**2*b4

    R = ( 9.0*a1*a2 -27.0*a3 -2.0*a1**3 ) / 54.0
    S = ( 3.0*a2 -a1**2 ) / 9.0
    T = R**2 + S**3

    x1 =  ( R + np.sqrt(T) )**(1/3) + ( R - np.sqrt(T) )**(1/3) - a1/3.0

    B = 0.5*(b1 + np.sqrt(b1**2 - 4.0*b2 + 4.0*x1) )
    C = 0.5*(x1 - np.sqrt(x1**2 - 4.0*b4) )

    v = 0.5* ( -B + np.sqrt(B**2 - 4.0*C) )
    gamma_rel = 1.0/np.sqrt(1.0 - v**2 )
    rho = D / gamma_rel
    vx  = mx*v/M
    vy  = my*v/M
    vz  = mz*v/M
  else:
    v = 0.0
    gamma_rel = 1.00
    rho = D / gamma_rel
    vx = 0.0
    vy = 0.0
    vz = 0.0

  Pth = (gamma_ad-1.0)*(E - mx*vx -my*vy -mz*vz - rho)

  return rho, vx, vy, vz, Pth

'''
   Computes the mass-loss/accretion rate of a spherical region
   of radius Rs centered at a position x0,y0,z0
   Position is measured with respect of a corner of the simulation
   I/O is in [ cgs ]
   nout,neq,path='',base='points'):

'''
def get_Mdot(nout, x0, y0, z0, R0, path='', base='points', verbose=False):
  #  get boxsize
  Nx, Ny, Nz = get_boxsize(nout, path+'BIN/', verbose=False)
  #  read scalings
  rsc, vsc, rhosc, Bsc = get_scalings(nout,path=path+'BIN/', verbose=False)
  #  from header read spacings
  header = read_header(path+'BIN/points'+str(0).zfill(3)+'.'+str(nout).zfill(3)+'.bin', verbose=False)
  dx, dy, dz  = header[3]
  dx *= rsc ; dy *= rsc ; dz *= rsc  # [ cgs ]

  #  locate bounds
  #  left
  iL = int( (x0 - 1.05*R0)/dx  + 0.5)
  iH = int( (x0 + 1.05*R0)/dx  + 0.5)
  jL = int( (y0 - 1.05*R0)/dy  + 0.5)
  jH = int( (y0 + 1.05*R0)/dy  + 0.5)
  kL = int( (z0 - 1.05*R0)/dz  + 0.5)
  kH = int( (z0 + 1.05*R0)/dz  + 0.5)

  Cells_in  = 0
  Cells_out = 0
  Fx = 0 ; Fy = 0 ; Fz = 0
  #  loop the domain
  for l in range(1,4):   #  l--> x, y, z directions

    flux   = readbin3d_all(nout=nout,neq=l,path=path+'BIN/',verbose=False,conserved=True,scale=False)

    for i in range(iL,iH) :
      for j in range(jL,jH) :
        for k in range(kL,kH) :

          x = float(i+0.5)*dx
          y = float(j+0.5)*dx
          z = float(k+0.5)*dx

          rad = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2 )

          if (rad <= R0) :    # inside the Gaussian Surface

            if l == 1 :       #  x direction
              Fx += ( flux[k,j,i+1] - flux[k,j,i-1])

            elif l==2 :       # y direction
              Fy += ( flux[k,j+1,i] - flux[k,j-1,i])

            elif l==3 :       # z direction
              Fz += ( flux[k+1,j,i] - flux[k-1,j,i])

            Cells_in  += 1
          else:
            Cells_out += 1

  Fsc = rhosc*vsc           # [scaling factor]
  dA  = dx*dy               # [ cgs ]

  Mdx = Fx * Fsc * dA       # [cgs]
  Mdy = Fy * Fsc * dA       # [cgs]
  Mdz = Fz * Fsc * dA       # [cgs]
  Mdot = Mdx + Mdy + Mdz
  print(f'{Cells_in} cells inside region, {Cells_out}, outside')

  return Mdot
