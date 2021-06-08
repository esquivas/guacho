!=======================================================================
!> @file lyman_alpha_tau.f90
!> @brief Lyman_alpha_utilities
!> @author M. Schneiter, Alejandro Esquivel
!> @date 4/May/2016

! Copyright (c) 2016 Guacho Co-Op
!
! This file is part of Guacho-3D.
!
! Guacho-3D is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief Lyman_alpha_utilities
!> @details Utilities to compute the Lyman-@f \alpha @f$ absorption

module polfracrad_utilities

contains

!> @brief Initializes data
!> @details Initializes data, MPI and other stuff

subroutine init_polfrac()

!  Initializes MPI, data arrays, etc
use parameters
use globals, only : u, up, dx, dy, dz, coords, rank, left, right   &
                     , top, bottom, out, in, rank, comm3d
implicit none
  integer :: nps, err
  integer, dimension(0:ndim-1) :: dims
  logical, dimension(0:ndim-1) :: period

  !initializes MPI
#ifdef MPIP
#ifdef PERIODX
  logical, parameter :: perx=.true.
#else
  logical, parameter :: perx=.false.
#endif
#ifdef PERIODY
  logical, parameter :: pery=.true.
#else
  logical, parameter :: pery=.false.
#endif
#ifdef PERIODZ
  logical, parameter :: perz=.true.
#else
  logical, parameter :: perz=.false.
#endif
  period(0)=perx
  period(1)=pery
  period(2)=perz
  dims(0)  =MPI_NBX
  dims(1)  =MPI_NBY
  dims(2)  =MPI_NBZ
  
  call mpi_init (err)
  call mpi_comm_rank (mpi_comm_world,rank,err)
  call mpi_comm_size (mpi_comm_world,nps,err)
  if (nps.ne.np) then
     print*, 'processor number (',nps,') is not equal to pre-defined number (',np,')'
     call mpi_finalize(err) 
     stop
  endif
#else
  rank=0
  coords(:)=0
#endif
  if(rank.eq.master) then
  print '(a)' ,"*******************************************"
     print '(a)' ,"                        _                 *"
     print '(a)' ,"  __   _   _  __ _  ___| |__   ___    3   *"
     print '(a)' ," / _ `| | | |/ _` |/ __| '_ \ / _ \    D  *"
     print '(a)' ,"| (_| | |_| | (_| | (__| | | | (_) |      *"
     print '(a)' ," \__, |\__,_|\__,_|\___|_| |_|\___/       *"
     print '(a)' ," |___/                                    *"
     print '(a)' ,"                                          *"
  endif
#ifdef MPIP
  if(rank.eq.master) then
     print '(a,i3,a)','*    running with mpi in', np , ' processors    *'
     print '(a)' ,'*******************************************'
     print '(a)', 'Calculating Lyman Alpha Tau'
  end if
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, 1            &
       , comm3d, err)
  call mpi_comm_rank(comm3d, rank, err)
  call mpi_cart_coords(comm3d, rank, ndim, coords, err)
  print '(a,i3,a,3i4)', 'processor ', rank                              &
       ,' ready w/coords',coords(0),coords(1),coords(2)   
  call mpi_cart_shift(comm3d, 0, 1, left  , right, err)
  call mpi_cart_shift(comm3d, 1, 1, bottom, top  , err)
  call mpi_cart_shift(comm3d, 2, 1, out   , in   , err)
  call mpi_barrier(mpi_comm_world, err)   
  !
#else
  print '(a)' ,'*******************************************'
  print '(a)' ,'*     running on a single processor       *'
  print '(a)' ,'*******************************************'
  print '(a)', 'Calculating Lyman Alpha Tau'
#endif

!   grid spacing
  dx=xmax/nxtot
  dy=ymax/nytot
  dz=zmax/nztot

!   allocate big arrays in memory
allocate( u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
allocate( up(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )

end subroutine init_polfrac

!=======================================================================

!> @brief reads data from file
!> @details reads data from file
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param integer [in] itprint : number of output
!> @param string [in] filepath : path where the output is

subroutine read_data(u,itprint,filepath)

  use parameters, only : &
                       np, neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals, only : rank, comm3d
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: itprint
  character (len=128), intent(in) :: filepath
  integer :: unitin, ip, err
  character (len=128) file1
  character           :: byte_read
  character, parameter  :: lf = char(10) 
  integer :: nxp, nyp, nzp, x0p, y0p, z0p, mpi_xp, mpi_yp, mpi_zp,neqp, neqdynp, nghostp
  real :: dxp, dyp, dzp, scal(3), cvp

  
  take_turns : do ip=0,np-1
    if (rank == ip) then

#ifdef MPIP
        write(file1,'(a,i3.3,a,i3.3,a)')  &
             trim(filepath)//'BIN/points',rank,'.',itprint,'.bin'
        unitin=rank+10
#else
         write(file1,'(a,i3.3,a)')         &
              trim(filepath)//'BIN/points',itprint,'.bin'
         unitin=10
#endif
         open(unit=unitin,file=file1,status='unknown', access='stream', &
              convert='LITTLE_ENDIAN')

         !   discard the ascii header
         do while (byte_read /= achar(255) )
            read(unitin) byte_read
            !print*, byte_read
         end do

    !  read bin header, sanity check to do
         read(unitin) byte_read
         read(unitin) byte_read
         read(unitin) nxp, nyp, nzp
         read(unitin) dxp, dyp, dzp
         read(unitin) x0p, y0p, z0p
         read(unitin) mpi_xp, mpi_yp, mpi_zp
         read(unitin) neqp, neqdynp
         read(unitin) nghostp
         read(unitin) scal(1:3)
         read(unitin) cvp
         read(unitin) u(:,:,:,:)
         close(unitin)

    print'(i3,a,a)',rank,' read: ',trim(file1)

    end if
    call mpi_barrier(comm3d,err)
  end do take_turns

end subroutine read_data

!=======================================================================

!> @brief gets position of a cell
!> @details Returns the position and spherical radius calculated with
!! respect to  the center of the grid
!> @param integer [in] i : cell index in the x direction
!> @param integer [in] j : cell index in the y direction
!> @param integer [in] k : cell index in the z direction
!> @param real    [in] x : x position in the grid 
!> @param real    [in] y : y position in the grid 
!> @param real    [in] z : z position in the grid 

  subroutine getXYZ(i,j,k,x,y,z)

    use globals,    only : dx, dy, dz, coords
    use parameters, only : nx, ny, nz, nxtot, nytot, nztot
    implicit none
    integer, intent(in)  :: i, j, k
    real,    intent(out) :: x, y, z
 
    x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
    y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
    z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
        
  end subroutine getXYZ

!=======================================================================

!> @brief Rotation around the X axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid

 subroutine rotation_x(theta,x,y,z,xn,yn,zn)

   ! rotation around the x axis by an angle theta

   implicit none
   real, intent(in ) :: theta, x, y, z
   reaL, intent(out) :: xn, yn, zn
   xn =   x
   yn =   y*cos(theta) - z*sin(theta)
   zn =   y*sin(theta) + z*cos(theta)
 end subroutine rotation_x

!=======================================================================

!> @brief Rotation around the Y axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid

 subroutine rotation_y(theta,x,y,z,xn,yn,zn)

   implicit none
   real, intent(in ) :: theta, x, y, z
   real, intent(out) :: xn, yn, zn
   xn =   x*cos(theta) + z*sin(theta)
   yn =   y
   zn = - x*sin(theta) + z*cos(theta)
 end subroutine rotation_y

!=======================================================================

!> @brief Rotation around the Z axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid

 subroutine rotation_z(theta,x,y,z,xn,yn,zn)

   implicit none
   real, intent(in ) :: theta, x, y, z
   real, intent(out) :: xn, yn, zn
   xn =   x*cos(theta) - y*sin(theta)
   yn =   x*sin(theta) + y*cos(theta)
   zn =   z
 end subroutine rotation_z

!=======================================================================

!> @brief Fill target map
!> @details Fills the target map of one MPI block
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param real [in] u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax) : 
!! conserved variables
!> @param real [out] map(nxmap,mymap) : Target map
!> @param real [in] dxT : target pixel width
!> @param real [in] dyT : target pixel height
!> @param real [in] thetax : Rotation around X
!> @param real [in] thetay : Rotation around Y
!> @param real [in] thetaz : Rotation around Z

subroutine init_azar()
  implicit none
  integer :: rand_size
  integer, allocatable, dimension(:) :: rand_seed
  character (len=10) :: system_time
  real :: rtime

  call random_seed(size=rand_size)
  allocate(rand_seed(1:rand_size))
  call date_and_time(time=system_time)
  read(system_time,*) rtime
  rand_seed=int(rtime*1000.)
  call random_seed(put=rand_seed)    
  deallocate(rand_seed)

end subroutine init_azar

 
 subroutine statb(u,up,sumbx,sumbx2,sumby,sumby2,ncontb,&
      sumb_2x,sumb_2x2,sumb_2y,sumb_2y2,ncontb_2,theta_x,theta_y,theta_z)
  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         neq, nx, ny, nz, vsc2, rsc, rhosc,nztot, neqdyn, &
                         bsc,psc
  use globals, only : dz, dx, dy
  use hydro_core, only : u2prim
  use ieee_arithmetic
  implicit none
  !
!  integer, intent(in) :: nxmap,nymap,nzmap
  real, intent(in) :: u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax)
  real, intent(in) :: up(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax)
  real , intent(in) :: theta_x, theta_y, theta_z
  real, intent(out) :: sumbx,sumbx2,sumby,sumby2
  real, intent(out) :: sumb_2x,sumb_2x2,sumb_2y,sumb_2y2
  integer,intent(out)::ncontb,ncontb_2
  real :: T, prim(neq)
  integer :: i,j,k
  real :: bxp,byp,bzp,bxn,byn,bzn,pres1,pres2,r1,r2,gg,deltaS
  real :: b2xp,b2yp,b2zp
  real :: vx, vy, vz,vsc,vt
!
  sumbx=0.
  sumbx2=0.
  sumby=0.
  sumby2=0.  
  ncontb=0
  sumb_2x=0.
  sumb_2x2=0.
  sumb_2y=0.
  sumb_2y2=0.  
  ncontb_2=0

  gg=5./3.
!  gg=1.3
  
  !
  vsc=sqrt(vsc2)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           call u2prim(u(:,i,j,k),prim,T)
           pres2=prim(5)*psc
           r2=prim(1)*rhosc
           b2xp=prim(6)*bsc
           b2yp=prim(7)*bsc
           b2zp=prim(8)*bsc
           call u2prim(up(:,i,j,k),prim,T)
           pres1=prim(5)*psc
           r1=prim(1)*rhosc
           bxp=prim(6)*bsc
           byp=prim(7)*bsc
           bzp=prim(8)*bsc
           deltaS=1.5*log(pres2*r1**gg/(pres1*r2**gg))
           !deltaS=3.333*log(pres2*r1**gg/(pres1*r2**gg))
           !
           if (deltaS.ge.15.) then
!           if (deltaS.ge.32.) then              
!           if (pres2.ge.7.e-8) then
              !
              ncontb=ncontb+1
              call rotation_x(theta_x,bxp,byp,bzp,bxn,byn,bzn)
              call rotation_y(theta_y,bxn,byn,bzn,bxp,byp,bzp)
              call rotation_z(theta_z,bxp,byp,bzp,bxn,byn,bzn)
              !
              sumbx=sumbx+bxn
              sumbx2=sumbx2+bxn**2.
              sumby=sumby+byn
              sumby2=sumby2+byn**2.
!
              ncontb_2=ncontb_2+1
              call rotation_x(theta_x,b2xp,b2yp,b2zp,bxn,byn,bzn)
              call rotation_y(theta_y,bxn,byn,bzn,b2xp,b2yp,b2zp)
              call rotation_z(theta_z,b2xp,b2yp,b2zp,bxn,byn,bzn)
              !
              sumb_2x=sumb_2x+bxn
              sumb_2x2=sumb_2x2+bxn**2.
              sumb_2y=sumb_2y+byn
              sumb_2y2=sumb_2y2+byn**2.
              !
           end if
        end do
     end do
  end do 
          !          
end subroutine statb

subroutine fill_map(nxmap,nymap,nvmap,vmin,vmax,u,up,map,dxT,dyT,&
                   theta_x,theta_y,theta_z,baver,dbaver,baver2,dbaver2)

  use constants, only : clight, pi
  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         neq, nx, ny, nz, vsc2, rsc, rhosc,nztot, neqdyn, &
                         bsc,psc
  use globals, only : dz, dx, dy
  use hydro_core, only : u2prim
  use ieee_arithmetic
  implicit none

  integer, intent(in) :: nxmap,nymap,nvmap
  real, intent(in) :: vmin, vmax 
  real, intent(in) :: u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax)
  real, intent(in) :: up(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax)  
  real , intent(in) :: dxT, dyT, theta_x, theta_y, theta_z
  real, intent(out) :: map(nxmap,nymap,nvmap)
  integer :: i,j,k, iobs, jobs, ibb!, kobs
  real :: x,y,z,xn,yn,zn, vx, vy, vz,vxn, vyn, vzn, velsc, pres1,pres2
  real :: T, prim(neq), profile(nvmap)
!  real, parameter :: sigmaLA = 0.01105, lambdaLA=1.215668e-5 !(c/nu0=lambda)
  real :: alpha, f0, rp, norx, nory, norz, bt, vt, dv, bper, alphaE, valf, bper1
  !
  real :: norbx, norby, norbz, coef, cos2xi, sin2xi, cos2xipi, sin2xipi
  real :: bx, by, bz, bxn, byn, bzn, bxp, byp, bzp
  real :: bxx, byy, bzz, exx, eyy, fluxc, r, bx1, by1, bz1, bt1
  real :: cos2_B, sin2_B, fmax, cos2_B1, sin2_B1
  real :: eta,crx, xp, aa, bb, r1,r2, crho, coefp, deltaS, gg
  real :: phirx(201)
  real :: bfn(7)
  integer :: qp, ip, np, ip1, itrx
  real :: xcenter, ycenter, bnor, beast, rdis, chi, sinr, cosr
  real :: farrot, fp, gp
  real :: p_ip1,p_im1,dpx,dpy,dpz,dptot,txr,baver,dbaver,fef,baver2,dbaver2,eta2
  real :: chiref(nxmap,nymap)
  character*80 :: filein  = 'coef0.1_2.4kev.dat'! Input File
  
  velsc=sqrt(vsc2)

  !some constants
  itrx    = 0
!  alpha   = 0.6
  !  eta     = 10.
  fef = 1.e-3
  farrot  = 0. !faraday rotation in degrees
  qp      = 3 ! 0 isotropic, 1 quasiparallel, 2 quasiperpendicular, 3 hibrido 
  dv      = sqrt(dz**2+dy**2+dx**2) ! CELL VOLUME
  alphaE  = (90.d0-farrot)*pi/180.d0 !
  gg=5./3.
!  gg=1.3
  !
!
  eta=(baver/(dbaver))**2.
  eta2=(baver2/(dbaver2))**2.
  do k=1,nz
     do j=1,ny
        do i=1,nx

          !  obtain original position
          call getXYZ(i,j,k, x,y,z)
          
          !distance to the centre of the grid
          rp =sqrt( x**2+ y**2+ z**2)

!!$          ! normalized position w.r.t c. of grid
!          norx = x/rp
!          nory = y/rp
!          norz = z/rp

          !only do to the side facing the observer
!          if(z < -30.*dz ) then
            !  do the rotation of the coordinates
            call rotation_x(theta_x,x,y,z,xn,yn,zn)
            call rotation_y(theta_y,xn,yn,zn,x,y,z)
            call rotation_z(theta_z,x,y,z,xn,yn,zn)
            ! This is the position projected on 
            ! the target (centered)
            ! Integration is along Z
            iobs=int(xn/dxT+nxmap/2)
            jobs=int(yn/dyT+nymap/2)
            !kobs=zn/dzT + nzmap/2
            
            !  get the velocity in cm/s
            call u2prim(u(:,i,j,k),prim,T)
            vx=prim(2)*velsc
            vy=prim(3)*velsc
            vz=prim(4)*velsc
            pres2=prim(5) *psc
            r2= prim(1)*rhosc
            bx=prim(6)*bsc ! Should we scale?
            by=prim(7)*bsc
            bz=prim(8)*bsc
            txr=T
            bt = sqrt(bx**2+by**2+bz**2) ! DEFINIR
            !
            vt = sqrt(vx**2+vy**2+vz**2) ! DEFINIR
            call u2prim(up(:,i,j,k),prim,T)
            r1 = prim(1)*rhosc
            pres1=prim(5)*psc
            bx1=prim(6)*bsc ! Should we scale?
            by1=prim(7)*bsc
            bz1=prim(8)*bsc
            bt1 = sqrt(bx1**2+by1**2+bz1**2) ! DEFINIR
            !            crho= r/r2
            !dpres_x
            call u2prim(u(:,i+1,j,k),prim,T)
            p_ip1=prim(5)
!
            call u2prim(u(:,i-1,j,k),prim,T)
            p_im1=prim(5)
!
            dpx=(p_ip1-p_im1)/(2.d0*dx)
            !dpres_y
            call u2prim(u(:,i,j+1,k),prim,T)
            p_ip1=prim(5)
!            
            call u2prim(u(:,i,j-1,k),prim,T)
            p_im1=prim(5)
            dpy=(p_ip1-p_im1)/(2.d0*dy)

            !dpres_z
            call u2prim(u(:,i,j,k+1),prim,T)
            p_ip1=prim(5)
!            
            call u2prim(u(:,i,j,k-1),prim,T)
            p_im1=prim(5)
            dpz=(p_ip1-p_im1)/(2.d0*dz)
            !
            dptot=sqrt(dpx*dpx+dpy*dpy+dpz*dpz)
!!$            if (dptot == 0.) then
!!$               dptot=1.e-12
!!$            end if
            if (vt == 0.) then
               vt=1.e-1
            end if
!!$            norx = dpx/dptot
!!$            nory = dpy/dptot
!!$            norz = dpz/dptot
            norx = vx/vt
            nory = vy/vt
            norz = vz/vt
!
            coefp=2.1
            alpha = (coefp-1.)/2.
            f0  = Gamma(coefp/4.-1./12.)*Gamma(coefp/4.+1.25)/Gamma(coefp/4.+1.75)
            fp  = 2.**(0.5*(coefp+1.))/(coefp+1)*Gamma(coefp/4.+19./12.)*f0
            gp = 2.**(0.5*(coefp-3.))*Gamma(coefp/4.+7./12.)*f0
            !            f0      = (alpha + 1) / (alpha + 5./3.)
            !
            deltaS=1.5*log(pres2*r1**gg/(pres1*r2**gg))
!            deltaS=3.333*log(pres2*r1**gg/(pres1*r2**gg))
            if (deltaS.ge.15.) then
!            if (deltaS.ge.32.) then                              
!               fluxc = vt**(4.*alpha)*r2
               fluxc = pres2**(2.*alpha)*r2**(1.-2.*alpha)
            else
               fluxc=0.
            end if
!            fluxc = pres**(2.*alpha)*r**(1.-2.*alpha)
            if ( ieee_is_nan(fluxc) ) then
               write(*,*) 'map is NaN', fluxc, vt, r, velsc, vx, vy, vz
               stop
            endif

            !  obtain the LOS velocity
            call rotation_x(theta_x,vx,vy,vz,vxn,vyn,vzn)
            call rotation_y(theta_y,vxn,vyn,vzn,vx,vy,vz)
            call rotation_z(theta_z,vx,vy,vz,vxn,vyn,vzn)
            !          
            ! rx
            !
            IF(ITRX.EQ.0) THEN
               ITRX=1
               open(unit=10,file=trim(filein),status='unknown')
               READ(10,*) NP
               DO IP=1,NP
                  READ(10,*) AA,BB,PHIRX(IP)
               END DO
               CLOSE(UNIT=10)
            END IF
            XP=(LOG10(TXR)-4.)/0.02+1.000001
            IP=MAX(1,INT(XP))
            IP=MIN(IP,200)
            IP1=IP+1
            CRX=PHIRX(IP)+(PHIRX(IP1)-PHIRX(IP))*(XP-FLOAT(IP))
            IF(TXR.GE.1.E8) CRX=PHIRX(200)*(TXR/1.E8)**0.5
            !
            !
            !

!            valf = bt/sqrt(4.*pi*r)
!            eta = valf/(fef*vt)


            cos2_B = ((bx*norx+by*nory+bz*norz)/bt)**2.
            sin2_B = 1.d0-cos2_B
            cos2_B1 = ((bx1*norx+by1*nory+bz1*norz)/bt1)**2.
            sin2_B1 = 1.d0-cos2_B1
            !
            !  obtain the rotated magnetid field
            call rotation_x(theta_x,bx,by,bz,bxn,byn,bzn)
            call rotation_y(theta_y,bxn,byn,bzn,bx,by,bz)
            call rotation_z(theta_z,bx,by,bz,bxn,byn,bzn)
            bper = sqrt(bxn**2+byn**2)
            ! This is the B-field on the target plane

!!$            if(dbaver2.gt.abs(baver2-bper))then
!!$               eta2=(bper/(bper-baver2))**2.
!!$            end if
!!$            if(ncontb.gt.1)then

!!$            else
!!$               eta=1000.
!!$            end if
            if (qp.eq.0.) coef = 1.d0
            if (qp.eq.1.) coef = cos2_B
            if (qp.eq.2.) coef = sin2_B
            if (qp.eq.3.) coef =eta2*(cos2_B*eta2**2.+1.)/(1.+eta2**2)!*&
!            if (qp.eq.3.) coef =(cos2_B1*eta**2.+1.)/(1.+eta**2)!*&         
                 !eta2*(cos2_B1*eta2**2+1.)/(1.+eta2**2)
            !explain what this is
            !magnetic fields unit vectors (plane of the sky)
            bxx  = byn/bper
            byy  = -bxn/bper

            !electric fields unit vectors (plane of the sky)
            exx  =  cos(alphaE)*bxx + sin(alphaE)*byy
            eyy  = -sin(alphaE)*bxx + cos(alphaE)*byy

            cos2xi    = bxx*bxx-byy*byy
            sin2xi    = 2.d0*bxx*byy
            cos2xipi  = exx*exx-eyy*eyy
            sin2xipi  = 2.d0*exx*eyy


            !  calculate the line profile function
!            call phigauss(T, vzn,vmin,vmax,nvmap,profile) 
            !  make sure the result lies in the map bounds
            if( (iobs >=1    ).and.(jobs >=1    ).and. &
                (iobs <=nxmap).and.(jobs <=nymap) ) then
!              if ((T < 1e5).and.(prim(7)<0)) then
!              if ((T < 1e5)) then

               !FLUXES, 1==Q, 2==U

               !FLUXQ
               map(iobs,jobs,1)= map(iobs,jobs,1) + &
                    fluxc*coef*gp*cos2xi*bper**(alpha+1.d0)*dv
               
               !FLUXQE
               map(iobs,jobs,2)= map(iobs,jobs,2) + &                
                    fluxc*coef*gp*cos2xipi*bper**(alpha+1.d0)*dv

               !FLUXU
               map(iobs,jobs,3)= map(iobs,jobs,3) + &                
                    fluxc*coef*gp*sin2xi*bper**(alpha+1.d0)*dv

               !FLUXUE
               map(iobs,jobs,4)= map(iobs,jobs,4) + &                
                    fluxc*coef*gp*sin2xipi*bper**(alpha+1.d0)*dv                 
               !FLUXI
               map(iobs,jobs,5)= map(iobs,jobs,5) + &                
                    fluxc*coef*fp*bper**(alpha+1.d0)*dv
               !if(map(iobs,jobs,5).lt.1e-1) map(iobs,jobs,5)=0.
               
               !EM
               map(iobs,jobs,6)= map(iobs,jobs,6) + &                
                   (r2/rhosc)*(r2/rhosc)*dv
               !Rx
               map(iobs,jobs,7)= map(iobs,jobs,7) + &                
                   (r2/rhosc)*(r2/rhosc)*crx*dv
            end if
!          end if
      end do
    end do
  end do


  
! map(:,:,1:5) are: fluxQ, fluxQE, FluxU, FluxUE, FluxI, FluxRX, posangle

  map(:,:,1)= map(:,:,1)/dxT/dyT!/maxval(map(:,:,1))/(dxT*dyT)
  map(:,:,2)= map(:,:,2)/dxT/dyT!/maxval(map(:,:,2))/(dxT*dyT)
  map(:,:,3)= map(:,:,3)/dxT/dyT!/maxval(map(:,:,3))/(dxT*dyT)
  map(:,:,4)= map(:,:,4)/dxT/dyT!/maxval(map(:,:,4))/(dxT*dyT)
  map(:,:,5)= map(:,:,5)/dxT/dyT!/maxval(map(:,:,5))/(dxT*dyT)!/fmax
  map(:,:,6)= map(:,:,6)/dxT/dyT!/maxval(map(:,:,5))/(dxT*dyT)!/fmax
  map(:,:,7)= map(:,:,7)*rsc/dxT/dyT!/maxval(map(:,:,5))/(dxT*dyT)!/fmax

        
!!  do i = 1, nxmap
!!     xcenter= float(i-nxmap/2)
!!     do j = 1, nymap
!!        ycenter = float(j-nymap/2)
!!        rdis= sqrt(xcenter**2+ycenter**2)
!!
!!        if(rdis.gt.70) then
!!           if(abs(map(i,j,1)).lt.1.e-5) map(i,j,1)=0.
!!           if(abs(map(i,j,2)).lt.1.e-5) map(i,j,2)=0.
!!           if(abs(map(i,j,3)).lt.1.e-5) map(i,j,3)=0.
!!           if(abs(map(i,j,4)).lt.1.e-5) map(i,j,4)=0.
!!           if(abs(map(i,j,5)).lt.1.e-5) map(i,j,5)=0.
!!           if(abs(map(i,j,6)).lt.1.e-5) map(i,j,6)=0.
!!           if(abs(map(i,j,7)).lt.1.e-5) map(i,j,7)=0.
!!           if(abs(map(i,j,8)).lt.1.e-5) map(i,j,8)=0.
!!        end if
!! 
!!!        if(rdis.gt.75) then
!!!           map(i,j,1)=0.
!!!           map(i,j,2)=0.
!!!           map(i,j,3)=0.
!!!           map(i,j,4)=0.
!!!           map(i,j,5)=0.
!!!           map(i,j,6)=0.
!!!           map(i,j,7)=0.
!!!           map(i,j,8)=0.
!!!        end if
!!        
!!     end do
!!  end do
end subroutine fill_map

!=======================================================================

!> @brief Writes projection to file
!> @details Writes projection to file
!> @param integer [in] itprint : number of output
!> @param string [in] filepath : path where to write
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param integer [in] nvmap : Number of velocity channels
!> @param real [in] map(nxmap,mymap) : Target map

subroutine  write_polfrac(itprint,filepath,nxmap,nymap,nvmap,map)
  use constants, only : pi
  use ieee_arithmetic
  implicit none
  integer, intent(in) :: nxmap, nymap,nvmap,itprint
  character (len=128), intent(in) :: filepath
  real, intent(inout) :: map(nxmap,nymap,nvmap)
  character (len=128) file1, file2
  integer ::  unitout, unitout2
  integer :: i, j, k
  real :: farrot
  real :: xcenter, ycenter, sinr, cosr, bnor, beast, chiref(nxmap,nymap)
  real :: chi, rdis

  write(file1,'(a,i3.3,a)')  '/storage2/pablov/guacho/arturo/sinx50y0z0hybpol_',itprint,'.bin'  
!  write(file1,'(a,i3.3,a)')  '/storage2/pablov/guacho/aroche/outpart/sinx0y90z0hybpol_',itprint,'.bin'
!  write(file1,'(a,i3.3,a)')  trim(filepath)//'BIN/sinx90y0z0hybtot_',itprint,'.bin'
!  write(file1,'(a,i3.3,a)')  trim(filepath)//'BIN/sinperx0y30z0noRM_',itprint,'.bin'
  unitout=11
  unitout2=12

  write(file2,'(a,i3.3,a)')  '/storage2/pablov/guacho/arturo/sincro_alfahigh_',itprint,'.dat'
!  write(file2,'(a,i3.3,a)')  trim(filepath)//'BIN/sincro_alfahigh_',itprint,'.dat'
  
  open(unit=unitout,file=file1,status='unknown',access='stream', &
       convert='LITTLE_ENDIAN')
  open(unit=unitout2,file=file2,status='unknown')
  farrot = 0.
  
! CHECK IF THERE ARE NANS
  do i = 1, nxmap
     do j = 1, nymap
        do k = 1,nvmap
          
          if ( ieee_is_nan(map(i,j,k)) ) then
             print*,map(i,j,k)
          end if

        end do
     end do
  end do

!!$  !  COMPUTE REFERENCE ANGLE
!!$  map(:,:,7)=90.*atan2(map(:,:,4),map(:,:,2))/pi+(90.-farrot)
!!$  do i = 1, nxmap
!!$    xcenter= real(i)-real(nxmap/2)
!!$    do j = 1, nymap
!!$      ycenter = real(j)-real(nymap/2)
!!$      !distance from centred of grid (plane of the sky)
!!$      rdis= sqrt(xcenter**2+ycenter**2)
!!$      sinr=sin(atan2(-xcenter,ycenter))
!!$      cosr=cos(atan2(-xcenter,ycenter))
!!$      chi          = 0.5*atan2(map(i,j,4),map(i,j,2))
!!$      bnor         = cos(chi+(90.-farrot)*pi/180.)
!!$      beast        = sin(chi+(90.-farrot)*pi/180.)
!!$      chiref(i,j)  = acos((cosr*bnor+sinr*beast))
!!$      map(i,j,8)  = chiref(i,j)*180./pi
!!$      map(i,j,9)  = atan2(-xcenter,ycenter)*180./pi
!!$    end do
!!$  end do

  write (unitout) map(:,:,:)
!  print*, map(:,:,1)
!  write (unitout2,*) sngl(map(:,:,:))
!  print* map(:,100,9)

  close(unitout)
  print'(a,a)'," wrote file:",trim(file1)
  
end subroutine write_polfrac

!!!!=======================================================================
!!!
!!!!> @brief This routine computes a gaussian line profile
!!!!> @details This routine computes a gaussian line profile
!!!
!!!subroutine phigauss(T,vzn,vmin,vmax,nvmap,profile) 
!!!
!!!  use constants, only: amh, pi, kB, clight
!!!  implicit none
!!!  real, intent(in) :: T, vzn, vmin, vmax
!!!  integer, intent(in) :: nvmap
!!!  real, intent(out) :: profile(nvmap)
!!!  integer :: i
!!!  real :: coef, dv, vr
!!!  
!!!  profile(:)=0.
!!!  dv=(vmax-vmin)/float(nvmap)
!!!  
!!!  coef=amh/(2.*kB*T)
!!!  
!!!  do i=1,nvmap
!!!     vr=(float(i)-0.5)*dv+vmin
!!!     profile(i)=sqrt(coef/pi)*exp(-coef*((vr-vzn)**2) )
!!!  end do
!!!  
!!!end subroutine phigauss
!!!
!!!!=======================================================================

end module polfracrad_utilities

!=======================================================================

!> @brief Computes the Ly-alpha apbsorption
!> @details Computes the Ly-alpha apbsorption
!! @n It rotates the data along each of the coordinates axis
!! by an amount @f$ \theta_x, \theta_y, \theta_z @f$, and the LOS
!! is along the Z axis

program polfracrad

  use constants, only : pi
  use parameters, only : xmax,master, mpi_real_kind
  use globals, only : u, up, dx, dy, rank, comm3d
  use polfracrad_utilities

  implicit none
#ifdef MPIP
  include "mpif.h"
#endif
  character (len=128) :: filepath
  integer :: err
  integer :: itprint
  ! 
  real, parameter :: theta_x = 50.1*pi/180.
  real, parameter :: theta_y = 0.1*pi/180.
  real, parameter :: theta_z = 0.1*pi/180.
  !   map and its dimensions
  integer, parameter :: nxmap=300, nymap=300, nzmap=300, nvmap=7
  real :: dxT, dyT, vmin,vmax
  real :: map(nxmap, nymap,nvmap), map1(nxmap, nymap,nvmap)
  real :: dbaver,sumbx,sumbx2,sumby,sumby2,ncontbtot,baver
  real :: dbaver2,sumb_2x,sumb_2x2,sumb_2y,sumb_2y2,baver2
  real :: sumtotx, sum2totx,sumtoty, sum2toty,varx,vary,varx2,vary2
!  real :: theta_x,theta_z,ran(2)
  integer:: ncontb,ncontb_2
  

  ! initializes program
  call init_polfrac()

  !  minumim and maximum of the velocity channels
  !  Target pixel size, relative to the simulation
  dxT= xmax/float(nxmap)
  dyT= dxT
  ! dzT= dyT
!  call init_azar()
!  call random_number(ran)
!  if (ran(1).le.0.5) ran(1)=1.-ran(1)
  !theta_x=ran(1)*90.*pi/180.


!  call init_azar()
!  call random_number(ran)
  !theta_z=ran(2)*90.*pi/180.
  
  !print*,'theta_x, ran: ',theta_x*180./pi,ran(1)
  !print*,'theta_z,ran:',theta_z*180./pi,ran(2)

  ! chose output (fix later to input form screen)
!  filepath='../aroche/outBnoise0.5Cour0.2/'
!  filepath='../aroche/outBnoise0.9Cour0.2/'
!  filepath='../aroche/outsinruido/'
!  filepath = '../aroche/outnoise0.5bgrad/'
!  filepath='../aroche/outsinruidoB1muGCour0.2/'
!  filepath='../aroche/outsinruidoCour0.2eta0.005/'
!  filepath='/storage2/esquivel/guacho-working/snr-pic/'
!  filepath='/storage2/aroche/guacho/snr/output30hyb_080319/'
!  filepath='/storage2/aroche/guacho/snr/output70hyb_110319/'
!   filepath ='../snr/output70hyb_110319/'
!  filepath ='/storage2/aroche/guacho/snr/output05hyb_240419/'
!  filepath ='/storage2/aroche/guacho/snr/output15hyb_240419/'
!  filepath='../aroche/outg1.3eta0.002rad25/'
!  filepath = '../aroche/outnoise0.5bgrad/'
!  filepath='../aroche/outsinruidoB1muGCour0.2/'
!  filepath='../aroche/outsinruidoCour0.2eta0.005/'
!  filepath='/storage2/esquivel/guacho-working/snr-pic/'
!  filepath='/storage2/aroche/guacho/snr/outputeta0.01/'
!  filepath='/storage2/esquivel/guacho-working/turb_snr-M10/' 
  filepath='/storage2/acruz/guacho/SGR19/output/'
  
  loop_over_outputs : do itprint=20,20
    
    !  read ph and u from file
    call read_data(u,itprint,filepath)
    call read_data(up,itprint-1,filepath)
    !
    !  resets map
    map(:,:,:)=0.
    map1(:,:,:)=0.
    sumtotx=0.
    sum2totx=0.
    sumtoty=0.
    sum2toty=0.
    ncontbtot=0.
    !
  !
    if (rank == master) then
       print'(a)', 'Calculating projection with angles of rotaton'
       print'(f6.2,a,f6.2,a,f6.2,a)',theta_x*180./pi,'° around X, ' &
                                    ,theta_y*180./pi,'° around Y, '&
                                    ,theta_z*180./pi,'° around Z, '
    end if
    ! obtain sum and sum^2
    call statb(u,up,sumbx,sumbx2,sumby,sumby2,ncontb,&
         sumb_2x,sumb_2x2,sumb_2y,sumb_2y2,ncontb_2,theta_x,theta_y,theta_z)
    !
    !print*,'sum: ',sumb,sumb2,ncontb*1.
    ! calculate sigma
    call mpi_allreduce(sumbx,sumtotx,1, mpi_real_kind, mpi_sum,comm3d, err)
    !
    call mpi_allreduce(sumbx2,sum2totx,1, mpi_real_kind, mpi_sum,comm3d, err)
    !
    call mpi_allreduce(sumby,sumtoty,1, mpi_real_kind, mpi_sum,comm3d, err)
    !
    call mpi_allreduce(sumby2,sum2toty,1, mpi_real_kind, mpi_sum,comm3d, err)
    !    
    call mpi_allreduce(ncontb*1.,ncontbtot,1, mpi_real_kind, mpi_sum,comm3d, err)
    !
    ! print*,'ncontbtot: ',ncontbtot, sumtot,sum2tot
!    stop
    if(ncontbtot.gt.1.)then
       baver=sqrt((sumtotx/ncontbtot)**2.+(sumtoty/ncontbtot)**2.)
       varx=(sum2totx-(sumtotx**2./ncontbtot))/(ncontbtot-1.)
       vary=(sum2toty-(sumtoty**2./ncontbtot))/(ncontbtot-1.)
       dbaver=sqrt(varx+vary)
       print*,'Bmean1;  dbaver1, eta1 = ',baver, dbaver, (baver/dbaver)**2.
!       stop
    else
       dbaver=1.d-10
       print*,'entre aqui'
    end if

    ! calculo el ruido del lado chocado
    !
    sumtotx=0.
    sum2totx=0.
    sumtoty=0.
    sum2toty=0.
    ncontbtot=0.
    !
    call mpi_allreduce(sumb_2x,sumtotx,1, mpi_real_kind, mpi_sum,comm3d, err)
    !
    call mpi_allreduce(sumb_2x2,sum2totx,1, mpi_real_kind, mpi_sum,comm3d, err)
    !
    call mpi_allreduce(sumb_2y,sumtoty,1, mpi_real_kind, mpi_sum,comm3d, err)
    !
    call mpi_allreduce(sumb_2y2,sum2toty,1, mpi_real_kind, mpi_sum,comm3d, err)
    !    
    call mpi_allreduce(ncontb_2*1.,ncontbtot,1, mpi_real_kind, mpi_sum,comm3d, err)
    !
    ! print*,'ncontbtot: ',ncontbtot, sumtot,sum2tot
!    stop
    if(ncontbtot.gt.1.)then
       baver2=sqrt((sumtotx/ncontbtot)**2.+(sumtoty/ncontbtot)**2.)
       varx2=(sum2totx-(sumtotx**2./ncontbtot))/(ncontbtot-1.)
       vary2=(sum2toty-(sumtoty**2./ncontbtot))/(ncontbtot-1.)
       dbaver2=sqrt(varx2+vary2)
       print*,'Bmean2;  dbaver2, eta2 = ',baver2,dbaver2, (baver2/dbaver2)**2. 
!       stop
    else
       dbaver=1.d-10
       print*,'entre aqui'
    end if

    !  add info to the map
    call fill_map(nxmap,nymap,nvmap,vmin,vmax,u,up,map,dxT,dyT, theta_x, theta_y, theta_z,baver,dbaver,baver2,dbaver2)
    !  sum all the partial sums
    call mpi_reduce(map,map1,nxmap*nymap*nvmap, mpi_real_kind, mpi_sum, master, comm3d, err)
    
    !  write result
    if (rank == master) then 
      call write_polfrac(itprint,filepath,nxmap,nymap,nvmap,map1)
    end if
    
  end do loop_over_outputs

  if (rank == master) print*, 'my work here is done, have a  nice day'
#ifdef MPIP
  call mpi_finalize(err)
#endif
  !
  stop

end program polfracrad

!=======================================================================





