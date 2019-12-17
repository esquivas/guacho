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
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, .true.            &
       , comm3d, err)
  call mpi_comm_rank(comm3d, rank, err)
  call mpi_cart_coords(comm3d, rank, ndim, coords, err)
  print '(a,i3,a,3i4)', 'processor ', rank                                  &
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

    x=(float(i+coords(0)*nx-nxtot/2)-0.5)*dx
    y=(float(j+coords(1)*ny-nytot/2)-0.5)*dy
    z=(float(k+coords(2)*nz-nztot/2)-0.5)*dz

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


subroutine fill_map(nxmap,nymap,nvmap,vmin,vmax,u,up,map,dxT,dyT,&
                   theta_x,theta_y,theta_z)

  use constants, only : clight, pi
  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         neq, nx, ny, nz, vsc2, rsc, rhosc,nztot, neqdyn, &
                         bsc, mhd
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
  real :: x,y,z,xn,yn,zn, vx, vy, vz,vxn, vyn, vzn, velsc, pres1, pres2,presc
  real :: T, prim(neq), profile(nvmap)
  real, parameter :: sigmaLA = 0.01105, lambdaLA=1.215668e-5 !(c/nu0=lambda)
  real :: alpha, f0, rp, norx, nory, norz, bt, vt, dv, bper, alphaE, valf
  real :: norbx, norby, norbz, coef, cos2xi, sin2xi, cos2xipi, sin2xipi
  real :: bx, by, bz, bxn, byn, bzn
  real :: bxx, byy, bzz, exx, eyy, fluxc,r, r1, r2, gg
  real :: cos2_B, sin2_B, fmax
  real :: eta,crx, xp, aa, bb, crho, coefp, deltaS
  real :: phirx(201)
  real :: bfn(7)
  integer :: qp, ip, np, ip1, itrx
  real :: xcenter, ycenter, bnor, beast, rdis, chi, sinr, cosr
  real :: farrot, cut_rho, fp, gp
  real :: p_ip1,p_im1,dpx,dpy,dpz,dptot,txr,baver,dbaver,fef
  real :: chiref(nxmap,nymap)
  character (len=128) :: filein  = './rx/coef0.1_2.4kev.dat'! Input File

  velsc=sqrt(vsc2)
  presc=rhosc*vsc2

  !some constants
  itrx    = 0
!
  do k=1,nz
     do j=1,ny
        do i=1,nx

          !  obtain original position
          call getXYZ(i,j,k, x,y,z)

          !distance to the centre of the grid
          rp =sqrt( x**2+ y**2+ z**2)

!!$          ! normalized position w.r.t c. of grid
!!$          norx = x/rp
!!$          nory = y/rp
!!$          norz = z/rp

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
          pres2=prim(5)*presc
          r2 = prim(1)*rhosc
          if(mhd) then
            bx=prim(6)*bsc ! Should we scale?
            by=prim(7)*bsc
            bz=prim(8)*bsc
            bt = sqrt(bx**2+by**2+bz**2) ! DEFINIR
            !            print*, 'Temperature',T
          endif
          txr=T
          vt = sqrt(vx**2+vy**2+vz**2) ! DEFINIR
          !
          if ( ieee_is_nan(fluxc) ) then
             write(*,*) 'map is NaN', fluxc, vt, r, velsc, vx, vy, vz
             stop
          endif
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

          !  calculate the line profile function
          !            call phigauss(T, vzn,vmin,vmax,nvmap,profile)
          !  make sure the result lies in the map bounds
          if( (iobs >=1    ).and.(jobs >=1    ).and. &
               (iobs <=nxmap).and.(jobs <=nymap) ) then
             !EM
             map(iobs,jobs,1)= map(iobs,jobs,1) + &
                  (r2/rhosc)*(r2/rhosc)*dv
             !Rx
             map(iobs,jobs,2)= map(iobs,jobs,2) + &
                  (r2/rhosc)*(r2/rhosc)*crx
          end if
          !          end if
       end do
    end do
 end do



! map(:,:,1:5) are: fluxQ, fluxQE, FluxU, FluxUE, FluxI, FluxRX, posangle


  map(:,:,1)= map(:,:,1)/dxT/dyT!/maxval(map(:,:,5))/(dxT*dyT)!/fmax
  map(:,:,2)= map(:,:,2)*rsc/dxT/dyT!/maxval(map(:,:,5))/(dxT*dyT)!/fmax



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
  character (len=128) file1
  integer ::  unitout
  integer :: i, j, k
  real :: farrot
  real :: xcenter, ycenter, sinr, cosr, bnor, beast, chiref(nxmap,nymap)
  real :: chi, rdis

  write(file1,'(a,i3.3,a)')  trim(filepath)//'BIN/xray_',itprint,'.bin'
  unitout=11


  open(unit=unitout,file=file1,status='unknown',access='stream', &
       convert='LITTLE_ENDIAN')
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
#ifdef MPIP
  use mpi
#endif

  implicit none
  character (len=128) :: filepath
  integer :: err
  integer :: itprint
  !
  real, parameter :: theta_x = 0.*pi/180.
  real, parameter :: theta_y = 0.*pi/180.
  real, parameter :: theta_z = 0.*pi/180.
  !   map and its dimensions
  integer, parameter :: nxmap=256, nymap=256, nvmap=2
  real :: dxT, dyT, vmin,vmax
  real :: map(nxmap, nymap,nvmap), map1(nxmap, nymap,nvmap)

  ! initializes program
  call init_polfrac()

  !  minumim and maximum of the velocity channels
  !  Target pixel size, relative to the simulation
  dxT= xmax/float(nxmap)
  dyT= dxT
 ! dzT= dyT

  ! chose output (fix later to input form screen)
!  filepath='/datos/esquivel/EXO-GUACHO/P1c/'
!  filepath='../alicia/outputRuidoby/'
!  filepath='../tycho/output_run_ek01_vbpar/'
!    filepath='../ran/output/M1/'
  filepath='./'

  loop_over_outputs : do itprint=45,50

    !  read ph and u from file
    call read_data(u,itprint,filepath)
!    call read_data(up,itprint-4,filepath)
    !
    !  resets map
    map(:,:,:)=0.
    map1(:,:,:)=0.
    !
    if (rank == master) then
       print'(a)', 'Calculating projection with angles of rotaton'
       print'(f6.2,a,f6.2,a,f6.2,a)',theta_x*180./pi,'° around X, ' &
                                    ,theta_y*180./pi,'° around Y, '&
                                    ,theta_z*180./pi,'° around Z, '
    end if

    !  add info to the map
    call fill_map(nxmap,nymap,nvmap,vmin,vmax,u,up,map,dxT,dyT, theta_x, theta_y, theta_z)
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
