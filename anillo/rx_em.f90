!=======================================================================
!> @file Rx_em.f90
!> @brief Rx Emission map utilities
!> @author P.  Velázquez M. Schneiter, Alejandro Esquivel
!> @date 4/May/2020

! Copyright (c) 2020 Guacho Co-Op
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

!> @brief Rx_utilities
!> @details Utilities to compute X ray maps

module rx_utilities

  real, allocatable :: phirx_soft(:), phirx_hard(:)  ! Rx emissivity read from table

contains

!> @brief Initializes data
!> @details Initializes data, MPI and other stuff
subroutine init_RX()

!  Initializes MPI, data arrays, etc
use parameters
use globals, only : u, dx, dy, dz, coords, rank, left, right   &
                     , top, bottom, out, in, rank, comm3d
implicit none
  integer :: nps, err
  integer, dimension(0:ndim-1) :: dims
  logical, dimension(0:ndim-1) :: period
  character (len=128) ::  filesoft = 'coef0.2_2kev.dat' ! Input File
  character (len=128) ::  filehard = 'coef2_10kev.dat'   ! Input File
  integer             :: ip, nd, proc
  real                :: aa, bb
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

!  master reads table and shares

!   soft Xrays
if (rank==master) then
  open(unit=10,file=trim(workdir)//'../../src/xrlib/'//trim(filesoft),status='unknown')
  read(10,*) nd
  allocate(phirx_soft(nd))
  do ip=1,nd
    read(10,*) aa, bb, phirx_soft(ip)
  end do
end if
call mpi_bcast(nd,  1,mpi_integer  ,0, comm3d, err)
if (rank /= master) allocate(phirx_soft(nd))
call mpi_bcast(phirx_soft, nd ,mpi_real_kind  ,0, comm3d, err)

! hard x rays
if (rank==master) then
  open(unit=10,file=trim(workdir)//'../../src/xrlib/'//trim(filehard),status='unknown')
  read(10,*) nd
  allocate(phirx_hard(nd))
  do ip=1,nd
    read(10,*) aa, bb, phirx_hard(ip)
  end do
end if
call mpi_bcast(nd,  1,mpi_integer  ,0, comm3d, err)
if (rank /= master) allocate(phirx_hard(nd))
call mpi_bcast(phirx_hard, nd ,mpi_real_kind  ,0, comm3d, err)

end subroutine init_RX

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
  character                       :: byte_read
  integer :: nxp, nyp, nzp, x0p, y0p, z0p, &
             mpi_xp, mpi_yp, mpi_zp,neqp, neqdynp, nghostp
  integer :: unitin, ip, err
    real :: dxp, dyp, dzp, scal(3), cvp
  character (len=128) file1

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
         open(unit=unitin,file=file1,status='unknown', access='stream')

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

         print'(i3,a,a)',rank,' read file:',trim(file1)

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

    x=( real(i+coords(0)*nx-nxtot/2) - 0.5 )*dx
    y=( real(j+coords(1)*ny-nytot/2) - 0.5 )*dy
    z=( real(k+coords(2)*nz-nztot/2) - 0.5 )*dz

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
   real, intent(out) :: xn, yn, zn
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
subroutine fill_map(nxmap,nymap,u,map,dxT,dyT,&
                   theta_x,theta_y,theta_z)

  use constants, only : clight
  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         neq, nx, ny, nz, vsc2, rsc,nztot, neqdyn
  use globals, only : dx, dy, dz
  use hydro_core, only : u2prim

  implicit none

  integer, intent(in) :: nxmap,nymap
  real, intent(in) :: u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax)
  real , intent(in) :: dxT, dyT, theta_x, theta_y, theta_z
  real, intent(out) :: map(nxmap,nymap,2)
  integer :: i,j,k, iobs, jobs, ip, ip1
  real :: x,y,z,xn,yn,zn,dv, xp
  real :: T, prim(neq), crx_soft, crx_hard

  dv = dx*dy*dz

  do k=1,nz
     do j=1,ny
        do i=1,nx

          !  obtain original position
          call getXYZ(i,j,k, x,y,z)
          !  do the rotation of the coordinates
          call rotation_x(theta_x,x,y,z,xn,yn,zn)
          call rotation_y(theta_y,xn,yn,zn,x,y,z)
          call rotation_z(theta_z,x,y,z,xn,yn,zn)
          ! This is the position on the target (centered)
          ! Integration is along Z
          iobs=xn/dxT + nxmap/2
          jobs=yn/dyT + nymap/2

          !  make sure the result lies within the map bounds
          if( (iobs >=1    ).and.(jobs >=1    ).and. &
             (iobs <=nxmap).and.(jobs <=nymap) ) then

            !  get the velocity in cm/s, and T
            call u2prim(u(:,i,j,k),prim,T)

            !  do the interpolation (needs to be made more general,
            !  for now it is assumed that both tables have 201 elements)
            xp=(log10(T)-4.)/0.02+1.000001
            ip=max(1,int(xp))
            ip=min(ip,200)
            ip1=ip+1
            crx_soft=phirx_soft(ip)+(phirx_soft(ip1)-phirx_soft(ip))*(xp-real(ip))
            crx_hard=phirx_hard(ip)+(phirx_hard(ip1)-phirx_hard(ip))*(xp-real(ip))
            if(T.ge.1.e8) then
              crx_soft=phirx_soft(200)*(T/1.e8)**0.5
              crx_hard=phirx_hard(200)*(T/1.e8)**0.5
            end if
            !soft Rx
            map(iobs,jobs,1)= map(iobs,jobs,1) + &
                                  crx_soft*prim(1)**2 *rsc* dz   !dv/dxT/dyT
            ! hard
            map(iobs,jobs,2)= map(iobs,jobs,2) + &
                                  crx_hard*prim(1)**2 *rsc* dz   !dv/dxT/dyT

          end if
        end do
    end do
  end do

end subroutine fill_map

!=======================================================================
!> @brief Writes projection to file
!> @details Writes projection to file
!> @param integer [in] itprint : number of output
!> @param string [in] filepath : path where to write
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param real [in] map(nxmap,mymap) : Target map
subroutine  write_RX(itprint,filepath,nxmap,nymap,map)
  implicit none
  integer, intent(in) :: nxmap,nymap,itprint
  character (len=128), intent(in) :: filepath
  real, intent(in) :: map(nxmap,nymap,2)
  character (len=128) file1
  integer ::  unitout

  write(file1,'(a,i3.3,a)')  trim(filepath)//'BIN/Rx-2-',itprint,'.bin'
  unitout=11
  open(unit=unitout,file=file1,status='unknown',access='stream')

  write (unitout) map(:,:,:)
  close(unitout)

  print'(a,a)'," wrote file:",trim(file1)

end subroutine write_RX


!=======================================================================

end module rx_utilities

!=======================================================================
!> @brief Computes the Ly-alpha apbsorption
!> @details Computes the Ly-alpha apbsorption
!! @n It rotates the data along each of the coordinates axis
!! by an amount @f$ \theta_x, \theta_y, \theta_z @f$, and the LOS
!! is along the Z axis
program Rx_em
  use constants, only : pi
  use parameters, only : xmax,master, mpi_real_kind, &
                         outputpath, nxtot, nytot
  use globals, only : u, rank, comm3d
  use rx_utilities
#ifdef MPIP
  use mpi
#endif
  implicit none
  character (len=128) :: filepath
  integer :: err
  integer :: itprint
  !
  real, parameter :: theta_x = -90.0 *pi/180.
  real, parameter :: theta_y =  0.00 *pi/180.
  real, parameter :: theta_z =  0.00 *pi/180.
  !   map and its dimensions
  integer            :: nxmap, nymap
  real :: dxT, dyT, vmin,vmax
  real, allocatable :: map(:,:,:), map1(:,:,:)

  nxmap = nxtot
  nymap = nytot
  allocate( map(nxmap, nymap, 2 ))
  allocate(map1(nxmap, nymap, 2 ))

  ! initializes program
  call init_RX()

  !  Target pixel size, relative to the simulation
  dxT= xmax/real(nxmap)
  dyT= dxT

  ! chose output (fix later to input form screen)
  filepath=trim(outputpath) !'/datos/esquivel/EXO-GUACHO/P1c/'

  loop_over_outputs : do itprint=1,100

    !  read ph and u from file
    call read_data(u,itprint,filepath)

    !  resets map
    map ( :,:, : )=0.0
    map1(: ,:, : )=0.0

    if (rank == master) then
       print'(a)', 'Calculating projection with angles of rotaton'
       print'(f6.2,a,f6.2,a,f6.2,a)',theta_x*180./pi,'° around X, ' &
                                    ,theta_y*180./pi,'° around Y, ' &
                                    ,theta_z*180./pi,'° around Z, '
    end if

    !  add info to the map
    call fill_map(nxmap,nymap,u,map, dxT,dyT, theta_x, theta_y, theta_z)
    !  sum all the partial sums
    call mpi_reduce(map,map1,nxmap*nymap*2, mpi_real_kind, mpi_sum, master, comm3d, err)

    !  write result
    if (rank == master) then
      call write_RX(itprint,filepath,nxmap,nymap,map1)
    end if

  end do loop_over_outputs

  if (rank == master) print*, 'my work here is done, have a  nice day'
#ifdef MPIP
  call mpi_finalize(err)
#endif

  stop

end program Rx_em

!=======================================================================
