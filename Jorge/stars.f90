!=======================================================================
!> @file exoplanet.f90
!> @brief Exoplanet problem module
!> @author M. Schneiter, C. Villarreal  D'Angelo, A. Esquivel
!> @date 2/Nov/2014

! Copyright (c) 2014 A. Esquivel, M. Schneiter, C. Villareal D'Angelo
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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief Exoplanet module
!> @details Problem Module for exoplanet

module stars

  use parameters
  implicit none

  integer :: Nstars
  real, allocatable :: xstar(:), ystar(:), zstar(:), vws(:), mdots(:),sstar(:)
  real :: Tsw, rw
contains

  !=======================================================================
  !> @brief Module initialization
  !> @details Here the parameters of the Star are initialized, and scaled
  !! to code units
  subroutine init_stars()

    use parameters, only : workdir, master, rsc, Tempsc
    use globals,    only : rank
    use constants,  only : Msun, yr, pc
#ifdef MPIP
    use mpi
#endif
    implicit none
    integer :: i, err
    real    :: datain(6)

    Tsw = 1E4/Tempsc
    rw = 0.7*pc/rsc

    !  Master reads data
    if(rank == master) then
      open(unit=10,file= trim(workdir)//'stars.dat',status='old')
      read(10,*) Nstars
      allocate( xstar(Nstars) )
      allocate( ystar(Nstars) )
      allocate( zstar(Nstars) )
      allocate(   vws(Nstars) )
      allocate( mdots(Nstars) )
      allocate( sstar(Nstars) )
  !    print*, 'Nstars:', Nstars
!Reads position (x, y, z), wind velocities,
! mass loss rate and photon emission from a table.
      do i=1,Nstars
        read(10,*) datain(1:6)
        xstar(i) = datain(1) *pc / rsc
        zstar(i) = datain(3) *pc / rsc
        ystar(i) = datain(2) *pc / rsc
        vws(i)   = datain(4) *1.0e5
        mdots(i) = 10**datain(5)*Msun/yr
        sstar(i) = 10**datain(6)
      end do
      close(unit=10)
    endif

    !  master distributes data
#ifdef MPIP
    call mpi_bcast(Nstars,1,mpi_integer,master,mpi_comm_world,err)
    if (rank /= master) then
      allocate( xstar(Nstars) )
      allocate( ystar(Nstars) )
      allocate( zstar(Nstars) )
      allocate(   vws(Nstars) )
      allocate( mdots(Nstars) )
      allocate( sstar(Nstars) )
    end if

    call mpi_bcast(zstar,Nstars,mpi_real_kind,master,mpi_comm_world,err)
    call mpi_bcast(xstar,Nstars,mpi_real_kind,master,mpi_comm_world,err)
    call mpi_bcast(ystar,Nstars,mpi_real_kind,master,mpi_comm_world,err)
    call mpi_bcast(  vws,Nstars,mpi_real_kind,master,mpi_comm_world,err)
    call mpi_bcast(mdots,Nstars,mpi_real_kind,master,mpi_comm_world,err)
    call mpi_bcast(sstar,Nstars,mpi_real_kind,master,mpi_comm_world,err)
#endif

  end subroutine init_stars

  !=======================================================================
  !> @brief Inject sources of wind
  !> @details Imposes the sources of wond from the star and planet
  !> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
  !! conserved variables
  !> @param real [time] time : current integration timr
  !--------------------------------------------------------------------
  subroutine impose_stars(u,time)

    use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax,&
                           rhosc, vsc, rsc
    use globals, only : coords, dx, dy, dz
    use constants, only : pi
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent (in) :: time
    real :: x,y,z,rad,densw,dens,velx,vely,velz
    integer :: i,j,k,l

    do k = nzmin,nzmax
      do j = nymin,nymax
        do i = nxmin,nxmax
          ! Position measured from the centre of the grid (star)
          x=(real(i+coords(0)*nx-nxtot/2)-0.5)*dx
          y=(real(j+coords(1)*ny-nytot/2)-0.5)*dy
          z=(real(k+coords(2)*nz-nztot/2)-0.5)*dz

          do l=1,Nstars

            rad=sqrt((x-xstar(l))**2+(y-ystar(l))**2+(z-zstar(l))**2)
            if (rad <= rw) then
              if(rad == 0.) rad=dx*0.10
              densw=((mdots(l)/rw)/(4*pi*rw*vws(l)))   ! stellar wind density
              densw=densw/(rhosc*rsc**2)

              velx = vws(l) * ( x - xstar(l) )/rad/vsc
              vely = vws(l) * ( y - ystar(l) )/rad/vsc
              velz = vws(l) * ( z - zstar(l) )/rad/vsc

              dens=densw!*rw**2/rad**2

              !   total density and momenta
              u(1,i,j,k) = dens
              u(2,i,j,k) = dens*velx
              u(3,i,j,k) = dens*vely
              u(4,i,j,k) = dens*velz
              ! total energy
              u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2) &
                        + cv*dens*1.9999*Tsw

              u(6,i,j,k) = 1.0E-4*dens
              u(7,i,j,k) = dens
            end if
          end do

        end do
      end do
    end do

  end subroutine impose_stars

  !=======================================================================

end module stars
