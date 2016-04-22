!=======================================================================
!> @file user_mod.f90
!> @brief User input module
!> @author Alejandro Esquivel
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
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief User imput module
!> @details  This is an attempt to have all input neede from user in a
!! single file
!!!  This module should load additional modules (i.e. star, jet, sn), to
!!  impose initial and boundary conditions (such as sources)

module user_mod

  use exoplanet
  ! load auxiliary modules
  implicit none
 
contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty 
subroutine init_user_mod()

  implicit none      
  !  initialize modules loaded by user
  call init_exo()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals,    only: coords, dx ,dy ,dz

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

  integer :: i,j,k
  real :: x,y,z, rads, velx, vely, velz, dens,cpi
  !  the star wind does not cover the entire domain, we fill here 
  !  as if the exoplanet is absent
  do i=nxmin,nxmax
    do j=nymin,nymax
      do k=nzmin,nzmax

        ! Position measured from the centre of the grid (star)
        x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
        y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
        z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz

        ! Distance from the centre of the star
        rads=sqrt(x**2+y**2+z**2)

        VelX=VSW*X/RADS
        VelY=VSW*Y/RADS
        VelZ=VSW*Z/RADS
        DENS=DSW*RSW**2/RADS**2
        !   total density and momena
        u(1,i,j,k) = dens
        u(2,i,j,k) = dens*velx
        u(3,i,j,k) = dens*vely
        u(4,i,j,k) = dens*velz
#if defined(PMHD) || defined(MHD)
        cpi = bsw*(RSW/rads)**3/(2.*rads**2)
        u(6,i,j,k) = 3.*y*x*cpi
        u(7,i,j,k) = (3.*y**2-rads**2)*cpi
        u(8,i,j,k) = 3.*y*z*cpi

#endif
#ifdef MHD
        ! total energy
        u(5,i,j,k)=0.5*dens*vsw**2         &
             + cv*dens*Tsw       & 
             + 0.5*(u(6,i,j,k)**2+u(7,i,j,k)**2+u(8,i,j,k)**2)
#else
              ! total energy
        u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2) &
             + cv*dens*1.9999*Tsw
#endif
#ifdef PASSIVES
        !  density of neutrals
        u(neqdyn+1,i,j,k)= 0.0001*dens
        !   passive scalar (h-hot, c-cold, i-ionized, n-neutral)
        u(neqdyn+2,i,j,k)= dens   ! passive scalar
#endif
      end do
    end do
  end do

  call impose_exo(u,0.)

end subroutine initial_conditions
  
!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

#ifdef OTHERB
subroutine impose_user_bc(u)

  use parameters, only:  neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals,    only: time

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

  call impose_exo(u,time)
 
end subroutine impose_user_bc

!=======================================================================

#endif

end module user_mod

!=======================================================================
