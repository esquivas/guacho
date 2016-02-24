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

  ! load auxiliary modules
  use jet
  implicit none
 
contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty 
subroutine init_user_mod()

  implicit none      
  !  if needed initialize modules loaded by user
  call init_jet()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax,Tempsc, cv, rsc
  use globals, only : coords, dx, dy, dz, rank, time
  use constants, only : pc
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, parameter :: T_ism = 10. 


  integer :: i,j,k
  
  !  ISM parameters
  do k=nzmin,nzmax
    do i=nxmin,nxmax
      do j=nymin,nymax

        !  dynamical equations
        u(1,i,j,k) = 1.e3
        u(2,i,j,k) = 0.
        u(3,i,j,k) = 0.
        u(4,i,j,k) = 0.
        u(5,i,j,k) = cv*u(1,i,j,k)*1.001*T_ism/Tempsc

        !  passive scalars needed for the chemistry network
        u( 6,i,j,k) = 0.0001    *  u(1,i,j,k)
        u( 7,i,j,k) = 0.9989/2. *  u(1,i,j,k)
        u( 8,i,j,k) = 0.001     *  u(1,i,j,k)
        u( 9,i,j,k) = 0.0001    *  u(1,i,j,k)
        u(10,i,j,k) = - u(1,i,j,k)

      end do
    end do
  end do

  !  place SN at the center of the grid
  call impose_jet(u, 0.)


end subroutine initial_conditions
  
!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

#ifdef OTHERB
subroutine impose_user_bc(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals   , only : time 
  implicit none
  real :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

  call impose_jet(u, time )
 
end subroutine impose_user_bc

!=======================================================================

#endif

end module user_mod

!=======================================================================
