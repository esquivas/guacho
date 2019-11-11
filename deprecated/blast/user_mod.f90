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
  use snr
  implicit none

contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty
subroutine init_user_mod()

  implicit none
  !  if needed initialize modules loaded by user

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax,Tempsc, cv, rsc, rhosc, vsc
  use globals, only : coords, dx, dy, dz, rank, time
  use constants, only : pc, amh, Kb
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, parameter :: T_ism = 1000., n_ism = 1.


  integer :: i,j,k

  !  ISM parameters
  do k=nzmin,nzmax
    do i=nxmin,nxmax
      do j=nymin,nymax
        !  conserved variables, scaled to code units, everithing in cgs, then divided by the corresponding scalings
        u(1,i,j,k) = n_ism*amh/rhosc  ! Density in g cm^{-3}/rhosc
        u(2,i,j,k) = 0.
        u(3,i,j,k) = 0.
        u(4,i,j,k) = 0.
        u(5,i,j,k) = cv*n_ism*T_ism/Tempsc
        ! this is the same as cv* nKT / Psc, to this we should
        ! the kinetic energy ek=rho*v^2/2, but since v=0 -> eK=0
      end do
    end do
  end do

  !  place SN at the center of the grid (15pc/rsc = 0.5)
  !  this is in the module snr in snr.f90
  call impose_snr(u, 0.5, 0.5, 0.5 )


end subroutine initial_conditions

!=====================================================================
!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)
!> @param integer [in] order : order (use to distinguish the 1st
!>  order half timestep, from the 2nd order full timestep)

subroutine impose_user_bc(u,order)

  use parameters, only:  neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals,    only: time
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: order

  !  In this case the boundary is the same for 1st and second order)
  !if (order >= 1) then
  !  here we can call an external module to impose e.g. a jet,
  !  a  wind or something like that
  !end if

end subroutine impose_user_bc


!=======================================================================

!> @brief User Defined source terms
!> This is a generic interrface to add a source term S in the equation
!> of the form:  dU/dt+dF/dx+dG/dy+dH/dz=S
!> @param real [in] pp(neq) : vector of primitive variables
!> @param real [inout] s(neq) : vector with source terms, has to add to
!>  whatever is there, as other modules can add their own sources
!> @param integer [in] i : cell index in the X direction
!> @param integer [in] j : cell index in the Y direction
!> @param integer [in] k : cell index in the Z direction

subroutine get_user_source_terms(pp,s, i, j , k)
  use parameters, only : neq
  implicit none
  integer, intent(in) :: i, j, k
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)

end subroutine get_user_source_terms



end module user_mod

!=======================================================================
