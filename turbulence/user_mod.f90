!=======================================================================
!> @file user_mod.f90
!> @brief User input module
!> @author C. Villareal, M. Shneiter, A. Esquivel
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

!> @brief User imput module
!> @details  This is an attempt to have all input neede from user in a
!! single file
!!!  This module should load additional modules (i.e. star, jet, sn), to
!!  impose initial and boundary conditions (such as sources)

module user_mod

  ! load auxiliary modules
  use msnr
  
  implicit none
 
contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty 
subroutine init_user_mod()

  implicit none      
  !  initialize modules loaded by user
  !call init_msnr()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
       cv, rsc, psc, rhosc, vsc, bsc
  
  use constants, only : pi, kb, amh

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  
  ! Defining parameters for SNR 
  real, parameter :: mu_ism =1.3
  real :: de, temp

  ! ENVIRONMENT
  de = 5.d-2           !NUMBER DENSITY
  temp = 1000d0        !TEMPERATURE
  
  u(:,:,:,:)=0.


  !START ENVIRONMENT
  u(1,:,:,:) = de*amh*mu_ism/rhosc
  u(6,:,:,:) = 0.
  u(7,:,:,:) = 2.e-6/bsc
  u(8,:,:,:) = 0.
  u(5,:,:,:) = cv*(de*kb*temp)/psc+ &
       0.5*(u(6,:,:,:)**2+u(7,:,:,:)**2+u(8,:,:,:)**2) 
  
  !   place SN at the center of the grid
  call impose_msnr(u)

 end subroutine initial_conditions
  
!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)
!> @param integer [in] order : order (mum of cells to be filled in case
!> domain boundaries are being set, valid values are 1 or 2)

subroutine impose_user_bc(u,order)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         bc_user, tsc
  use globals   , only : time 
  implicit none
  real :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: order
  real    :: time_sec

  time_sec = time*tsc  

  if (bc_user) then  
    if (order == 1) then 

    else if (order == 2) then
    
    end if
  end if

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

  ! in this example a constant gravity is added
  use constants,  only : Ggrav,Msun,Rsun
  use parameters, only : neq, nymin, nymax, rsc, vsc
  implicit none
  integer, intent(in) :: i,j,k
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)


end subroutine get_user_source_terms

!=====================================================================

end module user_mod

!=======================================================================
