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
!!!  This module should load additional modules (i.e. star, corona, sn), to
!!  impose initial and boundary conditions (such as sources)

module user_mod

  ! load auxiliary modules
  use corona
  implicit none
 
contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty 
subroutine init_user_mod()

  implicit none      
  !  if needed initialize modules loaded by user
  call init_corona()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, cv, rsc, vsc2, rhosc, Tempsc,ny, gamma
#ifdef MHD_BSPLIT
  use globals, only : coords, dx, dy, dz, rank, time, B0
#else
  use globals, only : coords, dx, dy, dz, rank, time
#endif  
  use constants, only : Ggrav, Msun, Rsun, pi, Rg
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) 
  real :: g, H 
  real :: y, rho_y 
  integer :: i, j
  
  g = Ggrav*Msun/Rsun/Rsun
  !H = mu*g/Rg/Tempc
  
!   do j = nymin,nymax
!     y = (float(j+coords(1)*ny) + 0.5)*dy*rsc
!     rho_y = rhoc*exp(-H*y) 
!     u(1,:,j,:) = rho_y/rhosc
!     u(2,:,j,:) = 0.
!     u(3,:,j,:) = 0.
!     u(4,:,j,:) = 0.    
!     u(5,:,j,:)= (cv*rho_y*Rg*Tempc/mu)/Psc     
!   end do
  
!   g = Ggrav*Msun/Rsun/Rsun
!   cc = mu*g/Rg
!   
  do j = nymin,nymax
    y = (float(j+coords(1)*ny) + 0.5)*dy*rsc
	if(y.le.2.E8) then
          P_y = P_0*exp(-c1*y/Temp1)
        else	
          P_y = P_0*exp(-c2*y/Temp1)
        endif	
    u(5,:,j,:) = cv*P_y/Psc 
    u(2,:,j,:) = 0.
    u(3,:,j,:) = 0.
    u(4,:,j,:) = 0.    
    u(1,:,j,:)= (mu*P_y/Rg/Tempc)/rhosc
      
  end do
        

end subroutine initial_conditions
  
!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)
!> @param integer [in] order : order (mum of cells to be filled in case
!> domain boundaries are being set, valid values are 1 or 2)

subroutine impose_user_bc(u,order)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, bc_other
  use globals   , only : time 
  implicit none
  real :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: order

  !  In this case the boundary is the same for 1st and second order)
  if (order >= 1) then 
    call driver(u,time)
  end if

end subroutine impose_user_bc
!=======================================================================

end module user_mod

!=======================================================================
