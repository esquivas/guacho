!=======================================================================
!> @file user_mod.f90
!> @brief User input module
!> @author C. Villarreal, M. Schneiter, A. Esquivel
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
  use jet 

  implicit none
 
contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty 
subroutine init_user_mod()

  implicit none      
  !  initialize modules loaded by user
  call init_jet()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, Tempsc
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, parameter :: T_ism = 100./Tempsc, dens_ism= 100.
  integer :: i,j,k
 
  do k=nzmin,nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax

        !   total density and momenta
        u(1,i,j,k) = dens_ism
        u(2,i,j,k) = 0.
        u(3,i,j,k) = 0.
        u(4,i,j,k) = 0.

        !  passive scalars needed for the chemistry network
        !  molecular environmet
        u( 6,i,j,k) = 0.0010   * dens_ism  ! nH0
        u( 7,i,j,k) = 0.0001   * dens_ism  ! nH+
        u( 8,i,j,k) = 0.9989/2.* dens_ism  ! nH2
        u( 9,i,j,k) = 0.0001   * dens_ism  ! ne (=nH+)
        u(10,i,j,k) = - u(1,i,j,k)

        !   energy
        u(5,i,j,k)=cv*0.50065*dens_ism*T_ism

      end do        
    end do
  end do

  call impose_jet(u,0.)



end subroutine initial_conditions
  
!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)
!> @param integer [in] order : order (mum of cells to be filled in case
!> domain boundaries are being set)

subroutine impose_user_bc(u,order)

  use parameters, only:  neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals,    only: time
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: order

  !  In this case the boundary is the same for 1st and second order)
  if (order >= 1) then 
    call impose_jet(u,time)
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
  !use constants,  only : Ggrav
  !use parameters, only : nx, ny, nz, nxtot, nytot, nztot, rsc, vsc2
  !use globals,    only : dx, dy, dz, coords
  !use jet
  implicit none
  integer, intent(in) :: i,j,k
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)
  !real :: xc,yc,zc, GM, rad
 
  !GM=Ggrav*MassS/rsc/vsc2

  !   get cell position
  !xc=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx
  !yc=(real(j+coords(1)*ny-nytot/2)+0.5)*dy
  !zc=(real(k+coords(2)*nz-nztot/2)+0.5)*dz

  ! calculate distance from the sources
  !rad = xc**2+yc**2+zc**2
  ! update source terms with gravity
  ! momenta
  !s(2)= s(2)-pp(1)*GM*xc/(rad**1.5)
  !s(3)= s(3)-pp(1)*GM*yc/(rad**1.5)
  !s(4)= s(4)-pp(1)*GM*zc/(rad**1.5)
  ! energy
  !s(5)= s(5)-pp(1)*GM*( pp(2)*xc +pp(3)*yc +pp(4)*zc )  &
  !       /(rad**1.5 )

end subroutine get_user_source_terms


!=======================================================================

end module user_mod

!=======================================================================
