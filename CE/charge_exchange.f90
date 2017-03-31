!=======================================================================
!> @file charge_exchange.f90
!> @brief Implicit & explicit scheme for hydrogen charge exchange
!> @author Krapp Leonardo
!> @date 2/Mar/2016

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

module charge_exchange


  implicit none
contains


!!!-------------------------------------

subroutine exchange(dt)

  use parameters
  use globals
  use constants
  use exoplanet
  implicit none

  real :: dt, xc_n, xc_i, xh_n, xh_i, b
  real :: beta, s1, s2    !rate coefficient cm^3/s
  real :: u_n     !keep the xhi's n-step value
  integer :: i, j, k

  ! beta*dt has units of cm^3 - u has units of cm^{-3} -> b is adimentional

  if(active.eq.1) then
    beta = 4.0E-08
  else
    beta = 1.0E-18
  endif

  do i=1,nx
    do j=1,ny
      do k=1,nz

        b = beta*dt*u(1,i,j,k)

        xc_i = u(neqdyn+2,i,j,k)/u(1,i,j,k)
        xc_n = u(neqdyn+3,i,j,k)/u(1,i,j,k)
        xh_i = u(neqdyn+4,i,j,k)/u(1,i,j,k)
        xh_n = u(neqdyn+5,i,j,k)/u(1,i,j,k)

        u_n = xh_i !! x hot (star)

        ! update xhi

        !!Implicit integration
        xh_i =  (u_n +  b*(u_n+xh_n)*(u_n+xc_i))/(1.0+b)

        !!Explicit integration
        !xh_i = u_n + b*(xh_n*xc_i - u_n*xc_n)

        !update xcn
        xc_n =  xc_n +  (xh_i-u_n)
        !update xhn
        xh_n =  xh_n -  (xh_i-u_n)

        !update xc_i
        xc_i = xc_i - (xh_i-u_n)


        !update fractions
        u(neqdyn+2,i,j,k) = xc_i*u(1,i,j,k)
        u(neqdyn+3,i,j,k) = xc_n*u(1,i,j,k)
        u(neqdyn+4,i,j,k) = xh_i*u(1,i,j,k)
        u(neqdyn+5,i,j,k) = xh_n*u(1,i,j,k)

        !update density of neutrals
        u(neqdyn+1,i,j,k) = (xc_n+xh_n)*u(1,i,j,k)

      end do
    end do
  end do


end subroutine exchange

!=======================================================================

end module charge_exchange
