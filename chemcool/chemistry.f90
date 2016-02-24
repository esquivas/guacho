!=======================================================================
!> @file chemistry.f90
!> @brief chemistry  module
!> @author A. Castellanos, A. Rodriguez, A. Raga  and A. Esquivel
!> @date 24/Nov/2014

! Copyright (c) 2014 A. Esquivel et al.
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

!> @brief chenistry module
!> @details module to solve the chemical/ionic network, and estimate
!> the cooling assiciated with it.

module chemistry

  use network
  implicit none

contains

!=======================================================================

subroutine chemstep(y,y0,T,deltt)
  !use coolingchem_module
  use linear_system
  use network, only : n_spec, n_reac, n_elem, get_reaction_coefficients,  &
                      derv, get_jacobian, n_nequ
  implicit none
  real (kind=8) :: dt,dtm
  real (kind=8) :: y1(n_spec),yt(n_spec),yin(n_spec),yp(n_spec),y(n_spec)
  real (kind=8) :: rate(n_reac),dydt(n_spec),jac(n_spec,n_spec),y0(n_elem)
  real (kind=8) :: T, deltt
  integer, parameter  :: niter=50        ! number of iterations
  integer :: n,iflag,i,iff,j
  
  
  n=0
  iflag=1
  dtm=1./deltt
  iff=1
  yin(1:n_spec) =y (1:n_spec)
  yp (1:n_spec) =y (1:n_spec)

  call get_reaction_coefficients(rate,T)

  !  call totalelem
  call inittest2(yp,y)

  do while (n <= niter. and. iflag == 1)
          
    call derv(y,rate,dydt,y0)
    call get_jacobian(y,jac,rate)
    
    do i=1,n_nequ
      jac(i,i)=jac(i,i)-dtm
      dydt(i)=dydt(i)-(y(i)-yin(i))*dtm
    end do
    y1(:)=-dydt(:)
    
    call linsys(jac,y1)

    y(1:n_spec)=y(1:n_spec) + y1(1:n_spec))

    yt(1:n_spec)=y1(1:n_spec)/y(1:n_spec)

    y(:)=max(y(:),1.e-100)

    if(all(abs(yt(:)) <= 0.0001)) iflag=0
    
    n=n+1

  end do

  return

end subroutine chemstep

!=======================================================================

  end module chemistry
