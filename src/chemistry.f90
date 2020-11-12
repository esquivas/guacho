!=======================================================================
!> @file chemistry.f90
!> @brief chemistry  module
!> @author A. Castellanos, P. Rivera A. Rodriguez, A. Raga  and A. Esquivel
!> @date 10/Mar/2016
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

!> @brief chemistry module
!> @details module to solve the chemical/ionic network.

module chemistry

  use network
  implicit none

contains

  !=======================================================================
  !> @brief Advances the chemistry network
  !> @details Advances the chemistry network on the entire domain
  !! (except ghost cells), updates primitives and conserved variables
  !! in globals
  subroutine update_chem()

    use parameters, only : neq, neqdyn, n_spec, nx, ny, nz, tsc, rhosc
    use globals, only : u, primit, dt_CFL
    use network, only : n_elem
    use hydro_core, only : u2prim
    implicit none
    real :: dt_seconds, T, y(n_spec), y0(n_elem)
    integer :: i, j, k, l

    dt_seconds = dt_CFL*tsc

    do k=1,nz
      do j=1,ny
        do i=1,nx

          !   get the primitives (and T)
          call u2prim(u(:,i,j,k),primit(:,i,j,k),T)

          y(1:n_spec) = primit(neqdyn+1:neqdyn+n_spec,i,j,k)
          y0(1) = primit(1,i,j,k)

          !  update the passive primitives (should not work in single precision)
          call chemstep(y, y0, T, dt_seconds)

          !  update the primitives and conserved variables
          do l = 1, n_spec
            primit(l+neqdyn,i,j,k) = y(l)
            u     (l+neqdyn,i,j,k) = y(l)
          end do

        end do
      end do
    end do

  end subroutine update_chem

  !=======================================================================
  !> @brief Advances the chemistry network in one cell
  !> @details Advances the chemistry network on the in one cell
  !> @param real [inout] y(n_spec) : number densities of the species
  !> to be updated by the chemistry
  !> @param real [in] y[n_elem] : total number density of each of the
  !! elements involved in the reactions
  !> @param real [in] T : Temperature [K]
  !> @param real [in] deltt : time interval (from the hydro, in seconds)
  subroutine chemstep(y,y0,T, deltt)

    use linear_system
    use network, only : n_spec, n_reac, n_elem, get_reaction_rates,            &
                        derv, get_jacobian, n_nequ, check_no_conservation
    implicit none
    real (kind=8), intent(inout) :: y(n_spec)
    real (kind=8), intent(in)    ::    y0(n_elem), T, deltt
    real (kind=8) :: dtm
    real (kind=8) :: y1(n_spec),yt(n_spec),yin(n_spec), y0_in(n_elem)
    real (kind=8) :: rate(n_reac),dydt(n_spec),jacobian(n_spec,n_spec)
    integer, parameter  :: niter=100       ! number of iterations
    integer :: n,i,iff

    n=0
    dtm=1./deltt
    iff=1
    yin(:) =y (:)
    y0_in(:) = y0(:)

    call get_reaction_rates(rate,min(T,3E4))

    !  initial guess for Newton-Raphson
    if ( check_no_conservation(y,y0_in) ) then
      !print*, '*****Reset Initial Guess ********'
      !print*, "T=", T
      call nr_init(y,y0_in)
    end if

    do while ( n <= niter )

      call derv(y,rate,dydt,y0)
      call get_jacobian(y,jacobian,rate)

      do i=1,n_nequ
        jacobian(i,i)=jacobian(i,i)-dtm
        dydt(i)=dydt(i)-(y(i)-yin(i))*dtm
      end do
      y1(:)=-dydt(:)

      call linsys(jacobian,y1, n_spec)

      y(:) = y(:) + y1(:)
      y(:)=max(y(:),1.e-40)

      yt(:)=y1(:)/y(:)

      !  exit the loop if converged
      if(all(abs(y1(:)) <= 0.0001)) exit

      n=n+1

    end do

    !if (n >= niter) then
    !  print*, "failed to converge after ", niter, " iterations"
    !else
    !  print*, 'converged after ', n+1, ' iterations'
    !end if

    return

  end subroutine chemstep

  !=======================================================================

end module chemistry
