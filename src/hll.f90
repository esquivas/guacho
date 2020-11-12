!=======================================================================
!> @file hll.f90
!> @brief HLL approximate Riemann solver module
!> @author Alejandro Esquivel
!> @date 4/May/2016
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

!> @brief HLL approximate Riemann solver module
!! @details The module contains the routines needed to Solve the Riemann
!! problem in the entire domain and return the physical fluxes in x,y,z
!! with the HLL solver

module hll

contains

  !=======================================================================
  !> @brief Solves the Riemann problem at the interface PL,PR
  !> using the HLL solver
  !> @details Solves the Riemann problem at the interface betweem
  !! PL and PR using the HLL solver
  !> @n The fluxes are computed in the X direction, to obtain the
  !! y and z directions a swap is performed
  !> @param real [in] primL : primitives at the Left state
  !> @param real [in] primR : primitives at the Right state
  !> @param real [out] ff : fluxes at the interface (@f$ F_{i+1/2} @f$)
  subroutine prim2fhll(priml,primr,ff)

    use parameters, only : neq
    use hydro_core, only : csound, prim2f, prim2u
    implicit none
    real, dimension(neq),intent(in   ) :: priml, primr
    real, dimension(neq),intent(inout) :: ff
    real, dimension(neq)               :: uR, uL, fL, fR
    real :: csl, csr, sl, sr

    call csound(priml(5),priml(1),csl)
    call csound(primr(5),primr(1),csr)

    sr=max(priml(2)+csl,primr(2)+csr)
    sl=min(priml(2)-csl,primr(2)-csr)

    if (sl > 0) then
       call prim2f(priml,ff)
       return
    endif

    if (sr < 0) then
       call prim2f(primr,ff)
       return
    endif

    call prim2f(priml,fL)
    call prim2f(primr,fR)
    call prim2u(priml,uL)
    call prim2u(primr,uR)

    ff(:)=(sr*fL(:)-sl*fR(:)+sl*sr*(uR(:)-uL(:)))/(sr-sl)

    return

  end subroutine prim2fhll

  !=======================================================================
  !> @brief Calculates HLL fluxes from the primitive variables
  !>   on all the domain
  !> @details Calculates HLL fluxes from the primitive variables
  !! on all the domain
  !> @param integer [in] choice : 1, uses primit for the 1st half of timestep
  !! (first order)
  !> @n 2 uses primit for second order timestep
  subroutine hllfluxes(choice)

    use parameters, only : neq, nx, ny, nz
    use globals, only : primit, f, g, h
    use hydro_core, only : swapy, swapz, limiter
    implicit none
    integer, intent(in) :: choice
    integer :: i, j, k
    real, dimension(neq) :: priml, primr, primll, primrr, ff

    select case(choice)

    case(1)        ! 1st half timestep

      do k=0,nz
        do j=0,ny
          do i=0,nx

            !------- x direction -------------------------------------
            priml(:)=primit(:,i  ,j ,k )
            primr(:)=primit(:,i+1,j ,k )

            call prim2fhll(priml,primr,ff)
            f(:,i,j,k)=ff(:)

            !------- y direction -------------------------------------
            priml(:)=primit(:,i ,j  ,k )
            primr(:)=primit(:,i, j+1,k )
            call swapy(priml,neq)
            call swapy(primr,neq)

            call prim2fhll(priml,primr,ff)
            call swapy(ff,neq)
            g(:,i,j,k)=ff(:)

            !------- z direction -------------------------------------
            priml(:)=primit(:,i ,j ,k  )
            primr(:)=primit(:,i, j, k+1)
            call swapz(priml,neq)
            call swapz(primr,neq)

            call prim2fhll(priml,primr,ff)
            call swapz(ff,neq)
            h(:,i,j,k)=ff(:)

          end do
        end do
      end do

    case (2)   !  2nd half timestep

      do k=0,nz
        do j=0,ny
          do i=0,nx

            !------- x direction ------------------------------------
            priml (:)=primit(:,i,  j,k )
            primr (:)=primit(:,i+1,j,k )
            primll(:)=primit(:,i-1,j,k )
            primrr(:)=primit(:,i+2,j,k )
            call limiter(primll,priml,primr,primrr,neq)

            call prim2fhll(priml,primr,ff)
            f(:,i,j,k)=ff(:)

            !------- y direction ------------------------------------
            priml (:)=primit(:,i,j  ,k )
            primr (:)=primit(:,i,j+1,k )
            primll(:)=primit(:,i,j-1,k )
            primrr(:)=primit(:,i,j+2,k )
            call swapy(priml,neq)
            call swapy(primr,neq)
            call swapy(primll,neq)
            call swapy(primrr,neq)
            call limiter(primll,priml,primr,primrr,neq)

            call prim2fhll(priml,primr,ff)
            call swapy(ff,neq)
            g(:,i,j,k)=ff(:)

            !------- z direction ------------------------------------
            priml (:)=primit(:,i,j,k  )
            primr (:)=primit(:,i,j,k+1)
            primll(:)=primit(:,i,j,k-1)
            primrr(:)=primit(:,i,j,k+2)
            call swapz(priml,neq)
            call swapz(primr,neq)
            call swapz(primll,neq)
            call swapz(primrr,neq)
            call limiter(primll,priml,primr,primrr,neq)

            call prim2fhll(priml,primr,ff)
            call swapz(ff,neq)
            h(:,i,j,k)=ff(:)

          end do
        end do
      end do

    end select

  end subroutine hllfluxes

  !=======================================================================

end module hll
