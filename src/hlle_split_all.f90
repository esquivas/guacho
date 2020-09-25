!=======================================================================
!> @file hlleSplitall.f90
!> @brief HLLE approximate Riemann solver module split All version
!> @author Valeria Sieyra, Matias Schneiter, Alejandro Esquivel
!> @date 04/May/2016

! Copyright (c) 2020 Guacho Co-Op.
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

!> @brief HLLE approximate Riemann solver module, slit Version
!! @details The module contains the routines needed to Solve the Riemann
!! problem in the entire domain and return the physical fluxes in x,y,z
!! with the HLLE solver

module hlleSplitAll

#ifdef BFIELD

contains

  !=======================================================================
  !> @brief Solves the Riemann problem at the interface PL,PR
  !! using the HLLE solver with split in all variables
  !> @details Solves the Riemann problem at the interface betweem
  !! PL and PR using the HLLE solver
  !> @n The fluxes are computed in the X direction, to obtain the
  !! y and z directions a swap is performed
  !> @param real [in] primL : primitives at the Left state (fluctuation)
  !> @param real [in] primR : primitives at the Right state (fluctuation)
  !> @param real [in] prim0L : primitives at the Left state (background)
  !> @param real [in] prim0R : primitives at the Right state (background)
  !> @param real [out] ff : fluxes at the interface (@f$ F_{i+1/2} @f$)
  subroutine prim2fhlleSplitAll(priml,primr,prim0l,prim0r,ff)

    use parameters, only : neq
    use hydro_core, only : cfastX, prim2f, prim2u
    implicit none
    real, dimension(neq),intent(in   ) :: priml, primr, prim0l, prim0r
    real, dimension(neq),intent(inout) :: ff
    real, dimension(neq)               :: uR, uL, fL, fR
    real :: csl, csr, sl, sr

    call cfastX(priml+prim0l,csl)
    call cfastX(primr+prim0r,csr)

    sr=max(priml(2)+prim0l(2)+csl,primr(2)+prim0r(2)+csr)
    sl=min(priml(2)+prim0l(2)-csl,primr(2)+prim0r(2)-csr)

    if (sl > 0) then
      call prim2f(priml,ff,prim0=prim0l)
      return
    endif

    if (sr < 0) then
      call prim2f(primr,ff,prim0=prim0r)
      return
    endif

    call prim2f(priml,fL,prim0=prim0l)
    call prim2f(primr,fR,prim0=prim0r)
    call prim2u(priml,uL, prim0=prim0L)
    call prim2u(primr,uR, prim0=prim0R)

    ff(:)=(sr*fL(:)-sl*fR(:)+sl*sr*(uR(:)-uL(:)))/(sr-sl)

  end subroutine prim2fhlleSplitAll

  !=======================================================================
  !> @brief Calculates HLLE fluxes from the primitive variables
  !! on all the domain
  !> @details Calculates HLLE fluxes from the primitive variables
  !! on all the domain, split version
  !> @param integer [in] choice : 1, uses primit for the 1st half of timestep
  !! (first order)
  !> @n 2 uses primit for second order timestep
  subroutine hllEfluxesSplitAll(choice)

    use parameters, only : neq, nx, ny, nz
    use globals, only : primit, f, g, h, primit0
    !  use hydro_core, only : swapy, swapz, swapy_bsplit, swapz_bsplit, limiter
    use hydro_core, only : swapy, swapz, limiter
    implicit none
    integer, intent(in) :: choice
    integer :: i, j, k
    real, dimension(neq) :: priml, primr, primll, primrr, ff, prim0l, prim0r,  &
                            prim0ll, prim0rr
    !  real, dimension(3)   :: 0l, B0r, B0ll, B0rr

    select case(choice)

    case(1)        ! 1st half timestep

      do k=0,nz
        do j=0,ny
          do i=0,nx

            !------- x direction -------------------------------------
            priml(:)=primit(:,i  ,j ,k )
            primr(:)=primit(:,i+1,j ,k )
            prim0l(:)=primit0(:,i ,j ,k)
            prim0r(:)=primit0(:,i+1 ,j ,k)

            call prim2fhlleSplitAll(priml,primr,prim0l,prim0r,ff)
            f(:,i,j,k)=ff(:)

            !------- y direction -------------------------------------
            priml(:)=primit(:,i ,j  ,k )
            primr(:)=primit(:,i, j+1,k )
            prim0l(:)=primit0(:,i ,j   ,k)
            prim0r(:)=primit0(:,i ,j+1 ,k)

            call swapy(priml,neq)
            call swapy(primr,neq)
            call swapy(prim0l,neq)
            call swapy(prim0r,neq)

            call prim2fhlleSplitAll(priml,primr,prim0l,prim0r,ff)
            call swapy(ff,neq)
            g(:,i,j,k)=ff(:)

            !------- z direction -------------------------------------
            priml(:)=primit(:,i ,j ,k  )
            primr(:)=primit(:,i, j ,k+1)
            prim0l(:)=primit0(:,i ,j ,k)
            prim0r(:)=primit0(:,i ,j ,k+1)

            call swapz(priml,neq)
            call swapz(primr,neq)
            call swapz(prim0l,neq)
            call swapz(prim0r,neq)

            call prim2fhlleSplitAll(priml,primr,prim0l,prim0r,ff)
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
            prim0l (:)=primit0(:,i   ,j ,k)
            prim0r (:)=primit0(:,i+1 ,j ,k)
            prim0ll(:)=primit0(:,i-1 ,j ,k)
            prim0rr(:)=primit0(:,i+2 ,j ,k)

            call limiter(primll,priml,primr,primrr,neq)
            call limiter(prim0ll,prim0l,prim0r,prim0rr,neq)

            call prim2fhlleSplitAll(priml,primr,prim0l,prim0r,ff)
            f(:,i,j,k)=ff(:)

            !------- y direction ------------------------------------
            priml (:)=primit(:,i,j  ,k )
            primr (:)=primit(:,i,j+1,k )
            primll(:)=primit(:,i,j-1,k )
            primrr(:)=primit(:,i,j+2,k )
            prim0l (:)=primit0(:,i,j  ,k )
            prim0r (:)=primit0(:,i,j+1,k )
            prim0ll(:)=primit0(:,i,j-1,k )
            prim0rr(:)=primit0(:,i,j+2,k )

            call swapy(priml,neq)
            call swapy(primr,neq)
            call swapy(primll,neq)
            call swapy(primrr,neq)

            call swapy(prim0l,neq)
            call swapy(prim0r,neq)
            call swapy(prim0ll,neq)
            call swapy(prim0rr,neq)

            call limiter(primll,priml,primr,primrr,neq)
            call limiter(prim0ll,prim0l,prim0r,prim0rr,neq)

            call prim2fhlleSplitAll(priml,primr,prim0l,prim0r,ff)
            call swapy(ff,neq)
            g(:,i,j,k)=ff(:)

            !------- z direction ------------------------------------
            priml (:)=primit(:,i,j,k  )
            primr (:)=primit(:,i,j,k+1)
            primll(:)=primit(:,i,j,k-1)
            primrr(:)=primit(:,i,j,k+2)
            prim0l (:)=primit0(:,i,j,k  )
            prim0r (:)=primit0(:,i,j,k+1)
            prim0ll(:)=primit0(:,i,j,k-1)
            prim0rr(:)=primit0(:,i,j,k+2)

            call swapz(priml,neq)
            call swapz(primr,neq)
            call swapz(primll,neq)
            call swapz(primrr,neq)

            call swapz(prim0l,neq)
            call swapz(prim0r,neq)
            call swapz(prim0ll,neq)
            call swapz(prim0rr,neq)

            call limiter(primll,priml,primr,primrr,neq)
            call limiter(prim0ll,prim0l,prim0r,prim0rr,neq)

            call prim2fhlleSplitAll(priml,primr,prim0l,prim0r,ff)
            call swapz(ff,neq)
            h(:,i,j,k)=ff(:)

          end do
        end do
      end do

    end select

  end subroutine hllEfluxesSplitAll

  !=======================================================================

#endif

end module hlleSplitAll
