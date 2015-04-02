!=======================================================================
!> @file hlle.f90
!> @brief HLLE approximate Riemann solver module
!> @author C. Villarreal  D'Angelo, A. Esquivel, M. Schneiter
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

!> @brief HLLE approximate Riemann solver module
!! @details The module contains the routines needed to Solve the Riemann
!! problem in the entire domain and return the physical fluxes in x,y,z
!! with the HLLE solver

module hlle

#ifdef HLLE

contains

contains

!=======================================================================

!> @brief Solves the Riemann problem at the interface PL,PR
!! using the HLLE solver
!> @details Solves the Riemann problem at the interface betweem 
!! PL and PR using the HLLE solver
!> @n The fluxes are computed in the X direction, to obtain the
!! y ans z directions a swap is performed
!> @param real [in] primL : primitives at the Left state
!> @param real [in] primR : primitives at the Right state
!> @param real [out] ff : fluxes at the interface (@f$ F_{i+1/2} @f$)

subroutine prim2fhlle(priml,primr,ff)

  use parameters, only : neq
  use hydro_core, only : cfastX, prim2f, prim2u
  implicit none
  real, dimension(neq),intent(in   ) :: priml, primr
  real, dimension(neq),intent(inout) :: ff
  real, dimension(neq)               :: uR, uL, fL, fR
  real :: csl, csr, sl, sr

  call cfastX(priml,csl)
  call cfastX(primr,csr)

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

  end subroutine prim2fhlle

!=======================================================================

!> @brief Calculates HLLE fluxes from the primitive variables 
!!   on all the domain
!> @details Calculates HLLE fluxes from the primitive variables 
!!   on all the domain
!> @param integer [in] choice : 1, uses primit for the 1st half of timestep
!! (first order)
!!                  @n 2 uses primit for second order timestep

subroutine hllEfluxes(choice)

  use parameters, only : neq, nx, ny, nz
  use globals, only : primit, f, g, h
  use hydro_core, only : swapy, swapz, limiter
  implicit none
  integer, intent(in) :: choice
  integer :: i, j, k, ip, jp, kp, im, jm, km, ip2, jp2, kp2
  real, dimension(neq) :: priml, primr, primll, primrr, ff, uu
  !
  select case(choice)
 
  case(1)        ! 1st half timestep
 
     do i=0,nx
        do j=0,ny
           do k=0,nz
 
              ip=i+1
              jp=j+1
              kp=k+1
 
              !------- x direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,ip,j ,k )
 
              call prim2fhlle(priml,primr,ff)
              f(:,i,j,k)=ff(:)
              !------- y direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, jp,k )
              call swapy(priml,neq)          !swaps primL for L state
              call swapy(primr,neq)          !swaps primR for R state 
 
              call prim2fhlle(priml,primr,ff)  !gets fluxes (swapped)
              call swapy(ff,neq)             !swaps back the fluxes
              g(:,i,j,k)=ff(:)
              !------- z direction -------------------------------------
              priml(:)=primit(:,i ,j ,k )
              primr(:)=primit(:,i, j, kp)
              call swapz(priml,neq)
              call swapz(primr,neq)
              !
              call prim2fhlle(priml,primr,ff)
              call swapz(ff,neq)
              h(:,i,j,k)=ff(:)
 
           end do
        end do
     end do
 
  case (2)   !  2nd half timestep
 
     do i=0,nx
        do j=0,ny
           do k=0,nz
 
              ip=i+1
              ip2=i+2
              im=i-1
              jp=j+1
              jp2=j+2
              jm=j-1
              kp=k+1
              kp2=k+2
              km=k-1
 
              !------- x direction ------------------------------------
              priml (:)=primit(:,i,  j,k )
              primr (:)=primit(:,ip, j,k )
              primll(:)=primit(:,im, j,k )
              primrr(:)=primit(:,ip2,j,k )
              call limiter(primll,priml,primr,primrr,neq)
 
              call prim2fhlle(priml,primr,ff)
              f(:,i,j,k)=ff(:)
              !------- y direction ------------------------------------
              priml (:)=primit(:,i,j  ,k )
              primr (:)=primit(:,i,jp ,k )
              primll(:)=primit(:,i,jm ,k )
              primrr(:)=primit(:,i,jp2,k )
              call swapy(priml,neq)
              call swapy(primr,neq)
              call swapy(primll,neq)
              call swapy(primrr,neq)
              call limiter(primll,priml,primr,primrr,neq)
 
              call prim2fhlle(priml,primr,ff)
              call swapy(ff,neq)
              g(:,i,j,k)=ff(:)
              !------- z direction ------------------------------------
              priml (:)=primit(:,i,j,k  )
              primr (:)=primit(:,i,j,kp )
              primll(:)=primit(:,i,j,km )
              primrr(:)=primit(:,i,j,kp2)
              call swapz(priml,neq)
              call swapz(primr,neq)
              call swapz(primll,neq)
              call swapz(primrr,neq)
              call limiter(primll,priml,primr,primrr,neq)
 
              call prim2fhlle(priml,primr,ff)
              call swapz(ff,neq)
              h(:,i,j,k)=ff(:)
              !
           end do
        end do
     end do

  end select

end subroutine hllEfluxes

#endif

end module hlle

!=======================================================================
