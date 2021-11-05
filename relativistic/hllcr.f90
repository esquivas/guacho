!=======================================================================
!> @file hllc.f90
!> @brief HLLC approximate Riemann solver module
!> @author Alejandro Esquivel
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

!> @brief HLLC approximate Riemann solver module
!! @details The module contains the routines needed to Solve the Riemann
!! problem in the entire domain and return the physical fluxes in x,y,z
!! with the HLLC solver

module hllcr

contains

!> @brief Solves the Riemann problem at the interface PL,PR
!! using the HLLC solver
!> @details Solves the Riemann problem at the interface betweem 
!! PL and PR using the HLLC solver
!> @n The fluxes are computed in the X direction, to obtain the
!! y ans z directions a swap is performed
!> @param real [in] primL : primitives at the Left state
!> @param real [in] primR : primitives at the Right state
!> @param real [out] ff : fluxes at the interface (@f$ F_{i+1/2} @f$)

subroutine prim2fhllcr(priml,primr,ff)

  use parameters, only : neq, neqdyn, cv, pmhd, passives
  use hydro_core, only : csound, prim2f, prim2u
  implicit none
  real, dimension(neq),intent(in   ) :: priml, primr
  real, dimension(neq),intent(inout) :: ff
  real, dimension(neq)               :: uu, uuk
  real :: csl, csr, sl, sr, slmul, srmur, rholul, rhorur, sst
  real :: rhost,ek
  
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
  
  slmul=sl-priml(2)
  srmur=sr-primr(2)
  rholul=priml(1)*priml(2)
  rhorur=primr(1)*primr(2)

  sst = (srmur*rhorur-slmul*rholul-primr(5)+priml(5) )        &  
        / (srmur*primr(1)-slmul*priml(1) )
  
  if (sst >= 0.) then
    rhost=priml(1)*(slmul)/(sl-sst)
    ek= 0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5)

    uuk(1)=rhost
    uuk(2)=rhost*sst
    uuk(3)=rhost*priml(3)
    uuk(4)=rhost*priml(4)
    uuk(5)=rhost*( ek/priml(1)+(sst-priml(2))*(sst+priml(5)/(priml(1)*slmul)) )

  if (pmhd) then
#ifdef BFIELD
    uuk(6:8)=rhost*priml(6:8)/priml(1)
#endif 
  end if
#ifdef PASSIVES
  if (passives) then
    uuk(neqdyn+1:neq)=rhost*priml(neqdyn+1:neq)/priml(1)
  end if
#endif

    call prim2f(priml,ff)
    call prim2u(priml,uu)
    ff(:)=ff(:) + sl*( uuk(:)-uu(:) )
    return
  endif

  if (sst <= 0.) then
    rhost=primr(1)*(srmur)/(sr-sst)
    ek= 0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5)

    uuk(1)=rhost
    uuk(2)=rhost*sst
    uuk(3)=rhost*primr(3)
    uuk(4)=rhost*primr(4)
    uuk(5)=rhost*( ek/primr(1)+(sst-primr(2))*(sst+primr(5)/(primr(1)*srmur)) )

  if (pmhd) then
#ifdef BFIELD
      !uuk(5)= 0.
    uuk(6:8)=rhost*primr(6:8)/primr(1)
#endif
  end if
  
#ifdef PASSIVES
  if (passives) then
    uuk(neqdyn+1:neq)=rhost*primr(neqdyn+1:neq)/primr(1)
  end if
#endif

    call prim2f(primr,ff)
    call prim2u(primr,uu)
    ff(:)=ff(:) + sr*( uuk(:)-uu(:) )
    return
  endif

  print*, 'Error in hllc'
  print*, 'primL: ',priml(:)
  print*, 'primR: ',primr(:)
  stop

end subroutine prim2fhllcr


!=======================================================================

!> @brief Calculates HLLC fluxes from the primitive variables 
!!   on all the domain
!> @details Calculates HLLC fluxes from the primitive variables 
!!   on all the domain
!> @param integer [in] choice : 1, uses primit for the 1st half of timestep (first order)
!!                  @n 2 uses primit for second order timestep

subroutine hllcrfluxes(choice)

  use parameters, only : neq, nx, ny, nz
  use globals, only : primit, f, g, h
  use hydro_core, only : swapy, swapz, limiter
  implicit none
  integer, intent(in) :: choice
  integer :: i, j, k
  real, dimension(neq) :: priml, primr, primll, primrr, ff
  !
  select case(choice)

  case(1)        ! 1st half timestep

     do k=0,nz
        do j=0,ny
           do i=0,nx

              !------- x direction -------------------------------------
              priml(:)=primit(:,i  ,j ,k )
              primr(:)=primit(:,i+1,j ,k )

              call prim2fhllcr(priml,primr,ff)
              f(:,i,j,k)=ff(:)

              !------- y direction -------------------------------------
              priml(:)=primit(:,i ,j  ,k )
              primr(:)=primit(:,i, j+1,k )
              call swapy(priml,neq)
              call swapy(primr,neq)

              call prim2fhllcr(priml,primr,ff)
              call swapy(ff,neq)
              g(:,i,j,k)=ff(:)

              !------- z direction -------------------------------------
              priml(:)=primit(:,i ,j ,k  )
              primr(:)=primit(:,i, j, k+1)
              call swapz(priml,neq)
              call swapz(primr,neq)

              call prim2fhllcr(priml,primr,ff)
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

              call prim2fhllcr(priml,primr,ff)
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

              call prim2fhllcr(priml,primr,ff)
              call swapy(ff,neq)
              g(:,i,j,k)=ff(:)

              !------- z direction ------------------------------------
              priml (:)=primit(:,i,j,k  )
              primr (:)=primit(:,i,j,k+1 )
              primll(:)=primit(:,i,j,k-1 )
              primrr(:)=primit(:,i,j,k+2)
              call swapz(priml,neq)
              call swapz(primr,neq)
              call swapz(primll,neq)
              call swapz(primrr,neq)
              call limiter(primll,priml,primr,primrr,neq)

              call prim2fhllcr(priml,primr,ff)
              call swapz(ff,neq)
              h(:,i,j,k)=ff(:)

           end do
        end do
     end do

  end select

end subroutine hllcrfluxes


end module hllcr

!=======================================================================
