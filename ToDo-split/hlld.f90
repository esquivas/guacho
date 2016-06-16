!=======================================================================
!> @file hlld.f90
!> @brief HLLD approximate Riemann solver module
!> @author  C. Villarreal  D'Angelo, A. Esquivel, M. Schneiter
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

!> @brief HLLD approximate Riemann solver module
!! @details The module contains the routines needed to Solve the Riemann
!! problem in the entire domain and return the physical fluxes in x,y,z
!! with the HLLD solver

module hlld

#ifdef HLLD

contains

!=======================================================================

!> @brief Solves the Riemann problem at the interface PL,PR
!! using the HLLD solver
!> @details Solves the Riemann problem at the interface betweem 
!! PL and PR using the HLLD solver
!> @n The fluxes are computed in the X direction, to obtain the
!! y ans z directions a swap is performed
!> @param real [in] primL : primitives at the Left state
!> @param real [in] primR : primitives at the Right state
!> @param real [out] ff : fluxes at the interface (@f$ F_{i+1/2} @f$)
#ifndef MHD_BSPLIT
subroutine prim2fhlld(priml,primr,ff)

  use parameters, only : neq, cv
  use hydro_core, only : cfastX, prim2f
  implicit none
  real, dimension(neq),intent(in   ) :: priml, primr   
  real, dimension(neq),intent(inout) :: ff
  real, dimension(neq)               :: pp
  real, dimension(neq)               ::fL, fR, uL, uR
  real, dimension(neq)               :: uu, ust
  real :: csl, csr, sl, sr, slmul, srmur, rholul, rhorur, sM
  real :: pTL, pTR, Bx, signBx
  real :: slmsM, srmsM, rhostl, rhostr, sstl, sstr
  real :: pst, el, er, denl, denr, sMmul, sMmur
  real :: vstl, wstl, bystl, bzstl, estl, vdotbl, vstdotbstl
  real :: vstr, wstr, bystr, bzstr, estr, vdotbr, vstdotbstr
  real :: sMmsstl, sMmsstr
  real :: dd, vstst, wstst, bystst, bzstst
  real ::  vststdotbstst, eststl, eststr
  integer :: err

  call cfastX(priml,csl)
  call cfastX(primr,csr)

  sr=max(priml(2)+csl,primr(2)+csr)
  sl=min(priml(2)-csl,primr(2)-csr)
!   print*,'sL',sl,'sR',sr

  ! UL region -----------------------------------
  if (sl > 0) then
     call prim2f(priml,ff)
     return
  endif

  ! UR region -----------------------------------
  if (sr < 0) then
     call prim2f(primr,ff)
     return
  endif

  Bx= 0.5* (primL(6)+primR(6) )
  signBx= sign(1.,Bx)

  !  Total pressure
  pTL=primL(5) + 0.5*( bx**2+primL(7)**2+primL(8)**2 )
  pTR=primR(5) + 0.5*( bx**2+primR(7)**2+primR(8)**2 )

  slmul=sl-priml(2)  ! sL-uL
  srmur=sr-primr(2)  ! sR-uR

  rholul=priml(1)*priml(2)  !rhoL*uL
  rhorur=primr(1)*primr(2)  !rhoR*uR

  sM = (srmur*rhorur-slmul*rholul-pTR+pTL)/( srmur*primr(1)-slmul*priml(1) )
!   print*,'sM',sM, 'uL',priml(2),'uR',primr(2)

  srmsM=sr-sM
  slmsM=sl-sM

  rhostl=priml(1)*slmul/slmsM      !rhoL*
  rhostr=primr(1)*srmur/srmsM      !rhoR*

  sstl=sM - abs(bx)/sqrt(rhostl)  !SL*
  sstr=sM + abs(bx)/sqrt(rhostr)  !SR*

  pst= (srmur*primr(1)*pTL - slmul*priml(1)*pTR               &  !pT* 
        + priml(1)*primr(1)*srmur*slmul*(primr(2)-priml(2)) ) &
        /( srmur*primr(1)-slmul*priml(1) )

  ! UL* region -----------------------------------
  if(sstl >= 0) then
    !
    el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &  !eL
               +0.5*(bx**2+priml(7)**2+priml(8)**2)

    sMmuL=sM - priml(2)
    denl=priml(1)*slmul*slmsM-bx**2

    if(denl == 0) then  !si sL=uL, sL=sM y Bx=0 (si pasa esto las cantidades vL*,wL*,ByL*,BzL*  quedarían 0/0.
      vstl = primL(3)
      wstl = primL(4)
      bystl= 0.! primL(7)
      bzstl= 0. !primL(8)
!       print*,'stopped @ HLLD','denl*','sL',sl,'ul',ul,'sM',sM,'Bx',bx
!       stop
    else

      vstl = priml(3) - bx*priml(7)*sMmul/denl                      !vL*
      wstl = priml(4) - bx*priml(8)*sMmul/denl                      !wL*
      bystl= priml(7)*( priml(1)*slmul**2 - bx**2 )/denl            !byL*
      bzstl= priml(8)*( priml(1)*slmul**2 - bx**2 )/denl            !bzL*
      
    endif
    
    vdotbl    = priml(2)*bx + priml(3)*priml(7) + priml(4)*priml(8) !vL dot BL
    vstdotbstl= sM*bx + vstl*bystl + wstl*bzstl                     !vL* dot BL*  

    estl= ( slmul*el -pTL*priml(2) +pst*sM +bx*(vdotbl-vstdotbstl) )/slmsM !eL*

    ff(1) = rhostl*sM
    ff(2) = rhostl*SM**2+pst-bx**2
    ff(3) = rhostl*sM*vstl-bx*bystl
    ff(4) = rhostl*sM*wstl-bx*bzstl
    ff(5) = sM*(estl+pst)-bx*(vstdotbstl)
    ff(6) = 0.
    ff(7) = bystl*sM-bx*vstl
    ff(8) = bzstl*sM-bx*wstl

#ifdef PASSIVES
      ff(neqdyn+1:neq)=sM*priml(neqdyn+1:neq)*slmul/slmsM
#endif   

    return
  endif

  ! UR* region -----------------------------------
  if(sstr <= 0) then

    er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &  !eR
               +0.5*(bx**2+primr(7)**2+primr(8)**2)

    sMmuR=sM - primr(2)
    denr=primr(1)*srmur*sRmsM-bx**2

    if(denr == 0) then
      vstr = primR(3)
      wstr = primR(4)
      bystr= 0.!primR(7)
      bzstr= 0.!primR(8)
!       print*,'stopped @ HLLD','denr*','sr',sr,'ur',ur,'sM',sM,'Bx',bx
! !       stop
    else

      vstr = primr(3) - bx*primr(7)*sMmuR/denr                      !vR*
      wstr = primr(4) - bx*primr(8)*sMmuR/denr                      !wR*
      bystr= primr(7)*( primr(1)*srmur**2 - bx**2 )/denr            !byR*
      bzstr= primr(8)*( primr(1)*srmur**2 - bx**2 )/denr            !bzR*

    end if
   
    vdotbr    = primr(2)*bx + primr(3)*primr(7) + primr(4)*primr(8) !vR dot BR
    vstdotbstr= sM*bx + vstr*bystr + wstr*bzstr                     !vR* dot BR*  

    estr= ( srmur*er -pTR*primr(2) +pst*sM +bx*(vdotbr-vstdotbstr) )/srmsM !eR*

    ff(1) = rhostr*sM
    ff(2) = rhostr*SM**2+pst-bx**2
    ff(3) = rhostr*sM*vstr-bx*bystr
    ff(4) = rhostr*sM*wstr-bx*bzstr
    ff(5) = sM*(estr+pst)-bx*(vstdotbstr)
    ff(6) = 0.
    ff(7) = bystr*sM-bx*vstr
    ff(8) = bzstr*sM-bx*wstr

#ifdef PASSIVES
    ff(neqdyn+1:neq)=sM*primr(neqdyn+1:neq)*srmur/srmsM
#endif

    return
  endif

  !   All this are needed on both the UL** and UR** regions
  sMmul= sM - priml(2)
  sMmur= sM - primr(2)

  denl=priml(1)*slmul*slmsM-bx**2 !+1.e-30
  denr=primr(1)*srmur*srmsM-bx**2 !+ 1.e-30

  if(denl == 0) then
    vstl =priml(3)
    wstl =priml(4)
    bystl=0.!priml(7)
    bzstl=0.!priml(8)
!     print*,'stopped @ HLLD','denl**',denl,'rhoL',priml(1),'sL',sl,'uL',priml(2),'sM',sM,'Bx',bx
!     stop
  else
    vstl = priml(3) - bx*priml(7)*sMmul/denl                      !vL*
    wstl = priml(4) - bx*priml(8)*sMmul/denl                      !wL*
    bystl= priml(7)*( priml(1)*slmul**2 - bx**2 )/denl            !byL*
    bzstl= priml(8)*( priml(1)*slmul**2 - bx**2 )/denl            !bzL*
  endif

  if(denr == 0) then
    vstr =primr(3)
    wstr =primr(4)
    bystr=0.!primr(7)
    bzstr=0.!primr(8)
!     print*,'stopped @ HLLD','denr**','sr',sr,'ur',primr(2),'sM',sM,'Bx',bx
!     stop
  else
    vstr = primr(3) - bx*primr(7)*sMmuR/denr                      !vR*
    wstr = primr(4) - bx*primr(8)*sMmuR/denr                      !wR*
    bystr= primr(7)*( primr(1)*srmur**2 - bx**2 )/denr            !byR*
    bzstr= primr(8)*( primr(1)*srmur**2 - bx**2 )/denr            !bzR*
  endif

  dd=sqrt(rhostl)+sqrt(rhostr)

  vstst =(sqrt(rhostl)*vstl + sqrt(rhostr)*vstr + (bystr-bystl)*signBx )/dd  !v** 
  wstst =(sqrt(rhostl)*wstl + sqrt(rhostr)*wstr + (bzstr-bzstl)*signBx )/dd  !w**

  bystst=(sqrt(rhostl)*bystr + sqrt(rhostr)*bystl +  &      !by**
          sqrt(rhostl*rhostr)*(vstr-vstl)*signBx )/dd                            
 
  bzstst=(sqrt(rhostl)*bzstr + sqrt(rhostr)*bzstl +  &      !bz**
          sqrt(rhostl*rhostr)*(wstr-wstl)*signBx )/dd

  vststdotbstst= sM*bx + vstst*bystst + wstst*bzstst        !v** dot B**

  ! UL** region -----------------------------------
  if(sM >= 0 ) then
    !
    el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &  !eL
                +0.5*(bx**2+priml(7)**2+priml(8)**2)

    vdotbl    = priml(2)*bx + priml(3)*priml(7) + priml(4)*priml(8) !vL dot BL
    vstdotbstl= sM*bx+ vstl*bystl + wstl*bzstl                      !vL* dot BL*  

    estl= ( slmul*el -pTL*priml(2) +pst*sM +bx*(vdotbl-vstdotbstl) )/slmsM !eL*

    eststl= estl - sqrt(rhostl)*(vstdotbstl-vststdotbstst)*signBx !eL**

    ff(1) = rhostl*sM
    ff(2) = rhostl*SM**2+pst-bx**2
    ff(3) = rhostl*sM*vstst-bx*bystst
    ff(4) = rhostl*sM*wstst-bx*bzstst
    ff(5) = sm*(eststl+pst)-bx*(vststdotbstst)
    ff(6) = 0.
    ff(7) = bystst*sM-bx*vstst
    ff(8) = bzstst*sM-bx*wstst

#ifdef PASSIVES
    ff(neqdyn+1:neq) = sM*priml(neqdyn+1:neq)*slmul/slmsM
#endif

  return
  endif

  ! UR** region -----------------------------------
  if(sM <= 0 ) then

    er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &  !eR
                +0.5*(bx**2+primr(7)**2+primr(8)**2)

    vdotbr    = primr(2)*bx + primr(3)*primr(7) + primr(4)*primr(8) !vR dot BR
    vstdotbstr= sM*bx+ vstr*bystr + wstr*bzstr                      !vR* dot BR*  

    estr= ( srmur*er -pTR*primr(2) +pst*sM +bx*(vdotbr-vstdotbstr) )/srmsM !eR*

    eststr= estr + sqrt(rhostr)*(vstdotbstr-vststdotbstst)*signBx !eR**

    ff(1) = rhostr*sM
    ff(2) = rhostr*SM**2+pst-bx**2
    ff(3) = rhostr*sM*vstst-bx*bystst
    ff(4) = rhostr*sM*wstst-bx*bzstst
    ff(5) = sm*(eststr+pst)-bx*(vststdotbstst)
    ff(6) = 0.
    ff(7) = bystst*sM-bx*vstst
    ff(8) = bzstst*sM-bx*wstst

#ifdef PASSIVES
    ff(neqdyn+1:neq) = sM*primr(neqdyn+1:neq)*srmur/srmsM
#endif

  return
  endif

  print'(a,5es12.3)', 'Error in HLLD routine', sM,sl,sr, csl, csr
  stop

end subroutine prim2fhlld

!=======================================================================

!> @brief Calculates HLLD fluxes from the primitive variables 
!!   on all the domain
!> @details Calculates HLLD fluxes from the primitive variables 
!!   on all the domain
!> @param integer [in] choice : 1, uses primit for the 1st half of timestep
!! (first order)
!!                  @n 2 uses primit for second order timestep

subroutine hlldfluxes(choice)

  use parameters, only : neq, nx, ny, nz
  use globals, only : primit, f, g, h
  use hydro_core, only : swapy, swapz, limiter
  implicit none
  integer, intent(in) :: choice
  integer :: i, j, k
  real, dimension(neq) :: priml, primr, primll, primrr, ff, uu

  select case(choice)

  case(1)        ! 1st half timestep
     !
     do k=0,nz
        do j=0,ny
           do i=0,nx

              !------- x direction -------------------------------------
              priml(:)=primit(:,i  ,j ,k )
              primr(:)=primit(:,i+1,j ,k )
              !
              call prim2fhlld(priml,primr,ff)
              f(:,i,j,k)=ff(:)

              !------- y direction -------------------------------------
              priml(:)=primit(:,i ,j  ,k )
              primr(:)=primit(:,i, j+1,k )
              call swapy(priml,neq)          !swaps primL for L state
              call swapy(primr,neq)          !swaps primR for R state 
              !
              call prim2fhlld(priml,primr,ff)  !gets fluxes (swapped)
              call swapy(ff,neq)             !swaps back the fluxes
              g(:,i,j,k)=ff(:)

              !------- z direction -------------------------------------
              priml(:)=primit(:,i ,j ,k  )
              primr(:)=primit(:,i, j, k+1)
              call swapz(priml,neq)
              call swapz(primr,neq)

              call prim2fhlld(priml,primr,ff)
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

              call prim2fhlld(priml,primr,ff)
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

              call prim2fhlld(priml,primr,ff)
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

              call prim2fhlld(priml,primr,ff)
              call swapz(ff,neq)
              h(:,i,j,k)=ff(:)

           end do
        end do
     end do

  end select

end subroutine hlldfluxes

#endif
!=======================================================================
#ifdef MHD_BSPLIT 
!Miyoshi 2010
subroutine prim2fhlld(priml,primr,B0l,B0r,ff)

  use parameters, only : neq, cv
  use hydro_core, only : cfastX_bsplit, prim2f_bsplit
  implicit none
  real, dimension(neq),intent(in   ) :: priml, primr
  real, dimension(3),intent(in   )   :: B0l, B0r
  real, dimension(neq),intent(inout) :: ff
  real, dimension(neq)               :: pp
  real, dimension(neq)               ::fL, fR, uL, uR
  real, dimension(neq)               :: uu, ust
  real :: csl, csr, sl, sr, slmul, srmur, rholul, rhorur, sM
  real :: pTL, pTR, Bx, Bxtot, signBx, signBxtot
  real :: slmsM, srmsM, rhostl, rhostr, sstl, sstr
  real :: pst, el, er, denl, denr, sMmul, sMmur
  real :: vstl, wstl, bystl, bzstl, estl, vdotbl, vstdotbstl
  real :: vstr, wstr, bystr, bzstr, estr, vdotbr, vstdotbstr
  real :: sMmsstl, sMmsstr
  real :: dd, vstst, wstst, bystst, bzstst
  real ::  vststdotbstst, eststl, eststr
  integer :: err

  call cfastX_bsplit(priml,B0l,csl) ! según el paper de Miyoshi et al 2010 esta se calculan 
  call cfastX_bsplit(primr,B0r,csr) ! solo con B1

  sr=max(priml(2)+csl,primr(2)+csr)
  sl=min(priml(2)-csl,primr(2)-csr)

  ! UL region -----------------------------------
  if (sl > 0) then
     call prim2f_bsplit(priml,B0l,ff)
     return
  endif

  ! UR region -----------------------------------
  if (sr < 0) then
     call prim2f_bsplit(primr,B0r,ff)
     return
  endif

  Bx= 0.5* (primL(6)+primR(6))  ! porque calculamos el promedio??? 19/04/16
  Bxtot= 0.5* (primL(6)+B0l(1)+primR(6)+B0r(1))
  signBx= sign(1.,Bx)
  signBxtot= sign(1.,Bxtot)

  !  Total pressure
  pTL=primL(5) + 0.5*( bx**2+primL(7)**2+primL(8)**2 + bx*B0l(1) + primL(7)*B0l(2) + primL(8)*B0l(3))
  pTR=primR(5) + 0.5*( bx**2+primR(7)**2+primR(8)**2 + bx*B0r(1) + primR(7)*B0r(2) + primR(8)*B0r(3))

  slmul=sl-priml(2)  ! sL-uL
  srmur=sr-primr(2)  ! sR-uR

  rholul=priml(1)*priml(2)  !rhoL*uL
  rhorur=primr(1)*primr(2)  !rhoR*uR

  sM = (srmur*rhorur-slmul*rholul-pTR+pTL)/( srmur*primr(1)-slmul*priml(1) )

  srmsM=sr-sM
  slmsM=sl-sM

  rhostl=priml(1)*slmul/slmsM      !rhoL*
  rhostr=primr(1)*srmur/srmsM      !rhoR*

  sstl=sM - abs(bx+B0l(1))/sqrt(rhostl)  !SL*
  sstr=sM + abs(bx+B0r(1))/sqrt(rhostr)  !SR*

  pst= (srmur*primr(1)*pTL - slmul*priml(1)*pTR               &  !pT* 
        + priml(1)*primr(1)*srmur*slmul*(primr(2)-priml(2)) ) &
        /( srmur*primr(1)-slmul*priml(1) )

  ! UL* region -----------------------------------
  if(sstl >= 0) then
    !
    el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &  !eL
               +0.5*(bx**2+priml(7)**2+priml(8)**2)

    sMmuL=sM - priml(2)
    denl=priml(1)*slmul*slmsM-bxtot**2

    if(denl == 0) then  !si sL=uL, sL=sM y Bx=0 (si pasa esto las cantidades vL*,wL*,ByL*,BzL*  quedarían 0/0.
      vstl = primL(3)
      wstl = primL(4)
      bystl= 0.!primL(7)
      bzstl= 0.!primL(8)
!       print*,'stopped @ HLLD','denl*','sL - uL',slmul,'sl - sM',slmsM,'Bxtot',bxtot
!       stop
    else

      vstl = priml(3) - bxtot*(priml(7)+B0l(2))*sMmul/denl                      !vL*
      wstl = priml(4) - bxtot*(priml(8)+B0l(3))*sMmul/denl                      !wL*
      bystl= (priml(7)+B0l(2))*( priml(1)*slmul**2 - bxtot**2 )/denl - B0l(2)            !byL*
      bzstl= (priml(8)+B0l(3))*( priml(1)*slmul**2 - bxtot**2 )/denl - B0l(3)           !bzL*

    endif
        
    vdotbl    = priml(2)*bx + priml(3)*priml(7) + priml(4)*priml(8) !vL dot BL
    vstdotbstl= sM*bx + vstl*bystl + wstl*bzstl                     !vL* dot BL*  

    estl= ( slmul*el -pTL*priml(2) +pst*sM +bxtot*(vdotbl-vstdotbstl) )/slmsM !eL*

    ff(1) = rhostl*sM
    ff(2) = rhostl*sM**2+pst-bxtot*bx-bx*B0l(1)
    ff(3) = rhostl*sM*vstl-bxtot*bystl-bx*B0l(2)
    ff(4) = rhostl*sM*wstl-bxtot*bzstl-bx*B0l(3)
    ff(5) = sM*(estl+pst)-bxtot*(vstdotbstl)
    ff(6) = 0.
    ff(7) = (bystl+B0l(2))*sM-bxtot*vstl
    ff(8) = (bzstl+B0l(3))*sM-bxtot*wstl

#ifdef PASSIVES
      ff(neqdyn+1:neq)=sM*priml(neqdyn+1:neq)*slmul/slmsM  
#endif

  return
  endif

  ! UR* region -----------------------------------
  if(sstr <= 0) then

    er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &  !eR
               +0.5*(bx**2+primr(7)**2+primr(8)**2)

    sMmuR=sM - primr(2)
    denr=primr(1)*srmur*sRmsM-bxtot**2

    if(denr == 0) then
      vstr = primR(3)
      wstr = primR(4)
      bystr= 0.!primR(7)
      bzstr= 0.!primR(8)
!       print*,'stopped @ HLLD','denr*','sr - ur',srmur,'sr - sM',sRmsM,'Bxtot',bxtot
!       stop
    else

      vstr = primr(3) - bxtot*(primr(7)+B0r(2))*sMmuR/denr                   !vR*
      wstr = primr(4) - bxtot*(primr(8)+B0r(3))*sMmuR/denr                      !wR*
      bystr= (primr(7)+B0r(2))*( primr(1)*srmur**2 - bxtot**2 )/denr - B0r(2)
      bzstr= (primr(8)+B0r(3))*( primr(1)*srmur**2 - bxtot**2 )/denr - B0r(3)           !bzR*

    end if
    
    vdotbr    = primr(2)*bx + primr(3)*primr(7) + primr(4)*primr(8) !vR dot BR
    vstdotbstr= sM*bx + vstr*bystr + wstr*bzstr                     !vR* dot BR*  

    estr= ( srmur*er -pTR*primr(2) +pst*sM +bxtot*(vdotbr-vstdotbstr) )/srmsM !eR*
    
    ff(1) = rhostr*sM
    ff(2) = rhostr*sM**2+pst-bxtot*bx-bx*B0r(1)
    ff(3) = rhostr*sM*vstr-bxtot*bystr-bx*B0r(2)
    ff(4) = rhostr*sM*wstr-bxtot*bzstr-bx*B0r(3)
    ff(5) = sM*(estr+pst)-bxtot*(vstdotbstr)
    ff(6) = 0.
    ff(7) = (bystr+B0r(2))*sM-bxtot*vstr
    ff(8) = (bzstr+B0r(3))*sM-bxtot*wstr

#ifdef PASSIVES
      ff(neqdyn+1:neq)=sM*primr(neqdyn+1:neq)*srmur/srmsM
#endif
    
    return
  endif

  !   All this are needed on both the UL** and UR** regions
  sMmul= sM - priml(2)
  sMmur= sM - primr(2)

  denl=priml(1)*slmul*slmsM-bxtot**2
  denr=primr(1)*srmur*srmsM-bxtot**2

  if(denl == 0) then
    vstl =priml(3)
    wstl =priml(4)
    bystl=0.!priml(7)
    bzstl=0.!priml(8)
!     print*,'stopped @ HLLD','denl**',denl,'sL - uL',slmul,'sL - sM',slmsM,'Bxtot',bxtot
!     stop
  else
    vstl = priml(3) - bxtot*(priml(7)+B0l(2))*sMmul/denl                      !vL*      
    wstl = priml(4) - bxtot*(priml(8)+B0l(3))*sMmul/denl                      !wL*
    bystl= (priml(7)+B0l(2))*( priml(1)*slmul**2 - bxtot**2 )/denl - B0l(2)            !byL*
    bzstl= (priml(8)+B0l(3))*( priml(1)*slmul**2 - bxtot**2 )/denl - B0l(3)           !bzL*

  endif

  if(denr == 0) then
    vstr =primr(3)
    wstr =primr(4)
    bystr=0.!primr(7)
    bzstr=0.!primr(8)
!     print*,'stopped @ HLLD','denr**','sr - ur',srmur,'sr - sM',sRmsM,'Bxtot',bxtot
!     stop
  else
    vstr = primr(3) - bxtot*(primr(7)+B0r(2))*sMmuR/denr                   !vR*
    wstr = primr(4) - bxtot*(primr(8)+B0r(3))*sMmuR/denr                      !wR*
    bystr= (primr(7)+B0r(2))*( primr(1)*srmur**2 - bxtot**2 )/denr - B0r(2)
    bzstr= (primr(8)+B0r(3))*( primr(1)*srmur**2 - bxtot**2 )/denr - B0r(3)           !bzR*
  endif

  dd=sqrt(rhostl)+sqrt(rhostr)

  vstst =(sqrt(rhostl)*vstl + sqrt(rhostr)*vstr + (bystr-bystl)*signBxtot )/dd  !v** 
  wstst =(sqrt(rhostl)*wstl + sqrt(rhostr)*wstr + (bzstr-bzstl)*signBxtot )/dd  !w**

  bystst=(sqrt(rhostl)*bystr + sqrt(rhostr)*bystl +  &      !by**
          sqrt(rhostl*rhostr)*(vstr-vstl)*signBxtot )/dd                            
 
  bzstst=(sqrt(rhostl)*bzstr + sqrt(rhostr)*bzstl +  &      !bz**
          sqrt(rhostl*rhostr)*(wstr-wstl)*signBxtot )/dd

  vststdotbstst= sM*bx + vstst*bystst + wstst*bzstst        !v** dot B**

  ! UL** region -----------------------------------
  if(sM >= 0 ) then
    !
    el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &  !eL
                +0.5*(bx**2+priml(7)**2+priml(8)**2)

    vdotbl    = priml(2)*bx + priml(3)*priml(7) + priml(4)*priml(8) !vL dot BL
    vstdotbstl= sM*bx+ vstl*bystl + wstl*bzstl                      !vL* dot BL*  

    estl= ( slmul*el -pTL*priml(2) +pst*sM +bxtot*(vdotbl-vstdotbstl) )/slmsM !eL*

    eststl= estl - sqrt(rhostl)*(vstdotbstl-vststdotbstst)*signBxtot !eL**
    
    ff(1) = rhostl*sM
    ff(2) = rhostl*sM**2+pst-bxtot*bx-bx*B0l(1)
    ff(3) = rhostl*sM*vstst-bxtot*bystst-bx*B0l(2)
    ff(4) = rhostl*sM*wstst-bxtot*bzstst-bx*B0l(3)
    ff(5) = sM*(eststl+pst)-bxtot*(vststdotbstst)
    ff(6) = 0.
    ff(7) = (bystst+B0l(2))*sM-bxtot*vstst
    ff(8) = (bzstst+B0l(3))*sM-bxtot*wstst

#ifdef PASSIVES
    ff(neqdyn+1:neq) = sM*priml(neqdyn+1:neq)*slmul/slmsM
#endif

  return
  endif

  ! UR** region -----------------------------------
  if(sM <= 0 ) then

    er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &  !eR
                +0.5*(bx**2+primr(7)**2+primr(8)**2)

    vdotbr    = primr(2)*bx + primr(3)*primr(7) + primr(4)*primr(8) !vR dot BR
    vstdotbstr= sM*bx+ vstr*bystr + wstr*bzstr                      !vR* dot BR*  

    estr= ( srmur*er -pTR*primr(2) +pst*sM +bxtot*(vdotbr-vstdotbstr) )/srmsM !eR*

    eststr= estr + sqrt(rhostr)*(vstdotbstr-vststdotbstst)*signBxtot !eR**
    
    ff(1) = rhostr*sM
    ff(2) = rhostr*sM**2+pst-bxtot*bx-bx*B0r(1)
    ff(3) = rhostr*sM*vstst-bxtot*bystst-bx*B0r(2)
    ff(4) = rhostr*sM*wstst-bxtot*bzstst-bx*B0r(3)
    ff(5) = sM*(eststr+pst)-bxtot*(vststdotbstst)
    ff(6) = 0.
    ff(7) = (bystst+B0r(2))*sM-bxtot*vstst
    ff(8) = (bzstst+B0r(3))*sM-bxtot*wstst

#ifdef PASSIVES
    ff(neqdyn+1:neq) = sM*primr(neqdyn+1:neq)*srmur/srmsM
#endif

  return
  endif

  print'(a,5es12.3)', 'Error in HLLD routine', sM,sl,sr, csl, csr
  stop

end subroutine prim2fhlld

!=======================================================================

!> @brief Calculates HLLD fluxes from the primitive variables 
!!   on all the domain
!> @details Calculates HLLD fluxes from the primitive variables 
!!   on all the domain
!> @param integer [in] choice : 1, uses primit for the 1st half of timestep
!! (first order)
!!                  @n 2 uses primit for second order timestep

subroutine hlldfluxes(choice)

  use parameters, only : neq, nx, ny, nz
  use globals, only : primit, f, g, h, B0
  use hydro_core, only : swapy, swapz, swapy_bsplit, swapz_bsplit, limiter
  implicit none
  integer, intent(in) :: choice
  integer :: i, j, k
  real, dimension(neq) :: priml, primr, primll, primrr, ff, uu
  real, dimension(3)   :: B0l, B0r, B0ll, B0rr

  select case(choice)

  case(1)        ! 1st half timestep
     !
     do k=0,nz
        do j=0,ny
           do i=0,nx

              !------- x direction -------------------------------------
              priml(:)=primit(:,i  ,j ,k )
              primr(:)=primit(:,i+1,j ,k )
              B0l(:)=B0(:,i ,j ,k)
              B0r(:)=B0(:,i+1 ,j ,k)
              !
              call prim2fhlld(priml,primr,B0l,B0r,ff)
              f(:,i,j,k)=ff(:)

              !------- y direction -------------------------------------
              priml(:)=primit(:,i ,j  ,k )
              primr(:)=primit(:,i, j+1,k )
              B0l(:)=B0(:,i ,j   ,k)
              B0r(:)=B0(:,i ,j+1 ,k)
              
              call swapy(priml,neq)          !swaps primL for L state
              call swapy(primr,neq)          !swaps primR for R state 
              call swapy_bsplit(B0l,3)
              call swapy_bsplit(B0r,3)
              !
              call prim2fhlld(priml,primr,B0l,B0r,ff)  !gets fluxes (swapped)
              call swapy(ff,neq)             !swaps back the fluxes
              g(:,i,j,k)=ff(:)

              !------- z direction -------------------------------------
              priml(:)=primit(:,i ,j ,k  )
              primr(:)=primit(:,i, j, k+1)
              B0l(:)=B0(:,i ,j ,k)
              B0r(:)=B0(:,i ,j ,k+1)
          
              call swapz(priml,neq)
              call swapz(primr,neq)
              call swapz_bsplit(B0l,3)
              call swapz_bsplit(B0r,3) 

              call prim2fhlld(priml,primr,B0l,B0r,ff)
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
              
              B0l (:)=B0(:,i   ,j ,k)
              B0r (:)=B0(:,i+1 ,j ,k)
              B0ll(:)=B0(:,i-1 ,j ,k)
              B0rr(:)=B0(:,i+2 ,j ,k)
           
              call limiter(primll,priml,primr,primrr,neq)
              call limiter(B0ll,B0l,B0r,B0rr,3)

              call prim2fhlld(priml,primr,B0l,B0r,ff)
              f(:,i,j,k)=ff(:)

              !------- y direction ------------------------------------
              priml (:)=primit(:,i,j  ,k )
              primr (:)=primit(:,i,j+1,k )
              primll(:)=primit(:,i,j-1,k )
              primrr(:)=primit(:,i,j+2,k )
              
              B0l (:)=B0(:,i,j  ,k )
              B0r (:)=B0(:,i,j+1,k )
              B0ll(:)=B0(:,i,j-1,k )
              B0rr(:)=B0(:,i,j+2,k )
              
              call swapy(priml,neq)
              call swapy(primr,neq)
              call swapy(primll,neq)
              call swapy(primrr,neq)
              
              call limiter(primll,priml,primr,primrr,neq)
              call limiter(B0ll,B0l,B0r,B0rr,3)
 
              call swapy_bsplit(B0l,3)
              call swapy_bsplit(B0r,3)


              call prim2fhlld(priml,primr,B0l,B0r,ff)
              call swapy(ff,neq)
              g(:,i,j,k)=ff(:)

              !------- z direction ------------------------------------
              priml (:)=primit(:,i,j,k  )
              primr (:)=primit(:,i,j,k+1)
              primll(:)=primit(:,i,j,k-1)
              primrr(:)=primit(:,i,j,k+2)
              
              B0l (:)=B0(:,i,j,k  )
              B0r (:)=B0(:,i,j,k+1)
              B0ll(:)=B0(:,i,j,k-1)
              B0rr(:)=B0(:,i,j,k+2)
              
              call swapz(priml,neq)
              call swapz(primr,neq)
              call swapz(primll,neq)
              call swapz(primrr,neq)
              
              call limiter(primll,priml,primr,primrr,neq)
              call limiter(B0ll,B0l,B0r,B0rr,3)
              
              call swapz_bsplit(B0l,3)
              call swapz_bsplit(B0r,3)
 

              call prim2fhlld(priml,primr,B0l,B0r,ff)
              call swapz(ff,neq)
              h(:,i,j,k)=ff(:)

           end do
        end do
     end do

  end select

end subroutine hlldfluxes

#endif

!=======================================================================

#endif

end module hlld
