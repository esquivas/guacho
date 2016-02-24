!=======================================================================
!> @file hydro_core.f90
!> @brief Hydrodynamical and Magnetohidrodynamocal bacic module
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

!> @brief Basic hydro (and MHD) subroutines utilities
!> @details This module contains subroutines and utilities that are the
!> core of the hydro (and MHD) that are common to most implementations 
!> and will be used for the different specific solvers

module hydro_core

#ifdef EOS_CHEM
  use network,  only : n_spec 
#endif
  implicit none

contains

!> @brief Computes the primitive variables and temperature from conserved
!!  variables on a single cell
!> @details Computes the primitive variables and temperature from conserved
!!  variables on a single cell
!> @param real [in] uu(neq) : conserved variables in one cell
!> @param real [out] prim(neq) : primitives in one cell
!> @param real [out] T : Temperature [K]

subroutine u2prim(uu, prim, T)

  use parameters, only : neq, neqdyn, Tempsc, vsc2, cv
  implicit none
  real,    intent(in),  dimension(neq)  :: uu
  real,    intent(out), dimension(neq)  :: prim
  real,    intent(out)                  :: T
  real :: r
#if defined(EOS_H_RATE) || defined(EOS_CHEM)
  real :: dentot
#endif

  r=max(uu(1),1e-15)
  prim(1)=r
  prim(2)=uu(2)/r 
  prim(3)=uu(3)/r
  prim(4)=uu(4)/r
  
#ifdef MHD
prim(5)=( uu(5)-0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2)   & 
               -0.5*  (  uu(6)**2+  uu(7)**2  +uu(8)**2) ) /cv  
#else
  prim(5)=( uu(5)-0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2) ) /cv
#endif
  
  prim(5)=max(prim(5),1e-16)
  
#if defined(PMHD) || defined(MHD) 
  prim(6:8) = uu(6:8)
#endif 

#ifdef PASSIVES
  prim(neqdyn+1:neq) = max( uu(neqdyn+1:neq), 0.)
#endif
  
  !   Temperature calculation

#ifdef EOS_ADIABATIC
  T=(prim(5)/r)*Tempsc
#endif

#ifdef EOS_SINGLE_SPECIE
  ! assumes it is fully ionized
  r=max(r,1e-15)
  T=max(1.,(prim(5)/r)*Tempsc)
  prim(5)=r*T/Tempsc
#endif

#ifdef EOS_H_RATE
  dentot=(2.*r-prim(neqdyn+1))
  dentot=max(dentot,1e-15)
  T=max(1.,(prim(5)/dentot)*Tempsc)
  prim(5)=dentot*T/Tempsc
#endif

#ifdef EOS_CHEM
  !  Assumes that rho scaling is mu*mh
  dentot= sum(prim( neqdyn+1 : neqdyn+n_spec ) )
  T=max(1.,(prim(5)/dentot)*Tempsc)
  prim(5) = dentot * T /Tempsc
#endif

end subroutine u2prim  

!=======================================================================

!> @brief Updated the primitives, using the conserved variables in the
!! entire domain
!> @details Updated the primitives, using the conserved variables in the
!! entire domain
!> @param real [in] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [out] prim(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! primitive variables

subroutine calcprim(u,primit)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
#ifdef THERMAL_COND
  use globals, only : Temp
#endif
  implicit none
  real,intent(in), dimension(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :: u
  real,intent(out),dimension(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :: primit
#ifndef THERMAL_COND
  real                 :: T
#endif
  integer :: i,j,k
  !
  do k=nzmin,nzmax
     do j=nymin,nymax
        do i=nxmin,nxmax
           
#ifdef THERMAL_COND
           call u2prim(u(:,i,j,k),primit(:,i,j,k),Temp(i,j,k) )
#else
           call u2prim(u(:,i,j,k),primit(:,i,j,k),T)
#endif
           
        end do
     end do
  end do

end subroutine calcprim

!=======================================================================

!> @brief Computes the conserved conserved variables from the
!! primitives in a single cell
!> @details Computes the conserved conserved variables from the
!! primitives in a single cell
!> @param real [in] prim(neq) : primitives in one cell
!> @param real [out] uu(neq) : conserved varibles in one cell

subroutine prim2u(prim,uu)

  use parameters
  implicit none
  real, dimension(neq), intent(in)  :: prim
  real, dimension(neq), intent(out) :: uu
  !
  uu(1) = prim(1)
  uu(2) = prim(1)*prim(2)
  uu(3) = prim(1)*prim(3)
  uu(4) = prim(1)*prim(4)

#ifdef MHD
  !   kinetic+thermal+magnetic energies
  uu(5) = 0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)+cv*prim(5) &
                 +0.5*(prim(6)**2+prim(7)**2+prim(8)**2)
#else
  !   kinetic+thermal energies
  uu(5) = 0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)+cv*prim(5)
#endif

#if defined(PMHD) || defined(MHD)
  uu(6:8)=prim(6:8)
#endif

#ifdef PASSIVES
  uu(neqdyn+1:neq) = prim(neqdyn+1:neq)
#endif

end subroutine prim2u

!=======================================================================

!> @brief Computes the Euler Fluxes in one cell
!> @details Computes the Euler Fluxes in one cell, using the primitices
!! @n It returns the flux in the x direction (i.e. F), the y and z fluxes
!! can be obtained swaping the respective entries (see swapy and swapz 
!!	subroutines)
!> @param real [in] prim(neq) : primitives in one cell
!> @param real [out] ff(neq) : Euler Fluxes (x direction)

subroutine prim2f(prim,ff)
  use parameters, only : neq, neqdyn, cv
  implicit none
  real,    dimension(neq), intent(in)  :: prim
  real,    dimension(neq), intent(out) :: ff
  real :: etot

  !  If MHD (active) not defined
#ifdef MHD
  ! MHD
  etot= 0.5*( prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)    &
                     + prim(6)**2+prim(7)**2+prim(8)**2  )  &
                     + cv*prim(5)
  
  ff(1) = prim(1)*prim(2)
  ff(2) = prim(1)*prim(2)*prim(2)+prim(5)+0.5*(prim(7)**2+prim(8)**2-prim(6)**2)
  ff(3) = prim(1)*prim(2)*prim(3)-prim(6)*prim(7)
  ff(4) = prim(1)*prim(2)*prim(4)-prim(6)*prim(8)
  ff(5) = prim(2)*(etot+prim(5)+0.5*(prim(6)**2+prim(7)**2+prim(8)**2) ) &
         -prim(6)*(prim(2)*prim(6)+prim(3)*prim(7)+prim(4)*prim(8))
#else
  ! HD or PMHD
  etot= 0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)+cv*prim(5)
  
  ff(1) = prim(1)*prim(2)
  ff(2) = prim(1)*prim(2)*prim(2)+prim(5)
  ff(3) = prim(1)*prim(2)*prim(3)
  ff(4) = prim(1)*prim(2)*prim(4)
  ff(5) = prim(2)*(etot+prim(5))
#endif
  
#if defined(PMHD) || defined(MHD)
  ff(6)=0.0
  ff(7)=prim(2)*prim(7)-prim(6)*prim(3)
  ff(8)=prim(2)*prim(8)-prim(6)*prim(4)
#endif

#ifdef PASSIVES
  ff(neqdyn+1:neq) = prim(neqdyn+1:neq)*prim(2)
#endif

end subroutine prim2f

!=======================================================================

!> @brief Swaps the x and y components in a cell.
!> @details Swaps the x and y components in a cell.
!> @param real [inout] var(neq) : variable to be swapped
!> @param real [in] neq : number of equations in the code

subroutine swapy(var,neq)

  implicit none
  real, intent(inout), dimension(neq) :: var
  integer, intent(in) :: neq
  real :: aux
  
  aux=var(2)
  var(2)=var(3)
  var(3)=aux
  
#if defined(PMHD) || defined(MHD)
  aux=var(6)
  var(6)=var(7)
  var(7)=aux
#endif 

end subroutine swapy

!=======================================================================

!> @brief Swaps the x and z components in a cell.
!> @details Swaps the x and z components in a cell.
!> @param real [inout] var(neq) : variable to be swapped
!> @param real [in] neq : number of equations in the code

subroutine swapz(var,neq)
  implicit none
  real, intent(inout), dimension(neq) :: var
  integer, intent(in) :: neq
  real :: aux
  
  aux=var(2)
  var(2)=var(4)
  var(4)=aux
  
#if defined(PMHD) || defined(MHD)
  aux=var(6)
  var(6)=var(8)
  var(8)=aux
#endif

end subroutine swapz

!=======================================================================

!> @brief Computes the sound speed
!> @details Computes the sound speed
!> @param real [in] p : value of pressure
!> @param real [in] d : value of density
!> @param real [out] cs : sound speed

subroutine csound(p,d,cs)

  use parameters, only : gamma
  implicit none
  real, intent(in) :: p, d
  real, intent(out) ::cs
  
  cs=sqrt(gamma*p/d)

end subroutine csound

!=======================================================================

#if defined(PMHD) || defined(MHD) 

!> @brief Computes the fast magnetosonic speeds  in the 3 coordinates
!> @details Computes the fast magnetosonic speeds  in the 3 coordinates
!> @param real [in] p  : value of pressure
!> @param real [in] d  : value of density
!> @param real [in] Bx : value of the x component of the magnetic field
!> @param real [in] By : value of the y component of the magnetic field
!> @param real [in] Bz : value of the z component of the magnetic field
!> @param real [out] csx : fast magnetisonic speed in x
!> @param real [out] csy : fast magnetisonic speed in y
!> @param real [out] csz : fast magnetisonic speed in z

subroutine cfast(p,d,bx,by,bz,cfx,cfy,cfz)

  use parameters, only : gamma
  implicit none
  real, intent(in) :: p, d, bx, by, bz
  real, intent(out) ::cfx,cfy,cfz
  real :: b2
  
	b2=bx*bx+by*by+bz*bz
  cfx=sqrt(0.5*((gamma*p+b2)/d+sqrt((gamma*p+b2)/d)**2-4.*gamma*p*bx*bx/d**2))
  cfy=sqrt(0.5*((gamma*p+b2)/d+sqrt((gamma*p+b2)/d)**2-4.*gamma*p*by*by/d**2))
  cfz=sqrt(0.5*((gamma*p+b2)/d+sqrt((gamma*p+b2)/d)**2-4.*gamma*p*bz*bz/d**2))                
  
end subroutine cfast

#endif

!=======================================================================

#if defined(PMHD) || defined(MHD) 

!> @brief Computes the fast magnetosonic speed in the x direction
!> @details Computes the fast magnetosonic speed in the x direction
!> @param real [in] prim(neq) : vector with the primitives in one cell

subroutine cfastX(prim,cfX)
  
  use parameters, only : neq, gamma
  implicit none
  real, intent(in) :: prim(neq)
  real, intent(out) ::cfX
  real :: b2, cs2va2
  
  b2=prim(6)**2+prim(7)**2+prim(8)**2
  cs2va2 = (gamma*prim(5)+b2)/prim(1)   ! cs^2 + ca^2

  cfx=sqrt(0.5*(cs2va2+sqrt(cs2va2**2-4.*gamma*prim(5)*prim(6)**2/prim(1)/prim(1) ) ) )
 
  end subroutine cfastX

#endif

!=======================================================================

!> @brief Otains the timestep allowed by the CFL condition in the entire
!> @details Otains the timestep allowed by the CFL condition in the entire
!> domain using the global primitives, and sets logical variable to dump 
!> output
!> @param integer [in] current_iter : Current iteration, it starts with a small
!! but increasing CFL in the first N_trans iterarions
!> @param integer [in] n_iter : Number of iterations to go from a small CFL to
!! the final CFL (in parameters.f90)
!> @param real [in] current_time : Current (global) simulation time
!> @param real [in] tprint : time for the next programed disk dump
!> @param real [out] : @f$ \Delta t@f$ allowed by the CFL condition
!> @param logical [out] dump_flag : Flag to write to disk

subroutine get_timestep(current_iter, n_iter, current_time, tprint, dt, dump_flag)

  use parameters, only : nx, ny, nz, cfl, mpi_real_kind
  use globals, only : primit, dx, dy, dz
  implicit none
#ifdef MPIP
  include "mpif.h"
#endif
  integer, intent(in)  :: current_iter, n_iter
  real,    intent(in)  :: current_time, tprint
  real,    intent(out) :: dt
  logical, intent(out) :: dump_flag
  real              :: dtp
#ifdef MHD
  real              :: cx, cy, cz
#else
  real              :: c
#endif
  integer :: i, j, k, err
  
  dtp=1.e30
  
  do k=1,nz
    do j=1,ny
      do i=1,nx
      
#ifdef MHD
        call cfast(primit(5,i,j,k),primit(1,i,j,k),&
            primit(6,i,j,k), primit(7,i,j,k), primit(8,i,j,k), &
            cx,cy,cz)
        dtp=min(dtp,dx/(abs(primit(2,i,j,k))+cx))  
        dtp=min(dtp,dy/(abs(primit(3,i,j,k))+cy))
        dtp=min(dtp,dz/(abs(primit(4,i,j,k))+cz))
#else
        call csound(primit(5,i,j,k),primit(1,i,j,k),c)
        dtp=min(dtp,dx/(abs(primit(2,i,j,k))+c))  
        dtp=min(dtp,dy/(abs(primit(3,i,j,k))+c))
        dtp=min(dtp,dz/(abs(primit(4,i,j,k))+c))
#endif
      end do
    end do
  end do

  if (current_iter <= n_iter ) then
    !   boots up simulation with small CFL
    dtp = cfl*(2.**(-real(n_iter+1-current_iter)))*dtp
  else
    dtp=cfl*dtp
  end if

#ifdef MPIP
  call mpi_allreduce(dtp, dt, 1, mpi_real_kind, mpi_min, mpi_comm_world,err)
#else
  dt=dtp
#endif

  !  Adjust dt so t+dt doesnt overshoot tprint
  if(( current_time + dt )>=tprint) then
    dt= tprint-current_time
    dump_flag=.true.
  endif


end subroutine get_timestep

!=======================================================================

!> @brief Performs a linear reconstruction of the primitive variables
!> @details returns a linear reconstruction of the variables at the
!> interface beteen the primitives PLL, PL, PR, PRR
!! @n The reconstruction is made with a slope limiter chosen at
!! compilation time (i.e. set on the Makefile)
!> @param real [in]    : primitives at the left of the left state
!> @param real [inout] : primitives at the left state
!> @param real [inout] : primitives at the right state
!> @param real [in]    : primitives at the right of the right state
!> @param real [in]    : number of equations

subroutine limiter(PLL,PL,PR,PRR,neq)

  implicit none
  real, dimension(neq), intent(inout) :: pl,  pr
  real, dimension(neq), intent(in)    :: pll, prr
  integer, intent (in)  :: neq
  real :: dl, dm, dr, al, ar
  integer :: ieq
  
  do ieq=1,neq
     dl=pl(ieq)-pll(ieq)
     dm=pr(ieq)-pl(ieq)
     dr=prr(ieq)-pr(ieq)
     al=average(dl,dm)
     ar=average(dm,dr)
     pl(ieq)=pl(ieq)+al*0.5
     pr(ieq)=pr(ieq)-ar*0.5
  end do
  
contains

  real function average(a,b)
    implicit none
    real, intent(in)    :: a, b

#if LIMITER==-1

    !   no average (reduces to 1st order)
    average=0.
#endif

#if LIMITER==0
    !   no limiter
    average=0.5*(a+b)
#endif

#if LIMITER==1
    !   Minmod - most diffusive
    real :: s
    s=sign(1.,a)
    average=s*max(0.,min(abs(a),s*b))
#endif

#if LIMITER==2
    !   Falle Limiter (Van Leer)
    if(a*b.le.0.) then
       average=0.
    else
       average=a*b*(a+b)/(a**2+b**2)
    end if
#endif

#if LIMITER==3
    !   Van Albada
    real, parameter :: delta=1.e-7
    average=(a*(b*b+delta)+b*(a*a+delta))/(a*a+b*b+delta)
#endif

#if LIMITER==4
    !   UMIST Limiter - less diffusive
    real :: s, c, d
    s=sign(1.,a)
    c=0.25*a+0.75*b
    d=0.75*a+0.25*b
    average=min(2.*abS(a),2.*s*b,s*c,s*d)
    average=s*max(0.,average)
#endif

#if LIMITER==5
    !    Woodward Limiter (o MC-limiter; monotonized central difference)
    real :: s, c
    s=sign(1.,a)
    c=0.5*(a+b)
    average=min(2.*abs(a),2.*s*b,s*c)
    average=s*max(0.,average)
#endif
#if LIMITER==6
    !   superbee Limiter (tends to flatten circular waves)
    real :: s, av1, av2
    s=sign(1.,b)
    av1=min(2.*abs(b),s*a)
    av2=min(abs(b),2.*s*a)
    average=s*max(0.,av1,av2)
#endif

  end function average

end subroutine limiter

!=======================================================================

end module hydro_core

!=======================================================================
