!=======================================================================
!> @file cooling_h.f90
!> @brief Cooling with hydrogen rate parametrized cooling
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

!> @brief Cooling with parametrized cooling and H rate equation
!> @details Cooling with parametrized cooling and H rate equation

module cooling_H

#ifdef PASSIVES

  implicit none

contains

!> @brief High level wrapper to apply cooling
!> @details High level wrapper to apply cooling
!! @n  parametrized cooling curve, uses the ionization state of
!! hydrogen and ties the O I and II to it
!> @param real [in] dt : timestep in seconds

subroutine coolingh()

  use parameters, only : neq, nx, ny, nz, tsc, charge_exchange
  use globals, only : u, coords, dt_CFL
  use difrad, only  : phCold, phHot, ph

  implicit none
  real :: dt_seconds
  integer :: i,j,k

  dt_seconds = dt_CFL*tsc

  do i=1,nx
    do j=1,ny
      do k=1,nz
        if (charge_exchange) then
          call atomic(dt_seconds,u(:,i,j,k),1.,phCold(i,j,k), phHot(i,j,k))
        else
          call atomic(dt,u(:,i,j,k),1.,ph(i,j,k), 1.)
        end if
      end do
    end do
  end do

end subroutine coolingh

!======================================================================

!> @brief calculates the recombination rate (case B)
!> @details calculates the recombination rate (case B)
!> @param real8 [in] T : Temperature K

function alpha(T)

  implicit none

  real (kind=8) :: alpha
  real (kind=8), intent(in) :: T

  alpha=2.55d-13*(1.d4/T)**0.79

end function alpha

!======================================================================

!> @brief calculates the recombination rate to level 1
!> @details calculates the recombination rate to level 1
!> @param real8 [in] T : Temperature K

function alpha1(T)
  implicit none

  real (kind=8) :: alpha1
  real (kind=8), intent(in) :: T

  alpha1=1.57d-13*(1.d4/T)**0.52

end function alpha1

!======================================================================

!> @brief calculates the collisional ionization rate
!> @details calculates the collisional ionization rate
!> @param real8[in] T : Temperature K

function colf(T)

  implicit none

  real (kind=8) :: colf
  real (kind=8), intent(in) :: T

  colf=5.83d-11*sqrt(T)*exp(-157828./T)

end function colf

!======================================================================

!> @brief betaH(T)
!> @details @f$ \beta_H(T) @f$
!> @param real 8[in] T : Temperature K

function betah(T)

  implicit none

  real (kind=8) ::  betah
  real (kind=8), intent(in) ::  T
  real (kind=8)             ::  a

  a=157890./T
  betah=1.133D-24/sqrt(a)*(-0.0713+0.5*log(a)+0.640*a**(-0.33333))

end function betah

!======================================================================


!> @brief Non equilibrium cooling
!> @details   Non-equilibrium energy loss for low temperatures
!!     considering the collisional excitation of [O I] and
!!   [O II] lines and radiative recombination of H. This
!!   cooling rate is multiplied by a factor of 7.033 so
!!   that it has the same value as the "coronal equilibrium"
!!   cooling rate at a temperature of 44770 K (at temperatures
!!   higher than this value, the equilibrium cooling rate is
!!   used). The collisional ionization of H and excitation
!!   of Lyman-alpha are computed separately, and added to
!!   the cooling rate.
!> @param real8 [in] x1  : initial H ionization fraction
!> @param real8 [in] x2  : final H ionization fraction
!> @param real [in] dt  : timestep
!> @param real8 [in] den : total density of hydrogen
!> @param real8 [in] dh0 : density of neutral hydrogen
!> @param real8 [in] Te0 : Temperature

FUNCTION ALOSS(X1,X2,DT,DEN,DH0,TE0)

  implicit none

  real (kind=8) :: ALOSS
  real, intent(in)          :: DT
  !real (kind=8), intent(in) :: X1,X2,DEN,DH0,TE0
  real (kind=8), intent(in) :: X1,X2,DEN,DH0,TE0
  real, parameter :: XION=2.179e-11, XH=0.9,XO=1.e-3
  real, parameter :: C0=0.5732,C1=1.8288e-5,C2=-1.15822e-10,C3=9.4288e-16
  real, parameter :: D0=0.5856,D1=1.55083e-5,D2=-9.669e-12, D3=5.716e-19
  real, parameter :: ENK=118409.,EN=1.634E-11
  real (kind=8) :: SHP
  !real (kind=8) :: TE, DH,DHP,DE,DOI,DOII,OMEGA,OMEGAL,OMEGAH,FRAC,QLA
  !real (kind=8) :: ECOLL,CION,EION,EREC,TM,T2,EOI,EOII,EQUIL,FR,EX2,TANH
  !real (kind=8) :: BETAH, BETAF, HIICOOL
  real (kind=8) :: TE, DH,DHP,DE,DOI,DOII,OMEGA,OMEGAL,OMEGAH,FRAC,QLA
  real (kind=8) :: ECOLL,CION,EION,EREC,TM,T2,EOI,EOII,EQUIL,FR,EX2,TANH
  real (kind=8) :: BETAF, HIICOOL!, BETAH

  Te=MAX(Te0,1000.)
  DH=DEN
  DHP=(1.-X1)*DH
  DE=DHP+1.E-4*DH
  DOI=XO*DH0
  DOII=XO*DHP

  SHP=-DH*(X1-X2)/DT  !  ?

  !   Collisionally excited Lyman alpha

  IF(TE.LE.55000.) OMEGA=C0+TE*(C1+TE*(C2+TE*C3))
  IF(TE.GE.72000.) OMEGA=D0+TE*(D1+TE*(D2+TE*D3))
  IF(TE.GT.55000..AND.TE.LT.72000.) THEN
     OMEGAL=C0+TE*(C1+TE*(C2+TE*C3))
     OMEGAH=D0+TE*(D1+TE*(D2+TE*D3))
     FRAC=(TE-55000.)/17000.
     OMEGA=(1.-FRAC)*OMEGAL+FRAC*OMEGAH
  END IF
  QLA=8.6287E-6/(2.*SQRT(TE))*OMEGA*EXP(-ENK/TE)
  ECOLL=DE*DH0*QLA*EN
  ECOLL=MAX(ECOLL,0.)

  !   Hydrogen recombination and collisional ionization

  CION=5.834E-11*SQRT(TE)*EXP(-1.579E5/TE)
  !     AREC=2.61E-10/TE**0.7
  !     EREC=DE*DHP*(BETAH(TE)+AREC*XION)
  !     EION=SHP*XION
  EION=DE*DH0*CION*XION
  !     EION=0.
  EREC=DE*DHP*(BETAH(TE))
  EREC=MAX(EREC,0.)

  !   [O I] and [O II] coll. excited lines

  TM=1./TE
  T2=TM*TM
  EOI=DE*DOI*10.**(1381465*T2-12328.69*TM-19.82621)
  EOII=DE*DOII*10.**(-2061075.*T2-14596.24*TM-19.01402)
  EOI=MAX(EOI,0.)
  EOII=MAX(EOII,0.)

  !   free-free cooling

  BETAF=1.3*1.42E-27*TE**0.5
  HIICOOL=DE*DHP*BETAF

  !   Equilibrium cooling (for high Te)

  EQUIL=(1.0455E-18/TE**0.63)*(1.-EXP(-(TE*1.E-5)**1.63))*DE*DEN+HIICOOL

  !   switch between non-equil. and equil. ionization cooling

  IF(TE.LE.44770.) FR=0.
  IF(TE.GE.54770.) FR=1.
  IF(TE.GT.44770..AND.TE.LT.54770.) THEN
     EX2=EXP(-2.*(TE-49770.)/500.)
     TANH=(1.-EX2)/(1.+EX2)
     FR=0.5*(1.+TANH)
  END IF

  ALOSS=ECOLL+EION+(EREC+7.033*(EOI+EOII))*(1.-FR)+EQUIL*FR
  !
END FUNCTION ALOSS

!=======================================================================

!> @brief Updates the ionization fraction and applpies cooling
!> @details Calculates the new ionization state and energy density
!!      using a time dependent ionization calculation and an
!!      approximate time dependent cooling calculation
!> @param real [in] dt      : timestep (seconds)
!> @param real [in] uu(neq) : conserved variablas in one cell
!> @param real [in] tau     : optical depth (not in use)
!> @param real [in] radphi  : photoionizing rate

subroutine atomic(dt,uu,tau,radphi1, radphi2)

  use parameters
  use globals, only : dx
  use hydro_core, only : u2prim
  implicit none

  integer                          :: nn, ll, lp
  real, intent(in)                 :: dt, tau
  real                             :: radphi, radphi1, radphi2
  real, intent(out),dimension(neq) :: uu
  real, dimension(neq)             :: prim
  real                             :: T, densn,f
  !real (kind=16) :: etau, psi, dh, y0, fpn, g0, e, y1, t1,dh0, al
  !real (kind=16) :: gain, tprime, ce, ALOSS
  real (kind=8) :: etau, dh, y0, fpn, g0, e, y1, t1,dh0, al, p1,p2,p3,p4,p5, r, dentot
  real (kind=8) :: gain, tprime, ce  !, ALOSS

  !    these need to be double precision in order for
  !      the ionization calculation to work

  !  real (kind=8) :: te, col,rec,a,b,c,d,colf,alpha
  real (kind=8) :: col,rec,a,b,c,d  !,colf,alpha
  !
  !    parameters
  !      xi - neutral carbon abundance (for non-zero electron density
  !      boltzm - Boltzmann's constant
  !
  !  real (kind=8), parameter ::  xi=1.d-4,boltzm=1.3807d-16,AMASS=2.158d-24
  real (kind=8), parameter ::  xi=1.d-4,boltzm=1.3807d-16,AMASS=2.158d-24

  !  the following is not used, and kept only to avoid compiler warnings
  !ph0=0.
  !psi0=0.
  !   atenuate photoionization with optical depth (already atenuated in radif)
  etau=exp(-tau)
  !!!   radphi=radphi0*etau
  !psi=psi0*etau

  !   solve for the ionization fraction and the internal energy

  !the charge exchange use planet and stellar neutrals
  !this loop allow to update both neutrals densities and compute
  !the energy at the end with the sum of both
nn = 0
lp = neqdyn+1
densn = 0.0
radphi = radphi1 !!total ph if exchange is turn off

if (charge_exchange) then
  lp = neqdyn+3
  radphi = radphi1
  nn = 1
end if

do ll=1, 1+nn

   call u2prim(uu,prim,T)     !# temperature

   col= colf(dble(t))  !# collisional ionization rate
   rec= alpha(dble(t)) !# rad. recombination rate

   y0=dble( uu(lp)/uu(1) )     !# neutral H fraction
   dh=dble( uu(1) )                  !# H density

   if (dif_rad) then
     fpn=dble(radphi)/dh               !# ionizing flux per nucleus
   end if

   !   solve for the new neutral fraction using the analytical
   !    solution (see notes)

   a=rec+col
if (dif_rad) then
   b=-((2.+xi)*rec+(1.+xi)*col+fpn)
else
   b=-((2.+xi)*rec+(1.+xi)*col)
end if
   c=(1.+xi)*rec
   d=sqrt(b**2-4.*a*c)
   g0=(2.*a*y0+b+d)/(2.*a*y0+b-d)
   e=exp( -d*dh*dble(dt) )

   y1=(-b-d*(1.+g0*e)/(1.-g0*e))/(2.*a) !# the new neutral fraction

   y1=min(y1,0.9999)
   y1=max(y1,0.)

   !    find the new total energy using the cooling and heating rates
   !    and assuming that the cooling goes approximately linear with
   !    temperature

   dh0=dble( uu(lp) )
   al=ALOSS(y0,y1,dt,dh,dh0,dble(t))/dh**2

   !if(al.lt.0.) write(*,*) 'que paso !'
   !if(al.lt.0.) al=0.
   !  al=al*(1.-(0.5e4/max(1.e4,t))**4)
   !if(t.le.1.e4) al=al*dble((t/1.e4)**4)

  if (dif_rad) then
    gain=dble(radphi)*dh0*boltzm*1.E4 !3.14d5
    tprime=max( gain*dble(T)/(dh**2*al),1000.)
  else
    tprime=10.
  end if
   !tprime=1000.

   ce=(2.*dh*al)/(3.*boltzm*dble(T))
   t1=tprime+(t-tprime)*exp(-ce*dt) !# new temperature

   t1=max(t1,0.1*dble(t) )
   t1=min(t1,10.*dble(t) )
   !  t1=max(t1,tprime)

!#if TWOTEMP
!   t1=1.E4-9990.*y1
!#endif
   !   update the uu array (also the ion for charge exchange)
   uu(lp)=sngl(y1)*uu(1)

if (charge_exchange) then
  uu(lp-1) = uu(lp-1) + (sngl(y0-y1))*uu(1)
end if

   densn = densn + uu(lp) !used for energy

   lp = lp + 2 !!neqdyn+3 -> neqdyn+5
   radphi = radphi2 !! hot ph
enddo

uu(5) = cv*(2.*uu(1)-densn)*sngl(t1)/Tempsc        &
     +0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)

end subroutine atomic

#endif

end module cooling_H

!======================================================================
