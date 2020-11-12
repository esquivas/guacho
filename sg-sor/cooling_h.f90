!=======================================================================
!> @file cooling_h.f90
!> @brief Cooling with hydrogen rate parametrized cooling
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

!> @brief Cooling with parametrized cooling and H rate equation
!> @details Cooling with parametrized cooling and H rate equation

module cooling_H

#ifdef PASSIVES

  implicit none

contains

  !=======================================================================
  !> @brief High level wrapper to apply cooling
  !> @details High level wrapper to apply cooling
  !! @n  parametrized cooling curve, uses the ionization state of
  !! hydrogen and ties the O I and II to it
  subroutine coolingh()
    use parameters, only : neq, nx, ny, nz, tsc, dif_rad, charge_exchange
    use globals, only : u, primit, dt_CFL
    use difrad, only : ph, phCold, phHot
    implicit none
    real    :: dt_seconds
    integer :: i,j,k

    dt_seconds = dt_CFL*tsc

    do k=1,nz
      do j=1,ny
        do i=1,nx

          if (dif_rad) then
            if (.not.charge_exchange) then
              call cooling_h_neq( primit(:,i,j,k), u(:,i,j,k), dt_seconds,     &
                                  ph(i,j,k) )
            else
              call cooling_h_neq( primit(:,i,j,k), u(:,i,j,k), dt_seconds,     &
                                  phCold(i,j,k)+phHot(i,j,k))
            end if
          else
            call cooling_h_neq( primit(:,i,j,k), u(:,i,j,k), dt_seconds, 0. )
          end if

        end do
      end do
    end do

  end subroutine coolingh

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
    betah=1.133e-24/sqrt(a)*(-0.0713+0.5*log(a)+0.640*a**(-0.33333))

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
    real (kind=8), intent(in) :: X1,X2,DEN,DH0,TE0
    real, parameter :: XION=2.179e-11, XH=0.9,XO=1.e-3
    real, parameter :: C0=0.5732,C1=1.8288e-5,C2=-1.15822e-10,C3=9.4288e-16
    real, parameter :: D0=0.5856,D1=1.55083e-5,D2=-9.669e-12, D3=5.716e-19
    real, parameter :: ENK=118409.,EN=1.634E-11
    real (kind=8) :: SHP
    real (kind=8) :: TE, DH,DHP,DE,DOI,DOII,OMEGA,OMEGAL,OMEGAH,FRAC,QLA
    real (kind=8) :: ECOLL,CION,EION,EREC,TM,T2,EOI,EOII,EQUIL,FR,EX2,TANH
    real (kind=8) :: BETAF, HIICOOL

    Te=MAX(Te0,10.)
    DH=DEN
    DHP=(1.-X1)*DH
    DE=DHP+1.E-4*DH
    DOI=XO*DH0
    DOII=XO*DHP

    SHP=-DH*(X1-X2)/DT  !  ?```
    if(TE <= 1e4 ) then
      ALOSS = 0.
      return
    end if

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
  !> @brief
  !> @details
  !> @param real [in] uu(neq) : primitive variablas in one cell
  !> @param real [in] uu(neq) : conserved variablas in one cell
  !> @param real [in] dt      : timestep (seconds)
  !> @param real [in] radphi  : photoionizing rate
  subroutine  cooling_h_neq(pp, uu, dt, radphi)
    use parameters, only : neqdyn, dif_rad, mhd, cv, neq, Tempsc, eq_of_state
    use constants
    use hydro_core, only : u2prim
    implicit none
    real, intent(inout) :: uu(neq), pp(neq)
    real, intent(in)    :: dt, radphi
    real                :: T, ch_factor
    real(kind = 8)      :: y0, y1, dh, dh0, gain, Tprime, al, ce, T1

    y0 =  real( pp(neqdyn+1)/pp(1), 8 )  !# neutral H fraction (t0)
    y1  = real( uu(neqdyn+1)/uu(1), 8 )  !# neutral H fraction (t0+dt)
    dh  = real( pp(1)             , 8)   !# total NH
    dh0 = real( pp(neqdyn+1)      , 8)   !# neutrals density

    call u2prim(uu,pp,T)               !# get temperature

    !  get the energy losses
    al=ALOSS(y0,y1,dt,dh,dh0,real(T,8))/dh**2

    if (dif_rad) then
      gain=real(radphi,8)*dh0*Kb*1.E4 !3.14d5
      Tprime=max( gain*real(T,8)/(dh**2*al),1000.)
    else
      tprime=10.
    end if
    !tprime=1000.

    ce=(2.*dh*al)/(3.*Kb*real(T,8))
    T1=Tprime+(T-Tprime)*exp(-ce*dt) !# new temperature

    T1=max(T1,0.5*real(T,8) )
    T1=min(T1,2.*real(T,8) )
    !  t1=max(t1,tprime)

    ch_factor = real(T1)/T

    !  update pressure
    pp(5) = pp(5) * ch_factor

    !  set pressure floor "a la mala"
    !pp(5) = max(pp(5),(pp(1)+pp(11))*10000./Tempsc)

    !  update total energy density
    if (mhd) then
#ifdef BFIELD

      uu(5) = cv*pp(5) + 0.5*pp(1)*(pp(2)**2+pp(3)**2+pp(4)**2)                &
                     + 0.5*      (pp(6)**2+pp(7)**2+pp(8)**2)

#endif
    else
      uu(5) = cv*pp(5) + 0.5*pp(1)*(pp(2)**2+pp(3)**2+pp(4)**2)
    end if

  end subroutine cooling_h_neq

#endif

end module cooling_H

!======================================================================
