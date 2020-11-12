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
  !> @n  parametrized cooling curve, uses the ionization state of
  !> hydrogen and ties the O I and II to it
  subroutine coolingh()

    use parameters, only : neq, nx, ny, nz, tsc, dif_rad
    use globals, only : u, coords, dt_CFL
    use difrad, only : ph

    implicit none
    real    :: dt_seconds
    integer :: i,j,k

    dt_seconds = dt_CFL*tsc

    do k=1,nz
      do j=1,ny
        do i=1,nx

          if (dif_rad) then
            call atomic(dt_seconds,u(:,i,j,k),1.,ph(i,j,k) )
          else
            call atomic(dt_seconds,u(:,i,j,k),1.,1.)
          endif

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

    alpha=2.55e-13*(1.e4/T)**0.79

  end function alpha

  !======================================================================
  !> @brief calculates the recombination rate to level 1
  !> @details calculates the recombination rate to level 1
  !> @param real8 [in] T : Temperature K
  function alpha1(T)

    implicit none
    real (kind=8) :: alpha1
    real (kind=8), intent(in) :: T

    alpha1=1.57e-13*(1.e4/T)**0.52

  end function alpha1

  !======================================================================
  !> @brief calculates the collisional ionization rate
  !> @details calculates the collisional ionization rate
  !> @param real8[in] T : Temperature K
  function colf(T)

    implicit none
    real (kind=8) :: colf
    real (kind=8), intent(in) :: T

    colf=5.83e-11*sqrt(T)*exp(-157828./T)

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
  function aloss(x1,x2,dt,den,dh0,te0)

    implicit none

    real (kind=8) :: aloss
    real, intent(in)          :: dt
    !real (kind=8), intent(in) :: x1,x2,den,dh0,te0
    real (kind=8), intent(in) :: x1,x2,den,dh0,te0
    real, parameter :: xion=2.179e-11, xh=0.9,xo=1.e-3
    real, parameter :: c0=0.5732,c1=1.8288e-5,c2=-1.15822e-10,c3=9.4288e-16
    real, parameter :: d0=0.5856,d1=1.55083e-5,d2=-9.669e-12, d3=5.716e-19
    real, parameter :: enk=118409.,en=1.634e-11
    real (kind=8) :: shp
    !real (kind=8) :: te, dh,dhp,de,doi,doii,omega,omegal,omegah,frac,qla
    !real (kind=8) :: ecoll,cion,eion,erec,tm,t2,eoi,eoii,equil,fr,ex2,tanh
    !real (kind=8) :: betah, betaf, hiicool
    real (kind=8) :: te, dh,dhp,de,doi,doii,omega,omegal,omegah,frac,qla
    real (kind=8) :: ecoll,cion,eion,erec,tm,t2,eoi,eoii,equil,fr,ex2,tanh
    real (kind=8) :: betaf, hiicool!, betah

    te=max(te0,10.)
    dh=den
    dhp=(1.-x1)*dh
    de=dhp+1.e-4*dh
    doi=xo*dh0
    doii=xo*dhp

    shp=-dh*(x1-x2)/dt  !  ?```
    if(te <= 1e4 ) then
      aloss = 0.
      return
    end if

    !   collisionally excited lyman alpha
    if(te.le.55000.) omega=c0+te*(c1+te*(c2+te*c3))
    if(te.ge.72000.) omega=d0+te*(d1+te*(d2+te*d3))
    if(te.gt.55000..and.te.lt.72000.) then
       omegal=c0+te*(c1+te*(c2+te*c3))
       omegah=d0+te*(d1+te*(d2+te*d3))
       frac=(te-55000.)/17000.
       omega=(1.-frac)*omegal+frac*omegah
    end if
    qla=8.6287e-6/(2.*sqrt(te))*omega*exp(-enk/te)
    ecoll=de*dh0*qla*en
    ecoll=max(ecoll,0.)

    !   hydrogen recombination and collisional ionization

    cion=5.834e-11*sqrt(te)*exp(-1.579e5/te)
    !     arec=2.61e-10/te**0.7
    !     erec=de*dhp*(betah(te)+arec*xion)
    !     eion=shp*xion
    eion=de*dh0*cion*xion
    !     eion=0.
    erec=de*dhp*(betah(te))
    erec=max(erec,0.)

    !   [o i] and [o ii] coll. excited lines

    tm=1./te
    t2=tm*tm
    eoi=de*doi*10.**(1381465*t2-12328.69*tm-19.82621)
    eoii=de*doii*10.**(-2061075.*t2-14596.24*tm-19.01402)
    eoi=max(eoi,0.)
    eoii=max(eoii,0.)

    !   free-free cooling

    betaf=1.3*1.42e-27*te**0.5
    hiicool=de*dhp*betaf

    !   equilibrium cooling (for high te)

    equil=(1.0455e-18/te**0.63)*(1.-exp(-(te*1.e-5)**1.63))*de*den+hiicool

    !   switch between non-equil. and equil. ionization cooling

    if(te.le.44770.) fr=0.
    if(te.ge.54770.) fr=1.
    if(te.gt.44770..and.te.lt.54770.) then
       ex2=exp(-2.*(te-49770.)/500.)
       tanh=(1.-ex2)/(1.+ex2)
       fr=0.5*(1.+tanh)
    end if

    aloss=ecoll+eion+(erec+7.033*(eoi+eoii))*(1.-fr)+equil*fr
    !
  end function aloss

  !=======================================================================
  !> @brief Updates the ionization fraction and applpies cooling
  !> @details Calculates the new ionization state and energy density
  !!      using a time dependent ionization calculation and an
  !!      approximate time dependent cooling calculation
  !> @param real [in] dt      : timestep (seconds)
  !> @param real [in] uu(neq) : conserved variablas in one cell
  !> @param real [in] tau     : optical depth (not in use)
  !> @param real [in] radphi  : photoionizing rate
  subroutine atomic(dt,uu,tau,radphi)

    use parameters
    use hydro_core, only : u2prim
    implicit none
    real, intent(in)                 :: dt, tau, radphi
    real, intent(out),dimension(neq) :: uu
    real, dimension(neq)             :: prim
    real                             :: T
    real (kind=8) :: etau, dh, y0, g0, e, y1, t1,dh0, al
    real (kind=8) :: tprime, ce  !, ALOSS
    real(kind=8) :: fpn, gain

    !    these need to be double precision in order for
    !      the ionization calculation to work
    real (kind=8) :: col,rec,a,b,c,d
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

    call u2prim(uu,prim,T)               !# temperature
    col=colf(real(t,8))                  !# collisional ionization rate
    rec=alpha(real(t,8))                 !# rad. recombination rate
    y0=real( uu(neqdyn+1)/uu(1), 8 )     !# neutral H fraction
    dh=real( uu(1), 8 )                  !# H density
    fpn=real(radphi, 8)/dh               !# ionizing flux per nucleus
    !print*,fpn
    !fpn=0.

    !    solve for the new neutral fraction using the analytical
    !    solution (see notes)

    a=rec+col

    if (dif_rad) then
      b=-((2.+xi)*rec+(1.+xi)*col+fpn)
    else
      b=-((2.+xi)*rec+(1.+xi)*col    )
    end if

    c=(1.+xi)*rec
    d=sqrt(b**2-4.*a*c)
    g0=(2.*a*y0+b+d)/(2.*a*y0+b-d)
    e=exp( -d*dh*real(dt,8) )

    y1=(-b-d*(1.+g0*e)/(1.-g0*e))/(2.*a) !# the new neutral fraction
    y1=min(y1,0.9999)
    y1=max(y1,0.)

    !    find the new total energy using the cooling and heating rates
    !    and assuming that the cooling goes approximately linear with
    !    temperature

    dh0=real( uu(neqdyn+1), 8 )
    al=ALOSS(y0,y1,dt,dh,dh0,real(t,8))/dh**2

    !if(al.lt.0.) write(*,*) 'que paso !'
    !if(al.lt.0.) al=0.
    !  al=al*(1.-(0.5e4/max(1.e4,t))**4)
    !if(t.le.1.e4) al=al*real((t/1.e4,8)**4)

    if (dif_rad) then
      gain=real(radphi,8)*dh0*boltzm*4.E4
      tprime=max( gain*real(T,8)/(dh**2*al),1000.)
    else
      tprime=10.
    end if
    !tprime=1000.

    ce=(2.*dh*al)/(3.*boltzm*real(T,8))
    t1=tprime+(t-tprime)*exp(-ce*dt) !# new temperature

    t1=max(t1,0.1*real(t,8) )
    t1=min(t1,10.*real(t,8) )
    !  t1=max(t1,tprime)

!#if TWOTEMP
!    t1=1.E4-9990.*y1
!#endif
    !   update the uu array
    uu(neqdyn+1)=real(y1)*uu(1)

    if (mhd) then
#ifdef BFIELD
      uu(5) = cv*(2.*uu(1)-uu(neqdyn+1))*real(t1)/Tempsc                       &
            + 0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)                   &
            + 0.5*        (prim(6)**2+prim(7)**2+prim(8)**2)
#endif
    else
      uu(5) = cv*(2.*uu(1)-uu(neqdyn+1))*real(t1)/Tempsc                       &
            + 0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)
    end if

  end subroutine atomic


  !======================================================================

#endif

end module cooling_H
