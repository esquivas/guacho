!=======================================================================
!> @file cooling_h.f90
!> @brief Hydrogen rate eqn.
!> @author Alejandro Esquivel & M. Schneiter
!> @date 7/Apr/2017

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

module hrate

#ifdef PASSIVES

implicit none

contains

  !======================================================================
  !> @brief High level wrapper to apply cooling
  !> @details High level wrapper to apply cooling
  !! @n  parametrized cooling curve, uses the ionization state of
  !! hydrogen and ties the O I and II to it
  subroutine update_neutral_fraction()
    use parameters, only : neq, nx, ny, nz, tsc, dif_rad, charge_exchange
    use globals, only : u, primit, coords, dt_CFL
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
              call solve_h_rate(dt_seconds,u(:,i,j,k),primit(:,i,j,k),ph(i,j,k))
            else
              call solve_h_rate(dt_seconds,u(:,i,j,k), primit(:,i,j,k),        &
                                phHot(i,j,k)+phCold(i,j,k) )
            end if
          else
            call solve_h_rate(dt_seconds,u(:,i,j,k), primit(:,i,j,k), 0.)
          end if

        end do
      end do
    end do

  end subroutine update_neutral_fraction

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
  !> @brief calculates the collisional ionization rate
  !> @details calculates the collisional ionization rate
  !> @param real8[in] T : Temperature K
  function colf(T)
    implicit none

    real (kind=8) :: colf
    real (kind=8), intent(in) :: T

    colf=5.83e-11*sqrt(T)*exp(-157828./T)

  end function colf

  !=======================================================================
  !> @brief Updates the ionization fraction using Hrate eqn.
  !> @param real [in] dt        : timestep (seconds)
  !> @param real [in] uu(neq)   : conserved variablas in one cell
  !> @param real [in] prim(neq) : primitive variablas in one cell
  !> @param real [in] radphi    : photoionizing rate
  subroutine solve_h_rate(dt,uu,prim,radphi)
    use parameters
    use hydro_core, only : u2prim
    implicit none
    real, intent(in)                   :: dt, radphi
    real, intent(inout),dimension(neq) :: uu, prim
    real                               :: T
    real (kind=8) :: dh, y0, g0, e, y1
    real (kind=8) :: fpn
    real (kind=8) :: col,rec,a,b,c,d
    !    parameters
    !      xi - neutral carbon abundance (for non-zero electron density
    real (kind=8), parameter ::  xi=1.e-4

    !   solve for the ionization fraction and the internal energy
    call u2prim(uu,prim,T)            !# temperature
    col=colf(real(t,8))               !# collisional ionization rate
    rec=alpha(real(t,8))              !# rad. recombination rate
    y0=real( uu(neqdyn+1)/uu(1), 8 )  !# neutral H fraction
    dh=real( uu(1), 8 )               !# H density
    fpn=real(radphi, 8)/dh            !# ionizing flux per nucleus

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
    y1=min(y1,1.-xi)
    y1=max(y1,0.)

    !   update the density of neutrals (conserved vars)
    uu(neqdyn+1)=real(y1)*uu(1)

  end subroutine solve_h_rate

#endif

end module hrate

!======================================================================
