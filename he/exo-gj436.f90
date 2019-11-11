!=======================================================================
!> @file exoplanet.f90
!> @brief Exoplanet problem module
!> @author M. Schneiter, C. Villarreal  D'Angelo, A. Esquivel
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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief Exoplanet module
!> @details Problem Module for exoplanet

module exoplanet

  use parameters
  implicit none

  real :: RSW     !< Stellar wind radius
  real :: TSW     !< Stellar wind temperature
  real :: VSW     !< Stellar wind velocity
  real :: dsw     !< Stellar Wind Density
  real :: RPW     !< Planetary radius
  real :: TPW     !< Planetary wind temperature
  real :: VPW     !< Planetary wind velocity
  real :: dpw     ! Planetary wind density
  real :: S0      !< Photons per sec from the stellar surface
  real :: Rstar   !< Stellar radius
  real :: Vesc, a, v0
  real :: fneutro_s  !< Stellar wind ionization fraction

  real :: torb    !< planet: orbital period
  real :: rorb     !<  orbital radius
  real :: MassS   !< Mass of the Star
  real :: MassP   !< Mass of the Planet
  real :: xp      !< X position of the planet
  real :: yp      !< Y position of the planet
  real :: zp      !< Z position of the planet
  real :: fneutro_p  !< Planet wind ionization fraction
  real :: phi

 contains

!=======================================================================

!> @brief Module initialization
!> @details Here the parameters of the Star are initialized, and scaled
!! to code units

subroutine init_exo()

  use constants, only :Rsun,Msun,yr,Ggrav,pi,mjup,au,day,Rjup,Kb
  use parameters

  implicit none
  real :: amdot  ! mdot_star (MSUN/yr)
  real :: ampdot ! mdot_planet (g/s)

  !----------------STAR PARAMETERS ------------------

  MassS = 0.452*msun
  Rstar = 0.493*Rsun
  AMDOT = 0.1*2.e-14*msun/yr            ! Stellar Mass Loss rate (g s^-1)
  Tsw   = 3.e6                       ! Stellar temperature (K)
  Vesc  = sqrt(2.*Ggrav*MassS/Rstar)  !*******      !escape velocity (cm/s)
  a     = sqrt(Kb*Tsw/(0.5*amh))
  v0    = a*(0.5*vesc/a)**2*exp(-0.5*(vesc/a)**2+3./2.)
  !  Stellar wind, imposed at the 2x sonic point (cm)
  RSW   = 1.5*Ggrav*MassS/2./(Rg*Tsw/0.5) !********************
  call get_parker(rsw,vsw)
  dsw   = ((AMDOT/RSW)/(4*pi*RSW*VSW))     ! Stellar numerical density @RS (cm^-3)
  S0    = 1.92e37/2.                       !Leuv*1.6242e11/13.6
  fneutro_s= 1.e-5

  !----------------PLANET PARAMETERS------------------

  MassP = 0.0737*Mjup
  AMPDOT= 9.76e9                  ! Planetary Mass Loss rate (g/s)
  TPW   = 4008.4                     ! Planets temperature
  RPW   = 5.*0.361*Rjup           ! Planetary wind radius (cm)
  vpw   = 12.e5                     ! Planets wind velocity (cm/s)
  dpw   = ((AMPDOT/RPW)/(4*pi*RPW*VPW))   ! Planetary wind density
  fneutro_p= 1.-0.52

  !ORBITAL PARAMETERS
  rorb= 0.0287*AU
  torb= 2.644*day

  ! change to code units
  dsw=dsw/rhosc
  vsw=vsw/vsc
  v0=v0/vsc
  Tsw=Tsw/Tempsc
  Rsw=Rsw/rsc
  Rstar=Rstar/rsc

  dpw=dpw/rhosc
  vpw=vpw/vsc
  Tpw=Tpw/Tempsc
  Rpw=Rpw/rsc

  rorb=rorb/rsc
  torb=torb/tsc

  !  initial position
  phi= 4.*pi/180.
  xp= Rorb*cos(phi)
  yp= 0.
  zp= Rorb*sin(phi)

end subroutine init_exo

!=======================================================================

!> @brief Inject sources of wind
!> @details Imposes the sources of wond from the star and planet
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserver variables
!> @param real [time] time : current integration timr
  !--------------------------------------------------------------------

subroutine impose_exo(u,time)

  use constants, only : pi
  use globals, only : coords, dx, dy, dz

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x, y, z, xpl, ypl, zpl
  real :: velx, vely, velz, rads, dens, radp
  real :: vxorb, vzorb, vyorb, omega! cpi
!  real :: vrparker
  integer ::  i,j,k!, iplanet, jplanet, kplanet
  
  omega=2.*pi/TORB

  XP=Rorb*COS(omega*TIME+phi)
  ZP=Rorb*SIN(omega*TIME+phi)

  !Orbital planetary velocity (moves in the xz-plane)
  vxorb= -omega*Rorb*sin(omega*TIME+phi)
  vzorb=  omega*Rorb*cos(omega*TIME+phi)
  vyorb=0.

  !
  ! the variables are return in physical units!!
  !vrparker=vrparker/vsc

  do i=nxmin,nxmax
    do j=nymin,nymax
      do k=nzmin,nzmax

        ! Position measured from the centre of the grid (star)
        x=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx
        y=(real(j+coords(1)*ny-nytot/2)+0.5)*dy
        z=(real(k+coords(2)*nz)+0.5)*dz

        ! Position measured from the centre of the planet
        xpl=x-xp
        ypl=y
        zpl=z-zp

        ! Distance from the centre of the star
        rads=sqrt(x**2+y**2+z**2)

        ! Distance from the centre of the planet
        radp=sqrt(xpl**2+ypl**2+zpl**2)

        ! IF INSIDE THE STAR
        if( rads <= rsw) then

          if(rads == 0.) rads=dx*0.10
          VelX=VSW*X/RADS
          VelY=VSW*Y/RADS
          VelZ=VSW*Z/RADS
          DENS=DSW*RSW**2/RADS**2

          !   total density and momenta
          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz
          !   Here the number density of the wind and planet
          !   components separately
          u(neqdyn+2,i,j,k) = (1.-fneutro_s)*dens   ! xhi*rho S ion
          u(neqdyn+3,i,j,k) = fneutro_s*dens        ! xhn*rho S neutro
          u(neqdyn+4,i,j,k) = 0.*dens               ! xci*rho P ion
          u(neqdyn+5,i,j,k) = 0.*dens               ! xcn*rho P neutro
          ! ne
          u(neqdyn+6,i,j,k) = u(neqdyn+2,i,j,k)+u(neqdyn+4,i,j,k)
          !density of neutrals
          u(neqdyn+1,i,j,k) = u(neqdyn+3,i,j,k)+u(neqdyn+5,i,j,k)
          !u(neqdyn+1,i,j,k) = fneutro_s*dens
          !   passive scalar (tag) for stellar material
          u(neqdyn+7,i,j,k) = 1000*dens

         !! if (u(neqdyn+3,i,j,k)<=0) then
         !!        print*, u(neqdyn+3,i,j,k)
         !!        stop
         !! endif

          ! total energy
          u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) &
                      + cv*(2.*dens-u(neqdyn+1,i,j,k))*Tsw
          if (u(5,i,j,k)<0) print*,'star pressure', u(5,i,j,k)

          ! IF INSIDE THE PLANET
        else if(radp <= rpw) then

          if(radp == 0.) radp=dx*0.10

          VelX=VXORB+VPW*XPL/RADP
          VelY=VYORB+VPW*YPL/RADP
          VelZ=VZORB+VPW*ZPL/RADP
          DENS=DPW*RPW**2/RADP**2

          !   total density and momenta
          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz
          !   Here the number density of the wind and planet
          !   components separately
          u(neqdyn+2,i,j,k) = 0.*dens      ! xhi*rho S ion
          u(neqdyn+3,i,j,k) = 0.*dens      ! xhn*rho S neutro
          u(neqdyn+4,i,j,k) = (1.-fneutro_p)*dens     ! xci*rho P ion
          u(neqdyn+5,i,j,k) = fneutro_p*dens     ! xcn*rho P neutro

          ! ne
          u(neqdyn+6,i,j,k) = u(neqdyn+2,i,j,k)+u(neqdyn+4,i,j,k)
          !density of neutrals
          u(neqdyn+1,i,j,k) = u(neqdyn+3,i,j,k)+u(neqdyn+5,i,j,k)
          !u(neqdyn+1,i,j,k) = fneutro_p*dens
          !   passive scalar (tag) for planetary material
          u(neqdyn+7,i,j,k)= -1000*dens

          !!if (u(neqdyn+5,i,j,k)<=0) then
          !!       print*, u(neqdyn+5,i,j,k)
          !!       stop
          !!endif

          ! total energy
          u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) &
                       + cv*(2.*dens-u(neqdyn+1,i,j,k))*Tpw
          if (u(5,i,j,k)<0) print*,'planet pressure', u(5,i,j,k)

        end if

      end do
    end do
  end do

end subroutine impose_exo

!=======================================================================
subroutine get_parker(r,Vrparker)

use constants, only: Ggrav,kB,amh

implicit none

real, intent(in)::r
real, intent(out)::vrparker
real:: vesc
real:: csS
real:: lambda
real:: Vr,Vr0, dVr,Rc,LHS, RHS, LHSold

!unscale the variables to work on physicals units
    vesc   = (2.*Ggrav*MassS/(Rstar))**0.5
    csS    = (kB*Tsw/(0.5*amh))**0.5
    lambda = 0.5*(vesc/csS)**2
! get sonic point radius
    Rc = Ggrav * MassS / (2.0 * csS**2.0)

    Vr0 = csS
    dVr = csS/10
    if (r < Rc) then
        dVr = - dVr
    end if

    LHS = Vr0 * exp(-Vr0**2.0 / (2.0 * csS**2.0) )
    RHS = csS * (Rc/r)**2.0 * exp(-2.0 * Rc/r + 3.0/2.0)

    Vr = Vr0

    do while (abs(LHS/RHS-1.0) > 1.e-8)
            ! save old LHS
            LHSold = LHS
            !update Vr
            Vr = Vr + dVr
            ! calculate new LHS
            LHS = Vr * exp(-Vr**2.0/(2.0 *csS**2.0) )
            ! see if crossed solution
            if((LHS/RHS-1.0)*(LHSold/RHS-1.0)<0)then
                dVr=-dVr/2.0
            end if
    end do
    Vrparker=Vr

end subroutine get_parker
!=======================================================================
end module exoplanet

!=======================================================================
