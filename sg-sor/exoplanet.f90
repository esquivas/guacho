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

  real :: torb    !< planet: orbital period
  real :: rorb    !<  orbital radius
  real :: MassS   !< Mass of the Star
  real :: MassP   !< Mass of the Planet
  real :: xp      !< X position of the planet
  real :: yp      !< Y position of the planet
  real :: zp      !< Z position of the planet
 contains

!=======================================================================

!> @brief Module initialization
!> @details Here the parameters of the Star are initialized, and scaled
!! to code units

subroutine init_exo()

  use constants, only : Msun, Yr, Ggrav, pi, mjup, au, day,Rjup
  use globals, only :rank, dx
  use parameters
  implicit none
  real :: amdot  ! mdot_star (MSUN/yr)
  real :: ampdot ! mdot_planet (g/s)

  !----------------STAR PARAMETERS ------------------
  MassS = 1.1*msun
  AMDOT = 2.E-14*msun/yr              ! Stellar Mass Loss rate (g s^-1)
  TSW   = 1.3E6 !3.0E6                ! Stellar temperature (K)
  !  Stellar wind, imposed at the 1.5x  sonic point (cm)
  RSW   = 1.1*Ggrav*MassS/2./(Rg*Tsw/0.6)
  vsw   = 200.e5                      ! Stellar wind velocity (cm/s)
  dsw   =((AMDOT/RSW)/(4*pi*RSW*VSW)) ! Stellar density @RS (g cm^-3)

  !----------------PLANET PARAMETERS------------------
  MassP =0.67*mjup
  AMPDOT= 5.e10                       ! Planetary Mass Loss rate (g/s)
  TPW   = 1.E4                        ! Planet temperature
  RPW   = 2.* 3.*1.38*Rjup            ! Planetary wind radius (cm)
  vpw   = 10.e5                       ! Planets wind velocity (cm/s)
  dpw=((AMPDOT/RPW)/(4*pi*RPW*VPW))   ! Planetary wind density

  !ORBITAL PARAMETERS
  rorb= 0.047*AU !0.47**AU
  torb= 3.52*day

  ! change to code units
  dsw=dsw/rhosc
  vsw=vsw/sqrt(vsc2)
  Tsw=Tsw/Tempsc
  Rsw=Rsw/rsc

  dpw=dpw/rhosc
  vpw=vpw/sqrt(vsc2)
  Tpw=Tpw/Tempsc
  Rpw=Rpw/rsc

  rorb=rorb/rsc
  torb=torb/tsc

  !  initial position
  xp= Rorb*cos(-25.*pi/180.)
  yp= 0.
  zp= Rorb*sin(-25.*pi/180.)

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
  use globals, only : coords, dx, dy, dz, rank
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x, y, z, xpl, ypl, zpl
  real :: velx, vely, velz, rads, dens, radp, phi!, a
  real :: vxorb, vzorb, vyorb, omega!, deltax ! cpi
  integer ::  i,j,k!, iplanet, jplanet, kplanet
  phi=-25.*pi/180.

  omega=2.*pi/TORB

  XP=Rorb*COS(omega*TIME+phi)
  ZP=Rorb*SIN(omega*TIME+phi)

  !Orbital planetary velocity (moves in the xz-plane)
  vxorb= -omega*Rorb*sin(omega*TIME+phi)
  vzorb=  omega*Rorb*cos(omega*TIME+phi)
  vyorb=0.

  do k=nzmin,nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax

        ! Position measured from the centre of the grid (star)
        x=(real(i+coords(0)*nx-nxtot/2)-0.5)*dx
        y=(real(j+coords(1)*ny-nytot/2)-0.5)*dy
        z=(real(k+coords(2)*nz-nztot/2)-0.5)*dz

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
          ! total energy
          u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2) &
          + cv*dens*1.9999*Tsw

          !  density of neutrals
          u(neqdyn+1,i,j,k) =  1.E-4*dens
          !  passive scalar (tag) for stellar material
          u(neqdyn+2,i,j,k)= 1000*dens

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

          ! total energy
          u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2) &
          + cv*dens*1.5*Tpw

          u(neqdyn+1,i,j,k) = 0.5*dens

          !   passive scalar (tag) for planetary material
          u(neqdyn+2,i,j,k)= -1000*dens

        end if

      end do
    end do
  end do

end subroutine impose_exo

!=======================================================================

end module exoplanet

!=======================================================================
