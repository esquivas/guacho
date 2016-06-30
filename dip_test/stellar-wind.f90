!=======================================================================
!> @file exoplanet.f90
!> @brief Exoplanet problem module
!> @author M. Schneiter, C. Villarreal  D'Angelo, A. Esquivel
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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief Exoplanet module
!> @details Problem Module for exoplanet

module stellarwind

  use parameters
  implicit none
  real :: RSW     !< Stellar wind radius
  real :: TSW     !< Stellar wind temperature
  real :: VSW     !< Stellar wind velocity
  real :: dsw     !< Stellar Wind Density
  real :: RsS
  real :: bsw
  real :: MassS   !< Mass of the Star

contains

!=======================================================================

!> @brief Module initialization
!> @details Here the parameters of the Star are initialized, and scaled 
!! to code units

subroutine init_wind()
  
  use constants, only : Msun, Yr, Ggrav, pi,  rsun

  implicit none
  real :: amdot  ! mdot_star (MSUN/yr)

  !----------------STAR PARAMETERS ------------------
  MassS = msun
  RsS   = rsun
  AMDOT = 2.E-14*msun/yr              ! Stellar Mass Loss rate (g s^-1)
  TSW   = 1.56E6     !*************   ! Stellar temperature (K)
  !  Stellar wind, imposed at the 1.5x  sonic point (cm)
  RSW   = rsun   !*************
  VSW   = 10.e5       !*************      ! Stellar wind velocity (cm/s)
  dsw   = 1.544e-16!((AMDOT/RSW)/(4*pi*RSW*VSW))   ! Stellar density @RS (g cm^-3)
  bsw   = 1.0                            ! Stellar magnetic field (g)

  ! change to code units
  dsw=dsw/rhosc
  vsw=vsw/vsc
  Tsw=Tsw/Tempsc
  Rsw=Rsw/rsc
  RsS=RsS/rsc
  bsw=bsw/bsc 

end subroutine init_wind

!=======================================================================

!> @brief Inject sources of wind
!> @details Imposes the sources of wond from the star and planet
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserver variables
!> @param real [time] time : current integration timr
  !--------------------------------------------------------------------

subroutine impose_wind(u,time)
  
  use constants, only : pi
  use globals, only : coords, dx, dy, dz
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x, y, z
  real :: velx, vely, velz, rads, dens
#ifdef BFIELD
  real :: cpi
#endif
  integer ::  i,j,k

  do i=nxmin,nxmax
    do j=nymin,nymax
      do k=nzmin,nzmax

        ! Position measured from the centre of the grid (star)
        x=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx
        y=(real(j+coords(1)*ny-nytot/2)+0.5)*dy
        z=(real(k+coords(2)*nz-nztot/2)+0.5)*dz

        ! Distance from the centre of the star
        rads=sqrt(x**2+y**2+z**2)

        ! IF INSIDE THE STAR
        if( rads <= rsw) then
          if(rads == 0.) rads=dx*0.10

          VelX=VSW*X/RADS
          VelY=VSW*Y/RADS
          VelZ=VSW*Z/RADS
          DENS=DSW!*RSW**2/RADS**2
          !   total density and momena
          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz
          !   magnetic field

          if (pmhd .or. mhd) then
#ifdef BFIELD
            if (rads <= 0.8*rsw) then
              cpi = bsw*(1./0.8 )**3/(2.*rads**2)
            else
              cpi = bsw*(RSW/rads)**3/(2.*rads**2)
            end if
            u(6,i,j,k) = 3.*y*x*cpi
            u(7,i,j,k) = (3.*y**2-rads**2)*cpi
            u(8,i,j,k) = 3.*y*z*cpi
#endif
          end if

          if (mhd) then
#ifdef BFIELD
            ! total energy
            u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2)        &
            + cv*dens*1.9999*Tsw       & 
            + 0.5*(u(6,i,j,k)**2+u(7,i,j,k)**2+u(8,i,j,k)**2)
#endif
          else
            ! total energy
            u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2) &
            + cv*dens*Tsw
          endif

          if (passives) then
#ifdef PASSIVES
            !  density of neutrals
            u(neqdyn+1,i,j,k)= 0.0001*dens
#endif
          end if

        end if

      end do
    end do
  end do

end subroutine impose_wind

!=======================================================================

end module stellarwind
!=======================================================================
