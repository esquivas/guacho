!=======================================================================
!> @file 2anilloW.f90
!> @2Wind problem module
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

module anilloW 

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

  real :: omega, phi0

contains

!=======================================================================
!> @brief Module initialization
!> @details Here the parameters of the Star are initialized, and scaled
!! to code units
subroutine init_anilloW()

  use constants, only : Msun, Yr, Ggrav, pi, mjup, au, day,Rjup,deg
  use globals, only :rank, dx
  use parameters
  implicit none
  real :: amsdot  ! mdot_star (MSUN/yr)

  !------------  (M1) COMPANION STAR PARAMETERS ------------
  ! Viento de M1
  MassS = 0.7*msun
  AMSDOT= 5.e-6*msun/yr                 ! M1 Mass Loss rate (g/s)
  TSW   = 100.0 ! 5.e5
  RSW   = 150.*AU                        ! M1 wind radius 10 cells
  VSW   = 15.0e5                         ! M1 wind velocity (cm/s)
  DSW   = ((AMSDOT/RSW)/(4*pi*RSW*VSW))                
  
  !ORBITAL PARAMETERS 
  rorb  = 100.*AU  
  torb  = 558.*yr
  phi0  = 0.

  ! change to code units
  dsw=dsw/rhosc
  vsw=vsw/vsc
  Tsw=Tsw/Tempsc
  Rsw=Rsw/rsc

  ! Orbital parameters in code units
  rorb  = rorb/rsc
  torb  = torb/tsc
  omega = 2.*pi / torb

  !  initial position
  xp= 0.0  !Rorb*cos(phi0*pi/180.)
  yp= 0.0  !Rorb*sin(phi0*pi/180.)
  zp= 0.0


end subroutine init_anilloW

!=======================================================================
!> @brief Inject sources of wind
!> @details Imposes the sources of wond from the star and planet
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserver variables
!> @param real [time] time : current integration timr
  !--------------------------------------------------------------------
subroutine impose_anilloW(u,time)

  use constants, only : pi
  use globals, only : coords, dx, dy, dz, rank
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x, y, z, xs, ys, zs
  real :: velx, vely, velz, rads, dens
  real :: vxorb, vzorb, vyorb
  integer ::  i, j, k

  !XP=Rorb*COS(omega*TIME+phi)
  !YP=Rorb*SIN(omega*TIME+phi)
  
  !Orbital velocity (moves in the xy-plane)
  vxorb= -omega*Rorb*sin(omega*TIME+phi0)
  vyorb=  omega*Rorb*cos(omega*TIME+phi0)
  vzorb=0.
        
  do k=nzmin,nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax

        ! Position measured from the centre of the grid (star)
        x=(real(i+coords(0)*nx-nxtot/2)-0.5)*dx
        y=(real(j+coords(1)*ny-nytot/2)-0.5)*dy
        z=(real(k+coords(2)*nz-nztot/2)-0.5)*dz

        ! Position measured from the centre of the wind source
        xs=x-xp
        ys=y-yp
        zs=z-zp

        ! Distance from the centre of the planet
        rads=sqrt(xs**2+ys**2+zs**2)

        ! IF INSIDE THE WIND SOURCE
        if( rads <= rsw) then

          if(rads == 0.) rads=dx*0.10
          VelX=VSW*X/RADS
          VelY=VSW*Y/RADS
          VelZ=VSW*Z/RADS
          DENS=DSW*RSW**2/RADS**2
          !   total density and momenta
          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*(velx + vxorb)
          u(3,i,j,k) = dens*(vely + vyorb)
          u(4,i,j,k) = dens*(velz + vzorb)
          ! total energy
          u(5,i,j,k)=0.5*dens*((velx+vxorb)**2+(vely+vyorb)**2+(velz+vzorb)**2) &
          + cv*dens*1.9999*Tsw

          !  density of neutrals
          u(neqdyn+1,i,j,k) =  1.E-4*dens
          !  passive scalar (tag) for stellar material
          u(neqdyn+2,i,j,k)= 1000.*dens

       endif
      end do
    end do
  end do

end subroutine impose_anilloW

!=======================================================================

end module anilloW

!=======================================================================
