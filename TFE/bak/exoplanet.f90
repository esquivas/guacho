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

module exoplanet

  use parameters
  implicit none

  !Stellar wind parameters
  real :: RSW     !< Radius at wich sw is imposed
  real :: TSW     !< Stellar wind temperature
  real :: VSW     !< Stellar wind velocity
  real :: dsw     !< Stellar Wind Density
  real :: bsw     !< Magnetic Field

  !Environment parameters
  real :: Tenv     !< Stellar wind temperature
  real :: Venv     !< Stellar wind velocity
  real :: denv     !< Stellar Wind Density
  real :: benv     !< Magnetic Field

  !Planetary wind parameters
  real :: bpw     !< Planetary Magnetic Field  
  real :: RPW     !< Planetary radius
  real :: TPW     !< Planetary wind temperature
  real :: VPW     !< Planetary wind velocity
  real :: dpw     ! Planetary wind density

  !Some other planetary parameters
  real :: MassP   !< Mass of the Planet
  real :: xp      !< X position of the planet
  real :: yp      !< Y position of the planet
  real :: zp      !< Z position of the planet
  real :: rhop    !< planetary density
  real :: volp    !< planetary volume
  real :: rorb    !<  orbital radius
  
  !Some other stellar parameters
  real :: MassS   !< Mass of the Star
  real :: Rss     !< star radius

  !  Flag to tag the presence of the planet
  logical, allocatable :: flagP(:,:,:)
contains

!=======================================================================

!> @brief Module initialization
!> @details Here the parameters of the Star are initialized, and scaled 
!! to code units

subroutine init_exo()
  use globals,    only: coords, dx ,dy ,dz
  use constants, only : Msun, Yr, Ggrav, pi, mjup, au, day,Rjup, rsun
  implicit none
  real :: amdot  ! mdot_star (MSUN/yr)
  real :: ampdot ! mdot_planet (g/s)
  integer :: i,j,k
  real    :: xpl, ypl, zpl, radp


  !ORBITAL PARAMETERS
  rorb=.047*AU!0.47**AU
  
  !----------------STAR PARAMETERS ------------------
  MassS = 1.1*msun
  RsS   = 1.2*rsun
  AMDOT = 2.E-14*msun/yr              ! Stellar Mass Loss rate (g s^-1)
  TSW   = 1.56E6     !*************   ! Stellar temperature (K)
  !  Stellar wind, imposed at the 1.5x  sonic point (cm)
  RSW   = rorb-xphys/2. !Distance from centre of star a.w.t sw is imposed
  vsw   = 200.e5     !*************      ! Stellar wind velocity (cm/s)
  dsw   = ((AMDOT/RSW)/(4*pi*RSW*VSW))   ! Stellar density @RS (g cm^-3)
  bsw   = 1.E-5                          ! Stellar magnetic field (g)


  !----------------IONIZED ENVIRONMENT PARAMETERS ------------------
!  denv   = 1.1*msun
!         = 
!         =      ! Stellar Mass Loss rate (g s^-1)
!         =      ! Stellar temperature (K)
!         = 
!         =      ! Stellar wind velocity (cm/s)
!         =      ! Stellar density @RS (g cm^-3)
!         =      ! Stellar magnetic field (g)
!

  
  !----------------PLANET PARAMETERS------------------
  MassP = 0.67*mjup
  AMPDOT= 1.E10   !***********         ! Planetary Mass Loss rate (g/s)
  TPW   = 1E4                          ! Planets temperature
  RPW   = 1.38*Rjup                    ! Planetary wind radius (cm) 
  VolP  = 4.*pi*RPW**3./3.             ! Planets density
  rhoP  = MassP/VolP                   ! Planets density
  vpw   = 60.e5 !Ves at 1R_e =42 km/s  ! Planets wind velocity (cm/s)
  dpw=((AMPDOT/RPW)/(4*pi*RPW*VPW))    ! Planetary wind density
  bpw   = 0.04                         ! Planetary magnetic field (g)

  
  ! change to code units
  dsw=dsw/rhosc
  vsw=vsw/sqrt(vsc2)
  Tsw=Tsw/Tempsc
  Rsw=Rsw/rsc
  RsS=RsS/rsc
  bsw=bsw/bsc 
  bpw=bpw/bsc
  dpw=dpw/rhosc
  vpw=vpw/sqrt(vsc2)
  Tpw=Tpw/Tempsc
  Rpw=Rpw/rsc

  !  tag planet
  allocate(flagP(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
  flagP(:,:,:) = .False.
  do k=nzmin,nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax

        ! Position measured from the centre of the grid (planet)
        xpl=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
        ypl=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
        zpl=(float(k+coords(2)*nz-nztot/2)+0.5)*dz

        ! Distance from the centre of the planet (centred)
        radp=sqrt(xpl**2+ypl**2+zpl**2)
        ! IF INSIDE THE PLANET
        if(radp <= rpw) flagP(i,j,k)=.true.

      end do
    end do
  end do

  
end subroutine init_exo

!=======================================================================

!> @brief Inject sources of wind
!> @details Imposes the sources of wond from the star and planet
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserver variables
!> @param real [time] time : current integration timr
  !--------------------------------------------------------------------

subroutine impose_exo(u,time)
!subroutine impose_exo(u,time)
  
  use constants, only : pi
  use globals, only : coords, dx, dy, dz
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x, y, z, xpl, ypl, zpl
  real :: velx, vely, velz, rads, dens, radp, phi
#ifdef BFIELD
  real :: cpi
#endif
  integer ::  i,j,k

  do k=nzmin,nzmax
     do j=nymin,nymax
          do i=nxmin,nxmax
           
           ! Position measured from the centre of the grid (planet)
           xpl=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
           ypl=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
           zpl=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
           
           ! Distance from the centre of the planet (centred)
           radp=sqrt(xpl**2+ypl**2+zpl**2)
           
           ! IF INSIDE THE PLANET
                   
           if(radp <= 1.1*rpw) then
              if(radp == 0.) radp=dx*0.10
              VelX=VPW*XPL/RADP
              VelY=VPW*YPL/RADP
              VelZ=VPW*ZPL/RADP
              DENS=DPW!*RPW**2/RADP**2
              
              if(radp <= 0.9*rpw) then
                 
                 VelX=0.!VPW*XPL/RADP
                 VelY=0.!VPW*YPL/RADP
                 VelZ=0.!VPW*ZPL/RADP
                 DENS=DPW!*RPW**2/RADP**2
                 !                 DENS=DPW*RPW**2/RADP**2
              end if
!           else if (radp > 1.1*rpw.and.time.eq.0) then
!              VelX=VPW*XPL/RADP
!              VelY=VPW*YPL/RADP
!              VelZ=VPW*ZPL/RADP
!              DENS=DPW*RPW**2/RADP**2
!           end if
              !   total density and momenta
                 
              u(1,i,j,k) = dens
              u(2,i,j,k) = dens*velx
              u(3,i,j,k) = dens*vely
              u(4,i,j,k) = dens*velz
           
              !  Magnetic fields (IF MHD or PMHD)
#ifdef BFIELD
!              cpi = bpw*(rpw/radp)**3/(2.*radp**2)
!              u(6,i,j,k) = 3.*ypl*xpl*cpi
!              u(7,i,j,k) = (3.*ypl**2-radp**2)*cpi
!              u(8,i,j,k) = 3.*ypl*zpl*cpi
              
              u(6,i,j,k) = 0.
              u(7,i,j,k) = 0.
              u(8,i,j,k) = 0.
              
#endif
              !   energy
              if (mhd) then
#ifdef BFIELD
                 u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2) &
                      + cv*dens*Tpw                !    & 
 !                     + 0.5*(u(6,i,j,k)**2+u(7,i,j,k)**2+u(8,i,j,k)**2)
#endif
              else
                 u(5,i,j,k)=0.5*dens*(velx**2+vely**2+velz**2) &
                      + cv*dens*Tpw
              end if
              
              if (passives) then
#ifdef PASSIVES
                 !  density of neutrals
                 u(neqdyn+1,i,j,k)=1.*dens                
                 !   passive scalar (h-hot, c-cold, i-ionized, n-neutral)
                 u(neqdyn+2,i,j,k)= -dens   ! passive scalar
              endif
#endif
              
           end if
           
        end do
     end do
  end do
  
end subroutine impose_exo

!=======================================================================

end module exoplanet

!=======================================================================
