!=======================================================================
!> @file constants.f90
!> @brief Constants module
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

!> @brief Module containing physical and asronomical constants

module constants
  implicit none

  real, parameter :: pi=acos(-1.)      !< @f$ \pi @f$
  real, parameter :: amh=1.66e-24      !< hydrogen mass
  real, parameter :: mu=0.5            !< mean atomic mass
  real, parameter :: Kb=1.38e-16       !< Boltzmann constant (cgs)  
  real, parameter :: Rg=8.3145e7       !< Gas constant (cgs)
  real, parameter :: Ggrav=6.67259e-8  !< Gravitational constant (cgs)
  real, parameter :: clight=2.99E10    !< speed of light in vacuum (cgs)

  real, parameter :: Msun=1.99E33      !< solar radius (cgs)
  real, parameter :: Rsun=6.955e10     !< solar mass (cgs)
  real, parameter :: Mjup=1.898E30     !< Jupiter mass (cgs)
  real, parameter :: Rjup=7.1492E9     !< Jupiter radius (cgs)
     
  real, parameter :: AU=1.496e13       !< 1AU in cm
  real, parameter :: pc=3.0857E18      !< 1pc in cm
  real, parameter :: kpc=3.0857E21     !< 1Kpc in cm
  real, parameter :: hr=3600.          !< 1hr in seconds
  real, parameter :: day=86400.        !< 1day in seconds
  real, parameter :: yr=3.1536E7       !< 1yr in seconds
  real, parameter :: Myr=3.1536E13     !< 1Myr in seconds

end module constants

!=======================================================================
