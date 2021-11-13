!=======================================================================
!> @file constants.f90
!> @brief Constants module
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

!> @brief Module containing physical, asronomical constants, and other
!> named constants

module constants
  implicit none

  real, parameter :: pi=acos(-1.)         !< @f$ \pi @f$
  real, parameter :: amh=1.66e-24         !< hydrogen mass
  real, parameter :: Kb=1.38e-16          !< Boltzmann constant (cgs)
  real, parameter :: Rg=8.3145e7          !< Gas constant (cgs)
  real, parameter :: Ggrav=6.67259e-8     !< Gravitational constant (cgs)
  real, parameter :: clight=2.99E10       !< speed of light in vacuum (cgs)
  real, parameter :: echarge=4.8032e-10   !< electron charge statcoulomb (cgs)
  real, parameter :: emass=9.10938e-28    !< electron mass (g)
  real, parameter :: sigma_SB=5.6704e-5   !< Stephan Boltzmann constant (cgs)
  real, parameter :: sigma_T =6.65245e-25 !< Thompson-scattering cross section

  real, parameter :: Msun=1.99E33         !< solar radius (cgs)
  real, parameter :: Rsun=6.955e10        !< solar mass (cgs)
  real, parameter :: gsun=274.e2          !< solar gravity (cgs)
  real, parameter :: Mjup=1.898E30        !< Jupiter mass (cgs)
  real, parameter :: Rjup=7.1492E9        !< Jupiter radius (cgs)

  real, parameter :: AU=1.496e13          !< 1AU in cm
  real, parameter :: pc=3.0857E18         !< 1pc in cm
  real, parameter :: kpc=3.0857E21        !< 1Kpc in cm
  real, parameter :: deg = pi/180.        !< conversion from deg to rad
  real, parameter :: hr=3600.             !< 1hr in seconds
  real, parameter :: day=86400.           !< 1day in seconds
  real, parameter :: yr=3.1536E7          !< 1yr in seconds
  real, parameter :: Myr=3.1536E13        !< 1Myr in seconds
  real, parameter :: eV=1.60218E-12       !< 1 ev in ergs

  !  Named constants

  !  Approximate Riemann solvers
  integer, parameter :: SOLVER_HLL            =  1
  integer, parameter :: SOLVER_HLLC           =  2
  integer, parameter :: SOLVER_HLLE           =  3
  integer, parameter :: SOLVER_HLLD           =  4
  integer, parameter :: SOLVER_HLLE_SPLIT_B   =  5
  integer, parameter :: SOLVER_HLLD_SPLIT_B   =  6
  integer, parameter :: SOLVER_HLLE_SPLIT_ALL =  7
  integer, parameter :: SOLVER_HLLD_SPLIT_ALL =  8
  integer, parameter :: SOLVER_RHLL           =  9
  integer, parameter :: SOLVER_RHLLC          = 10

  !  Equations of state
  integer, parameter :: EOS_ADIABATIC     = 1
  integer, parameter :: EOS_SINGLE_SPECIE = 2
  integer, parameter :: EOS_H_RATE        = 3
  integer, parameter :: EOS_CHEM          = 4
  integer, parameter :: EOS_REL_IDEAL     = 5
  integer, parameter :: EOS_REL_TM        = 6

  !  Cooling Schemes
  integer, parameter :: COOL_NONE  = 0
  integer, parameter :: COOL_H     = 1
  integer, parameter :: COOL_BBC   = 2
  integer, parameter :: COOL_DMC   = 3
  integer, parameter :: COOL_CHI   = 4
  integer, parameter :: COOL_SKKKV = 5
  integer, parameter :: COOL_CHEM  = 6

  !  Boundary conditions
  integer, parameter :: BC_OUTFLOW  = 1
  integer, parameter :: BC_CLOSED   = 2
  integer, parameter :: BC_PERIODIC = 3
  integer, parameter :: BC_OTHER    = 4

  !  Slope limiters
  integer, parameter :: LIMITER_NO_AVERAGE = -1
  integer, parameter :: LIMITER_NO_LIMIT   =  0
  integer, parameter :: LIMITER_MINMOD     =  1
  integer, parameter :: LIMITER_VAN_LEER   =  2
  integer, parameter :: LIMITER_VAN_ALBADA =  3
  integer, parameter :: LIMITER_UMIST      =  4
  integer, parameter :: LIMITER_WOODWARD   =  5
  integer, parameter :: LIMITER_SUPERBEE   =  6

  !  Thermal conduction
  integer, parameter :: TC_OFF         = 0
  integer, parameter :: TC_ISOTROPIC   = 1
  integer, parameter :: TC_ANISOTROPIC = 2

end module constants
