!=======================================================================
!> @file globals.f90
!> @brief Global variables
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

!> @brief Module containing global variables
!> @details This mudules contains variables that are treated as global
!! in the code

module globals
  implicit none

  real, allocatable ::      u(:,:,:,:) !< conserved varibles
  real, allocatable ::     up(:,:,:,:) !< conserved varibles after 1/2 timestep
  real, allocatable :: primit(:,:,:,:) !< primitive varibles
  real, allocatable ::      f(:,:,:,:) !< X fluxes
  real, allocatable ::      g(:,:,:,:) !< Y fluxes
  real, allocatable ::      h(:,:,:,:) !< Z fluxes
  real, allocatable ::   Temp  (:,:,:) !< Temperature array [K]

  real :: dx  !< grid spacing in X
  real :: dy  !< grid spacing in Y
  real :: dz  !< grid spacing in Z

  !> position of neighboring MPI blocks
  integer, dimension(0:2) :: coords

  integer :: left   !< MPI neighbor in the -x direction
  integer :: right  !< MPI neighbor in the +x direction
  integer :: top    !< MPI neighbor in the -y direction
  integer :: bottom !< MPI neighbor in the +y direction
  integer :: out    !< MPI neighbor in the -z direction
  integer :: in     !< MPI neighbor in the +z direction

  integer :: rank     !< MPI rank
  integer :: comm3d   !< Cartessian MPI comunicator

  !> Current time
  real :: time
  !> Current CFL $\Delta t$
  real :: dt_CFL
  !> Current iteration
  integer :: currentIteration

end module globals

!=======================================================================
