!=======================================================================
!> @file self_gravity.f90
!> @brief Guacho-3D main program
!> @author Veronica Lora & Alejandro Esquivel
!> @date 9/March/2020
!
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

!> @brief Self gravity Module
!> @details Solves the Poisson equation with a successive over-relaxation (SOR)
!> method to get the gravitational potential.
!> The sources are added in 'sources.f90'
module self_gravity
  use constants, only : pi, Ggrav
  implicit none

  real, allocatable :: phi_grav(:,:,:)
  real, parameter :: four_pi_G = 4.*pi*Ggrav

contains

  !================================================================
  !> @brief Initialization of module
  !> @details Allocates memory for all global variables that correspond
  !> to the module
  subroutine init_self_gravity()

    implicit none

  end subroutine init_self_gravity

  !================================================================
  !> @brief Solve Poisson equation
  !> @details Compute the gravitational potential with a SOR mehtod
  subroutine solve_poisson()

    implicit none

  end subroutine solve_poisson

  !================================================================
  !> @brief Add self gravity sources
  !> @details Adds the sources due to the gravitationl potential, the gradient
  !> of the potential (gravityationa force), and the work done by it.
  subroutine add_self_gravity()

    implicit none

  end subroutine add_self_gravity

  !================================================================

end module self_gravity

!================================================================
