!=======================================================================
!> @file cooling_chem.f90
!> @brief Cooling with hydrogen rate parametrized cooling
!> @author Alejandro Esquivel & Carolina Villarreal
!> @date 11/Nov/2019

! Copyright (c) 2019 Guacho Co-Op
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

!> @brief Cooling with chemistry module
!> @details Cooling with chemistry module including Helium 10380 line

module cooling_chem

#ifdef PASSIVES

  implicit none

contains

  !=======================================================================
  subroutine init_cooling_chem()

    use cooling_schure, only : init_cooling_schure
    implicit none

    call init_cooling_schure()

  end subroutine init_cooling_chem

  !=======================================================================
  !> @brief Cooling module using the chemistry
  !> @details Cooling using the ionization states computed with the ionic/
  !> chemistry module
  subroutine chem_cool()

    use cooling_schure, only : coolingschure
    implicit none

    call coolingschure()

  end subroutine chem_cool

  !=======================================================================

#endif

end module cooling_chem
