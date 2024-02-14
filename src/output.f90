!=======================================================================
!> @file output.f90
!> @brief Writes Output
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

!> @brief Writes output
!> @details This module writes the ouput in the formats specified in
!! the makefile

module output

contains

  !=======================================================================
  !> @brief Writes output
  !> @details Writes output, the format is chosen in makefile
  !! @n Supported formats are
  !! *.bin and VTK (both BINARY), Silo (+hdf5)
  !> @param integer [in] itprint : number of output
  subroutine write_output(itprint)

    use parameters, only : out_bin, out_vtk, out_silo, enable_lmp, master
    use Out_BIN_Module
    use Out_Silo_Module
    use Out_VTK_Module
    use lmp_module
    implicit none
    integer, intent(in) :: itprint

    if (out_bin )   call write_BIN(itprint)
    if (out_vtk )   call write_VTK(itprint)
    if (out_silo)   call write_silo(itprint)
    if (enable_lmp) call write_LMP(itprint)

  end subroutine write_output

  !=======================================================================

end module output

!=======================================================================
