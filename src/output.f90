!=======================================================================
!> @file output.f90
!> @brief Writes Output
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

#ifdef OUTBIN
  use Out_BIN_Module
#endif
#ifdef OUTSILO
  use Out_Silo_Module
#endif
#ifdef OUTVTK
  use Out_VTK_Module
#endif
implicit none 
integer, intent(in) :: itprint

#ifdef OUTBIN
  call write_BIN(itprint)
#endif
#ifdef OUTVTK
  call write_VTK(itprint)
#endif
#ifdef OUTSILO     					      
  call outputsilo(itprint)
#endif							      

end subroutine write_output

!=======================================================================

end module output

!=======================================================================
