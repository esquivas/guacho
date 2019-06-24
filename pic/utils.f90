!=======================================================================
!> @file globals.f90
!> @brief Utilities module
!> @author Alejandro Esquivel
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
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief Module containing general purpose utilities
!> @details This module contains utilities that could be called from then
!> user module (and any other as well)

module utilities

  implicit none

contains

  !================================================================
  ! @brief In domain function (logical)
  ! @details Determines if the position of a given point lies within the
  !> procesor domain, returns true if it is, false if it isn't
  ! @param real [in] : 3D position with respect to a cornet of the domain
  function isInDomain(pos)

    use parameters, only : nx, ny, nz
    use globals,    only : coords, dx, dy, dz
    implicit none
    logical          :: isInDomain
    real, intent(in) :: pos(3)
    integer          :: ind(3)

    ! shift to the local processor
    ind(1) = int(pos(1)/dx) - coords(0)*nx
    ind(2) = int(pos(2)/dy) - coords(1)*ny
    ind(3) = int(pos(3)/dz) - coords(2)*nz

    if ( ind(1)<0  .or. ind(2)<0  .or. ind(3)<0 .or. &
         ind(1)>=nx .or. ind(2)>=ny .or. ind(3)>=nz ) then

      isInDomain = .false.

    else

      isInDomain = .true.

    end if

  end function isInDomain

    !================================================================
  ! @brief In Which domain function
  ! @details Determines if the rank that has the position of a given point,
  !> return -1 if point lies outside domain
  ! @param real [in] : 3D position with respect to a cornet of the domain
  function inWhichDomain(pos)

    use parameters, only : nx, ny, nz
    use globals,    only : dx, dy, dz, comm3d
#ifdef MPIP
    use mpi
#endif
    implicit none
    integer          :: inWhichDomain
    real, intent(in) :: pos(3)
    integer          :: ind(0:2),err

    ! get coord of rank
    ind(0) = int(pos(1)/dx)/nx
    ind(1) = int(pos(2)/dy)/ny
    ind(2) = int(pos(3)/dz)/nz

#ifdef MPIP
    call mpi_cart_rank(comm3d,ind,inWhichDomain,err)
#else
     if(.not.isInDomain(pos)) inWhichDomain = -1
#endif

  end function inWhichDomain


end module utilities

!=======================================================================
