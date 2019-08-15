!=======================================================================
!> @file utils.f90
!> @brief Utilities module
!> @author Alejandro Esquivel
!> @date 1/jul/2019
!
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

!> @brief Module containing general purpose utilities
!> @details This module contains utilities that could be called from then
!> user module as they are compiled before it

module utilities

  implicit none

contains

  !================================================================
  !> @brief Is in domain? function (logical)
  !> @details Determines if the position of a given point lies within the
  !> procesor domain, returns true if it is, false if it isn't
  !> @param real [in] : 3D position WITH RESPECT TO A CORNER of the domain
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
         ind(1)>nx .or. ind(2)>ny .or. ind(3)>nz ) then

      isInDomain = .false.

    else

      isInDomain = .true.

    end if

  end function isInDomain

  !================================================================
  !> @brief is in shock? function (logical)
  !> @details Determines if the position of a given point is inside a
  !> shocked region (as tracked by the global array 'shockF', see also
  !> flag_shock subroutine). Returns 'true' if it is, 'false' if it isn't
  !> @param real [in] : 3D position with respect to a cornet of the domain
  function isInShock(pos)
    use parameters, only : nx, ny, nz
    use globals,    only : coords, dx, dy, dz, shockF
    implicit none
    logical          :: isInShock
    real, intent(in) :: pos(3)
    integer          :: i, j, k

    ! shift to the local processor
    i = int( pos(1)/dx ) - coords(0)*nx + 1
    j = int( pos(2)/dy ) - coords(1)*ny + 1
    k = int( pos(3)/dz ) - coords(2)*nz + 1

    isInShock = .false.

    if ( i < 1  .or. j < 1  .or. k < 1 .or. &
         i > nx .or. j > ny .or. k > nz  ) then

      return

    elseif ( shockF(i,j,k) == 1 ) then

      isInShock = .true.
      return

    end if

  end function isInShock

  !================================================================
  !> @brief In Which domain function
  !> @details Determines if the rank that has the position of a given point,
  !> returns -1 if point lies outside domain
  ! @param real [in] : 3D position with respect to a corner of the domain
  function inWhichDomain(pos)
    use parameters, only : nx, ny, nz
    use globals,    only : dx, dy, dz, comm3d
#ifdef MPIP
    use mpi
#endif
    implicit none
    integer          :: inWhichDomain
    real, intent(in) :: pos(3)
    integer          :: ind(0:2), err

    ! get coords of rank
    ind(0) = int(pos(1)/dx)/nx
    ind(1) = int(pos(2)/dy)/ny
    ind(2) = int(pos(3)/dz)/nz

#ifdef MPIP
    if(isInDomain(pos)) then
      call mpi_cart_rank(comm3d,ind,inWhichDomain,err)
    else
      inWhichDomain = -1
    endif
#else
    if(.not.isInDomain(pos)) inWhichDomain = -1
#endif

  end function inWhichDomain

  !=======================================================================
  !> @brief shock detector
  !> @details Schock detection, similar to Mignone et al 2012
  !> @param integer [out] shock :: shock (one if shocked material)
  subroutine flag_shock()
    use parameters
    use globals,    only : dx, dy, dz, primit, shockF
    implicit none
    integer :: i,j,k
    !>  threshold for shock detection in thermal pressure gradient
    real, parameter :: epsilon_sh = 3.  !  changed from 3.
    real    :: gradP, divV

    do k =1,nz
      do j =1,ny
        do i =1,nx

          shockF(i,j,k) = 0

          gradP = abs(primit(5,i-1,j,k)-primit(5,i+1,j,k))/                    &
                  min(primit(5,i-1,j,k),primit(5,i+1,j,k)) +                   &
                  abs(primit(5,i,j-1,k)-primit(5,i,j+1,k))/                    &
                  min(primit(5,i,j-1,k),primit(5,i,j+1,k)) +                   &
                  abs(primit(5,i,j,k-1)-primit(5,i,j,k+1))/                    &
                  min(primit(5,i,j,k-1),primit(5,i,j,k+1))

          if (gradP >= epsilon_sh) then
          divV =  (primit(2,i+1,j,k)-primit(2,i-1,j,k))/(2.*dx)                &
                + (primit(3,i,j+1,k)-primit(3,i,j-1,k))/(2.*dy)                &
                + (primit(4,i,j,k+1)-primit(4,i,j,k-1))/(2.*dz)

          if (divV < 0.) shockF(i,j,k) = 1

          end if

        end do
      end do
    end do

end subroutine flag_shock

  !=======================================================================

end module utilities
