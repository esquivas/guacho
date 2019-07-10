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

!    #ifdef MPIP
!        call mpi_cart_rank(comm3d,ind,inWhichDomain,err)
!    #else
!         if(.not.isInDomain(pos)) inWhichDomain = -1
!    #endif

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
  !> @brief Computes div(V)
  !> @details Computes div(V) in all domain
  !> @param real [out] d :: div(V)

  subroutine divergence_V()
    use parameters
    use globals,    only : dx, dy, dz, primit, divV, coords, comm3d,   &
                           left, right, top, bottom, out, in, divV
    implicit none
    integer :: i,j,k
    integer, parameter :: nxm1=nx-1 ,nxp1=nx+1
    integer, parameter :: nym1=ny-1, nyp1=ny+1
    integer, parameter :: nzm1=nz-1, nzp1=nz+1
    integer:: status(MPI_STATUS_SIZE), err
    real, dimension(1,0:nyp1,0:nzp1)::sendr,recvr,sendl,recvl
    real, dimension(0:nxp1,1,0:nzp1)::sendt,recvt,sendb,recvb
    real, dimension(0:nxp1,0:nyp1,1)::sendi,recvi,sendo,recvo
    integer, parameter :: bxsize=(ny+2)*(nz+2)
    integer, parameter :: bysize=(nx+2)*(nz+2)
    integer, parameter :: bzsize=(nx+2)*(ny+2)

    do k =1,nz
      do j =1,ny
        do i =1,nx
          divV(i,j,k) =  (primit(2,i+1,j,k)-primit(2,i-1,j,k))/(2.*dx)  &
                       + (primit(3,i,j+1,k)-primit(3,i,j-1,k))/(2.*dy)  &
                       + (primit(4,i,j,k+1)-primit(4,i,j,k-1))/(2.*dz)
        end do
      end do
    end do

    !  pass ghost cells to 1st order
#ifdef MPIP

  !   Exchange boundaries between processors
  !   -------------------------------------------------------------

  !   boundaries to procs: right, left, top, bottom, in and out
  sendr(1,:,:)=divV(nx    ,0:nyp1,0:nzp1)
  sendl(1,:,:)=divV(1     ,0:nyp1,0:nzp1)
  sendt(:,1,:)=divV(0:nxp1,ny    ,0:nzp1)
  sendb(:,1,:)=divV(0:nxp1,1     ,0:nzp1)
  sendi(:,:,1)=divV(0:nxp1,0:nyp1,nz    )
  sendo(:,:,1)=divV(0:nxp1,0:nyp1,1     )

  call mpi_sendrecv(sendr, bxsize, mpi_real_kind, right, 0,           &
                   recvl, bxsize, mpi_real_kind, left,   0,           &
                   comm3d, status , err)

  call mpi_sendrecv(sendt, bysize, mpi_real_kind, top,   0,           &
                   recvb, bysize, mpi_real_kind, bottom ,0,           &
                   comm3d, status , err)

  call mpi_sendrecv(sendi, bzsize, mpi_real_kind, in,    0,           &
                   recvo, bzsize, mpi_real_kind, out,    0,           &
                   comm3d, status , err)

  call mpi_sendrecv(sendl, bxsize, mpi_real_kind, left,  0,           &
                   recvr, bxsize, mpi_real_kind, right,  0,           &
                   comm3d, status , err)

  call mpi_sendrecv(sendb, bysize, mpi_real_kind, bottom, 0,           &
                   recvt, bysize, mpi_real_kind, top,     0,           &
                   comm3d, status , err)

  call mpi_sendrecv(sendo, bzsize, mpi_real_kind, out,    0,           &
                   recvi, bzsize, mpi_real_kind, in,      0,           &
                   comm3d, status , err)

  if (left  .ne. -1) divV(0,     0:nyp1,0:nzp1)=recvl(1,:,:)
  if (right .ne. -1) divV(nxp1,  0:nyp1,0:nzp1)=recvr(1,:,:)
  if (bottom.ne. -1) divV(0:nxp1,0,     0:nzp1)=recvb(:,1,:)
  if (top   .ne. -1) divV(0:nxp1,nyp1,  0:nzp1)=recvt(:,1,:)
  if (out   .ne. -1) divV(0:nxp1,0:nyp1,0     )=recvo(:,:,1)
  if (in    .ne. -1) divV(0:nxp1,0:nyp1,nzp1  )=recvi(:,:,1)

#else

  !   periodic BCs
  if (bc_left == BC_PERIODIC .and. bc_right == BC_PERIODIC) then
    !   Left BC
    if (coords(0).eq.0) then
      divV(0,:,:)=divV(nx,:,:)
    end if
    !   Right BC
    if (coords(0).eq.MPI_NBX-1) then
      divV(nxp1,:,:)=divV(1,:,:)
    end if
  end if

  if (bc_bottom == BC_PERIODIC .and. bc_top == BC_PERIODIC) then
    !   bottom BC
    if (coords(1).eq.0) then
      divV(:,0,:)= divV(:,ny,:)
     end if
     !   top BC
     if (coords(1).eq.MPI_NBY-1) then
      divV(:,nyp1,:)= divV(:,1,:)
  endif

  if (bc_out == BC_PERIODIC .and. bc_in == BC_PERIODIC) then
    !   out BC
    if (coords(2).eq.0) then
      divV(:,:,0)= divV(:,:,nz)
    end if
    !   in BC
    if (coords(2).eq.MPI_NBZ-1) then
      divV(:,:,nzp1)= divV(:,:,1)
    end if
  end if

#endif

  !   Reflecting BCs
  !     left
  if (bc_left == BC_CLOSED) then
    if (coords(0).eq.0) then
      primit(2  ,nxmin,:,:) = -primit(2  ,nghost,:,:)
      primit(3:4,nxmin,:,:) =  primit(3:4,nghost,:,:)
      do k=0,nzp1
        do j=0,nyp1
          divV(0,j,k) = (primit(2,1,j,  k  )-primit(2,-1,j  ,k  ))/(2.*dx)  &
                      + (primit(3,0,j+1,k  )-primit(3, 0,j-1,k  ))/(2.*dy)  &
                      + (primit(4,0,j  ,k+1)-primit(4, 0,j,  k-1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   right
  if (bc_right == BC_CLOSED) then
    if (coords(0).eq.(MPI_NBX-1)) then
      primit(2  ,nxmax,:,:) = -primit(2  ,nx-1,:,:)
      primit(3:4,nxmax,:,:) =  primit(3:4,nx-1,:,:)
      do k=0,nzp1
        do j=0,nyp1
          divV(nxp1,j,k) = (primit(2,nxmax,j  ,k  )-primit(2,nx,  j,  k  ))/(2.*dx)  &
                         + (primit(3,nxp1, j+1,k  )-primit(3,nxp1,j-1,k  ))/(2.*dy)  &
                         + (primit(4,nxp1, j  ,k+1)-primit(4,nxp1,j,  k-1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   bottom
  if (bc_bottom == BC_CLOSED) then
    if (coords(1).eq.0) then
      primit(2,:,nymin,:) =  primit(2,:,nghost,:)
      primit(3,:,nymin,:) = -primit(3,:,nghost,:)
      primit(4,:,nymin,:) =  primit(4,:,nghost,:)
      do k=0,nzp1
        do i=0,nxp1
          divV(i,0,k) = (primit(2,i+1,0,k  )-primit(2,i-1,0, k  ))/(2.*dx)  &
                      + (primit(3,i  ,1,k  )-primit(3,i  ,-1,k  ))/(2.*dy)  &
                      + (primit(4,i  ,0,k+1)-primit(4,i  ,0, k-1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   top
  if (bc_top == BC_CLOSED) then
    if (coords(1).eq.(MPI_NBY-1)) then
      primit(2,:,nymax,:) =  primit(2,:,ny-1,:)
      primit(3,:,nymax,:) = -primit(3,:,ny-1,:)
      primit(4,:,nymax,:) =  primit(4,:,ny-1,:)
      do k=0,nzp1
        do i=0,nxp1
          divV(i,nyp1,k) = (primit(2,i+1,nyp1, k  )-primit(2,i-1,nyp1,  k  ))/(2.*dx)  &
                         + (primit(3,i  ,nymax,k  )-primit(3,i,  ny,    k  ))/(2.*dy)  &
                         + (primit(4,i  ,nyp1, k+1)-primit(4,i,  nyp1,  k-1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   out
  if (bc_out == BC_CLOSED) then
    if (coords(2).eq.0) then
      primit(2:3,:,:,nzmin) =  primit(2:3,:,:,nghost)
      primit(4  ,:,:,nzmin) = -primit(4  ,:,:,nghost)
      do j=0,nyp1
        do i=0,nxp1
          divV(i,j,0) = (primit(2,i+1,j,  0)-primit(2,i-1,j,   0))/(2.*dx)  &
                      + (primit(3,i,  j+1,0)-primit(3,i,  j-1, 0))/(2.*dy)  &
                      + (primit(4,i,  j,  1)-primit(4,i,  j,  -1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   in
  if (bc_in == BC_CLOSED) then
    if (coords(2).eq.MPI_NBZ-1) then
      primit(2:3,:,:,nzmax) =  primit(2:3,:,:,nz-1)
      primit(4  ,:,:,nzmax) = -primit(4  ,:,:,nz-1)
      do j=0,nyp1
        do i=0,nxp1
          divV(i,j,nzp1) = (primit(2,i+1,j,  nzp1 )-primit(2,i-1,j,  nzp1))/(2.*dx)  &
                         + (primit(3,i  ,j+1,nzp1 )-primit(3,i,  j-1,nzp1))/(2.*dy)  &
                         + (primit(4,i  ,j,  nzmax)-primit(4,i,  j,  nz  ))/(2.*dz)
        end do
      end do
    end if
  end if

  !   outflow BCs
  !   left
  if (bc_left == BC_OUTFLOW) then
    if (coords(0).eq.0) then
      primit(2:4,nxmin,:,:) =  primit(2:4,nghost,:,:)
      do k=0,nzp1
        do j=0,nyp1
          divV(0,j,k) = (primit(2,1,j,  k  )-primit(2,-1,j,  k  ))/(2.*dx)  &
                      + (primit(3,0,j+1,k  )-primit(3,0, j-1,k  ))/(2.*dy)  &
                      + (primit(4,0,j,  k+1)-primit(4,0, j  ,k-1))/(2.*dz)
        end do
      end do
     end if
  end if

  !   right
  if (bc_right == BC_OUTFLOW) then
    if (coords(0).eq.MPI_NBX-1) then
      primit(2:4,nxmax,:,:) =  primit(2:4,nx-1,:,:)
      do k=0,nzp1
        do j=0,nyp1
          divV(nxp1,j,k) = (primit(2,nxmax,j,  k  )-primit(2,nx,  j,  k  ))/(2.*dx)  &
                         + (primit(3,nxp1, j+1,k  )-primit(3,nxp1,j-1,k  ))/(2.*dy)  &
                         + (primit(4,nxp1, j,  k+1)-primit(4,nxp1,j,  k-1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   bottom
  if (bc_bottom == BC_OUTFLOW) then
    if (coords(1).eq.0) then
      primit(2:4,:,nymin,:) =  primit(2:4,:,nghost,:)
      do k=0,nzp1
        do i=0,nxp1
          divV(i,0,k) = (primit(2,i+1,0,k  )-primit(2,i-1, 0, k  ))/(2.*dx)  &
                      + (primit(3,i  ,1,k  )-primit(3,i  ,-1, k  ))/(2.*dy)  &
                      + (primit(4,i  ,0,k+1)-primit(4,i  , 0, k-1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   top
  if (bc_top == BC_OUTFLOW) then
    if (coords(1).eq.MPI_NBY-1) then
      primit(2:4,:,nymax,:) =  primit(2:4,:,ny-1,:)
      do k=0,nzp1
        do i=0,nxp1
          divV(i,nyp1,k) = (primit(2,i+1,nyp1, k  )-primit(2,i-1,nyp1,k  ))/(2.*dx)  &
                         + (primit(3,i,  nymax,k  )-primit(3,i  ,ny  ,k  ))/(2.*dy)  &
                         + (primit(4,i,  nyp1, k+1)-primit(4,i  ,nyp1,k-1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   out
  if (bc_out == BC_OUTFLOW) then
    if (coords(2).eq.0) then
      primit(2:4,:,:,nzmin) =  primit(2:4,:,:,nghost)
      do j=0,nyp1
        do i=0,nxp1
          divV(i,j,0) = (primit(2,i+1,j,  0)-primit(2,i-1,j,   0))/(2.*dx)  &
                      + (primit(3,i  ,j+1,0)-primit(3,i  ,j-1, 0))/(2.*dy)  &
                      + (primit(4,i  ,j,  1)-primit(4,i  ,j,  -1))/(2.*dz)
        end do
      end do
    end if
  end if

  !   in
  if (bc_in == BC_OUTFLOW) then
    if (coords(2).eq.MPI_NBZ-1) then
      primit(2:4,:,:,nzmax) =  primit(2:4,:,:,nz-1)
        do j=0,nyp1
        do i=0,nxp1
          divV(i,j,nzp1) = (primit(2,i+1,j,  nzp1 )-primit(2,i-1,j,  nzp1))/(2.*dx)  &
                         + (primit(3,i,  j+1,nzp1 )-primit(3,i  ,j-1,nzp1))/(2.*dy)  &
                         + (primit(4,i,  j,  nzmax)-primit(4,i,  j,  nz  ))/(2.*dz)
        end do
      end do
    end if
  end if


  end subroutine divergence_V

  !=======================================================================

end module utilities
