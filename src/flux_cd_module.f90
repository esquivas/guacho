!=======================================================================
!> @file flux_cd_module.f90
!> @brief Flux CD module
!> @author C. Villareal D'Angelo, M. Schneiter, A. Esquivel
!> @date 26/Apr/2016
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

!> @brief Module to computes the flux-CD div B correction
!> @details This module corrects the div B with a flux interpolated
!> central difference scheme
!> See. Sect. 4.5 of Toth 2000, Journal of Computational Physics 161, 605
!>
!> Changes to the E boundary treatment:
!>  - Only the tangential components of E are exchanged, and only over
!>    the interior range of each face (edges/corners of E are never read
!>    by flux_cd_update, and neither is the normal component).
!>  - The six exchanges are posted at once (irecv/isend + waitall), with
!>    tags that encode the direction of travel, so that runs with only
!>    two ranks along a periodic direction are safe.
!>  - Neighbours are compared against MPI_PROC_NULL instead of -1
!>    (MPI_PROC_NULL is -2 in Open MPI).
!>  - Physical boundaries are applied by a single routine, with a
!>    consistent parity at closed walls and a zero-gradient fallback for
!>    any boundary type not handled explicitly (e.g. BC_OTHER).
!>  - Serial branch: fixed the missing "end if" in the periodic BCs.

module flux_cd_module

#ifdef BFIELD

  implicit none
  real, allocatable :: e(:,:,:,:) !< electric field

  !> Parity of the tangential components of E at BC_CLOSED walls
  !> (the normal component gets the opposite parity).
  !>  +1. : tangential E even  <->  v_n odd, B_n even, B_t odd
  !>  -1. : tangential E odd   <->  v_n odd, B_n odd,  B_t even
  !>        (perfectly conducting wall)
  !> This MUST be consistent with the parity applied to B in the
  !> hydro boundaries for BC_CLOSED. Only the tangential sign affects
  !> the B update; the normal component in the ghost layer is unused.
  real, parameter :: closed_tang_sign = +1.

contains

  !=======================================================================
  !>@brief Boundary conditions (one cell) for flux-CD
  !>@details Fills the first ghost layer of E, used in the flux-CD update.
  !> At rank interfaces E must be the value computed by the neighbour
  !> in its own interior cell, otherwise div B is not preserved there.
  subroutine boundaryI_ef()

    use parameters
    use globals
    implicit none

    integer, parameter :: nxp1=nx+1, nyp1=ny+1, nzp1=nz+1

#ifdef MPIP

    integer, parameter :: bxsize=2*ny*nz
    integer, parameter :: bysize=2*nx*nz
    integer, parameter :: bzsize=2*nx*ny
    !  tags = direction of travel of the message
    integer, parameter :: tag_px=1, tag_mx=2, tag_py=3, tag_my=4,             &
                          tag_pz=5, tag_mz=6
    !  x faces carry (Ey,Ez), y faces (Ex,Ez), z faces (Ex,Ey)
    real, dimension(2,ny,nz) :: sendr, recvr, sendl, recvl
    real, dimension(2,nx,nz) :: sendt, recvt, sendb, recvb
    real, dimension(2,nx,ny) :: sendi, recvi, sendo, recvo
    integer :: req(12), stats(MPI_STATUS_SIZE,12), err

    !   post all receives
    call mpi_irecv(recvl, bxsize, mpi_real_kind, left  , tag_px, comm3d, req(1), err)
    call mpi_irecv(recvr, bxsize, mpi_real_kind, right , tag_mx, comm3d, req(2), err)
    call mpi_irecv(recvb, bysize, mpi_real_kind, bottom, tag_py, comm3d, req(3), err)
    call mpi_irecv(recvt, bysize, mpi_real_kind, top   , tag_my, comm3d, req(4), err)
    call mpi_irecv(recvo, bzsize, mpi_real_kind, out   , tag_pz, comm3d, req(5), err)
    call mpi_irecv(recvi, bzsize, mpi_real_kind, in    , tag_mz, comm3d, req(6), err)

    !   pack tangential components of the outermost interior layer
    sendr(:,:,:) = e(2:3, nx  , 1:ny, 1:nz)
    sendl(:,:,:) = e(2:3, 1   , 1:ny, 1:nz)

    sendt(1,:,:) = e(1  , 1:nx, ny  , 1:nz)
    sendt(2,:,:) = e(3  , 1:nx, ny  , 1:nz)
    sendb(1,:,:) = e(1  , 1:nx, 1   , 1:nz)
    sendb(2,:,:) = e(3  , 1:nx, 1   , 1:nz)

    sendi(:,:,:) = e(1:2, 1:nx, 1:ny, nz  )
    sendo(:,:,:) = e(1:2, 1:nx, 1:ny, 1   )

    call mpi_isend(sendr, bxsize, mpi_real_kind, right , tag_px, comm3d, req(7) , err)
    call mpi_isend(sendl, bxsize, mpi_real_kind, left  , tag_mx, comm3d, req(8) , err)
    call mpi_isend(sendt, bysize, mpi_real_kind, top   , tag_py, comm3d, req(9) , err)
    call mpi_isend(sendb, bysize, mpi_real_kind, bottom, tag_my, comm3d, req(10), err)
    call mpi_isend(sendi, bzsize, mpi_real_kind, in    , tag_pz, comm3d, req(11), err)
    call mpi_isend(sendo, bzsize, mpi_real_kind, out   , tag_mz, comm3d, req(12), err)

    call mpi_waitall(12, req, stats, err)

    !  unpack, receives from MPI_PROC_NULL (domain bounds) leave the buffer
    !  untouched
    if (left   /= MPI_PROC_NULL) e(2:3, 0   , 1:ny, 1:nz) = recvl(:,:,:)
    if (right  /= MPI_PROC_NULL) e(2:3, nxp1, 1:ny, 1:nz) = recvr(:,:,:)

    if (bottom /= MPI_PROC_NULL) then
      e(1, 1:nx, 0   , 1:nz) = recvb(1,:,:)
      e(3, 1:nx, 0   , 1:nz) = recvb(2,:,:)
    end if
    if (top    /= MPI_PROC_NULL) then
      e(1, 1:nx, nyp1, 1:nz) = recvt(1,:,:)
      e(3, 1:nx, nyp1, 1:nz) = recvt(2,:,:)
    end if

    if (out    /= MPI_PROC_NULL) e(1:2, 1:nx, 1:ny, 0   ) = recvo(:,:,:)
    if (in     /= MPI_PROC_NULL) e(1:2, 1:nx, 1:ny, nzp1) = recvi(:,:,:)

#else

    !   periodic BCs (single process)
    if (bc_left == BC_PERIODIC .and. bc_right == BC_PERIODIC) then
      e(:, 0   , 1:ny, 1:nz) = e(:, nx, 1:ny, 1:nz)
      e(:, nxp1, 1:ny, 1:nz) = e(:, 1 , 1:ny, 1:nz)
    end if

    if (bc_bottom == BC_PERIODIC .and. bc_top == BC_PERIODIC) then
      e(:, 1:nx, 0   , 1:nz) = e(:, 1:nx, ny, 1:nz)
      e(:, 1:nx, nyp1, 1:nz) = e(:, 1:nx, 1 , 1:nz)
    end if

    if (bc_out == BC_PERIODIC .and. bc_in == BC_PERIODIC) then
      e(:, 1:nx, 1:ny, 0   ) = e(:, 1:nx, 1:ny, nz)
      e(:, 1:nx, 1:ny, nzp1) = e(:, 1:nx, 1:ny, 1 )
    end if

#endif

    !   physical (non-periodic) boundaries
    if (coords(0) == 0        ) call phys_bc_face(1, 0   , 1 , bc_left  )
    if (coords(0) == MPI_NBX-1) call phys_bc_face(1, nxp1, nx, bc_right )
    if (coords(1) == 0        ) call phys_bc_face(2, 0   , 1 , bc_bottom)
    if (coords(1) == MPI_NBY-1) call phys_bc_face(2, nyp1, ny, bc_top   )
    if (coords(2) == 0        ) call phys_bc_face(3, 0   , 1 , bc_out   )
    if (coords(2) == MPI_NBZ-1) call phys_bc_face(3, nzp1, nz, bc_in    )

  end subroutine boundaryI_ef

  !=======================================================================
  !>@brief Fills one ghost layer of E at a physical boundary
  !>@param integer [in] axis : 1, 2 or 3 (normal direction of the face)
  !>@param integer [in] ig   : index of the ghost layer along axis
  !>@param integer [in] ii   : index of the adjacent interior layer
  !>@param integer [in] bc   : boundary type of this face
  subroutine phys_bc_face(axis, ig, ii, bc)

    use parameters
    use globals
    implicit none
    integer, intent(in) :: axis, ig, ii, bc
    real :: sn, st   ! signs for normal and tangential components

    select case (bc)
    case (BC_PERIODIC)
      return                        ! filled by the exchange/periodic copy
    case (BC_CLOSED)
      st =  closed_tang_sign
      sn = -closed_tang_sign
    case default
      !  BC_OUTFLOW, and fallback (zero gradient) for any other type,
      !  so that the ghost layer never keeps a stale value
      st = 1.
      sn = 1.
    end select

    select case (axis)
    case (1)
      e(1  , ig  , 1:ny, 1:nz) = sn * e(1  , ii  , 1:ny, 1:nz)
      e(2:3, ig  , 1:ny, 1:nz) = st * e(2:3, ii  , 1:ny, 1:nz)
    case (2)
      e(2  , 1:nx, ig  , 1:nz) = sn * e(2  , 1:nx, ii  , 1:nz)
      e(1  , 1:nx, ig  , 1:nz) = st * e(1  , 1:nx, ii  , 1:nz)
      e(3  , 1:nx, ig  , 1:nz) = st * e(3  , 1:nx, ii  , 1:nz)
    case (3)
      e(3  , 1:nx, 1:ny, ig  ) = sn * e(3  , 1:nx, 1:ny, ii  )
      e(1:2, 1:nx, 1:ny, ig  ) = st * e(1:2, 1:nx, 1:ny, ii  )
    end select

  end subroutine phys_bc_face

  !=======================================================================
  !>@brief Computes E
  !>@details Obtains the electric field from the fluxes
  !> (eq. 31 of Toth 2000). Requires f, g, h on faces 0:nx, 0:ny, 0:nz.
  subroutine get_efield()

    use parameters, only : nx, ny, nz
    use globals, only :  f, g, h
    implicit none
    integer :: i, j, k

    do k=1,nz
      do j=1,ny
        do i=1,nx

          e(1,i,j,k)=0.25*( -g(8,i,j-1,k) - g(8,i,j,k)                         &
                            +h(7,i,j,k-1) + h(7,i,j,k) )

          e(2,i,j,k)=0.25*( +f(8,i-1,j,k) + f(8,i,j,k)                         &
                            -h(6,i,j,k-1) - h(6,i,j,k) )

          e(3,i,j,k)=0.25*( -f(7,i-1,j,k) -f(7,i,j,k)                          &
                            +g(6,i,j-1,k) +g(6,i,j,k) )

        end do
      end do
    end do

    call boundaryI_ef()

  end subroutine get_efield

  !=======================================================================
  !> @brief Upper level wrapper for flux-CD update
  !> @details Upper level wrapper for flux-CD, updates the
  !> hydro variables with upwind scheme and the field as flux-CD
  !> @param integer [in] i : cell index in the X direction
  !> @param integer [in] j : cell index in the Y direction
  !> @param integer [in] k : cell index in the Z direction
  !> @param real [in] dt : timestep
  subroutine flux_cd_update(i,j,k,dt)

    use parameters, only : passives, neqdyn
    use globals, only : dx, dy, dz, up, u, f, g, h
    implicit none
    integer, intent(in)  :: i, j, k
    real, intent (in)    :: dt
    real :: dtdx, dtdy, dtdz

    dtdx=dt/dx
    dtdy=dt/dy
    dtdz=dt/dz

    !   hydro variables (and passive eqs.)
    up(:5,i,j,k)=u(:5,i,j,k)-dtdx*(f(:5,i,j,k)-f(:5,i-1,j,k))                  &
                            -dtdy*(g(:5,i,j,k)-g(:5,i,j-1,k))                  &
                            -dtdz*(h(:5,i,j,k)-h(:5,i,j,k-1))

#ifdef PASSIVES
    if (passives) &
      up(neqdyn+1:,i,j,k) = u(neqdyn+1:,i,j,k)                                 &
                          - dtdx*(f(neqdyn+1:,i,j,k)-f(neqdyn+1:,i-1,j,k))     &
                          - dtdy*(g(neqdyn+1:,i,j,k)-g(neqdyn+1:,i,j-1,k))     &
                          - dtdz*(h(neqdyn+1:,i,j,k)-h(neqdyn+1:,i,j,k-1))
#endif
    ! evolution of B with flux-CD
    ! (uses only tangential E in the first ghost layer)
    up(6,i,j,k) = u(6,i,j,k)                                                   &
                - 0.5*dtdy*(e(3,i,j+1,k)-e(3,i,j-1,k))                         &
                + 0.5*dtdz*(e(2,i,j,k+1)-e(2,i,j,k-1))

   up(7,i,j,k) = u(7,i,j,k)                                                    &
               + 0.5*dtdx*(e(3,i+1,j,k)-e(3,i-1,j,k))                          &
               - 0.5*dtdz*(e(1,i,j,k+1)-e(1,i,j,k-1))

   up(8,i,j,k) = u(8,i,j,k)                                                    &
               - 0.5*dtdx*(e(2,i+1,j,k)-e(2,i-1,j,k))                          &
               + 0.5*dtdy*(e(1,i,j+1,k)-e(1,i,j-1,k))

  end subroutine flux_cd_update

  !=======================================================================

#endif

end module flux_cd_module
