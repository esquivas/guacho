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

module flux_cd_module

#ifdef BFIELD

  implicit none
  real, allocatable :: e(:,:,:,:) !< electric field

contains

  !=======================================================================
  !>@brief Boundary conditions (one cell) for flux-CD
  !>@details Boundary conditions applied to E, used
  !> in the flux-CD calculation
  subroutine boundaryI_ef()

    use parameters
    use globals
    implicit none

    integer, parameter :: nxm1=nx-1 ,nxp1=nx+1
    integer, parameter :: nym1=ny-1, nyp1=ny+1
    integer, parameter :: nzm1=nz-1, nzp1=nz+1
    integer:: status(MPI_STATUS_SIZE), err
    real, dimension(3,1,0:nyp1,0:nzp1)::sendr,recvr,sendl,recvl
    real, dimension(3,0:nxp1,1,0:nzp1)::sendt,recvt,sendb,recvb
    real, dimension(3,0:nxp1,0:nyp1,1)::sendi,recvi,sendo,recvo
    integer, parameter :: bxsize=3*(ny+2)*(nz+2)
    integer, parameter :: bysize=3*(nx+2)*(nz+2)
    integer, parameter :: bzsize=3*(nx+2)*(ny+2)

#ifdef MPIP

    !   Exchange boundaries between processors
    !   -------------------------------------------------------------

    !   boundaries to procs: right, left, top, bottom, in and out
    sendr(:,1,:,:)=e(:,nx    ,0:nyp1,0:nzp1)
    sendl(:,1,:,:)=e(:,1     ,0:nyp1,0:nzp1)
    sendt(:,:,1,:)=e(:,0:nxp1,ny    ,0:nzp1)
    sendb(:,:,1,:)=e(:,0:nxp1,1     ,0:nzp1)
    sendi(:,:,:,1)=e(:,0:nxp1,0:nyp1,nz    )
    sendo(:,:,:,1)=e(:,0:nxp1,0:nyp1,1     )

    call mpi_sendrecv(sendr, bxsize, mpi_real_kind, right  ,0,                 &
                      recvl, bxsize, mpi_real_kind, left   ,0,                 &
                      comm3d, status , err)

    call mpi_sendrecv(sendt, bysize, mpi_real_kind, top    ,0,                 &
                      recvb, bysize, mpi_real_kind, bottom ,0,                 &
                      comm3d, status , err)

    call mpi_sendrecv(sendi, bzsize, mpi_real_kind, in     ,0,                 &
                      recvo, bzsize, mpi_real_kind, out    ,0,                 &
                      comm3d, status , err)

    call mpi_sendrecv(sendl, bxsize, mpi_real_kind, left  , 0,                 &
                      recvr, bxsize, mpi_real_kind, right , 0,                 &
                      comm3d, status , err)

    call mpi_sendrecv(sendb, bysize, mpi_real_kind, bottom, 0,                 &
                      recvt, bysize, mpi_real_kind, top   , 0,                 &
                      comm3d, status , err)

    call mpi_sendrecv(sendo, bzsize, mpi_real_kind, out   , 0,                 &
                      recvi, bzsize, mpi_real_kind, in    , 0,                 &
                      comm3d, status , err)

    if (left  .ne. -1) e(:,0     ,0:nyp1,0:nzp1)=recvl(:,1,:,:)
    if (right .ne. -1) e(:,nxp1  ,0:nyp1,0:nzp1)=recvr(:,1,:,:)
    if (bottom.ne. -1) e(:,0:nxp1,0     ,0:nzp1)=recvb(:,:,1,:)
    if (top   .ne. -1) e(:,0:nxp1,nyp1  ,0:nzp1)=recvt(:,:,1,:)
    if (out   .ne. -1) e(:,0:nxp1,0:nyp1,0     )=recvo(:,:,:,1)
    if (in    .ne. -1) e(:,0:nxp1,0:nyp1,nzp1  )=recvi(:,:,:,1)

#else

    !   periodic BCs
    if (bc_left == BC_PERIODIC .and. bc_right == BC_PERIODIC) then
      !   Left BC
      if (coords(0).eq.0) then
        e(:,0,:,:)=e(:,nx,:,:)
      end if
      !   Right BC
      if (coords(0).eq.MPI_NBX-1) then
        e(:,nxp1,:,:)=e(:,1,:,:)
      end if
    end if

    if ( bc_bottom == BC_PERIODIC .and. bc_top == BC_PERIODIC) then
      !   bottom BC
      if (coords(1).eq.0) then
        e(:,:,0,:)= e(:,:,ny,:)
      end if
      !   top BC
      if (coords(1).eq.MPI_NBY-1) then
        e(:,:,nyp1,:)= e(:,:,1,:)
      end if

      if (bc_out == BC_PERIODIC .and. bc_in == BC_PERIODIC) then
        !   out BC
        if (coords(2).eq.0) then
          e(:,:,:,0)= e(:,:,:,nz)
        end if
        !   in BC
        if (coords(2).eq.MPI_NBZ-1) then
          e(:,:,:,nzp1)= e(:,:,:,1)
        end if
      end if

#endif

    !   Reflecting BCs (not tested)
    !     left
    if (bc_left == BC_CLOSED) then
      if (coords(0).eq.0) then
        e(1  ,0,0:nyp1,0:nzp1) =-e(1  ,1,0:nyp1,0:nzp1)
        e(2:3,0,0:nyp1,0:nzp1) = e(2:3,1,0:nyp1,0:nzp1)
      end if
    end if

    !   right
    if (bc_right == BC_CLOSED) then
      if (coords(0).eq.(MPI_NBX-1)) then
        e(1  ,nxp1,0:nyp1,0:nzp1) =-e(1  ,nx,0:nyp1,0:nzp1)
        e(2:3,nxp1,0:nyp1,0:nzp1) = e(2:3,nx,0:nyp1,0:nzp1)
      end if
    end if

    !   bottom
    if (bc_bottom == BC_CLOSED) then
      if (coords(1).eq.0) then
        e(1,0:nxp1,0,0:nzp1) = e(1,0:nxp1,1,0:nzp1)
        e(2,0:nxp1,0,0:nzp1) =-e(2,0:nxp1,1,0:nzp1)
        e(3,0:nxp1,0,0:nzp1) = e(3,0:nxp1,1,0:nzp1)
      end if
    end if

    !   top
    if (bc_top == BC_CLOSED) then
      if (coords(1).eq.(MPI_NBY-1)) then
        e(1,0:nxp1,nyp1,0:nzp1) = e(1,0:nxp1,ny,0:nzp1)
        e(2,0:nxp1,nyp1,0:nzp1) =-e(2,0:nxp1,ny,0:nzp1)
        e(3,0:nxp1,nyp1,0:nzp1) = e(3,0:nxp1,ny,0:nzp1)
      end if
    end if

    !   out
    if (bc_out == BC_CLOSED) then
      if (coords(2).eq.0) then
        e(1:2,0:nxp1,0:nyp1,0) = e(1:2,0:nxp1,0:nyp1,1)
        e(3  ,0:nxp1,0:nyp1,0) = e(3  ,0:nxp1,0:nyp1,1)
      end if
    end if

    !   in
    if (bc_in == BC_CLOSED) then
      if (coords(2).eq.MPI_NBZ-1) then
        e(1:2,0:nxp1,0:nyp1,nzp1) = e(1:2,0:nxp1,0:nyp1,nz)
        e(3  ,0:nxp1,0:nyp1,nzp1) = e(3  ,0:nxp1,0:nyp1,nz)
      end if
    end if

    !   outflow BCs
    !   left
    if (bc_left == BC_OUTFLOW) then
      if (coords(0).eq.0) then
        e(:,0,0:nyp1,0:nzp1)=e(:,1 ,0:nyp1,0:nzp1)
      end if
    end if

    !   right
    if (bc_right == BC_OUTFLOW) then
      if (coords(0).eq.MPI_NBX-1) then
        e(:,nxp1,0:nyp1,0:nzp1)=e(:,nx,0:nyp1,0:nzp1)
      end if
    end if

    !   bottom
    if (bc_bottom == BC_OUTFLOW) then
      if (coords(1).eq.0) then
        e(:,0:nxp1,0,0:nzp1)=e(:,0:nxp1,1 ,0:nzp1)
      end if
    end if

    !   top
    if (bc_top == BC_OUTFLOW) then
      if (coords(1).eq.MPI_NBY-1) then
        e(:,0:nxp1,nyp1,0:nzp1)=e(:,0:nxp1,ny,0:nzp1)
      end if
    end if

    !   out
    if (bc_out == BC_OUTFLOW) then
      if (coords(2).eq.0) then
        e(:,0:nxp1,0:nyp1,0)=e(:,0:nxp1,0:nyp1,1 )
      end if
    end if

    !   in
    if (bc_in == BC_OUTFLOW) then
      if (coords(2).eq.MPI_NBZ-1) then
        e(:,0:nxp1,0:nyp1,nzp1)=e(:,0:nxp1,0:nyp1,nz)
      end if
    end if

  end subroutine boundaryI_ef

  !=======================================================================
  !>@brief Computes E
  !>@details Obtains the electric field from the fluxes
  !> (eq. 31 of Toth 2000)
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
