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
  real, parameter   :: four_pi_G = 4.*pi*Ggrav

contains

  !================================================================
  !> @brief Initialization of module
  !> @details Allocates memory for all global variables that correspond
  !> to the module
  subroutine init_self_gravity()
    use parameters, only : nx, ny, nz
    implicit none
    !  allocate one ghost cell (needed for gradient/laplacian)
    allocate( phi_grav(0:nx+1,0:ny+1,0:nz+1) )

  end subroutine init_self_gravity

  !=======================================================================
  !>@brief Boundary conditions (one cell) for flux-CD
  !>@details Boundary conditions applied to E, used
  !> in the flux-CD calculation
  subroutine phi_grav_boundaries()
    use parameters
    use globals
    implicit none

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

#ifdef MPIP

    !   Exchange boundaries between processors
    !   -------------------------------------------------------------

    !   boundaries to procs: right, left, top, bottom, in and out
    sendr(1,:,:)=phi_grav(nx    ,0:nyp1,0:nzp1)
    sendl(1,:,:)=phi_grav(1     ,0:nyp1,0:nzp1)
    sendt(:,1,:)=phi_grav(0:nxp1,ny    ,0:nzp1)
    sendb(:,1,:)=phi_grav(0:nxp1,1     ,0:nzp1)
    sendi(:,:,1)=phi_grav(0:nxp1,0:nyp1,nz    )
    sendo(:,:,1)=phi_grav(0:nxp1,0:nyp1,1     )

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

    if (left  .ne. -1) phi_grav(0     ,0:nyp1,0:nzp1)=recvl(1,:,:)
    if (right .ne. -1) phi_grav(nxp1  ,0:nyp1,0:nzp1)=recvr(1,:,:)
    if (bottom.ne. -1) phi_grav(0:nxp1,0     ,0:nzp1)=recvb(:,1,:)
    if (top   .ne. -1) phi_grav(0:nxp1,nyp1  ,0:nzp1)=recvt(:,1,:)
    if (out   .ne. -1) phi_grav(0:nxp1,0:nyp1,0     )=recvo(:,:,1)
    if (in    .ne. -1) phi_grav(0:nxp1,0:nyp1,nzp1  )=recvi(:,:,1)

#else

    !   periodic BCs
    if (bc_left == BC_PERIODIC .and. bc_right == BC_PERIODIC) then
      !   Left BC
      if (coords(0).eq.0) then
        phi_grav(0,:,:)=phi_grav(nx,:,:)
      end if
      !   Right BC
      if (coords(0).eq.MPI_NBX-1) then
        phi_grav(nxp1,:,:)=phi_grav(1,:,:)
      end if
    end if

    if ( bc_bottom == BC_PERIODIC .and. bc_top == BC_PERIODIC) then
      !   bottom BC
      if (coords(1).eq.0) then
        phi_grav(:,0,:)= phi_grav(:,ny,:)
      end if
      !   top BC
      if (coords(1).eq.MPI_NBY-1) then
        phi_grav(:,nyp1,:)= phi_grav(:,1,:)
      end if

      if (bc_out == BC_PERIODIC .and. bc_in == BC_PERIODIC) then
        !   out BC
        if (coords(2).eq.0) then
          phi_grav(:,:,0)= phi_grav(:,:,nz)
        end if
        !   in BC
        if (coords(2).eq.MPI_NBZ-1) then
          phi_grav(:,:,nzp1)= phi_grav(:,:,1)
        end if
      end if

#endif

    !   outflow BCs
    !   left
    if (bc_left == BC_OUTFLOW) then
      if (coords(0).eq.0) then
        phi_grav(0,0:nyp1,0:nzp1)=phi_grav(1 ,0:nyp1,0:nzp1)
      end if
    end if

    !   right
    if (bc_right == BC_OUTFLOW) then
      if (coords(0).eq.MPI_NBX-1) then
        phi_grav(nxp1,0:nyp1,0:nzp1)=phi_grav(nx,0:nyp1,0:nzp1)
      end if
    end if

    !   bottom
    if (bc_bottom == BC_OUTFLOW) then
      if (coords(1).eq.0) then
        phi_grav(0:nxp1,0,0:nzp1)=phi_grav(0:nxp1,1 ,0:nzp1)
      end if
    end if

    !   top
    if (bc_top == BC_OUTFLOW) then
      if (coords(1).eq.MPI_NBY-1) then
        phi_grav(0:nxp1,nyp1,0:nzp1)=phi_grav(0:nxp1,ny,0:nzp1)
      end if
    end if

    !   out
    if (bc_out == BC_OUTFLOW) then
      if (coords(2).eq.0) then
        phi_grav(0:nxp1,0:nyp1,0)=phi_grav(0:nxp1,0:nyp1,1 )
      end if
    end if

    !   in
    if (bc_in == BC_OUTFLOW) then
      if (coords(2).eq.MPI_NBZ-1) then
        phi_grav(0:nxp1,0:nyp1,nzp1)=phi_grav(0:nxp1,0:nyp1,nz)
      end if
    end if

    print*, 'pase fronteras'

  end subroutine phi_grav_boundaries


  ! ======================================================================
  ! Computes the residue $\chi= \nabla^2 \phi - 4\pi G \rho$
  ! at a given cell
  subroutine get_residue(i,j,k,residue)
    use globals, only : dx, dy, dz, primit
    implicit none
    integer, intent(in)  :: i,j,k
    real,    intent(out) :: residue

    residue= ( phi_grav(i+1, j,  k   )+phi_grav(i-1,j,  k  ) -2.*phi_grav(i,j,k) )/dx**2      &
            +( phi_grav(i  , j+1,k   )+phi_grav(i  ,j-1,k  ) -2.*phi_grav(i,j,k) )/dy**2      &
            +( phi_grav(i  , j , k+1 )+phi_grav(i  ,j  ,k-1) -2.*phi_grav(i,j,k) )/dz**2      &
            - four_pi_G*primit(1,i,j,k)

  end subroutine get_residue

  !================================================================
  !> @brief Solve Poisson equation
  !> @details Compute the gravitational potential with a SOR mehtod
  subroutine solve_poisson()
    use parameters, only : nx, ny, nz
    use globals,    only : dx, dy, dz, primit, rank
    implicit none
    integer, parameter :: max_iterations=200
    real, parameter    :: Tol = 1E-4  !  Relative error tolerance
    real               :: omega, relative_error
    real    :: residue, max_error, e_ijk, ph0
    logical :: need_more=.false.
    integer :: iter, ipass, i,j,k, ksw, jsw, kpass

    omega=1.99

    e_ijk = -2.*( 1./dx**2 +1./dy**2 + 1./dz**2 )

    main_loop : do iter=1, max_iterations
      !
      ! do k = 1,nz
      !   do j = 1, ny
      !     do i = 1, nx
      !       call get_residuphi_grav(i,j,k,residue)
      !       !relative_error= omega*abs(residue)/abs(phi(i,j,k)*e_ijk)
      !       ph0         = phi_grav(i,j,k)
      !       phip(i,j,k) = ph0 - omega*residue/e_ijk
      !
      !       !call get_phi_star(rho,phi,i,j,k,dx,residue)
      !       !phi(i,j,k)=omega*residue + (1.-omega)*ph0
      !
      !       relative_error= abs(phip(i,j,k)-ph0)/abs(ph0)
      !       max_error     = max(max_error,relative_error)
      !
      !       !if(relative_error > Tol)
      !       need_more=.true.
      !     end do
      !   end do
      ! end do

      ksw=1

      black_red: do kpass=1,2

        max_error=-10.
        do ipass=1,2
          jsw=ksw
          do k=1, nz
            do j=jsw,ny,2
              do i=ipass,nx,2

                call get_residue(i,j,k,residue)
                !relative_error= omega*abs(residue)/abs(phi(i,j,k)*e_ijk)
                ph0             = phi_grav(i,j,k)
                phi_grav(i,j,k) = ph0 - omega*residue/e_ijk

                !call get_phi_star(rho,phi,i,j,k,dx,residue)
                !ph0=phi(i,j,k)
                !phi(i,j,k)=omega*residue + (1.-omega)*ph0

                relative_error= abs(phi_grav(i,j,k)-ph0)/abs(ph0)
                max_error     = max(max_error,relative_error)

                !if(relative_error > Tol)
                need_more=.true.

              end do
            end do
            jsw=2-jsw
          end do
          ksw=2-ksw
        end do

      end do black_red

      call phi_grav_boundaries()

      !print*,rank, 'errror:',max_error

      if(.not.need_more) then
        print*, 'Converged in ', iter, 'iterations',max_error
        return
      end if
      !  reset convergence flag
      need_more=.false.

    end do main_loop

    print'(a)', 'SOR exceeded maximum number of iterations'


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
