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
  !>@brief Boundary conditions (one cell) for SOR itrations
  !>@details Boundary conditions applied to phi_grav, used in the SOR method to
  !>to solve the Poisson eq. requiered to add self-gravity.
  !>@nAt the moment works only with open and perdiodic boundaries
  subroutine phi_grav_boundaries()
    use parameters
    use globals
    implicit none

    integer, parameter :: nxm1=nx-1 ,nxp1=nx+1
    integer, parameter :: nym1=ny-1, nyp1=ny+1
    integer, parameter :: nzm1=nz-1, nzp1=nz+1
    integer:: status(MPI_STATUS_SIZE), err
    integer, parameter :: bxsize=(ny+2)*(nz+2)
    integer, parameter :: bysize=(nx+2)*(nz+2)
    integer, parameter :: bzsize=(nx+2)*(ny+2)

#ifdef MPIP

    !   Exchange boundaries between processors
    !   -------------------------------------------------------------
    call mpi_sendrecv(phi_grav(nx, :, :), bxsize, mpi_real_kind, right  ,0,    &
                      phi_grav( 0, :, :), bxsize, mpi_real_kind, left   ,0,    &
                      comm3d, status , err)

    call mpi_sendrecv(phi_grav( :,ny, :), bysize, mpi_real_kind, top    ,1,    &
                      phi_grav( :, 0, :), bysize, mpi_real_kind, bottom ,1,    &
                      comm3d, status , err)

    call mpi_sendrecv(phi_grav( :, :,nz), bzsize, mpi_real_kind, in     ,2,    &
                      phi_grav( :, :, 0), bzsize, mpi_real_kind, out    ,2,    &
                      comm3d, status , err)

    call mpi_sendrecv(phi_grav(   1, :, :), bxsize, mpi_real_kind, left  , 3,  &
                      phi_grav(nxp1, :, :), bxsize, mpi_real_kind, right , 3,  &
                      comm3d, status , err)

    call mpi_sendrecv(phi_grav( :,   1, :), bysize, mpi_real_kind, bottom, 4,  &
                      phi_grav( :,nyp1, :), bysize, mpi_real_kind, top   , 4,  &
                      comm3d, status , err)

    call mpi_sendrecv(phi_grav( :, :,    1), bzsize, mpi_real_kind, out  , 5,  &
                      phi_grav( :, :, nzp1), bzsize, mpi_real_kind, in   , 5,  &
                      comm3d, status , err)

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
        phi_grav(0,:,:)=phi_grav(1,:,:)
      end if
    end if

    !   right
    if (bc_right == BC_OUTFLOW) then
      if (coords(0).eq.MPI_NBX-1) then
        phi_grav(nxp1,:,:)=phi_grav(nx,:,:)
      end if
    end if

    !   bottom
    if (bc_bottom == BC_OUTFLOW) then
      if (coords(1).eq.0) then
        phi_grav(:,0,:)=phi_grav(:,1 ,:)
      end if
    end if

    !   top
    if (bc_top == BC_OUTFLOW) then
      if (coords(1).eq.MPI_NBY-1) then
        phi_grav(:,nyp1,:)=phi_grav(:,ny,:)
      end if
    end if

    !   out
    if (bc_out == BC_OUTFLOW) then
      if (coords(2).eq.0) then
        phi_grav(:,:,0)=phi_grav(:,:,1 )
      end if
    end if

    !   in
    if (bc_in == BC_OUTFLOW) then
      if (coords(2).eq.MPI_NBZ-1) then
        phi_grav(:,:,nzp1)=phi_grav(:,:,nz)
      end if
    end if

  end subroutine phi_grav_boundaries


  ! ======================================================================
  !> @brief Computes residue in a given cell
  !> @details Computes the residue at a given cell from the density taken from
  !> the global variable in primit
  !> @param integer [in] i : index in the x direction
  !> @param integer [in] j : index in the y direction
  !> @param integer [in] k : index in the z direction
  !> @param real [out] xi  : residue $\xi= \nabla^2 \phi - 4\pi G \rho$
  subroutine get_residue(i,j,k,xi)
    use globals, only : dx, dy, dz, primit
    implicit none
    integer, intent(in)  :: i,j,k
    real,    intent(out) :: xi

    xi = ( phi_grav(i+1,j,k) + phi_grav(i-1,j,k) - 2.*phi_grav(i,j,k) )/dx**2  &
       + ( phi_grav(i,j+1,k) + phi_grav(i,j-1,k) - 2.*phi_grav(i,j,k) )/dy**2  &
       + ( phi_grav(i,j,k+1) + phi_grav(i,j,k-1) - 2.*phi_grav(i,j,k) )/dz**2  &
       - four_pi_G*primit(1,i,j,k)

  end subroutine get_residue

  !================================================================
  !> @brief Solve Poisson equation
  !> @details Compute the gravitational potential with a SOR mehtod
  subroutine solve_poisson()
    use parameters, only : nx, ny, nz, mpi_real_kind, master
#ifdef MPIP
    use mpi
#endif
    use globals,    only : dx, dy, dz, primit, comm3d, rank,time
    use constants,  only : pi
    implicit none
    integer, parameter :: max_iterations=10000
    real, parameter    :: Tol = 1E-4   !  Relative error tolerance
    real               :: omega, relative_error
    real               :: xi , max_error, e_ijk, ph0, max_error_local
    logical, parameter :: enable_chebyshev_accel = .true.
    logical            :: converged
    integer            :: iter, err
    integer            :: i_rb, isw, jsw, i, j,k, ksw, kpass, ipass
    real, allocatable  :: phiP(:,:,:)


    converged = .false.
    e_ijk = ( 2.0/dx**2 + 2.0/dy**2 + 2.0/dz**2 )

    if(.not.enable_chebyshev_accel) allocate( phiP(0:nx+1,0:ny+1,0:nz+1) )

    main_loop : do iter=1, max_iterations

      omega = 1.5
      !omega = 2.0 / ( 1.0 + sin( pi/real(iter+1) ) )

      max_error_local = -1.0
      max_error       = 1e20

      if( enable_chebyshev_accel ) then

      black_red: do i_rb=1,2

    		jsw = i_rb
    		do k=1,nz
    			isw =jsw
    			do j=1,ny
    				do i=isw,nx,2

              call get_residue(i,j,k,xi)

              !ph0             = phi_grav(i,j,k)
              phi_grav(i,j,k) = xi/e_ijk !ph0 - omega*xi/e_ijk

              !relative_error  = abs(phi_grav(i,j,k)-ph0)/abs(ph0+1e-30)
              relative_error  = abs( omega*xi/e_ijk ) /                        &
                                max( abs( phi_grav(i,j,k) ),  1e-30 )

              max_error_local = max(max_error_local,relative_error)

            end do
            isw = 3 - isw
          end do
          jsw = 3 - jsw
        end do
      end do black_red

        !  this shouldn't work but the code commented above should, and it
        !  doesn't hav to check this later
        ! ksw=1
        ! black_red: do kpass=1,2
        !
        !   do ipass=1,2
        !     jsw=ksw
        !     do k=1,nz
        !       do j=jsw,ny,2
        !         do i=ipass,nx,2
        !
        !           call get_residue(i,j,k,xi)
        !
        !           ph0             = phi_grav(i,j,k)
        !           phi_grav(i,j,k) = ph0 - omega*xi/e_ijk
        !
        !           relative_error  = abs(phi_grav(i,j,k)-ph0)/abs(ph0 + 1e-30)
        !           max_error_local = max(max_error_local,relative_error)
        !
        !         end do
        !       end do
        !       jsw=3-jsw
        !     end do
        !     ksw=3-ksw
        !   end do
        !
        ! end do black_red
      else

        do k=1,nz
          do j=1,ny
            do i=1,nx

              call get_residue(i,j,k,xi)

              phiP(i,j,k)     = phi_grav(i,j,k) - omega*xi/e_ijk

              relative_error  = abs( omega*xi/e_ijk ) /                        &
                              max(  abs(phi_grav(i,j,k)),  1e-30 )

              max_error_local = max(max_error_local,relative_error)

            end do
          end do
        end do

        phi_grav(:,:,:) = phiP(:,:,:)

      end if

      !  need to share the error among the different cores
      call mpi_allreduce(max_error_local, max_error, 1, mpi_real_kind, mpi_max,&
                         comm3d, err)

      if(max_error < Tol) converged = .true.

      print*, max_error_local, max_error, omega

      call phi_grav_boundaries()

      if(converged) then
        if (rank == master) print'(a,i0,a,2es12.5)', 'SOR converged in ', iter,&
                                    ' iterations with an error of', max_error, &
                                    max_error_local
        if (.not.enable_chebyshev_accel) deallocate (phiP)
        return
      end if

    end do main_loop

    print'(a)', 'SOR exceeded maximum number of iterations'
    if (.not.enable_chebyshev_accel) deallocate(phiP)

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
