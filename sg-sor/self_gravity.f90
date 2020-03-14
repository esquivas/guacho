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
    integer, parameter :: max_iterations=100
    real, parameter    :: Tol = 1E-4  !  Relative error tolerance
    real               :: omega, relative_error
    real    :: residue, max_error, e_ijk, ph0
    logical :: need_more=.false.
    integer :: iter, ipass, i,j,k, ksw, jsw, kpass
    real    :: phip(0:nx+1,0:ny+1,0:nz+1)

    omega=1.99

    e_ijk = -2.*( 1./dx**2 +1./dy**2 + 1./dz**2 )

    main_loop : do iter=1, max_iterations

      do k = 1,nz
        do j = 1, ny
          do i = 1, nx
            call get_residue(i,j,k,residue)
            !relative_error= omega*abs(residue)/abs(phi(i,j,k)*e_ijk)
            ph0             = phi_grav(i,j,k)
            phip(i,j,k) = ph0 - omega*residue/e_ijk

            !call get_phi_star(rho,phi,i,j,k,dx,residue)
            !ph0=phi(i,j,k)
            !phi(i,j,k)=omega*residue + (1.-omega)*ph0

            relative_error= abs(phi_grav(i,j,k)-ph0)/abs(ph0)
            max_error     = max(max_error,relative_error)

            !if(relative_error > Tol)
            need_more=.true.
          end do
        end do
      end do

      ! ksw=1

      ! black_red: do kpass=1,2
      !
      !   max_error=-10.
      !   do ipass=1,2
      !     jsw=ksw
      !     do k=1, nz
      !       do j=jsw,ny,2
      !         do i=ipass,nx,2
      !
      !           call get_residue(i,j,k,residue)
      !           !relative_error= omega*abs(residue)/abs(phi(i,j,k)*e_ijk)
      !           ph0             = phi_grav(i,j,k)
      !           phi_grav(i,j,k) = ph0 - omega*residue/e_ijk
      !
      !           !call get_phi_star(rho,phi,i,j,k,dx,residue)
      !           !ph0=phi(i,j,k)
      !           !phi(i,j,k)=omega*residue + (1.-omega)*ph0
      !
      !           relative_error= abs(phi_grav(i,j,k)-ph0)/abs(ph0)
      !           max_error     = max(max_error,relative_error)
      !
      !           !if(relative_error > Tol)
      !           need_more=.true.
      !
      !         end do
      !       end do
      !       jsw=3-jsw
      !     end do
      !     ksw=3-ksw
      !   end do
      !
      ! end do black_red

      phi_grav = phip
      phi_grav(0   ,:   , :  ) = phi_grav(1 ,: , : )
      phi_grav(nx+1,:   , :  ) = phi_grav(nx,: , : )
      phi_grav(:   ,0   , :  ) = phi_grav(: ,1 , : )
      phi_grav(:   ,ny+1, :  ) = phi_grav(: ,ny, : )
      phi_grav(:   , :  ,0   ) = phi_grav(: , :,1  )
      phi_grav(:   , :  ,nz+1) = phi_grav(: , :,nz )

      print*,rank, 'errror:',max_error

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
