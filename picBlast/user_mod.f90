!=======================================================================
!> @file user_mod.f90
!> @brief User input module
!> @author C. Villarreal, M. Schneiter, A. Esquivel
!> @date 4/May/2016
!
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

!> @brief User imput module
!> @details  This is an attempt to have all input neede from user in a
!> single file
!> This module should load additional modules (i.e. star, jet, sn), to
!> impose initial and boundary conditions (such as sources)

module user_mod

  ! load auxiliary modules

  implicit none

contains

  !=====================================================================
  !> @brief Initializes variables in the module, as well as other
  !! modules loaded by user.
  !! @n It has to be present, even if empty
  subroutine init_user_mod()

    implicit none
    !  initialize modules loaded by user
    !  call init_vortex()

  end subroutine init_user_mod

  !=====================================================================
  !> @brief Here the domain is initialized at t=0
  !> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
  !> conserved variables
  !> @param real [in] time : time in the simulation (code units)
  subroutine initial_conditions(u)

    use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax,      &
          pmhd, mhd, passives, rsc,rhosc, vsc, psc, cv, Tempsc, neqdyn, tsc,   &
          gamma, nx, ny, nz, nxtot, nytot, nztot, N_MP, NBinsSEDMP, np

    use globals,    only : coords, dx ,dy ,dz, rank,                           &
                           Q_MP0, partID, partOwner, n_activeMP, MP_SED
    use constants,  only : pi
    use utilities,  only : isInDomain

    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    !logical ::  isInDomain
    integer :: i,j,k
    real    :: dens, temp, rad, x, y, z, radSN, pressSN, eSN, nu
    integer :: yj,xi
    real    :: pos(3), E0
    real, parameter :: gamma_pic=3.,de=6./real(NBinsSEDMP)
    real    :: emin, emax !lower and upper energy limits of mp spectra
    integer :: m, ntot_mp !m index of initial spectra, ntot_mp number of mp
    !-----------------------------------------------------------------------------
    !       HIDRODINAMICA : MEDIO AMBIENTE
    !       BLAST PROBLEM
    !       (high Order Finite Difference and Finite Volume WENO Schemes
    !       and Discontinuous Galerkin Methodsfor CFDChi-Wang Shu)
    !-----------------------------------------------------------------------------

    !ENVIRONMENT

    u(1,:,:,:) = 1.
    u(2,:,:,:) = 0.
    u(3,:,:,:) = 0.
    u(4,:,:,:) = 0.
    u(5,:,:,:) = cv*1e-5
    u(6,:,:,:) = 0.
    u(7,:,:,:) = .25
    u(8,:,:,:) = 0.

    !We have to impose the blast according to .....FLASH CODE?
    eSN = 2. !energy of the SN
    nu  = 3. !3 is for spherical and 2 for cylindrical

    !blast
    do k=nzmin,nzmax
      do j=nymin,nymax
        do i=nxmin,nxmax

          !  this is the position with respect of the grid center
          x= ( real(i+coords(0)*nx-nxtot/2) - 0.5) *dx
          y= ( real(j+coords(1)*ny-nytot/2) - 0.5) *dy
          z= ( real(k+coords(2)*nz-nztot/2) - 0.5) *dz
          rad=sqrt(x**2+y**2)


          if (rad.le.dx*4.) then
            pressSN=3.*(gamma-1)*eSN/((nu+1)*3.14*rad)
            !  total energy (kinetic + thermal)
            u(5,i,j,k) = cv*pressSN
          endif

        end do
      end do
    end do

    !  TRACER PARTICLES
    !  initialize Owners (-1 means no body has claimed the particle)
    partOwner(:) = -1
    !  initialize Particles ID, not active is ID 0
    partID(:)    =  0
    n_activeMP   =  0

    !Insert homogenously distributed particles
    do yj=2,ny,8
      do xi=2,nx,8

        !  position of particles (respect to a corner --needed by isInDomain--)
        pos(1)= real(xi+ coords(0)*nx + 0.5) * dx
        pos(2)= real(yj+ coords(1)*ny + 0.5) * dy
        pos(3)= real( 1+ coords(2)*nz + 0.5) * dz

        if(isInDomain(pos) ) then

          n_activeMP            = n_activeMP + 1
          partOwner(n_activeMP) = rank
          partID   (n_activeMP) = n_activeMP + rank*N_MP
          Q_MP0(n_activeMP,1:3) = pos(:)
          E0 =  10**( -gamma_pic*(-4 + 1.5*de) )
          Q_MP0(n_activeMP,4:10) = 0.
          do i = 1,NBinsSEDMP
            MP_SED(1,i,n_activeMP)=10**(-4+(0.5+real(i))*de)
            MP_SED(2,i,n_activeMP)= MP_SED(1,i,n_activeMP)**(-gamma_pic)/E0
          end do

        endif

      end do
    end do

    print*, rank, 'has ', n_activeMP, ' active MPs'

  end subroutine initial_conditions

  !=====================================================================
  !> @brief User Defined Boundary conditions
  !> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
  !! conserved variables
  !> @param real [in] time : time in the simulation (code units)
  !> @param integer [in] order : order (mum of cells to be filled in case
  !> domain boundaries are being set)
  subroutine impose_user_bc(u,order)
    use parameters, only:  neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, tsc
    use globals,    only: time, dt_CFL
    implicit none
    real, intent(out)    :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, save           :: w(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    integer, intent(in)  :: order
    integer              :: i, j, k

  end subroutine impose_user_bc

  !=======================================================================
  !> @brief User Defined source terms
  !> This is a generic interrface to add a source term S in the equation
  !> of the form:  dU/dt+dF/dx+dG/dy+dH/dz=S
  !> @param real [in] pp(neq) : vector of primitive variables
  !> @param real [inout] s(neq) : vector with source terms, has to add to
  !>  whatever is there, as other modules can add their own sources
  !> @param integer [in] i : cell index in the X direction
  !> @param integer [in] j : cell index in the Y direction
  !> @param integer [in] k : cell index in the Z direction
  subroutine get_user_source_terms(pp,s, i, j , k)
    use parameters, only : neq, NBinsSEDMP
    implicit none
    real, intent(in)   :: pp(neq)
    real, intent(out)  :: s(neq)
    integer :: i, j, k

  end subroutine get_user_source_terms

  !=======================================================================

end module user_mod

!=======================================================================
