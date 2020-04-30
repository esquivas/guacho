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
  use snr

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
          gamma, nx, ny, nz, nxtot, nytot, nztot, N_MP, NBinsSEDMP,            &
          np, xmax, ymax, zmax, vsc2, vsc, bsc, np

    use globals,    only : coords, dx ,dy ,dz, rank,                           &
                           Q_MP0, partID, partOwner, n_activeMP, MP_SED, comm3d
    use constants,  only : pi, eV
    use lmp_module, only : interpBD
    use utilities,  only : isInDomain
    implicit none
    real, intent(out) :: u(neq,nxmin,nymin:nymax,nzmin:nzmax)
    !logical ::  isInDomain
    integer :: i,j,k,ip, err
    real, allocatable   :: data(:,:,:,:)
    real, parameter     :: rho_env=0.005, T_env=1000.0/Tempsc, B_env=2.e-6/bsc
    ! real    :: dens, temp, rad, x, y, z, radSN, pressSN, eSN, nu
    integer :: yj,xi
    real    :: pos(3)
    !  initial MP spectra parameters
    real    :: emin, emax, deltaE, gamma_lmp, N0, chi0
    logical, parameter :: uniform = .false. ! place the MPs uniformly
    real    :: xc, yc, zc
    integer :: ind(3), l
    real    :: weights(8), rhoI

    !ENVIRONMENT  (Read data anf fill only physical domain)
    allocate (data(7,nxtot,nytot,nztot))
    do ip=0, np-1
      if(rank == ip) then

        ! Read entire input
        open(unit=10,file='/storage2/esquivel/turb-mhd-sims/b1p.01/dens.bin',  &
             access='stream')
        open(unit=11,file='/storage2/esquivel/turb-mhd-sims/b1p.01/velx.bin',  &
             access='stream')
        open(unit=12,file='/storage2/esquivel/turb-mhd-sims/b1p.01/vely.bin',  &
             access='stream')
        open(unit=13,file='/storage2/esquivel/turb-mhd-sims/b1p.01/velz.bin',  &
             access='stream')
        open(unit=14,file='/storage2/esquivel/turb-mhd-sims/b1p.01/magx.bin',  &
             access='stream')
        open(unit=15,file='/storage2/esquivel/turb-mhd-sims/b1p.01/magy.bin',  &
             access='stream')
        open(unit=16,file='/storage2/esquivel/turb-mhd-sims/b1p.01/magz.bin',  &
             access='stream')

        do k=1,nztot
          do j=1,nytot
            do i=1,nxtot
              read(10) data(1,i,j,k)
              read(11) data(2,i,j,k)
              read(12) data(3,i,j,k)
              read(13) data(4,i,j,k)
              read(14) data(5,i,j,k)
              read(15) data(6,i,j,k)
              read(16) data(7,i,j,k)

            end do
          end do
        end do

        close(10)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)

        !   Fill physical domain
        do k=1,nz
          do j=1,ny
            do i=1,nx
              !  density
              u(1,i,j,k) = data(1,i+coords(0)*nx,j+coords(1)*ny,k+coords(2)*nz) * rho_env
              !  momenta
              u(2,i,j,k) = data(2,i+coords(0)*nx,j+coords(1)*ny,k+coords(2)*nz) * u(1,i,j,k)
              u(3,i,j,k) = data(3,i+coords(0)*nx,j+coords(1)*ny,k+coords(2)*nz) * u(1,i,j,k)
              u(4,i,j,k) = data(4,i+coords(0)*nx,j+coords(1)*ny,k+coords(2)*nz) * u(1,i,j,k)
              !  B field
              u(6,i,j,k) = data(5,i+coords(0)*nx,j+coords(1)*ny,k+coords(2)*nz) * Benv
              u(7,i,j,k) = data(6,i+coords(0)*nx,j+coords(1)*ny,k+coords(2)*nz) * Benv
              u(8,i,j,k) = data(7,i+coords(0)*nx,j+coords(1)*ny,k+coords(2)*nz) * Benv
              !  Total energy
              u(5,i,j,k) = 0.5*( u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2 )/ u(1,i,j,k) &                                    &
                         + 0.5*( u(6,i,j,k)**2+u(7,i,j,k)**2+u(8,i,j,k)**2 )  &
                         + cv * u(1,i,j,k) * T_env

            end do
          end do
        end do

      end if

      call mpi_barrier(comm3d, err)

    end do


    xc = 12.* pc/rsc
    yc = 12.* pc/rsc
    zc = 12.* pc/rsc

    call impose_snr(u,xc,yc,zc)


    !  TRACER PARTICLES
    !  initialize Owners (-1 means no body has claimed the particle)
    partOwner(:) = -1
    !  initialize Particles ID, not active is ID 0
    partID(:)    =  0
    n_activeMP   =  0

    ! delta E in log bins (this is in ergs)
    Emin = 0.63*1e6 *eV !/(rhosc*rsc**3*vsc2)
    Emax = 0.31*1e12*eV !/(rhosc*rsc**3*vsc2)
    gamma_lmp = 3.
    deltaE = (log10(Emax)-log10(Emin))/ (real(NBinsSEDMP)-1.)

    if (uniform) then
      !Insert homogenously distributed particles
      do yj=4,ny,8
        do xi=4,nx,8

          !  position of MPs (respect to a corner --needed by isInDomain--)
          pos(1)= real(xi+ coords(0)*nx + 0.5) * dx
          pos(2)= real(yj+ coords(1)*ny + 0.5) * dy
          pos(3)= real( 1+ coords(2)*nz + 0.5) * dz

          if(isInDomain(pos) ) then
            n_activeMP            = n_activeMP + 1
            partOwner(n_activeMP) = rank
            partID   (n_activeMP) = n_activeMP + rank*N_MP
            Q_MP0(n_activeMP,:) = 0.
            Q_MP0(n_activeMP,1:3) = pos(:)

            !  Interpolate density to normalize SED
            call interpBD(Q_MP0(n_activeMP,1:3),ind,weights)
            rhoI = 0.0
            l    = 1
            do k= ind(3),ind(3)+1
              do j=ind(2),ind(2)+1
                do i=ind(1),ind(1)+1
                  rhoI = rhoI + u(1,i,j,k) * weights(l)
                  l  = l + 1
                end do
              end do
            end do
            N0     =  1.0e-6 / rhoI
            chi0   = N0 * (1.-gamma_lmp)/                                      &
                          ( Emax**(1.-gamma_lmp)-Emin**(1.-gamma_lmp) )

            do i = 1,NBinsSEDMP
              MP_SED(1,i,n_activeMP)=10.**(log10(Emin)+real(i-1)*deltaE)
              MP_SED(2,i,n_activeMP)= chi0*MP_SED(1,i,n_activeMP)**(-gamma_lmp)
            end do

          endif

        end do
      end do
    else
      !  this will insert N_MP/2 (N_MP is set in parameters) randomly
      !  distributed within each processor domain
      do while (n_activeMP < N_MP/2)
        call random_number(pos(1:3))
        pos(1) = pos(1) * xmax
        pos(2) = pos(2) * ymax
        pos(3) = pos(3) * zmax

        if(isInDomain(pos) ) then
            n_activeMP            = n_activeMP + 1
            partOwner(n_activeMP) = rank
            partID   (n_activeMP) = n_activeMP + rank*N_MP
            Q_MP0(n_activeMP,:) = 0.
            Q_MP0(n_activeMP,1:3) = pos(:)

            !  Interpolate density to normalize SED
            call interpBD(Q_MP0(n_activeMP,1:3),ind,weights)
            rhoI = 0.0
            l    = 1
            do k= ind(3),ind(3)+1
              do j=ind(2),ind(2)+1
                do i=ind(1),ind(1)+1
                  rhoI = rhoI + u(1,i,j,k) * weights(l)
                  l  = l + 1
                end do
              end do
            end do
            N0     =  1.0e-6 / rhoI
            chi0   = N0 * (1.-gamma_lmp)/                                      &
                          ( Emax**(1.-gamma_lmp)-Emin**(1.-gamma_lmp) )

            do i = 1,NBinsSEDMP
              MP_SED(1,i,n_activeMP)=10.**(log10(Emin)+real(i-1)*deltaE)
              MP_SED(2,i,n_activeMP)= chi0*MP_SED(1,i,n_activeMP)**(-gamma_lmp)
            end do

          endif

      end do
    end if

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
