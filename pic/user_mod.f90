!=======================================================================
!> @file user_mod.f90
!> @brief User input module
!> @author C. Villarreal, M. Schneiter, A. Esquivel
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

!> @brief User imput module
!> @details  This is an attempt to have all input neede from user in a
!! single file
!!!  This module should load additional modules (i.e. star, jet, sn), to
!!  impose initial and boundary conditions (such as sources)

module user_mod

  ! load auxiliary modules
  use constants, only : msun
  use starwind
  use nbody
  implicit none

contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty
subroutine init_user_mod()

  implicit none
  !  initialize modules loaded by user
  call init_starwind()
  call init_nbody()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         pmhd, mhd, passives, rhosc, vsc, psc, cv, Tempsc, neqdyn, tsc
  use globals,    only: coords, dx ,dy ,dz, time, dt_CFL

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

  integer :: i,j,k
  real :: velx, vely, velz, dens, y0, Temp,densa,rc,Tempa,x,y,z,radp
  !--------------------------------------------------------------------------------------------
  !             MEDIO AMBIENTE
  !             con dos zonas: (dens,Temp) para r<rc y (densa,Tempa) para r>rc
  !--------------------------------------------------------------------------------------------
  dens = 1.*1.3*1.66e-24 /rhosc
  velx = 0.                /vsc
  vely = 0.                /vsc
  velz = 0.                /vsc
  Temp = 10000.              /Tempsc
  y0   = .9999    !  neutral fraction

  do i=nxmin,nxmax
     do j=nymin,nymax
        do k=nzmin,nzmax
           x=(float(i+coords(0)*nx-nxtot/2) - 0.5)*dx
           y=(float(j+coords(1)*ny-nytot/2) - 0.5)*dy
           z=(float(k+coords(2)*nz-nztot/2) - 0.5)*dz
           radp=sqrt(x**2+y**2+z**2)
        !   total density and momenta
           u(1,i,j,k) = dens
           u(2,i,j,k) = dens*velx
           u(3,i,j,k) = dens*vely
           u(4,i,j,k) = dens*velz

        !  density of neutrals
           u(neqdyn+1,i,j,k)= y0*dens
             !   passive scalar (h-hot, c-cold, i-ionized, n-neutral)
           u(neqdyn+2,i,j,k)= dens   ! passive scalar

        !  total energy (kinetic + thermal)
           u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) + cv*Temp*(2.*dens-y0*dens)

        end do
    end do
  end do

  call impose_starwind(u,0.)

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

  !  Take advantage of the use attribute to calculate the winds only
  !  once per timestep

  if (order == 1) then
    call do_nbody(time*tsc, dt_CFL*tsc)
    w(:,:,:,:) = 0.
    call impose_starwind(w,time)

    do k=nzmin,nzmax
      do j=nymin,nymax
        do i=nxmin,nxmax
          if(w(1,i,j,k) /= 0.) u(:,i,j,k) = w(:,i,j,k)
        end do
      end do
    end do

  else
    !call impose_starwind(u,time)
    do k=nzmin,nzmax
      do j=nymin,nymax
        do i=nxmin,nxmax
          if (w(1,i,j,k) /= 0.) u(:,i,j,k) = w(:,i,j,k)
        end do
      end do
    end do

  end if

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

  ! in this example a constant gravity is added
  use constants,  only : Ggrav
  use parameters, only : nx, ny, nz, nxtot, nytot, nztot, rsc, vsc2
  use globals,    only : dx, dy, dz, coords
  use starwind
  use nbody,      only : MassBH, N_star, x_star, y_star, z_star, massp
  implicit none
  integer, intent(in) :: i,j,k
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)
  real :: x,y,z, GM, rad2,xc ,yc, zc
  integer :: i_star

  !   get cell position respect to the center (BH)
  xc=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
  yc=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
  zc=(float(k+coords(2)*nz-nztot/2)+0.5)*dz

  !  Add the gravity due to the BH
  GM = Ggrav*MassBH/rsc/vsc2
  rad2 = xc**2 +yc**2 + zc**2
  ! momenta
  s(2)= s(2)-pp(1)*GM*xc/(rad2**1.5)
  s(3)= s(3)-pp(1)*GM*yc/(rad2**1.5)
  s(4)= s(4)-pp(1)*GM*zc/(rad2**1.5)
  ! energy
  s(5)= s(5)-pp(1)*GM*( pp(2)*xc +pp(3)*yc +pp(4)*zc )/(rad2**1.5 )

  do i_star=2, N_star

    x  = xc - x_star(i_star) / rsc
    y  = yc - y_star(i_star) / rsc
    z  = zc - z_star(i_star) / rsc

    GM = Ggrav*massp(i_star)/rsc/vsc2
    rad2 = x**2 +y**2 + z**2
    ! momenta
    s(2)= s(2)-pp(1)*GM*x/(rad2**1.5)
    s(3)= s(3)-pp(1)*GM*y/(rad2**1.5)
    s(4)= s(4)-pp(1)*GM*z/(rad2**1.5)
    ! energy
    s(5)= s(5)-pp(1)*GM*( pp(2)*x +pp(3)*y +pp(4)*z )/(rad2**1.5 )

  enddo

end subroutine get_user_source_terms


!=======================================================================

end module user_mod

!=======================================================================
