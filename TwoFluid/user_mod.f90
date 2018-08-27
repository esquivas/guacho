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
  use exoplanet

  implicit none

contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty
subroutine init_user_mod()

  implicit none
  !  initialize modules loaded by user
  call init_exo()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         pmhd, mhd, passives, mu, twofluid, Tempsc, rhosc, vsc, &
                         nxtot,nytot,nztot
  use globals,    only: coords, dx ,dy ,dz, rank

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
!  real, intent(out) :: un(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

  integer :: i,j,k
  real :: x,y,z, rads, velx, vely, velz, dens,cpi,xpl,ypl,zpl,radp

  !  the star wind does not cover the entire domain, we fill here
  !  as if the exoplanet is absent
  ! if twofluids, u (prim) is used for ions an un (primn) for neutrals

  do k=nzmin,nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax

        ! Position measured from the centre of the grid (planet)
        xpl=(float(i+coords(0)*nx-nxtot/2) - 0.5)*dx
        ypl=(float(j+coords(1)*ny-nytot/2) - 0.5)*dy
        zpl=(float(k+coords(2)*nz-nztot/2) - 0.5)*dz

        ! Distance from the centre of the planet (centred)
        radp=sqrt(xpl**2+ypl**2+zpl**2)

        if(radp > Rpw) then

          u(1,i,j,k) = dsw
          u(2,i,j,k) = 0.!dsw*vsw
          u(3,i,j,k) = 0.
          u(4,i,j,k) = 0.
          u(6,i,j,k) = 0.
          u(7,i,j,k) = Bsw
          u(8,i,j,k) = 0.
          u(5,i,j,k) = cv*(dsw/0.63)*Tsw  + 0.5*Bsw**2 ! &
                       !  + 0.5*dsw*vsw**2

        else
          u(1,i,j,k) = 1E-5*dsw
          u(2,i,j,k) = 0.
          u(3,i,j,k) = 0.
          u(4,i,j,k) = 0.
          u(6,i,j,k) = 0.
          u(7,i,j,k) = 0.
          u(8,i,j,k) = 0.
          u(5,i,j,k) = cv*(1e-5*dsw/0.63)*Tsw

        end if
      end do
    enddo
  enddo
end subroutine initial_conditions


!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditionsN(u)

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

  !  outside planet wind
  u(1,:,:,:) = 0.1*dpw
  u(2,:,:,:) = 0.
  u(3,:,:,:) = 0.
  u(4,:,:,:) = 0.
  u(6,:,:,:) = 0.
  u(7,:,:,:) = 0.
  u(8,:,:,:) = 0.
  u(5,:,:,:) = cv*(0.1*dpw/1.3)*Tpw

  ! rewrite inside planet
  call impose_exo(u, 0.)

end subroutine initial_conditionsN

!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)
!> @param integer [in] order : order (mum of cells to be filled in case
!> domain boundaries are being set)

subroutine impose_user_bc(u,order,neutral)

  use parameters, only:  neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         rhosc, vsc, Tempsc
  use globals,    only: time, coords
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: order
  logical, optional, intent(in) :: neutral
  !  real :: temp, de
  integer :: i, j, k, im, ip, jm, jp

  !  In this case the boundary is the same for 1st and second order)
  select case(order)

  case(1)

    if(present(neutral).and.neutral) then
      call impose_exo(u,time)
    else

      !  primer orden
      !  viento plano paralelo
      do k=0,nz+1
        do j = 0,ny+1
          if (coords(0) == 0) then
            u(1,0,j,k) = dsw
            u(2,0,j,k) = dsw*vsw
            u(3,0,j,k) = 0.
            u(4,0,j,k) = 0.
            u(6,0,j,k) = 0.
            u(7,0,j,k) = Bsw
            u(8,0,j,k) = 0.
            u(5,0,j,k) = cv*(dsw/0.63)*Tsw  + &
                 0.5*dsw*vsw**2 + 0.5*bsw**2
          end if

          do i=0,nx+1

             if(flagP(i,j,k)) then

              u(1,i,j,k) = 1E-5*dsw
              u(2:8,i,j,k) = 0.

              !   reflect x
              if ( flagP(i-1,j,k).and.(.not.flagP(i+1,j,k)).and.(i>0)) then
                u(1,i,j,k) =  u(1,i+1,j,k)
                u(2,i,j,k) = -u(2,i+1,j,k)
              endif

              if (flagP(i+1,j,k).and.(.not.flagP(i-1,j,k)).and.(i<nx+1)) then
                u(1,i,j,k) =  u(1,i-1,j,k)
                u(2,i,j,k) = -u(2,i-1,j,k)
              endif

             !   reflect y
             if (flagP(i,j-1,k).and.(.not.flagP(i,j+1,k)).and.(j>0)) then
               u(1,i,j,k) = u(1,i,j+1,k)
               u(3,i,j,k) =-u(3,i,j+1,k)
             endif
             if (flagP(i,j+1,k).and.(.not.flagP(i,j-1,k)).and.(j<ny+1)) then
               u(1,i,j,k) = u(1,i,j-1,k)
               u(3,i,j,k) =-u(3,i,j-1,k)
             endif

             !   reflect z
             if (flagP(i,j,k-1).and.(.not.flagP(i,j,k+1)).and.(k>0)) then
               u(1,i,j,k) = u(1,i,j,k+1)
               u(4,i,j,k) =-u(4,i,j,k+1)
             endif
             if (flagP(i,j,k+1).and.(.not.flagP(i,j,k-1)).and.(k<nz+1)) then
               u(1,i,j,k) = u(1,i,j,k-1)
               u(4,i,j,k) =-u(4,i,j,k-1)
             endif

              u(5,i,j,k) = 0.5*(u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)/u(1,i,j,k) + &
                           cv*(u(1,i,j,k)/0.63)*Tsw     +  &
                            + 0.5*bsw**2

            end if

          end do
        end do
      end do

    end if

  case(2)
    ! segundo orden
    if(present(neutral).and.neutral) then
      call impose_exo(u,time)
    else

      ! viento plano paralelo
      do k=-1,nz+2
        do j = -1,ny+2
          if (coords(0) == 0) then
            u(1,-1:0,j,k) = dsw
            u(2,-1:0,j,k) = dsw*vsw
            u(3,-1:0,j,k) = 0.
            u(4,-1:0,j,k) = 0.
            u(6,-1:0,j,k) = 0.
            u(7,-1:0,j,k) = Bsw
            u(8,-1:0,j,k) = 0.
            u(5,-1:0,j,k) = cv*(dsw/0.63)*Tsw  + &
                 0.5*dsw*vsw**2 + 0.5*bsw**2
          end if

          do i=-1,nx+2

            if(flagP(i,j,k)) then

              u(1,i,j,k) = 1E-5*dsw
              u(2:8,i,j,k) = 0.

              !   reflect x
              if (flagP(i-1,j,k).and.(.not.flagP(i+1,j,k)).and.(i>-1)) then
                u(1,i  ,j,k) = u(1,i+1,j,k)
                u(2,i  ,j,k) =-u(2,i+1,j,k)
                !u(1,i-1,j,k) = u(1,i+2,j,k)
                !u(2,i-1,j,k) =-u(2,i+2,j,k)
              endif
              if (flagP(i+1,j,k).and.(.not.flagP(i-1,j,k)).and.(i<nx+2)) then
                u(1,i  ,j,k) = u(1,i-1,j,k)
                u(2,i  ,j,k) =-u(2,i-1,j,k)
                !u(1,i+1,j,k) = u(1,i-2,j,k)
                !u(2,i+1,j,k) =-u(2,i-2,j,k)
              endif

              !   reflect y
              if (flagP(i,j-1,k).and.(.not.flagP(i,j+1,k)).and.(j>-1)) then
                u(1,i,j  ,k) = u(1,i,j+1,k)
                u(3,i,j  ,k) =-u(3,i,j+1,k)
                !u(1,i,j-1,k) = u(1,i,j+2,k)
                !u(3,i,j-1,k) =-u(3,i,j+2,k)
              endif
              if (flagP(i,j+1,k).and.(.not.flagP(i,j-1,k)).and.(j<ny+2)) then
                u(1,i,j  ,k) = u(1,i,j-1,k)
                u(3,i,j  ,k) =-u(3,i,j-1,k)
                !u(1,i,j+1,k) = u(1,i,j-2,k)
                !u(3,i,j+1,k) =-u(3,i,j-2,k)
              endif

              !   reflect z
              if (flagP(i,j,k-1).and.(.not.flagP(i,j,k+1)).and.(k>-1)) then
                u(1,i,j,k  ) = u(1,i,j,k+1)
                u(4,i,j,k  ) =-u(4,i,j,k+1)
                !u(1,i,j,k-1) = u(1,i,j,k+2)
                !u(4,i,j,k-1) =-u(4,i,j,k+2)
              endif
              if (flagP(i,j,k+1).and.(.not.flagP(i,j,k-1)).and.(k<nz+2)) then
                u(1,i,j,k  ) = u(1,i,j,k-1)
                u(4,i,j,k  ) =-u(4,i,j,k-1)
                !u(1,i,j,k+1) = u(1,i,j,k-2)
                !u(4,i,j,k+1) =-u(4,i,j,k-2)
              endif
              u(5,i,j,k) = 0.5*(u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)/u(1,i,j,k) + &                             cv*(u(1,i,j,k)/0.63)*Tsw     +  &
                             + 0.5*bsw**2 + cv*(u(1,i,j,k)/0.63)*Tsw
            end if

          end do
        end do
      end do

    end if

  end select

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
  use exoplanet
  implicit none
  integer, intent(in) :: i,j,k
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)
  integer, parameter  :: nb=2   ! 2 particles
  real :: x(nb),y(nb),z(nb), GM(nb), rad2(nb)
  integer :: l
  real    :: xc ,yc, zc

  ! Gravity felt by the planet at the centrer of the grid
  GM(1)=Ggrav*MassP/rsc/vsc2
    ! Gravity felt by the planet at the centrer of the grid
!  GM(2)=0.3*Ggrav*MassS/rsc/vsc2


  !   get cell position
  xc=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
  yc=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
  zc=(float(k+coords(2)*nz-nztot/2)+0.5)*dz

  ! calculate distance from the sources
  ! planet
  x(1)=xc
  y(1)=yc
  z(1)=zc
  rad2(1) = x(1)**2 +y(1)**2 + z(1)**2

  ! update source terms with gravity

!     ! momenta
     s(2)= s(2)-pp(1)*GM(1)*x(1)/(rad2(1)**1.5)
     s(3)= s(3)-pp(1)*GM(1)*y(1)/(rad2(1)**1.5)
     s(4)= s(4)-pp(1)*GM(1)*z(1)/(rad2(1)**1.5)
     ! energy
     s(5)= s(5)-pp(1)*GM(1)*( pp(2)*x(1) +pp(3)*y(1) +pp(4)*z(1) )  &
          /(rad2(1)**1.5 )

end subroutine get_user_source_terms

!=======================================================================

end module user_mod

!=======================================================================
