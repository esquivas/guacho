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

  use parameters, only: neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         neqdyn,rsc,vsc,rhosc,tempsc
  use globals,    only: coords, dx ,dy ,dz
  use constants,  only: Rjup,pi

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

  integer :: i,j,k
  real :: x,y,z, rads, velx, vely, velz, dens, sum_n
  real :: xx,tp,vp,fion,rhop, xpl,ypl,zpl,radp
  real :: Vr,Vr0, dVr,Rc,LHS, RHS, LHSold
  real :: omega,vxorb,vyorb,vzorb

  omega=2.*pi/TORB

  !Orbital planetary velocity (moves in the xz-plane)
  vxorb= -omega*Rorb*sin(phi)
  vzorb=  omega*Rorb*cos(phi)
  vyorb= 0.

  !  the star wind does not cover the entire domain, we fill here
  !  as if the exoplanet is absent
  do i=nxmin,nxmax
    do j=nymin,nymax
      do k=nzmin,nzmax
 
        ! Position measured from the centre of the grid (star)
        x=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx
        y=(real(j+coords(1)*ny-nytot/2)+0.5)*dy
        z=(real(k+coords(2)*nz)+0.5)*dz

        ! Distance from the centre of the star
        rads=sqrt(x**2+y**2+z**2)
   
        xpl=x-xp
        ypl=y
        zpl=z-zp
        radp=sqrt(xpl**2+ypl**2+zpl**2)

        !impose the planetary wind initially in a 10Rp radius region
        if(radp<3*rpw)then
            !print*,'entre'

             xx=radp*rsc/(0.361*Rjup)
             !velocity function above 5Rp
             ! y=10^b*x^a
             vp=(10**5.615*xx**0.682)/vsc
             tp=(10**3.727*xx**(-0.176))/tempsc
             rhop=10**(-15.560)*xx**(-2.682)
             fion=10**(-0.509)*xx**0.327

             VelX=vxorb+vp*xpl/radp
             VelY=vyorb+vp*ypl/radp
             VelZ=vzorb+vp*zpl/radp
             dens=rhop/rhosc

             !   total density and momenta
             u(1,i,j,k) = dens
             u(2,i,j,k) = dens*velx
             u(3,i,j,k) = dens*vely
             u(4,i,j,k) = dens*velz

             !   Here the number density of the wind and planet
             !   components separately
             !!u(neqdyn+2,i,j,k) = 0.*dens   ! xhi*rho S ion
             !!u(neqdyn+3,i,j,k) = 0.*dens   ! xhn*rho S neutro
             !!u(neqdyn+4,i,j,k) = fion*dens   ! xci*rho P ion
             !!u(neqdyn+5,i,j,k) = (1.-fion)*dens   ! xcn*rho P neutro

             !if (u(neqdyn+5,i,j,k)<=0) then
             !       print*,'planet neutralse=', u(neqdyn+5,i,j,k)
             !       stop
             !endif
             !!! ne
             !!u(neqdyn+6,i,j,k) = u(neqdyn+2,i,j,k)+u(neqdyn+4,i,j,k)
             !!!density of neutrals
             !!u(neqdyn+1,i,j,k) = u(neqdyn+3,i,j,k)+u(neqdyn+5,i,j,k)

             u(neqdyn+1,i,j,k) = (1.-fion)*dens
             !   passive scalar (tag) for stellar material
             u(neqdyn+2,i,j,k)= -1000*dens

             ! total energy
             u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) &
                          + cv*(2.*dens-u(neqdyn+1,i,j,k))*tp

       else

            if (rads<= 1.1*Rstar) then
              if(rads == 0.) rads=dx*0.1
              vr= v0
            else
              ! fill with parker solution all the mesh
              ! get sonic point radius
                  Rc = Ggrav * MassS / (2.0 * a**2.0)
    
                  Vr0 = a
                  dVr = a/10
                  if (rads*rsc < Rc) then
                      dVr = - dVr
                  end if
    
                  LHS = Vr0 * exp(-Vr0**2.0 / (2.0 * a**2.0) )
                  RHS = a * (Rc/(rads*rsc))**2.0 * exp(-2.0 * Rc/(rads*rsc) + 3.0/2.0)
    
                  Vr = Vr0
    
                  do while (abs(LHS/RHS-1.0) > 1.e-8)
                          ! save old LHS
                          LHSold = LHS
                          !update Vr
                          Vr = Vr + dVr
                          ! calculate new LHS
                          LHS = Vr * exp(-Vr**2.0/(2.0 *a**2.0) )
                          ! see if crossed solution
                          if((LHS/RHS-1.0)*(LHSold/RHS-1.0)<0)then
                              dVr=-dVr/2.0
                          end if
                  end do
                  vr=vr/vsc
            endif
    
            VelX=vr*X/RADS
            VelY=vr*Y/RADS
            VelZ=vr*Z/RADS
            DENS=DSW*RSW**2/RADS**2
            !   total density and momenta
            u(1,i,j,k) = dens
            u(2,i,j,k) = dens*velx
            u(3,i,j,k) = dens*vely
            u(4,i,j,k) = dens*velz
    
            !   Here the number density of the wind and planet
            !   components separately
            !!u(neqdyn+2,i,j,k) = 0.999*dens   ! xhi*rho S ion
            !!u(neqdyn+3,i,j,k) = 1.e-4*dens   ! xhn*rho S neutro
            !!u(neqdyn+4,i,j,k) = 0.*dens   ! xci*rho P ion
            !!u(neqdyn+5,i,j,k) = 0.*dens   ! xcn*rho P neutro
            !!if (u(neqdyn+3,i,j,k)<=0) then
            !!       print*,'star neutrals=', u(neqdyn+3,i,j,k)
            !!       stop
            !!endif
            !!! ne
            !!u(neqdyn+6,i,j,k) = u(neqdyn+2,i,j,k)+u(neqdyn+4,i,j,k)
            !density of neutrals
            !!u(neqdyn+1,i,j,k) = u(neqdyn+3,i,j,k)+u(neqdyn+5,i,j,k)
            u(neqdyn+1,i,j,k) = fneutro_s*dens
            !   passive scalar (tag) for stellar material
            u(neqdyn+2,i,j,k)= 1000*dens
    
            ! total energy
            u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) &
                          + cv*(2.*dens-u(neqdyn+1,i,j,k))*Tsw
            
        endif

      end do
    end do
  end do

  call impose_exo(u,0.)

end subroutine initial_conditions

!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)
!> @param integer [in] order : order (mum of cells to be filled in case
!> domain boundaries are being set)

subroutine impose_user_bc(u,order)

  use parameters, only:  neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals,    only: time
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: order

  !  In this case the boundary is the same for 1st and second order)
  if (order >= 1) then
    call impose_exo(u,time)
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

  ! Adds the Rad Pressure according to the Beta profile of Bourrier
  use constants,  only : Ggrav
  use parameters, only : nx, ny, nz, nxtot, nytot, nztot, rsc, vsc2,&
                         beta_pressure, vsc,neqdyn
  use globals,    only : dx, dy, dz, coords
  use exoplanet
  use radpress
 
  implicit none
  integer, intent(in) :: i, j, k
  integer             :: l, index
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)
  integer, parameter  :: nb=2   ! 2 particles
  real :: x(nb),y(nb),z(nb), GM(nb), rad2(nb)
  real :: xc ,yc, zc
  real :: v, fracv, frac_neutro !, a, b, c

  GM(1)= Ggrav*MassS/rsc/vsc2
  GM(2)= Ggrav*MassP/rsc/vsc2

  !   get cell position
  xc=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx
  yc=(real(j+coords(1)*ny-nytot/2)+0.5)*dy
  zc=(real(k+coords(2)*nz)+0.5)*dz

  ! calculate distance from the sources
  ! star
  x(1)=xc
  y(1)=yc
  z(1)=zc
  rad2(1) = x(1)**2 +y(1)**2 + z(1)**2
  ! planet
  x(2)=xc-xp
  y(2)=yc
  z(2)=zc-zp
  rad2(2) = x(2)**2 +y(2)**2 + z(2)**2


if ( beta_pressure ) then

    beta(i,j,k) = 0.
    !  do only outside BC
    if( (rad2(1) >= rsw**2) .and. (rad2(2) >= rpw**2) ) then

      frac_neutro = pp(neqdyn+1)/pp(1)        !!Each cell feels a given pressure proporcional to the neutrals fraction
      !  Radial velocity in km s^-1
      v =  ( (pp(2)*xc + pp(3)*yc + pp(4)*zc)/sqrt(rad2(1)) ) * (vsc)

      fracv = (v-vr(1))/(vr(Nr)-vr(1))*Nr
      index = int(fracv)+1

      if (index < 1) then
        index = 1
      else if ( index > Nr-1 ) then
        index = Nr-1
      end if

      Beta(i,j,k) = ( Br(index)  + (v-vr(index))*( Br(index+1)-Br(index) ) / ( vr(index+1)-vr(index) ) ) *frac_neutro!*active

      !!Linear interpolation for Beta, active allows turn on the Beta term.
      GM(1)=GM(1)*(1.-Beta(i,j,k)) !!Update scale factor GM

    end if

  endif

    ! update source terms with gravity
    do l=1, nb
      ! momenta
      s(2)= s(2)-pp(1)*GM(l)*x(l)/(rad2(l)**1.5)
      s(3)= s(3)-pp(1)*GM(l)*y(l)/(rad2(l)**1.5)
      s(4)= s(4)-pp(1)*GM(l)*z(l)/(rad2(l)**1.5)
      ! energy
      s(5)= s(5)-pp(1)*GM(l)*( pp(2)*x(l) +pp(3)*y(l) +pp(4)*z(l) )  &
      /(rad2(l)**1.5 )
    end do

end subroutine get_user_source_terms


!=======================================================================

end module user_mod

!=======================================================================
