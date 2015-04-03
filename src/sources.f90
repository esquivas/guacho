!=======================================================================
!> @file sources.f90
!> @brief Adds source terms
!> @author Alejandro Esquivel
!> @date 2/Nov/2014

! Copyright (c) 2014 A. Esquivel, M. Schneiter, C. Villareal D'Angelo
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

!> @brief Adds source terms
!> @details This module adds the source terms from gravity, radiation
!! pressure (not fully tested), and div(B) cleaning if the 8 wave
!! scheme is used

#if defined(GRAV) || defined(RADPRES) || defined(EIGHT_WAVE)

  module sources
  use parameters, only : neq, neqdyn, nxtot, nytot, nztot, &
                         rsc, rhosc, vsc2, nx, ny, nz
  use globals,    only : dx, dy, dz, coords
  implicit none
  
contains

!=======================================================================

!> @brief Gets position in the grid
!> @details Gets the position and spherical radius calculated with
!! respect to the center of the grid
!> @param integer [in] i : index in the X direction
!> @param integer [in] j : index in the Y direction
!> @param integer [in] k : index in the Z direction
!> @param real [out] x : X position form the center of the grid (code units)
!> @param real [out] y : Y position form the center of the grid (code units)
!> @param real [out] z : Z position form the center of the grid (code units)
!> @param real [out] r : Spherical radius form the center of the grid 
!! (code units)

subroutine getpos(i,j,k,x,y,z,r)

  implicit none
  integer, intent(in)  :: i, j, k
  real,    intent(out) :: x, y, z, r

  x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
  y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
  z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
  
  r  = sqrt(x**2 +y**2 +z**2 )
  
end subroutine getpos

!=======================================================================

#ifdef GRAV

!> @brief Gravity due to point sources
!> @details Adds the gravitational force due to point particles, at this
!!  moment is fixed to two point sources (exoplanet)
!> @param real [in] xc : X position of the cell
!> @param real [in] yc : Y position of the cell
!> @param real [in] zc : Z position of the cell
!> @param real [in] pp(neq) : vector of primitive variables
!> @param real [out] s(neq) : vector with source terms

subroutine grav_source(xc,yc,zc,pp,s)
  use constants, only : Ggrav
  use exoplanet  ! this module contains the position of the planet
  implicit none
  real, intent(in)    :: xc, yc, zc
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)
  integer, parameter  :: nb=2   ! 2 particles
  real :: x(nb),y(nb),z(nb), GM(nb), rad2(nb)
  integer :: i

  GM(1)=0.3*Ggrav*MassS/rsc/vsc2
  GM(2)=Ggrav*MassP/rsc/vsc2

  !calculate distance from the sources
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

  ! update source terms
  do i=1, nb
    ! momenta
    s(2)= s(2)-pp(1)*GM(i)*x(i)/(rad2(i)**1.5)
    s(3)= s(3)-pp(1)*GM(i)*y(i)/(rad2(i)**1.5)
    s(4)= s(4)-pp(1)*GM(i)*z(i)/(rad2(i)**1.5)
    ! energy
    s(5)= s(5)-pp(1)*GM(i)*( pp(2)*x(i) +pp(3)*y(i) +pp(4)*z(i) )  &
           /(rad2(i)**1.5 )
  end do


end subroutine grav_source

#endif

!=======================================================================

!> @brief Radiation pressure force
!> @details Adds the radiaiton pressure force due to photo-ionization
!> @param integer [in] i : cell index in the X direction
!> @param integer [in] j : cell index in the Y direction
!> @param integer [in] k : cell index in the Z direction
!> @param real [in] xc : X position of the cell
!> @param real [in] yc : Y position of the cell
!> @param real [in] zc : Z position of the cell
!> @param reak [in] rc : @f$ \sqrt{x^2+y^2+z^2} @f$
!> @param real [in] pp(neq) : vector of primitive variables
!> @param real [out] s(neq) : vector with source terms

#ifdef RADPRES

  subroutine radpress_source(i,j,k,xc,yc,zc,rc,pp,s)
#ifdef RADDIFF
  use difrad
#endif
    implicit none
    integer, intent(in)  :: i,j,k
    real,    intent(in)  :: xc, yc, zc, rc, pp(neq)
    real,    intent(inout) :: s(neq)
    real :: radphi
    !  the following is h/912Angstroms = h/lambda
    real :: hlambda= 7.265e-22, Frad

  radphi= ph(i,j,k) 

  !  update the source terms
  Frad=radphi*pp(6)*hlambda*rsc/rhosc/vsc2
  !  momenta
  s(2) = s(2) + Frad*xc/rc
  s(3) = s(3) + Frad*yc/rc
  s(4) = s(4) + Frad*zc/rc
  !  energy
  s(5) = s(5)+  Frad*( xc*pp(2) + yc*pp(3) + zc*pp(4) )/rc

  end subroutine radpress_source

#endif

!=======================================================================

#ifdef EIGHT_WAVE

!> @brief Computes div(B)
!> @details Computes div(B)
!> @param integer [in] i : cell index in the X direction
!> @param integer [in] j : cell index in the Y direction
!> @param integer [in] k : cell index in the Z direction
!> @param real [out] d :: div(B)

subroutine divergence_B(i,j,k,d)
use globals
implicit none
integer, intent(in) :: i,j,k
real, intent(out)   :: d
  
d=  (primit(6,i+1,j,k)-primit(6,i-1,j,k))/(2.*dx)  &
  + (primit(7,i,j+1,k)-primit(7,i,j-1,k))/(2.*dy)  &
  + (primit(8,i,j,k+1)-primit(8,i,j,k-1))/(2.*dz)  

end subroutine divergence_B

!=======================================================================

!> @brief 8 Wave source terms for div(B) correction
!> @details  Adds terms proportional to div B in Faraday's Law,
!! momentum equationand energy equation as propoes in Powell et al. 1999
!> @param integer [in] i : cell index in the X direction
!> @param integer [in] j : cell index in the Y direction
!> @param integer [in] k : cell index in the Z direction
!> @param real [in] pp(neq) : vector of primitive variables
!> @param real [out] s(neq) : vector with source terms

subroutine divbcorr_source(i,j,k,pp,s)
 
  implicit none
  integer, intent(in) :: i, j, k
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)
  real                :: divB
  
  call divergence_B(i,j,k,divB)

  ! update source terms
    ! momenta
    s(2)= s(2)-divB*pp(6) 
    s(3)= s(3)-divB*pp(7) 
    s(4)= s(4)-divB*pp(8) 

    ! energy
    s(5)= s(5)-divB*(pp(2)*pp(6)+pp(3)+pp(7)+pp(4)*pp(8))      

    ! Faraday law
    s(6)=s(6)-divB*pp(2)
    s(7)=s(7)-divB*pp(3)
    s(8)=s(8)-divB*pp(4)

end subroutine divbcorr_source

#endif

!=======================================================================

!> @brief Upper level wrapper for sources
!> @details Upper level wrapper for sources
!! @n Main driver, this is called from the upwind stepping
!> @param integer [in] i : cell index in the X direction
!> @param integer [in] j : cell index in the Y direction
!> @param integer [in] k : cell index in the Z direction
!> @param real [in] prim(neq) : vector of primitive variables
!> @param real [out] s(neq) : vector with source terms'

subroutine source(i,j,k,prim,s)

  implicit none
  integer, intent(in)  :: i, j, k
  real, intent(in)     :: prim(neq)
  real, intent(out)    :: s(neq)
  real :: x, y, z, r
  
  ! resets the source terms
  s(:) = 0.
  
  ! position with respect to the center of the grid
  call getpos( i, j, k, x, y ,z, r) 

#ifdef GRAV
  !  point source(s) gravity
  call grav_source(x,y,z,prim,s)
#endif

#ifdef RADPRES
  !  photoionization radiation pressure
  call radpress_source(i,j,k,x,y,z,r,prim,s)
#endif
#ifdef EIGHT_WAVE
 
  !  divergence correction Powell et al. 1999
  call divbcorr_source(i,j,k,prim,s)
#endif

  
  return
end subroutine source

!=======================================================================
  
end module sources

#endif

!=======================================================================
