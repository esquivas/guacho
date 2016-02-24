!=======================================================================
!> @file jet.f90
!> @brief jet module
!> @author A. Esquivel
!> @date 24/Nov/2014

! Copyright (c) 2014 A. Esquivel et al.
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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief jet module
!> @details Module to impose a jet with precesion and variability

module jet

  use parameters
  implicit none
  real, save :: Rj, Lj, denj, Tempj, vj0, dvj, tau, omega 
  real, save :: posj(3)
  !  the direction can be obtained with the following parameters
  !  alpha is the angle with respect to z at t=0
  !  the angle can be adjusted by rotating it by a precesion
  !  angle (omegaP x t)
  real, save :: alpha, omegaP

contains

  !--------------------------------------------------------------------
  !   Here the parameters of the jet are initialized, and scaled to
  !   code units
  !--------------------------------------------------------------------

  subroutine init_jet()

    use constants, only : au, pi
    implicit none
    
    Rj    = 400.*au/rsc   !  jet radius
    Lj    = 400.*au/rsc   !  jet length
    
    !  jet position
    posj(1)= 7.5e3*au /rsc
    posj(2)= 7.5e3*au /rsc
    posj(3)= 0.!e3*au /rsc
    
    !  jet orientation parameters
    alpha =6.*pi/180.
    omegaP=2.*pi/(2142.*yr/tsc)

    denj  = 3.e3                     !  density
    Tempj = 100.                      !  jet temperature
    vj0   = 200.e5/vsc                !  mean velocity
    dVj   = (200./3.)*1e5/vsc         !  amplitude of variability
    tau   = 535.*yr/tsc               !  period of variability
    omega = 2.*pi/tau                 !  initial position
    

  end subroutine init_jet
 
  !--------------------------------------------------------------------
 
  subroutine impose_jet(u,time)
    use globals, only : coords, dx, dy, dz, rank
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent (in)   :: time
    real :: omegat, x, y, z, rad, xp, yp, zp, rx, ry, rz, vjet
    !  precesion opening angle (or initial direction, repect to the z axis)
    real :: sina, cosa
    real :: coso, sino
    integer ::  i,j,k


    sina= sin(alpha)
    cosa= cos(alpha)
    
    sino= sin(omegaP*time)
    coso= cos(omegaP*time)

    omegat = omega*time   ! for jet variability

    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax
           
             !   measured from the corner of the computational mesh
             x=(float(i+coords(0)*nx) - 0.5)*dx
             y=(float(j+coords(1)*ny) - 0.5)*dy
             z=(float(k+coords(2)*nz) - 0.5)*dz

             xp=x-posj(1)
             yp=y-posj(2)
             zp=z-posj(3)

             rx= xp*coso      - yp*sino
             ry= xp*cosa*sino + yp*cosa*coso - zp*sina
             rz= xp*sina*sino + yp*sina*coso + zp*cosa
             
             rad=sqrt(rx**2+ry**2)          

             !if( (j.eq.0).and.(i.eq.0).and.(rank.eq.0)) print*,k,z,zp

             if( (abs(rz) <= Lj).and.(rad <= Rj) ) then

                !  inside the jet source
                vjet= vj0 + dvj*sin(omegat)
                !vjet=sign(vjet,rz)
                !
                !   total density and momenta
                u(1,i,j,k) = denj
                u(2,i,j,k) = denj*vjet*sina*coso
                u(3,i,j,k) = denj*vjet*sina*sino
                u(4,i,j,k) = denj*vjet*cosa
                !   energy
                u(5,i,j,k)=0.5*denj*vjet**2+cv*denj*1.1*Tempj/Tempsc
                
                !  passive scalars needed for the chemistry network
                u( 6,i,j,k) = 0.0001 * u(1,i,j,k)
                u( 7,i,j,k) = 0.1/2. * u(1,i,j,k)
                u( 8,i,j,k) = 0.8999 * u(1,i,j,k)
                u( 9,i,j,k) = 0.0001 * u(1,i,j,k)
                u(10,i,j,k) = + u(1,i,j,k)

             endif
             
          end do
       end do
    end do

  end subroutine impose_jet
  !--------------------------------------------------------------------
end module jet
