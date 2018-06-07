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
!> (modified to impose two jets at once)

module jet

  use parameters
  implicit none
  integer, parameter :: njets = 2
  real :: Rj(njets), Lj(njets), denj(njets), Tempj(njets), vj0(njets), &
         dvj(njets), tau(njets), omega(njets)
  real :: posj(njets,3)
  !  the direction can be obtained with the following parameters
  !  alpha is the angle with respect to z at t=0
  !  the angle can be adjusted by rotating it by a precesion
  !  angle (omegaP x t)
  real, save :: alpha(njets), omegaP(njets)

contains

  !--------------------------------------------------------------------
  !   Here the parameters of the jet are initialized, and scaled to
  !   code units
  !--------------------------------------------------------------------

  subroutine init_jet()

    use constants, only : au, pi, deg
    implicit none

    Rj(:)    = 400.*au/rsc   !  jet radius
    Lj(:)    = 400.*au/rsc   !  jet length

    !  jet position(s)
    posj(1,1)= -5.e3*au /rsc
    posj(1,2)=  5.e3*au /rsc
    posj(1,3)= 0.!e3*au /rsc

    posj(2,1)=  5e3*au /rsc
    posj(2,2)=  5e3*au /rsc
    posj(2,3)= 0.!e3*au /rsc

    !  jet orientation parameters
    alpha(1) = 30.*deg
    alpha(2) =-10.*deg
    omegaP(1)=2.*pi/(2000.*yr/tsc)
    omegaP(2)=2.*pi/(2000.*yr/tsc)

    !  outflow parameters
    denj(:)  = 300.                      !  density
    Tempj(:) = 1000./Tempsc              !  jet temperature
    vj0(:)   = 200.e5/vsc                !  mean velocity
    dVj(:)   = 0.!(200./3.)*1e5/vsc      !  amplitude of variability
    tau(:)   = 500.*yr/tsc               !  period of variability
    omega(:) = 2.*pi/tau

  end subroutine init_jet

  !--------------------------------------------------------------------

  subroutine impose_jet(u,time)
    use globals, only : coords, dx, dy, dz
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent (in)   :: time
    real :: x, y, z, rad, xp, yp, zp, rx, ry, rz, vjet
    !  precesion opening angle (or initial direction, repect to the z axis)
    real :: omegat(njets)
    real :: sina(njets), cosa(njets)
    real :: coso(njets), sino(njets)
    integer ::  i,j,k, nj

    do nj=1,njets
      sina(nj)= sin(alpha(nj))
      cosa(nj)= cos(alpha(nj))

      sino(nj)= sin(omegaP(nj)*time)
      coso(nj)= cos(omegaP(nj)*time)

      omegat(nj) = omega(nj)*time   ! for jet variability
    end do

    do i=nxmin,nxmax
      do j=nymin,nymax
        do k=nzmin,nzmax
          !   measured from the center of the computational mesh
          x=(float(i+coords(0)*nx-nxtot/2) - 0.5)*dx
          y=(float(j+coords(1)*ny-nytot/2) - 0.5)*dy
          z=(float(k+coords(2)*nz-nztot/2) - 0.5)*dz

          do nj = 1, njets

            xp=x-posj(nj,1)
            yp=y-posj(nj,2)
            zp=z-posj(nj,3)

            rx= xp*coso(nj)          - yp*sino(nj)
            ry= xp*cosa(nj)*sino(nj) + yp*cosa(nj)*coso(nj) - zp*sina(nj)
            rz= xp*sina(nj)*sino(nj) + yp*sina(nj)*coso(nj) + zp*cosa(nj)

            rad=sqrt(rx**2+ry**2)

            !if( (j.eq.0).and.(i.eq.0).and.(rank.eq.0)) print*,k,z,zp

            if( (abs(rz) <= Lj(nj)).and.(rad <= Rj(nj)) ) then
              !  inside the jet source
              vjet= vj0(nj) + dvj(nj)*sin(omegat(nj))
              vjet=sign(vjet,rz)
              !
              !   total density and momenta
              u(1,i,j,k) = denj(nj)
              u(2,i,j,k) = denj(nj)*vjet*sina(nj)*coso(nj)
              u(3,i,j,k) = denj(nj)*vjet*sina(nj)*sino(nj)
              u(4,i,j,k) = denj(nj)*vjet*cosa(nj)
              !   energy
              u(5,i,j,k)=0.5*denj(nj)*vjet**2 + cv*denj(nj)*Tempj(nj)
              !  passive scalars needed for the chemistry network
              u(1,i,j,k) = denj(nj)
            end if

          end do

        end do
      end do
    end do

  end subroutine impose_jet
  !--------------------------------------------------------------------
end module jet
