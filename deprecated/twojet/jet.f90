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
  real :: alpha(njets), omegaP(njets), phiJ(njets), theta(njets)

contains

  !--------------------------------------------------------------------
  !   Here the parameters of the jet are initialized, and scaled to
  !   code units
  !--------------------------------------------------------------------

  subroutine init_jet()

    use constants, only : au, pi, deg
    implicit none

    Rj(1)    = 2*300.*au/rsc   !  jet radius
    Rj(2)    = 2*300.*au/rsc   !  jet radius
    Lj(1)    = 2*500.*au/rsc   !  jet length
    Lj(2)    = 2*500.*au/rsc   !  jet length
    !  jet position(s)
    posj(1,1)= -1560*au /rsc
    posj(1,2)=  0.!*au /rsc
    posj(1,3)=  -140*au /rsc

    posj(2,1)= 1560*au /rsc
    posj(2,2)=-1000*au /rsc
    posj(2,3)=  140*au /rsc

    !  jet orientation parameters
    alpha(1) =   7.*deg
    alpha(2) =  29.*deg
    omegaP(1)=2.*pi/(200.*yr/tsc)
    omegaP(2)=2.*pi/(200.*yr/tsc)
    phiJ(1)   = 30.*deg
    phiJ(2)   =-60.*deg
    theta(1)  = 15.*deg
    theta(2)  = 15.*deg
    !  outflow parameters
    denj(:)  = 10.                      !  density
    Tempj(1) = 1000./Tempsc              !  jet temperature
    Tempj(2) = 1000./Tempsc              !  jet temperature
    vj0(:)   = 300.e5/vsc                !  mean velocity
    dVj(:)   = 0.!(200./3.)*1e5/vsc      !  amplitude of variability
    tau(:)   = 500.*yr/tsc               !  period of variability
    omega(:) = 2.*pi/tau

  end subroutine init_jet

  !--------------------------------------------------------------------

  subroutine impose_jet(u,time)
    use globals, only : coords, dx, dy, dz
    use constants, only : deg
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent (in)   :: time
    real :: x, y, z, rad, xp, yp, zp, rx, ry, rz, vjet, radC
    !  precesion opening angle (or initial direction, repect to the z axis)
    real :: omegat(njets)
    real :: sina(njets), cosa(njets)
    real :: coso(njets), sino(njets)
    integer ::  i,j,k,nj

    do nj=1,njets
      sina(nj)= sin( alpha(nj) )
      cosa(nj)= cos( alpha(nj) )
      sino(nj)= sin( phiJ(nj) )!sin(omegaP(nj)*time + phij(nj) )
      coso(nj)= cos( phiJ(nj) )!cos(omegaP(nj)*time + phij(nj) )

      omegat(nj) = omega(nj)*time   ! for jet variability
    end do

    do i=nxmin,nxmax
      do j=nymin,nymax
        do k=nzmin,nzmax
          !   measured from the center of the computational mesh
          x=(real(i+coords(0)*nx-nxtot/2) - 0.5)*dx
          y=(real(j+coords(1)*ny-nytot/2) - 0.5)*dy
          z=(real(k+coords(2)*nz-nztot/2) - 0.5)*dz

          do nj =1,2

            xp=x-posj(nj,1)
            yp=y-posj(nj,2)
            zp=z-posj(nj,3)

            rx =   xp*cosa(nj)*coso(nj) + yp*cosa(nj)*sino(nj) - zp*sina(nj)
            ry = - xp*sino(nj)          + yp*coso(nj)
            rz =   xp*sina(nj)*coso(nj) + yp*sina(nj)*sino(nj) + zp*cosa(nj)

            rad=sqrt( rx**2 + ry**2 )

            !if( (j.eq.0).and.(i.eq.0).and.(rank.eq.0)) print*,k,z,zp

            if( (abs(rz) <= Lj(nj)).and.(rad <= Rj(nj)) ) then

            !  !  inside the jet source
            !  vjet= vj0(nj) + dvj(nj)*sin(omegat(nj))
            !  vjet=sign(vjet,rz)

            !  if (nj <= 2) then
            !   !   total density and momenta
            !    u(1,i,j,k) =   denj(nj)
            !    u(2,i,j,k) =   denj(nj)*vjet*sina(nj)*coso(nj)
            !    u(3,i,j,k) =   denj(nj)*vjet*sina(nj)*sino(nj)
            !    u(4,i,j,k) =   denj(nj)*vjet*cosa(nj)
            !    !   energy
            !    u(5,i,j,k)=0.5*denj(nj)*vjet**2 + cv*denj(nj)*Tempj(nj)
            !    !  passive scalars needed for the chemistry network
            !    u(6,i,j,k) = denj(nj)

            !  else

                if ( atan(rad/abs(rz)) > theta(nj) ) then
                  vjet=0.
                else
                  vjet= vj0(nj) + dvj(nj)*sin(omegat(nj))
                end if
                radC = sqrt( xp**2 + yp**2 + zp**2 )
                !   total density and momenta
                u(1,i,j,k) = denj(nj)
                u(2,i,j,k) = denj(nj)*vjet*xp/radC
                u(3,i,j,k) = denj(nj)*vjet*yp/radC
                u(4,i,j,k) = denj(nj)*vjet*zp/radC
                !   energy
                u(5,i,j,k)=0.5*denj(nj)*vjet**2 + cv*denj(nj)*Tempj(nj)
                !  passive scalars needed for the chemistry network
                u(6,i,j,k) = denj(nj)

              endif

            !endif

          end do

        end do
      end do
    end do

  end subroutine impose_jet

  !--------------------------------------------------------------------
end module jet
