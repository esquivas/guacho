!=======================================================================
!> @file snr.f90
!> @brief Exoplanet problem module
!> @author M. Schneiter, C. Villarreal  D'Angelo, A. Esquivel
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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief super nova remnant module
!> @details Module to impose a SN explosion

module snr

  use constants, only : pc, Msun
  implicit none
  !   SN parameters
  real :: Esn = 1.e51    !<  Energy in the SN
  real :: Rsn = 1.0*pc   !<  Initial radius of the SN
  real :: Msn = 2.0*Msun !<  Mass inside Rsn
  real, parameter :: chi =0.5       !<  Fraction of kinetic to total energy
  real, parameter :: mu_sn = 0.61   !<  Mean atomic mass of the SN

contains

!=======================================================================

!> @brief Module initialization
!> @details Here the parameters of the SNR are initialized, and scaled 
!! to code units

subroutine init_snr()
  use parameters, only: rhosc, psc, rsc
  implicit none

  Esn = Esn / Psc
  Rsn = Rsn / rsc
  Msn = Msn * rsc**3 / rhosc

end subroutine init_snr

!=======================================================================

!> @brief  detonate SN
!> @details Imposes a SN explosion centered at entered at (xc, yc, zc)
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserver variables
!> @param real [in] xc, yc, zc : coordinates of the center of the SN
!--------------------------------------------------------------------

subroutine impose_snr(u, xc, yc, zc)
  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         nxtot, nytot, nztot, neq, nx, ny, nz
  use globals,    only : coords, dx, dy, dz
  use constants,  only : pi 
  implicit none
  real, intent(inout) :: u(neq,nxmin:nxmax,nymin:nymax)
  real, intent(in)    :: xc, yc, zc
  integer :: i, j, k
  real    :: x, y, z, r,  dens, vr, Eth, Mtot
  
  do k=nzmin, nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax 

        ! Position measured from the centre of the grid
        x = ( float( i+coords(0)*nx-nxtot/2 )+0.5 ) * dx
        y = ( float( j+coords(1)*ny-nytot/2 )+0.5 ) * dy
        z = ( float( k+coords(2)*nz-nztot/2 )+0.5 ) * dz

        !  radius from the SN center
        r = sqrt( (x-xc)**2 + (y-yc)**2 + (z-zc)**2 )

        if (r <= Rsn) then

        !   inside SN
        dens=(3./4./pi)*Msn/(Rsn**3)* mu_sn
        vr=(R/Rsn)*Sqrt(10.*chi*Esn/(3*Msn))
        Eth=(3./4./pi)*Esn*(1.-chi)/(Rsn**3)

        u(1,i,j)= dens
        u(2,i,j)= dens*vr*(x-xc)/(r)
        u(3,i,j)= dens*vr*(y-yc)/(r)
        u(4,i,j)= dens*vr*(z-yc)/(r)
        u(5,i,j)= 0.5*dens*vr**2+Eth

        endif
      end do
    end do
  end do

end subroutine impose_snr

!=======================================================================

end module snr


!=======================================================================