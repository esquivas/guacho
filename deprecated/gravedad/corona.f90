!=======================================================================
!> @file corona.f90
!> @brief corona module
!> @author V. Sieyra, M. Schneiter & A. Esquivel
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

!> @brief corona module
!> @details Module to impose a corona with precesion and variability

module corona

  use parameters
  implicit none
  real, save :: nc, rhoc, P_0, Temp1, Temp2
!   real, save :: bc1x, bc1y, bc1z, ec1
!   real, save :: bc0x, bc0y, bc0z, ec0

contains

  !--------------------------------------------------------------------
  !   Here the parameters of the corona are initialized, and scaled to
  !   code units
  !--------------------------------------------------------------------

  subroutine init_corona()

    use constants, only : amh, pi, kb
    implicit none
    
      nc = 1.e20
   !   rhoc = mu_1*amh*nc
      Temp1 = 1.e4 ! cuando mu=1 no puede valer 1.e6 porque la raiz de la densidad se hace negativa
      P_0 = nc*kb*Temp1
      Temp2 = 1.e6

      !       print*,'rhoc corona',rhoc
  end subroutine init_corona
 
  !--------------------------------------------------------------------
 
  subroutine driver(u,time) ! driver Murawski, 2014
  
  use constants, only : pi
  use globals, only : coords, dx, dy, dz, rank
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: Av, w, xc, yc, x, y, Pd
  integer :: i, j

!     Av=3.e5/vsc ! 3km/s vsc = 1664782 cm/s = 16 km/s
!     Pd=300./tsc
!     w=1.e8/rsc  
!     xc = dx*(nxtot*0.5-0.5) 
!     yc = 0.
! 
!       do i = nxmin,nxmax
!         x = (real(i+coords(0)*nx) - 0.5)*dx
!         y = (real(j+coords(1)*ny) - 0.5)*dy
!         if (coords(1).eq.0) then
!           u(3,i,0:10,:)=u(1,i,0:10,:)*Av*exp(-((x-xc)**2.+(y-yc)**2.)/w**2.)!*sin(2.*pi*time/Pd)
!           u(5,i,0:10,:) = cv*u(1,i,0:10,:)*Tempc + 0.5*u(3,i,0:10,:)*u(3,i,0:10,:)/u(1,i,0:10,:) &
!           + 0.5*(u(6,i,0:10,:)*u(6,i,0:10,:)+u(7,i,0:10,:)*u(7,i,0:10,:)+u(8,i,0:10,:)*u(8,i,0:10,:))
!         end if
!       end do 
  end subroutine driver
  
!   subroutine impose_jet(u,time)
!     use globals, only : coords, dx, dy, dz, rank
!     implicit none
!     real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
!     real, intent (in)   :: time
!     real :: omegat, x, y, z, rad, xp, yp, zp, rx, ry, rz, vjet
!     !  precesion opening angle (or initial direction, repect to the z axis)
!     real :: sina, cosa
!     real :: coso, sino
!     integer ::  i,j,k
! 

!     sina= sin(alpha)
!     cosa= cos(alpha)
!     
!     sino= sin(omegaP*time)
!     coso= cos(omegaP*time)
! 
!     omegat = omega*time   ! for corona variability
! 
!     do i=nxmin,nxmax
!        do j=nymin,nymax
!           do k=nzmin,nzmax
!            
!              !   measured from the corner of the computational mesh
!              x=(real(i+coords(0)*nx) - 0.5)*dx
!              y=(real(j+coords(1)*ny) - 0.5)*dy
!              z=(real(k+coords(2)*nz) - 0.5)*dz
! 
!              xp=x-posj(1)
!              yp=y-posj(2)
!              zp=z-posj(3)
! 
!              rx= xp*coso      - yp*sino
!              ry= xp*cosa*sino + yp*cosa*coso - zp*sina
!              rz= xp*sina*sino + yp*sina*coso + zp*cosa
!              
!              rad=sqrt(rx**2+ry**2)          
! 
!              !if( (j.eq.0).and.(i.eq.0).and.(rank.eq.0)) print*,k,z,zp
! 
!              if( (abs(rz) <= Lj).and.(rad <= Rj) ) then
! 
!                 !  inside the corona source
!                 vjet= vj0 + dvj*sin(omegat)
!                 !vjet=sign(vjet,rz)
!                 !
!                 !   total density and momenta
!                 u(1,i,j,k) = denj
!                 u(2,i,j,k) = denj*vjet*sina*coso
!                 u(3,i,j,k) = denj*vjet*sina*sino
!                 u(4,i,j,k) = denj*vjet*cosa
!                 !   energy
!                 u(5,i,j,k)=0.5*denj*vjet**2+cv*denj*0.9501*Tempj/Tempsc
!                 
!                 !  passive scalars needed for the chemistry network
!                 u( 6,i,j,k) = 0.9989 * u(1,i,j,k)
!                 u( 7,i,j,k) = 0.0001 * u(1,i,j,k)
!                 u( 8,i,j,k) = 0.001/2. * u(1,i,j,k)
!                 u( 9,i,j,k) = 0.0001 * u(1,i,j,k)
!                 u(10,i,j,k) = + u(1,i,j,k)
! 
!              endif
!              
!           end do
!        end do
!     end do

!   end subroutine impose_jet
  !--------------------------------------------------------------------

end module corona

!=======================================================================
