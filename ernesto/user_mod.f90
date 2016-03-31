!=======================================================================
!> @file user_mod.f90
!> @brief User input module
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

!> @brief User imput module
!> @details  This is an attempt to have all input neede from user in a
!! single file
!!!  This module should load additional modules (i.e. star, jet, sn), to
!!  impose initial and boundary conditions (such as sources)

module user_mod

  ! load auxiliary modules
  implicit none
 
contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty 
subroutine init_user_mod()

  implicit none      
  !  initialize modules loaded by user

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables

subroutine initial_conditions(u)

  use parameters,only: neq,neqdyn,nxmin,nxmax,nymin,nymax,nzmin,nzmax,nx,ny,nz, & 
                       nxtot,nytot,nztot,gamma,mu,rhosc,tempsc,vsc,Psc,Esc,Bsc
  use constants,only: amh,Rg
  use globals,only: coords,dx,dy,dz,rank, time

  implicit none
  real,intent(out):: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer:: i,j,k
  real:: x,y,z, t
  real:: dens,pres,velx,vely,velz, B0
  real:: dens_phys,temp_phys,pres_phys,velx_phys,vely_phys,velz_phys,B0_phys
  real:: pulso=1.1,r_pulso=0.1,dist
  
  t = time

  velx_phys = 0.   
  vely_phys = 0.
  velz_phys = 0.

  velx=velx_phys/vsc
  vely=vely_phys/vsc
  velz=velz_phys/vsc

  dens_phys = 2.1e-15
  temp_phys = 3e+6
  B0_phys = 5.
  B0 = B0_phys/Bsc
  
  do k=nzmin,nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax

        ! Position measured from the centre of the grid
        x = (float(i+coords(0)*nx-nxtot/2)+0.5)*dx
        y = (float(j+coords(1)*ny-nytot/2)+0.5)*dy
        z = (float(k+coords(2)*nz-nztot/2)+0.5)*dz

        dist=sqrt(x*x+y*y+z*z) 
        if (dist <= r_pulso) then 

          dens = dens_phys*pulso/rhosc 

          pres_phys = Rg/mu*dens_phys*temp_phys
          pres = pres_phys/Psc

          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz
          u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) &
           + pres/(gamma-1.0) + 0.5*B0**2
          u(6,i,j,k) = 0.
          u(7,i,j,k) = B0
          u(8,i,j,k) = 0.

        else
     
          dens = dens_phys/rhosc 

          pres_phys = Rg/mu*dens_phys*temp_phys
          pres = pres_phys/Psc

          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz
          u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) &
          + pres/(gamma-1.0) + 0.5*B0**2
          u(6,i,j,k) = 0.
          u(7,i,j,k) = B0
          u(8,i,j,k) = 0.

        endif

      enddo
    enddo
  enddo


end subroutine initial_conditions
  
!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables

#ifdef OTHERB
subroutine impose_user_bc(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals,    only : time   
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  
end subroutine impose_user_bc
#endif

!=======================================================================


end module user_mod

!=======================================================================
