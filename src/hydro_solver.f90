!=======================================================================
!> @file hydro_solver.f90
!> @brief Hydrodynamical and Magnetohidrodynamocal solver module
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

!> @brief Advances the simulation one timestep
!> @details Advances the solution from @f$ t @f$ to @f$ t + \Delta t @f$

module hydro_solver

#ifdef HLL
 use hll
#endif
#ifdef HLLC
  use hllc
#endif
#ifdef HLLE
  use hllE
#endif
#ifdef HLLD
  use hlld
#endif
#ifdef EOS_CHEM
  use chemistry
#endif
  implicit none

contains

!> @brief Adds artificial viscosity to the conserved variables
!> @details Adds artificial viscosity to the conserved variables
!! @n Takes the variables from the globals module and it assumes
!! that the up are the stepped variables, while u are unstepped
subroutine viscosity()

  use parameters, only : nx, ny, nz, eta
  use globals, only: u, up
  implicit none
  integer :: i, j, k
  
  do k=1,nz
     do j=1,ny
        do i=1,nx
           up(:,i,j,k)=up(:,i,j,k)+eta*( u(:,i+1,j,k)+u(:,i-1,j,k)       &
                                        +u(:,i,j+1,k)+u(:,i,j-1,k)       &
                                        +u(:,i,j,k+1)+u(:,i,j,k-1)       &
                                     -6.*u(:,i,j,k) )
        end do
     end do
  end do
  
end subroutine viscosity
!=======================================================================
#ifdef CT
subroutine current()

  use parameters, only : nx, ny, nz
  use globals, only :  f, g, h, e
  use boundaries, only: boundaryI_ct

  implicit none

  integer :: i, j, k

  do k=1,nz!+1
     do j=1,ny!+1
        do i=1,nx!+1
           !
! Determination of electric field (E= v x B) 
!         i=max(i,0)            
!         j=max(j,0)            
!         k=max(k,0)            
!         i=min(i,nx+1)
!         j=min(j,ny+1)
!         k=min(k,nz+1)


           e(1,i,j,k)=0.25*(-g(8,i,j-1,k)-g(8,i,j,k) &
                              +h(7,i,j,k-1)+h(7,i,j,k))
           e(2,i,j,k)=0.25*(+f(8,i-1,j,k)+f(8,i,j,k) &
                              -h(6,i,j,k-1)-h(6,i,j,k))
           e(3,i,j,k)=0.25*(-f(7,i-1,j,k)-f(7,i,j,k) &
                              +g(6,i,j-1,k)+g(6,i,j,k)) 

        end do
     end do
  end do

  call boundaryI_ct()  

end subroutine current
#endif
!=======================================================================
!> @brief Upwind timestep
!> @details Performs the upwind timestep according to
!! @f[ U^{n+1}_i= U^n_i -\frac{\Delta t}{\Delta x} 
!!\left[F^{n+1/2}_{i+1/2}-F^{n+1/2}_{i-1/2} \right] @f]
!! (in 3D), it takes @f$ U^{n+1} @f$=up from the global variables
!! and @f$ U^{n} @f$=u
!> @param real [in] dt : timestep

subroutine step(dt)

  use parameters, only : nx, ny, nz, neqdyn
#ifdef CT
  use globals, only:  up, u, primit, f, g, h, dx, dy, dz, e
#else
  use globals, only : up, u, primit, f, g, h, dx, dy, dz
#endif 
#if defined(GRAV) || defined(RADPRES) || defined(EIGHT_WAVE)
  use sources
#endif
  implicit none
#if defined(GRAV) || defined(RADPRES) || defined(EIGHT_WAVE)
  real :: s(neq)
#endif
  real, intent(in) :: dt
  integer :: i, j, k
  real :: dtdx, dtdy, dtdz

  dtdx=dt/dx
  dtdy=dt/dy
  dtdz=dt/dz

  do k=1,nz
     do j=1,ny
        do i=1,nx

           up(:5,i,j,k)=u(:5,i,j,k)-dtdx*(f(:5,i,j,k)-f(:5,i-1,j,k))    &
                                 -dtdy*(g(:5,i,j,k)-g(:5,i,j-1,k))    &
                                 -dtdz*(h(:5,i,j,k)-h(:5,i,j,k-1))
#ifdef MHD
#ifdef CT                                  
           up(6,i,j,k)=u(6,i,j,k)-0.5*dtdy*(e(3,i,j+1,k)-e(3,i,j-1,k))    &
                +0.5*dtdz*(e(2,i,j,k+1)-e(2,i,j,k-1))
           
           up(7,i,j,k)=u(7,i,j,k)+0.5*dtdx*(e(3,i+1,j,k)-e(3,i-1,j,k))    &
                -0.5*dtdz*(e(1,i,j,k+1)-e(1,i,j,k-1))
           
           up(8,i,j,k)=u(8,i,j,k)-0.5*dtdx*(e(2,i+1,j,k)-e(2,i-1,j,k))    &
                              +0.5*dtdy*(e(1,i,j+1,k)-e(1,i,j-1,k))

#else
           up(6:8,i,j,k)=u(6:8,i,j,k)-dtdx*(f(6:8,i,j,k)-f(6:8,i-1,j,k))    &
                                 -dtdy*(g(6:8,i,j,k)-g(6:8,i,j-1,k))    &
                                 -dtdz*(h(6:8,i,j,k)-h(6:8,i,j,k-1))
#endif
#endif

#ifdef PASSIVES
           up(neqdyn+1:,i,j,k)=u(neqdyn+1:,i,j,k)-dtdx*(f(neqdyn+1:,i,j,k)-f(neqdyn+1:,i-1,j,k))    &
                                 -dtdy*(g(neqdyn+1:,i,j,k)-g(neqdyn+1:,i,j-1,k))    &
                                 -dtdz*(h(neqdyn+1:,i,j,k)-h(neqdyn+1:,i,j,k-1))
#endif
#if defined(GRAV) || defined(RADPRES) || defined(EIGHT_WAVE)
           call source(i,j,k,primit(:,i,j,k),s)
           up(:,i,j,k)= up(:,i,j,k)+dt*s(:)
#endif

        end do
     end do
  end do

end subroutine step


!> @brief High level wrapper to advancce the simulation
!> @details High level wrapper to advancce the simulation
!! @n The variables are taken from the globals module.

subroutine tstep()

  use parameters, only : tsc
  use globals
  use hydro_core, only : calcprim
  use boundaries
#ifdef C2ray
  use C2Ray_RT, only: C2Ray_radiative_transfer
#endif
#ifdef COOLINGH
  use cooling_H
#endif
#ifdef COOLINGDMC
  use cooling_DMC
#endif
#ifdef COOLINGCHI
  use cooling_CHI
#endif
#ifdef RADDIFF
  use difrad
#endif

#ifdef THERMAL_COND
  use thermal_cond
#endif
  implicit none
  real :: dtm
   
  !  1st half timestep ========================   
  dtm=dt_CFL/2.
  !   calculate the fluxes using the primitives
  !   (piecewise constant)
#ifdef HLL
  call hllfluxes(1)
#endif
#ifdef HLLC
  call hllcfluxes(1)
#endif
#ifdef HLLE
  call hllEfluxes(1)
#endif
#ifdef HLLD
  call hlldfluxes(1)
#endif

  !calculates the electric current for CT
#ifdef CT
  call current()
#endif  
  !   upwind timestep
  call step(dtm)
  
  !   add viscosity
  !call viscosity()
  
  !  2nd half timestep ========================
  !  boundaries in up and  primitives up ---> primit
  call boundaryII()
  call calcprim(up,primit)

  !   calculate the fluxes using the primitives
  !   with linear reconstruction (piecewise linear)
#ifdef HLL
  call hllfluxes(2)
#endif
#ifdef HLLC
  call hllcfluxes(2)
#endif
#ifdef HLLE
  call hllEfluxes(2)
#endif
#ifdef HLLD
  call hlldfluxes(2)
#endif

  !calculates the electric current for CT
#ifdef CT
  call current()
#endif
  !  upwind timestep
  call step(dt_CFL)

  !  add viscosity
  call viscosity()
  
  !  copy the up's on the u's
  u=up

  ! update the chemistry network
  ! at this point is in cgs
#ifdef EOS_CHEM
  !  the primitives in the physical celles are upated
  call update_chem()
#endif

 !  Do the Radiaiton transfer (Monte Carlo type)
#ifdef RADDIFF
  call diffuse_rad()
#endif

#ifdef CEXCHANGE
 !****************************************************
  !   apply charge exchange TO BE ADDED IN COLLING?
  !   not fully implemented
  !****************************************************
  call cxchange(dt*tsc)
#endif

  !   apply cooling/heating
#ifdef COOLINGH
  !   add cooling to the conserved variables
  call coolingh()
   !  update the primitives with u
  call calcprim(u, primit)
#endif
#ifdef COOLINGDMC
  !   the primitives are updated in the cooling routine
  call coolingdmc()
#endif
#ifdef COOLINGCHI
  !   the primitives are updated in the cooling routine
  call coolingchi()
#endif
#ifdef COOLINGCHEM
  !the primitives are already updated in update_chem
  call cooling_chem()
#endif

#ifdef C2ray
  !  Apply Rad transfer, and heating & Cooling w/C2-Ray
  call C2ray_radiative_transfer(dt_CFL*tsc)
#endif

  !   boundary contiditions on u
  call boundaryI()

  !  update primitives on the boundaries
  call calcprim(u,primit,only_ghost=.true.)


#ifdef THERMAL_COND
  !  Thermal conduction
  call thermal_conduction()
#endif

end subroutine tstep

end module hydro_solver

!=======================================================================
