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
  
  do i=1,nx
     do j=1,ny
        do k=1,nz
           up(:,i,j,k)=up(:,i,j,k)+eta*( u(:,i+1,j,k)+u(:,i-1,j,k)       &
                                        +u(:,i,j+1,k)+u(:,i,j-1,k)       &
                                        +u(:,i,j,k+1)+u(:,i,j,k-1)       &
                                     -6.*u(:,i,j,k) )
        end do
     end do
  end do
  
end subroutine viscosity

!=======================================================================

!> @brief Upwind timestep
!> @details Performs the upwind timestep according to
!! @f[ U^{n+1}_i= U^n_i -\frac{\Delta t}{\Delta x} 
!!\left[F^{n+1/2}_{i+1/2}-F^{n+1/2}_{i-1/2} \right] @f]
!! (in 3D), it takes @f$ U^{n+1} @f$=up from the global variables
!! and @f$ U^{n} @f$=u
!> @param real [in] dt : timestep

subroutine step(dt)
  use parameters, only : nx, ny, nz
  use globals, only : up, u, primit, f, g, h, dx, dy, dz
#ifdef SOURCE
  use sources
#endif
  implicit none
#ifdef SOURCE
  real :: s(neq)
#endif
  real, intent(in) :: dt
  integer :: i, j, k
  real :: dtdx, dtdy, dtdz

  dtdx=dt/dx
  dtdy=dt/dy
  dtdz=dt/dz

  do i=1,nx
     do j=1,ny
        do k=1,nz
           !
           up(:,i,j,k)=u(:,i,j,k)-dtdx*(f(:,i,j,k)-f(:,i-1,j,k))    &
                                 -dtdy*(g(:,i,j,k)-g(:,i,j-1,k))    &
                                 -dtdz*(h(:,i,j,k)-h(:,i,j,k-1))

#ifdef SOURCE
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
!> @param real [in] time :  integration time
!> @param real [in] dt   :  timestep

subroutine tstep(time,dt)

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
  real, intent(in):: dt, time
  real            :: dtm
   
  !  1st half timestep ========================   
  dtm=dt/2.
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

  !   upwind timestep
  call step(dtm)
  
  !   add viscosity
  !call viscosity()
  
  !  2nd half timestep ========================
  !  boundaries in up and  primitives up ---> primit
  call boundaryII(time,dt)
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

  !  upwind timestep
   call step(dt)

  !  add viscosity
  call viscosity()
  
  !  copy the up's on the u's
  u=up
 
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
  call coolingh(dt*tsc)
#endif
#ifdef COOLINGDMC
  !   the primitives are updated in the cooling routine
  call coolingdmc(dt*tsc)
#endif
#ifdef COOLINGCHI
  !   the primitives are updated in the cooling routine
  call coolingchi(dt*tsc)
#endif
  !   BBC cooling not implemented yet

  !  update the primitives with u
  call calcprim(u, primit)
  

#ifdef C2ray
  !  Apply Rad transfer, and heating & Cooling w/C2-Ray
  call C2ray_radiative_transfer(dt*tsc)
#endif

  !   boundary contiditions on u
  call boundaryI(time,dt)

#ifdef THERMAL_COND
  !  Thermal conduction
  call thermal_conduction(dt*tsc)
#endif

end subroutine tstep

end module hydro_solver

!=======================================================================
