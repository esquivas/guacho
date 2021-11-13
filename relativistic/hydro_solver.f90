!=======================================================================
!> @file hydro_solver.f90
!> @brief Hydrodynamical and Magnetohidrodynamocal solver module
!> @author Alejandro Esquivel
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

!> @brief Advances the simulation one timestep
!> @details Advances the solution from @f$ t @f$ to @f$ t + \Delta t @f$

module hydro_solver

  use hll
  use hllr
  use hllc
  use hllcr
  use hllE
  use hlld
  use hlleSplitAll
  use chemistry
  implicit none

contains

!> @brief Adds artificial viscosity to the conserved variables
!> @details Adds artificial viscosity to the conserved variables
!> using the advanced solution and copying the result in the U array
!! @n Takes the variables from the globals module and it assumes
!! that the up are the stepped variables, while u are unstepped
!! It only performs the copy in the physical domain

  subroutine viscous_copy()

  use parameters, only : nx, ny, nz, eta
  use globals, only: u, up
  implicit none
  integer :: i, j, k

  do k=1,nz
     do j=1,ny
        do i=1,nx
           u(:,i,j,k)=up(:,i,j,k)+eta*(  up(:,i+1,j,k)+up(:,i-1,j,k)       &
                                        +up(:,i,j+1,k)+up(:,i,j-1,k)       &
                                        +up(:,i,j,k+1)+up(:,i,j,k-1)       &
                                     -6.*up(:,i,j,k) )
        end do
     end do
  end do

end subroutine viscous_copy

!=======================================================================

!> @brief Upwind timestep
!> @details Performs the upwind timestep according to
!! @f[ U^{n+1}_i= U^n_i -\frac{\Delta t}{\Delta x}
!!\left[F^{n+1/2}_{i+1/2}-F^{n+1/2}_{i-1/2} \right] @f]
!! (in 3D), it takes @f$ U^{n+1} @f$=up from the global variables
!! and @f$ U^{n} @f$=u
!> @param real [in] dt : timestep

subroutine step(dt)

  use parameters, only : nx, ny, nz, neqdyn, &
                         user_source_terms, radiation_pressure, &
                         eight_wave, enable_flux_cd

  use globals, only : up, u, primit,f, g, h, dx, dy, dz
  use flux_cd_module
  use sources
  implicit none
  real :: s(neq)
  real, intent(in) :: dt
  integer :: i, j, k
  real :: dtdx, dtdy, dtdz

  dtdx=dt/dx
  dtdy=dt/dy
  dtdz=dt/dz

#ifdef BFIELD
  if (enable_flux_cd) call get_efield()
#endif

  do k=1,nz
    do j=1,ny
      do i=1,nx

        if (.not.enable_flux_cd) then
          !  upwind step for all variables
          up(:,i,j,k)=u(:,i,j,k)-dtdx*(f(:,i,j,k)-f(:,i-1,j,k))     &
                                -dtdy*(g(:,i,j,k)-g(:,i,j-1,k))     &
                                -dtdz*(h(:,i,j,k)-h(:,i,j,k-1))
        else

#ifdef BFIELD
          call flux_cd_update(i,j,k,dt)
#endif
        endif

        if (user_source_terms     .or. &
            radiation_pressure    .or. &
            eight_wave) then
          call source(i,j,k,primit(:,i,j,k),s)
          up(:,i,j,k)= up(:,i,j,k)+dt*s(:)

        end if

        end do
     end do
  end do

end subroutine step

!> @brief High level wrapper to advancce the simulation
!> @details High level wrapper to advancce the simulation
!! @n The variables are taken from the globals module.

subroutine tstep()

  use parameters, only : tsc, riemann_solver, eq_of_state, &
                         dif_rad, cooling, &
                         th_cond
  use constants
  use globals
  use hydro_core, only : calcprim
  use boundaries
  use cooling_H
  use cooling_DMC
  use cooling_CHI
  use difrad
  use thermal_cond
  implicit none
  real :: dtm

  !  1st half timestep ========================
  dtm=dt_CFL/2.
  !   calculate the fluxes using the primitives
  !   (piecewise constant)
  if (riemann_solver == SOLVER_HLL )  call hllfluxes(1)
  if (riemann_solver == SOLVER_HLLC)  call hllcfluxes(1)
  if (riemann_solver == SOLVER_HLLE)  call hllefluxes(1)
  if (riemann_solver == SOLVER_HLLD)  call hlldfluxes(1)
  !if (riemann_solver == SOLVER_HLLE_SPLIT_B) call hllefluxes(1)
  !if (riemann_solver == SOLVER_HLLD_SPLIT_B) call hllefluxes(1)
  if (riemann_solver == SOLVER_HLLE_SPLIT_ALL) call hllefluxesSplitAll(1)
  !if (riemann_solver == SOLVER_HLLD_SPLIT_ALL) call hllefluxes(1)
  if (riemann_solver == SOLVER_RHLL ) call hllrfluxes(1)
  if (riemann_solver == SOLVER_RHLLC) call hllcrfluxes(1)

  !   upwind timestep
  call step(dtm)

  !  2nd half timestep ========================
  !  boundaries in up and  primitives up ---> primit
  call boundaryII()
  call calcprim(up,primit)

  !   calculate the fluxes using the primitives
  !   with linear reconstruction (piecewise linear)
  if (riemann_solver == SOLVER_HLL )  call hllfluxes(2)
  if (riemann_solver == SOLVER_HLLC)  call hllcfluxes(2)
  if (riemann_solver == SOLVER_HLLE)  call hllefluxes(2)
  if (riemann_solver == SOLVER_HLLD)  call hlldfluxes(2)
  !if (riemann_solver == SOLVER_HLLE_SPLIT_B) call hllefluxes(2)
  !if (riemann_solver == SOLVER_HLLD_SPLIT_B) call hllefluxes(2)
  if (riemann_solver == SOLVER_HLLE_SPLIT_ALL) call hllefluxesSplitAll(2)
  !if (riemann_solver == SOLVER_HLLD_SPLIT_ALL) call hllefluxes(2)
  if (riemann_solver == SOLVER_RHLL ) call hllrfluxes(2)
  if (riemann_solver == SOLVER_RHLLC) call hllcrfluxes(2)

  !  upwind timestep
  call step(dt_CFL)

  !  copy the up's on the u's and include viscosity
  !  warning: The copy does not do the ghost cells
  call viscous_copy()

  !  Do the Radiation transfer (Monte Carlo type)
  if (dif_rad) call diffuse_rad()

  ! update the chemistry network
  ! at this point is in cgs
  !  the primitives in the physical domain are upated
  if (eq_of_state == EOS_CHEM) call update_chem()

  !----------------------------------
  !   apply cooling/heating terms

  !   add cooling (H rat e)to the conserved variables
  if (cooling == COOL_H) then
    call coolingh()
  end if

  ! DMC cooling (the primitives are updated in the cooling routine)
  if (cooling == COOL_DMC) call coolingdmc()

  ! Chianti cooling (the primitives are updated in the cooling routine)
  if (cooling == COOL_CHI) call coolingchi()

  ! Chemistry network cooling (primitives are already updated in update_chem)
  if (cooling == COOL_CHEM) call cooling_chem()

  !   boundary contiditions on u
  call boundaryI()

  if ( (cooling == COOL_H).or.(cooling == COOL_NONE) ) then
    !  must update primitives in all the domain
        call calcprim(u,primit)
  else
    !  update primitives on the boundaries
    call calcprim(u,primit,only_ghost=.true.)
  end if

  !  Thermal conduction
  if (th_cond /= 0 ) call thermal_conduction()

end subroutine tstep

end module hydro_solver

!=======================================================================
