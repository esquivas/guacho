!=======================================================================
!> @file main.f90
!> @brief Guacho-3D main program
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

!> @brief Guacho-3D Main Program
!! This is the main program unit of the Guacho-3D code.
!! @n The code itegrates Euler equations in three dimensions, 
!! the choice of the integration method is set in the makefile.
!! @n The flow (conserved) variables are taken to be:
!! @n ieq=
!! @n 1 : rho  (total)
!! @n      2 : rho u
!! @n      3 : rho v
!> @n      4 : rho w
!> @n      5 : Internal energy (thermal+kinetic)
!> @n      6 : bx  (optional, if MHD or PMHD)
!> @n      7 : by  (optional, if MHD or PMHD)
!> @n      8 : bz  (optional, if MHD or PMHD)
!> @n      additional variables advected into the flow, e.g.:
!> @n      9   (6):  n_HI
!> @n      10  (7):  n_HII
!> @n      11  (8):  n_HeI
!> @n      12  (9):  n_HeII
!> @n      13  (10): n_HeIII
!> @n      14  (11): rho*zbar
!> @n      15  (12): ne
!! @n     This can be changed bu the user according to cooling
!! function for instance

program guacho

  use constants, only : day
  use parameters
  use globals
  use init
  use hydro_core, only : calcprim, get_timestep
  use output
  use hydro_solver
  use boundaries

#ifdef RADDIFF
  use difrad
#endif
#ifdef C2ray
  use C2Ray_init
#endif

  implicit none
  integer :: err
  integer :: itprint
  real    :: time, dt, tprint

  !   initializes mpi, and global variables
  call initmain(time, tprint, itprint)

  !   initialize u's
  call initflow(itprint)

  !   impose  boundaries (needed if imposing special BCs)
  call boundaryI(time,0.0)

  !   update primitives with u
  call calcprim(u,primit)

#ifdef RADDIFF
  call diffuse_rad()
#endif
#ifdef C2ray
  call C2ray_initialization()
#endif

#ifdef MPIP
  call mpi_barrier(mpi_comm_world, err)
#endif

  !   time integration
  do while (time.le.tmax)
  
     !   output at intervals tprint
     if(time.ge.tprint) then
        !call diffuse_rad()
        call write_output(itprint)
        if (rank.eq.0) then
           print'(a,i4)', &
           '****************** wrote output *************** &
            & :' , itprint
        end if
        tprint=tprint+dtprint
        itprint=itprint+1
     end if

     !   computes the timestep
     call get_timestep(dt)
     time = time + dt
     if (rank.eq.0) print'(a,es12.3,a,es12.3,a,es12.3,a)',            &
     !'time=',time,'  dt=', dt,' tprint=',tprint
     'time=',time*tsc/day,' day  dt=', dt*tsc/day,' day tprint=',tprint*tsc/day,' day'
     !'time=',time*tsc/myr,' Myr  dt=', dt*tsc/myr,' Myr tprint=',tprint*tsc/myr,' yr'
     
     !   advances the solution
     call tstep(time,dt)
     
  end do

  !   finishes
  if (rank.eq.0) print'(a)',"--- My work here is done, have a nice day ---"

#ifdef MPIP
  call mpi_finalize(err)
#endif
  stop

end program  guacho

!=======================================================================
