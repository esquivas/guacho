!=======================================================================
!> @file init.f90
!> @brief Guacho-3D initialization module
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

!> @brief Guacho-3D initialization
!> @details This module contains the routines needed to initializa the
!! code, it also initiaizes all the modules set by the user.

module init

contains

!> @brief Main initialization routine
!> @details This subsroutine initializes all the variables in the globals
!! module, MPI, cooling and user_mod routines;
!! and outputs to screen the main parameters used in the run
!> @param real [out] time   : time (code units)
!> @param real [out] tprint : time of next output
!> @param integer [out] itprint : number of next output

subroutine initmain(time, tprint, itprint)

  use constants, only : yr
  use parameters
  use globals
#ifdef COOLINGDMC
  use cooling_dmc, only : read_table, cooltab
#endif
#ifdef COOLINGCHI
  use cooling_chi, only : read_table, cooltab
#endif
#ifdef RADDIFF
  use difrad
#endif
#ifdef THERMAL_COND
  use thermal_cond
#endif
  use user_mod

  implicit none
  real,    intent(out) :: time
  real,    intent(out) ::tprint
  integer, intent(out) :: itprint

#ifdef MPIP
  integer :: err, nps
  integer, dimension(0:ndim-1) :: dims
  logical, dimension(0:ndim-1) :: period
#endif  
  !initializes MPI

#ifdef MPIP
#ifdef PERIODX
  logical, parameter :: perx=.true.
#else
  logical, parameter :: perx=.false.
#endif
#ifdef PERIODY
  logical, parameter :: pery=.true.
#else
  logical, parameter :: pery=.false.
#endif
#ifdef PERIODZ
  logical, parameter :: perz=.true.
#else
  logical, parameter :: perz=.false.
#endif
  period(0)=perx
  period(1)=pery
  period(2)=perz
  dims(0)  =mpicol
  dims(1)  =mpirow
  dims(2)  =mpirowz
  !
  call mpi_init (err)
  call mpi_comm_rank (mpi_comm_world,rank,err)
  call mpi_comm_size (mpi_comm_world,nps,err)
  if (nps.ne.np) then
     print*, 'processor number (',nps,') is not equal to pre-defined number (',np,')'
     call mpi_finalize(err) 
     stop
  endif
#else
  rank=0
  coords(:)=0
#endif
  if(rank.eq.master) then
     print '(a)' ,"*******************************************"
     print '(a)' ,"                        _                 *"
     print '(a)' ,"  __   _   _  __ _  ___| |__   ___    3   *"
     print '(a)' ," / _ `| | | |/ _` |/ __| '_ \ / _ \    D  *"
     print '(a)' ,"| (_| | |_| | (_| | (__| | | | (_) |      *"
     print '(a)' ," \__, |\__,_|\__,_|\___|_| |_|\___/       *"
     print '(a)' ," |___/                                    *"
  endif
#ifdef MPIP
  if(rank.eq.master) then
     print '(a,i3,a)','*    running with mpi in', np , ' processors    *'
     print '(a)' ,'*******************************************'
  end if
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, 1            &
       , comm3d, err)
  call mpi_comm_rank(comm3d, rank, err)
  call mpi_cart_coords(comm3d, rank, ndim, coords, err)
  print '(a,i3,a,3i4)', 'processor ', rank                              &
       ,' ready w/coords',coords(0),coords(1),coords(2)   
  call mpi_cart_shift(comm3d, 0, 1, left  , right, err)
  call mpi_cart_shift(comm3d, 1, 1, bottom, top  , err)
  call mpi_cart_shift(comm3d, 2, 1, out   , in   , err)
  call mpi_barrier(mpi_comm_world, err)   
  !
#else
  print '(a)' ,'*******************************************'
  print '(a)' ,'*     running on a single processor       *'
  print '(a)' ,'*******************************************'
#endif
  
  !   grid spacing
  dx=xmax/nxtot
  dy=ymax/nytot
  dz=zmax/nztot
  
  !   initialize time integration
  if (.not.iwarm) then
     if(rank.eq.master) then
        print'(a)', 'Starting cold'
        print'(a)',' ' 
     endif
     itprint=0
     time=0.
     tprint=0.
  else
     itprint=itprint0
     time=real(itprint)*dtprint
     if(rank.eq.master) then
        print'(a,i,a,es12.3,a)', 'Warm start , from output ',itprint,' at a tim!e ',time*tsc/yr,' yr'
        print'(a)',' ' 
     end if
     tprint=time+dtprint
  end if
  
  !   allocate big arrays in memory
  allocate (     u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (    up(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (primit(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (     f(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (     g(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  allocate (     h(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
  
  !   DMC cooling
#ifdef COOLINGDMC
  call read_table()
#endif

!   CHIANTI COOLING
#ifdef COOLINGCHI
  call read_table()
#endif

!   BBC COOLING
!#ifdef COOLINGBBC
!  do ii=0,(nps-1)
!     if (rank.eq.ii) then
!        call bbcrd
!        print'(a,i4,a)','rank:',rank,' Just read the tables'
!     endif
!#ifdef MPIP
!     call mpi_barrier (mpi_comm_world, err)
!#endif
!  end do
!#endif

#ifdef THERMAL_COND
  call init_thermal_cond()
#endif
#ifdef RADDIFF
  call init_rand()
#endif

  !  User input initialization, it is called always, 
  !  it has to be there, even empty
  call init_user_mod()
  
  !   write report of compilation parameters

#ifdef MPIP
  call mpi_barrier(mpi_comm_world, err)
  if(rank.eq.master) then
    print'(a)',''
#endif
    print'(a,i0,a)', 'Running with ',neq,' total equations' 
    print'(a,i0,a,i0,a,i0)','Resolution is (nxtot, nytot, nztot) ', nxtot,' ',nytot,' ',nztot
    print'(a)',''
#ifdef MHD
     print'(a)', 'Full MHD enabled'
     print'(a)', ''
#endif
#ifdef PMHD
     print'(a)', 'Passive MHD enabled'
     print'(a)', ''
#endif
#ifdef DOUBLEP
     print'(a)', 'Double precision used (reals are 8 bytes long)'
     print'(a)', ''
#else
     print'(a)', 'Single precision used (reals are 4 bytes long)'
     print'(a)', ''
#endif
#ifdef HLLC
     print'(a)', 'The Riemann solver is HLLC'
     print'(a)', ''
#endif
#ifdef HLL
     print'(a)', 'The Riemann solver is HLL'
     print'(a)', ''
#endif
#ifdef HLL_HLLC
     print'(a)', 'The Riemann solver is HLL-HLLC (hybrid)'
     print'(a)', ''
#endif
#ifdef HLLE
     print'(a)', 'The Riemann solver is HLLE)'
     print'(a)', ''
#endif
#ifdef HLLD
     print'(a)', 'The Riemann solver is HLLD'
     print'(a)', ''
#endif
#ifdef EOS_ADIABATIC
     print'(a)', 'The code uses an AIABATIC EOS'
     print'(a)', ''
#endif
#ifdef EOS_SINGLE_SPECIE
     print'(a)', 'The EOS considers only one specie of H'
     print'(a)', ''
#endif
#ifdef EOS_H_RATE
     print'(a)', 'The EOS considers a rate equation for H'
     print'(a)', ''
#endif
#ifdef EOS_MULTI_SPECIES
     print'(a)', 'The EOS considers multiple species'
     print'(a)', ''
#endif
#ifdef NO_COOL
     print'(a)', 'Cooling is turned off'
     print'(a)', ''
#endif
#ifdef COOLINGH
     print'(a)', 'Radiative cooling ON (w/parametrized cooling curve)'
     print'(a)', ''
#endif
#ifdef COOLINGBBG
     print'(a)', 'Radiative cooling ON (w/ Benjamin Benson & Cox 2003 prescription)'
     print'(a)', ''
#endif
#ifdef COOLINGDMC
     print'(a)', 'Radiative cooling ON (w/ Dalgarno & Mc Cray, coronal eq.)'
     print'(a)', ''
#endif
#ifdef COOLINGCHI
     print'(a)', 'Radiative cooling ON (Uses table from CHIANTI)'
     print'(a)', ''
#endif
#ifdef RADDIFF
     print'(a)','Diffuse radiative transfer enabled, local'
#endif
     print'(a)', '-----  OUTPUT -----------------------'
     print'(a)', 'path: '//trim(outputpath)
     print'(a)', 'in the following format(s):'
#ifdef OUTBIN
     print'(a)', '*.bin (raw unformatted)'
#endif
#ifdef OUTDAT
     print'(a)', '*.dat (formatted, beware of big files)'
#endif
#ifdef OUTVTK
     print'(a)', '*.vtk (binary VTK)'
#endif
     print'(a)', ''
     print'(a)', '----- BOUNDARY CONDITIONS -----------'
#ifdef PERIODX
     print'(a)', 'LEFT & RIGHT: PERIODIC'
#endif
#ifdef PERIODY
     print'(a)', 'BOTTOM & TOP: PERIODIC'
#endif
#ifdef PERIODZ
     print'(a)', 'IN & OUT: PERIODIC'
#endif
#ifdef REFXL
     print'(a)', 'LEFT:   REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFXL
     print'(a)', 'LEFT:   OUTFLOW    (OPEN)'
#endif
#ifdef INFXL
     print'(a)', 'LEFT:   INFLOW     (user defined)'
#endif
#ifdef REFXR
     print'(a)', 'RIGHT:  REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFXR
     print'(a)', 'RIGHT:  OUTFLOW    (OPEN)'
#endif
#ifdef INFXR
     print'(a)', 'RIGHT:  INFLOW     (user defined)'
#endif
#ifdef REFYB
     print'(a)', 'BOTTOM: REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFYB
     print'(a)', 'BOTTOM: OUTFLOW    (OPEN)'
#endif
#ifdef INFYB
     print'(a)', 'BOTTOM: INFLOW     (user defined)'
#endif
#ifdef REFYT
     print'(a)', 'TOP:    REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFYT
     print'(a)', 'TOP:    OUTFLOW    (OPEN)'
#endif
#ifdef INFYT
     print'(a)', 'TOP:    INFLOW     (user defined)'
#endif
#ifdef REFZO
     print'(a)', 'OUT:    REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFZO
     print'(a)', 'OUT:    OUTFLOW    (OPEN)'
#endif
#ifdef INFZO
     print'(a)', 'OUT:    INFLOW     (user defined)'
#endif
#ifdef REFZI
     print'(a)', 'IN: REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFZI
     print'(a)', 'IN:     OUTFLOW    (OPEN)'
#endif
#ifdef INFZI
     print'(a)', 'IN: INFLOW     (user defined)'
#endif
    print'(a)', ''
     print'(a)', '----- OTHER STUFF -----------'
#ifdef RADDIFF
     print'(a)', 'Diffuse radiation (local+MPI) enabled'
#endif
#ifdef THERMAL_COND
     print'(a)', 'Thermal conduction included (isotropic)'
#endif
#ifdef OTHERB
     print'(a)', 'Other boundaries enabled (otherbounds.f90)'
#endif
     print'(a)', ''
#if LIMITER==-1
     print'(a)', 'No average in the limiter (reduces to 1st order)'
#endif
#if LIMITER==0
     print'(a)', 'No limiter'
#endif
#if LIMITER==1
     print'(a)', 'MINMOD limiter -most diffusive-'
#endif
#if LIMITER==2
     print'(a)', 'Falle Limiter (Van Leer)'
#endif
#if LIMITER==3
     print'(a)', 'Van Albada Limiter'
#endif
#if LIMITER==4
     print'(a)', 'UMIST limiter -least diffusive-'
#endif
#if LIMITER==5
     print'(a)', 'Woodward Limiter (MC-limiter; monotonized central difference)'
#endif
#if LIMITER==6
     print'(a)', 'SUPERBEE limiter (tends to flatten circular waves)'
#endif
     print'(a)', ''
     print'(a)','***********************************************'
#ifdef MPIP
  end if
  call mpi_barrier(mpi_comm_world, err)
#endif

end subroutine initmain

!====================================================================

!> @brief Initializes the conserved variables, in the globals module
!> @details Initializes the conserved variables, in the globals module
!> @param real [inout] itprint : number of current output

subroutine initflow(itprint)

  use parameters, only : outputpath, iwarm, itprint
  use globals, only : u, rank
  use user_mod, only : initial_conditions
  implicit none

#ifdef MPIP
  include "mpif.h"
#endif

  integer , intent(inout) :: itprint
  integer ::  unitin,err
  character (len=128) file1

  if (.not.iwarm) then

    call initial_conditions(u,0.)

  else

     !   read from previous (.bin) output
#ifdef MPIP
     write(file1,'(a,i3.3,a,i3.3,a)')  &
          trim(outputpath)//'BIN/points',rank,'.',itprint,'.bin'
     unitin=rank+10
#else
     write(file1,'(a,i3.3,a)')         &
          trim(outputpath)//'BIN/points',itprint,'.bin'
     unitin=10
#endif
     open(unit=unitin,file=file1,status='unknown',form='unformatted', &
          convert='LITTLE_ENDIAN')
     
     read(unitin) u(:,:,:,:)
     close(unitin)
     
     print'(i3,a,a)',rank,'read',trim(file1)
     itprint=itprint+1
     
     
#ifdef MPIP
     call mpi_barrier(mpi_comm_world,err)
#endif

  end if

end subroutine initflow

!====================================================================

end module init

!====================================================================





