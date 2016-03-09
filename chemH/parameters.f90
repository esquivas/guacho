!=======================================================================
!> @file parameters.f90
!> @brief parameters module
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

!> @brief Parameters module
!> @details This module contains parameters of the run, some of this
!! can be moved later to a runtime input file

module parameters
  use constants, only : Rg, amh, au, yr
  implicit none
#ifdef MPIP
  include "mpif.h"
#endif
  !> Path used to write the output
 character (len=128),parameter ::  outputpath='/datos/esquivel/Guacho-1.2/chemcol/'
 !> working directory
 character (len=128),parameter ::  workdir='./'
  
#if defined(MHD) || defined(PMHD)
  integer, parameter :: neqdyn=8        !< num. of eqs  (+scal)
#else
  integer, parameter :: neqdyn=5        !< num. of eqs  (+scal)
#endif
  integer, parameter :: ndim=3          !< num. of dimensions
#ifdef PASSIVES
  integer, parameter :: npas=5          !< num. of passive scalars 
#else
  integer, parameter :: npas=0          !< num. of passive scalars 
#endif
  integer, parameter :: nghost=2        !< num. of ghost cells

  integer, parameter :: nxtot=192       !< Total grid size in X
  integer, parameter :: nytot=192       !< Total grid size in Y
  integer, parameter :: nztot=512       !< Total grid size in Z

#ifdef MPIP
  !   mpi array of processors
  integer, parameter :: MPI_NBX=2        !< number of MPI blocks in X
  integer, parameter :: MPI_NBY=2        !< number of MPI blocks in Y
  integer, parameter :: MPI_NBZ=8        !< number of MPI blocks in Z   
  !> total number of MPI processes
  integer, parameter :: np=MPI_NBX*MPI_NBY*MPI_NBZ
#endif

  !   some parameters for coupling with C2-ray
  !real, parameter :: nphot=5.10e48, ts=1.e5 !Number of Phot. source temp.

  !  set box size   
  real, parameter :: xmax=1.         !< grid extent in X (code units)
  real, parameter :: ymax=1.         !< grid extent in Y (code units)
  real, parameter :: zmax=8./3.      !< grid extent in Z (code units)
  real, parameter :: xphys=1.5E4*au  !< grid extent in X (pohysical units, cgs)

  !  For the equation of state
  real, parameter :: cv=1.5            !< Specific heat at constant volume (/R)
  real, parameter :: gamma=(cv+1.)/cv  !< Cp/Cv
  real, parameter :: mu = 1.           !< mean atomic mass
  
  !  scaling factors to physical (cgs) units
  real, parameter :: T0=10.                 !<  reference temperature (to set cs)
  real, parameter :: rsc=xphys/xmax         !<  distance scaling
  real, parameter :: rhosc=amh*mu           !<  mass density scaling
  real, parameter :: Tempsc=T0*gamma        !<  Temperature scaling
  real, parameter :: vsc2 = gamma*Rg*T0/mu  !<  Velocity scaling
  real, parameter :: vsc = sqrt(vsc2)
  real, parameter :: Psc = rhosc*vsc2       !<  Pressure scaling
  real, parameter :: tsc =rsc/sqrt(vsc2)    !<  time scaling
#ifdef PMHD
!> magnetic fiewld scaling
   real, parameter :: bsc = sqrt(4.0*pi*Psc)
#endif
#ifdef MHD
!> magnetic fiewld scaling
   real, parameter :: bsc = sqrt(4.0*pi*Psc)     
#endif
  !> Maximum integration time
  real, parameter :: tmax    = 2.E3*yr/tsc
  !> interval between consecutive outputs
  real, parameter :: dtprint =  0.1E2*yr/tsc
  real, parameter :: cfl=0.1        !< Courant-Friedrichs-Lewy number
  real, parameter :: eta=0.01       !< artificial viscosity

  !> Warm start flag, if true restarts the code from previous output
  logical, parameter :: iwarm=.true.
  integer            :: itprint0=141  !< number of output to do warm start

  !-------------------------------------------------------------------------
  !  some derived parameters (no need of user's input below this line)

  integer, parameter :: neq=neqdyn + npas  !< number of equations

#ifdef MPIP
  !>  number of physical cells in x in each MPI block
  integer, parameter :: nx=nxtot/MPI_NBX
  !>  number of physical cells in y in each MPI block
  integer, parameter :: ny=nytot/MPI_NBY
  !>  number of physical cells in z in each MPI block
  integer, parameter :: nz=nztot/MPI_NBZ
#else
  integer, parameter :: nx=nxtot, ny=nytot, nz=nztot
  integer, parameter :: np=1, MPI_NBY=1, MPI_NBX=1, MPI_NBZ=1
#endif

  integer, parameter :: nxmin = 1  - nghost   !< lower bound of hydro arrays in x
  integer, parameter :: nxmax = nx + nghost   !< upper bound of hydro arrays in x
  integer, parameter :: nymin = 1  - nghost   !< lower bound of hydro arrays in y
  integer, parameter :: nymax = ny + nghost   !< upper bound of hydro arrays in y
  integer, parameter :: nzmin = 1  - nghost   !< lower bound of hydro arrays in z
  integer, parameter :: nzmax = nz + nghost   !< upper bound of hydro arrays in z
  
  !  more mpi stuff
  integer, parameter ::master=0  !<  rank of master of MPI processes
  
  !   set floating point precision (kind) for MPI messages
#ifdef MPIP
#ifdef DOUBLEP
  integer, parameter :: mpi_real_kind=mpi_real8  !< MPI double precision
#else
  integer, parameter :: mpi_real_kind=mpi_real4  !< MPI single precision
#endif
#endif

end module parameters

!=======================================================================

