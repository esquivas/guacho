!=======================================================================
!> @file init.f90
!> @brief Guacho-3D initialization module
!> @author Alejandro Esquivel
!> @date 4/May/2016
! Copyright (c) 2020 Guacho Co-Op
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
!> code, it also initiaizes all the modules set by the user.

module init

contains

  !> @brief Main initialization routine
  !> @details This subsroutine initializes all the variables in the globals
  !> module, MPI, cooling and user_mod routines;
  !> and outputs to screen the main parameters used in the run
  !> @param real [out] tprint : time of next output
  !> @param integer [out] itprint : number of next output
  subroutine initmain(tprint, itprint)

    use constants
    use parameters
    use globals
    use lmp_module
    use cooling_dmc
    use cooling_chi
    use cooling_schure
    use difrad
    use thermal_cond
    use flux_cd_module
    use user_mod
    implicit none
    real,    intent(out) ::tprint
    integer, intent(out) :: itprint

#ifdef MPIP
    integer :: err, nps
    integer, dimension(0:ndim-1) :: dims
    logical, dimension(0:ndim-1) :: period
    logical :: perx=.false., pery=.false., perz=.false.
#endif

    !initializes MPI
#ifdef MPIP

    if (bc_left   == BC_PERIODIC .and. bc_right == BC_PERIODIC) perx=.true.
    if (bc_bottom == BC_PERIODIC .and. bc_top   == BC_PERIODIC) pery=.true.
    if (bc_out    == BC_PERIODIC .and. bc_in    == BC_PERIODIC) perz=.true.

    period(0)=perx
    period(1)=pery
    period(2)=perz
    dims(0)  =MPI_NBX
    dims(1)  =MPI_NBY
    dims(2)  =MPI_NBZ

    call mpi_init (err)
    call mpi_comm_rank (mpi_comm_world,rank,err)
    call mpi_comm_size (mpi_comm_world,nps,err)
    if (nps.ne.np) then
      print*,                                                                  &
        'processor number (',nps,') is not equal to pre-defined number (',np,')'
      call mpi_finalize(err)
      stop
    endif
#else
    rank=0
    coords(:)=0
#endif
    if(rank.eq.master) then
      print '(a)' ,"*********************************************"
      print '(a)' ,"*                         _                 *"
      print '(a)' ,"*   __   _   _  __ _  ___| |__   ___    3   *"
      print '(a)' ,"*  / _ `| | | |/ _` |/ __| '_ \ / _ \    D  *"
      print '(a)' ,"* | (_| | |_| | (_| | (__| | | | (_) |      *"
      print '(a)' ,"*  \__, |\__,_|\__,_|\___|_| |_|\___/       *"
      print '(a)' ,"*  |___/                                    *"
    endif
#ifdef MPIP
    if(rank.eq.master) then
      print '(a,i3,a)','*    running with mpi in ', np, ' processors     *'
      print '(a)', '*********************************************'
    end if
    call mpi_cart_create(mpi_comm_world, ndim, dims, period,.true., comm3d, err)
    call mpi_comm_rank(comm3d, rank, err)
    call mpi_cart_coords(comm3d, rank, ndim, coords, err)
    print '(a,i3,a,3i4)', 'processor ', rank                                   &
    ,' ready w/coords',coords(0),coords(1),coords(2)
    call mpi_cart_shift(comm3d, 0, 1, left  , right, err)
    call mpi_cart_shift(comm3d, 1, 1, bottom, top  , err)
    call mpi_cart_shift(comm3d, 2, 1, out   , in   , err)
    call mpi_barrier(mpi_comm_world, err)
#else
    print '(a)' ,'*********************************************'
    print '(a)' ,'*      running on a single processor        *'
    print '(a)' ,'*********************************************'
#endif

    !   grid spacing
    dx=xmax/nxtot
    dy=ymax/nytot
    dz=zmax/nztot

    !   initialize time integration
    currentIteration = 1
    if (.not.iwarm) then
      if(rank.eq.master) then
        print'(a)', 'Starting cold'
        print'(a)',' '
      endif
      itprint=0
      time= time_0
      tprint=dtprint
    else
      itprint=itprint0
      time=time_0
      if(rank == master) then
        print'(a,i0,a,es12.3,a)', 'Warm start , from output ',itprint,         &
                                  ' at a time ',time*tsc/yr,' yr'
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
    allocate (Temp(nxmin:nxmax,nymin:nymax,nzmin:nzmax) )

    if (riemann_solver == SOLVER_HLLE_SPLIT_ALL )                              &
    allocate (primit0(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax))

#ifdef BFIELD
    if (enable_flux_cd)                                                        &
    allocate ( e(3,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )
#endif

    if (enable_lmp) call init_lmp()

    !   DMC cooling
    if (cooling == COOL_DMC) call init_cooling_dmc()

    !   CHIANTI COOLING
    if (cooling == COOL_CHI) call init_cooling_chianti()

    !   SKKKV cooling
    if (cooling == COOL_SKKKV) call init_cooling_schure()

    !  Thermal conduction
    if (th_cond /= TC_OFF) call init_thermal_cond()

    !  create directories to write the outputs
    if (rank == master) then
      if (out_bin) then
        call system('if [ ! -e '//trim(outputpath)//'BIN ]; then mkdir -p '    &
                                //trim(outputpath)//'BIN ; fi')
      end if
      if (out_vtk) then
        call system('if [ ! -e '//trim(outputpath)//'VTK ]; then mkdir -p '    &
                                //trim(outputpath)//'VTK ; fi')
      end if
      if (out_silo) then
        call system('if [! -e'//trim(outputpath)//'SILO/BLOCKS];then mkdir -p '&
                              //trim(outputpath)//'SILO/BLOCKS ; fi')
      end if
    end if

    !  User input initialization, it is called always,
    !  it has to be there, even if empty
    call init_user_mod()

    !  Diffuse radiation transfer module required random numbers
    if (dif_rad) call init_difrad()

    !   write report of compilation parameters
#ifdef MPIP
    call mpi_barrier(mpi_comm_world, err)
    if(rank.eq.master) then
      print'(a)',''
#endif
      print'(a,i0,a)', 'Running with ',neq,' total equations'
      print'(a,i0,a,i0,a)', '(',neqdyn,' dynamical, and ',npas,' passives)'
      print'(a,i0,a,i0,a,i0)','Resolution is (nxtot, nytot, nztot) ',          &
            nxtot,' ', nytot,' ',nztot
      print'(a)',''

      if (mhd) then
        print'(a)', 'Full MHD enabled'
        print'(a)', ''
      end if
      if (pmhd) then
        print'(a)', 'Passive MHD enabled'
        print'(a)', ''
      end if
      if (mhd .and. pmhd) then
        print'(a)', "Error, select only one of the options, 'mhd' or 'pmhd'"
        print'(a)', ''
        stop
      end if

#ifdef DOUBLEP
      print'(a)', 'Double precision used (reals are 8 bytes long)'
      print'(a)', ''
#else
      print'(a)', 'Single precision used (reals are 4 bytes long)'
      print'(a)', ''
#endif

      if (riemann_solver == SOLVER_HLL) then
        print'(a)', 'The Riemann solver is HLL'
        print'(a)', ''
      else if (riemann_solver == SOLVER_HLLC) then
        print'(a)', 'The Riemann solver is HLLC'
        print'(a)', ''
      else if (riemann_solver == SOLVER_HLLE) then
        print'(a)', 'The Riemann solver is HLLE'
        print'(a)', ''
      else if (riemann_solver == SOLVER_HLLD) then
        print'(a)', 'The Riemann solver is HLLD'
        print'(a)', ''
      else if (riemann_solver == SOLVER_HLLE_SPLIT_B) then
        print'(a)', 'The Riemann solver is HLLE with split B field'
        print'(a)', ''
      else if (riemann_solver == SOLVER_HLLD_SPLIT_B) then
        print'(a)', 'The Riemann solver is HLLD with split B field'
        print'(a)', ''
      else if (riemann_solver == SOLVER_HLLE_SPLIT_ALL) then
        print'(a)', 'The Riemann solver is HLLE with split in All Variables'
        print'(a)', ''
      else if (riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
        print'(a)', 'The Riemann solver is HLLD with split in All Variables'
        print'(a)', ''
      else if (riemann_solver == SOLVER_RHLL) then
        print'(a)', 'The Riemann solver is the relativistic HLL'
        print'(a)', ''
      else if (riemann_solver == SOLVER_RHLLC) then
        print'(a)', 'The Riemann solver is the relativistic HLLC'
        print'(a)', ''
      else
        print'(a)', 'Unrecognized Riemann Solver'
        print'(a)', ''
        stop
      end if

      if (enable_flux_cd) then
        print'(a)', 'div(B) constrained with flux-CD method'
        print'(a)', ''
      end if

      if (eight_wave) then
        print'(a)', 'div(B) constrained with 8 wave method'
        print'(a)', ''
      end if

      if (eq_of_state == EOS_ADIABATIC) then
        print'(a)', 'The code uses an AIABATIC EOS'
        print'(a)', ''
      else if (eq_of_state == EOS_SINGLE_SPECIE) then
        print'(a)', 'The EOS considers only one specie of H'
        print'(a)', ''
      else if (eq_of_state == EOS_H_RATE) then
        print'(a)', 'The EOS considers a rate equation for H'
        print'(a)', ''
      else if (eq_of_state == EOS_CHEM) then
        print'(a)', 'The EOS considers multiple species (chemical network)'
        print'(a)', ''
      else if (eq_of_state == EOS_REL_IDEAL) then
        print'(a)', 'Relativistic EOS with constant Gamma'
        print'(a)', ''
      else if (eq_of_state == EOS_REL_TM) then
        print'(a)', 'Relativistic TM EOS (Mignone et al 2005)'
        print'(a)', ''

      else
        print'(a)', 'Unrecognized equation of state'
        print'(a)', ''
        stop
      end if

      if (cooling == COOL_NONE) then
        print'(a)', 'Cooling is turned off'
        print'(a)', ''
      else if (cooling == COOL_H) then
        print'(a)', 'Radiative cooling ON (w/parametrized cooling curve)'
        print'(a)', ''
      else if (cooling == COOL_BBC) then
        print'(a)', 'Rad. cooling ON (Benjamin Benson & Cox 2003 prescription)'
        print'(a)', ''
      else if (cooling == COOL_DMC) then
        print'(a)', 'Radiative cooling ON (w/ Dalgarno & Mc Cray, coronal eq.)'
        print'(a)', ''
      else if (cooling == COOL_CHI) then
        print'(a)', 'Radiative cooling ON (Uses table from CHIANTI)'
        print'(a)', ''
      else if (cooling == COOL_SKKKV) then
        print'(a)', 'Radiative cooling ON (Uses tables from Schure et al. 2009)'
        print'(a)', ''
      else
        print'(a)', 'Unrecognized cooling scheme'
        print'(a)', ''
        stop
      end if

      if (dif_rad) print'(a)','Diffuse radiative transfer enabled, local'

      print'(a)', '-----  OUTPUT -----------------------'
      print'(a)', 'path: '//trim(outputpath)
      print'(a)', 'in the following format(s):'
      if (out_bin) print'(a)', '*.bin (binary, with a small header)'
      if (out_vtk) print'(a)', '*.vtk (binary VTK)'
      print'(a)', ''

      print'(a)', '----- BOUNDARY CONDITIONS -----------'
      if (bc_left == BC_PERIODIC .and. bc_right == BC_PERIODIC) then
        print'(a)', 'LEFT & RIGHT: PERIODIC'
      else if (bc_left == BC_PERIODIC .and. bc_right /= bc_left) then
        print'(a)', 'Invalid periodic BCs'
        stop
      end if
      if (bc_bottom == BC_PERIODIC .and. bc_top == BC_PERIODIC) then
        print'(a)', 'BOTTOM & TOP: PERIODIC'
      else if (bc_bottom == BC_PERIODIC .and. bc_top /= bc_bottom) then
        print'(a)', 'Invalid periodic BCs'
        stop
      end if
      if (bc_out == BC_PERIODIC .and. bc_in == BC_PERIODIC) then
        print'(a)', 'IN & OUT: PERIODIC'
      else if (bc_out == BC_PERIODIC .and. bc_in /= bc_out) then
        print'(a)', 'Invalid periodic BCs'
        stop
      end if
      if (bc_left == BC_OUTFLOW  ) print'(a)', 'LEFT:   OUTFLOW    (OPEN)'
      if (bc_left == BC_CLOSED   ) print'(a)', 'LEFT:   REFLECTIVE (CLOSED)'
      if (bc_left == BC_OTHER    ) print'(a)', 'LEFT:   OTHER      (user set)'
      if (bc_right == BC_OUTFLOW ) print'(a)', 'RIGHT:  OUTFLOW    (OPEN)'
      if (bc_right == BC_CLOSED  ) print'(a)', 'RIGHT:  REFLECTIVE (CLOSED)'
      if (bc_right == BC_OTHER   ) print'(a)', 'RIGHT:  OTHER      (user set)'
      if (bc_bottom == BC_OUTFLOW) print'(a)', 'BOTTOM: OUTFLOW    (OPEN)'
      if (bc_bottom == BC_CLOSED ) print'(a)', 'BOTTOM: REFLECTIVE (CLOSED)'
      if (bc_bottom == BC_OTHER  ) print'(a)', 'BOTTOM: OTHER      (user set)'
      if (bc_top == BC_OUTFLOW   ) print'(a)', 'TOP:    OUTFLOW    (OPEN)'
      if (bc_top == BC_CLOSED    ) print'(a)', 'TOP:    REFLECTIVE (CLOSED)'
      if (bc_top == BC_OTHER     ) print'(a)', 'TOP:    OTHER      (user set)'
      if (bc_out == BC_OUTFLOW   ) print'(a)', 'OUT:    OUTFLOW    (OPEN)'
      if (bc_out == BC_CLOSED    ) print'(a)', 'OUT:    REFLECTIVE (CLOSED)'
      if (bc_out == BC_OTHER     ) print'(a)', 'OUT:    OTHER      (user set)'
      if (bc_in == BC_OUTFLOW    ) print'(a)', 'IN:     OUTFLOW    (OPEN)'
      if (bc_in == BC_CLOSED     ) print'(a)', 'IN:     REFLECTIVE (CLOSED)'
      if (bc_in == BC_OTHER      ) print'(a)', 'IN:     OTHER      (user set)'
      if (bc_user)  print'(a)', 'Other boundaries enabled (user_mod.f90)'
      print'(a)', ''

      if(enable_lmp) then
        print'(a)', 'Lagrangian tracer particles module enabled'
        print'(a)', ''
        if (lmp_distf) then
          print'(a)', 'Evolution of LMPs SED enabled'
        print'(a)', ''
        end if
      end if

      print'(a)', '----- OTHER STUFF -----------'
      if (dif_rad) print'(a)', 'Diffuse radiation (local+MPI) enabled'
      print'(a)', ''
      if (th_cond == TC_ISOTROPIC) then
        print'(a)', 'Thermal conduction enabled (isotropic)'
      else if (th_cond == TC_ANISOTROPIC) then
        print'(a)', 'Thermal conduction enabled (Anisotropic)'
      endif
      print'(a)', ''
      if (slope_limiter == LIMITER_NO_AVERAGE) then
        print'(a)', 'No average in the limiter (reduces to 1st order)'
      else if (slope_limiter == LIMITER_NO_LIMIT) then
        print'(a)', 'No limiter'
      else if (slope_limiter == LIMITER_MINMOD) then
        print'(a)', 'MINMOD limiter -most diffusive-'
      else if (slope_limiter == LIMITER_VAN_LEER) then
        print'(a)', 'Falle Limiter (Van Leer)'
      else if (slope_limiter == LIMITER_VAN_ALBADA) then
        print'(a)', 'Van Albada Limiter'
      else if (slope_limiter == LIMITER_UMIST) then
        print'(a)', 'UMIST limiter -least diffusive-'
      else if (slope_limiter == LIMITER_WOODWARD) then
        print'(a)', 'Woodward Limiter (MC-limiter)'
      else if (slope_limiter == LIMITER_SUPERBEE) then
        print'(a)', 'SUPERBEE limiter (tends to flatten circular waves)'
      else
        print'(a)', 'Unrecognized limiter'
        stop
      end if
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

    use parameters, only : outputpath, iwarm, enable_lmp, lmp_distf !, itprint0
    use globals, only : u, rank, Q_MP0, MP_SED, partID, P_DSA, partID,         &
                        partOwner, n_activeMP
    use user_mod, only : initial_conditions
#ifdef MPIP
    use mpi
#endif
    implicit none

    integer , intent(inout) :: itprint
    integer ::  unitin,err
    character (len=128) :: file1
    character           :: byte_read
    character, parameter  :: lf = char(10)
    integer :: nxp, nyp, nzp, x0p, y0p, z0p, mpi_xp, mpi_yp, mpi_zp,neqp,      &
    neqdynp, nghostp
    real :: dxp, dyp, dzp, scal(3), cvp
    integer :: npp, n_mpp, NBinsSEDMPP, i_mp

    if (.not.iwarm) then

      call initial_conditions(u)

    else

      !   read from previous (.bin) output
#ifdef MPIP
    write(file1,'(a,i3.3,a,i3.3,a)')                                           &
          trim(outputpath)//'BIN/points',rank,'.',itprint,'.bin'
    unitin=10  !*rank  (not needed)
#else
    write(file1,'(a,i3.3,a)')                                                  &
          trim(outputpath)//'BIN/points',itprint,'.bin'
    unitin=10
#endif
      open(unit=unitin,file=file1,status='old', access='stream' )

      !   discard the ascii header
      do while (byte_read /= achar(255) )
        read(unitin) byte_read
        !print*, byte_read
      end do
      !  read bin header, sanity check to do
      read(unitin) byte_read
      read(unitin) byte_read
      read(unitin) nxp, nyp, nzp
      read(unitin) dxp, dyp, dzp
      read(unitin) x0p, y0p, z0p
      read(unitin) mpi_xp, mpi_yp, mpi_zp
      read(unitin) neqp, neqdynp
      read(unitin) nghostp
      read(unitin) scal(1:3)
      read(unitin) cvp
      read(unitin) u(:,:,:,:)
      close(unitin)

      print'(i3,a,a)',rank,' read: ',trim(file1)

      if(enable_lmp) then
        !   Read Lagrangian Particles info if they are enabled
#ifdef MPIP
        write(file1,'(a,i3.3,a,i3.3,a)')                                       &
              trim(outputpath)//'BIN/lmp',rank,'.',itprint,'.bin'
#else
        write(file1,'(a,i3.3,a)') trim(outputpath)//'BIN/lmp',itprint,'.bin'
#endif

        unitin=10
        open(unit=unitin,file=file1,status='unknown',access='stream')

        read(unitin) npp, n_mpp, n_activeMP, NBinsSEDMPP

        do i_mp=1,n_activeMP

          partOwner(i_mp) = rank
          read(unitin) partID(i_mp)
          read(unitin) Q_MP0(i_mp,1:3)

          if(lmp_distf) then
            read(unitin) Q_MP0(i_mp,11:12)
            read(unitin) MP_SED(1,:,i_mp)
            read(unitin) MP_SED(2,:,i_mp)
            read(unitin) P_DSA(i_mp,:,:)
          end if

        end do

        close(unitin)

        print'(i3,a,a,a,i0,a)',rank,' read: ',trim(file1), ' (',n_activeMP,             &
                              ' active particles)'


      end if

      itprint=itprint+1

#ifdef MPIP
      call mpi_barrier(mpi_comm_world,err)
#endif

    end if

  end subroutine initflow

  !====================================================================

end module init

!====================================================================
