!=======================================================================
!> @file cooling_schure.f90
!> @brief Cooling module with cooling curves of Schure et al. 2009
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

!> @brief Cooling module with cooling curves of Schure et al. 2009
!> @details Cooling module cooling curves of Schure et al. 2009, A&A, 508,751
!> @n The location of the tables is assumed to be in
!> src/cool_lib/coolingSKKKV.tab

module cooling_schure

  implicit none
  real (kind=8), allocatable :: cooltab(:,:)

  !> The following selects the column of the ionization fraction used for low T
  !> dmc_f = 1 : f_i = 1e-4
  !>         2 : f_i = 1e-3
  !>         3 : f_i = 1e-2
  !>         4 : f_i = 1e-1
  integer, parameter         :: dmc_f = 4


contains

  !=======================================================================
  !> @brief Initializes the DMC cooling
  !> @details Declares variables and reads table
  subroutine init_cooling_schure()

    implicit none

    allocate(cooltab(2,180))
    call read_table_schure()

  end subroutine init_cooling_schure

  !=======================================================================
  !> @brief Reads the cooling curve table
  !> @details Reads the cooling curve table,
  !! the location is assumed in /src/cool_lib/coolingSKKKV.tab
  subroutine read_table_schure()

    use parameters, only : workdir, master
    use globals, only : rank
#ifdef MPIP
    use mpi
#endif
    implicit none
    integer :: i, err
    real (kind=8) :: data(5)  ! Table has five columns

    if(rank == master) then
      open(unit=10,file= trim(workdir)//'../src/cool_lib/coolingSKKKV.tab',&
           status='old')
      do i=1,180
        read(10,*) data(:)
        cooltab(1,i)= 10.0**data(1)
        cooltab(2,i)= 10.0**(- data(1+dmc_f) )
      end do
      close(unit=10)
    endif
#ifdef MPIP
    call mpi_bcast(cooltab,360,mpi_double_precision,0,mpi_comm_world,err)
#endif

end subroutine read_table_schure

  !=======================================================================
  !> @brief Returns the cooling coefficient interpolating the table
  !> @param real [in] T : Temperature K
  function get_lambda(T)

    implicit none
    real , intent(in) :: T
    integer           :: if1
    real, parameter   :: deltaTemp=0.04 !  spacing of T in tables
    real, parameter   :: logTmin = 1.0 !  log of minimum value of temp
    real (kind=8)     :: get_lambda, T0, T1, C0, C1


    if(T.gt.1e8) then
      get_lambda=0.21e-26*sqrt(T)
    else
      if1=int(( log10(T)- logTmin) /deltaTemp) + 1
      if (if1 < 1) then
        get_lambda = cooltab(2,1)
        return
      end if
      T0=cooltab(1,if1)
      c0=cooltab(2,if1)
      T1=cooltab(1,if1+1)
      c1=cooltab(2,if1+1)
      get_lambda=(c1-c0)*(T-T0)/(T1-T0)+c0
    end if

  end function get_lambda

  !=======================================================================
  !> @brief High level wrapper to apply cooling with CHIANTI tables
  !> @details High level wrapper to apply cooling with CHIANTI tables
  !> @n cooling is applied in the entire domain and updates both the
  !! conserved and primitive variables
  subroutine coolingschure()

    use parameters, only : nx, ny, nz, cv, Psc, tsc, mhd, n1_chem
    use constants,  only : Kb
    use globals, only : u, primit, dt_CFL
    use hydro_core, only : u2prim
    use network
    implicit none
    real                 :: T , dens
    real, parameter      :: Tmin=10.
    real (kind=8)        :: Lambda0, emtauC, gain
    integer              :: i, j, k
    real                 :: dt_seconds, ch_factor

    dt_seconds = dt_CFL*tsc

    do k=1,nz
      do j=1,ny
        do i=1,nx

          !   get the primitives (and T)
          call u2prim(u(:,i,j,k),primit(:,i,j,k),T)

          if(T > Tmin) then

            Lambda0=get_lambda(T)

            !   here we should add support for ionization heating
            gain=0.0

            dens=primit(1,i,j,k)
            !  e^{-dt/Tau}=e^{-2. L0 dt/(3 n K T)}
            emtauC = exp( -2.0*dt_seconds*dens*Lambda0/(3.0*Kb*T) )
            !  this is the Temperature factor of change
            ch_factor = (gain/(dens**2*Lambda0))*(1.0-emtauC)+emtauC

            !  limit changes to avoid catastrophic cooling
            ch_factor = max(ch_factor, 0.1)
            ch_factor = min(ch_factor,10.0)

            !  apply cooling to primitive and conserved variables
            primit(5,i,j,k)=primit(5,i,j,k)*ch_factor
            !  update total energy density

            u(5,i,j,k) = cv*primit(5,i,j,k)                                    &
                         + 0.5*primit(1,i,j,k)*(  primit(2,i,j,k)**2           &
                                                + primit(3,i,j,k)**2           &
                                                + primit(4,i,j,k)**2  )
#ifdef BFIELD
          if (mhd) then
            u(5,i,j,k) = u(5,i,j,k) + 0.5*(  primit(6,i,j,k)**2                &
                                           + primit(7,i,j,k)**2                &
                                           + primit(8,i,j,k)**2  )
          end if
#endif

          end if
        end do
      end do
    end do


  end subroutine coolingschure

  !======================================================================

end module cooling_schure
