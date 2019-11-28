!=======================================================================
!> @file cooling_schure.f90
!> @brief Cooling module with CHIANTI generated cooling curves
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

!> @brief Cooling module with CHIANTI generated cooling curves
!> @details Cooling module with CHIANTI generated cooling curves
!> @n The location of the tables is assumed to be in
!! src/cool_lib/coolingCHIANTI.tab

module cooling_schure

  implicit none
  real (kind=8), allocatable :: cooltab_chianti(:,:)

contains

  !=======================================================================
  !> @brief Initializes the DMC cooling
  !> @details Declares variables and reads table
  subroutine init_cooling_schure()

    implicit none

    allocate(cooltab_chianti(2,41))
    call read_table_chianti()

  end subroutine init_cooling_schure

  !=======================================================================
  !> @brief Reads the cooling curve table
  !> @details Reads the cooling curve table generated by CHUANTI,
  !! the location is assumed in /src/cool_lib/coolingCHIANTI.tab
  subroutine read_table_chianti()

    use parameters, only : workdir, master
    use globals, only : rank
#ifdef MPIP
    use mpi
#endif
    implicit none
    integer :: i, err
    real (kind=8) :: a, b

    if(rank == master) then
      open(unit=10,file= trim(workdir)//'../src/cool_lib/coolingCHIANTI.tab',&
           status='old')
      do i=1,41
        read(10,*) a, b
        cooltab_chianti(1,i)=a
        cooltab_chianti(2,i)=b
      end do
      close(unit=10)
    endif
#ifdef MPIP
    call mpi_bcast(cooltab_chianti,82,mpi_double_precision,0,mpi_comm_world,err)
#endif

  end subroutine read_table_chianti

  !=======================================================================
  !> @brief Returns the cooling coefficient interpolating the table
  !> @param real [in] T : Temperature K
  function get_lambda(T)

    implicit none
    real , intent(in) :: T
    integer           :: if1
    real (kind=8)     :: get_lambda, T0, T1, C0, C1

    if(T.gt.1e8) then
      get_lambda=0.21e-26*Sqrt(real(T,8))
    else
      if1=int(log10(T)*10)-39
      T0=cooltab_chianti(1,if1)
      c0=cooltab_chianti(2,if1)
      T1=cooltab_chianti(1,if1+1)
      c1=cooltab_chianti(2,if1+1)
      get_lambda=(c1-c0)*(real(T,8)-T0)/(T1-T0)+c0
    end if

  end function get_lambda

  !=======================================================================
  !> @brief High level wrapper to apply cooling with CHIANTI tables
  !> @details High level wrapper to apply cooling with CHIANTI tables
  !> @n cooling is applied in the entire domain and updates both the
  !! conserved and primitive variables
  subroutine coolingschure()

    use parameters, only : nx, ny, nz, cv, Psc, tsc, dif_rad, mhd, n1_chem
    use globals, only : u, primit, dt_CFL
    use hydro_core, only : u2prim
    use difrad
    use network
    implicit none
    real                 :: T ,Eth0, dens
    real, parameter      :: Tmin=1000.
    real (kind=8)        :: ALOSS, Ce, gain
    integer              :: i, j, k, iiHI, iiHeI, iiHeII
    real                 :: dt_seconds, Tprime, T1, ch_factor

    dt_seconds = dt_CFL*tsc
    iiHI   = n1_chem - 1 + iHI
    iiHeI  = n1_chem - 1 + iHeI
    iiHeII = n1_chem - 1 + iHeII

    do k=1,nz
      do j=1,ny
        do i=1,nx

          !   get the primitives (and T)
          call u2prim(u(:,i,j,k),primit(:,i,j,k),T)

          if(T > Tmin) then

            Eth0=cv*primit(5,i,j,k)

            Aloss=get_lambda(T)

            if (dif_rad) then
              !  energy per photo ionization from Black 1981 (in erg)
              gain = phHI(i,j,k)   * u(iiHI  ,i,j,k) * 7.75e-12                &
                   + phHeI(i,j,k)  * u(iiHeI ,i,j,k) * 2.19e-11                &
                   + phHeII(i,j,k) * u(iiHeII,i,j,k) * 3.10e-11

              Tprime=max( gain*T/aloss,Tmin)
              if(Tprime < Tmin) print*, 'Tprime=', Tprime

            else

              Tprime=Tmin

            end if

            dens=primit(1,i,j,k)
            Ce=(Aloss*dens**2)/(Eth0*Psc)  ! cgs

            T1 = Tprime+(T-Tprime)*exp(-Ce*dt_seconds) !# new temperature
            T1 = max( T1, 0.2*T )
            T1 = min( T1, 5.0*T )

            ch_factor = real(T1)/T

            !  apply cooling to primitive and conserved variables
            primit(5,i,j,k)=primit(5,i,j,k)*ch_factor
            !  update total energy density

            u(5,i,j,k) = cv*primit(5,i,j,k)                                    &
                         + 0.5*primit(1,i,j,k)*(  primit(2,i,j,k)**2           &
                                                + primit(3,i,j,k)**2           &
                                                + primit(4,i,j,k)**2  )
#ifdef BFIELD
          if (mhd) then
            u(5,i,j,k) = u(5,i,j,k) + 0.5*(  primit(6,i,j,k)**2                  &
                                           + primit(7,i,j,k)**2                  &
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
