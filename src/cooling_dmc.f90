!=======================================================================
!> @file cooling_dmc.f90
!> @brief Cooling module with Dlgarno Mac Cray coronal cooling curve
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

!> @brief Cooling module with Dalgarno McCray coronal cooling curve
!> @details Cooling module with Dalgarno McCray coronal cooling curve
!> @n The location of the tables is assumed to be in 
!! src/DMClib/coolingDMC.tab, it is read by init subroutine

module cooling_dmc

#ifdef COOLINGDMC

  implicit none
  real (kind=8) :: cooltab(2,41)

contains

!> @brief Reads the cooling curve table
!> @details Reads the Dalgarno McCray cooling courve
!! the location is assumed in src/DMClib/coolingDMC.tab, 
!! it is read by init subroutine

subroutine read_table()

  use parameters, only : workdir, master
  use globals, only : rank
  implicit none
#ifdef MPIP
  include "mpif.h"
#endif
  integer :: i, err
  real (kind=8) :: a, b


  if(rank.eq.master) then
     open(unit=10,file=trim(workdir)//'../src/DMClib/coolingDMC.tab',status='old')
     do i=1,41
        read(10,*) a, b
        cooltab(1,i)=10.d0**(a)
        cooltab(2,i)=10.d0**(-b)
     end do
     close(unit=10)
  endif
#ifdef MPIP
  call mpi_bcast(cooltab,82,mpi_real8,0,mpi_comm_world,err)
#endif

end subroutine read_table

!=======================================================================

!> @brief Returns the cooling coefficient interpolating the table
!> @param real [in] T : Temperature K

function cooldmc(T)
  
  real , intent(in) :: T
  integer           :: if1
  real (kind=8)     :: cooldmc, T0, T1, C0, C1

  if(T.gt.1e8) then
    cooldmc=0.27D-26*Sqrt(dble(T))
  else
    if1=int(log10(T)*10)-39
    T0=cooltab(1,if1)
    c0=cooltab(2,if1)
    T1=cooltab(1,if1+1)
    c1=cooltab(2,if1+1)
    cooldmc=(c1-c0)*(dble(T)-T0)/(T1-T0)+c0
  end if

end function cooldmc

!=======================================================================

!> @brief High level wrapper to apply cooling with DMC table
!> @details High level wrapper to apply cooling with DMC table
!> @n cooling is applied in the entire domain and updates both the 
!! conserved and primitive variables

subroutine coolingdmc()

  use parameters, only : nx, ny, nz, cv, Psc, tsc
  use globals, only : u, primit, dt_CFL
  use hydro_core, only : u2prim
  implicit none
  real                 :: T ,Eth0, dens
  real, parameter :: Tmin=10000.
  real (kind=8)        :: ALOSS, Ce
  integer :: i, j, k
  real :: dt_seconds

  dt_seconds = dt_CFL*tsc

  do k=1,nz
     do j=1,ny
        do i=1,nx

           !   get the primitives (and T)
           call u2prim(u(:,i,j,k),primit(:,i,j,k),T)

           if(T > Tmin) then

              Eth0=cv*primit(5,i,j,k)

              Aloss=cooldmc(T)
              dens=primit(1,i,j,k)
              Ce=(Aloss*dble(dens)**2)/(Eth0*Psc)  ! cgs

              !  apply cooling to primitive and conserved variables
              primit(5,i,j,k)=primit(5,i,j,k)*exp(-ce*dt_seconds)

              !   u(neqdyn,new)=Ekin0+Eth_new
              u(5,i,j,k)=u(5,i,j,k)-Eth0+cv*primit(5,i,j,k)

           end if

        end do
     end do
  end do

end subroutine coolingdmc

!======================================================================

#endif

end module cooling_dmc

!======================================================================