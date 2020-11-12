!=======================================================================
!> @file Out_BIN_Module.f90
!> @brief Output in BIN Format
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

!> @brief Output in BIN format
!> @details This module writes the ouput in BIN format

module  Out_BIN_Module

  use parameters
  use globals
  use constants
contains

!=======================================================================
!> @brief Writes header
!> @details Writes header for binary input
!> @param integer [in] unit : number of logical unit

subroutine write_header(unit, neq_out, nghost_out)
  implicit none
  integer, intent(in) :: unit, neq_out, nghost_out
  character, parameter  :: lf = char(10)
  character (len=128) :: cbuffer

  !  Write ASCII header
  write(unit) "**************** Output for Guacho v1.3****************",lf

  write(cbuffer,'("Dimensions    : ", i0,1x,i0,1x,i0)') NX, NY, NZ
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Spacings      : ", 3(es10.3))') dX, dY, dZ
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Block Origin, cells    : ", i0,1x,i0,1x,i0)') &
        coords(0)*nx, coords(1)*ny, coords(2)*nz
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("MPI blocks (X, Y, Z)   : ", i0,1x,i0,1x,i0)') &
        MPI_NBX, MPI_NBY, MPI_NBZ
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Number of Equations/dynamical ones  ", i0,"/",i0)') &
        neq_out, neqdyn
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Number of Ghost Cells  ", i0)') nghost_out
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Scalings ", a )')
  write(unit) trim(cbuffer), lf

  write(cbuffer, '("r_sc: ",es10.3," v_sc: ",es10.3," rho_sc: ",es10.3)') &
    rsc, vsc, rhosc
  write(unit) trim(cbuffer), lf

   write(cbuffer, '("Specfic heat at constant volume Cv: ",f7.2)') cv
  write(unit) trim(cbuffer), lf

#ifdef DOUBLEP
  write(unit) "Double precision 8 byte floats",lf
#else
  write(unit) "Double precision 4 byte floats",lf
#endif

  write(unit) "*******************************************************",lf
  write(unit) achar(255),lf
#ifdef DOUBLEP
  write(unit) 'd'
#else
  write(unit) 'f'
#endif
  write(unit) nx, ny, nz
  write(unit) dx, dy, dz
  write(unit) coords(0)*nx, coords(1)*ny, coords(2)*nz
  write(unit) MPI_NBX, MPI_NBY, MPI_NBZ
  write(unit) neq_out, neqdyn
  write(unit) nghost_out
  write(unit) rsc, vsc, rhosc
  write(unit) cv

end subroutine write_header

!=======================================================================

!> @brief  Writes Data, one file per processor
!> @details  Writes Data in BIN format one file per processor
!> @param integer [in] itprint : number of output
subroutine write_BIN(itprint)

  use difrad
  use self_gravity, only : phi_grav, four_pi_G
  use radpress, only : beta
  implicit none
  integer, intent(in) :: itprint
  character (len=128) :: file1

#ifdef MPIP
  integer :: err
#endif
  integer :: unitout
  integer :: ip
  integer ::  i, j, k
#ifdef BFIELD
  real, allocatable :: divB(:,:,:)
#ENDIF
  real, allocatable :: rho_grav(:,:,:)


#ifdef MPIP
  write(file1,'(a,i3.3,a,i3.3,a)')  &
        trim(outputpath)//'BIN/points',rank,'.',itprint,'.bin'
  unitout=rank+10
#else
  write(file1,'(a,i3.3,a)')  trim(outputpath)//'BIN/points',itprint,'.bin'
  unitout=10
#endif

  ! take turns
  do ip=0, np-1
   if(rank == ip) then
    open(unit=unitout,file=file1,status='replace',access='stream')

    ! write header, then data
    call write_header(unitout,neq,nghost)
    if (riemann_solver == SOLVER_HLLE_SPLIT_ALL) then
       do k =nzmin,nzmax
          do j = nymin,nymax
             do i = nxmin, nxmax
                write(unitout) u(1,i,j,k)+primit0(1,i,j,k)
                write(unitout) u(2,i,j,k)
                write(unitout) u(3,i,j,k)
                write(unitout) u(4,i,j,k)
                write(unitout) u(5,i,j,k)+cv*primit0(5,i,j,k)+&
              0.5*(primit0(6,i,j,k)**2+primit0(7,i,j,k)**2+primit0(8,i,j,k)**2)
                write(unitout) u(6,i,j,k)+primit0(6,i,j,k)
                write(unitout) u(7,i,j,k)+primit0(7,i,j,k)
                write(unitout) u(8,i,j,k)+primit0(8,i,j,k)
             end do
          end do
       end do

    else
       write(unitout) u(:,:,:,:)
    endif

    close(unitout)
    print'(i3,a,a)',rank," wrote file:",trim(file1)

   end if
#ifdef MPIP
     call mpi_barrier(mpi_comm_world, err)
#endif
  end do

     !   write the emissvity and photoionizing rate
     !   if diffuse radiation enabled
  if (dif_rad) then

       ! take turns
    do ip=0, np-1
      if(rank == ip) then

        write(file1,'(a,i3.3,a,i3.3,a)') &
              trim(outputpath)//'BIN/em-',rank,'.',itprint,'.bin'
        unitout=10+rank
        open(unit=unitout,file=file1,status='replace', access='stream')
        call write_header(unitout,1,0)
        write (unitout) em(:,:,:)
        close(unitout)
        print'(i3,a,a)',rank," wrote file:",trim(file1)

        write(file1,'(a,i3.3,a,i3.3,a)') &
              trim(outputpath)//'BIN/ph-',rank,'.',itprint,'.bin'
        unitout=10+rank
        open(unit=unitout,file=file1,status='replace',access='stream')
        call write_header(unitout,1,0)
        write (unitout) ph(:,:,:)
        close(unitout)
        print'(i3,a,a)',rank," wrote file:",trim(file1)

      end if
#ifdef MPIP
        call mpi_barrier(mpi_comm_world, err)
#endif
    end do

  end if


  if (enable_self_gravity) then

       allocate ( rho_grav(nx,ny,nz) )
       do i=1,nx
         do j=1,ny
           do k=1,nz
              rho_grav(i,j,k) =                                                &
          (   ( phi_grav(i+1,j  ,k  ) + phi_grav(i-1,j  ,k  ) - 2.*phi_grav(i,j,k) ) / dx**2   &
           +  ( phi_grav(i  ,j+1,k  ) + phi_grav(i  ,j-1,k  ) - 2.*phi_grav(i,j,k) ) / dy**2   &
           +  ( phi_grav(i  ,j  ,k+1) + phi_grav(i  ,j  ,k-1) - 2.*phi_grav(i,j,k) ) / dz**2  )&
           /  four_pi_G

           end do
         end do
       end do

           ! take turns
        do ip=0, np-1
          if(rank == ip) then

            write(file1,'(a,i3.3,a,i3.3,a)') &
                  trim(outputpath)//'BIN/rho_grav-',rank,'.',itprint,'.bin'
            unitout=10+rank
            open(unit=unitout,file=file1,status='replace',access='stream')
            call write_header(unitout,1,0)
            write (unitout) rho_grav(1:nx,1:ny,1:nz)
            close(unitout)
            print'(i3,a,a)',rank," wrote file:",trim(file1)
            deallocate ( rho_grav )

            write(file1,'(a,i3.3,a,i3.3,a)') &
                  trim(outputpath)//'BIN/phi_grav-',rank,'.',itprint,'.bin'
            unitout=10+rank
            open(unit=unitout,file=file1,status='replace',access='stream')
            call write_header(unitout,1,0)
            write (unitout) phi_grav(1:nx,1:ny,1:nz)
            close(unitout)
            print'(i3,a,a)',rank," wrote file:",trim(file1)

          end if
#ifdef MPIP
            call mpi_barrier(mpi_comm_world, err)
#endif
        end do

  end if

#ifdef BFIELD
  if (dump_divb) then
    !   This is a hack to write div(B) to plot it easily
    !  compute div(B)
    allocate(divB(nx,ny,nz))

    do k=1,nz
      do j=1,ny
        do i=1,nx
          divB(i,j,k) = (u(6,i+1,j,k)-u(6,i-1,j,k))/(2.*dx) + &
                        (u(7,i,j+1,k)-u(7,i,j-1,k))/(2.*dy) + &
                        (u(8,i,j,k+1)-u(8,i,j,k-1))/(2.*dz)
        end do
      end do
    end do

    ! take turns to write to disk
    do ip=0, np-1
      if(rank == ip) then
        write(file1,'(a,i3.3,a,i3.3,a)') &
              trim(outputpath)//'BIN/divB-',rank,'.',itprint,'.bin'
        unitout=10+rank

        open(unit=unitout,file=file1,status='replace',access='stream')

        call write_header(unitout,1,0)
        write (unitout) divB(:,:,:)
        close(unitout)
        print'(i3,a,a)',rank," wrote file:",trim(file1)

      end if
#ifdef MPIP
        call mpi_barrier(mpi_comm_world, err)
#endif
    end do

    deallocate(divB)

  end if
#endif

  if (beta_pressure) then
    ! take turns to write to disk
    do ip=0, np-1
      if(rank == ip) then
        write(file1,'(a,i3.3,a,i3.3,a)') &
              trim(outputpath)//'BIN/betaB-',rank,'.',itprint,'.bin'
        unitout=10+rank

        open(unit=unitout,file=file1,status='replace',access='stream')

        call write_header(unitout,1,0)
        write (unitout) beta(:,:,:)
        close(unitout)
        print'(i3,a,a)',rank," wrote file:",trim(file1)
      end if
#ifdef MPIP
        call mpi_barrier(mpi_comm_world, err)
#endif
    end do
  end if

end subroutine write_BIN

!=======================================================================

end module Out_BIN_Module

!=======================================================================
