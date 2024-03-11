!=======================================================================
!> @file Out_VTK_Module.f90
!> @brief Output in VTK Format
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

!> @brief Output in VTK format
!> @details This module writes the ouput in VTK format

module  Out_VTK_Module

  use parameters
  use globals

contains

!=======================================================================

!> @brief  Writes Data, one file per processor
!> @details  Writes Data in VTK format one file per processor
!> @param integer [in] itprint : number of output

subroutine write_VTK(itprint)
  use hydro_core, only : u2prim
  implicit none

  integer, intent(in) :: itprint
  character (len=128) :: file1
  character (len=128) :: cbuffer
  real  :: t, x0,y0, z0
  real, dimension(neq) :: prim
  integer :: nCells, unitout, i, j, k
  character, parameter  :: lf = char(10)
#ifdef MPIP
  character (len=128) :: file2
#endif

  !   write to .visit file to include several subsets
#ifdef MPIP
  if (rank.eq.0) then
     write(file2,'(a)') trim(outputpath)//'VTK/master.visit'
     if (itprint.eq.itprint0) then
        open(unit=7,file=file2,status='replace',form='formatted')
        write(7,'(a,i0)') '!NBLOCKS ',np
        do i=0,np-1
           write(7,'(a,i3.3,a,i3.3,a)')  'out-',i,'.',itprint,'.vtk'
        end do
        close(7)
     else
        open(unit=7,file=file2,status='old',form='formatted',position='append')
        do i=0,np-1
           write(7,'(a,i3.3,a,i3.3,a)')  'out-',i,'.',itprint,'.vtk'
        end do
        close(7)
     end if
  end if

  write(file1,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'VTK/out-',rank,'.',itprint,'.vtk'
  unitout=rank+10
#else
  write(file1,'(a,i3.3,a)') trim(outputpath)//'VTK/out-',itprint,'.vtk'
  unitout=10
#endif

  open(unit=unitout,file=file1,status='replace',access='stream', convert='BIG_ENDIAN')

  !   write the header
  x0=( real(coords(0)*nx) )*dx
  y0=( real(coords(1)*ny) )*dy
  z0=( real(coords(2)*nz) )*dz
  nCells = nx*ny*nz

  write(cbuffer,'(a)') '# vtk DataFile Version 3.0 '
  write(unitout) trim(cbuffer),lf

  write(cbuffer,'(a)') 'output from Guacho-3D'
  write(unitout) trim(cbuffer),lf

  write(cbuffer,'(a)') 'BINARY'
  write(unitout) trim(cbuffer),lf

  write(cbuffer,'(a)') 'DATASET STRUCTURED_POINTS'
  write(unitout) trim(cbuffer),lf

  write(cbuffer, '("DIMENSIONS ",(i0,1x,i0,1x,i0))') nx+1,ny+1,nz+1
  write(unitout) trim(cbuffer),lf

  write(cbuffer, '("ORIGIN "    ,3e15.7)') x0*rsc,y0*rsc,z0*rsc
  write(unitout) trim(cbuffer),lf

  write(cbuffer, '("SPACING",3e15.7)') dx*rsc,dy*rsc,dz*rsc
  write(unitout) trim(cbuffer),lf

  !write(cbuffer,'(a)') 'TIME 1 1 double'
  !write(unitout) trim(cbuffer),lf
  !write(unitout) real(itprint*time*tsc,kind=8)
  !write(unitout) lf

  !  writes the variables, scalars first then vectors
  write(cbuffer,'(a,i0)') 'CELL_DATA ',nCells
  write(unitout) trim(cbuffer),lf

  !  Density
  write(cbuffer,'(a)') 'FIELD FieldData 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a,1x,i0,1x,i0,a)') 'Density',1,nCells,' float'
  write(unitout) trim(cbuffer),lf

  do k=1,nz
     do j=1,ny
        do i=1, nx
           write(unitout) real(primit(1,i,j,k)*rhosc,4)
        end do
     end do
  end do
  write(unitout) lf

 !  Gas pressure
 write(cbuffer,'(a)') 'FIELD FieldData 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a,1x,i0,1x,i0,a)') 'Thermal_Pressure',1,nCells,' float'
  write(unitout) trim(cbuffer),lf

  do k=1,nz
    do j=1,ny
     do i=1,nx
           write(unitout) real(primit(5,i,j,k)*Psc,4)
        end do
     end do
  end do
  write(unitout) lf

  !  Temperature
  write(cbuffer,'(a)') 'FIELD FieldData 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a,1x,i0,1x,i0,a)') 'Temperature',1,nCells,' float'
  write(unitout) trim(cbuffer),lf

  do k=1,nz
     do j=1,ny
        do i=1,nx
           call u2prim (u(:,i,j,k),prim,T)
           write(unitout) real(T,4)
        end do
     end do
  end do
  write(unitout) lf

  !  Velocity
  write(cbuffer,'(a)') 'VECTORS Velocity float'
  write(unitout) trim(cbuffer),lf
   do k=1,nz
      do j=1,ny
         do i=1,nx
            write(unitout)  real(primit(2,i,j,k)*vsc,4),    &
                            real(primit(3,i,j,k)*vsc,4),    &
                            real(primit(4,i,j,k)*vsc,4)
         end do
      end do
   end do
   write(unitout) lf

  if (pmhd .or. mhd) then
  !  Magnetic field
    write(cbuffer,'(a)') 'VECTORS BField float'
    write(unitout) trim(cbuffer),lf
     do k=1,nz
        do j=1,ny
           do i=1,nx
              write(unitout)  real(primit(6,i,j,k)*bsc,4),           &
                              real(primit(7,i,j,k)*bsc,4),           &
                              real(primit(8,i,j,k)*bsc,4)
           end do
        end do
     end do
     write(unitout) lf
  end if

  close(unitout)

  print'(i3,a,a)',rank," wrote file : ",trim(file1)

end subroutine write_VTK

!=======================================================================

end module Out_VTK_Module

!=======================================================================
