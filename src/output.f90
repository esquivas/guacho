!=======================================================================
!> @file output.f90
!> @brief Writes Output
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

!> @brief Writes output
!> @details This module writes the ouput in the formats specified in 
!! the makefile

module output

contains

!=======================================================================
!> @brief Writes output
!> @details Writes output, the format is chosen in makefile
!! @n Supported formats are
!! *.dat (ASCII), *.bin and VTK (both BINARY), Silo (+hdf5) 
!> @param integer [in] itprint : number of output

subroutine write_output(itprint)

  use parameters
  use globals
#ifdef OUTSILO
  use Out_Silo_Module
#endif
#ifdef RADDIFF
  use difrad
#endif
#ifdef OUTVTK
  use hydro_core, only : u2prim
#endif
  implicit none
  integer, intent(in) :: itprint
  character (len=128) :: file1
#ifdef OUTVTK
  character (len=50) cbuffer
  real  :: t, x0,y0, z0
  real, dimension(neq) :: prim
  integer :: npoints
#ifdef MPIP
  character (len=128) :: file2
#endif
#endif
#ifdef MPIP
  integer :: err
#endif
  !real, dimension(neq+1,nxmin:nxmax,nymin:nymax,nzmin:nzmax):: uout
  character(1), parameter  :: lf = char(10)
  integer :: unitout
  
  integer :: ip
#if defined (OUTDAT) || defined(DIVBCORR) || defined(OUTVTK)
  integer ::  i, j, k
#endif
#ifdef DIVBCORR
 real :: divB(nx,ny,nz)
#endif

  !   output *.dat
#ifdef OUTDAT
#ifdef MPIP
  write(file1,'(a,i3.3,a,i3.3,a)') trim(outputpath)//'DAT/points',rank,'.',itprint,'.dat'
  unitout=rank+10
#else
  write(file1,'(a,i3.3,a)')  trim(outputpath)//'DAT/points',itprint,'.dat'
  unitout=10
#endif

  open(unit=unitout,file=file1,status='unknown')
  !
  do i=1,nx
     do j=1,ny
        do k=1,nz
           !uu(:)=u(:,i,j,k)
           !   calculate the primitive variables (which are the output)
           !   and write them to file           call u2prim(uu,prim,T)
           !write(unitout,'(9(e16.8))') float(i)*dx,float(j)*dy,float(k)*dz &
           !     ,prim,t
           write(unitout,'(14(e16.8))') float(i)*dx,float(j)*dy,float(k)*dz &
                ,primit(:,i,j,k)!,t
        end do
     end do
  end do

  close(unitout)
#endif

     !   output *.bin ( used for warm start )
#ifdef OUTBIN

#ifdef MPIP
    write(file1,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'BIN/points',rank,'.',itprint,'.bin'
  unitout=rank+10
#else
  write(file1,'(a,i3.3,a)')  trim(outputpath)//'BIN/points',itprint,'.bin'
  unitout=10
#endif
     ! take turns
     do ip=0, np-1
       if(rank == ip) then
         open(unit=unitout,file=file1,status='unknown',form='unformatted', &
                       convert='LITTLE_ENDIAN')
          write(unitout) u(:,:,:,:)

         close(unitout)

       end if
#ifdef MPIP
         call mpi_barrier(mpi_comm_world, err)
#endif
     end do

#endif
#ifdef MPIP
     print'(i3,a,a)',rank," wrote file:",trim(file1)
     call mpi_barrier(mpi_comm_world, err)
#endif

     !   write the emissvity and photoionizing rate
     !   if diffuse radiation enabled
#ifdef RADDIFF

     ! take turns
     do ip=0, np-1
       if(rank == ip) then

          write(file1,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'BIN/em-',rank,'.',itprint,'.bin'
          unitout=10+rank
          open(unit=unitout,file=file1,status='unknown',form='unformatted', &
               convert='LITTLE_ENDIAN')
          write (unitout) em(:,:,:)
          close(unitout)
          print'(i3,a,a)',rank," wrote file:",trim(file1)

          write(file1,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'BIN/ph-',rank,'.',itprint,'.bin'
          unitout=10+rank
          open(unit=unitout,file=file1,status='unknown',form='unformatted', &
                              convert='LITTLE_ENDIAN')
          write (unitout) ph(:,:,:)
          close(unitout)
          print'(i3,a,a)',rank," wrote file:",trim(file1)
          
       end if
       call mpi_barrier(mpi_comm_world, err)
    end do
             
#endif

#ifdef DIVBCORR     
    !   This is a hack to write div(B) to plot it easily
    !  compute div(B)
    do i=1,nx
      do j=1,ny
        do k=1,nz
          divB(i,j,k) = (u(6,i+1,j,k)-u(6,i-1,j,k))/(2.*dx) + &
                        (u(7,i,j+1,k)-u(7,i,j-1,k))/(2.*dy) + &
                        (u(8,i,j,k+1)-u(8,i,j,k-1))/(2.*dz)
        end do
      end do
    end do

    ! take turns to write to disk
     do ip=0, np-1
       if(rank == ip) then
          write(file1,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'BIN/divB-',rank,'.',itprint,'.bin'
          unitout=10+rank

          open(unit=unitout,file=file1,status='unknown',form='unformatted', &
               convert='LITTLE_ENDIAN')

          write (unitout) divB(:,:,:)
          close(unitout)
          
       end if
       call mpi_barrier(mpi_comm_world, err)
    end do
#endif


     !   output *.vtk
#ifdef OUTVTK

  !   write to .visit file to include several subsets
#ifdef MPIP
  if (rank.eq.0) then
     write(file2,'(a)') trim(outputpath)//'VTK/.visit'
     if (itprint.eq.0) then
        open(unit=7,file=file2,status='unknown',form='formatted')
        write(7,'(a,i)') '!NBLOCKS ',np
        do i=0,np-1
           write(7,'(a,i3.3,a,i3.3,a,i3.3,a)')  trim(outputpath)//'VTK/',nxtot,'-',i,'.',itprint,'.vtk'
        end do
        close(7)
     else
        open(unit=7,file=file2,status='unknown',form='formatted',position='append')
        do i=0,np-1
           write(7,'(a,i3.3,a,i3.3,a,i3.3,a)')  trim(outputpath)//'VTK/',nxtot,'-',i,'.',itprint,'.vtk'
        end do
        close(7)
     end if
  end if

  write(file1,'(a,i3.3,a,i3.3,a,i3.3,a)')  trim(outputpath)//'VTK/',nxtot,'-',rank,'.',itprint,'.vtk'
  unitout=rank+10
#else
  write(file1,'(a,i3.3,a,i3.3,a)') trim(outputpath)//'VTK/out-',nxtot,'-',itprint,'.vtk'
  unitout=10
#endif
  open(unit=unitout,file=file1,status='unknown',form='formatted')!,   &
!!$       form='formatted',recordtype='STREAM', &
!!$       action='write',convert='LITTLE_ENDIAN',    &
!!$       access='sequential' )

  !   write the header
  x0=( float(coords(0)*nx) )*dx
  y0=( float(coords(1)*ny) )*dy
  z0=( float(coords(2)*nz) )*dz
  npoints=(nx+1)*(ny+1)*(nz+1)

  write(unitout,'(a)') '# vtk DataFile Version 2.0 '
  write(unitout,'(a)') 'output from diable'
  write(unitout,'(a)') 'ASCII'
  write(unitout,'(a)') 'DATASET STRUCTURED_POINTS'
  write(unitout, '("DIMENSIONS ",(3i6,1x))') ,nx+1,ny+1,nz+1
 ! write(unitout,'(a)') trim(cbuffer),lf
  write(unitout, '("ORIGIN "    ,3e15.7)'), x0,y0,z0
  !write(unitout,'(a)') trim(cbuffer),lf
  write(unitout, '("SPACING",3e15.7)'), dx,dy,dz
 ! write(unitout,'(a)') trim(cbuffer),lf

  !   writes the variables, scalars first then vectors

  write(unitout,'(a,i10)') 'POINT_DATA ',npoints
 ! write(unitout) trim(cbuffer),lf
  !   DENSITY
  write(unitout,'(a)') 'SCALARS Density double 1'
  !write(unitout) trim(cbuffer),lf
  write(unitout,'(a)') 'LOOKUP_TABLE default'
  !write(unitout) trim(cbuffer),lf

  do k=0,nz
     do j=0,ny
           write(unitout,*) (primit(1,i,j,k)*rhosc, i=0,nx)
     end do
  end do

!!$  write(unitout) lf  

  !   GAS PRESSURE
  write(unitout,'(a)') 'SCALARS Pressure double 1'
  !write(unitout) trim(cbuffer),lf
  write(unitout,'(a)') 'LOOKUP_TABLE default'
  !write(unitout) trim(cbuffer),lf

  do k=0,nz
     do j=0,ny
           write(unitout,*) (primit(5,i,j,k)*psc, i=0,nx)
     end do
  end do
  !write(unitout) lf  

  !   TEMP
  write(unitout,'(a)') 'SCALARS Temperature double 1'
 ! write(unitout) trim(cbuffer),lf
  write(unitout,'(a)') 'LOOKUP_TABLE default'
  !write(unitout) trim(cbuffer),lf

  do k=0,nz
     do j=0,ny
        do i=0,nx
           call u2prim (u(:,i,j,k),prim,T)
           write(unitout,*) T
        end do
     end do
  end do
  !write(unitout) lf

  !   VELOCITY
  write(unitout,'(a)') 'VECTORS Velocity double'
  !write(unitout) trim(cbuffer),lf
     do k=0,nz
        do j=0,ny
           do i=0,nx
              write(unitout,*) primit(2,i,j,k)*sqrt(vsc2),    &
                               primit(3,i,j,k)*sqrt(vsc2),    &
                               primit(4,i,j,k)*sqrt(vsc2)
           end do
        end do
     end do
     !write(unitout) lf  

#ifdef PMHD
  !   MAGNETIC FIELD
  write(unitout,'(a)') 'VECTORS BField double'
  !write(unitout) trim(cbuffer),lf
     do k=0,nz
        do j=0,ny
           do i=0,nx
              write(unitout,*) primit(6,i,j,k)*bsc,           &
                               primit(7,i,j,k)*bsc,           &
                               primit(8,i,j,k)*bsc
           end do
        end do
     end do
     !write(unitout) lf  

#endif
#ifdef MHD
  !   MAGNETIC FIELD
  write(unitout,'(a)') 'VECTORS BField double'
  !write(unitout) trim(cbuffer),lf
     do k=0,nz
        do j=0,ny
           do i=0,nx
              write(unitout,*) primit(6,i,j,k)*bsc,           &
                               primit(7,i,j,k)*bsc,           &
                               primit(8,i,j,k)*bsc
           end do
        end do
     end do
     !write(unitout) lf  

#endif

  close(unitout)
#endif

  !  output *.silo					      
#ifdef OUTSILO						     
							      
  call outputsilo(itprint)				      
							      
#endif							      

end subroutine write_output

!=======================================================================

end module output

!=======================================================================
