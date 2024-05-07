!=======================================================================
!> @file Out_Silo_Module.f90
!> @brief Output in Silo Format
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

!> @brief Output in Silo (+HDF5) Format
!> @details This module writes the ouput in SILO (HDF5) format

module  Out_Silo_Module

#ifdef OUT_SILO

  use parameters
  use globals
  include "silof90.inc"

contains

  !=======================================================================
  !> @brief  Writes Data, one file per processor
  !> @details  Writes Data in silo format one file per processor
  !> @param integer [in] itprint : number of output
  subroutine writeblocks(itprint)

    use globals
    use parameters
    implicit none
    integer, intent(in) :: itprint
    character (len=128) file
    integer ::  unitout
    integer :: err, ierr, meshdims(3), vardims(3), optlist
    integer :: i, j, k
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real, allocatable :: xax(:), yax(:), zax(:)
    integer :: lo_offset(3), hi_offset(3)

    !  open the silo file
#ifdef MPIP
    write(file,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'SILO/BLOCKS/out',rank,&
                                     '.',itprint,'.silo'
    unitout=rank+10
#else
    write(file,'(a,i3.3,a)') trim(outputpath)//'SILO/BLOCKS/out',itprint,'.silo'
    unitout=10
#endif

    !  determine the size of the blocks to be dumped, external ghost zones
    !  are not written, only 1 ghost cell

    !  in X
    if (coords(0).eq.0) then
      imin=1
      lo_offset(1)=0
    else
      imin=0
      lo_offset(1)=1
    end if
    if ( coords(0).eq.(MPI_NBX-1) ) then
        imax=nx
        hi_offset(1)=0
      else
        imax=nx+1
        hi_offset(1)=1
    end if

    !  in Y
    if (coords(1).eq.0) then
      jmin=1
      lo_offset(2)=0
    else
      jmin=0
      lo_offset(2)=1
    end if
    if ( coords(1).eq.(MPI_NBY-1) ) then
      jmax=ny
      hi_offset(2)=0
    else
      jmax=ny+1
      hi_offset(2)=1
    end if

    !  in Z
    if (coords(2).eq.0) then
      kmin=1
      lo_offset(3)=0
    else
      kmin=0
      lo_offset(3)=1
    end if
    if ( coords(2).eq.(MPI_NBZ-1) ) then
      kmax=nz
      hi_offset(3)=0
    else
      kmax=nz+1
      hi_offset(3)=1
    end if

    !  number of cells
    vardims(1)=imax-imin+1
    vardims(2)=jmax-jmin+1
    vardims(3)=kmax-kmin+1
    !  number of cell interfaces
    meshdims(:)=vardims(:)+1

    !   generate the axes
    allocate ( xax(imin:imax+1) )
    allocate ( yax(jmin:jmax+1) )
    allocate ( zax(kmin:kmax+1) )

    do i=imin,imax+1
      xax(i)=(real( i +coords(0)*nx -1 ) )*dx
    end do

    do j=jmin,jmax+1
      yax(j)=(real( j +coords(1)*ny -1 ) )*dy
    end do

    do k=kmin,kmax+1
      zax(k)=(real( k +coords(2)*nz -1 ) )*dz
    end do

    !  open the silo file
    ierr = dbcreate(file, len_trim(file), DB_CLOBBER, DB_LOCAL,                &
                     "Comment about the data", 22, DB_HDF5, unitout)

    !  put mesh
    !   options to flag ghost zones
    err = dbmkoptlist(2, optlist)
    err = dbaddiopt(optlist, DBOPT_HI_OFFSET, hi_offset)
    err = dbaddiopt(optlist, DBOPT_LO_OFFSET, lo_offset)
    !   actual mesh

    err = dbputqm (unitout, "quadmesh", 8, "xc", 2, "yc", 2, "zc", 2,          &
                   xax, yax, zax, meshdims, ndim, DB_DOUBLE, DB_COLLINEAR,     &
                   optlist, ierr)
    err = dbfreeoptlist(optlist)

    !  write density
    err = dbputqv1(unitout, "rho", 3, "quadmesh", 8,                           &
                   primit(1,imin:imax,jmin:jmax,kmin:kmax), vardims, ndim,     &
                  DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)

    !  velocity components
    err = dbputqv1(unitout, "vx", 2, "quadmesh", 8,                            &
                   primit(2,imin:imax,jmin:jmax,kmin:kmax), vardims, ndim,     &
                   DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)

    err = dbputqv1(unitout, "vy", 2, "quadmesh", 8,                            &
                   primit(3,imin:imax,jmin:jmax,kmin:kmax), vardims, ndim,     &
                   DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)

    err = dbputqv1(unitout, "vz", 2, "quadmesh", 8,                            &
                   primit(4,imin:imax,jmin:jmax,kmin:kmax), vardims, ndim,     &
                   DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)

    !  write Gas pressure
    err = dbputqv1(unitout, "Pth", 3, "quadmesh", 8,                           &
                  primit(5,imin:imax,jmin:jmax,kmin:kmax), vardims, ndim,      &
                  DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)

#ifdef BFIELD
    !  B components
    err = dbputqv1(unitout, "bx", 2, "quadmesh", 8,                            &
                  primit(6,imin:imax,jmin:jmax,kmin:kmax), vardims, ndim,      &
                  DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)

    err = dbputqv1(unitout, "by", 2, "quadmesh", 8,                            &
                  primit(7,imin:imax,jmin:jmax,kmin:kmax), vardims, ndim,      &
                  DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)

    err = dbputqv1(unitout, "bz", 2, "quadmesh", 8,                            &
                  primit(8,imin:imax,jmin:jmax,kmin:kmax), vardims, ndim,      &
                  DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)
#endif

    !  close the silo file
    ierr = dbclose(unitout)

    print'(i3,a,a)',rank," wrote file : ",trim(file)

  end subroutine writeblocks

  !=======================================================================
  !> @brief Writes the Master File
  !> @details Writes the master file with the metadata and multivars
  !> @param integer [in] itprint : number of output
  subroutine writemaster(itprint)

    use globals
    use parameters
    implicit none
    integer, intent(in) :: itprint
    character (len=128) :: buf
    character (len=4096) :: names
    integer ::  unitout, err, i, Lstring, ierr, oldlen
    integer ::  lnames(np), types(np)

    !  create file name and open it
    write(buf,'(a,i3.3,a)') trim(outputpath)//'SILO/master',itprint,'.root'
    err = dbcreate(buf, len_trim(buf) , DB_CLOBBER, DB_LOCAL, "multimesh root",&
                   14, DB_HDF5, unitout)
    if(unitout.eq.-1) then
      write (6,*) "Could not create Silo file!"
      return
    endif

    !  add multimesh
    do i=0, np-1
      write(buf,'(a,i3.3,a,i3.3,a)') './BLOCKS/out',i,'.',    &
            itprint,'.silo:quadmesh'
      Lstring=len_trim(buf)
      names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
      lnames(i+1) =  Lstring
      types (i+1) = DB_QUAD_RECT
    end do
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(lnames(1))
    err = dbputmmesh(unitout, "mesh", 4, np, names, lnames, types, DB_F77NULL, &
                     ierr)
    err = dbset2dstrlen(oldlen)

    !  add multivars
    !  density
    do i=0, np-1
      write(buf,'(a,i3.3,a,i3.3,a)') './BLOCKS/out',i,'.',    &
            itprint,'.silo:rho'
      Lstring=len_trim(buf)
      names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
      lnames(i+1)=  Lstring
      types (i+1) = DB_QUADVAR
    end do
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(lnames(1))
    err = dbputmvar(unitout, "rho", 3, np, names, lnames,types, DB_F77NULL,    &
                   ierr)
    err = dbset2dstrlen(oldlen)

    !  velocity components
    !V_x
    do i=0, np-1
      write(buf,'(a,i3.3,a,i3.3,a)') './BLOCKS/out',i,'.',     &
           itprint,'.silo:vx'
      Lstring=len_trim(buf)
      names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
      lnames(i+1)=  Lstring
      types (i+1) = DB_QUADVAR
    end do
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(lnames(1))
    err = dbputmvar(unitout, "vx", 2, np, names, lnames,types, DB_F77NULL, ierr)
    err = dbset2dstrlen(oldlen)
    !V_y
    do i=0, np-1
      write(buf,'(a,i3.3,a,i3.3,a)') './BLOCKS/out',i,'.',    &
            itprint,'.silo:vy'
      Lstring=len_trim(buf)
      names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
      lnames(i+1)=  Lstring
      types (i+1) = DB_QUADVAR
    end do
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(lnames(1))
    err = dbputmvar(unitout, "vy", 2, np, names, lnames,types, DB_F77NULL, ierr)
    err = dbset2dstrlen(oldlen)
    !V_z
    do i=0, np-1
      write(buf,'(a,i3.3,a,i3.3,a)') './BLOCKS/out',i,'.',    &
            itprint,'.silo:vz'
      Lstring=len_trim(buf)
      names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
      lnames(i+1)=  Lstring
      types (i+1) = DB_QUADVAR
    end do
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(lnames(1))
    err = dbputmvar(unitout, "vz", 2, np, names, lnames,types, DB_F77NULL, ierr)
    err = dbset2dstrlen(oldlen)

    !  Thermal pressure
    do i=0, np-1
      write(buf,'(a,i3.3,a,i3.3,a)')  './BLOCKS/out',i,'.',    &
            itprint,'.silo:Pth'
      Lstring=len_trim(buf)
      names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
      lnames(i+1)=  Lstring
      types (i+1) = DB_QUADVAR
    end do
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(lnames(1))
    err = dbputmvar(unitout, "Pth", 3, np, names, lnames,types, DB_F77NULL,ierr)
    err = dbset2dstrlen(oldlen)

    if (pmhd .or. mhd) then
      !  B components
      !B_x
      do i=0, np-1
        write(buf,'(a,i3.3,a,i3.3,a)') './BLOCKS/out',i,'.',  &
              itprint,'.silo:bx'
        Lstring=len_trim(buf)
        names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
        lnames(i+1)=  Lstring
        types (i+1) = DB_QUADVAR
      end do
      oldlen = dbget2dstrlen()
      err = dbset2dstrlen(lnames(1))
      err = dbputmvar(unitout,"bx",2, np, names, lnames, types, DB_F77NULL,ierr)
      err = dbset2dstrlen(oldlen)
      !B_y
      do i=0, np-1
        write(buf,'(a,i3.3,a,i3.3,a)') './BLOCKS/out',i,'.',  &
              itprint,'.silo:by'
        Lstring=len_trim(buf)
        names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
        lnames(i+1)=  Lstring
        types (i+1) = DB_QUADVAR
      end do
      oldlen = dbget2dstrlen()
      err = dbset2dstrlen(lnames(1))
      err = dbputmvar(unitout,"by",2, np, names, lnames, types, DB_F77NULL,ierr)
      err = dbset2dstrlen(oldlen)
      !B_z
      do i=0, np-1
        write(buf,'(a,i3.3,a,i3.3,a)') './BLOCKS/out',i,'.',  &
              itprint,'.silo:bz'
        Lstring=len_trim(buf)
        names(i*Lstring+1:(i+1)*Lstring ) = buf(1:Lstring)
        lnames(i+1)=  Lstring
        types (i+1) = DB_QUADVAR
      end do
      oldlen = dbget2dstrlen()
      err = dbset2dstrlen(lnames(1))
      err = dbputmvar(unitout,"bz",2, np, names, lnames, types, DB_F77NULL,ierr)
      err = dbset2dstrlen(oldlen)
    endif

    !  close master file
    err = dbclose(unitout)

    print*, 'Finished writing master file'

  end subroutine writemaster

  !=======================================================================
  !> @brief Upper level wrapper
  !> @details Upper level wrapper for the SILO output
  !> @param integer [in] itprint : number of output
  subroutine write_silo(itprint)

    implicit none
    integer, intent(in) :: itprint
    integer :: ip, err
    !   al processors dumps their output in a separate file
    ! take turns
    do ip=0, np-1
      if(rank == ip) then
        call writeblocks(itprint)
      end if
      call mpi_barrier(mpi_comm_world, err)
    end do

    !  if MPI enabled, the master node wirtes a master file
#ifdef MPIP
    if (master.eq.rank) call writemaster(itprint)
#endif

    end subroutine write_silo

!=======================================================================

#endif

end module Out_Silo_Module
