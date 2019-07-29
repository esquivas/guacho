module pic_module
  use parameters, only : N_MP, pic_distF, NBinsSEDMP
  use globals   , only : Q_MP0, Q_MP1, partOwner, n_activeMP, partID, &
                         divV, MP_SED
  implicit none

contains

  !================================================================
  ! @brief initialization of module
  subroutine init_pic()
    use parameters, only : nx, ny, nz
    implicit none

    if(pic_distF) then
      allocate( Q_MP0(N_MP,8) )
    else
      allocate( Q_MP0(N_MP,6) )
    end if
    allocate( Q_MP1(N_MP,3) )
    allocate( partID   (N_MP) )
    allocate( partOwner(N_MP) )
    if (pic_distF) then
      allocate( divV(0:nx+1,0:ny+1,0:nz+1) )
      allocate( MP_SED(N_MP, NBinsSEDMP) )
    endif

  end subroutine init_pic

  !================================================================
  subroutine deactivateMP(i_mp)
    use globals,    only : partID, n_activeMP
    implicit none
    integer, intent(in) :: i_mp

    !  deactivate particle
    partID(i_mp) = 0

    !  if necessary move tail index
    if (i_mp == n_activeMP) n_activeMP = n_activeMP - 1

  end subroutine deactivateMP

  !================================================================
  subroutine addMP(ID, ndata, Qdata, i_mp)
    use parameters, only : N_MP
    use globals,    only : partID,Q_MP0
    implicit none
    integer, intent(in)  :: ID,ndata
    real,    intent(in)  :: Qdata(ndata)
    integer, intent(out) :: i_mp

    ! sweep list and add element in empty slot
    do i_mp=1, n_mp
      if(partID(i_mp)==0) then
        partID(i_mp)  = ID
        Q_MP0(i_mp,1:ndata) = Qdata(1:ndata)
        if (i_mp > n_activeMP) n_activeMP = i_mp
        return
      end if
    end do

  end subroutine addMP
  !================================================================
  subroutine PICpredictor

    use globals,   only : primit, dt_CFL, rank, dx, dy, dz,comm3d
    use parameters
    use utilities, only : isInDomain, inWhichDomain
    implicit none
    integer :: i_mp, i, j, k, l, ind(3), dest, nLocSend, &
               sendLoc(0:np-1), sendList(0:np-1,0:np-1), &
               dataLoc(N_MP), iS, iR
    real    :: weights(8), SingleRec(6)
    integer:: status(MPI_STATUS_SIZE), err

    ! initialize send and recv lists
    dataLoc(:)    =  0
    sendLoc(:)    =  0
    sendlist(:,:) =  0
    nLocSend      =  0

    do i_mp=1, n_MP
      ! execute only if particle i_mp is in the active list
      if (partID(i_mp)/=0) then

        ! proceed further only if paricle is in domain
        if ( isInDomain( Q_MP0(i_mp,1:3) ) ) then

          ! Calculate interpolation reference and weights
          call interpBD(Q_MP0(i_mp,1:3),ind,weights)

          !  Interpolate the velocity field to particle position
          l=1
          Q_MP0(i_mp,4:6) = 0
          do k= ind(3),ind(3)+1
            do j=ind(2),ind(2)+1
              do i=ind(1),ind(1)+1
                Q_MP0(i_mp,4) = Q_MP0(i_mp,4) + primit(2,i,j,k)*weights(l)
                Q_MP0(i_mp,5) = Q_MP0(i_mp,5) + primit(3,i,j,k)*weights(l)
                Q_MP0(i_mp,6) = Q_MP0(i_mp,6) + primit(4,i,j,k)*weights(l)
                l = l + 1
              end do
            end do
          end do

          !  predictor step
          Q_MP1(i_mp,1:3) = Q_MP0(i_mp,1:3)+ dt_CFL*Q_MP0(i_mp,4:6)

        end if

        !  check if particle is leaving the domain
        dest = inWhichDomain(Q_MP1(i_mp,1:3))
        if( dest /= rank ) then
          !  count for MPI exchange
          sendLoc(dest) = sendLoc(dest) + 1
          nLocSend      = nLocSend      + 1
          dataLoc(nLocSend) = i_mp
          !print'(i2,a,i4,a,i4)', rank, ' *** particle ',partID(i_mp), ' is going to ',DEST
        end if

      end if
    end do

    !   consolidate list to have info of all send/receive operations
    call mpi_allgather(sendLoc(:),  np, mpi_integer, &
                       sendList, np, mpi_integer, comm3d,err)

    !  exchange particles
    do iR=0,np-1
      do iS=0,np-1

        if(sendList(iR,iS) /= 0) then
          if(iS == rank) then
            !print'(i0,a,i0,a,i0)', rank,'-->', IR, ':',sendList(iR,iS)
            do i=1,sendlist(iR,iS)
              !print*,'>>>',rank,partID(dataLoc(i)),dataLoc(i)
              call mpi_send( Q_MP0(dataLoc(i),1:6) , 6, mpi_real_kind ,IR, &
                            partID(dataLoc(i)), comm3d,err)
              !  deactivate particle from current processor
              call deactivateMP(dataLoc(i))
            end do

          endif
          if(iR == rank) then
            !print'(i0,a,i0,a,i0)', rank,'<--', IS, ':',sendList(iR,iS)
            do i=1,sendList(iR,iS)
              !print*,'<<<',rank,' will recv here from ',IS

              call mpi_recv(SingleRec(1:6), 6, mpi_real_kind, IS, mpi_any_tag, &
                            comm3d, status, err)
              !print*,'received successfuly', status(MPI_TAG)

              !  add current particle in list and data in new processor
              call addMP( status(MPI_TAG), 6, SingleRec(1:6), i_mp )

              !  recalculate predictor step for newcomer
              Q_MP1(i_mp,1:3) = Q_MP0(i_mp,1:3)+ dt_CFL*Q_MP0(i_mp,4:6)

            end do
          end if
        end if
      end do
    end do

  end subroutine PICpredictor

  !================================================================
  subroutine PICcorrector

    use globals,   only : primit, dt_CFL, rank, comm3d
    use parameters
    use utilities, only : inWhichDomain
    implicit none
    integer :: i_mp, i, j, k, l, ind(3)
    real    :: weights(8), vel1(3), SingleRec(3)
    integer :: dest, nLocSend, dataLoc(N_mp),  sendLoc(0:np-1), &
               sendList(0:np-1,0:np-1)
    integer :: status(MPI_STATUS_SIZE), err, iR, iS

    ! initialize send and recv lists
    dataLoc(:)    =  0
    sendLoc(:)    =  0
    sendlist(:,:) =  0
    nLocSend      =  0


    do i_mp=1, n_MP
      ! execute only if particle i_mp is in the active list
      if (partID(i_mp)/=0) then
        ! Calculate interpolation reference and weights
        call interpBD(Q_MP1(i_mp,1:3),ind,weights)
        !  Interpolate the velocity field to particle position, and add
        !  to the velocity from the corrector step
        l=1
        vel1(:) = 0
        do k= ind(3),ind(3)+1
          do j=ind(2),ind(2)+1
            do i=ind(1),ind(1)+1
              vel1(1) = vel1(1) + primit(2,i,j,k)*weights(l)
              vel1(2) = vel1(2) + primit(3,i,j,k)*weights(l)
              vel1(3) = vel1(3) + primit(4,i,j,k)*weights(l)
              l = l + 1
            end do
          end do
        end do

        !  corrector step
        Q_MP0(i_mp,1:3) = Q_MP0(i_mp,1:3) &
        + 0.5*dt_CFL*( Q_MP0(i_mp,4:6) + vel1(1:3) )

        !  check if particle is leaving the domain
        dest = inWhichDomain( Q_MP0(i_mp,1:3) )
        if (dest /= rank) then
          !  count for MPI exchange
          sendLoc(dest) = sendLoc(dest) + 1
          nLocSend      = nLocSend      + 1
          dataLoc(nLocSend) = i_mp
        end if

      end if
    end do

      !   consolidate list to have info of all send/receive operations
    call mpi_allgather(sendLoc(:),  np, mpi_integer, &
                       sendList, np, mpi_integer, comm3d,err)


    !  exchange particles
    do iR=0,np-1
      do iS=0,np-1

        if(sendList(iR,iS) /= 0) then
          if(iS == rank) then
            do i=1,sendlist(iR,iS)
              call mpi_send( Q_MP0(dataLoc(i),1:3) , 3, mpi_real_kind ,IR, &
                            partID(dataLoc(i)), comm3d,err)
              !  deactivate particle from current processor
              call deactivateMP(dataLoc(i))
            end do

          endif
          if(iR == rank) then
            do i=1,sendList(iR,iS)
              call mpi_recv(SingleRec(1:3), 3, mpi_real_kind, IS, mpi_any_tag, &
                            comm3d, status, err)
              !  add current particle in list and data in new processor
              call addMP( status(MPI_TAG), 3, SingleRec(1:3), i_mp )

            end do
          end if
        end if
      end do
    end do

  end subroutine PICcorrector

  !================================================================
  !  @brief bilineal interpolator
  !  @ details: returns weights and indices corresponding to a bilineal
  !> interpolation at for the particle position.
  !> @param real [in] pos(3) : Three dimensional position of particle,
  !> @param integer   ind(3) : Reference corner i0,j0,k0 indices of the cube of
  !> cells that are used in the interpolation
  !> @param real [out] weights(8) : weights associated with each corner of the
  !> interopolation in the following order
  !> 1-> i0,j0,k0,  3-> i0,j1,k0,  5-> i0,j0,k1,  7-> i0,j1,k1
  !> 2-> i1,j0,k0,  4-> i1,j1,k0,  6-> i1,j0,k1,  8-> i1,j1,k1
    subroutine interpBD(pos,ind,weights)

    use parameters, only : nx, ny, nz
    use globals,    only : coords, dx, dy, dz
    implicit none
    real,    intent(in)  :: pos(3)
    integer, intent(out) :: ind(3)
    real,    intent(out) :: weights(8)
    real                 :: x0, y0, z0

    ! get the index in the whole domain (integer part)
    x0 = pos(1)/dx
    y0 = pos(2)/dy
    z0 = pos(3)/dz

    ! shift to the local processor
    ind(1) = int(x0) - coords(0)*nx
    ind(2) = int(y0) - coords(1)*ny
    ind(3) = int(z0) - coords(2)*nz

    ! get the reminder
    x0 = x0 - int(x0)
    y0 = y0 - int(y0)
    z0 = z0 - int(z0)

    weights(1) = (1-x0)*(1-y0)*(1-z0)
    weights(2) =   x0  *(1-y0)*(1-z0)
    weights(3) = (1-x0)*  y0  *(1-z0)
    weights(4) =   x0  *  y0  *(1-z0)
    weights(5) = (1-x0)*(1-y0)*  z0
    weights(6) =   x0  *(1-y0)*  z0
    weights(7) = (1-x0)*  y0  *  z0
    weights(8) =   x0  *  y0  *  z0

    return

  end subroutine interpBD

  !================================================================
  ! @brief Writes pic output
  ! @details Writes position and particles velocities in a single file
  !> format is binary, one integer with the number of points in domain, and then
  !> all the positions and velocities ( in the default real precision)
  subroutine write_pic(itprint)

    use parameters, only : outputpath, np
    use globals,    only : rank
    implicit none
    integer, intent(in) :: itprint
    character(len = 128) :: fileout
    integer              :: unitout, i_mp, i_active

#ifdef MPIP
    write(fileout,'(a,i3.3,a,i3.3,a)')  &
        trim(outputpath)//'BIN/pic',rank,'.',itprint,'.bin'
    unitout=rank+10
#else
    write(fileout,'(a,i3.3,a)')  trim(outputpath)//'BIN/pic',itprint,'.bin'
    unitout=10
#endif

  open(unit=unitout,file=fileout,status='unknown',access='stream')

  i_active = 0
  do i_mp=1,N_MP
    if(partID(i_mp)/=0) i_active = i_active + 1
  end do

  !  write how many particles are in rank domain
  write(unitout) np, N_MP, i_active !n_activeMP


  !  loop over particles owned by processor and write the active ones
  do i_mp=1,N_MP
    !if (partOwner(i_mp)==rank) then
    if (partID(i_mp) /=0) then
      write(unitout) partID(i_mp)
      write(unitout) Q_MP0(i_mp,1:6)

    endif
  end do

end subroutine write_pic

end module pic_module

!================================================================