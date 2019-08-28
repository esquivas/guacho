!=======================================================================
!> @file pic_module.f90
!> @brief PIC module
!> @author Alejandro Esquivel & Matias Schneiter
!> @date 7/Ago/2019
!
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

!> @brief PIC module
!> @details Implementation of a prticle module, based in
!> Vaidya et al. 2918, ApJ, 865, 144

module pic_module

  implicit none

contains

  !================================================================
  !> @brief Initialization of module
  !> @details Allocates memory for all global variables that correspond
  !> to the particle module
  subroutine init_pic()
    use parameters, only : nx, ny, nz, pic_distF, N_MP, NBinsSEDMP
    use globals, only : Q_MP0, Q_MP1, MP_SED, P_DSA, shockF, partID, partOwner
    implicit none

    if(pic_distF) then
      allocate( Q_MP0(N_MP,12) )
      !Q_MP0(i, eq) has the following info:
      ! eq = 1-3 : x, y, z
      ! eq = 4-6 : vx, vy, vz
      ! eq = 7   : b**2 = bx**2+by**2+bz**2
      ! eq = 8   : rho
      ! eq = 9   : P
      ! eq = 10  : shock flag (1 if shocked)
      ! eq = 11  : compression ratio (does not reset)
      ! eq = 12  : angle between the shock normal and the preshock field

      allocate( shockF(nx,ny,nz) )
      !  used to mark in the MHD grid shocked regions (shockF(i,j,k)=1)

      allocate( MP_SED(2,NBinsSEDMP,N_MP) )
      !MP_SED(1,:,i) :  Energy (Lagrangian) bins
      !MP_SED(1,:,i) :  Number of MP with Energy E_i +- Delta E

      allocate( P_DSA(N_MP,2,8))
      !P_DSA(i, 1, :) : Pre  shock MHD info (U1 in Vaidya et al 2018)
      !P_DSA(i, 2, :) : Post shock MHD info (U2 in Vaidya et al 2018)
    else
      allocate( Q_MP0(N_MP,6) )
      !Q_MP0(i, eq) has the following info:
      ! eq = 1-3 : x, y, z
      ! eq = 4-6 : vx, vy, vz
    end if

    allocate( Q_MP1(N_MP,3) )    ! x,y,z position advanced by the predictor
    allocate( partID   (N_MP) )  ! Individual particle identifier
    allocate( partOwner(N_MP) )  ! Rank of the processor that owns said particle

  end subroutine init_pic

  !================================================================
  !> @brief Deactivation of particle
  !> @details Takes out particle from the active list, and modifies then
  !> updates the n_activeMP variable if needed
  !> @param integer [in] i_mp : local position of the particle to be deactivated
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
  !> @brief Insertion of of new particle
  !> @details Add a new particle, and its data to the local arrays
  !> partID(i_mp), Q_MP0(i_mp)
  !> updates the n_activeMP variable if needed
  !> @param integer [in ] ID    : ID of particle to be inserted
  !> @param integer [in ] ndata : number of data fields to be inserted
  !> @param integer [in ] Qdata(ndata) : Data to be loaded into Q_MP0
  !> @param integer [out] i_mp  : local position of the particle added
  subroutine addMP(ID, ndata, Qdata, i_mp)
    use parameters, only : N_MP
    use globals,    only : partID,Q_MP0, n_activeMP
    implicit none
    integer, intent(in)  :: ID,ndata
    real,    intent(in)  :: Qdata(ndata)
    integer, intent(out) :: i_mp

    ! sweep list and add element in empty slot
    do i_mp=1, n_mp
      if(partID(i_mp) == 0) then
        partID(i_mp)  = ID
        Q_MP0(i_mp,1:ndata) = Qdata(1:ndata)
        if (i_mp > n_activeMP) n_activeMP = i_mp
        return
      end if
    end do

  end subroutine addMP

  !================================================================
  !> @brief Predictor step subroutine
  !> @details Advances the position of the particle for the predictor step
  !> as described in Vaidya et. al (2018)
  !> It also implements the required update of the SED of each MP, including
  !> the Diffuse Shock Acceleration treatment
  subroutine PICpredictor()

    use globals,   only : primit, dt_CFL, rank, comm3d, &
                          Q_MP0, Q_MP1, P_DSA, MP_SED, partID, currentIteration
    use parameters
    use constants, only : pi
    use utilities, only : isInDomain, inWhichDomain, isInShock
    implicit none
    integer :: i_mp, i, j, k, l, ind(3), dest, nLocSend, &
               sendLoc(0:np-1), sendList(0:np-1,0:np-1), &
               dataLoc(N_MP), iS, iR
    real    :: weights(8)
    real    :: fullSend(2*NBinsSEDMP+28), fullRecv(2*NBinsSEDMP+28)
    !          that is 2*NBinsSEDMP of the SED, 12 of Q_MP0 and 2*8 from P_DSA
    integer :: status(MPI_STATUS_SIZE), err
    real    :: normal(3), comp, thB1, thB2
    ! initialize send and recv lists
    dataLoc(:)    =  0
    sendLoc(:)    =  0
    sendlist(:,:) =  0
    nLocSend      =  0

    do i_mp=1, n_MP
      ! execute only if particle i_mp is in the active list
      if (partID(i_mp)/=0) then

        ! do only if paricle is in domain
        if ( isInDomain( Q_MP0(i_mp,1:3) ) ) then

          ! Calculate interpolation reference and weights
          call interpBD(Q_MP0(i_mp,1:3),ind,weights)

          !  Clear arrays for interpolation
          if (pic_distF) then
            Q_MP0(i_mp,4:9) = 0.
          else
            Q_MP0(i_mp,4:6) = 0.
          end if

          !  Interpolate u, [B^2 & rho if needed] to particle position
          l = 1
          do k= ind(3),ind(3)+1
            do j=ind(2),ind(2)+1
              do i=ind(1),ind(1)+1
                !   interpolate velocity to the particle position
                Q_MP0(i_mp,4) = Q_MP0(i_mp,4) + primit(2,i,j,k) * weights(l)
                Q_MP0(i_mp,5) = Q_MP0(i_mp,5) + primit(3,i,j,k) * weights(l)
                Q_MP0(i_mp,6) = Q_MP0(i_mp,6) + primit(4,i,j,k) * weights(l)

                !  if we're following the MP SED, interpolate alpha and beta
                if(pic_distF) then
                  !  computed following Vaidya et al, 2018 ApJ
                  !  B**2
                  Q_MP0(i_mp,7) = Q_MP0(i_mp,7) + weights(l)**2 *              &
                                                (  primit(6,i,j,k)**2 +        &
                                                   primit(7,i,j,k)**2 +        &
                                                   primit(8,i,j,k)**2  )
                  !  aiabatic expansion term rho^n
                  Q_MP0(i_mp,8) = Q_MP0(i_mp,8) + primit(1,i,j,k)*weights(l)
                end if

                l = l + 1
              end do
            end do
          end do

          !   DSA calculation
          if (pic_distF) then

          !   If particle was already inside shock
            if (Q_MP0(i_mp,10) /= 0.) then

              !print*, 'particle ', partID(i_mp),                                 &
              !        'marked inside the shock region', currentIteration
              !  interpolate Pressure
              !Q_MP0(i_mp,9) = 0.
              l = 1
              do k= ind(3),ind(3)+1
                do j=ind(2),ind(2)+1
                  do i=ind(1),ind(1)+1
                    Q_MP0(i_mp,9) = Q_MP0(i_mp,9) + primit(5,i,j,k)*weights(l)
                    l = l + 1
                  end do
                end do
              end do

              if(Q_MP0(i_mp,9) <= P_DSA(i_mp,1,5)) then  ! (P < Pmin)
                ! interpolate primitives and load them to P_DSA(i_mp,1,:)
                !print*, 'Set P1 of particle ', partID(i_mp)
                P_DSA(i_mp,1,1) = Q_MP0(i_mp,8)  !  Density
                P_DSA(i_mp,1,2) = Q_MP0(i_mp,4)  !  vx
                P_DSA(i_mp,1,3) = Q_MP0(i_mp,5)  !  vy
                P_DSA(i_mp,1,4) = Q_MP0(i_mp,6)  !  vz
                P_DSA(i_mp,1,5) = Q_MP0(i_mp,9)  !  P
                P_DSA(i_mp,1,6:8) = 0.
                l = 1
                do k= ind(3),ind(3)+1
                  do j=ind(2),ind(2)+1
                    do i=ind(1),ind(1)+1
                      !  Interpolate B field components
                      P_DSA(i_mp,1,6)=P_DSA(i_mp,1,6)+primit(6,i,j,k)*weights(l)
                      P_DSA(i_mp,1,7)=P_DSA(i_mp,1,7)+primit(7,i,j,k)*weights(l)
                      P_DSA(i_mp,1,8)=P_DSA(i_mp,1,8)+primit(8,i,j,k)*weights(l)
                      l = l + 1
                    end do
                  end do
                end do

              else if (Q_MP0(i_mp,9) >= P_DSA(i_mp,2,5)) then  ! (P > Pmax)
                ! interpolate primitives and load them to P_DSA(i_mp,2,:)
                !print*, 'Set P2 of particle ', partID(i_mp)
                P_DSA(i_mp,2,1) = Q_MP0(i_mp,8)  !  Density
                P_DSA(i_mp,2,2) = Q_MP0(i_mp,4)  !  vx
                P_DSA(i_mp,2,3) = Q_MP0(i_mp,5)  !  vy
                P_DSA(i_mp,2,4) = Q_MP0(i_mp,6)  !  vz
                P_DSA(i_mp,2,5) = Q_MP0(i_mp,9)  !  P            .
                P_DSA(i_mp,2,6:8) = 0.
                l = 1
                do k= ind(3),ind(3)+1
                  do j=ind(2),ind(2)+1
                    do i=ind(1),ind(1)+1
                      !  Interpolate B field components
                      P_DSA(i_mp,2,6)=P_DSA(i_mp,2,6)+primit(6,i,j,k)*weights(l)
                      P_DSA(i_mp,2,7)=P_DSA(i_mp,2,7)+primit(7,i,j,k)*weights(l)
                      P_DSA(i_mp,2,8)=P_DSA(i_mp,2,8)+primit(8,i,j,k)*weights(l)
                      l = l + 1
                    end do
                  end do
                end do

              end if

              if (.not.isInShock(Q_MP0(i_mp,1:3))) then
                !print*, 'particle ', partID(i_mp),                               &
                !        'has left the shock region', currentIteration

                !***  Here, the SED should be updated with te DSA prescription***
                call get_NRth(P_DSA(i_mp,1,:),P_DSA(i_mp,2,:), normal, comp,   &
                              thB1,thB2)

                Q_MP0(i_mp,11) = max(comp, Q_MP0(i_mp,11))
                Q_MP0(i_mp,12) = thB2*180./pi

                !   Mark particle with -1 as "just left shock region"
                Q_MP0(i_mp,10) = -1.

              end if

            else
              !  particle not marked as in shock, it is either entering or
              !  happily living its life
              if (isInShock(Q_MP0(i_mp,1:3))) then

                !print*, 'particle ', partID(i_mp),                               &
                !        ' has just entered shock', currentIteration

                !  Mark it as shocked for future Reference
                Q_MP0(i_mp,10) = 1.

                !  interpolate primitives and load them to both P_DSA(i_mp,1:2,:)
                P_DSA(i_mp,1,1) = Q_MP0(i_mp,8)  !  Density
                P_DSA(i_mp,1,2) = Q_MP0(i_mp,4)  !  vx
                P_DSA(i_mp,1,3) = Q_MP0(i_mp,5)  !  vy
                P_DSA(i_mp,1,4) = Q_MP0(i_mp,6)  !  vz
                P_DSA(i_mp,1,6:8) = 0.
                Q_MP0(i_mp,9) = 0.
                l = 1
                do k= ind(3),ind(3)+1
                  do j=ind(2),ind(2)+1
                    do i=ind(1),ind(1)+1
                      !  Interpolate Pressure
                      Q_MP0(i_mp,9) = Q_MP0(i_mp,9) + primit(5,i,j,k)*weights(l)
                      !  Interpolate B field components
                      P_DSA(i_mp,1,6)=P_DSA(i_mp,1,6)+primit(6,i,j,k)*weights(l)
                      P_DSA(i_mp,1,7)=P_DSA(i_mp,1,7)+primit(7,i,j,k)*weights(l)
                      P_DSA(i_mp,1,8)=P_DSA(i_mp,1,8)+primit(8,i,j,k)*weights(l)
                      l = l + 1
                    end do
                  end do
                end do
                P_DSA(i_mp,1,5) = Q_MP0(i_mp,9)  !  P

                !  Copy to P_DSA(i_mp,2,:)
                P_DSA(i_mp,2,:) = P_DSA(i_mp,1,:)

              endif

            end if

          end if

          !  predictor step
          Q_MP1(i_mp,1:3) = Q_MP0(i_mp,1:3)+ dt_CFL*Q_MP0(i_mp,4:6)

        end if  ! isIndomain (l 158)

        !  check if particle is leaving the domain
        dest = inWhichDomain(Q_MP1(i_mp,1:3))
        if (dest == -1) call deactivateMP(i_mp)
        if( dest /= rank .and. dest /= -1 ) then
          !  count for MPI exchange
          sendLoc(dest) = sendLoc(dest) + 1
          nLocSend      = nLocSend      + 1
          dataLoc(nLocSend) = i_mp
          !print'(i2,a,i4,a,i4)', rank, ' *** particle ',partID(i_mp),  &
          !                     ' is going to ',DEST
        end if

      end if  ! if part was active
    end do    ! loop over all particles

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
              !    if (iR /= -1) then
              if (pic_distF) then

                !print*,'>>>',rank,partID(dataLoc(i)),dataLoc(i)
                !  pack info if we are solving the SED
                fullSend( 1:12) = Q_MP0(dataLoc(i),1:12)
                fullSend(13:20) = P_DSA(dataLoc(i),1,:)
                fullSend(21:28) = P_DSA(dataLoc(i),2,:)
                fullSend(29:           28  +NBinsSEDMP)=                       &
                                              MP_SED(1,1:NBinsSEDMP,dataLoc(i) )
                fullSend(29+NBinsSEDMP:28+2*NBinsSEDMP)=                       &
                                              MP_SED(2,1:NBinsSEDMP,dataLoc(i) )
                !  send the whole thing
                call mpi_send( fullSend ,28+2*NBinsSEDMP, mpi_real_kind ,IR,   &
                               partID(dataLoc(i)), comm3d,err )

              else

                !  in case we are only passing the particles and not their SED
                call mpi_send( Q_MP0(dataLoc(i),1:6) , 6, mpi_real_kind ,IR,   &
                               partID(dataLoc(i)), comm3d,err )

              endif

              !  deactivate particle from current processor
              call deactivateMP(dataLoc(i))

            end do

          end if

          if(iR == rank) then

            !print'(i0,a,i0,a,i0)', rank,'<--', IS, ':',sendList(iR,iS)
            do i=1,sendList(iR,iS)

              if (pic_distF) then

                !print*,'<<<',rank,' will recv here from ',IS
                !  in case we are only the particles w/their SED
                call mpi_recv(fullRecv,28+2*NBinsSEDMP, mpi_real_kind, IS,     &
                              mpi_any_tag,comm3d, status, err)

                !  add current particle in list and data in new processor
                call addMP( status(MPI_TAG), 12, fullRecv(1:12), i_mp )
                P_DSA(i_mp,1,1:8) = fullRecv(13:20)
                P_DSA(i_mp,2,1:8) = fullRecv(21:28)
                ! unpack the SED
                MP_SED(1,1:100,i_mp) = fullRecv(29           :28  +NBinsSEDMP)
                MP_SED(2,1:100,i_mp) = fullRecv(29+NBinsSEDMP:28+2*NBinsSEDMP)

              else

                !  in case we are only passing the particles and not their SED
                call mpi_recv(fullRecv(1:6), 6, mpi_real_kind, IS, mpi_any_tag,&
                              comm3d, status, err)
                !print*,'received successfuly', status(MPI_TAG)

                !  add current particle in list and data in new processor
                call addMP( status(MPI_TAG), 6, fullRecv(1:6), i_mp )

              end if

              !  recalculate predictor step for newcomer
              Q_MP1(i_mp,1:3) = Q_MP0(i_mp,1:3)+ dt_CFL*Q_MP0(i_mp,4:6)

            end do

          end if

        end if
      end do
    end do

  end subroutine PICpredictor

  !================================================================
  !> @brief Corrector step subroutine
  !> @details Advances the position of the particle for the corrector step
  !> as described in Vaidya et. al (2018)
  !> It also implements the required update of the SED of each MP, including
  !> the Diffuse Shock Acceleration treatment
  subroutine PICcorrector()

    use globals,   only : primit, dt_CFL, rank, comm3d,                        &
                          MP_SED, Q_MP0, Q_MP1, partID, P_DSA
    use parameters
    use utilities, only : inWhichDomain, isInDomain, isInShock
    implicit none
    integer :: i_mp, i, j, k, l, ind(3)
    real    :: weights(8), vel1(3)
    integer :: dest, nLocSend, dataLoc(N_mp),  sendLoc(0:np-1),                &
               sendList(0:np-1,0:np-1)
    real    :: fullSend(2*NBinsSEDMP+5), fullRecv(2*NBinsSEDMP+5)
    integer :: status(MPI_STATUS_SIZE), err, iR, iS, ib
    real    :: ema, bdist, rhoNP1, crNP1, dataIn(12)
    ! initialize send and recv lists
    dataLoc(:)    =  0
    sendLoc(:)    =  0
    sendlist(:,:) =  0
    nLocSend      =  0

    do i_mp=1, n_MP
      ! execute only if particle i_mp is in the active list
        if (partID(i_mp)/=0 .and. isInDomain(Q_MP1(i_mp,1:3)) ) then

          !  clear come variables
          if(pic_distF) then
            ema = 0.
            bdist = 0.
            rhoNP1 = 0.
            crNP1  = 0.
          end if

          !  Interpolate the velocity field to particle position, and add
          !  to the velocity from the corrector step
          !  Calculate interpolation reference and weights
          call interpBD(Q_MP1(i_mp,1:3),ind,weights)
          !  Interpolates the magnetic field to calculate beta.
          vel1(:) = 0.
          l = 1
          do k=ind(3),ind(3)+1
            do j=ind(2),ind(2)+1
              do i=ind(1),ind(1)+1
                !  interpolate velocity / density and B**2
                vel1(1) = vel1(1) + primit(2,i,j,k)*weights(l)
                vel1(2) = vel1(2) + primit(3,i,j,k)*weights(l)
                vel1(3) = vel1(3) + primit(4,i,j,k)*weights(l)

                !  if we're following the MP SED, interpolate alpha and beta
                if (pic_distF) then
                  !   source terms
                  rhoNP1 = rhoNP1 + primit(1,i,j,k)*weights(l)
                  crNP1  = crNP1  + weights(l)**2 *                            &
                   ( primit(6,i,j,k)**2 +primit(7,i,j,k)**2+primit(8,i,j,k)**2 )

                end if
                l = l + 1
              end do
            end do
          end do

          !  corrector step (only position is updated)
          Q_MP0(i_mp,1:3) = Q_MP0(i_mp,1:3)                                    &
                          + 0.5*dt_CFL*( Q_MP0(i_mp,4:6) + vel1(1:3) )

          if (pic_distF) then
            !  exp(-a) and cr
            ema    = (rhoNP1/Q_MP0(i_mp,8))**(1./3.)
            ! eq. (23) Vaidya et al. 2018
            bdist  = 0.5*dt_CFL*( Q_MP0(i_mp,7) + ema*crNP1 )

            !  update only if *not* currently marked as inside shock
            if (Q_MP0(i_mp,10) == 0.) then
              do ib=1,NBinsSEDMP

                MP_SED(2,ib,i_mp)=MP_SED(2,ib,i_mp)*ema*                       &
                                  (1.+bdist*MP_SED(1,ib,i_mp))**2

                MP_SED(1,ib,i_mp)=MP_SED(1,ib,i_mp)*ema/                       &
                                  (1.+bdist*MP_SED(1,ib,i_mp))

              end do
            else if (Q_MP0(i_mp,10) == -1.) then  ! inject spectra after shock

              call inject_spectrum(i_mp)

              !  Clear primit P1/P2 arrays and mark as no longer in shock
              P_DSA(i_mp,:,:)= 0.
              Q_MP0(i_mp,10) = 0.

            end if

        end if

        !  check if particle is leaving the domain
        dest = inWhichDomain( Q_MP0(i_mp,1:3) )
        if (dest == -1) call deactivateMP(i_mp)
        if (dest /= rank .and. dest /= -1) then
          !  count for MPI exchange
          sendLoc(dest) = sendLoc(dest) + 1
          nLocSend      = nLocSend      + 1
          dataLoc(nLocSend) = i_mp
        end if

      end if
    end do

    !   consolidate list to have info of all send/receive operations
    call mpi_allgather(sendLoc(:),  np, mpi_integer,                           &
                       sendList, np, mpi_integer, comm3d,err)

    !  exchange particles
    do iR=0,np-1
      do iS=0,np-1
        if(sendList(iR,iS) /= 0) then

          if(iS == rank) then
            do i=1,sendlist(iR,iS)

              !      if (iR /= -1) then
              if (pic_distF) then

                !  pack info if we are solving the SED
                fullSend(1:3) = Q_MP0(dataLoc(i),1:3)
                fullSend(4:5) = Q_MP0(dataLoc(i),4:5)
                fullSend(6:             5+NBinsSEDMP)=                         &
                                              MP_SED(1,1:NBinsSEDMP,dataLoc(i) )
                fullSend(6+NBinsSEDMP:5+2*NBinsSEDMP)=                         &
                                              MP_SED(2,1:NBinsSEDMP,dataLoc(i) )

                  !  send the whole thing
                call mpi_send( fullSend ,5+2*NBinsSEDMP, mpi_real_kind ,IR,    &
                               partID(dataLoc(i)), comm3d,err )

              else

                !   if not solving the SED, send only X
                call mpi_send( Q_MP0(dataLoc(i),1:3) , 3, mpi_real_kind ,IR,   &
                               partID(dataLoc(i)), comm3d,err)

              end if

              !  deactivate particle from current processor
              call deactivateMP(dataLoc(i))

            end do
          end if

          if(iR == rank) then
            do i=1,sendList(iR,iS)

              ! Sets de distribution fuction if enabled
              !equation 3 in Vaidya et al. 2016.
              if (pic_distF) then
                !  in case we are only the particles w/their SED
                call mpi_recv(fullRecv,5+2*NBinsSEDMP, mpi_real_kind, IS,      &
                                       mpi_any_tag,comm3d, status, err)

                 !  add current particle in list and data in new processor
                 dataIn(1:3) = fullRecv(1:3)
                 dataIn(4:10) = 0.
                 dataIn(11:12) = fullRecv(4:5)
                 call addMP( status(MPI_TAG), 3, dataIn(1:12), i_mp )

                 ! unpack the SED
                 MP_SED(1,1:100,i_mp) = fullRecv(6:             5+NBinsSEDMP)
                 MP_SED(2,1:100,i_mp) = fullRecv(6+NBinsSEDMP:5+2*NBinsSEDMP)
              else

                call mpi_recv( fullRecv(1:3), 3, mpi_real_kind, IS,            &
                               mpi_any_tag, comm3d, status, err)
                !  add current particle in list and data in new processor
                call addMP( status(MPI_TAG), 3, fullRecv(1:3), i_mp )

              end if

            end do
          end if

        end if
      end do
    end do

  end subroutine PICcorrector

  !================================================================
  !> @brief bilineal interpolator
  !> @ details: returns weights and indices corresponding to a bilineal
  !> interpolation at for the particle position.
  !> @param real [in] pos(3) : Three dimensional position of particle,
  !> @param integer   ind(3) : Reference corner i0,j0,k0 indices of the
  !> cube of cells that are used in the interpolation
  !> @param real [out] weights(8) : weights associated with each corner
  !> of the interopolation in the following order
  !> 1-> i0,j0,k0,  3-> i0,j1,k0,  5-> i0,j0,k1,  7-> i0,j1,k1
  !> 2-> i1,j0,k0,  4-> i1,j1,k0,  6-> i1,j0,k1,  8-> i1,j1,k1
  subroutine interpBD(pos,ind,weights)

    use parameters, only : nx, ny, nz
    use globals,    only : coords, dx, dy, dz
    implicit none
    real,    intent(in)  :: pos(3)
    integer, intent(out) :: ind(3)
    real,    intent(out) :: weights(8)
    real                 :: x0, y0, z0, remx, remy, remz, distx, disty, distz

    ! get the index in the whole domain (integer part)
    ! this is the particle position in "cell units" from 1 to nxmax (real)
    x0 = pos(1)/dx
    y0 = pos(2)/dy
    z0 = pos(3)/dz

    !  get reminder
    remx = x0 - int(x0)
    remy = y0 - int(y0)
    remz = z0 - int(z0)

    !  get bounds
    if (remx < 0.5) then
      ind(1)= int(x0) - coords(0)*nx
    else
      ind(1)= int(x0) - coords(0)*nx + 1
    end if

    if (remy < 0.5) then
      ind(2)= int(y0)- coords(1)*ny
    else
      ind(2)= int(y0)- coords(1)*ny + 1
    end if

    if (remz < 0.5) then
      ind(3)= int(z0)- coords(2)*nz
    else
      ind(3)= int(z0)- coords(2)*nz + 1
    end if

    !  get distances
    distx = x0 - ( real(ind(1) + nx*coords(0) ) -0.5 )
    disty = y0 - ( real(ind(2) + ny*coords(1) ) -0.5 )
    distz = z0 - ( real(ind(3) + nz*coords(2) ) -0.5 )

    weights(1) = (1.-distx) * (1.-disty) * (1.-distz)
    weights(2) =     distx  * (1.-disty) * (1.-distz)
    weights(3) = (1.-distx) *     disty  * (1.-distz)
    weights(4) =     distx  *     disty  * (1.-distz)
    weights(5) = (1.-distx) * (1.-disty) *     distz
    weights(6) =     distx  * (1.-disty) *     distz
    weights(7) = (1.-distx) *     disty  *     distz
    weights(8) =     distx  *     disty  *     distz

    return

  end subroutine interpBD

  !================================================================
  !> @brief Writes pic output
  !> @details Writes position and particles velocities in a single file
  !> format is binary, one integer with the number of points in domain,
  !> and then all the positions and velocities (in the default real precision)
  !> param integer [in] itprint : number of itration (coincides with the [M]HD
  !> output)
  subroutine write_pic(itprint)

    use parameters, only : outputpath, np, pic_distF, N_MP, NBinsSEDMP
    use utilities
    use globals,    only : rank, Q_MP0, MP_SED, partID, P_DSA
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
    if (.not.pic_distF) then
      write(unitout) np, N_MP, i_active, 0  !n_activeMP
    else
      write(unitout) np, N_MP, i_active, NBinsSEDMP
    end if

    !  loop over particles owned by processor and write the active ones
    do i_mp=1,N_MP
      if (partID(i_mp) /=0) then

        write(unitout) partID(i_mp)
        write(unitout) Q_MP0(i_mp,1:3)
        if(pic_distF) then
          write(unitout) Q_MP0(i_mp,11:12)
          write(unitout) MP_SED(1:2,:,i_mp)
          write(unitout) P_DSA(i_mp,:,:)
        end if

      end if
    end do

    close(unitout)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  tHIS IS ONLY FOR DEBUGGING PURPOSES
    !  repeat only for shocked PARTICLES
    write(fileout,'(a,i3.3,a,i3.3,a)')  &
    trim(outputpath)//'BIN/pic-shocked-',rank,'.',itprint,'.bin'
    unitout=rank+10

    open(unit=unitout,file=fileout,status='unknown',access='stream')

    i_active = 0
    do i_mp=1,N_MP
      if (partID(i_mp)/=0) then
        if(isInShock(Q_MP0(i_mp,1:3)) ) i_active = i_active + 1
      endif
    end do
    !print*, i_active, ' particles in shock'

    write(unitout) np, N_MP, i_active, 0
    do i_mp=1,N_MP
      if (partID(i_mp) /=0) then
        if (isInShock(Q_MP0(i_mp,1:3))) then
          write(unitout) Q_MP0(i_mp,1:3)
        end if
      end if
    end do

    close(unitout)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine write_pic

  !=======================================================================
  !> @brief Get n and r
  !> @details Compute the shock normal and conpression ratio
  !> param real [in]  prim1(8) : primitives in pre-shock region
  !> param real [in]  prim2(8) : primitives in post-shock region
  !> param real [out] n(3)     : unitary vector normal to the shock
  !> param real [out] r        : compression ratio (rho2/rho1)
  !> param real [out] thB1     : angle between shock normal and B1
  !> param real [out] thB2     : angle between shock normal and B2
  subroutine get_NRth(prim1,prim2,nsh,r,thB1,thb2)
    !see equation 26 to 27 Vaidya 2018
    implicit none
    real, intent(in) :: prim1(8), prim2(8)
    real, intent(out) :: nsh(3), r, thB1, thB2
    real :: delV(3), delB(3), BdotN, magN, magB

    !  the compression ratio
    r = prim2(1)/prim1(1)

    delB(1)=prim2(6)-prim1(6)  ! deltaBx
    delB(2)=prim2(7)-prim1(7)  ! deltaBy
    delB(3)=prim2(8)-prim1(8)  ! deltaBz

    delV(1)=prim2(2)-prim1(2)  ! deltaVx
    delV(2)=prim2(3)-prim1(3)  ! deltaVy
    delV(3)=prim2(4)-prim1(4)  ! deltaVz

    !  get normal vector for th /= 0 or 90  (magnetic coplanarity)
    nsh(1) = (prim1(8)*delV(1)-prim1(6)*delV(3))*delB(3)                       &
           - (prim1(6)*delV(2)-prim1(7)*delV(1))*delB(2)

    nsh(2) = (prim1(6)*delV(2)-prim1(7)*delV(1))*delB(1)                       &
           - (prim1(7)*delV(3)-prim1(8)*delV(2))*delB(3)

    nsh(3) = (prim1(7)*delV(3)-prim1(8)*delV(2))*delB(2)                       &
           - (prim1(8)*delV(1)-prim1(6)*delV(3))*delB(1)

    magN = sqrt(nsh(1)**2 + nsh(2)**2 + nsh(3)**2)
    if (magN == 0) then
      nsh(:) = 0.
    else
      nsh(:) = nsh(:)/magN
    end if

    !  Compute thB1
    magB = sqrt( prim1(6)**2 + prim1(7)**2 + prim1(8)**2 )
    if (magB /= 0.) then
      BdotN = abs( nsh(1)*prim1(6) + nsh(2)*prim1(7) + nsh(3)*prim1(8) ) / magB
    else
      BdotN =  0.
    end if
    !  Recompute normal if thB1 will be close to 0 or 90 degrees
    if ( (BdotN > 0.996 ).or.(BdotN < 0.087 ) )   then
      nsh(1:3) = delV(1:3)
      magN = sqrt(delV(1)**2 + delV(2)**2 + delV(3)**2)
      if (magN  == 0.) then
        nsh(:) = 0.
      else
        nsh(1:3)=nsh(1:3)/magN
      end if
      magB = sqrt( prim1(6)**2 + prim1(7)**2 + prim1(8)**2 )
      if (magB /= 0.) then
        BdotN = abs( nsh(1)*prim1(6) + nsh(2)*prim1(7) + nsh(3)*prim1(8) ) /magB
      else
        BdotN =  0.
      end if
    end if
    BdotN = min(BdotN,1.0)
    BdotN = max(BdotN,0.0)
    thB1 = acos(BdotN)

    !  Compute thB2
    magB = sqrt( prim2(6)**2 + prim2(7)**2 + prim2(8)**2 )
    if (magB /= 0.) then
      BdotN = abs( nsh(1)*prim2(6) + nsh(2)*prim2(7) + nsh(3)*prim2(8) ) /magB
    else
      BdotN =  0.
    end if
    !  Recompute normal if thB1 will be close to 0 or 90 degrees
    if ( (BdotN > 0.996 ).or.(BdotN < 0.087 ) )   then
    !if ( (BdotN > 0.98 ).or.(BdotN < 0.16 ) )   then
      nsh(1:3) = delV(1:3)
      magN = sqrt(delV(1)**2 + delV(2)**2 + delV(3)**2)
      if (magN  == 0.) then
        nsh(:) = 0.
      else
        nsh(1:3)=nsh(1:3)/magN
      end if
      magB = sqrt( prim2(6)**2 + prim2(7)**2 + prim2(8)**2 )
      if (magB /= 0.) then
        BdotN = abs( nsh(1)*prim2(6) + nsh(2)*prim2(7) + nsh(3)*prim2(8) ) /magB
      else
        BdotN =  0.
      end if
    end if
    BdotN = min(BdotN,1.0)
    BdotN = max(BdotN,0.0)
    thB2  = acos(BdotN)

  end subroutine get_NRth

  !=======================================================================
  !> @brief Inject inject_spectrum
  !> @details Inject new spectrum, after DSA subgrid calculation
  !> @param integer [in] i_mp : index of the MP to which the SED is updated
  subroutine inject_spectrum(i_mp)
    use globals,    only : MP_SED, Q_MP0, P_DSA
    use parameters, only : NBinsSEDMP
    implicit none
    integer, intent(in) :: i_mp
    integer :: i
    real    :: q_index, r, E0

    r = max(1.5,Q_MP0(i_mp,11))

    q_index = 3.*r/(r-1.)
    E0 = MP_SED(2,1,i_mp)*MP_SED(1,1,i_mp)

    do i = 1,NBinsSEDMP
      MP_SED(2,i,i_mp)= MP_SED(1,i,i_mp)**(-q_index)/E0
    end do

  end subroutine inject_spectrum

end module pic_module

!================================================================
