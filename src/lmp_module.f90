!=======================================================================
!> @file lmp_module.f90
!> @brief LMP module
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

!> @brief LMP module
!> @details Implementation of a particle module, based in the Lagrangian tracer
!> particles in Vaidya et al. 2918, ApJ, 865, 144

module lmp_module

  implicit none

contains

  !================================================================
  !> @brief Initialization of module
  !> @details Allocates memory for all global variables that correspond
  !> to the particle module
  subroutine init_lmp()
    use parameters, only : nx, ny, nz, lmp_distf, N_MP, NBinsSEDMP
    use globals, only : Q_MP0, Q_MP1, MP_SED, P_DSA, shockF, partID, partOwner
    implicit none

    if(lmp_distf) then
      allocate( Q_MP0(N_MP,12) )
      ! Q_MP0(i, eq) has the following info:
      ! eq = 1-3 : x, y, z
      ! eq = 4-6 : vx, vy, vz
      ! eq = 7   : b**2/2 = (bx**2+by**2+bz**2)/2
      ! eq = 8   : rho
      ! eq = 9   : P
      ! eq = 10  : shock flag (1 if shocked)
      ! eq = 11  : compression ratio (does not reset)
      ! eq = 12  : angle between the shock normal and the preshock field

      allocate( shockF(nx,ny,nz) )
      !  used to mark in the MHD grid shocked regions (shockF(i,j,k)=1)

      allocate( MP_SED(2,NBinsSEDMP,N_MP) )
      !  MP_SED(1,:,i) :  Energy (Lagrangian) bins
      !  MP_SED(2,:,i) :  Number of MP with Energy E_i +- Delta E

      allocate( P_DSA(N_MP,2,8))
      !  P_DSA(i, 1, :) : Pre  shock MHD info (U1 in Vaidya et al 2018)
      !  P_DSA(i, 2, :) : Post shock MHD info (U2 in Vaidya et al 2018)
    else
      allocate( Q_MP0(N_MP,6) )
      !  Q_MP0(i, eq) has the following info:
      !  eq = 1-3 : x, y, z
      !  eq = 4-6 : vx, vy, vz
    end if

    allocate( Q_MP1(N_MP,3) )    ! x,y,z position advanced by the predictor
    allocate( partID   (N_MP) )  ! Individual particle identifier
    allocate( partOwner(N_MP) )  ! Rank of the processor that owns said particle

  end subroutine init_lmp

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

    ! loop over list and add element in empty slot
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
  subroutine LMPpredictor()
    use globals,   only : primit, dt_CFL, rank, comm3d, &
                          Q_MP0, Q_MP1, P_DSA, MP_SED, partID
    use parameters
    use constants, only : pi
    use utilities, only : isInDomain, inWhichDomain, isInShock
    implicit none
    integer :: i_mp, i, j, k, l, ind(3)
    real    :: weights(8)
    integer :: dest, nLocSend, sendLoc(0:np-1), sendList(0:np-1,0:np-1),       &
               dataLoc(N_MP), iS, iR, status(MPI_STATUS_SIZE), err
    real    :: fullSend(2*NBinsSEDMP+28), fullRecv(2*NBinsSEDMP+28)
    !          above is 2*NBinsSEDMP of the SED, 12 of Q_MP0 and 2*8 from P_DSA
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
          if (lmp_distf) then
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
                if(lmp_distf) then
                  !  computed following Vaidya et al, 2018 ApJ
                  !  B**2
                  Q_MP0(i_mp,7) = Q_MP0(i_mp,7) + weights(l)**2 *              &
                                                (  primit(6,i,j,k)**2 +        &
                                                   primit(7,i,j,k)**2 +        &
                                                   primit(8,i,j,k)**2  ) /2.
                  !  aiabatic expansion term rho^n
                  Q_MP0(i_mp,8) = Q_MP0(i_mp,8) + primit(1,i,j,k)*weights(l)
                end if

                l = l + 1
              end do
            end do
          end do

          !   DSA calculation
          if (lmp_distf) then

          !   If particle was already inside shock
            if (Q_MP0(i_mp,10) /= 0.) then

              !print*, 'particle ', partID(i_mp),                              &
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
                !print*, 'particle ', partID(i_mp),                            &
                !        'has left the shock region', currentIteration

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

                !print*, 'particle ', partID(i_mp),                            &
                !        ' has just entered shock', currentIteration

                !  Mark it as shocked for future Reference
                Q_MP0(i_mp,10) = 1.

                ! interpolate primitives and load them to both P_DSA(i_mp,1:2,:)
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
          !print'(i2,a,i4,a,i4)', rank, ' *** particle ',partID(i_mp),        &
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
              if (lmp_distf) then

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

              if (lmp_distf) then

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

  end subroutine LMPpredictor

  !================================================================
  !> @brief Corrector step subroutine
  !> @details Advances the position of the particle for the corrector step
  !> as described in Vaidya et. al (2018)
  !> It also implements the required update of the SED of each MP, including
  !> the Diffuse Shock Acceleration treatment
  subroutine LMPcorrector()
    use globals,   only : primit, dt_CFL, rank, comm3d,                        &
                          MP_SED, Q_MP0, Q_MP1, partID, P_DSA
    use parameters
    use constants, only : sigma_SB, sigma_T,clight,emass
    use utilities, only : inWhichDomain, isInDomain, isInShock
    implicit none
    integer :: i_mp, i, j, k, l, ib, ind(3)
    real    :: weights(8)
    integer :: dest, nLocSend, dataLoc(N_mp),  sendLoc(0:np-1),                &
               sendList(0:np-1,0:np-1),status(MPI_STATUS_SIZE), err, iR, iS
    real    :: fullSend(2*NBinsSEDMP+5), fullRecv(2*NBinsSEDMP+5), dataIn(12)
    real    :: rhoNP1, vel1(3), pNP1, B_2Np1, BI, EI
    real    :: ema, bNP1
    real    :: q_NR, Emin, Emax, chi0
    !> RH term eq (7) Vaidya +
    real, parameter :: Tcmb = 2.278
    real, parameter :: Urad = sigma_SB*(Tcmb**4)/clight/Psc  !~1.05e-13
    !> constant in front of eq(7) in cgs
    real, parameter ::Cr0= ( 4.*sigma_T )/(3.* emass**2 * clight**3 )

    ! initialize send and recv lists
    dataLoc(:)    =  0
    sendLoc(:)    =  0
    sendlist(:,:) =  0
    nLocSend      =  0

    do i_mp=1, n_MP
      ! execute only if particle i_mp is in the active list
        if (partID(i_mp)/=0 .and. isInDomain(Q_MP1(i_mp,1:3)) ) then

          !  clear come variables
          if(lmp_distf) then
            ema    = 0.
            bNP1   = 0.
            rhoNP1 = 0.
            B_2NP1 = 0.
            pNP1   = 0.
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
                if (lmp_distf) then
                  !   source terms
                  rhoNP1 = rhoNP1 + primit(1,i,j,k)*weights(l)
                  B_2NP1  = B_2NP1  + 0.5 * weights(l)**2 *                    &
                   ( primit(6,i,j,k)**2 +primit(7,i,j,k)**2+primit(8,i,j,k)**2 )
                   pNP1  = pNP1   + primit(5,i,j,k)*weights(l)
                end if
                l = l + 1
              end do
            end do
          end do

          !  corrector step (only position is updated)
          Q_MP0(i_mp,1:3) = Q_MP0(i_mp,1:3)                                    &
                          + 0.5*dt_CFL*( Q_MP0(i_mp,4:6) + vel1(1:3) )

          if (lmp_distf) then
            !  exp(-a) and cr
            ema    = (rhoNP1/Q_MP0(i_mp,8))**(1./3.)
            ! eq. (23) Vaidya et al. 2018
            bNP1  = 0.5*dt_CFL*Cr0*( (Q_MP0(i_mp,7)+Urad) + ema*(B_2NP1+Urad) )
            !  convert to cgs to update SED
            bNP1  = bNP1 * tsc * Psc

            !  update only if *not* currently marked as inside shock
            if (Q_MP0(i_mp,10) == 0.) then
              do ib=1,NBinsSEDMP

                MP_SED(2,ib,i_mp)=MP_SED(2,ib,i_mp)*ema*                       &
                                  (1.+bNP1*MP_SED(1,ib,i_mp))**2

                MP_SED(1,ib,i_mp)=MP_SED(1,ib,i_mp)*ema/                       &
                                  (1.+bNP1*MP_SED(1,ib,i_mp))

              end do
            else if (Q_MP0(i_mp,10) == -1.) then  ! inject spectra after shock

              ! Call routine that calculates energy limits
              ! and power law parameters
              BI = 2. * B_2NP1
              EI = cv * pNP1

              call get_PL_parameters(i_mp,P_DSA(i_mp,1,:),P_DSA(i_mp,2,:)      &
                                     ,rhoNP1, EI, BI, chi0, q_NR, Emin, Emax)

              call inject_PL_spectrum(i_mp,chi0,q_NR,Emin,Emax)

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
              if (lmp_distf) then

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

              !  Sets de distribution fuction if enabled
              !  equation 3 in Vaidya et al. 2016.
              if (lmp_distf) then
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

  end subroutine LMPcorrector

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
  !> @brief Writes LMP output
  !> @details Writes position and particles velocities in a single file
  !> format is binary, one integer with the number of points in domain,
  !> and then all the positions and velocities (in the default real precision)
  !> param integer [in] itprint : number of itration (coincides with the [M]HD
  !> output)
  subroutine write_LMP(itprint)

    use parameters, only : outputpath, np, lmp_distf, N_MP, NBinsSEDMP,        &
                           rhosc, rsc, vsc2
    use utilities
    use globals,    only : rank, Q_MP0, MP_SED, partID, P_DSA
    implicit none
    integer, intent(in) :: itprint
    character(len = 128) :: fileout
    integer              :: unitout, i_mp, i_active

#ifdef MPIP
    write(fileout,'(a,i3.3,a,i3.3,a)')                                         &
        trim(outputpath)//'BIN/lmp',rank,'.',itprint,'.bin'
    unitout=rank+10
#else
    write(fileout,'(a,i3.3,a)')  trim(outputpath)//'BIN/lmp',itprint,'.bin'
    unitout=10
#endif

    open(unit=unitout,file=fileout,status='unknown',access='stream')

    i_active = 0
    do i_mp=1,N_MP
      if(partID(i_mp)/=0) then
         if( isInDomain(Q_MP0(i_mp,1:3)) ) i_active = i_active + 1
       end if
    end do

    !  write how many particles are in rank domain
    if (.not.lmp_distf) then
      write(unitout) np, N_MP, i_active, 0  !n_activeMP
    else
      write(unitout) np, N_MP, i_active, NBinsSEDMP
    end if

    !  loop over particles owned by processor and write the active ones
    do i_mp=1,N_MP
      if (partID(i_mp) /=0) then
        if( isInDomain(Q_MP0(i_mp,1:3)) ) then
          write(unitout) partID(i_mp)
          write(unitout) Q_MP0(i_mp,1:3)
          if(lmp_distf) then
            write(unitout) Q_MP0(i_mp,11:12)
            write(unitout) MP_SED(1,:,i_mp)
            write(unitout) MP_SED(2,:,i_mp)
            write(unitout) P_DSA(i_mp,:,:)
          end if
        endif
      end if
    end do

    close(unitout)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  tHIS IS ONLY FOR DEBUGGING PURPOSES
    !  repeat only for shocked PARTICLES
    write(fileout,'(a,i3.3,a,i3.3,a)')  &
    trim(outputpath)//'BIN/lmp-shocked-',rank,'.',itprint,'.bin'
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

  end subroutine write_LMP

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
  !> @details Inject new power law spectrum (for DSA subgrid calculation)
  !> of the form N \propto A0 E^(-m)
  !> @param integer [in] i_mp : index of the MP to which the SED is updated
  !> @param real    [in] chi0 : Amplitude
  !> @param real    [in] q    : spectral index
  !> @param real    [in] Emin : lower end energy in the spectrum
  !> @param real    [in] Emax : Upper end energy in the spectrun
  subroutine inject_PL_spectrum(i_mp, chi0, q, Emin, Emax)
    use globals,    only : MP_SED
    use parameters, only : NBinsSEDMP
    implicit none
    integer, intent(in) :: i_mp
    real,    intent(in) :: chi0, q, Emin, Emax
    integer :: i
    real    :: deltaE, logE0, logE1, m, comp

    !    print*, 'A0', A0,'m', m,'Emin', Emin,'Emax', Emax
    logE0 = LOG10(Emin)
    logE1 = LOG10(Emax)

    comp = q/(q-3.)
    if (comp > 2.) then

      m = q - 2.

      !  delta E in logerithmic bins
      deltaE = ( logE1 - logE0 ) / (real(NBinsSEDMP)-1.)

      do i = 1,NBinsSEDMP
        MP_SED(1,i,i_mp) = 10.**(logE0+real(i-1)*deltaE)
        MP_SED(2,i,i_mp) = chi0*MP_SED(1,i,i_mp)**(-m)
      end do

    end if

  end subroutine inject_PL_spectrum

  !================================================================
  !> @brief Obtain parameters for PL spectrum to be injected
  !> @details Get the parameters of new power-law spectrum to be imposed
  !> as a result of the subgrid recipe for the Diffuse Shock Acceleration.
  !> (Vaidya 2018 et al)
  !> @param real [in]  prim1(8) : primitives in the pre-shock region
  !> @param real [in]  prim2(8) : primitives in the post-shock region
  !> @param real [in]  rhoI     : density interpolated at the MP position
  !> @param real [in]  EI       : energy density interpolated at the MP position
  !> @param real [in]  BI       : B field interpolated at the MP position
  !> @param real [out] chi0     : Amplitude for the new PL
  !> @param real [out] qNR      : q index (Drury et al. 1983)
  !> @param real [out] Emin     : minimum energy (e0, Vaidya et al 2018)
  !> @param real [out] Emax     : maximum energy (e1, Vaidya et al 2018)
  subroutine get_PL_parameters(i_mp,prim1,prim2,rhoI,EI,BI,chi0,qNR,Emin,Emax)
    use parameters, only : rhosc, rsc, vsc2, NBinsSEDMP, Bsc, vsc2, Psc
    use constants,  only : eV, clight, echarge
    use globals, only : MP_SED, dx
    implicit none
    integer, intent(in) :: i_mp
    real, intent(in)    :: prim1(8), prim2(8), rhoI, EI, BI
    real, intent(out)   :: chi0, qNR, Emin, Emax
    real                :: normal(3), r, thB1, thB2, v1, v2, v_shock, Brat,    &
                           beta1sq, lambda_eff, Elimit
    real, parameter     :: lmp_eta = 4.25        !  for eq(32) in Vaidya et al.
    real, parameter     :: e1const = 37.702      ! m^2c^4 sqrt(9/(8pie^3))
    real, parameter     :: deltaN  = 0.01        ! Mimica et al. 2009
    real, parameter     :: deltaE  = 0.1         ! Mimica et al. 2009
    real :: e_old, n_old, c1, c2

    call get_NRth(prim1,prim2,normal,r,thB1,thB2)
    r = max(r,2.0)      !  this should be <= comp in if in injectcion routine
    r = min(r,4.0)      !  Test
    qNR = 3.*r/(r -1.)

    v1 = normal(1)*prim1(2)+normal(2)*prim1(3)+normal(3)*prim1(4)
    v2 = normal(1)*prim2(2)+normal(2)*prim2(3)+normal(3)*prim2(4)

    v_shock = v1-r*v2 /(1.-r)

    Brat = sqrt(prim1(6)**2+prim1(7)**2+prim1(8)**2)/                          &
           sqrt(prim2(6)**2+prim2(7)**2+prim2(8)**2)

    beta1sq = (v1-v_shock)**2*vsc2 / (clight**2)
    beta1sq = max(beta1sq, 1e-30)

    lambda_eff = lmp_eta*r / ( beta1sq*(r-1.) )                                &
               * (            cos(thB1)**2 + sin(thB1)**2/(1.+lmp_eta**2)      &
                   + r*Brat*( cos(thB2)**2 + sin(thB2)**2/(1.+lmp_eta**2) )  )

    !>  Upper limit to energy to remain inside one cell
    Elimit = 0.5*echarge*BI*Bsc*dx*rsc
    !  get E1
    Emax   = e1const / sqrt(BI*Bsc*lambda_eff)
    Emax = min(Emax,Elimit)

    call get_nE_old(NBinsSEDMP, MP_SED(:, :, i_mp), n_old, E_old)

    !  do not update enegy bounds (testing purposes)
    !Emax = MP_SED(1,100,i_mp)
    Emin = MP_SED(1,1,i_mp)

    !calculate Emin (see Esquivas notes)
    c1   = n_old + deltaN*rhoI
    c2   = E_old + deltaE*EI*Psc

    Emin = c2/c1 * (4.-qNR)/(3.-qNR)
    Emin = max(Emin,1e-8)
    !print*, '------', E_old, c2/c1 *(4.0-qNR)/(3.0-qNR)
    !    print*,"n_old, e_old", n_old,E_old
    chi0   =c1*(3.0-qNR)/(Emax**(3.0-qNR)-Emin**(3.0-qNR) ) /rhoI

  end subroutine get_PL_parameters

  !=======================================================================
  !> @brief Get nold and Eold
  !> @details Obtains n and E intrgrating the spectra using the trapezoidal
  !> rule (assuming logarithmic spaced bins)
  !> param integer [in]  nbins        : number of (log spaced) bins in SED
  !> param real    [in]  SED(2,nbins) : spectra, 1st index 0 is E, 1 is chi=N/n
  !> param real    [out]  nold        : total # of particles: int N(E) dE
  !> param real    [out]  Eold        : total Energy: int N(E) E dE
  subroutine get_nE_old(nbins,SED,nold,Eold)
    implicit none
    integer, intent(in)  :: nbins
    real,    intent(in)  :: SED(2,nbins)
    real,    intent(out) :: nold, Eold
    integer :: i
    real    :: slopeN, slopeE, Fn0, Fn1, FE0, FE1, x0, x1

    nold = 0.0
    Eold = 0.0

    do i = 1, nbins-1

      x0  = SED(1,i  )
      x1  = SED(1,i+1)
      Fn0 = SED(2,i  )
      Fn1 = SED(2,i+1)
      FE0 = SED(2,i  )*SED(1,i  )
      FE1 = SED(2,i+1)*SED(1,i+1)

      if (x1 /= x0) then
        slopeN = ( log(Fn1/Fn0) )/( log(x1/x0) )
        slopeE = ( log(FE1/FE0) )/( log(x1/x0) )
      else
        nold = -1.
        Eold = -1.
        return
      end if

      if (slopeN /= -1.) then
        nold = nold + Fn0/(slopeN+1.)*(x1*(x1/x0)**slopeN-x0)
      else
        nold = nold + Fn0*x0*log(x1/x0)
      end if

      if (slopeE /= -1.) then
        Eold = Eold + FE0/(slopeE+1.)*(x1*(x1/x0)**slopeE-x0)
      else
        Eold = Eold + FE0 * x0 * log(x1/x0)
      end if

    end do

  end subroutine get_nE_old

  !================================================================

end module lmp_module

!================================================================
