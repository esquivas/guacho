!=======================================================================
!> @file thermal_cond.f90
!> @brief Thermal conduction module
!> @author Alejandro Esquivel & Ernesto Zurbiggen
!> @date 07/Sep/2015

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
! along with this program.  If not, see httƒSATp://www.gnu.org/licenses/.
!=======================================================================

!> @brief Adds thermal conducion
!> @details Adds a thermal conduction term, affects both the primitive
!! and conserved variables

module thermal_cond
  use globals
  use parameters
  implicit none

  !> Parameter for the sturated regime in McKee
  real, parameter :: ph=0.4
  real, parameter :: nu=0.01            !< Super-stepping damping factor
  real, parameter :: snu=sqrt(nu)       !< Sqrt of damping factor
  integer, parameter ::  Max_iter = 100 !< Maximum number of iterations
  !> timestep reduction factor for the conduction
  real, parameter :: tstep_red_factor=0.25
  real :: dt_cond                       !< conduction timestep
  integer :: tc_log !< loical unit to write TC log

contains

!=======================================================================

!> @brief Intializes Temperature array
!> @details Intializes Temperature array
!> (to resolve dependencies it was moved to the globals module)

subroutine init_thermal_cond()
  implicit none

  !  create log dir if not present
  if (rank == master) then

  !  call execute_command_line('if [ ! -e '//trim(workdir)&
  !    //'logs ]; then mkdir '//trim(workdir)//'logs ; fi')

  call system('if [ ! -e '//trim(workdir)&
      //'logs ]; then mkdir '//trim(workdir)//'logs ; fi')

  open(newunit=tc_log,file=trim(workdir)//'logs/thermal_conduction.log')
  write(tc_log,'(a)') '***** Thermal conduction logfile ********'
  write(tc_log,'(a)') 'iteration    |    dt_hydro    |   dt_cond      | Nsteps  '

  end if

end subroutine init_thermal_cond

!=======================================================================

!> @brief computes conduction timescale
!> @details computes conduction timescale (in seconds)
!> @param real [out] dt :: conduction timescale

subroutine get_dt_cond(dt)

  use hydro_core, only : csound
  implicit none
  real, intent(out) :: dt
  real              :: dtp, ddx
  integer :: i, j, k, err
  !
  dtp=huge(1.)
  ddx=min(dx ,dy)
  ddx=min(ddx,dz)
  !
  do  k=1,nz
     do j=1,ny
        do i=1,nx

           !  spitzer timescale
           dtp = min( dtp, primit(1,i,j,k)/Ksp(Temp(i,j,k)) )

        end do
     end do
  end do

  dtp=tstep_red_factor*0.5*(ddx*rsc)**2*cv*Rg*dtp*rhosc/mu

#ifdef MPIP
  call mpi_allreduce(dtp, dt, 1, mpi_real_kind, mpi_min, mpi_comm_world,err)
#else
  dt=dtp
#endif


end subroutine get_dt_cond

!=======================================================================

!> @brief Progress bar
!> @details Progress bar
!! takes a number between 1 and tot
!> @param integer [in] j   : current iteration
!> @param integer [in] tot : total number of iterartions

subroutine progress(j,tot)
  implicit none
  integer(kind=4)::j,k
  integer(kind=4), intent(in) :: tot
  character(len=57)::bar="???% |                                                  |"
  open (unit=6)
  write(unit=bar(1:3),fmt="(i3)") 100*j/tot
  bar(7:56)="."
  do k=1, 50*j/tot
     bar(6+k:6+k)="="
  enddo
  ! print the progress bar.
  write(unit=6,fmt="(a1,a1,a57)",advance="no") '+',char(13), bar
  return
end subroutine progress

!=======================================================================

!> @brief Spitzer conductivity
!> @details Computes the Spitzer conductivity
!> @param real [in] T : temperature [K]

real function KSp(T)
  implicit none
  real, intent(in) :: T
  real, parameter :: beta=6.e-7

  Ksp= beta*T**(2.5)

end function KSp

!=======================================================================

!> @brief Spitzer parallel conductivity
!> @details Computes the Spitzer conductivity parallel to B
!> @param real [in] T : temperature [K]

real function KSp_parl(xtemp)
  implicit none
  real,intent(in):: xtemp
  real,parameter:: K0_parl = 9.2181e-7

  Ksp_parl= K0_parl*xtemp**(2.5)

end function KSp_parl

!=======================================================================

!> @brief Spitzer perpendicular conductivity
!> @details Computes the Spitzer conductivity perpendicular to B
!> @param real [in] T : temperature [K]

real function KSp_perp(xtemp,xdens,B2)
  implicit none
  real,intent(in):: xtemp,xdens,B2
  real,parameter:: K0_perp = 0.30089e+33
  Ksp_perp= K0_perp*xdens/(B2*sqrt(xtemp))*xdens

end function KSp_perp


!=======================================================================

!> @brief Returns Heat Fluxes
!> @details Heat flux, if saturation enabled it takes minimum of the
!! Spitzer and the saturated value
!! @n  The result is stored in the 5th component of global the
!! F,G,H fluxes (in cgs, conversion is done in dt product)

subroutine heatfluxes()
  use hydro_core, only : csound
  implicit none
  integer :: i, j, k
  real, parameter :: clight=3.E10, phi=0.3
  real:: cs, coef, dTx, dTy, dTz, meanT, meanP, meanDens, yhp

  do k=0,nz
     do j=0,ny
        do i=0,nx

           yhp=1.!-primit(neqdyn+1,i,j,k)/primit(1,i,j,k)

           !  get the flux in the X direction
           if (Temp(i,j,k) == Temp(i+1,j,k) ) then
              F(5,i,j,k)=0.
           else
              meanP   = 0.5*(primit(5,i,j,k)+primit(5,i+1,j,k))
              meanDens= 0.5*(primit(1,i,j,k)+primit(1,i+1,j,k))
              meanT   = 0.5*(    Temp(i,j,k)+    Temp(i+1,j,k))
              dTx=(Temp(i+1,j,k)-Temp(i,j,k))/(dx*rsc)

          if (tc_saturation) then
              call csound(meanP,meanDens,cs)
              cs=min(cs*sqrt(vsc2),clight)
              coef=min( Ksp(meanT) , 5.*ph*cs*meanP*Psc/abs(dTx) )
          else
              coef = Ksp(meanT)
          end if

              F(5,i,j,k)=-coef*dTx*yhp

           end if
           !  get the flux in the Y direction
           if (Temp(i,j,k) == Temp(i,j+1,k) ) then
              G(5,i,j,k)=0.
           else
              meanP   = 0.5*(primit(5,i,j,k)+primit(5,i,j+1,k))
              meanDens= 0.5*(primit(1,i,j,k)+primit(1,i,j+1,k))
              meanT   = 0.5*(    Temp(i,j,k)+    Temp(i,j+1,k))
              dTy=(Temp(i,j+1,k)-Temp(i,j,k))/(dy*rsc)

              if (tc_saturation) then
                call csound(meanP,meanDens,cs)
                cs=min(cs*sqrt(vsc2),clight)
                coef=min( Ksp(meanT) , 5.*ph*cs*meanP*Psc/abs(dTy) )
              else
                coef = Ksp(meanT)
              endif

              G(5,i,j,k)=-coef*dTy*yhp

           end if
           !  get the flux in the Z direction
           if (Temp(i,j,k) == Temp(i,j,k+1) ) then
              H(5,i,j,k)=0.
           else
              meanP   = 0.5*(primit(5,i,j,k)+primit(5,i,j,k+1))
              meanDens= 0.5*(primit(1,i,j,k)+primit(1,i,j,k+1))
              meanT   = 0.5*(  Temp(i,j,k)  +    Temp(i,j,k+1))
              dTz=(Temp(i,j,k+1)-Temp(i,j,k))/(dz*rsc)

              if (tc_saturation) then
                call csound(meanP,meanDens,cs)
                cs=min(cs*sqrt(vsc2),clight)
                coef=min( Ksp(meanT) , 5.*ph*cs*meanP*Psc/abs(dTz) )
              else
                coef = Ksp(meanT)
              endif

              H(5,i,j,k)=-coef*dTz*yhp

           end if

        end do
     end do
  end do
  !
end subroutine heatfluxes

!=======================================================================

!> @brief Returns Heat Fluxes with anisotropic thermal conduction
!> @details Heat flux, if sturation enabled takes minimum of the
!! Spitzer and the saturated value
!! @n  The result is stored in the 5th component of global the
!! F,G,H fluxes (in cgs, conversion is done in dt product)

subroutine MHD_heatfluxes()

  use hydro_core, only : csound
  implicit none
  integer :: i,j,k
  real,parameter :: clight=3.E10,phi=0.3,alpha=5.0*phi
  real :: cs,meanTemp,meanPres,meanDens
  !Gradient Temperature
  real :: gradTx,gradTy,gradTz
  !real :: gradT
  real :: gradT_parl_x,gradT_parl_y,gradT_parl_z,gradT_parl
  real :: gradT_perp_x,gradT_perp_y,gradT_perp_z,gradT_perp
  !Magnetic fiel unitary vector components
  real :: bx,by,bz,modB,B2
  !Internal product of vectors: B and gradT
  real :: bgradT
  !real :: ,bgradTx,bgradTy,bgradTz
  !Coeficientes de conducción
  real :: K_parl_x,K_parl_y,K_parl_z,K_perp_x,K_perp_y,K_perp_z
  real :: coefSatx,coefSaty,coefSatz

  F(5,:,:,:)=0.0  ; G(5,:,:,:)=0.0 ; H(5,:,:,:)=0.0

  do k=0,nz
     do j=0,ny
        do i=0,nx

          bx = primit(6,i,j,k)
          by = primit(7,i,j,k)
          bz = primit(8,i,j,k)
          B2 = bx*bx+by*by+bz*bz
          modB = sqrt(B2)
          bx = bx/modB
          by = by/modB
          bz = bz/modB

          !  not saturated conduction
          if(.not.tc_saturation) then
            !  get the flux in the X direction
            if ( abs(Temp(i,j,k)-Temp(i+1,j,k)) < 1.0e-14 ) then

              gradTx = 0.0
              K_parl_x = 0.0
              K_perp_x = 0.0

            else

              meanDens = 0.5*(primit(1,i,j,k)+primit(1,i+1,j,k))
              meanTemp = 0.5*(    Temp(i,j,k)+    Temp(i+1,j,k))

              gradTx = (Temp(i+1,j,k)-Temp(i,j,k))/(dx*rsc)
              K_parl_x = Ksp_parl(meanTemp)
              K_perp_x = Ksp_perp(meanTemp,meanDens*rhosc,B2*bsc**2)

            end if

            !  get the flux in the Y direction
            if ( abs(Temp(i,j,k)-Temp(i,j+1,k)) < 1.0e-14 ) then

              gradTy = 0.0
              K_parl_y = 0.0
              K_perp_y = 0.0

            else

              meanDens = 0.5*(primit(1,i,j,k)+primit(1,i,j+1,k))
              meanTemp = 0.5*(    Temp(i,j,k)+    Temp(i,j+1,k))

              gradTy = (Temp(i,j+1,k)-Temp(i,j,k))/(dy*rsc)
              K_parl_y = Ksp_parl(meanTemp)
              K_perp_y = Ksp_perp(meanTemp,meanDens*rhosc,B2*bsc**2)

            end if

            !  get the flux in the Z direction
            if ( abs(Temp(i,j,k)-Temp(i,j,k+1)) < 1.0e-14 ) then

              gradTz = 0.0
              K_parl_z = 0.0
              K_perp_z = 0.0

            else

              meanDens = 0.5*(primit(1,i,j,k)+primit(1,i,j,k+1))
              meanTemp = 0.5*(  Temp(i,j,k)  +    Temp(i,j,k+1))

              gradTz = (Temp(i,j,k+1)-Temp(i,j,k))/(dz*rsc)
              K_parl_z = Ksp_parl(meanTemp)
              K_perp_z = Ksp_perp(meanTemp,meanDens*rhosc,B2*bsc**2)

            end if

            !internal product of b.gradT
            bgradT = bx*gradTx+by*gradTy+bz*gradTz

            gradT_parl_x = bgradT*bx
            gradT_parl_y = bgradT*by
            gradT_parl_z = bgradT*bz

            gradT_perp_x = gradTx-gradT_parl_x
            gradT_perp_y = gradTy-gradT_parl_y
            gradT_perp_z = gradTz-gradT_parl_z

            F(5,i,j,k) = -K_parl_x*gradT_parl_x - K_perp_x*gradT_perp_x

            G(5,i,j,k) = -K_parl_y*gradT_parl_y - K_perp_y*gradT_perp_y

            H(5,i,j,k) = -K_parl_z*gradT_parl_z - K_perp_z*gradT_perp_z

          else
            !   Saturated conduction
            !  get the flux in the X direction
            if ( abs(Temp(i,j,k)-Temp(i+1,j,k)) < 1.0e-14 ) then

              gradTx = 0.0
              K_parl_x = 0.0
              K_perp_x = 0.0
              coefSatx = 0.0

            else

              meanPres = 0.5*(primit(5,i,j,k)+primit(5,i+1,j,k))
              meanDens = 0.5*(primit(1,i,j,k)+primit(1,i+1,j,k))
              meanTemp = 0.5*(    Temp(i,j,k)+    Temp(i+1,j,k))
              call csound(meanPres,meanDens,cs)
              cs = min(cs*vsc,clight)
              coefSatx = alpha*meanDens*cs**3

              gradTx = (Temp(i+1,j,k)-Temp(i,j,k))/(dx*rsc)
              K_parl_x = Ksp_parl(meanTemp)
              K_perp_x = Ksp_perp(meanTemp,meanDens*rhosc,B2*bsc**2)

            end if

             !  get the flux in the Y direction
            if ( abs(Temp(i,j,k)-Temp(i,j+1,k)) < 1.0e-14 ) then

              gradTy = 0.0
              K_parl_y = 0.0
              K_perp_y = 0.0
              coefSaty = 0.0

            else

              meanPres = 0.5*(primit(5,i,j,k)+primit(5,i,j+1,k))
              meanDens = 0.5*(primit(1,i,j,k)+primit(1,i,j+1,k))
              meanTemp = 0.5*(    Temp(i,j,k)+    Temp(i,j+1,k))
              call csound(meanPres,meanDens,cs)
              cs = min(cs*vsc,clight)
              coefSaty = alpha*meanDens*cs**3

              gradTy = (Temp(i,j+1,k)-Temp(i,j,k))/(dy*rsc)
              K_parl_y = Ksp_parl(meanTemp)
              K_perp_y = Ksp_perp(meanTemp,meanDens*rhosc,B2*bsc**2)

            end if

            !  get the flux in the Z direction
            if ( abs(Temp(i,j,k)-Temp(i,j,k+1)) < 1.0e-14 ) then

              gradTz = 0.0
              K_parl_z = 0.0
              K_perp_z = 0.0
              coefSatz = 0.0

            else

              meanPres = 0.5*(primit(5,i,j,k)+primit(5,i,j,k+1))
              meanDens = 0.5*(primit(1,i,j,k)+primit(1,i,j,k+1))
              meanTemp = 0.5*(  Temp(i,j,k)  +    Temp(i,j,k+1))
              call csound(meanPres,meanDens,cs)
              cs = min(cs*vsc,clight)
              coefSatz = alpha*meanDens*cs**3

              gradTz = (Temp(i,j,k+1)-Temp(i,j,k))/(dz*rsc)
              K_parl_z = Ksp_parl(meanTemp)
              K_perp_z = Ksp_perp(meanTemp,meanDens*rhosc,B2*bsc**2)

            end if

            !internal product of b.gradT
            bgradT = bx*gradTx+by*gradTy+bz*gradTz

            gradT_parl_x = bgradT*bx
            gradT_parl_y = bgradT*by
            gradT_parl_z = bgradT*bz
            ! |gradT_parl| == |(b.gradT)b| == |(b.gradT)|
            gradT_parl = bgradT

            gradT_perp_x = gradTx-gradT_parl_x
            gradT_perp_y = gradTy-gradT_parl_y
            gradT_perp_z = gradTz-gradT_parl_z
            ! |gradT_perp| == |gradT-gradT_parl|
            gradT_perp = sqrt( gradT_perp_x*gradT_perp_x + gradT_perp_y*gradT_perp_y + gradT_perp_z*gradT_perp_z )

            F(5,i,j,k) = - 1./( 1./(K_parl_x + 1.e-14) + gradT_parl/(coefSatx + 1.e-14) ) * gradT_parl_x &
                         - 1./( 1./(K_perp_x + 1.e-14) + gradT_perp/(coefSatx + 1.e-14) ) * gradT_perp_x

            G(5,i,j,k) = - 1./( 1./(K_parl_y + 1.e-14) + gradT_parl/(coefSaty + 1.e-14) ) * gradT_parl_y &
                         - 1./( 1./(K_perp_y + 1.e-14) + gradT_perp/(coefSaty + 1.e-14) ) * gradT_perp_y

            H(5,i,j,k) = - 1./( 1./(K_parl_z + 1.e-14) + gradT_parl/(coefSatz + 1.e-14) ) * gradT_parl_z &
                         - 1./( 1./(K_perp_z + 1.e-14) + gradT_perp/(coefSatz + 1.e-14) ) * gradT_perp_z

          end if

        end do
     end do
  end do

end subroutine MHD_heatfluxes


!=======================================================================

!> @brief Exchanges ghost cells for energy only
!> @details Exchanges one layer of boundaries, only the equation that
!!  corresponds to the energy

  subroutine thermal_bounds()

    implicit none
    integer, parameter :: nxp1=nx+1
    integer, parameter :: nyp1=ny+1
    integer, parameter :: nzp1=nz+1
#ifdef MPIP
    integer:: status(MPI_STATUS_SIZE), err
    real, dimension(1,1,0:nyp1,0:nzp1)::sendr,recvr,sendl,recvl
    real, dimension(1,0:nxp1,1,0:nzp1)::sendt,recvt,sendb,recvb
    real, dimension(1,0:nxp1,0:nyp1,1)::sendi,recvi,sendo,recvo
    integer, parameter :: bxsize=(ny+2)*(nz+2)
    integer, parameter :: bysize=(nx+2)*(nz+2)
    integer, parameter :: bzsize=(nx+2)*(ny+2)

    !   Exchange boundaries between processors

    !   boundaries to procs: right, left, top, bottom, in and out
    sendr(1,1,:,:)=u(5,nx    ,0:nyp1,0:nzp1)
    sendl(1,1,:,:)=u(5,1     ,0:nyp1,0:nzp1)
    sendt(1,:,1,:)=u(5,0:nxp1,ny    ,0:nzp1)
    sendb(1,:,1,:)=u(5,0:nxp1,1     ,0:nzp1)
    sendi(1,:,:,1)=u(5,0:nxp1,0:nyp1,nz    )
    sendo(1,:,:,1)=u(5,0:nxp1,0:nyp1,1     )
    !
    call mpi_sendrecv(sendr, bxsize, mpi_real_kind, right  ,0,          &
                      recvl, bxsize, mpi_real_kind, left   ,0,          &
                      comm3d, status , err)

    call mpi_sendrecv(sendt, bysize, mpi_real_kind, top    ,0,          &
                      recvb, bysize, mpi_real_kind, bottom ,0,          &
                      comm3d, status , err)

    call mpi_sendrecv(sendi, bzsize, mpi_real_kind, in     ,0,          &
                      recvo, bzsize, mpi_real_kind, out    ,0,          &
                      comm3d, status , err)

    call mpi_sendrecv(sendl, bxsize, mpi_real_kind, left  , 0,          &
                      recvr, bxsize, mpi_real_kind, right , 0,          &
                      comm3d, status , err)

    call mpi_sendrecv(sendb, bysize, mpi_real_kind, bottom, 0,          &
                      recvt, bysize, mpi_real_kind, top   , 0,          &
                      comm3d, status , err)

    call mpi_sendrecv(sendo, bzsize, mpi_real_kind, out   , 0,          &
                      recvi, bzsize, mpi_real_kind, in    , 0,          &
                      comm3d, status , err)

    if (left  .ne. -1) u(5,0     ,0:nyp1,0:nzp1)=recvl(1,1,:,:)
    if (right .ne. -1) u(5,nxp1  ,0:nyp1,0:nzp1)=recvr(1,1,:,:)
    if (bottom.ne. -1) u(5,0:nxp1,0     ,0:nzp1)=recvb(1,:,1,:)
    if (top   .ne. -1) u(5,0:nxp1,nyp1  ,0:nzp1)=recvt(1,:,1,:)
    if (out   .ne. -1) u(5,0:nxp1,0:nyp1,0     )=recvo(1,:,:,1)
    if (in    .ne. -1) u(5,0:nxp1,0:nyp1,nzp1  )=recvi(1,:,:,1)
    !
#else

    !   periodic BCs
    if (bc_left == BC_PERIODIC .and. bc_right == BC_PERIODIC) then
      !   Left BC
      if (coords(0).eq.0) then
         u(5,0,:,:)=u(5,nx,:,:)
      endif
      !   Right BC
      if (coords(0).eq.MPI_NBX-1) then
         u(5,nxp1,:,:)=u(5,1,:,:)
      endif
    end if

    if (bc_bottom == BC_PERIODIC .and. bc_top == BC_PERIODIC) then
      !   bottom BC
      if (coords(1).eq.0) then
         u(5,:,0,:)= u(5,:,ny,:)
      endif
      !   top BC
      if (coords(1).eq.MPI_NBY-1) then
         u(5,:,nyp1,:)= u(5,:,1,:)
      endif
    end if

    if (bc_out == BC_PERIODIC .and. bc_in == BC_PERIODIC) then
      !   out BC
      if (coords(2).eq.0) then
         u(5,:,:,0)= u(5,:,:,nz)
      endif
      !   in BC
      if (coords(2).eq.MPI_NBZ-1) then
         u(5,:,:,nzp1)= u(5,:,:,1)
      endif
    endif

#endif /* !MPIP */
    !   reflecting and outflow BCs

    !   left
    if (coords(0).eq.0) then
       u(5,0,   0:nyp1,0:nzp1)=u(5,1 ,0:nyp1,0:nzp1)
    endif
    !   right
    if (coords(0).eq.MPI_NBX-1) then
       u(5,nxp1,0:nyp1,0:nzp1)=u(5,nx,0:nyp1,0:nzp1)
    endif
    !   bottom
    if (coords(1).eq.0) then
       u(5,0:nxp1,0   ,0:nzp1)=u(5,0:nxp1,1 ,0:nzp1)
    endif
    !   top
    if (coords(1).eq.MPI_NBY-1) then
       u(5,0:nxp1,nyp1,0:nzp1)=u(5,0:nxp1,ny,0:nzp1)
    endif
    !   out
    if (coords(2).eq.0) then
       u(5,0:nxp1,0:nyp1,0   )=u(5,0:nxp1,0:nyp1,1 )
    endif
    !   in
    if (coords(2).eq.MPI_NBZ-1) then
       u(5,0:nxp1,0:nyp1,nzp1)=u(5,0:nxp1,0:nyp1,nz)
    endif

  end subroutine thermal_bounds

!=======================================================================

!> @brief  Length of superstep
!> @details Returns the length of the superstep with N inner substeps
!> @param integer [in] N : Nunber of inner substeps
!> @param real [in] snu : sqrt of daMPI_NBg factor

real function superstep(N,snu)

  implicit none
  integer :: N
  real,    intent(in) :: snu

  superstep=real(N)/(2.*snu) * ( (1+snu)**(2*N) - (1-snu)**(2*N) ) / &
       ( (1+snu)**(2*N) + (1-snu)**(2*N) )

  !1/( (nu-1.)*Cos(pi*(2*real(j)-1.)/(2.*real(N)) )+nu+1. )

end function superstep

!=======================================================================

!> @brief Size of substep j
!> @details Returns the size of substep j of N
!> @param  integer [in] j : index of current step
!> @param  integer [in] N : Total number of substeps
!> @param  real [in] nu : daMPI_NBg factor

real function substep(j,N,nu)

  implicit none
  integer, intent(in) :: j, N
  real,    intent(in) :: nu

  substep=1./( (nu-1.)*Cos(pi*real(2*j-1)/(2.*real(N)) )+nu+1. )

end function substep

!=======================================================================

!> @brief Returns the number of Supersteps
!> @details Returns the number of Supersteps
!> @param real fs    : ratio of dtcond/dthydro
!> @param integer Ns : Number of Supersteps
!> @param real fstep : Number of supersteps (float)

  subroutine ST_steps(fs,Ns,fstep)

    implicit none

    real ,   intent(in) :: fs
    integer, intent(out):: Ns
    real,    intent(out):: fstep
    integer, parameter :: jmax=199
    integer :: j
    !
    do j=1,jmax
       if (superstep(j,snu) > fs) exit
    end do

    Ns = j
    fstep = fs/superstep(Ns,snu)

  end subroutine ST_steps

!=======================================================================

!> @brief Upper level wrapper for thermal conduction
!> @details This routine adds the heat conduction, receives the hydro
!!  timestep in seconds, and assumes the primitives and Temp(i,j,k)
!!  arrays are updated

subroutine thermal_conduction()

  use hydro_core, only : calcprim
  implicit none
  real :: dt_hydro
  real :: dts, fstep
  integer :: n,i,j,k, nsteps
  logical :: SuperStep

  dt_hydro = dt_CFL*tsc

  !  get the conduction timescale
  call get_dt_cond(dt_cond)

  SuperStep = .True.

  if (dt_cond < dt_hydro) then

    if (SuperStep) then
       call ST_steps(dt_hydro/dt_cond,Nsteps,fstep)
    else
      Nsteps = min( ceiling(dt_hydro/dt_cond), Max_iter )
    end if

  else

    !   this oprevents use of superstep if Nsteps =1
    SuperStep = .False.
    fstep  = dt_hydro/dt_cond
    Nsteps = 1

  end if

  if (rank == master) then
    print*, 'Calculating thermal conduction'
    write(tc_log,'(i0,a,es15.7,a,es15.7,a,i4)') currentIteration,' |',dt_hydro,' |', dt_cond,' |', Nsteps
  end if

  steps : do n=1,Nsteps

    !  here i take care of the transformation from cgs to code units
    if (SuperStep) then
      dts=dt_cond*fstep*substep(n,Nsteps,nu)/Psc/rsc
    else
      dts=dt_hydro/real(Nsteps)/Psc/rsc
    end if

    !  show progress bar (comment to run in batch)
    if (rank == master ) call progress(n,nsteps)

    !  get the heat fluxes
    if (th_cond == TC_ANISOTROPIC) then
      call MHD_heatfluxes()
    end if
    if (th_cond == TC_ISOTROPIC) then
      call heatfluxes()
    end if

    !  update the conserved and primitives
    do k=1,nz
      do j=1,ny
        do i=1,nx
        u(5,i,j,k)=u(5,i,j,k)-dts*( ( f(5,i,j,k) - f(5,i-1,j,k) )/dx &
                                  + ( g(5,i,j,k) - g(5,i,j-1,k) )/dy &
                                  + ( h(5,i,j,k) - h(5,i,j,k-1) )/dz )
        end do
      end do
    end do

    !  boundary conditions
    !  (only one layer of u(5,:,:,:) is exchanged )
    call thermal_bounds()

    !  update primitives and Temperature
    call calcprim(u, primit)

  end do steps

end subroutine thermal_conduction

!=======================================================================

end module thermal_cond

!=======================================================================
