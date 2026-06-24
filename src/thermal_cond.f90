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
  use constants, only : clight, pi
  implicit none

  !> Parameter for the sturated regime in McKee
  real, parameter :: eps_perp = 1.0e-4  !< K_perp/K_par (field wandering turb)
  real, parameter :: phi=0.3
  real, parameter :: nu=0.01            !< ST damping factor ~[0.001-0.1]
  real, parameter :: beta_sp = 6.4E-8   !< Effective conductivity spitzer=6.4E-7

  !> Maximum number of iterations ST block, or subcycles if a ST not used
  integer, parameter :: Max_iter = 200

  !> timestep reduction factor for the conduction
  real, parameter :: tstep_red_factor=0.50   ! [0.1-0.5]
  real    :: dt_cond    !< conduction timestep
  real    :: dt_sat     !< conduction timestep including saturation
  logical :: kappa_eff  = .true. !< compute dt including saturation
  !>  use if too saturated or if dt_cond becomes too restrictive
  integer :: tc_log !< logical unit to write TC log

contains

!=======================================================================
!> @brief Intializes Temperature array
!> @details Intializes Temperature array
!> (to resolve dependencies it was moved to the globals module)
subroutine init_thermal_cond()
  implicit none

  !  create log dir if not present
  if (rank == master) then

  call system('if [ ! -e '//trim(outputpath)&
      //'logs ]; then mkdir '//trim(outputpath)//'logs ; fi')

  open(newunit=tc_log,file=trim(outputpath)//'logs/thermal_conduction-test.log')
  write(tc_log,'(a)') '************* Thermal conduction logfile ****************'
  write(tc_log,'(a)') '# iter  | dt_hydro/left |    dt_cond   | Nsteps | ST block'

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

          !  spitzer timescale (parabolic dt)
           dtp = min( dtp, primit(1,i,j,k)/Ksp(Temp(i,j,k)) )

        end do
     end do
  end do

  dtp=tstep_red_factor*(ddx*rsc)**2*cv*Rg*dtp*rhosc*mu/3.0  !  [seconds]

#ifdef MPIP
  call mpi_allreduce(dtp, dt, 1, mpi_real_kind, mpi_min, mpi_comm_world,err)
#else
  dt=dtp
#endif

! ==============================================================================
! Uses an effective kappa by including the saturated flux timestep
if ( tc_saturation .and. kappa_eff ) then
  if (TC_ISOTROPIC) then
    call heatfluxes()
  else if (TC_ANISOTROPIC) then
    call MHD_heatfluxes()
  end if

#ifdef MPIP
  call mpi_allreduce(dt_sat, dtp, 1, mpi_real_kind, mpi_min, mpi_comm_world,err)
#else
  dtp = dt_sat
#endif
  !if (rank == master) then
  !   print('(a,es15.7,a,es15.7,a,f7.3)'), &
  !   'dt_spitzer: ', dt, ' dt_saturated: ', dtp, ' speedup factor: ', dtp/dt
  !end if
  dt = max(dt, dtp)

  end if
  !=============================================================================

end subroutine get_dt_cond

!=======================================================================
!> @brief Progress bar
!> @details Progress bar
!! takes a number between 1 and tot
!> @param integer [in] j   : current iteration
!> @param integer [in] tot : total number of iterartions
subroutine progress(j, tot, done)
  implicit none

  integer, intent(in) :: j, tot
  logical, intent(in), optional :: done

  integer :: k, nfill, pct
  integer, parameter :: nbar = 37
  character(len=53) :: bar
  character(len=53) :: blank
  logical :: finished

  finished = .false.
  if (present(done)) finished = done

  blank = repeat(" ", len(blank))

  if (finished) then
     ! Clear the progress-bar line and return to beginning of line
     write(6, "(a1,a53,a1)", advance="no") char(13), blank, char(13)
     flush(6)
     return
  endif

  if (tot <= 0) return

  pct = int(100.0d0 * dble(j) / dble(tot))
  pct = max(0, min(100, pct))

  nfill = int(dble(nbar) * dble(j) / dble(tot))
  nfill = max(0, min(nbar, nfill))

  bar = "th cond: ???% |                                     |"

  write(bar(10:12), "(i3)") pct

  bar(16:52) = "."

  do k = 1, nfill
     bar(15+k:15+k) = "="
  end do

  ! Print progress bar without newline
  write(6, "(a1,a53)", advance="no") char(13), bar
  flush(6)

end subroutine progress

!=======================================================================
!> @brief Spitzer conductivity
!> @details Computes the Spitzer conductivity
!> @param real [in] T : temperature [K]
real function KSp(T)
  implicit none
  real, intent(in) :: T

  Ksp = beta_sp*T**(2.5)

end function KSp

!=======================================================================
!> @brief Spitzer parallel conductivity
!> @details Computes the Spitzer conductivity parallel to B
!> @param real [in] T : temperature [K]
real function KSp_par(Temp)
  implicit none
  real,intent(in):: Temp

  Ksp_par = beta_sp*Temp**(2.5)

end function KSp_par

!=======================================================================
!> @brief Spitzer perpendicular conductivity
!> @details Computes the Spitzer conductivity perpendicular to B
!> @param real [in] T : temperature [K]
real function KSp_perp(Temp)
  implicit none
  real,intent(in):: Temp

  Ksp_perp = beta_sp*eps_perp*Temp**(2.5)

end function KSp_perp

!=======================================================================

!> @brief Returns Heat Fluxes
!> @details Heat flux, if saturation enabled the flux is limited with it
!! @n  The result is stored in the 5th component of global the
!! F,G,H fluxes (in cgs, conversion is done in dt product)
!! also compute a dt_sat in case the saturation is enabled
subroutine heatfluxes()
  use hydro_core, only : csound
  implicit none
  integer :: i, j, k
  real, parameter :: phi=0.3
  real :: cs, meanP, meanDens, kap_x, kap_y, kap_z
  real :: qx_sat, qy_sat, qz_sat, qx_cl, qy_cl, qz_cl, dTdx, dTdy, dTdz, fac

  F(5,:,:,:)=0.0  ; G(5,:,:,:)=0.0 ; H(5,:,:,:)=0.0
  dt_sat = huge(1.)

  do k=0,nz
     do j=0,ny
        do i=0,nx

          !--------------- X direction -----------------------------------------
          ! Temperature gradient in x
          dTdx = (Temp(i+1,j,k) - Temp(i,j,k)) / (dx*rsc)
          ! Kappa at face center
          kap_x = 2.0*Ksp(Temp(i,j,k))*Ksp(Temp(i+1,j,k)) / &
                  ( Ksp(Temp(i,j,k)) + Ksp(Temp(i+1,j,k)))
          ! Classic (spitzer) heatflux
          qx_cl = -kap_x * dTdx
          if (tc_saturation) then
            ! Cowley & McKee saturation
            meanDens= 0.5*(primit(1,i,j,k)+primit(1,i+1,j,k))
            meanP   = 0.5*(primit(5,i,j,k)+primit(5,i+1,j,k))
            call csound(meanP,meanDens,cs)
            cs=min(cs*vsc,clight)         ! scale and limit to < clight
            qx_sat = 5.0 * phi * meanDens*rhosc * cs**3
            ! harmonic average
            fac = 1.0 / ( 1.0 + abs(qx_cl)/qx_sat )
            F(5,i,j,k) = qx_cl *fac
            dt_sat = min(dt_sat, meanDens*rhosc/fac/kap_x )
          else
            F(5,i,j,k) = qx_cl
          end if

          !--------------- Y direction -----------------------------------------
          ! Temperature gradient in y
          dTdy = (Temp(i,j+1,k) - Temp(i,j,k)) / (dy*rsc)
          ! Kappa at face center
          kap_y = 2.0*Ksp(Temp(i,j,k))*Ksp(Temp(i,j+1,k)) / &
                  ( Ksp(Temp(i,j,k)) + Ksp(Temp(i,j+1,k)) )
          ! Classic (spitzer) heatflux
          qy_cl = -kap_y * dTdy
          if (tc_saturation) then
            ! Cowley & McKee saturation
            meanDens= 0.5*(primit(1,i,j,k)+primit(1,i,j+1,k))
            meanP   = 0.5*(primit(5,i,j,k)+primit(5,i,j+1,k))
            call csound(meanP,meanDens,cs)
            cs=min(cs*vsc,clight)         ! scale and limit to < clight
            qy_sat = 5.0 * phi * meanDens*rhosc * cs**3
            ! harmonic average'
            fac = 1./ ( 1.0 + abs(qy_cl)/qy_sat )
            G(5,i,j,k) = qy_cl * fac
            dt_sat = min(dt_sat, meanDens*rhosc/fac/kap_y )
          else
            G(5,i,j,k) = qy_cl
          end if

          !--------------- Z direction -----------------------------------------
          ! Temperature gradient in z
          dTdz = (Temp(i,j,k+1) - Temp(i,j,k)) / (dz*rsc)
          ! Kappa at face center
          kap_z = 2.0*Ksp(Temp(i,j,k))*Ksp(Temp(i,j,k+1)) / &
                  ( Ksp(Temp(i,j,k)) + Ksp(Temp(i,j,k+1)) )
          ! Classic (spitzer) heatflux
          qz_cl = -kap_z * dTdz
          if (tc_saturation) then
            ! Cowley & McKee saturation
            meanDens= 0.5*(primit(1,i,j,k)+primit(1,i,j,k+1))
            meanP   = 0.5*(primit(5,i,j,k)+primit(5,i,j,k+1))
            call csound(meanP,meanDens,cs)
            cs=min(cs*vsc,clight)         ! scale and limit to < clight
            qz_sat = 5.0 * phi * meanDens*rhosc * cs**3
            ! harmonic average
            fac = 1.0 / ( 1.0 + abs(qz_cl)/qz_sat )
            H(5,i,j,k) = qz_cl * fac
            dt_sat = min(dt_sat, meanDens*rhosc/fac/kap_z )
          else
            H(5,i,j,k) = qz_cl
          end if

        end do
     end do
  end do

  !   here I multiply rho/kappa_eff by the prefactors to get dt_sat
  dt_sat = dt_sat * tstep_red_factor*(dx*rsc)**2*cv*Rg*mu/3.0

end subroutine heatfluxes

!=======================================================================
!> @brief Returns Heat Fluxes with anisotropic thermal conduction
!> @details Heat flux, if sturation enabled takes minimum of the
!! Spitzer and the saturated value
!! @n  The result is stored in the 5th component of global the
!! F,G,H fluxes (in cgs, conversion is done in dt product)
!! also compute a dt_sat in case the saturation is enabled
subroutine MHD_heatfluxes()

  use hydro_core, only : csound
  implicit none
  integer :: i,j,k
  real, parameter :: phi=0.3
  real            :: Bx_f, By_f, Bz_f, Bmag, bxh, byh, bzh, dTdx, dTdy, dTdz
  real            :: dTdy_i, dTdz_i, dTdy_ip, dTdz_ip,                         &
                     dTdx_j, dTdz_j, dTdx_jp, dTdz_jp,                         &
                     dTdx_k, dTdy_k, dTdx_kp, dTdy_kp, Kpar_f, Kperp_f,        &
                     bdotgradT, qx_cl, qy_cl, qz_cl, qmag_cl, rho_f, p_f, cs,  &
                     fac, qsat

  F(5,:,:,:)=0.0  ; G(5,:,:,:)=0.0 ; H(5,:,:,:)=0.0

  dt_sat = huge(1.)

  do k=0,nz
    do j=0,ny
      do i=0,nx

        !--------------- X direction -------------------------------------------
        !  Interpolate magnetic field to i+1/2 interface
        Bx_f = 0.5*(primit(6,i,j,k) + primit(6,i+1,j,k))
        By_f = 0.5*(primit(7,i,j,k) + primit(7,i+1,j,k))
        Bz_f = 0.5*(primit(8,i,j,k) + primit(8,i+1,j,k))
        Bmag = sqrt( Bx_f**2 + By_f**2 + Bz_f**2 )+ tiny(1.0)
        !  (scaling units not needed for the following)
        bxh = Bx_f / Bmag
        byh = By_f / Bmag
        bzh = Bz_f / Bmag
        !  Temperature gradient at interface
        dTdx    = (Temp(i+1,j,k) - Temp(i,j,k) ) / ( dx*rsc)
        dTdy_i  = ( Temp( i ,j+1,k) - Temp( i ,j-1,k) ) / (2.0*dy*rsc)
        dTdy_ip = ( Temp(i+1,j+1,k) - Temp(i+1,j-1,k) ) / (2.0*dy*rsc)
        dTdy    = 0.5*(dTdy_i + dTdy_ip)
        dTdz_i  = ( Temp( i ,j,k+1) - Temp( i ,j,k-1) ) / (2.0*dz*rsc)
        dTdz_ip = ( Temp(i+1,j,k+1) - Temp(i+1,j,k-1) ) / (2.0*dz*rsc)
        dTdz    = 0.5*(dTdz_i + dTdz_ip)
        !  Interpolate conductivities to inteface
        Kpar_f  = 2.0*KSp_par(Temp(i,j,k)) * KSp_par(Temp(i+1,j,k)) /          &
                  (   Ksp_par(Temp(i,j,k)) + Ksp_par(Temp(i+1,j,k)) )
        Kperp_f = 2.0*KSp_perp(Temp(i,j,k)) * KSp_perp(Temp(i+1,j,k)) /        &
                  (   KSp_perp(Temp(i,j,k)) + Ksp_perp(Temp(i+1,j,k)) )
        !  b dot grad(T)
        bdotgradT = bxh*dTdx + byh*dTdy + bzh*dTdz
        !  compute classical flux at interface
        qx_cl = -Kperp_f * dTdx -( Kpar_f - Kperp_f ) * bxh * bdotgradT
        qy_cl = -Kperp_f * dTdy -( Kpar_f - Kperp_f ) * byh * bdotgradT
        qz_cl = -Kperp_f * dTdz -( Kpar_f - Kperp_f ) * bzh * bdotgradT
        !  Cowie & McKee saturation
        if (tc_saturation) then
          rho_f = 0.5*(primit(1,i,j,k)+primit(1,i+1,j,k))
          p_f   = 0.5*(primit(5,i,j,k)+primit(5,i+1,j,k))
          call csound(p_f,rho_f,cs)
          cs=min(cs*vsc,clight)         ! scale and limit to < clight
          qsat = 5.0 * phi * rho_f*rhosc * cs**3
          qmag_cl = sqrt(qx_cl**2 + qy_cl**2 + qz_cl**2)
          fac = 1.0 / (1.0 + qmag_cl/max(qsat, tiny(1.0)))
          dt_sat = min(dt_sat, rho_f*rhosc/fac/Kpar_f )
        else
          fac = 1.0
        end if
        F(5,i,j,k) = qx_cl * fac


        !--------------- Y direction -------------------------------------------
        !  Interpolate magnetic field to j+1/2 interface
        Bx_f = 0.5*(primit(6,i,j,k) + primit(6,i,j+1,k))
        By_f = 0.5*(primit(7,i,j,k) + primit(7,i,j+1,k))
        Bz_f = 0.5*(primit(8,i,j,k) + primit(8,i,j+1,k))
        Bmag = sqrt( Bx_f**2 + By_f**2 + Bz_f**2 ) + tiny(1.0)
        !  (scaling units not needed for the following)
        bxh = Bx_f / Bmag
        byh = By_f / Bmag
        bzh = Bz_f / Bmag
        !  Temperature gradient at interface
        dTdx_j  = ( Temp(i+1, j ,k) - Temp(i-1, j ,k) ) / (2.0*dx*rsc)
        dTdx_jp = ( Temp(i+1,j+1,k) - Temp(i-1,j+1,k) ) / (2.0*dx*rsc)
        dTdx    = 0.5*(dTdx_j + dTdx_jp)
        dTdy    = ( Temp(i,j+1,k) - Temp(i,j,k) ) / ( dy*rsc)
        dTdz_j  = ( Temp(i, j ,k+1) - Temp(i, j ,k-1) ) / (2.0*dz*rsc)
        dTdz_jp = ( Temp(i,j+1,k+1) - Temp(i,j+1,k-1) ) / (2.0*dz*rsc)
        dTdz    = 0.5*(dTdz_j + dTdz_jp)
        !  Interpolate conductivities to inteface
        Kpar_f  = 2.0*KSp_par(Temp(i,j,k)) * KSp_par(Temp(i,j+1,k)) /          &
                  (   KSp_par(Temp(i,j,k)) + KSp_par(Temp(i,j+1,k)))
        Kperp_f = 2.0*KSp_perp(Temp(i,j,k)) * KSp_perp(Temp(i,j+1,k)) /        &
                  (   KSp_perp(Temp(i,j,k)) + Ksp_perp(Temp(i,j+1,k)))
        !  b dot grad(T)
        bdotgradT = bxh*dTdx + byh*dTdy + bzh*dTdz
        !  compute classical flux at interface
        qx_cl = -Kperp_f * dTdx -( Kpar_f - Kperp_f ) * bxh * bdotgradT
        qy_cl = -Kperp_f * dTdy -( Kpar_f - Kperp_f ) * byh * bdotgradT
        qz_cl = -Kperp_f * dTdz -( Kpar_f - Kperp_f ) * bzh * bdotgradT
        !  Cowie & McKee saturation
        if (tc_saturation) then
          rho_f = 0.5*(primit(1,i,j,k)+primit(1,i,j+1,k))
          p_f   = 0.5*(primit(5,i,j,k)+primit(5,i,j+1,k))
          call csound(p_f,rho_f,cs)
          cs=min(cs*vsc,clight)         ! scale and limit to < clight
          qsat = 5.0 * phi * rho_f*rhosc * cs**3
          qmag_cl = sqrt(qx_cl**2 + qy_cl**2 + qz_cl**2)
          fac = 1.0 / (1.0 + qmag_cl/max(qsat, tiny(1.0)))
          dt_sat = min(dt_sat, rho_f*rhosc/fac/Kpar_f )
        else
          fac = 1.0
        end if
        G(5,i,j,k) = qy_cl * fac

        !--------------- Z direction -------------------------------------------
        !  Interpolate magnetic field to k+1/2 interface
        Bx_f = 0.5*(primit(6,i,j,k) + primit(6,i,j,k+1))
        By_f = 0.5*(primit(7,i,j,k) + primit(7,i,j,k+1))
        Bz_f = 0.5*(primit(8,i,j,k) + primit(8,i,j,k+1))
        Bmag = sqrt( Bx_f**2 + By_f**2 + Bz_f**2 )+ tiny(1.0)
        !  (scaling units not needed for the following)
        bxh = Bx_f / Bmag
        byh = By_f / Bmag
        bzh = Bz_f / Bmag
        !  Temperature gradient at interface
        dTdx_k  = ( Temp(i+1,j, k ) - Temp(i-1,j, k ) ) / (2.0*dx*rsc)
        dTdx_kp = ( Temp(i+1,j,k+1) - Temp(i-1,j,k+1) ) / (2.0*dx*rsc)
        dTdx    = 0.5*(dTdx_k + dTdx_kp)
        dTdy_k  = ( Temp(i,j+1, k ) - Temp(i,j-1, k) ) / (2.0*dy*rsc)
        dTdy_kp = ( Temp(i,j+1,k+1) - Temp(i,j-1,k+1) ) / (2.0*dy*rsc)
        dTdy    = 0.5*(dTdy_k + dTdy_kp)
        dTdz    = ( Temp(i,j,k+1) - Temp(i,j,k) ) / ( dz*rsc)
        !  Interpolate conductivities to inteface
        Kpar_f  = 2.0*KSp_par(Temp(i,j,k)) * KSp_par(Temp(i,j,k+1)) /          &
                  (   KSp_par(Temp(i,j,k)) + KSp_par(Temp(i,j,k+1)))
        Kperp_f = 2.0*KSp_perp(Temp(i,j,k)) * KSp_perp(Temp(i,j,k+1)) /        &
                  (   KSp_perp(Temp(i,j,k)) + Ksp_perp(Temp(i,j,k+1)))
        !  b dot grad(T)
        bdotgradT = bxh*dTdx + byh*dTdy + bzh*dTdz
        !  compute classical flux at interface
        qx_cl = -Kperp_f * dTdx -( Kpar_f - Kperp_f ) * bxh * bdotgradT
        qy_cl = -Kperp_f * dTdy -( Kpar_f - Kperp_f ) * byh * bdotgradT
        qz_cl = -Kperp_f * dTdz -( Kpar_f - Kperp_f ) * bzh * bdotgradT
        !  Cowie & McKee saturation
        if (tc_saturation) then
          rho_f = 0.5*(primit(1,i,j,k)+primit(1,i,j,k+1))
          p_f   = 0.5*(primit(5,i,j,k)+primit(5,i,j,k+1))
          call csound(p_f,rho_f,cs)
          cs=min(cs*vsc,clight)         ! scale and limit to < clight
          qsat = 5.0 * phi * rho_f*rhosc * cs**3
          qmag_cl = sqrt(qx_cl**2 + qy_cl**2 + qz_cl**2)
          fac = 1.0 / (1.0 + qmag_cl/max(qsat, tiny(1.0)))
          dt_sat = min(dt_sat, rho_f*rhosc/fac/Kpar_f )
        else
          fac = 1.0
        end if
        H(5,i,j,k) = qz_cl * fac
        end do
     end do
  end do

  !   here I multiply rho/kappa_eff for the constants in front of dt_sat
  dt_sat = dt_sat * tstep_red_factor*(dx*rsc)**2*cv*Rg*mu/3.0

end subroutine MHD_heatfluxes

!=======================================================================
!> @brief Exchanges ghost cells for energy only
!> @details Exchanges one layer of boundaries, only the equation that
!>  corresponds to the energy
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
! ********   Thius function only exact in th elimit of N large ****************
!real function superstep(N,snu)
!
!  implicit none
!  integer :: N
!  real,    intent(in) :: snu
!
!  superstep=real(N)/(2.*snu) * ( (1+snu)**(2*N) - (1-snu)**(2*N) ) / &
!       ( (1+snu)**(2*N) + (1-snu)**(2*N) )
!
!  !1/( (nu-1.)*Cos(pi*(2*real(j)-1.)/(2.*real(N)) )+nu+1. )
!
!end function superstep

!=======================================================================
!> @brief Size of substep j
!> @details Returns the size of substep j of N
!> @param  integer [in] j : index of current step
!> @param  integer [in] N : Total number of substeps
!> @param  real [in] nu : damping factor
real function substep(j,N,nu)

  implicit none
  integer, intent(in) :: j, N
  real,    intent(in) :: nu

  substep=1.0/( (nu-1.0)*Cos(pi*real(2*j-1)/(2.0*real(N)) ) + nu + 1.0 )

end function substep

!=======================================================================
!> @brief Returns the number of Supersteps
!> @details Returns the number of Supersteps
!> @param real dt_hydro   : dt hydro [seconds]
!> @param real dt_cond    : dt of classical conduction
!> @param integer Ns      : Number of Supersteps required
!> @param real scale      : Scaling factor to match up dt_cfl with Ns supersteps
subroutine ST_steps(dt_hydro,dt_cond, Ns,scale)

  implicit none
  real ,   intent(in)  :: dt_hydro, dt_cond
  integer, intent(out) :: Ns
  real,    intent(out) :: scale
  real                 :: dt_tot
  integer :: j

  Ns = 1
  do
    dt_tot = 0.0
    do j = 1, Ns
      dt_tot = dt_tot + dt_cond * substep(j, Ns, nu)
    end do

    if (dt_tot >= dt_hydro) exit

    Ns = Ns + 1

    if (Ns > Max_iter) then
      !print*, "Error: number of supersteps exceeded, increase Max_iter or" //&
      !        "reduce the damping factor nu ", Ns, Max_iter, dt_hydro, dt_cond
      !stop
      Ns = Max_iter
      scale = 1.0
      return !exit
    end if

  end do

  scale = dt_hydro / dt_tot

end subroutine ST_steps

!=======================================================================
!> @brief updates pressure and temperature
!> @details This routine makes a minimal update of only pressure and
!> temperature, takes the indices i, j, k of a cell, updates the global arrays
!> 'primit' and 'Temp'
!> @param integer [in] i : index in the X direction
!> @param integer [in] j : index in the Y direction
!> @param integer [in] k : index in the Z direction
!!**********  support for different EOS is under construction **************
subroutine update_PT(i, j, k)

  implicit none
  integer, intent(in) :: i, j, k
  real, parameter :: Temp_floor = 10.0

  if (mhd) then
    primit(5,i,j,k) = ( u(5,i,j,k)                                             &
            - 0.5 * ( u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)/u(1,i,j,k) )  &
            - 0.5 * ( u(6,i,j,k)**2+u(7,i,j,k)**2+u(8,i,j,k) ) /Cv
  else
    primit(5,i,j,k) = ( u(5,i,j,k)                                             &
            - 0.5 * ( u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)/u(1,i,j,k) )
  end if

  if (eq_of_state == EOS_SINGLE_SPECIE) then
    Temp(i,j,k) = max(Temp_floor,(primit(5,i,j,k)/primit(1,i,j,k))*Tempsc)
  else
    print*, 'Unsupported EOS'
    stop
  end if

end subroutine

!=======================================================================
!> @brief Upper level wrapper for thermal conduction
!> @details This routine adds the heat conduction, receives the hydro
!>  timestep in seconds, and assumes the primitives and Temp(i,j,k)
!>  arrays are updated
subroutine thermal_conduction()

  use hydro_core, only : calcprim
  implicit none
  real    :: dt_hydro
  real    :: dts, scale, dt_left
  integer :: n,i,j,k, nsteps, nSTb
  logical :: SuperStep = .true.  !  Enable superstepping
  logical :: progressbar = .false.  ! print progress bar within ST blocks

  dt_hydro = dt_CFL*tsc          !  [seconds]
  dt_left  = dt_hydro
  nSTb     = 1

 STblocks : do

  !  get the conduction timescale
    call get_dt_cond(dt_cond)

    !  compute the number of (super) steps needed to asdvance dt_cfl
    if (SuperStep) then
      call ST_steps(dt_left, dt_cond, Nsteps, scale )
    else
      Nsteps = ceiling(dt_left/dt_cond)
      scale = 1.0
      if (NSteps > Max_iter) then
        print*, "Not enough sub-cycles, try increasing MAX_iter, or enabling"//&
                "Super-steps"
      end if
    end if

    if (rank == master) then
      write(tc_log,'(i7,a,es12.4,a,es12.4,a,i6,a,i3)')                         &
                         currentIteration,' | ',dt_left,'  | ', dt_cond,' | ', &
                         Nsteps, ' | ', nSTb
      flush(tc_log)
    end if

    steps : do n=1,Nsteps

      if (SuperStep) then
        dts= dt_cond*scale*substep(n,Nsteps,nu) /Psc/rsc
      else
        dts=dt_hydro/real(Nsteps)               /Psc/rsc
      end if

      !  remaining time to cover (seconds)
      dt_left = dt_left - dts *Psc*rsc

      !  show progress bar
      if (rank == master .and. progressbar) call progress(n,nsteps)

      !  get the heat fluxes
      if (th_cond == TC_ANISOTROPIC) then
        call MHD_heatfluxes()
      end if
      if (th_cond == TC_ISOTROPIC) then
        call heatfluxes()
      end if

      !  update the conserved and primitive vars
      dts = dts /Psc/rsc  ! cgs to code units in the loop
      do k=1,nz
        do j=1,ny
          do i=1,nx
          u(5,i,j,k)=u(5,i,j,k)-dts*( ( f(5,i,j,k) - f(5,i-1,j,k) )/dx &
                                    + ( g(5,i,j,k) - g(5,i,j-1,k) )/dy &
                                    + ( h(5,i,j,k) - h(5,i,j,k-1) )/dz )

          call update_PT(i,j,k)

          end do
        end do
      end do

      !  boundary conditions
      !  (only one layer of u(5,:,:,:) is exchanged )
      call thermal_bounds()

      !  update primitives and Temperature
      !call calcprim(u, primit)

      !=========================================================================
      !  minimal update, only one ghost cell
      !  X (i=0 & i = nx+1)
      do k=0,nz+1
        do j=0,ny+1
          call update_PT(  0  ,j ,k )
          call update_PT(nx+1 ,j ,k )
        end do
      end do
      !  Y (j=0 & j = ny+1)
      do k=0,nz+1
        do i=0,nx+1
          call update_PT(i,   0  ,k )
          call update_PT(i, ny+1 ,k )
        end do
      end do
      !  Z (k=0 & k= nz+1)
      do j=0,ny+1
        do i=0,nx+1
          call update_PT(i ,j , 0    )
          call update_PT(i ,j , nz+1 )
        end do
      end do
      !=========================================================================

    end do steps

    if (rank == master .and. progressbar) call progress(n,nsteps, done=.true.)

    if (rank==master)  print('(a,i4,a,2es12.4,f10.1,a)'),                      &
        ' Finished block of: ', nsteps, ' STs, dt_left/dt_cond ', dt_left,     &
        dt_cond, dt_left/dt_cond, ' remain'

    !  if have finished
    if (abs( dt_left )  <= 0.01*dt_cond ) exit

    nSTb  = nSTb + 1

  end do STblocks

end subroutine thermal_conduction

!=======================================================================

end module thermal_cond
!
!=======================================================================
