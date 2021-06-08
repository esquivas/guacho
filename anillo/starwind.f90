!=======================================================================
!> @file exoplanet.f90
!> @brief Exoplanet problem module
!> @author M. Schneiter, C. Villarreal  D'Angelo, A. Esquivel
!> @date 4/May/2016

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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief starwind module
!> @details Problem Module for starwind (multiple stars with wind)

module starwind

  use parameters
  implicit none

  real :: dsw     !< Stellar Wind Density
  real :: RPW     !< source(s) radius
  real :: TPW     !< source(s) wind temperature
  real :: VPW     !< source(s) wind velocity
  real :: dpw     ! Planetary wind density
  real :: y0p

contains

!=======================================================================
!> @brief Module initialization
!> @details Here the parameters of the Star are initialized, and scaled
!> to code units
subroutine init_starwind()

  use constants, only : Msun, Yr, Ggrav, pi, mjup, au, day,Rjup, rsun
  implicit none
  real                   :: amdot        ! mdot_star (MSUN/yr)
  real                   :: ampdot       ! mdot_planet (g/s)
  integer                :: i_star


  !------------------       STAR PARAMETERS       --------------------
  !               ALL N STARS  with the same parameters
  !------------------------------------------------------------------------
  !
  AMPDOT= 1.e-6*msun/yr                  ! Star Mass Loss rate (g/s)
  TPW   = 1.e5                           ! Temperatura de la estrella
  RPW   = 1.e16                          ! Radio del viento de la estrella (cm)
  vpw   = 1000.e5                        ! Velocidad del viento de la estrella (cm/s)
  dpw=((AMPDOT/RPW)/(4*pi*RPW*VPW))      ! Densidad del viento de la estrella
  y0p = 1.e-4                            ! Densidad de neutros del viento

  ! -------------     ORBITAL PARAMETERS     ---------------
  !      For a star orbiting a BH with 4.e6 M_sun  Mass
  !         at a distance d=1.5e17 cm = 10026.88*AU
  ! rorb = 10026.88*AU
  ! torb = 502.736*365.*day

  ! CHANGE TO CODE UNITS
  dpw=dpw/rhosc
  vpw=vpw/sqrt(vsc2)
  Tpw=Tpw/Tempsc
  Rpw=Rpw/rsc

  !rorb=rorb/rsc
  !torb=torb/tsc
  !omegap=2.*pi/torb
  !
  !

end subroutine init_starwind

!=======================================================================
!> @brief Inject sources of wind
!> @details Imposes the sources of wond from the star and planet
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!> conserver variables
!> @param real [time] time : current integration timr
!--------------------------------------------------------------------
subroutine impose_starwind(u,time)
  use constants, only : pi
  use globals, only : coords, dx, dy, dz
  use nbody,   only : x_star,y_star,z_star,vx_star,vy_star,vz_star,N_star
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real              :: x, y, z
  real  :: radp, xp, yp, zp, vxorb, vyorb, vzorb
  real  :: xpl, ypl, zpl, velx,vely,velz
  !real ,dimension(N_star)  ::  rorb,torb,omegap ,phi
  integer :: i_star
  real :: dens
  real :: time_seconds

#ifdef BFIELD

  real :: cpi

#endif

  integer ::  i,j,k
  integer :: ii, jj,kk
  real*8  :: radio1, radio2, dif_radio

  !  integer :: matrix(2,3) !two dimensional real array
  time_seconds = time * tsc

!-------------------------------!
!  Para saber si se tienen dos  !
!  estrellas en la misma celda. !
! Si es asi, entonces se impone !
!     un viento del "doble".    !
!-------------------------------!
!!!kk =0
!!!  do ii=2,N_star
!!!     radio1 = sqrt(x_star(ii)**2 + y_star(ii)**2 + z_star(ii)**2)
!!!     do jj=2,N_star
!!!        if (ii.ne.jj) then
!!!        radio2 = sqrt(x_star(jj)**2 + y_star(jj)**2 + z_star(jj)**2)
!!!        dif_radio = abs(radio1-radio2)
!!!        if (dif_radio.le.(2.*RPW*rsc)) then
!!!           print *, 'Doble Viento...', ii, jj, dif_radio, x_star(ii),rpw*rsc
!!!           kk = kk +1
!!!        endif
!!!        endif
!!!     enddo
!!!  enddo
!!!  print *, 'aqui acabo ...', kk
!!!
!!!stop
!

  do k=nzmin,nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax

        do i_star=2,N_star

          xp = x_star(i_star) /rsc  !  posicion en unidades de codigo
          yp = y_star(i_star) /rsc  !  posicion en unidades de codigo
          zp = z_star(i_star) /rsc  !  posicion en unidades de codigo

          vxorb = vx_star(i_star) / vsc  !  orbital velocity code units
          vyorb = vy_star(i_star) / vsc  !  orbital velocity code units
          vzorb = vz_star(i_star) / vsc  !  orbital velocity code units

          ! Position measured from the centre of the grid (BH)
          x=(float(i+coords(0)*nx-nxtot/2) - 0.5)*dx
          y=(float(j+coords(1)*ny-nytot/2) - 0.5)*dy
          z=(float(k+coords(2)*nz-nztot/2) - 0.5)*dz

          ! distance measured from the centre of the wind source
          xpl = x - xp
          ypl = y - yp
          zpl = z - zp
          !
          radp = sqrt( xpl**2+ypl**2+zpl**2 )

          ! IF inside the star
          if(radp <= rpw ) then
            if(radp == 0. ) radp=dx*0.10
            !
            velx = vxorb + vpw * xpl/radp
            vely = vyorb + vpw * ypl/radp
            velz = vzorb + vpw * zpl/radp
            dens=dpw
            !
            !   total density and momenta
            u(1,i,j,k) = dens      + u(1,i,j,k)
            u(2,i,j,k) = dens*velx + u(2,i,j,k)
            u(3,i,j,k) = dens*vely + u(3,i,j,k)
            u(4,i,j,k) = dens*velz + u(4,i,j,k)
            !
            !  density of neutrals
            u(neqdyn+1,i,j,k)=y0p*dens + u(neqdyn+1,i,j,k)
            !
            !   passive scalar (h-hot, c-cold, i-ionized, n-neutral)
            u(neqdyn+2,i,j,k)= -dens   + u(neqdyn+2,i,j,k) ! passive scalar

            u(5,i,j,k)=0.5*dens*( velx**2+vely**2+velz**2 ) &
            + cv*Tpw*(2.*dens-y0p*dens) + u(5,i,j,k)

          end if
        
        enddo

      end do
    end do
  end do


end subroutine impose_starwind

!=======================================================================

end module starwind

!=======================================================================
