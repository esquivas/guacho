module two_winds

  use parameters
  implicit none

  integer, parameter :: Nsources = 2  ! Number of sources
  real :: dw(Nsources)     !< Stellar Wind Density
  real :: RW(Nsources)     !< source radius
  real :: TW(Nsources)     !< source wind temperature
  real :: VW(Nsources)     !< source wind velocity
  real :: y0(Nsources)     !< ionization fraction

  real :: mass(Nsources)   !< Mass of the sources
  real :: xp(Nsources)     !< X position of the sources
  real :: yp(Nsources)     !< Y position of the sources
  real :: zp(Nsources)     !< Z position of the sources

  real :: torb    !<  orbital period
  real :: rorb    !<  orbital radius
 contains

!=======================================================================
subroutine init_two_winds()

  use constants, only : Msun, Yr, Ggrav, pi, Mjup, au, hr, Rsun,yr
  use globals, only :rank, dx
  use parameters
  implicit none
  real :: mdot(Nsources)

  !  First source (star centered in domain)
  mass(1) = 90.0*Msun
  mdot(1) = 4.8e-4*Msun/yr   !  Mass Loss rate (g s^-1)
  TW(1)   = 3.5e4            !  wind Temperature (K)
  RW(1)   = 60.0*Rsun
  vw(1)   = 420.0e5          !  wind velocity
  dw(1)   =  (mdot(1)/rw(1))/(4.0*pi*rw(1)*vw(1))
  y0(1)   = 1e-4

  !  Second source (jupiter like planet in close orbit)
  mass(2) = 30.0*Msun
  mdot(2) = 1.45e-5          !  Mass Loss rate (g s^-1)
  TW(2)   = 3.5E4             !  wind Temperature (K)
  RW(2)   = 30.0*Rsun
  vw(2)   = 3000.0E5            !  wind velocity
  dw(2)   = (mdot(2)/rw(2))/(4.0*pi*rw(2)*vw(2))
  y0(2)   = 0.1


  print*, RW(:)

  !  ORBITAL PARAMETERS
  rorb= 1.45*AU
  torb= 5.54*yr

  ! change to code units
  dw(:) = dw(:)/rhosc
  vw(:) = vw(:)/vsc
  Tw(:) = Tw(:)/Tempsc
  Rw(:) = Rw(:)/rsc

  rorb=rorb/rsc
  torb=torb/tsc

  !  initial positions respect to the grid center
  xp(1) = 0.0  ; xp(2) = Rorb
  yp(1) = 0.0  ; yp(2) = 0.0
  yp(1) = 0.0  ; zp(2) = 0.0

end subroutine init_two_winds

!=======================================================================
subroutine impose_two_winds(u,time)

  use constants, only : pi
  use globals, only : coords, dx, dy, dz, rank
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x, y, z
  real :: velx, vely, velz, rads, dens, radp
  real :: vxorb, vzorb, vyorb, omega
  integer ::  i, j, k

  omega=2.*pi/torb

  !  Planet orbit
  xp(2) = Rorb*COS(omega*time)
  zp(2) = Rorb*SIN(omega*time)

  !Orbital planetary velocity (moves in the xz-plane)
  vxorb= -omega*Rorb*sin(omega*time)
  vzorb=  omega*Rorb*cos(omega*time)
  vyorb=0.

  do k=nzmin,nzmax
    do j=nymin,nymax
      do i=nxmin,nxmax

        ! Position measured from the centre of the grid
        x=(real(i+coords(0)*nx-nxtot/2)-0.5)*dx
        y=(real(j+coords(1)*ny-nytot/2)-0.5)*dy
        z=(real(k+coords(2)*nz-nztot/2)-0.5)*dz

        ! Distance from the centre of the first source
        rads=sqrt( (x-xp(1))**2 + (y-yp(1))**2 + (z-zp(1))**2 )

        ! Distance from the centre of the second source
        radp=sqrt( (x-xp(2))**2 + (y-yp(2))**2 + (z-zp(2))**2 )

        ! IF INSIDE THE STAR
        if( rads <= rw(1) ) then

          if(rads == 0.) rads=dx*0.10
          velx = vw(1)* (x-xp(1)) /rads
          vely = vw(1)* (y-yp(1)) /rads
          velz = vw(1)* (z-zp(1)) /rads
          dens = dw(1)*rw(1)**2/rads**2

          !   total density and momenta
          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz

          !  density of neutrals
          u(neqdyn+1,i,j,k) =  y0(1)*dens
          !  passive scalar (tag) for stellar material
          u(neqdyn+2,i,j,k)= dens

          ! total energy
          u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) +             &
                       cv * (2.0*u(1,i,j,k) - u(neqdyn+1,i,j,k) ) * Tw(1)


          ! IF INSIDE THE PLANET
        else if(radp <= rw(2) ) then

          if(radp == 0.) radp=dx*0.10

          velx = vw(2)* (x-xp(2)) /radp + vxorb
          vely = vw(2)* (y-yp(2)) /radp + vyorb
          velz = vw(2)* (z-zp(2)) /radp + vzorb
          dens = dw(2)*rw(2)**2/radp**2

          !   total density and momenta
          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz

          !  density of neutrals
          u(neqdyn+1,i,j,k) = y0(2)*dens

          !   passive scalar (tag) for planetary material
          u(neqdyn+2,i,j,k)= -dens

           ! total energy
          u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) +             &
                       cv * (2.0*u(1,i,j,k) - u(neqdyn+1,i,j,k) ) * Tw(2)

        end if

      end do
    end do
  end do

end subroutine impose_two_winds

!=======================================================================

end module two_winds
