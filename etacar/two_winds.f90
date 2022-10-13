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
   
  real :: vxorb(Nsources)   !< X orbital velocity elliptical orbit
  real :: vyorb(Nsources)   !< Y orbital velocity elliptical orbit
  real :: vzorb(Nsources)   !< Z orbital velocity elliptical orbit

  real :: e                 !<  eccentricity of the orbit
  real :: a                 !<  distance to compute position vector r 
  real :: r                 !<  position vector
  real :: rm(Nsources)      !<  reduced mass
  real :: torb              !<  orbital period
  real :: rorbT             !<  max separation between stars elliptical orbit      
  real :: rorb(Nsources)    !<  orbital radius
  real :: r0(Nsources)      !<  difference between the center of the ellipse
                            !<  and the center of mass of the system
  real :: omega             !<  orbital angular velocity 

 
 contains

!=======================================================================
subroutine init_twowinds()

  use constants, only : Msun, Yr, Ggrav, pi, Mjup, au, day,Rjup
  use globals, only :rank, dx
  use parameters
  implicit none
  real :: mdot(Nsources)  

  !  ORBITAL PARAMETERS
  e      = 0.9
  a      = 8.13*AU/rsc        ! Para garantizar una distancia entre estrellas = rorbT
  rorbT  = (15.45/2)*AU
  torb   =  5.54*yr/tsc
  omega  = 2.*pi/torb 

  !  First source (star centered in domain)
  mass(1)  = 90.*Msun 
  mdot(1)  = 4.8E-4*Msun/yr   !  Mass Loss rate (g s^-1)
  TW(1)    = 3.5E4            !  wind Temperature (K)
  RW(1)    = 0.5*AU 
  vw(1)    = 420.e5           !  wind velocity
  dw(1)    = (mdot(1)/rw(1))/(4.0*pi*rw(1)*vw(1))
  y0(1)    = 1e-4
  r0(1)    = -4.355*AU
  rorb(1)  = 3.86*AU
  

  !  Second source (jupiter like planet in close orbit)
  mass(2)  = 30.*Msun
  mdot(2)  = 1.4e-5*Msun/yr   !  Mass Loss rate (g s^-1)
  TW(2)    = 3.5E4            !  wind Temperature (K)
  RW(2)    = 0.5*AU
  vw(2)    = 3000.e5          !  wind velocity
  dw(2)    = (mdot(2)/rw(2))/(4.0*pi*rw(2)*vw(2))
  y0(2)    = 0.5
  r0(2)    = 3.36*AU
  rorb(2)  = 11.58*AU

  ! Massas reducidas
  rm(1) = -mass(2)/(mass(1)+mass(2))
  rm(2) =  mass(1)/(mass(1)+mass(2))

  ! change to code units
  dw(:) = dw(:)/rhosc
  vw(:) = vw(:)/vsc
  Tw(:) = Tw(:)/Tempsc
  Rw(:) = Rw(:)/rsc
  R0(:) = R0(:)/rsc
  Rorb(:) = Rorb(:)/rsc
  
    
  !  initial positions respect to the grid center i.e. center of mass of
  !                                                    elliptical orbit
  !xp(1) = -3.86*AU/rsc  ; xp(2) = 11.58*AU/rsc
  !yp(1) = 0.0           ; yp(2) = 0.0  
  !zp(1) = 9.0           ; zp(2) = 0.0

end subroutine init_twowinds

!=======================================================================
subroutine impose_winds(u,time)

  use constants, only : pi
  use globals, only : coords, dx, dy, dz, rank
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x, y, z
  real :: velx, vely, velz, rads, dens, radp
  integer ::  i, j, k

  
 r  = a*(1.-e**2) / (1.-(e*cos(omega*time)))

! Orbita excentrica de estrella masiva
  xp(1) = rm(1) * r *(cos(omega*time))
  yp(1) = 0.
  zp(1) = rm(1) * r *(sin(omega*time))

  xp(2) = rm(2) * r *(cos(omega*time))
  yp(2) = 0.
  zp(2) = rm(2) * r *(sin(omega*time))



  !Orbital planetary velocity (moves in the xz-plane)
  !vxorb= -omega*Rorb*sin(omega*time)
  !vzorb=  omega*Rorb*cos(omega*time)
 
  vxorb(1) =  - omega*rm(1)*r*sin(omega*time)
  vzorb(1) =    omega*rm(1)*r*cos(omega*time)
  vyorb(1) = 0.

  vxorb(2) =  - omega*rm(2)*r*sin(omega*time)
  vzorb(2) =    omega*rm(2)*r*cos(omega*time)
  vyorb(2) = 0.

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
          velx = vw(1)* (x-xp(1)) /rads + vxorb(1)
          vely = vw(1)* (y-yp(1)) /rads + vyorb(1)
          velz = vw(1)* (z-zp(1)) /rads + vzorb(1)
          dens = dw(1)*rw(1)**2/rads**2

          !   total density and momenta
          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz

          !  density of neutrals
          u(neqdyn+1,i,j,k) =  y0(1)*dens
          !  passive scalar (tag) for stellar material
          u(neqdyn+2,i,j,k)= 1000*dens

          ! total energy
          u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) +             &
                       cv * (2.0*u(1,i,j,k) - u(neqdyn+1,i,j,k) ) * Tw(1)
          

          ! IF INSIDE THE PLANET
        else if(radp <= rw(2) ) then

          if(radp == 0.) radp=dx*0.10

          velx = vw(2)* (x-xp(2)) /radp + vxorb(2)
          vely = vw(2)* (y-yp(2)) /radp + vyorb(2)
          velz = vw(2)* (z-zp(2)) /radp + vzorb(2)
          dens = dw(2)*rw(2)**2/radp**2

          !   total density and momenta
          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz

          !  density of neutrals
          u(neqdyn+1,i,j,k) = y0(2)*dens

          !   passive scalar (tag) for planetary material
          u(neqdyn+2,i,j,k)= -1000*dens

           ! total energy
          u(5,i,j,k) = 0.5*dens*(velx**2+vely**2+velz**2) +             &
                       cv * (2.0*u(1,i,j,k) - u(neqdyn+1,i,j,k) ) * Tw(2)

        end if

      end do
    end do
  end do

end subroutine impose_winds

!=======================================================================

end module two_winds

