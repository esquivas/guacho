!=======================================================================
!  Implements the KELVIN-HELMHOLTZ INSTABILITY test
!=======================================================================
module RayleighTaylor
  use parameters
  use constants
  implicit none
  real :: rho, vx, vy, vz, p
  real :: Bx, By, Bz
  real , parameter :: twopi=2.*pi

contains

  subroutine init_rt()
    implicit none
    
        
  end subroutine init_rt

  !--------------------------------------------------------------------
  ! Initial conditions for the rayleigh_taylor
  subroutine impose_rt(u)
    use globals, only : coords, dx, dy
    use parameters, only : nxtot,gamma,xmax,ymax
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real :: x, y
    integer ::  i,j,k
    integer,parameter :: seed = 86456
    real :: r
    call srand(seed)
    !CALL init_random_seed()         ! see example of RANDOM_SEED
    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax

             ! Position measured from the centre of the grid (star)
             x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
             y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy

             rho = 1.0

             vx = 0.0
             vy = 0.01*(1.0+cos(2.*pi*x/xmax))*(1.0+cos(2.0*pi*y/ymax))/4.0
             vz = 0.0

#ifdef MHD 
             bx = 0.5/sqrt(4*pi)
             by = 0.
             bz = 0.
#endif                              

             if (y > 0.0) then
               rho = 2.0
             endif

             !   total density and momenta
             u(1,i,j,k) = rho
             u(2,i,j,k) = rho*vx
             u(3,i,j,k) = rho*vy
             u(4,i,j,k) = rho*vz
             !   total energy

#ifdef MHD 
             u(5,i,j,k) =  0.5*rho*(vx**2+vy**2+vz**2)+cv*p+0.5*(Bx**2+By**2+Bz**2)
             u(6,i,j,k) =  Bx
             u(7,i,j,k) =  By
             u(8,i,j,k) =  BZ

#else
             if (y > 0.0) then
               u(5,i,j,k) = (1.0/gamma - 0.2*y)/(gamma - 1.0)
             else
               u(5,i,j,k) = (1.0/gamma - 0.1*y)/(gamma - 1.0)
             endif
             u(5,i,j,k) = u(5,i,j,k) + 0.5*rho*vy**2

#endif
          end do
       end do
    end do

  end subroutine impose_rt
  !--------------------------------------------------------------------

end module RayleighTaylor
