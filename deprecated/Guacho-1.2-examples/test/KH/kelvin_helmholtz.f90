!=======================================================================
!  Implements the KELVIN-HELMHOLTZ INSTABILITY test
!=======================================================================
module KevinHelmholtz
  use parameters
  use constants
  implicit none
  real :: rho, vx, vy, vz, p
  real :: Bx, By, Bz
  real , parameter :: twopi=2.*pi

contains

  subroutine init_kh()
    implicit none
    
    p   = 2.5
        
  end subroutine init_kh

  !--------------------------------------------------------------------
  ! Initial conditions for the kelvin_helmholtz
  subroutine impose_kh(u)
    use globals, only : coords, dx, dy
    use parameters, only : nxtot,gamma
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real :: x, y
    integer ::  i,j,k
    integer,parameter :: seed = 86456
    real :: vflow = 0.5
    real :: r
    call srand(seed)
    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax

             ! Position measured from the centre of the grid (star)
             x=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx
             y=(real(j+coords(1)*ny-nytot/2)+0.5)*dy

             rho = 2.0

             vx = vflow + 0.01*(rand() - 0.5)
             vy = 0.01*(rand() - 0.5)
             vz = 0.0

#ifdef MHD 
             bx = 0.5*sqrt(4*pi)
             by = 0.
             bz = 0.
#endif                              

             if (abs(y) > 0.25) then
               vx  = -vx
               rho = 0.5*rho
             endif

             !   total density and momenta
             u(1,i,j,k) = rho
             u(2,i,j,k) = rho*vx
             u(3,i,j,k) = rho*vy
             u(4,i,j,k) = rho*vz

             !   total energy
             u(5,i,j,k) = 0.5*rho*(vx**2+vy**2+vz**2) + p/(gamma-1.0)
#ifdef MHD 
             u(5,i,j,k) = u(5,i,j,k) + 0.5*(Bx**2+By**2+Bz**2)
             u(6,i,j,k) = bx
             u(7,i,j,k) = by
             u(8,i,j,k) = bz
#endif
          end do
       end do
    end do

  end subroutine impose_kh
  !--------------------------------------------------------------------

end module KevinHelmholtz
