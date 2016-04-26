!=======================================================================
!  Implements the Current Sheet test
!=======================================================================
module CurrentSheet
  use parameters
  use constants
  implicit none
  real :: rho, vx, vy, vz, p
  real :: Bx, By, Bz
  real , parameter :: twopi=2.*pi

contains

  subroutine init_currentsheet()
    implicit none
  end subroutine init_currentsheet

  !--------------------------------------------------------------------
  ! Initial conditions for the current sheet
  subroutine impose_currentsheet(u)

    use globals, only : coords, dx, dy, dz
    use parameters, only : nxtot,gamma,xmax,ymax,rsc
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real :: beta, A
    real :: x,y
    integer ::  i,j,k

    beta = 0.1
    A    = 0.1

    do k=nzmin,nzmax
       do j=nymin,nymax
          do i=nxmin,nxmax

             ! Position measured from the centre of the grid (star)
             x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx*rsc
             y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy*rsc

             rho = 1.0
             p = beta*0.5

             if(abs(x) > 0.25)then
               by = sqrt(4.0*pi)
             else
               by = -sqrt(4.0*pi)
             endif

             vx = A*sin(2.0*pi*y)
             vy = 0.0
             vz = 0.0

             !   total density and momenta
             u(1,i,j,k) = rho
             u(2,i,j,k) = rho*vx
             u(3,i,j,k) = rho*vy
             u(4,i,j,k) = rho*vz
             !   total energy
             u(5,i,j,k) =  0.5*rho*(vx**2+vy**2+vz**2) &
                           + p/(gamma - 1.0)

             u(5,i,j,k) =  u(5,i,j,k) + 0.5*(bx**2+by**2+bz**2)
             u(6,i,j,k) =  Bx
             u(7,i,j,k) =  By
             u(8,i,j,k) =  BZ

          end do
       end do
    end do

  end subroutine impose_currentsheet
  !--------------------------------------------------------------------

end module CurrentSheet


