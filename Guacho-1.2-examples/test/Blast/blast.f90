!=======================================================================
!  Implements the BLAST test
!=======================================================================
module Blast
  use parameters
  use constants
  implicit none
  real :: rho, vx, vy, vz, p
  real :: Bx, By, Bz
  real , parameter :: twopi=2.*pi

contains

  subroutine init_blast()
    implicit none
  end subroutine init_blast

  !--------------------------------------------------------------------
  ! Initial conditions for the blast
  subroutine impose_blast(u)

    use globals, only : coords, dx, dy, dz
    use parameters, only : nxtot,gamma,xmax,ymax
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real :: x, y, z, r, rin
    real :: theta,phi,b0
    integer ::  i,j,k

    theta = 45.0
    phi   = 45.0
    b0    = 1.0

    rin = 0.1
    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax

             ! Position measured from the centre of the grid (star)
             x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
             y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
             z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz

             r = sqrt(x*x+y*y+z*z)

             rho = 1.0
             p = 0.1

             if(r < rin)then
               p = 100.0*p
             endif

             vx = 0.0
             vy = 0.0
             vz = 0.0

#ifdef MHD 
             bx = b0*sin(theta*pi/180.0)*cos(phi*pi/180.0)
             by = b0*sin(theta*pi/180.0)*sin(phi*pi/180.0)
             bz = b0*cos(theta*pi/180.0)
#endif                              

             !   total density and momenta
             u(1,i,j,k) = rho
             u(2,i,j,k) = rho*vx
             u(3,i,j,k) = rho*vy
             u(4,i,j,k) = rho*vz
             !   total energy
             u(5,i,j,k) =  0.5*rho*(vx**2+vy**2+vz**2) &
                           + p/(gamma - 1.0)

#ifdef MHD 
             u(5,i,j,k) =  u(5,i,j,k) + 0.5*(Bx**2+By**2+Bz**2)
             u(6,i,j,k) =  Bx
             u(7,i,j,k) =  By
             u(8,i,j,k) =  BZ
#endif

          end do
       end do
    end do

  end subroutine impose_blast
  !--------------------------------------------------------------------

end module Blast


