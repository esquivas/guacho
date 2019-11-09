!=======================================================================
!  Implements the ORZAG-TANG VORTEX test
!=======================================================================
module OrzagTang2
  use parameters
  use constants
  implicit none
  real :: rho, vx, vy, vz, p
  real :: Bx, By, Bz
  real , parameter :: twopi=2.*pi

contains

  subroutine init_ot()
    implicit none
    
    rho = 25./(36.*pi)
    p   = 5. /(12.*pi)
        
  end subroutine init_ot

  !--------------------------------------------------------------------
  ! Initial conditions for the orzag_tang
  subroutine impose_ot(u)
    use globals, only : coords, dx, dy
    use parameters, only : nxtot
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real :: x, y
    integer ::  i,j,k
    
    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax

             ! Position measured from the bottom corner of the grid
             x=(float(i+coords(0)*nx)+0.5)*dx*rsc !WHY +0.5?
             y=(float(j+coords(1)*ny)+0.5)*dy*rsc !ANS: Position of cell center

                vx = -sin( y*twopi )
                vy =  sin( x*twopi )
                vz = 0.
#ifdef MHD 
                bx = -sin( y*twopi  )/sqrt(4*pi)
                by =  sin(2.*x*twopi)/sqrt(4*pi)
                bz = 0.
#endif                              

                !   total density and momenta
                u(1,i,j,k) = rho
                u(2,i,j,k) = rho*vx
                u(3,i,j,k) = rho*vy
                u(4,i,j,k) = rho*vz
                !   total energy

#ifdef MHD 
                u(5,i,j,k)=0.5*rho*(vx**2+vy**2+vz**2)+cv*p  &
                               +0.5*(Bx**2+By**2+Bz**2)
                u(6,i,j,k) =  Bx
                u(7,i,j,k) =  By
                u(8,i,j,k) =  BZ

#else
                u(5,i,j,k)=0.5*rho*(vx**2+vy**2+vz**2)+cv*p
#endif
          end do
       end do
    end do

  end subroutine impose_ot
  !--------------------------------------------------------------------

end module OrzagTang2
