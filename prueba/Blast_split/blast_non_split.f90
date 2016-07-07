!=======================================================================
!  Implements the STRONG BLAST WAVE test
!=======================================================================
module blast
  use parameters
  use constants
  implicit none
  real :: rho, vx, vy, vz, p, p1
  real :: Bx, By, Bz
  real, parameter :: r0 = 0.125

contains

  subroutine init_blast()
    implicit none
    
    rho = 1.
    p   = 0.1!1.
    p1  = 10.!100. !100 veces mas grande
        
  end subroutine init_blast

  !--------------------------------------------------------------------
  ! Initial conditions for the blast 
  subroutine impose_blast(u)
    use globals, only : coords, dx, dy
    use parameters, only : nxtot, nytot
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real :: xp, yp, r
    integer ::  i,j,k
    
    do i=nxmin,nxmax
      do j=nymin,nymax
        do k=nzmin,nzmax

          ! Position measured from the center of the grid
          xp=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
          yp=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
          r = sqrt(xp**2 +yp**2)
          
          vx = 0.
          vy = 0.
          vz = 0.
#ifdef BFIELD
          bx = sqrt(2*pi)!10
          by = sqrt(2*pi)!0.
          bz = 0.
#endif                              

          !   total density and momenta
           
           u(1,i,j,k) = rho
           u(2,i,j,k) = 0.
           u(3,i,j,k) = 0.
           u(4,i,j,k) = 0.
           
          !   total energy

#ifdef BFIELD

          u(6,i,j,k) = bx
          u(7,i,j,k) = by
          u(8,i,j,k) = bz
          
          if (r.lt.0.1) then
            u(5,i,j,k) = cv*(p1+p) + 0.5*(bx**2+by**2+bz**2) + 0.5*rho*(vx**2+vy**2+vz**2)
          else
            u(5,i,j,k) = cv*p + 0.5*(bx**2+by**2+bz**2) + 0.5*rho*(vx**2+vy**2+vz**2)
          end if

#else
          if (r.lt.0.1) then
            u(5,i,j,k)=0.5*rho*(vx**2+vy**2+vz**2)+cv*(p+p1)
          else
            u(5,i,j,k)=0.5*rho*(vx**2+vy**2+vz**2)+cv*p
          end if
#endif
        end do
       end do
    end do

  end subroutine impose_blast
  !--------------------------------------------------------------------

end module blast
