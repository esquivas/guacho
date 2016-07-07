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
    use globals, only : coords, dx, dy, primit0
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
           primit0(1,i,j,k) = rho
           primit0(2,i,j,k) = 0.
           primit0(3,i,j,k) = 0.
           primit0(4,i,j,k) = 0.
           primit0(5,i,j,k) = p
           
           u(1,i,j,k) = 0.
           u(2,i,j,k) = 0.
           u(3,i,j,k) = 0.
           u(4,i,j,k) = 0.
           
!            if (r.le.0.125) then
           if (r.lt.0.1) then
             u(5,i,j,k) = cv*p1
           else
             u(5,i,j,k) = 0.
           end if
           
          !   total energy

#ifdef BFIELD

          primit0(6,i,j,k) =  bx
          primit0(7,i,j,k) =  by
          primit0(8,i,j,k) =  bz

          u(6,i,j,k) = 0.
          u(7,i,j,k) = 0.
          u(8,i,j,k) = 0.

! #else
!           u(5,i,j,k)=0.5*rho*(vx**2+vy**2+vz**2)+cv*p
#endif
        end do
       end do
    end do

  end subroutine impose_blast
  !--------------------------------------------------------------------

end module blast
