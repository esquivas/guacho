module sod_tube
  implicit none
  real :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, pL, pR
contains

  subroutine init_sod()
    implicit none

    rhoL = 1.0    ; rhoR =  0.125
    vxL  = 0.0    ; vxR  =  0.0
    vyL  = 0.0    ; vyR  =  0.0
    vzL  = 0.0    ; vzR  =  0.0
    pL   = 1.0    ; pR   =  0.1

  end subroutine init_sod

  !--------------------------------------------------------------------
  ! Initial conditions for the Sod tube test
  subroutine impose_sod(u)

    use globals,    only : coords, dx
    use parameters, only : neq, nx, nxmin, nxmax, nymin, nymax, nzmin, nzmax,  &
                           nxtot, gamma
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real :: x, rho, px, py, pz, Etot, pas
    integer ::  i,j,k

    do k=nzmin,nzmax
      do j=nymin,nymax
        do i=nxmin,nxmax

          ! Position measured from the centre of the grid
          x=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx

          ! IF LEFT STATE
          if( x <= 0.) then

            rho  = rhoL
            px   = rhoL * vxL
            py   = rhoL * vyL
            pz   = rhoL * vzL
            Etot = 0.5*rhoL*(vxL**2 + vyL**2 +vzL**2) + pL/(gamma-1.0)
            pas  = rhoL

          else !RIGHT STATE

            rho  = rhoR
            px   = rhoR * vxR
            py   = rhoR * vyR
            pz   = rhoR * vzR
            Etot = 0.5*rhoR*(vxR**2 + vyR**2 +vzR**2) + pR/(gamma-1.0)
            pas  = - rhoR
          end if

          !  total density and momenta
          u(1,i,j,k) = rho
          u(2,i,j,k) = px
          u(3,i,j,k) = py
          u(4,i,j,k) = pz
          !  energy
          u(5,i,j,k) = Etot
          !  passive scalar
          u(6,i,j,k) = pas

        end do
      end do
    end do

  end subroutine impose_sod
  !--------------------------------------------------------------------

end module sod_tube
