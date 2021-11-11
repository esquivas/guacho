module shock_tube_rel
  implicit none
  real :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, pL, pR
contains

  subroutine init_rs_p1()
    implicit none

    rhoL = 1.0    ; rhoR =  1.0
    vxL  = 0.9    ; vxR  =  0.0
    vyL  = 0.0    ; vyR  =  0.0
    vzL  = 0.0    ; vzR  =  0.0
    pL   = 1.0    ; pR   =  10.0

  end subroutine init_rs_p1

  !--------------------------------------------------------------------
  ! Initial conditions for the relativistic shock tube (P1 of Mignone & Bodo 05)
  subroutine impose_rs(u)

    use globals,    only : coords, dx
    use constants,  only : EOS_REL_IDEAL, EOS_REL_TM
    use parameters, only : neq, nx, nxmin, nxmax, nymin, nymax, nzmin, nzmax,  &
                           nxtot, eq_of_state, gamma
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real :: x, gamma_rel, D, mx, my, mz, Etot, h
    integer ::  i,j,k

    do k=nzmin,nzmax
      do j=nymin,nymax
        do i=nxmin,nxmax

          ! Position measured from the centre of the grid
          x=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx

          ! IF LEFT STATE
          if( x <= 0.) then

            if (eq_of_state == EOS_rel_ideal ) then
              h = 1.0 + (pL/rhoL)*gamma/( gamma- 1.0 )
            end if

            gamma_rel = 1.0 / sqrt( 1.0 - vxL**2 - vyL**2 - vzL**2 )

            D    = gamma_rel * rhoL
            mx   = D * h * gamma_rel * vxL
            my   = D * h * gamma_rel * vyL
            mz   = D * h * gamma_rel * vzL
            Etot = D * h * gamma_rel - pL

          else !RIGHT STATE

            if (eq_of_state == EOS_rel_ideal ) then
              h = 1.0 + (pR/rhoR)*gamma/( gamma- 1.0 )
            end if

            gamma_rel = 1.0 / sqrt( 1.0 - vxR**2 - vyR**2 - vzR**2 )

            D    = gamma_rel * rhoR
            mx   = D * h *gamma_rel * vxR
            my   = D * h *gamma_rel * vyR
            mz   = D * h *gamma_rel * vzR
            Etot = D * h *gamma_rel - pR

          end if

          !  total density and momenta
          u(1,i,j,k) = D
          u(2,i,j,k) = mx
          u(3,i,j,k) = my
          u(4,i,j,k) = mz
          !  energy
          u(5,i,j,k) = Etot
          !  passive scalar
          u(6,i,j,k) = D

        end do
      end do
    end do

  end subroutine impose_rs
  !--------------------------------------------------------------------

end module shocK_tube_rel
