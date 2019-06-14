module two_fluid

  use parameters,  only : neq,  nxmin, nymin, nzmin, &
  nxmax, nymax, nzmax, nx, ny, nz , tsc, rhosc

  implicit none

  ! Look for the value of alpha (where did we take this from?)
  real,parameter :: alphac = 6.02e14*rhosc*tsc! 5e-2

contains

  !> @brief Twofluid interaction term
  !> @brief from Smith & Sakai (2008)

  subroutine get_TF_sources(pp, pn, tp, tn, s)

    use cooling_H, only: alpha, colf
    implicit none
    real, intent(in)  :: pp(neq), pn(neq), tp, tn
    real, intent(out) :: s(neq)
    real :: driftx, drifty, driftz,ci,ar
    !  update the source terms

    ci=colf(tn)
    ar=alpha(tp)
    
    driftx = pn(2)-pp(2)
    drifty = pn(3)-pp(3)
    driftz = pn(4)-pp(4)

    ! reset sources
    s(:) = 0.

    !  mass 14/06/2019
    s(1) = - pp(1)* (ar * pp(1) - ci * pn(1) )
    
    !  momenta
    s(2) = alphac * pp(1)*pn(1)*driftx - pp(1)*(ar*pp(1)*pp(2)-ci*pn(1)*pn(2))
    s(3) = alphac * pp(1)*pn(1)*drifty - pp(1)*(ar*pp(1)*pp(3)-ci*pn(1)*pn(3))
    s(4) = alphac * pp(1)*pn(1)*driftz - pp(1)*(ar*pp(1)*pp(4)-ci*pn(1)*pn(4))
    s(5) = alphac * pp(1)*pn(2)* (driftx*pp(2)+drifty*pp(3)+driftz*pp(4) )

  end subroutine get_TF_sources

  ! obtains R
  subroutine get_TF_R(pp, pn,tp,tn,dt, R)
    use cooling_H, only : alpha, colf
    implicit none
    real, intent(in)  :: pp(neq), pn(neq), tp, tn, dt
    real, intent(out) :: R(3)
    real :: den, ci, ar

    ci=colf(tn)
    ar=alpha(tp)
    
    den =  1+pp(1)*(ar*dt+ci*dt) +alphac*dt*(pp(1)+pn(1)) 
    
    R(1) =(pp(1)*(ci*dt*pn(1)*pn(2)-ar*dt*pp(1)*pp(2) ) + alphac*dt * pp(1)*pn(1) * ( pn(2)-pp(2) )) / den
    R(2) =(pp(1)*(ci*dt*pn(1)*pn(3)-ar*dt*pp(1)*pp(3) ) + alphac*dt * pp(1)*pn(1) * ( pn(3)-pp(3) )) / den
    R(3) =(pp(1)*(ci*dt*pn(1)*pn(4)-ar*dt*pp(1)*pp(4) ) + alphac*dt * pp(1)*pn(1) * ( pn(4)-pp(4) )) / den


  end subroutine get_TF_R


  subroutine update_TF(u,un,up,upn,pp,pn,dt)

    use globals, only : temp, tempn

    implicit none
    real, intent(in)     :: u  (neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(in)     :: un (neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(inout)  :: up (neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(inout)  :: upn(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(in)     :: pp (neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(in)     :: pn (neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(in)     :: dt
    real                 :: R(3)
    integer  :: i, j, k

    do k=1,nz
      do j=1,ny
        do i=1,nx

          !  Gets the source terms
          call get_TF_R(pp(:,i,j,k), pn(:,i,j,k),temp(i,j,k), tempn(i,j,k),dt,R )

          !  update momenta
          up(2,i,j,k) = up(2,i,j,k) + R(1)
          up(3,i,j,k) = up(3,i,j,k) + R(2)
          up(4,i,j,k) = up(4,i,j,k) + R(3)

          upn(2,i,j,k) = upn(2,i,j,k) - R(1)
          upn(3,i,j,k) = upn(3,i,j,k) - R(2)
          upn(4,i,j,k) = upn(4,i,j,k) - R(3)

          !  update enery
          up(5,i,j,k) = up(5,i,j,k)   + pp(2,i,j,k)*R(1) &
                                      + pp(3,i,j,k)*R(2) &
                                      + pp(4,i,j,k)*R(3)

          upn(5,i,j,k) = upn(5,i,j,k) - pn(2,i,j,k)*R(1) &
                                      - pn(3,i,j,k)*R(2) &
                                      - pn(4,i,j,k)*R(3)

        end do
      end do
    end do


  end subroutine update_TF


  !  Add sources to fluids
  subroutine add_TF_sources(upp,upn,primp,primn,dt)
    use globals, only: Temp, Tempn
    implicit none
    real, intent(inout)  ::   upp(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(inout)  ::   upn(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(in)     :: primp(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(in)     :: primn(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(in)     :: dt
    real                 :: s(neq)
    integer :: i, j, k

    do k=1,nz
      do j=1,ny
        do i=1,nx

          !  Gets the source terms
          call get_TF_sources(primp(:,i,j,k), primn(:,i,j,k),temp(i,j,k), tempn(i,j,k),s )

          !  Add the sources
          upp(:,i,j,k)=upp(:,i,j,k) + dt* s(:)
          upn(:,i,j,k)=upn(:,i,j,k) - dt* s(:)

        end do
      end do
    end do

    return

  end subroutine add_TF_sources

end module two_fluid
