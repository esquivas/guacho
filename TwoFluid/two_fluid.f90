module two_fluid

  use parameters,  only : neq,  nxmin, nymin, nzmin, &
                         nxmax, nymax, nzmax, nx, ny, nz

  implicit none

  real :: alpha = 1e-2


contains

!> @brief Twofluid interaction term
!> @brief from Draine 1980

subroutine get_TF_sources(pp, pn,s)
 
  implicit none
  real, intent(in)  :: pp(neq), pn(neq)
  real, intent(out) :: s(neq)
  real :: driftx, drifty, driftz
  !  update the source terms
  
  driftx = pn(2)-pp(2)
  drifty = pn(3)-pp(3)
  driftz = pn(4)-pp(4)
  
  ! reset sources
  s(:) = 0.

  !  momenta
  s(2) = alpha * pp(1)*pn(1)*driftx
  s(3) = alpha * pp(1)*pn(1)*drifty
  s(4) = alpha * pp(1)*pn(1)*driftz
  s(5) = alpha * pp(1)*pn(2)* (driftx*pp(2)+drifty*pp(3)+driftz*pp(4) )

end subroutine get_TF_sources

! obtains R
subroutine get_TF_R(pp, pn,dt, R)
 
  implicit none
  real, intent(in)  :: pp(neq), pn(neq), dt
  real, intent(out) :: R(3)
  
  R(1) = alpha * dt * pp(1) * pn(1) *( pn(2)-pp(2) )/(1+alpha*dt*(pp(1)+pn(1)))
  R(2) = alpha * dt * pp(1) * pn(1) *( pn(3)-pp(3) )/(1+alpha*dt*(pp(1)+pn(1)))
  R(3) = alpha * dt * pp(1) * pn(1) *( pn(4)-pp(4) )/(1+alpha*dt*(pp(1)+pn(1)))

end subroutine get_TF_R


subroutine update_TF(u,un,up,upn,pp,pn,dt)
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
          call get_TF_R(pp(:,i,j,k), pn(:,i,j,k),dt,R )
          
          !  update momenta
          up(2,i,j,k) = up(2,i,j,k) + R(1)
          up(3,i,j,k) = up(3,i,j,k) + R(2)
          up(4,i,j,k) = up(4,i,j,k) + R(3)

          upn(2,i,j,k) = upn(2,i,j,k) - R(1)
          upn(3,i,j,k) = upn(3,i,j,k) - R(2)
          upn(4,i,j,k) = upn(4,i,j,k) - R(3)
          
          !  update enery
          !up(5,i,j,k) = u(5,i,j,k)   + pp(2,i,j,k)*R(1) &
          !                           + pp(3,i,j,k)*R(2) &
          !                           + pp(4,i,j,k)*R(3)

          !upn(5,i,j,k) = un(5,i,j,k) - pn(2,i,j,k)*R(1) &
          !                           - pn(3,i,j,k)*R(2) &
          !                           - pn(4,i,j,k)*R(3)

        end do
     end do
  end do
  

end subroutine update_TF


!  Add sources to fluids
subroutine add_TF_sources(upp,upn,primp,primn,dt)
  implicit none
  real, intent(inout)  ::   upp(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent(inout)  ::   upn(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent(in)     :: primp(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent(in)     :: primn(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent(in)    :: dt
  real                 :: s(neq)
  integer :: i, j, k

  do k=1,nz
     do j=1,ny
        do i=1,nx

          !  Gets the source terms
          call get_TF_sources(primp(:,i,j,k), primn(:,i,j,k),s )
          
          !  Add the sources
          upp(:,i,j,k)=upp(:,i,j,k) + dt* s(:)
          upn(:,i,j,k)=upn(:,i,j,k) - dt* s(:)

        end do
     end do
  end do
  
  return

end subroutine add_TF_sources

end module two_fluid