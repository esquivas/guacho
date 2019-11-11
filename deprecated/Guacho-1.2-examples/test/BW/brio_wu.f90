module BrioWu
  use parameters
  implicit none
  real :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, pL, pR
  real :: BxL, BxR, ByL, ByR, BzL, BzR
contains

subroutine init_bw()
  implicit none
  
  rhoL = 1.   ; rhoR =  0.125
  vxL  = 0.   ; vxR  =  0.
  vyL  = 0.   ; vyR  =  0.
  vzL  = 0.   ; vzR  =  0.
  pL   = 1.   ; pR   =  0.1
  BxL  = 0.75 ; BxR  =  0.75
  ByL  = 1.   ; ByR  = -1.
  BzL  = 0.   ; BzR  =  0.

  !rhoL = 1.08             ; rhoR =  1.0
  !pL   = 0.95             ; pR   =  1.0
  !vxL  = 1.2              ; vxR  =  0.
  !vyL  = 0.01             ; vyR  =  0.
  !vzL  = 0.5              ; vzR  =  0.
  !BxL  = 2. /sqrt(4.*pi)  ; BxR  = 2. /sqrt(4.*pi)
  !ByL  = 3.6/sqrt(4.*pi)  ; ByR  = 4. /sqrt(4.*pi)
  !BzL  = 2. /sqrt(4.*pi)  ; BzR  = 2. /sqrt(4.*pi)


end subroutine init_bw

!--------------------------------------------------------------------
! Initial conditions for the BRIO-WU shock tube
subroutine impose_bw(u)

use globals, only : coords, dx
use parameters, only : nxtot
implicit none
real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
real :: x
integer ::  i,j,k

do i=nxmin,nxmax
   do j=nymin,nymax
      do k=nzmin,nzmax

         ! Position measured from the centre of the grid
         x=(real(i+coords(0)*nx-nxtot/2)+0.5)*dx

         ! IF LEFFT STATE
         if( x <= 0.) then

            !   total density and momena
            u(1,i,j,k) = rhoL
            u(2,i,j,k) = rhoL*vxL
            u(3,i,j,k) = rhoL*vyL
            u(4,i,j,k) = rhoL*vzL
            !   total energy
            u(5,i,j,k)=0.5*rhoL*(vxL**2+vyL**2+vzL**2)+cv*pL  &
                           +0.5*(BxL**2+ByL**2+BzL**2)
            u(6,i,j,k) =  BxL
            u(7,i,j,k) =  ByL
            u(8,i,j,k) =  BzL

         else !RIGHT STATE
            !   total density and momenta
            u(1,i,j,k) = rhoR
            u(2,i,j,k) = rhoR*vxR
            u(3,i,j,k) = rhoR*vyR
            u(4,i,j,k) = rhoR*vzR
            !   energy
            u(5,i,j,k)=0.5*rhoR*(vxR**2+vyR**2+vzR**2)+cv*pR   &
                           +0.5*(BxR**2+ByR**2+BzR**2)
            u(6,i,j,k) =  BxR
            u(7,i,j,k) =  ByR
            u(8,i,j,k) =  BzR

         end if
            
      end do
   end do
end do

end subroutine impose_bw
!--------------------------------------------------------------------

end module  BrioWu
