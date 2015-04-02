!=======================================================================
!   Cooling routine using Benjamin, Bensonm, & Cox (2001) ApJL
!   non-equilibrium ionization routines.
!=======================================================================

module cooling_bbc

subroutine coolingbbc(dt)
#ifdef COOLINGBBC
  use parameters
  use globals
  implicit none
  real, dimension(neq) ::  uu, prim
  real, dimension(npas) :: y, yi
  real    ::  T, dt, Eth, Emin, Efloor, dE
  real, parameter :: Tmin=10.**(3.5) 
  integer :: i,j,k
  !
    do i=0,nx+1
     do j=0,ny+1
        do k=0,nz+1
           !
           !  only outside the obstacle 
           !if(obs(i,j).eq.0) then
           !
           uu(:)=u(:,i,j,k)
           call u2prim(uu,prim,T)
           !
           if(T.ge.1.1e4) then

           yi(:) = prim(6:12)/prim(1)*amol
           !yi(6) = 0.183318
           !
           !  update ionization vector and passive scalars
           call bbcev(dt,prim(1)/amol,T, yi,y)           
           uu(6:neqpas)=y(:)*prim(1)/amol
           !
           Eth= uu(5)-0.5*prim(1)*(prim(2)**2.+prim(3)**2.+prim(4)**2)
           !
           Efloor=0.5*Eth
           !
           !  compute delta E w/cooling (in units of ergs s-1 cm-3)
           dE= bbccool(1.,prim(1)/amol,T,y)*dt
           !
           !  change to normalized units
           dE= dE/rhosc/vsc/vsc
           dE= min(dE, Efloor)
           !
           !  update thermal energy
           Eth=Eth*exp(-dE/Eth)
           !
           !  update Etot
           uu(5)=Eth+0.5*prim(1)*(prim(2)**2.+prim(3)**2.+prim(4)**2.)
           !
           endif

           u(:,i,j,k)=uu(:)
           !
           !end if
           !
        end do
     end do
  end do
  !
  return
#endif
end subroutine coolingbbc

!=======================================================================

end module

!=======================================================================