! ======================================================================
!    Non equilibrium cooling routine 
!   And the primitives are updated in this routine as well
!======================================================================
subroutine coolingne(dt)
#ifdef COOLINGNE
  use parameters
  use globals
  use atomic
  use NEcooling
  implicit none
  real, dimension(nbmin:nbmax,0:nx,0:ny) :: n_tot 
  real,    intent(in)  :: dt
  real                 :: T ,Eth0
  real, parameter :: Tmin=100.
  real (kind=8)        :: L_rad, Ce
  integer :: i, j,nb,k
  !
  do nb=nbmin,nbmax
     if ( leaf(nb) ) then
        n_tot(nb,0,0)=0.0
        do i=1,nx	
           do j=1,ny
              !	
              !   get the primitives (and T)
              call uprim(primit(nb,:,i,j),u(nb,:,i,j),T)	
              !	
!#ifdef PASSIVES	
!              primit(nb,5:neq,i,j) = primit(nb,5:neq,i,j)/primit(nb,1,i,j)
!#endif	
              if(T.gt.Tmin) then	
                 !
                 !
                 !                 print*, 'primit=', primit(nb,5,i,j), primit(nb,6,i,j), primit(nb,7,i,j), primit(nb,1,i,j),  'T_i=', T

!                 print*,  neq, u(nb,5:neq,i,j)
                 
                 !   solve rate equations
                 call atomstep(dt,u(nb,5:neq,i,j),u(nb,1,i,j),T)

                 !  calculate cooling
                 call cooling(T,U(nb,5:neq,i,j),L_rad)

!                 do k=5,neq
!                    n_tot(nb,i,j)= n_tot(nb,0,0)+U(nb,k,i,j)
!                 end do
!                 print *, U(nb,1,i,j), n_tot(nb,i,j)
!                 print *, U(nb,1,i,j), U(nb,5:neq,i,j)
                 
                 !  update the passive scalars in the primitives array
                 primit(nb,5:neq,i,j)= u(nb,5:neq,i,j)

                 !print *,  'T_f=', T
                 Eth0=cv*primit(nb,4,i,j)                 
                 Ce=L_rad/(Eth0*Psc)  ! cgs
                 !
                 !  apply cooling to primitive and conserved variables
                 primit(nb,4,i,j)=primit(nb,4,i,j)*exp(-ce*dt)

                 u(nb,4,i,j)=u(nb,4,i,j)-Eth0+cv*primit(nb,4,i,j)
                 !
              end if
!!              !
           end do
        end do
     endif
  end do

  !
  !--------------------------------------------------------------------
  !
#endif
end subroutine coolingne
!======================================================================
