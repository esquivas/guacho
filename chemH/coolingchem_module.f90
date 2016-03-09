! ==============================================!
!!!!!============ User Module =============!!!!!!
! ==============================================!
Module coolingchem_module
Contains
  subroutine initconds
    Use network
    implicit none
    Integer :: i
    real :: denn
!
denn=20000
y(iHII)=0.d0
y(iH2)=0.5d0*denn
y(iHI)=5.d-5*denn
y(ie)=y(iHII)
!!$y(CH2)=1.d0
!!$y(CO)=1.62e-4*1e7
!!$y(CO2)=1.d0
!!$y(ie)=1.d0
!!$y(H)=1.d0
!!$y(H2)=0.9998*1e7
!!$y(H2O)=1.d0
!!$y(HCO)=1.d0
!!$y(Hp)=1.d0
!!$y(O)=1.d0
!!$y(O2)=1.d0
!!$y(OH)=2.22e-6*1.e7
!
!-----------------------------------------------!
!
!!$    y(iM )=0.
!!$    Do i=1,nspec
!!$       y(iM )=y(i)+y(iM)
!!$    EndDo
    Return
  End Subroutine initconds
  !========================================
  
  Subroutine output(T,time)
    Use network
    Implicit None
    Real (kind=8), Intent(in) :: T,time
    Integer :: i
    !   
!    do i = 1, nsteps
       write(10,*) T,(y(i),i=1,n_spec+1)!,(y0(i),i=1,nelem)
       !
!    End do
!    close(10)  
    !
    Return
  End Subroutine output
  !
  Subroutine openfiles
    Implicit None
    Open(unit=10,file='yden.tab')
  End Subroutine openfiles
End Module coolingchem_module
