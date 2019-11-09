!=================================

 module network

   use parameters, only : n_spec
   implicit none

   ! number of equilibrium equations
   integer, parameter :: n_equi = 3

   ! number of non-equilibrium equations
   integer, parameter :: n_nequ = n_spec - n_equi

   ! number of total elements
   integer, parameter :: n_elem = 2

   ! indexes of the different species
   integer, parameter :: iHI    = 1
   integer, parameter :: iHII   = 2
   integer, parameter :: iHeI   = 3
   integer, parameter :: iHeII  = 4
   integer, parameter :: iHeIII = 5
   integer, parameter :: ie     = 6

   ! indexes of the equilibrium species
   integer, parameter :: iH  = 1
   integer, parameter :: iHe = 2

   ! number of reaction rates
   integer, parameter :: n_reac = 7

   ! indexes of the different rates
   integer, parameter :: ichi    = 1
   integer, parameter :: ichei   = 2
   integer, parameter :: icheii  = 3
   integer, parameter :: iahii   = 4
   integer, parameter :: iaheii  = 5
   integer, parameter :: iaheiii = 6
   integer, parameter :: iphiH   = 7

   ! names of the equation to be solved
   character(len=10)  :: name_eqs(n_spec)

   ! Size of the table
   real (kind=8) :: rate_table(141,n_reac)

   ! number of different densities
   integer, parameter     :: ndens = 21

   ! This is the table when we store the cooling rates
   real (kind=8) ::  lambda_table(200,21,n_spec)

   !=======

 contains

   !=======

   subroutine derv(y,rate,dydt,y0)

     implicit none
     real (kind=8), intent(in)  ::   y0(n_elem)
     real (kind=8), intent(in)  ::    y(n_spec)
     real (kind=8), intent(out) :: dydt(n_spec)
     real (kind=8), intent(in)  :: rate(n_reac)

     dydt(iHI) = - rate(ichi )*y(iHI )*y(ie)                                   &
                 + rate(iahii)*y(iHII)*y(ie)                                   &
                 - rate(iphiH)*y(iHI )

     dydt(iHeI) = - rate(ichei)*y(ie)*y(iHeI)                                  &
                  + rate(iaheii)*y(ie)*y(iHeII)

     dydt(iHeII) =  rate(ichei )*y(ie)*y(iHeI )                                &
                  - rate(icheii )*y(ie)*y(iHeII )                              &
                  - rate(iaheii)*y(ie)*y(iHeII)                                &
                  + rate(iaheiii)*y(ie)*y(iHeIII)

     !  "conservation" equations
     dydt(iHII  ) = - y0(iH) + y(iHI) + y(iHII)

     dydt(iHeIII) = - y0(iHe) + y(iHeI) + y(iHeII) + y(iHeIII)

     dydt(ie)     = y(ie) - y(iHII) - y(iHeII) - 2.*y(iHeIII)

   end subroutine derv

   !=======

   subroutine get_jacobian(y,jacobian,rate)

     implicit none
     real (kind=8), intent(in)  :: y(n_spec)
     real (kind=8), intent(out) ::jacobian(n_spec,n_spec)
     real (kind=8), intent(in)  :: rate(n_reac)

    jacobian(iHI,iHI)       = - rate(ichi)*y(ie) - rate(iphiH)
    jacobian(iHI,iHII)      = + rate(iahii)*y(ie)
    jacobian(iHI,iHeI)      =   0.
    jacobian(iHI,iHeII)     =   0.
    jacobian(iHI,iHeIII)    =   0.
    jacobian(iHI,ie)        = - rate(ichi)*y(iHI) + rate(iahii)*y(iHII)

    jacobian(iHII,iHI)      =   1.
    jacobian(iHII,iHII)     =   1.
    jacobian(iHII,iHeI)     =   0.
    jacobian(iHII,iHeII)    =   0.
    jacobian(iHII,iHeIII)   =   0.
    jacobian(iHII,ie)       =   0.

    jacobian(iHeI,iHI)      =   0.
    jacobian(iHeI,iHII)     =   0.
    jacobian(iHeI,iHeI)     = - rate(ichei)*y(ie)
    jacobian(iHeI,iHeII)    = + rate(iaheii)*y(ie)
    jacobian(iHeI,iHeIII)   =   0.
    jacobian(iHeI,ie)       = - rate(ichei)*y(iHeI) + rate(iaheii)*y(iHeII)

    jacobian(iHeII,iHI)     =  0.
    jacobian(iHeII,iHII)    =  0.
    jacobian(iHeII,iHeI)    = + rate(ichei)*y(ie)
    jacobian(iHeII,iHeII)   = - rate(icheii)*y(ie) - rate(iaheii)*y(ie)
    jacobian(iHeII,iHeIII)  = + rate(iaheiii)*y(ie)
    jacobian(iHeII,ie)      = + rate(ichei)*y(iHeI) - rate(icheii)*y(iHeII)    &
                          - rate(iaheii)*y(iHeII) + rate(iaheiii)*y(iHeIII)

    jacobian(iHeIII,iHI)    =   0.
    jacobian(iHeIII,iHII)   =   0.
    jacobian(iHeIII,iHeI)   =   1.
    jacobian(iHeIII,iHeII)  =   1.
    jacobian(iHeIII,iHeIII) =   1.
    jacobian(iHeIII,ie)     =   0.

    jacobian(ie,iHI)        =   0.
    jacobian(ie,iHII)       = - 1.
    jacobian(ie,iHeI)       =   0.
    jacobian(ie,iHeII)      = - 1.
    jacobian(ie,iHeIII)     = - 2.
    jacobian(ie,ie)         =   1.

   end subroutine get_jacobian

   !=======================================================================

   subroutine get_reaction_rates(rate,T,phiH)
     implicit none
     real (kind=8), intent(in) :: T, phiH
     real (kind=8),intent(out) :: rate(n_reac)

     ! indexes of the different rates
     rate(ichi   ) = 5.830e-11*sqrt(T)*exp(-157828.0/T)
     rate(ichei  ) = 2.707e-11*sqrt(T)*exp(-285400.0/T)
     rate(icheii ) = 5.711e-12*sqrt(T)*exp(-631000.0/T)
     rate(iahii  ) = 2.55e-13*(1.0e4/T)**0.79
     rate(iaheii ) = 4.30e-13*(1.0e4/T)**0.672                                 &
                   + 0.0019*T**(-1.5)*exp(-4.7e5/T)*(1.0+0.3*exp(-94000.0/T) )
     rate(iaheiii) = 2.21e-9*T**0.79
     rate(iphiH  ) = phiH

   end subroutine get_reaction_rates

   !=======================================================================

   subroutine nr_init(y,y0)
     implicit none
     real, intent(out) :: y(n_spec)
     real, intent(in ) :: y0(n_elem)

     y(iHI    ) = 0.5     * y0(iH )
     y(iHII   ) = 0.5     * y0(iH )
     y(iHeI   ) = 1.0/3.0 * y0(iHe)
     y(iHeII  ) = 1.0/3.0 * y0(iHe)
     y(iHeIII ) = 1.0/3.0 * y0(iHe)
     y(ie     ) = y(iHII) + y(iHeII) + 2.*y(ieHII)

     return
   end subroutine nr_init

   !=======================================================================

   logical function check_no_conservation(y,y0_in)
     implicit none
     real, intent(in)  :: y(n_spec)
     real, intent(in ) :: y0_in  (n_elem)
     real              :: y0_calc(n_elem)
     integer           :: i

     check_no_conservation = .false.

     y0_calc(iH )= y(iHI ) + y(iHII )
     y0_calc(iHe)= y(iHeI) + y(iHeII) + y(iHeIII)

     do i = 1, n_elem
       if ( y0_calc(i) > 1.001*y0_in(i) ) check_no_conservation = .true.
     end do

   end function check_no_conservation

   !=======================================================================


 end module network

!=================================
