!=======================================================================
!> @file hydro_core.f90
!> @brief Hydrodynamical and Magnetohidrodynamocal bacic module
!> @author Alejandro Esquivel, V. Sieyra, M. Shneiter
!> @date 4/May/2016
! Copyright (c) 2020 Guacho Co-Op
!
! This file is part of Guacho-3D.
!
! Guacho-3D is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief Basic hydro (and MHD) subroutines utilities
!> @details This module contains subroutines and utilities that are the
!> core of the hydro (and MHD) that are common to most implementations
!> and will be used for the different specific solvers

module hydro_core

  implicit none

contains

  !=======================================================================
  !> @brief Computes the primitive variables and temperature from conserved
  !!  variables on a single cell
  !> @details Computes the primitive variables and temperature from conserved
  !!  variables on a single cell
  !> @param real [in] uu(neq) : conserved variables in one cell
  !> @param real [out] prim(neq) : primitives in one cell
  !> @param real [out] T : Temperature [K]
  
  subroutine u2prim(uu, prim, T)

    use parameters, only : neq, neqdyn, Tempsc, vsc2, cv, passives,   &
                           pmhd, mhd, eq_of_state, n1_chem, n_spec,   &
                           riemann_solver, gamma
    use constants
    implicit none
    real,    intent(in),  dimension(neq)  :: uu
    real,    intent(out), dimension(neq)  :: prim
    real,    intent(out)                  :: T
    real :: r

    real :: a1, a2, a3
    real :: b1, b2, b3, b4
    real :: B, C, M, RR, SS, TT
    real :: x1
    real :: gamma_rel
    real :: v
    real :: denom

#ifdef PASSIVES
    integer :: i
    real :: dentot
#endif
    real, parameter :: dens_floor = 1e-8
    real, parameter :: Temp_floor = 10.0

    if (riemann_solver == SOLVER_RHLL .or. &
        riemann_solver == SOLVER_RHLLC) then
        
      M = sqrt(uu(2)**2 + uu(3)**2 + uu(4)**2)

      if (M > 0.0) then

        denom = (M**2 + uu(1)**2)*(gamma - 1.0)**2

        b1 = -(2*gamma*(gamma-1.0)*M*uu(5))/denom
        b2 = ((gamma*uu(5))**2 + 2*(gamma-1.0)*M**2 - ((gamma -1.0)*uu(1))**2)/denom
        b3 = -2*gamma*M*uu(5)/denom
        b4 = M**2/denom

        a1 = -b2
        a2 = b1*b3 - 4*b4
        a3 = 4*b2*b4 - b3**2 - b4*b1**2

        RR = (9*a1*a2 - 27*a3 - 2*a1**3)/54.0
        SS = (3*a2 - a1**2)/9.0
        TT = RR**2 + SS**3

        x1 = (RR + sqrt(TT))**(1.0/3.0) + (RR - sqrt(TT))**(1.0/3.0)-a1/3.0

        B = 0.5*(b1 + sqrt(b1**2 - 4*b2 + 4*x1))
        C = 0.5*(x1 - sqrt(x1**2 - 4*b4))

        v = 0.5*(-B + sqrt(B**2 - 4*C))
        gamma_rel = 1.0/sqrt(1 - v**2)
        prim(1) = uu(1)/gamma_rel
        prim(2) = uu(2)*v/M
        prim(3) = uu(3)*v/M
        prim(4) = uu(4)*v/M

      else

        v = 0.0
        gamma_rel = 1.0
        prim(1) = uu(1)/gamma_rel
        prim(2) = 0.0
        prim(3) = 0.0
        prim(4) = 0.0

      end if

      prim(5)=(gamma-1)*(uu(5) - uu(2)*prim(2) - uu(3)*prim(3) - uu(4)*prim(4) - prim(1))
  
    else
  
      r=max(uu(1),dens_floor)
      prim(1)=r
  
      prim(2)=uu(2)/r
      prim(3)=uu(3)/r
      prim(4)=uu(4)/r


      if (mhd) then
#ifdef  BFIELD
        prim(5)=( uu(5)-0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2)                 &
                        -0.5*  (  uu(6)**2+  uu(7)**2  +uu(8)**2) ) /cv
#endif
      else
        prim(5)=( uu(5)-0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2) ) /cv
      end if

      prim(5)=max(prim(5),1e-16)    

    end if

    if (mhd .or. pmhd) then
#ifdef  BFIELD
      prim(6:8) = uu(6:8)
#endif
    end if

    if (passives) then
#ifdef PASSIVES
      prim(neqdyn+1:neq) = uu(neqdyn+1:neq)
#endif
    end if

    !   Temperature calculation
    if (eq_of_state == EOS_ADIABATIC) then
      T=max( (prim(5)/r)*Tempsc, Temp_floor)
      prim(5)= r*T/Tempsc
    end if

    if (eq_of_state == EOS_SINGLE_SPECIE) then
      ! assumes it is fully ionized
      T=max(Temp_floor,(prim(5)/r)*Tempsc)
      prim(5)=r*T/Tempsc
    end if

#ifdef PASSIVES

    if (eq_of_state == EOS_H_RATE) then
      dentot=(2.*r-prim(neqdyn+1))
      dentot=max(dentot, dens_floor)
      T=max(Temp_floor,(prim(5)/dentot)*Tempsc)
      prim(5)=dentot*T/Tempsc
    end if

    if (eq_of_state == EOS_CHEM) then
      !  Assumes that rho scaling is mu*mh
      dentot = 0.
      do i = n1_chem, n1_chem+n_spec-1
        dentot = prim(i) + dentot
      end do
      dentot = max(dentot, dens_floor)
      T=max(Temp_floor,(prim(5)/dentot)*Tempsc )

      !T=(prim(5)/r)*Tempsc
      prim(5) = dentot * T /Tempsc
    end if

#endif

  end subroutine u2prim

  !=======================================================================
  !> @brief Computes the primitive variables and temperature from conserved
  !!  variables on a single cell
  !> @details Computes the primitive variables and temperature from conserved
  !!  variables on a single cell for the split method in all the variables
  !> @param real [in] uu(neq) : conserved variables in one cell (fluctuations)
  !> @param real [out] prim(neq) : primitives in one cell (fluctuations)
  !> @param real [in] prim0(neq) : primitives in one cell (mean value,
  !> initial conds)
  !> @param real [out] T : Temperature [K]
  subroutine u2primSplitAll(uu, prim, prim0, T)

    use parameters, only : neq, neqdyn, Tempsc, vsc2, cv, passives,            &
                           pmhd, mhd, eq_of_state, n_spec
    use constants
    implicit none
    real,    intent(in),  dimension(neq)  :: uu, prim0
    real,    intent(out), dimension(neq)  :: prim
    real,    intent(out)                  :: T
    real :: r
#ifdef PASSIVES
    integer :: i
    real :: dentot
#endif

    prim(1)=uu(1)
    !  r=max(uu(1),1e-15)
    r = prim(1)+prim0(1)

    prim(2)=uu(2)/r
    prim(3)=uu(3)/r
    prim(4)=uu(4)/r

#ifdef  BFIELD
    if (mhd .or. pmhd) then
      prim(6:8) = uu(6:8)
    end if
#endif

    if (mhd) then
#ifdef BFIELD
      prim(5)=( uu(5) - 0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2)               &
                      - 0.5*  (uu(6)**2+uu(7)**2+uu(8)**2)                     &
                      - prim0(6)*uu(6)-prim0(7)*uu(7)-prim0(8)*uu(8) ) /cv
#endif
    else
      prim(5)=( uu(5)-0.5*r*(prim(2)**2+prim(3)**2+prim(4)**2) ) /cv
    end if

  !  prim(5)=max(prim(5),1e-16)

#ifdef PASSIVES
    if (passives) then
       prim(neqdyn+1:neq) = uu(neqdyn+1:neq)
    end if
#endif

    !   Temperature calculation

    if (eq_of_state == EOS_ADIABATIC) then
      T=((prim(5)+prim0(5))/r)*Tempsc
    end if

    if (eq_of_state == EOS_SINGLE_SPECIE) then
      ! assumes it is fully ionized
      r=max(r,1e-15)
      T=max(1.,((prim(5)+prim0(5))/r)*Tempsc)
      !prim(5)=r*T/Tempsc ! REPENSAR
    end if

#ifdef PASSIVES

    if (eq_of_state == EOS_H_RATE) then
      dentot=(2.*r-prim(neqdyn+1)-prim0(neqdyn+1))
      dentot=max(dentot,1e-15)
      T=max(1.,((prim(5)+prim0(5))/dentot)*Tempsc)
      !prim(5)=dentot*T/Tempsc !REVISAR
    end if

    if (eq_of_state == EOS_CHEM) then
      !  Assumes that rho scaling is mu*mh
      dentot = 0.
      do i = neqdyn+1, neqdyn+n_spec
        dentot = prim(i) + prim0(i) + dentot
      end do
      dentot = max(dentot, 1e-15)
      T=max(1.,((prim(5)+prim0(5))/dentot)*Tempsc )

      !T=(prim(5)/r)*Tempsc
      !prim(5) = dentot * T /Tempsc
    end if

#endif

  end subroutine u2primSplitAll

  !=======================================================================
  !> @brief Updated the primitives, using the conserved variables in the
  !> entire domain
  !> @details Updated the primitives, using the conserved variables in the
  !> entire domain
  !> @param real [in] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
  !> conserved variables
  !> @param real [out] prim(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
  !> primitive variables
  !> @param logical [in] only_ghost : if set to true then updates the primitives
  !> only on the ghost cells, it defaults to false (the entire domain)
  subroutine calcprim(u,primit, only_ghost)

    use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                           nx, ny, nz, riemann_solver
    use constants
    use globals, only : Temp, primit0
    implicit none
    real,intent(in)  :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real,intent(out) :: primit(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    logical, optional, intent(in) :: only_ghost
    integer :: i,j,k

    if (present(only_ghost) .and. only_ghost) then
      !-----------------------
      !   k = 0, and, nz
      do j=0,ny+1
        do i=0,nx+1
          if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or. &
              riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
            call u2primSplitAll(u(:,i,j,0   ),    primit(:,i,j,0 ),            &
                                primit0(:,i,j,0 ),Temp(i,j,0   ) )
            call u2primSplitAll(u(:,i,j,nz+1),       primit(:,i,j,nz+1),       &
                                primit0(:,i,j,nz+1), Temp(i,j,nz+1) )
          else
            call u2prim(u(:,i,j,0   ), primit(:,i,j,0   ), Temp(i,j,0   ) )
            call u2prim(u(:,i,j,nz+1), primit(:,i,j,nz+1), Temp(i,j,nz+1) )
          end if
        end do
      end do

      !   j = 0, and, ny
      do k=0,nz+1
        do i=0,nx+1
          if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or. &
              riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
            call u2primSplitAll(u(:,i,0   ,k),       primit(:,i,0   ,k),       &
                                primit0(:,i,0   ,k), Temp(i,0   ,k) )
            call u2primSplitAll(u(:,i,ny+1,k),          primit(:,i,ny+1,k),    &
                                primit0(:,i,ny+1   ,k), Temp(i,ny+1,k) )
          else
            call u2prim(u(:,i,0   ,k),primit(:,i,0   ,k),Temp(i,0   ,k) )
            call u2prim(u(:,i,ny+1,k),primit(:,i,ny+1,k),Temp(i,ny+1,k) )
          end if
        end do
      end do

      !   i = 0, and, nx
      do k=0,nz+1
        do j=0,ny+1
          if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or. &
              riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
            call u2primSplitAll(u(:,0   ,j,k),      primit(:,0   ,j,k),        &
                               primit0(:,0   ,j,k), Temp(0,j,k) )
            call u2primSplitAll(u(:,nx+1,j,k),      primit(:,nx+1,j,k),        &
                               primit0(:,nx+1,j,k), Temp(nx+1,j,k) )
          else
            call u2prim(u(:,0   ,j,k), primit(:,0   ,j,k), Temp(0   ,j,k) )
            call u2prim(u(:,nx+1,j,k), primit(:,nx+1,j,k), Temp(nx+1,j,k) )
          end if
        end do
      end do

    else
      !   entire domain
      do k=nzmin,nzmax
        do j=nymin,nymax
          do i=nxmin,nxmax
            if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or. &
                riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
              call u2primSplitAll(u(:,i,j,k),primit(:,i,j,k), primit0(:,i,j,k),&
                                 Temp(i,j,k) )
            else
              call u2prim(u(:,i,j,k), primit(:,i,j,k), Temp(i,j,k) )
            end if

          end do
        end do
      end do

    end if

  end subroutine calcprim

  !=======================================================================
  !> @brief Computes the conserved conserved variables from the
  !> primitives in a single cell
  !> @details Computes the conserved conserved variables from the
  !> primitives in a single cell
  !> @param real [in] prim(neq) : primitives in one cell
  !> @param real [out] uu(neq) : conserved varibles in one cell

  ! AGREGAR LA CONVERSION RELATIVISTA ENTRE PRIMITVAS A CONSERVADAS

  subroutine prim2u(prim,uu, prim0)

    use parameters
    implicit none
    real, dimension(neq), intent(in)  :: prim
    real, dimension(neq), intent(in), optional :: prim0
    real, dimension(neq), intent(out) :: uu

    if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or.   &
        riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
      if (present(prim0)) then
        uu(1) = prim(1)
        uu(2) = prim(2)*(prim(1)+prim0(1))
        uu(3) = prim(3)*(prim(1)+prim0(1))
        uu(4) = prim(4)*(prim(1)+prim0(1))
      end if
    else if (riemann_solver == SOLVER_RHLL .or.   &
            riemann_solver == SOLVER_RHLLC) then

            ! CONVERSION

    end if

    ! energy for hydro and passive mhd
    uu(5) = 0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)+cv*prim(5)

#ifdef BFIELD
    if (mhd) then
      if (present(prim0)) then
        !   kinetic+thermal+magnetic energies
        uu(5) = uu(5) + 0.5*(prim(6)**2+prim(7)**2+prim(8)**2)                 &
              + 0.5*prim0(1)*(prim(2)**2+prim(3)**2+prim(4)**2)                &
              + prim0(6)*prim(6)+prim0(7)*prim(7)+prim0(8)*prim(8)
      else
        uu(5) = uu(5) + 0.5*(prim(6)**2+prim(7)**2+prim(8)**2)
      end if
    end if
#endif

#ifdef BFIELD
    if (mhd .or. pmhd) then
      uu(6:8)=prim(6:8)
    end if
#endif

#ifdef PASSIVES
    if (passives) then
      uu(neqdyn+1:neq) = prim(neqdyn+1:neq)
    end if
#endif

  end subroutine prim2u

  !=======================================================================
  !> @brief Computes the Euler Fluxes in one cell
  !> @details Computes the Euler Fluxes in one cell, using the primitices
  !> @n It returns the flux in the x direction (i.e. F), the y and z fluxes
  !> can be obtained swaping the respective entries (see swapy and swapz
  !>	subroutines)
  !> @param real [in] prim(neq) : primitives in one cell
  !> @param real [out] ff(neq) : Euler Fluxes (x direction)
  subroutine prim2f(prim,ff,prim0)

    use parameters, only : neq, neqdyn, cv, pmhd, mhd, passives, riemann_solver
    use constants, only : SOLVER_HLLE_SPLIT_ALL, SOLVER_HLLD_SPLIT_ALL,        &
                          SOLVER_RHLL, SOLVER_RHLLC
    implicit none
    real,    dimension(neq), intent(in)  :: prim
    real, dimension(neq), intent(in), optional :: prim0
    real,    dimension(neq), intent(out) :: ff
    real :: etot

    if (mhd) then
      if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or. &
          riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
        if (present(prim0)) then
#ifdef BFIELD
          etot = 0.5*( (prim(1)+prim0(1))*(prim(2)**2+prim(3)**2+prim(4)**2)   &
          + prim(6)**2+prim(7)**2+prim(8)**2  )                                &
          + cv*prim(5)                                                         &
          + prim0(6)*prim(6)+prim0(7)*prim(7)+prim0(8)*prim(8)

          ff(1) = (prim(1)+prim0(1))*prim(2)
          ff(2) = (prim(1)+prim0(1))*prim(2)*prim(2)+prim(5)+0.5*(prim(7)**2   &
                +prim(8)**2-prim(6)**2) &
                - prim0(6)*prim(6)+prim0(7)*prim(7)+prim0(8)*prim(8)
          ff(3) = (prim(1)+prim0(1))*prim(2)*prim(3)-prim(6)*prim(7)           &
                - prim0(7)*prim(6)-prim0(6)*prim(7)
          ff(4) = (prim(1)+prim0(1))*prim(2)*prim(4)-prim(6)*prim(8)           &
                - prim0(8)*prim(6)-prim0(6)*prim(8)

          ff(5) = prim(2)*(etot+prim(5)+0.5*((prim(6)+prim0(6))**2+(prim(7)    &
                +prim0(7))**2+(prim(8)+prim0(8))**2)                           &
                +cv*prim0(5)+prim0(5)+0.5*(prim0(6)**2                         &
                +prim0(7)**2+prim0(8)**2) )                                    &
                -(prim(6)+prim0(6))*(prim(2)*(prim(6)+prim0(6))                &
                +prim(3)*(prim(7)+prim0(7))+prim(4)*(prim(8)+prim0(8))) !REVISAR
#endif
        end if
        else
#ifdef BFIELD
          !  MHD (active)
          etot = 0.5*( prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)              &
               + prim(6)**2+prim(7)**2+prim(8)**2  ) + cv*prim(5)

          ff(1) = prim(1)*prim(2)
          ff(2) = prim(1)*prim(2)*prim(2)+prim(5)+0.5*(prim(7)**2              &
                +  prim(8)**2-prim(6)**2)
          ff(3) = prim(1)*prim(2)*prim(3)-prim(6)*prim(7)
          ff(4) = prim(1)*prim(2)*prim(4)-prim(6)*prim(8)
          ff(5) = prim(2)*(etot+prim(5)+0.5*(prim(6)**2                        &
                + prim(7)**2+prim(8)**2) ) &
                - prim(6)*(prim(2)*prim(6)+prim(3)*prim(7)+prim(4)*prim(8))
#endif
      end if
    else if (riemann_solver == SOLVER_RHLL .or. &
             riemann_solver == SOLVER_RHLLC) then

             ! CALCULAR LOS FLUJOS RELATIVISTAS

    else
      ! HD or PMHD
      etot= 0.5*prim(1)*(prim(2)**2+prim(3)**2+prim(4)**2)+cv*prim(5)

      ff(1) = prim(1)*prim(2)
      ff(2) = prim(1)*prim(2)*prim(2)+prim(5)
      ff(3) = prim(1)*prim(2)*prim(3)
      ff(4) = prim(1)*prim(2)*prim(4)
      ff(5) = prim(2)*(etot+prim(5))
    end if

#ifdef BFIELD
    if (mhd .or. pmhd) then
      if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or. &
          riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
        ff(6)=0.
        ff(7)=prim(2)*(prim0(7)+prim(7))-prim(3)*(prim0(6)+prim(6))
        ff(8)=prim(2)*(prim0(8)+prim(8))-prim(4)*(prim0(6)+prim(6))
      else
        ff(6)=0.0
        ff(7)=prim(2)*prim(7)-prim(6)*prim(3)
        ff(8)=prim(2)*prim(8)-prim(6)*prim(4)
      endif
    endif
#endif

#ifdef PASSIVES
    if (passives) then
      ff(neqdyn+1:neq) = prim(neqdyn+1:neq)*prim(2)
    end if
#endif

  end subroutine prim2f

  !=======================================================================
  !> @brief Swaps the x and y components in a cell.
  !> @details Swaps the x and y components in a cell.
  !> @param real [inout] var(neq) : variable to be swapped
  !> @param real [in] neq : number of equations in the code
  subroutine swapy(var,neq)

    use parameters, only : pmhd, mhd
    implicit none
    integer, intent(in) :: neq
    real, intent(inout), dimension(neq) :: var
    real :: aux

    aux=var(2)
    var(2)=var(3)
    var(3)=aux

    if (mhd .or. pmhd) then
#ifdef BFIELD
      aux=var(6)
      var(6)=var(7)
      var(7)=aux
#endif
    end if

  end subroutine swapy

  !=======================================================================
  !> @brief Swaps the x and z components in a cell.
  !> @details Swaps the x and z components in a cell.
  !> @param real [inout] var(neq) : variable to be swapped
  !> @param real [in] neq : number of equations in the code
  subroutine swapz(var,neq)

    use parameters, only : pmhd, mhd
    implicit none
    integer, intent(in) :: neq
    real, intent(inout), dimension(neq) :: var
    real :: aux

    aux=var(2)
    var(2)=var(4)
    var(4)=aux

    if (mhd .or. pmhd) then
#ifdef BFIELD
      aux=var(6)
      var(6)=var(8)
      var(8)=aux
#endif
    end if

  end subroutine swapz

  !=======================================================================
  !> @brief Computes the sound speed
  !> @details Computes the sound speed
  !> @param real [in] p : value of pressure
  !> @param real [in] d : value of density
  !> @param real [out] cs : sound speed
  subroutine csound(p,d,cs)

    use parameters, only : gamma
    implicit none
    real, intent(in) :: p, d
    real, intent(out) ::cs

    cs=sqrt(gamma*p/d)

  end subroutine csound

  !=======================================================================
  !> @brief Computes the fast magnetosonic speeds  in the 3 coordinates
  !> @details Computes the fast magnetosonic speeds  in the 3 coordinates
  !> @param real [in] p  : value of pressure
  !> @param real [in] d  : value of density
  !> @param real [in] Bx : value of the x component of the magnetic field
  !> @param real [in] By : value of the y component of the magnetic field
  !> @param real [in] Bz : value of the z component of the magnetic field
  !> @param real [out] csx : fast magnetosonic speed in x
  !> @param real [out] csy : fast magnetosonic speed in y
  !> @param real [out] csz : fast magnetosonic speed in z
  subroutine cfast(p,d,bx,by,bz,cfx,cfy,cfz)

    use parameters, only : gamma
    implicit none
    real, intent(in) :: p, d, bx, by, bz
    real, intent(out) ::cfx,cfy,cfz
    real :: b2

    b2=bx*bx+by*by+bz*bz
    cfx=sqrt(0.5*((gamma*p+b2)+sqrt((gamma*p+b2)**2-4.*gamma*p*bx*bx))/d)
    cfy=sqrt(0.5*((gamma*p+b2)+sqrt((gamma*p+b2)**2-4.*gamma*p*by*by))/d)
    cfz=sqrt(0.5*((gamma*p+b2)+sqrt((gamma*p+b2)**2-4.*gamma*p*bz*bz))/d)

  end subroutine cfast

#ifdef BFIELD
  !=======================================================================
  !> @brief Computes the fast magnetosonic speed in the x direction
  !> @details Computes the fast magnetosonic speed in the x direction
  !> @param real [in] prim(neq) : vector with the primitives in one cell
  subroutine cfastX(prim,cfX)

    use parameters, only : neq, gamma
    implicit none
    real, intent(in) :: prim(neq)
    real, intent(out) ::cfX
    real :: b2, cs2va2

    b2=prim(6)**2+prim(7)**2+prim(8)**2
    cs2va2 = (gamma*prim(5)+b2)/prim(1)   ! cs^2 + ca^2

    cfx=sqrt(0.5*(cs2va2+sqrt(cs2va2**2                                        &
             - 4.*gamma*prim(5)*prim(6)**2/prim(1)/prim(1) ) ) )

  end subroutine cfastX

#endif

  !=======================================================================
  !> @brief Otains the timestep allowed by the CFL condition in the entire
  !> @details Otains the timestep allowed by the CFL condition in the entire
  !> domain using the global primitives, and sets logical variable to dump
  !> output
  !> @param integer [in] current_iter : Current iteration, it starts with a small
  !> but increasing CFL in the first N_trans iterarions
  !> @param integer [in] n_iter : Number of iterations to go from a small CFL to
  !> the final CFL (in parameters.f90)
  !> @param real [in] current_time : Current (global) simulation time
  !> @param real [in] tprint : time for the next programed disk dump
  !> @param real [out] : @f$ \Delta t@f$ allowed by the CFL condition
  !> @param logical [out] dump_flag : Flag to write to disk
  subroutine get_timestep(current_iter, n_iter, current_time, tprint, dt,      &
                          dump_flag)

    use constants, only : SOLVER_HLLE_SPLIT_ALL, SOLVER_HLLD_SPLIT_ALL
    use parameters, only : nx, ny, nz, cfl, mpi_real_kind, mhd, riemann_solver
    use globals, only : primit, dx, dy, dz, primit0
#ifdef MPIP
    use mpi
#endif
    implicit none
    integer, intent(in)  :: current_iter, n_iter
    real,    intent(in)  :: current_time, tprint
    real,    intent(out) :: dt
    logical, intent(out) :: dump_flag
    real              :: dtp
    real              :: c
#ifdef BFIELD
    real              :: cx, cy, cz
#endif
    integer :: i, j, k, err

    dtp=1.e30

    do k=1,nz
      do j=1,ny
        do i=1,nx

          if (mhd) then
#ifdef BFIELD
            if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or. &
                riemann_solver == SOLVER_HLLD_SPLIT_ALL) then
              call cfast(primit(5,i,j,k)+primit0(5,i,j,k),primit(1,i,j,k)      &
                        + primit0(1,i,j,k),&
                 primit(6,i,j,k)+primit0(6,i,j,k), primit(7,i,j,k)             &
                +primit0(7,i,j,k), primit(8,i,j,k)+primit0(8,i,j,k),           &
                 cx,cy,cz)
            else
              call cfast(primit(5,i,j,k),primit(1,i,j,k),                      &
                   primit(6,i,j,k), primit(7,i,j,k), primit(8,i,j,k),          &
                   cx,cy,cz)
            end if

            dtp=min(dtp,dx/(abs(primit(2,i,j,k))+cx))
            dtp=min(dtp,dy/(abs(primit(3,i,j,k))+cy))
            dtp=min(dtp,dz/(abs(primit(4,i,j,k))+cz))
#endif
          else
            ! plain old hydro
            call csound(primit(5,i,j,k),primit(1,i,j,k),c)
            dtp=min(dtp,dx/(abs(primit(2,i,j,k))+c))
            dtp=min(dtp,dy/(abs(primit(3,i,j,k))+c))
            dtp=min(dtp,dz/(abs(primit(4,i,j,k))+c))
          end if

        end do
      end do
    end do

    if (current_iter <= n_iter ) then
      !   boots up simulation with small CFL
      dtp = cfl*(2.**(-real(n_iter+1-current_iter)))*dtp
    else
      dtp=cfl*dtp
    end if

#ifdef MPIP
    call mpi_allreduce(dtp, dt, 1, mpi_real_kind, mpi_min, mpi_comm_world,err)
#else
    dt=dtp
#endif

    !  Adjust dt so t+dt doesnt overshoot tprint
    if(( current_time + dt )>=tprint) then
      dt= tprint-current_time
      dump_flag=.true.
    endif

  end subroutine get_timestep

  !=======================================================================
  !> @brief Performs a linear reconstruction of the primitive variables
  !> @details returns a linear reconstruction of the variables at the
  !> interface beteen the primitives PLL, PL, PR, PRR
  !> @n The reconstruction is made with a slope limiter chosen at
  !> compilation time (i.e. set on the Makefile)
  !> @param real [in]    : primitives at the left of the left state
  !> @param real [inout] : primitives at the left state
  !> @param real [inout] : primitives at the right state
  !> @param real [in]    : primitives at the right of the right state
  !> @param real [in]    : number of equations
  subroutine limiter(PLL,PL,PR,PRR,neq)

    implicit none
    integer, intent (in)  :: neq
    real, dimension(neq), intent(inout) :: pl,  pr
    real, dimension(neq), intent(in)    :: pll, prr
    real :: dl, dm, dr, al, ar
    integer :: ieq
    real :: s, c, d, av1, av2
    real, parameter :: delta=1.e-7

    do ieq=1,neq
      dl=pl(ieq)-pll(ieq)
      dm=pr(ieq)-pl(ieq)
      dr=prr(ieq)-pr(ieq)
      al=average(dl,dm)
      ar=average(dm,dr)
      pl(ieq)=pl(ieq)+al*0.5
      pr(ieq)=pr(ieq)-ar*0.5
    end do

  contains

    real function average(a,b)
      use constants
      use parameters, only : slope_limiter
      implicit none
      real, intent(in)    :: a, b

      If (slope_limiter == LIMITER_NO_AVERAGE) then
        !   no average (reduces to 1st order)
        average=0.
      end if

      if (slope_limiter == LIMITER_NO_LIMIT) then
        !   no limiter
        average=0.5*(a+b)
      end if

      if (slope_limiter == LIMITER_MINMOD) then
        !   Minmod - most diffusive
        s=sign(1.,a)
        average=s*max(0.,min(abs(a),s*b))
      end if

      if (slope_limiter == LIMITER_VAN_LEER) then
        !   Falle Limiter (Van Leer)
        if(a*b.le.0.) then
          average=0.
        else
          average=a*b*(a+b)/(a**2+b**2)
        end if
      end if

      if (slope_limiter == LIMITER_VAN_ALBADA) then
        !   Van Albada
        average=(a*(b*b+delta)+b*(a*a+delta))/(a*a+b*b+delta)
      end if

      if (slope_limiter == LIMITER_UMIST) then
        !   UMIST Limiter - less diffusive
        s=sign(1.,a)
        c=0.25*a+0.75*b
        d=0.75*a+0.25*b
        average=min(2.*abS(a),2.*s*b,s*c,s*d)
        average=s*max(0.,average)
      end if

      if (slope_limiter == LIMITER_WOODWARD) then
        !    Woodward Limiter (o MC-limiter; monotonized central difference)
        s=sign(1.,a)
        c=0.5*(a+b)
        average=min(2.*abs(a),2.*s*b,s*c)
        average=s*max(0.,average)
      end if

      if (slope_limiter == LIMITER_SUPERBEE) then
        !   superbee Limiter (tends to flatten circular waves)
        s=sign(1.,b)
        av1=min(2.*abs(b),s*a)
        av2=min(abs(b),2.*s*a)
        average=s*max(0.,av1,av2)
      end if

    end function average

  end subroutine limiter

  !=======================================================================

end module hydro_core
