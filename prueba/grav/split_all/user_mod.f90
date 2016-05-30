!=======================================================================
!> @file user_mod.f90
!> @brief User input module
!> @author Alejandro Esquivel
!> @date 2/Nov/2014

! Copyright (c) 2014 A. Esquivel, M. Schneiter, C. Villareal D'Angelo
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

!> @brief User imput module
!> @details  This is an attempt to have all input neede from user in a
!! single file
!!!  This module should load additional modules (i.e. star, corona, sn), to
!!  impose initial and boundary conditions (such as sources)

module user_mod

  ! load auxiliary modules
  use corona   ! Note that parameters is included here
!   implicit none
!   real :: equi(neq,nymin:nymax)
!   public :: integral
contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty 
subroutine init_user_mod()

  implicit none      
  !  if needed initialize modules loaded by user
  call init_corona()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, &
                        nzmax, cv, rsc, vsc2, rhosc, Tempsc, ny, &
                        gamma, mu, mu_1, mu_2, nghost

  use globals, only : coords, dy, dx, primit0, rank
  use constants, only : Rg, amh, Ggrav, Msun, Rsun
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) 
  integer :: i, j, jj, k, err, ndim
  real :: Temp_ch, Temp_c,mu_var, rho, P_0, P_y, ym, xm, rad, g, y, bb
  real :: Tempi(1-nghost:nytot+nghost)

  mu_var = mu
  g = Ggrav*Msun/Rsun/Rsun
  nc = 1.e19
  Temp_ch = 1.e4
  Temp_c = 1.e6
  rho = nc*mu_var*amh
  P_0 = nc*amh*Rg*Temp_ch
  bb = 1.e-10
  
! VARIABLES NO PERTURBADAS
  primit0(2,:,:,:)=0.
  primit0(3,:,:,:)=0.
  primit0(4,:,:,:)=0.
  primit0(6,:,:,:)=0.
  primit0(7,:,:,:)=bb/Bsc
  primit0(8,:,:,:)=0.  
  
  if (rank.eq.master) then 
    do j= 1-nghost,nytot+nghost
      y = (float(j) + 0.5)*dy*rsc
      if (y.le.2.e8) then
        Tempi(j) = 1.e4
      else
        Tempi(j) = 1.e6
      end if
    end do
  end if
  
  ndim= nytot+2*nghost-1
#ifdef MPIP
  call mpi_bcast(Tempi,ndim,mpi_double_precision,0,mpi_comm_world,err)
#endif

    
  do j = nymin,nymax
     jj = j + coords(1)*ny
     y = (float(jj) + 0.5)*dy*rsc
     P_y = P_0*exp(-integral(jj,Tempi,mu_var))
     do k =nzmin, nzmax
        do i = nxmin, nxmax
           primit0(5,i,j,k) = P_y/Psc
           primit0(1,i,j,k) = (mu_var*P_y/Rg/Tempi(jj))/rhosc
        end do
     end do
  end do

! VARIABLES PERTURBADAS
  u(:,:,:,:)=0.
  
! PULSO
!   do j=nymin,nymax
!     ym = (float(j+coords(1)*ny-nytot/2) + 0.5)*dy*rsc
!     do i=nxmin,nxmax
!       xm = (float(i+coords(0)*nx-nxtot/2) + 0.5)*dx*rsc
!       rad = sqrt(xm*xm+ym*ym)
!       if(rad.le.0.5e8) then
! !         print*, 'entre'
! 	u(5,i,j,:) = cv*0.12666724448792355*3./Psc  ! E-cv*primit0(5)
!       else
! 	u(5,i,j,:) = 0.
!       end if
!     end do
!   end do
  
 
end subroutine initial_conditions
 
real function integral(jj,Tempi,mu_var)
  use globals, only:dy
  implicit none
  integer, intent(in) :: jj 
  real,intent(in) ::Tempi(nymin:nymax),mu_var
  integer :: ii
  real :: y, g
  
  g = Ggrav*Msun/Rsun/Rsun
  integral=0.
  do ii=0,jj
     y = (float(ii)+0.5)*dy*rsc
     integral = integral + dy*rsc*mu_var*g/Rg/Tempi(ii)  
  end do
end function integral
  
!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) : 
!! conserved variables
!> @param real [in] time : time in the simulation (code units)
!> @param integer [in] order : order (mum of cells to be filled in case
!> domain boundaries are being set, valid values are 1 or 2)

! subroutine impose_user_bc(u,order)
! 
!   use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
!                          bc_user, nx, ny, nz, tsc
!   use globals   , only : time 
!   implicit none
!   real :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
!   integer, intent(in) :: order
!   integer :: i, k, ieq
!   real    :: time_sec

!   time_sec = time*tsc  
! 
!   !  In this case the boundary is the same for 1st and second order)
!   
!   if (order == 1) then 
!      !  BC in X
! !!$     do k=0, nz+1
! !!$        do j=0, ny+1
! !!$           do ieq = 1, neq
! !!$              u(ieq,   0,j,k) = equi(ieq,j)
! !!$              u(ieq,nx+1,j,k) = equi(ieq,j)
! !!$           end do
! !!$        end do
! !!$     end do
!      !  BC in Y
!      do k=0, nz+1
!         do i=0, nx+1
!            do ieq = 1, neq
!            ! porque son cero???
!            !--------------------------------------
!               u(ieq,i,0   ,k) = 0.
!               u(ieq,i,ny+1,k) = 0.
! !               u(2:4,i,0   ,k) = 0. ! velocidad cero
! !               u(2:4,i,ny+1,k) = 0. ! velocidad cero
! !               u(1,i,0   ,k) = u(1,i,1 ,k) ! outflow para las demas 
! !               u(1,i,ny+1,k) = u(1,i,ny,k)
! !               do ieq= 5,neq
! !                 u(ieq,i,0   ,k) = u(ieq,i,1 ,k) 
! !                 u(ieq,i,ny+1,k) = u(ieq,i,ny,k)
! !               end do
!            end do
!         end do
!      end do
! !!     !  BC in Z
! !!     do j=0, ny+1
! !!        do i=0, nx+1
! !!           do ieq = 1, neq
! !!              u(ieq,i,j,0  ) =  equi(ieq,j)
! !!              u(ieq,i,j,nz+1) = equi(ieq,j)
! !!           end do
! !!        end do
! !!     end do
!   else if (order == 2) then
!           !  BC in X
! !     do k=nzmin, nzmax
! !        do j=nymin, nymax
! !           do ieq = 1, neq
! !              u(ieq,   -1,j,k) = equi(ieq,j)
! !              u(ieq,    0,j,k) = equi(ieq,j)
! !              u(ieq, nx+1,j,k) = equi(ieq,j)
! !              u(ieq, nx+2,j,k) = equi(ieq,j)
! !           end do
! !        end do
! !     end do
!      !  BC in Y
!      do k=nzmin, nzmax
!         do i=nxmin, nxmax
!            do ieq = 1, neq
!            ! porque son cero??
!            !---------------------------------------
!               u(ieq,i,  -1,k) = 0.
!               u(ieq,i,   0,k) = 0.
!               u(ieq,i,ny+1,k) = 0.
!               u(ieq,i,ny+2,k) = 0.
! !               u(2:4,i,  -1,k) = 0.
! !               u(2:4,i,0   ,k) = 0. ! velocidad cero
! !               u(2:4,i,ny+1,k) = 0. ! velocidad cero
! !               u(2:4,i,ny+2,k) = 0. ! velocidad cero
! !               u(1,i,-1  ,k) = u(1,i,1 ,k)
! !               u(1,i,0   ,k) = u(1,i,1 ,k) ! outflow para las demas 
! !               u(1,i,ny+1,k) = u(1,i,ny,k)
! !               u(1,i,ny+2,k) = u(1,i,ny,k)
! !               do ieq= 5,neq
! !                 u(ieq,i,-1  ,k) = u(ieq,i,1 ,k) 
! !                 u(ieq,i,0   ,k) = u(ieq,i,1 ,k) 
! !                 u(ieq,i,ny+1,k) = u(ieq,i,ny,k)
! !                 u(ieq,i,ny+2,k) = u(ieq,i,ny,k)
! !               end do
!            end do
!         end do
!      end do
!      !  BC in Z
! !!$     do j=nymin, nymax
! !!$        do i=nxmin, nxmax
! !!$           do ieq = 1, neq
! !!$              u(ieq,i,j,-1  ) =  equi(ieq,j)
! !!$              u(ieq,i,j,0   ) =  equi(ieq,j)
! !!$              u(ieq,i,j,nz+1) =  equi(ieq,j)
! !!$              u(ieq,i,j,nz+2) =  equi(ieq,j)
! !!$           end do
! !!$        end do
! !!$     end do
!      end if
  
! end subroutine impose_user_bc

!=======================================================================

!> @brief User Defined source terms
!> This is a generic interrface to add a source term S in the equation
!> of the form:  dU/dt+dF/dx+dG/dy+dH/dz=S
!> @param real [in] pp(neq) : vector of primitive variables
!> @param real [inout] s(neq) : vector with source terms, has to add to
!>  whatever is there, as other modules can add their own sources
!> @param integer [in] iin : cell index in the X direction
!> @param integer [in] jin : cell index in the Y direction
!> @param integer [in] kin : cell index in the Z direction

subroutine get_user_source_terms(pp,s, iin, jin , kin,prim0)

  ! in this example a constant gravity is added
  use constants,  only : Ggrav,Msun,Rsun
  use parameters, only : neq, nymin, nymax, rsc, vsc
  implicit none
  integer, intent(in), optional    :: iin, jin, kin
  real, dimension(neq), intent(in), optional :: prim0
  integer  :: i,j,k
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)
  real :: g

  if (present(iin)) i=iin
  if (present(jin)) j=jin
  if (present(kin)) k=kin

  g = Ggrav*Msun/Rsun/Rsun
  g = g*rsc/vsc2

  ! momento
  s(2)= s(2)
  s(4)= s(4)
  s(3)= s(3)-pp(1)*g!GM/((R)**2.)
  ! energy
  if (riemann_solver == SOLVER_HLLE_SPLIT_ALL .or.   &
      riemann_solver == SOLVER_HLLD_SPLIT_ALL .or.   &
      riemann_solver == SOLVER_HLL_SPLIT_ALL  .or.   &
      riemann_solver == SOLVER_HLLC_SPLIT_ALL) then
     if (present(prim0)) then
       s(5)= s(5)-(pp(1)+prim0(1))*pp(3)*g!GM/((R)**2. ) 
     end if
  else
    s(5)= s(5)-pp(1)*pp(3)*g
  end if

end subroutine get_user_source_terms
!=====================================================================

end module user_mod



!=======================================================================
