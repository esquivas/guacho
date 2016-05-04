!=======================================================================
!> @file linear_system.f90
!> @brief linear system inversion  module
!> @author A. Castellanos, A. Rodriguez, A. Raga and A. Esquivel
!> @date 4/May/2016

! Copyright (c) 2016 A. Esquivel et al.
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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief linear system inversion  module
!> @details Inversion of a system of linear equations with an LU
!> decomposition method (these routines are from Numerical
!> Methods by Press et al.)

module linear_system
  implicit none

  contains

!======================================================================  

!> @brief LU decomposition
!> @details LU decomposition of a row-wise permutation
!> @param real [inout] a(n,n) : matrix to be decomposed result is done in place
!> @param integer [in] n : size of the matrix
!> @param real [out] index(n) : vector that contains the row permutation
!> affected by the partial pivoting
!> @param integer [inout] d : +/- 1 depending if the row intergarches is
!> even or odd

subroutine ludcmp(a, n, indx,d)

  implicit none
  integer, intent(in) :: n
  real (kind=8), intent(inout) :: a(n,n), d
  integer, intent(out) :: indx(n)

  real (kind=8) :: vv(n*n), aamax, dum, sum
  real (kind=8), parameter  :: tiny=1.0e-25
  integer :: i, imax, j, k

  d=1.
  do i=1,n
    aamax=0.
    do j=1,n
      if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
    end do
    !        if(aamax.eq.0.) pause ' singular matrix.'
    vv(i)=1./aamax
  end do

  do j=1,n
    if(j.gt.1) then
      do i=1,j-1
        sum=a(i,j)
        if(i.gt.1) then
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          end do
          a(i,j)=sum
        end if
      end do
    end if
    aamax=0.
    do i=j,n
      sum=a(i,j)
      if(j.gt.1) then
        do k=1,j-1
          sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
      end if
      dum=vv(i)*abs(sum)
      if(dum.ge.aamax) then
        imax=i
        aamax=dum
      end if
    end do
    if(j.ne.imax) then
      do k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
      end do
      d=-d
      vv(imax)=vv(j)
    end if
    indx(j)=imax
    if(j.ne.n) then
      if(a(j,j).eq.0.) a(j,j)=tiny
      dum=1./a(j,j)
      do i=j+1,n
        a(i,j)=a(i,j)*dum
      end do
    end if
  end do
  if(a(n,n).eq.0.) a(n,n)=tiny
  return

end subroutine ludcmp

!=======================================================================

!> @brief Solves a set of linear equations
!> @details Solves a linear set of equations of the form @f$A cdot X = B$
!> with a LU decomposition method
!> @param real [inout] a(n,n) : LU decomposition of the matrix A
!> @param integer [in] n : size of the matrix
!> @param real [out] index(n) : vector that contains the row permutation
!> affected by the partial pivoting (from ludcmp)
!> @param real [inout] b : right hand side vector, the result is returned
!> in this same vector

subroutine lubksb(a,n,indx,b)
  
  implicit none
  integer, intent(in) :: n
  real (kind=8), intent(in)    :: a(n,n)
  real (kind=8), intent(inout) :: b(n)
  integer, intent(in) :: indx(n)
  integer :: i, ii, j, ll
  real (kind = 8) :: sum

  ii=0
  do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if(ii.ne.0) then
      do j=ii,i-1
        sum=sum-a(i,j)*b(j)
      end do
    else if(sum.ne.0.) then
      ii=i
    end if
    b(i)=sum
  end do
  do i=n,1,-1
    sum=b(i)
    if(i.lt.n) then
      do j=i+1,n
        sum=sum-a(i,j)*b(j)
      end do
    end if
    b(i)=sum/a(i,i)
  end do

  return

end subroutine lubksb

!=======================================================================

!> @brief Driver to solves a set of linear equations
!> @details Solves a linear set of equations   @f$A cdot X = B$
!> with an LU decomposition 
!> mehtod
!> @param real [inout] a(n,n) : the matrix A
!> @param real [inout] b : right hand side vector, the result is returned
!> in this same vector
!> @param integer [in] n : size of the system

subroutine linsys(a,b,n)
  implicit none
  integer, intent(in) :: n
  real (kind=8) :: a(n,n),b(n), d
  integer :: indx(n*n)
        
  call ludcmp(a,n,indx,d)
  call lubksb(a,n,indx,b)
        
  return

end subroutine linsys


end module linear_system

!======================================================================  