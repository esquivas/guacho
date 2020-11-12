!=======================================================================
!> @file difrad.f90
!> @brief Diffuse radiation module
!> @author Alejandro Esquivel
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

!> @brief Ray tracing Radiative Trasnport
!> @details Ray tracing Radiative Trasnport
module difrad

  use parameters
  use globals
  implicit none
  real, parameter    :: a0=6.3e-18    !< Fotoionization cross section
  integer, parameter :: nrays=1000000 !< Number of rays
  real, allocatable  :: ph(:,:,:)     !< Photoionizing rate
  real, allocatable  :: em(:,:,:)     !< Photoionizing emissivity
  !  auxiliary MPI arrays
  real, allocatable  :: photL(:,:,:)  !< Auxiliary buffer for MPI
  real, allocatable  :: photR(:,:,:)  !< Auxiliary buffer for MPI
  real, allocatable  :: photB(:,:,:)  !< Auxiliary buffer for MPI
  real, allocatable  :: photT(:,:,:)  !< Auxiliary buffer for MPI
  real, allocatable  :: photO(:,:,:)  !< Auxiliary buffer for MPI
  real, allocatable  :: photI(:,:,:)  !< Auxiliary buffer for MPI
  integer            :: buffersize(6) !< Auxiliary buffer for MPI

contains

  !=======================================================================
  !> @brief initializes random number generation
  !> @details initializes random number generation
  subroutine init_difrad()

    implicit none
    integer :: rand_size
    integer, allocatable, dimension(:) :: rand_seed
    character (len=10) :: system_time
    real :: rtime

    allocate( ph(nx,ny,nz) )
    allocate( em(nx,ny,nz) )

    !allocate buffers for transmission of photons
    !  1st index 1:send 2:recv
    !  2nd index 3xposition + 3xdirection + 1for the F
    !  3rd index a large number can be changed if the code exits
    !  on error
    allocate (photL(2,7,nrays/5))
    allocate (photR(2,7,nrays/5))
    allocate (photB(2,7,nrays/5))
    allocate (photT(2,7,nrays/5))
    allocate (photO(2,7,nrays/5))
    allocate (photI(2,7,nrays/5))

    call random_seed(size=rand_size)
    allocate(rand_seed(1:rand_size))
    call date_and_time(time=system_time)
    read(system_time,*) rtime
    rand_seed=int(rtime*1000.0)
#ifdef MPIP
    rand_seed=rand_seed*rank
#endif
    call random_seed(put=rand_seed)
    deallocate(rand_seed)

  end subroutine init_difrad

  !=======================================================================
  !> @brief calculates the diffuse fotoionization emissivity
  !> @brief calculates the diffuse fotoionization emissivity in the
  !! entire domain
  !> @param real [out] emax : maximum emissivity in the entire grid
  subroutine emdiff(emax)

    use hydro_core, only : u2prim
    implicit none
    real, intent(out) :: emax
    real :: prim(neq)
    real :: T, emaxp, de, emLym, emHeII
    integer :: i ,j , k, err
    !real :: x,y,z,rad,rstar,vol,emstar

    !   clear  variables
    emaxp=0.

    do k=1,nz
      do j=1,ny
        do i=1,nx
          call u2prim(u(:,i,j,k), prim, T)
          !  electron density
          de=max(u(1,i,j,k)-u(neqdyn+1,i,j,k),0.)
          !  Lyman cont. emission
          emLym= de*de*1.57e-13*(1.e4/T)**0.52

          !  He II emission (turned off)
          emHeII= 0. ! 2.*0.1*de*de*0.101*8.629e-6*exp(-4.72e5/T)/sqrt(T)

          !  Total emissivity
          em(i,j,k) = emLym +emHeII
          !
          emaxp=max(emaxp, emLym+emHeII)

        end do
      end do
    end do

#ifdef MPIP
    call mpi_allreduce(emaxp, emax, 1, mpi_real_kind, mpi_max, comm3d,err)
#else
    emax=emaxp
#endif

  end subroutine emdiff

  !=======================================================================
  !> @brief returns the 3 components of a random versor
  !> @details returns the 3 components of a random versor (unit magnitude)
  !> @param real [out] xd : x component
  !> @param real [out] yd : y component
  !> @param real [out] zd : z component
  subroutine random_versor(xd,yd,zd)

    implicit none
    real, intent (out) :: xd, yd, zd
    real :: r2, w
    real :: ran(3)

    r2=2.
    do while (r2.gt.1.)
       call random_number(ran)
       xd=2.*(ran(1)-0.5)
       yd=2.*(ran(2)-0.5)
       zd=2.*(ran(3)-0.5)
       r2=xd**2+yd**2+zd**2
    end do
    w=1./sqrt(r2)
    xd=xd*w
    yd=yd*w
    zd=zd*w

  end subroutine random_versor

  !=======================================================================
  !> @brief Place photon packets at a "star" surface
  !> @details returns the random location and direction at a star surface,
  !! if the direction goes into the star, the direction is inverted
  !> @param real [in] Srad : radius of the "star"
  !> @param real [in] x0   : X position of the center of the star
  !> @param real [in] y0   : Y position of the center of the star
  !> @param real [in] y0   : Z position of the center of the star
  !> @param real [out] x   : random X position at the star surface
  !> @param real [out] y   : random Y position at the star surface
  !> @param real [out] z   : random Z position at the star surface
  !> @param real [out] xd  : random X direction
  !> @param real [out] yd  : random Y direction
  !> @param real [out] zd  : random Z direction
  subroutine starsource(srad,x0,y0,z0,x,y,z,xd,yd,zd)

    implicit none
    real, intent (in)  :: Srad, x0, y0, z0
    real, intent (out) :: x, y, z, xd, yd, zd
    real :: xs, ys,zs
    real :: rd, rs, w, prod
    real :: ran(3)

    !   get unit vector from star center
    rs=2.
    do while(rs > 1.)
       call random_number(ran)
       xs=2.*(ran(1)-0.5)
       ys=2.*(ran(2)-0.5)
       zs=2.*(ran(3)-0.5)
       rs= xs**2 + ys**2 + zs**2
    end do
    w=Srad/sqrt(rs)
    xs=xs*w
    ys=ys*w
    zs=zs*w

    !  the position of the source
    x=x0+xs
    y=y0+ys
    z=z0+zs

    !  get the direction of the ray
    rd=2.
    do while(rd > 1.)
       call random_number(ran)
       xd=2.*(ran(1)-0.5)
       yd=2.*(ran(2)-0.5)
       zd=2.*(ran(3)-0.5)
       rd = xd**2 + yd**2 + zd**2
    end do
    w=1./sqrt(rd)
    xd=xd*w
    yd=yd*w
    zd=zd*w

    !  reflect rays that point to the star
    prod=( xs*xd+ ys*yd + zs*zd )
    if (prod < 0. ) then
       xd=-xd
       yd=-yd
       zd=-zd
    endif

  end subroutine starsource

  !=======================================================================
  !> @brief Photon trajectories
  !> @details Launches a photon from cell (xc,yc,zc) in the (xd,yd,zd)
  !! direction, with f and ionizing photons, and updates the
  !! photoionizing rate
  !> @param real [in] xl0 : Initial X position
  !> @param real [in] yl0 : Initial Y position
  !> @param real [in] zl0 : Initial Z position
  !> @param real [in] xd :  Direction in X
  !> @param real [in] yd :  Direction in Y
  !> @param real [in] zd :  Direction in Z
  !> @param real [in] f  : NUmber of photoionizong photons
  subroutine photons(xl0,yl0,zl0,xd,yd,zd,f)

    implicit none
    real, intent(in) :: xl0, yl0, zl0, xd,yd,zd
    real, intent(inout) :: f
    real :: dl, dxl, dyl, dzl, xl, yl, zl, dtau
    real :: fmin
    integer :: i, j ,k

    !     photon trajectory
    fmin=1.e-10*f
    dl=0.5
    dxl=dl*xd
    dyl=dl*yd
    dzl=dl*zd

    xl=xl0+dxl
    yl=yl0+dyl
    zl=zl0+dzl

    i=int(xl + 0.5)
    j=int(yl + 0.5)
    k=int(zl + 0.5)

    if (i < 1) then
      buffersize(1)=buffersize(1)+1
      photL(1,1, buffersize(1))= xl0+real(nx)
      photL(1,2, buffersize(1))= yl0
      photL(1,3, buffersize(1))= zl0

      photL(1,4, buffersize(1))= xd
      photL(1,5, buffersize(1))= yd
      photL(1,6, buffersize(1))= zd
      photL(1,7, buffersize(1))= f
      !return
    else if(i > nx) then
      buffersize(2)=buffersize(2)+1
      photR(1,1, buffersize(2))= xl0-real(nx)
      photR(1,2, buffersize(2))= yl0
      photR(1,3, buffersize(2))= zl0

      photR(1,4, buffersize(2))= xd
      photR(1,5, buffersize(2))= yd
      photR(1,6, buffersize(2))= zd
      photR(1,7, buffersize(2))= f
      !return
    endif
    if (j < 1) then
      buffersize(3)=buffersize(3)+1
      photB(1,1, buffersize(3))= xl0
      photB(1,2, buffersize(3))= yl0+real(ny)
      photB(1,3, buffersize(3))= zl0

      photB(1,4, buffersize(3))= xd
      photB(1,5, buffersize(3))= yd
      photB(1,6, buffersize(3))= zd
      photB(1,7, buffersize(3))= f
      !return
    else if(j > ny) then
      buffersize(4)=buffersize(4)+1
      photT(1,1, buffersize(4))= xl0
      photT(1,2, buffersize(4))= yl0-real(ny)
      photT(1,3, buffersize(4))= zl0

      photT(1,4, buffersize(4))= xd
      photT(1,5, buffersize(4))= yd
      photT(1,6, buffersize(4))= zd
      photT(1,7, buffersize(4))= f
      !return
    endif
    if (k < 1) then
      buffersize(5)=buffersize(5)+1
      photO(1,1, buffersize(5))= xl0
      photO(1,2, buffersize(5))= yl0
      photO(1,3, buffersize(5))= zl0+real(nz)

      photO(1,4, buffersize(5))= xd
      photO(1,5, buffersize(5))= yd
      photO(1,6, buffersize(5))= zd
      photO(1,7, buffersize(5))= f
      !return
    else if(k > nz) then
      buffersize(6)=buffersize(6)+1
      photI(1,1, buffersize(6))= xl0
      photI(1,2, buffersize(6))= yl0
      photI(1,3, buffersize(6))= zl0-real(nz)

      photI(1,4, buffersize(6))= xd
      photI(1,5, buffersize(6))= yd
      photI(1,6, buffersize(6))= zd
      photI(1,7, buffersize(6))= f
      !return
    end if

    !  the photon will cross the domain boundaries

    do while( f.ge.fmin  &
        .and. (i <= nx).and.(i >= 1) &
        .and. (j <= ny).and.(j >= 1) &
        .and. (k <= nz).and.(k >= 1) &
        )

       dtau=u(neqdyn+1,i,j,k)*a0*dl*dx*rsc
       if (dtau < 1E-5) then
          ph(i,j,k)=f*a0*dl/((dx*rsc)**2)!*dx*rsc*/dx**3/rsc**3
          f=(1.-dtau)*f
       else
          ph(i,j,k)=ph(i,j,k)+f*(1.-exp(-dtau) )/(u(neqdyn+1,i,j,k)*(dx*rsc)**3)
          f=f*exp(-dtau)
       end if

       xl=xl+dxl
       yl=yl+dyl
       zl=zl+dzl

       i=int(xl+0.5)
       j=int(yl+0.5)
       k=int(zl+0.5)

       !  the photon reaches the domain boundaries
       if (i < 1) then
         buffersize(1)=buffersize(1)+1
         photL(1,1, buffersize(1))= xl-dxl+real(nx)
         photL(1,2, buffersize(1))= yl-dyl
         photL(1,3, buffersize(1))= zl-dzl

         photL(1,4, buffersize(1))= xd
         photL(1,5, buffersize(1))= yd
         photL(1,6, buffersize(1))= zd
         photL(1,7, buffersize(1))= f
         return
       else if(i > nx) then
         buffersize(2)=buffersize(2)+1
         photR(1,1, buffersize(2))= xl-dxl-real(nx)
         photR(1,2, buffersize(2))= yl-dyl
         photR(1,3, buffersize(2))= zl-dzl

         photR(1,4, buffersize(2))= xd
         photR(1,5, buffersize(2))= yd
         photR(1,6, buffersize(2))= zd
         photR(1,7, buffersize(2))= f
         return
       endif
       if (j < 1) then
         buffersize(3)=buffersize(3)+1
         photB(1,1, buffersize(3))= xl-dxl
         photB(1,2, buffersize(3))= yl-dyl+real(ny)
         photB(1,3, buffersize(3))= zl-dzl

         photB(1,4, buffersize(3))= xd
         photB(1,5, buffersize(3))= yd
         photB(1,6, buffersize(3))= zd
         photB(1,7, buffersize(3))= f
         return
       else if(j > ny) then
         buffersize(4)=buffersize(4)+1
         photT(1,1, buffersize(4))= xl-dxl
         photT(1,2, buffersize(4))= yl-dyl-real(ny)
         photT(1,3, buffersize(4))= zl-dzl

         photT(1,4, buffersize(4))= xd
         photT(1,5, buffersize(4))= yd
         photT(1,6, buffersize(4))= zd
         photT(1,7, buffersize(4))= f
         return
       endif
       if (k < 1) then
         buffersize(5)=buffersize(5)+1
         photO(1,1, buffersize(5))= xl-dxl
         photO(1,2, buffersize(5))= yl-dyl
         photO(1,3, buffersize(5))= zl-dzl+real(nz)

         photO(1,4, buffersize(5))= xd
         photO(1,5, buffersize(5))= yd
         photO(1,6, buffersize(5))= zd
         photO(1,7, buffersize(5))= f
         return
       else if(k > nz) then
         buffersize(6)=buffersize(6)+1
         photI(1,1, buffersize(6))= xl-dxl
         photI(1,2, buffersize(6))= yl-dyl
         photI(1,3, buffersize(6))= zl-dzl-real(nz)

         photI(1,4, buffersize(6))= xd
         photI(1,5, buffersize(6))= yd
         photI(1,6, buffersize(6))= zd
         photI(1,7, buffersize(6))= f
         return
       end if

     end do

   end subroutine photons

   !=======================================================================
   !> @brief follows the rays across MPI boundaries
   !> @details follows the rays across MPI boundaries
   subroutine radbounds()

#ifdef MPIP
    implicit none
    integer :: ip
    integer :: AllBufferSize(np*6), sizeSend, sizeRecv, niter
    integer,  parameter :: length=np*6
    integer:: status(MPI_STATUS_SIZE), err

    !   loop over MPI blocks to ensure the rays are followed to the
    !   entire domain
    do ip=1,int( sqrt(real(MPI_NBX)**2+real(MPI_NBY)**2+real(MPI_NBZ)**2) )

      !print'(a,2i3, 6i7)','**',rank,np,buffersize
      call mpi_allgather(buffersize(:)   , 6, mpi_integer,                     &
                         Allbuffersize(:), 6, mpi_integer, comm3d ,err)

      !  to and from the left
      if (left /= -1) then
        sizeSend=Allbuffersize(6*rank+1)
        if(sizeSend > 0) then
          call mpi_send(photL(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind,      &
                        left,100,comm3d,err)
          !print*,rank,'send -> L:',SizeSend,'to   ',left
        end if
        sizeRecv=Allbuffersize(6*left+2)
        if(sizeRecv > 0) then
          call mpi_recv(photL(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind,     &
                        left,102, comm3d, status,err)
           !print*,rank,'recv <- L:',SizeRecv,'from ',left
         end if
       end if

       !  from and to the right
       if (right /= -1) then
         sizeRecv=Allbuffersize(6*right+1)
         if(sizeRecv > 0 ) then
           call mpi_recv(photR(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind,    &
                         right,100,comm3d, status,err)
           !print*,rank, 'recv <- R:',SizeRecv,'from ',right
         end if
         sizeSend=Allbuffersize(6*rank +2)
         if(sizeSend > 0) then
           call mpi_send(photR(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind,     &
                         right,102,comm3d,err)
           !print*,rank,'send -> R:',SizeSend,'to   ',right
         end if
       end if

       !  to and from the bottom
       if (bottom /= -1) then
         sizeSend=Allbuffersize(6*rank+3)
         if(sizeSend > 0) then
           call mpi_send(photB(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind,     &
                         bottom,200,comm3d,err)
           !print*,rank,'send -> B:',SizeSend,'to   ',bottom
        end if
        sizeRecv=Allbuffersize(6*bottom+4)
        if(sizeRecv > 0) then
          call mpi_recv(photB(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind,     &
                        bottom,202, comm3d, status,err)
          !print*,rank,'recv <- B:',SizeRecv,'from ',bottom
        end if
      end if

      !  from and to the top
      if (top /= -1) then
        sizeRecv=Allbuffersize(6*top+3)
        if(sizeRecv > 0 ) then
          call mpi_recv(photT(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind,     &
                        top,200,comm3d, status,err)
          !print*,rank, 'recv <- T:',SizeRecv,'from ',top
        end if
        sizeSend=Allbuffersize(6*rank +4)
        if(sizeSend > 0) then
          call mpi_send(photT(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind,      &
                        top,202,comm3d,err)
          !print*,rank,'send -> T:',SizeSend,'to   ',top
        end if
      end if

      !  to and from the out
      if (out /= -1) then
        sizeSend=Allbuffersize(6*rank+5)
        if(sizeSend > 0) then
          call mpi_send(photO(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind,      &
                        out,300,comm3d,err)
          !print*,rank,'send -> O:',SizeSend,'to   ',out
        end if
        sizeRecv=Allbuffersize(6*out+6)
        if(sizeRecv > 0) then
          call mpi_recv(photO(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind,   &
                        out,302, comm3d, status,err)
          !print*,rank,'recv <- O:',SizeRecv,'from ',out
        end if
      end if

      !  from and to the in
      if (in /= -1) then
        sizeRecv=Allbuffersize(6*in+5)
        if(sizeRecv > 0 ) then
          call mpi_recv(photI(2,1:7,1:sizeRecv),7*sizeRecv, mpi_real_kind,     &
                        in,300,comm3d, status,err)
          !print*,rank, 'recv <- I:',SizeRecv,'from ',in
        end if
        sizeSend=Allbuffersize(6*rank +6)
        if(sizeSend > 0) then
          call mpi_send(photI(1,1:7,1:sizeSend),7*sizeSend,mpi_real_kind,      &
                        in,302,comm3d,err)
          !print*,rank,'send -> I:',SizeSend,'to   ',in
        end if
      end if

      !  Continue with photon tracing
      !  Reset the out buffers of individual MPI blocks
      buffersize(:)=0

      !  left face
      if(left /= -1) then
        do niter=1,Allbuffersize(6*left+2)
          call photons(photL(2,1,niter),photL(2,2,niter),photL(2,3,niter), &
          photL(2,4,niter),photL(2,5,niter),photL(2,6,niter),photL(2,7,niter) )
        end do
      endif

      !  right face
      if(right /= -1) then
        do niter=1,Allbuffersize(6*right+1)
          call photons(photR(2,1,niter),photR(2,2,niter),photR(2,3,niter), &
          photR(2,4,niter),photR(2,5,niter),photR(2,6,niter),photR(2,7,niter) )
        end do
      end if

      !  bottom face
      if(bottom /= -1) then
        do niter=1,Allbuffersize(6*bottom+4)
          call photons(photB(2,1,niter),photB(2,2,niter),photB(2,3,niter), &
          photB(2,4,niter),photB(2,5,niter),photB(2,6,niter),photB(2,7,niter) )
        end do
      end if

      !  top face
      if(top /= -1) then
        do niter=1,Allbuffersize(6*top+3)
          call photons(photT(2,1,niter),photT(2,2,niter),photT(2,3,niter), &
          photT(2,4,niter),photT(2,5,niter),photT(2,6,niter),photT(2,7,niter) )
        end do
      end if

      !  out face
      if(out /= -1) then
        do niter=1,Allbuffersize(6*out+6)
          call photons(photO(2,1,niter),photO(2,2,niter),photO(2,3,niter), &
          photO(2,4,niter),photO(2,5,niter),photO(2,6,niter),photO(2,7,niter) )
        end do
      end if

      !  in face
      if(in /= -1) then
        do niter=1,Allbuffersize(6*in+5)
          call photons(photI(2,1,niter),photI(2,2,niter),photI(2,3,niter), &
          photI(2,4,niter),photI(2,5,niter),photI(2,6,niter),photI(2,7,niter) )
        end do
      end if

      !  reset the remaining out buffers
      allBuffersize(:)=0

    end do
#endif

  end subroutine radbounds

  !=======================================================================
  !> @brief Progress bar
  !> @details Progress bar (only tested with Fortran conmpiler)
  !! takes a number between 1 and tot
  !> @param integer [in] j   : current iteration
  !> @param integer [in] tot : total number of iterartions
  subroutine progress(j,tot)

    implicit none
    integer(kind=4)::j,k
    integer(kind=4), intent(in) :: tot
    character(len=57)::bar="???% |                                                  |"
    open (unit=6)
    write(unit=bar(1:3),fmt="(i3)") 100*j/tot
    bar(7:56)="."
    do k=1, 50*j/tot
      bar(6+k:6+k)="="
    enddo
    ! print the progress bar.
    write(unit=6,fmt="(a1,a1,a57)",advance="no") '+',char(13), bar

  end subroutine progress

  !=======================================================================
  !> @brief  Diffuse radiation driver
  !> @details Upper level wrapper to compute the diffuse photoionization
  !>rate
  subroutine diffuse_rad()

    use constants, only : Rsun
    implicit none
    real :: f, dirx,diry,dirz
    real :: Srad, xc,yc,zc, xp, yp, zp
    integer :: i,j,k, err, niter, ii, jj, kk
    integer :: in, nmax

    !nmax = nxtot*nytot*nztot/10000/np
    !resets the photoionization rate
    ph (:,:,:)=0.
    em (:,:,:)=0.
    buffersize(:)=0

    !   computes the emissivity at each cell
    !call emdiff(emax)

    !   fire the photon torpedoes!

    !  Source position
    xc=real(nxtot/2)*dx
    yc=real(nytot/2)*dy
    zc=real(nztot/2)*dz

    ! To  impose the photoionizing field of a star
    ! this is the radius of the actual star [cm]
    srad=1.1*Rsun/rsc

    in=0  ! number or rays successfully injected
    do niter=1, nrays
      !  get the location and direction of the photon to be traced
      call starsource(srad,xc,yc,zc,xp,yp,zp,dirx,diry,dirz)
      !  obtain in which proc will be put
      !
      f=2.46E38   !  flujo de fotones ionizantes:not divided by np
      !f=nphot
      !f=5.E48  !  Iliev 2009 Test 5

      i=int(xp/dx - 0.5)
      j=int(yp/dy - 0.5)
      k=int(zp/dz - 0.5)

      ii=i/nx
      jj=j/ny
      kk=k/nz

      if( (ii==coords(0)).and.(jj==coords(1)).and.(kk==coords(2))) then
        in=in+1
        ! trace the photon
        call photons(real(i-coords(0)*nx)+0.5, &
                     real(j-coords(1)*ny)+0.5, &
                     real(k-coords(2)*nz)+0.5, &
                     dirx,diry,dirz,f)
        !call progress(niter,nrays)
      end if

    end do

    !determine the actual number of photons injected,
    !and divide ph by it
#ifdef MPIP
    call mpi_allreduce(in, nmax, 1, mpi_integer, mpi_sum, mpi_comm_world,err)
#else
    nmax=in
#endif
    ! trace photons across boundaries)
    call radbounds()

    em=ph
    ph(:,:,:)=ph(:,:,:)/real(nmax)

    return
  end subroutine diffuse_rad

  !=======================================================================

end module difrad
