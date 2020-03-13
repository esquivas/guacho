! ======================================================================
!  This module solves the Poisson equation to obtian the gravitational
!  potential and with it the force due to self-gravity
!  The Poisson equation is solved with a Successive Overrelaxation (SOR)
!  method (adapted from Numerical recipes, Press et al.)
module self_gravity
  implicit none

  !  gravitational potential before and after iterartion
  real, allocatable, save :: phi(:,:,:)
  integer :: nx, ny, nz
  real    :: dx, dy, dz

  real, parameter :: pi = acos(-1.)
  real, parameter :: Ggrav = 6.67259e-8   !  Gravitational const. [cgs]
  real, parameter :: four_pi_G = 4.*pi*Ggrav

contains

  ! ======================================================================
  !  Initializes variables
  subroutine init_SOR(nx_ext,ny_ext,nz_ext,dx_ext, dy_ext, dz_ext)
    implicit none
    integer, intent(in) :: nx_ext, ny_ext, nz_ext
    real,    intent(in) :: dx_ext, dy_ext, dz_ext

    nx=nx_ext
    ny=ny_ext
    nz=nz_ext

    dx=dx_ext
    dy=dy_ext
    dz=dz_ext

    allocate(phi (nx,ny,nz))
    phi (:,:,:)=0

  end subroutine init_SOR

  ! ======================================================================
  ! Computes the residue $\chi= \nabla^2 \phi - 4\pi G \rho$
  ! at a given cell
  subroutine get_residue(rho,phi,i,j,k,dx,dy,dz,residue)
    implicit none
    real,    intent(in)  :: rho(nx,ny,nz),phi(nx,ny,nz)
    integer, intent(in)  :: i,j,k
    real,    intent(in)  :: dx, dy, dz
    real,    intent(out) :: residue

    residue= ( phi(i+1, j,  k   )+phi(i-1,j,  k  ) -2.*phi(i,j,k) )/dx**2      &
            +( phi(i  , j+1,k   )+phi(i  ,j-1,k  ) -2.*phi(i,j,k) )/dy**2      &
            +( phi(i  , j , k+1 )+phi(i  ,j  ,k-1) -2.*phi(i,j,k) )/dz**2      &
            - four_pi_G*rho(i,j,k)

  end subroutine get_residue

  ! ======================================================================
  ! Computes phi*  assumes that dx=dy=dz
  subroutine get_phi_star(rho,phi,i,j,k,dx,phi_star)
    implicit none
    real,    intent(in)  :: rho(nx,ny,nz),phi(nx,ny,nz)
    integer, intent(in)  :: i,j,k
    real,    intent(in)  :: dx
    real,    intent(out) :: phi_star

    phi_star = ( phi(i+1,j,  k  )+phi(i-1,j  ,k  )                             &
               + phi(i  ,j+1,k  )+phi(i  ,j-1,k  )                             &
               + phi(i  ,j  ,k+1)+phi(i  ,j  ,k-1)                             &
               -four_pi_G*rho(i,j,k)*dx**2        )      / 6.

  end subroutine get_phi_star

  ! ======================================================================
  ! SOR method

  subroutine SOR(nx,ny,nz,dx,dy,dz,rho)
    implicit none
    integer, intent(in) :: nx, ny, nz  ! Number of physical cells
    real,    intent(in) :: dx, dy, dz
    real,    intent(in) :: rho (nx,ny,nz)   ! density
    ! Maximum number of iterations
    integer, parameter :: max_iterations=100
    real, parameter    :: Tol = 1E-4  !  Relative error tolerance
    real               :: omega, relative_error
    real    :: residue, max_error, e_ijk, ph0
    logical :: need_more=.false.
    integer :: iter, ipass, i,j,k, isw, jsw, kpass
    !real :: phiP(nx,ny,nz)

    omega=1.99

    e_ijk = -2.*( 1./dx**2 +1./dy**2 + 1./dz**2 )

    main_loop : do iter=1, max_iterations


      isw=1
      black_red: do ipass=1,2

        max_error=-10.

        do kpass=1,2
          jsw=isw
          do i=2, nx-1
            do j=jsw+1,ny-1,2
              do k=kpass+1,nz-1,2

                call get_residue(rho,phi,i,j,k,dx,dy,dz,residue)
                !relative_error= omega*abs(residue)/abs(phi(i,j,k)*e_ijk)
                ph0=phi(i,j,k)
                phi(i,j,k)=ph0-omega*residue/e_ijk

                !call get_phi_star(rho,phi,i,j,k,dx,residue)
                !ph0=phi(i,j,k)
                !phi(i,j,k)=omega*residue + (1.-omega)*ph0

                relative_error=abs(phi(i,j,k)-ph0)/abs(ph0)
                max_error=max(max_error,relative_error)

                !if(relative_error > Tol)
                need_more=.true.

              end do
            end do
            jsw=3-jsw
          end do
          isw=3-isw
        end do

      end do black_red
      phi(: ,1 ,: ) = phi(:   ,ny-1,:   )
      phi(: ,ny,: ) = phi(:   ,2   ,:   )
      phi(: ,: ,1 ) = phi(:   ,:   ,nz-1)
      phi(: , :,nz) = phi(:   ,:   ,2   )

      if(.not.need_more) then
        print*, 'Converged in ', iter, 'iterations',max_error
        return
      end if
      !  reset convergence flag
      need_more=.false.

    end do main_loop

    print'(a)', 'SOR exceeded maximum number of iterations'

  end subroutine SOR

  ! ======================================================================

end module self_gravity

! ======================================================================
!  utilities to testdrive the module

module utilities

  implicit none

contains

  subroutine read_data(nx,ny,nz,output_number,rho)
    implicit none
    integer, intent(in) :: nx,ny,nz, output_number
    real, intent (out)  :: rho(nx,ny,nz)
    real :: xpos, ypos, zpos
    integer :: i,j,k
    character (len=128) file_in

    write(file_in,'(a,i4.4,a)') 'DAT/out-',output_number,'.dat'
    open(unit=10,file=file_in,status='unknown')

    do i=1,nx
      do j=1,ny
        do k=1,nz
          read(10,'(4e16.7e3)') xpos, ypos, zpos, rho(i,j,k)
        end do
      end do
    end do

    close(10)

  end subroutine read_data

  ! ======================================================================

  subroutine  write_data(nx,ny,nz,itprint,data,name)

    implicit none
    integer, intent(in) :: nx, ny,nz,itprint
    real, intent(in) :: data(nx,ny,nz)
    character (len=128) :: file_out
    character (len=3)  :: name

    write(file_out,'(a,a,a,i4.4,a)') 'DAT/',trim(name),'-',itprint,'.bin'
    open(unit=11,file=file_out,status='unknown',form='unformatted', &
    convert='LITTLE_ENDIAN')

    write (11) data(:,:,:)
    close(11)
    print'(a,a)'," wrote file:",trim(file_out)

  end subroutine write_data

  ! =====================================================================
  !  computes laplacian
  subroutine get_laplacian(nx,ny,nz,dx,dy,dz,phi,laplacian)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real,    intent(in) :: dx, dy, dz
    real,    intent(in) :: phi(nx,ny,nz)
    real,    intent(out):: laplacian(nx,ny,nz)
    integer :: i,j,k

    do i=2,nx-1
      do j=2,ny-1
        do k=2,nz-1
          laplacian(i,j,k) =                                                   &
              ( phi(i+1,j  ,k  ) + phi(i-1,j  ,k  ) - 2.*phi(i,j,k) ) / dx**2  &
            + ( phi(i  ,j+1,k  ) + phi(i  ,j-1,k  ) - 2.*phi(i,j,k) ) / dy**2  &
            + ( phi(i  ,j  ,k+1) + phi(i  ,j  ,k-1) - 2.*phi(i,j,k) ) / dz**2
        end do
      end do
    end do

    laplacian(1 ,: ,: ) = laplacian(2   ,:   ,:   )
    laplacian(nx,: ,: ) = laplacian(nx-1,:   ,:   )
    laplacian(: ,1, : ) = laplacian(:   ,2   ,:   )
    laplacian(: ,ny,: ) = laplacian(:   ,ny-1,:   )
    laplacian(: ,: ,1 ) = laplacian(:   ,:   ,2   )
    laplacian(: ,: ,nz) = laplacian(:   ,:   ,nz-1)

  end subroutine get_laplacian

  ! =====================================================================
  !  computes gradient

  subroutine get_gradient(nx,ny,nz,dx,dy,dz,phi,gradient)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real,    intent(in) :: dx, dy, dz
    real,    intent(in) :: phi(nx,ny,nz)
    real,    intent(out):: gradient(nx,ny,nz)
    integer :: i,j,k

    do i=2,nx-1
      do j=2,ny-1
        do k=2,nz-1
          gradient(i,j,k) =  ( phi(i+1 ,j ,k  )-phi(i-1,j  ,k  ) )/ (2.*dx)    &
                           + ( phi(i  ,j+1,k  )-phi(i  ,j-1,k  ) )/ (2.*dy)    &
                           + ( phi(i  ,j  ,k+1)-phi(i  ,j  ,k-1) )/ (2.*dz)
        end do
      end do
    end do

    gradient(1 ,: ,: ) = gradient(2   ,:   ,:   )
    gradient(nx,: ,: ) = gradient(nx-1,:   ,:   )
    gradient(: ,1, : ) = gradient(:   ,2   ,:   )
    gradient(: ,ny,: ) = gradient(:   ,ny-1,:   )
    gradient(: ,: ,1 ) = gradient(:   ,:   ,2   )
    gradient(: ,: ,nz) = gradient(:   ,:   ,nz-1)


  end subroutine get_gradient

  ! =====================================================================

end module utilities

! =====================================================================
! Testdrive Self-gravity module

program test_SOR
  use self_gravity, only : phi, init_SOR, SOR, four_pi_G
  use utilities
  implicit none

  integer, parameter :: nx=128, ny=128, nz=128
  real :: density(nx,ny,nz)  ! assumes 1 guard cell
  real, parameter :: xphys=4.*3.1E18     ! 4pc
  real  :: dx,dy,dz
  real  :: rho1(nx,ny,nz), f_g(nx,ny,nz)

  dx=xphys/float(nx)
  dy=dx
  dz=dx

  !  initializes SOR
  call init_SOR(nx,ny,nz,dx,dy,dz)

  !  read data from file
  call read_data(nx,ny,nz,0,density)

  !density(:,:,:)=1E-18
  !phi(:,:,:)=density(:,:,:)
  !  obtain phi
  call SOR(nx,ny,nz,dx,dy,dz,density)

  !  write rho and phi
  call write_data(nx,ny,nz,0,density(1:nx,1:ny,1:nz),'rho')

  call get_laplacian(nx,ny,nz,dx,dy,dz,phi,rho1)
  rho1=rho1/four_pi_G
  call write_data(nx,ny,nz,0,rho1,'rh1')

  call get_gradient(nx,ny,nz,dx,dy,dz,phi,f_g)
  f_g(:,:,:)=-f_g(:,:,:)
  call write_data(nx,ny,nz,0,f_g,'f_g')

  call write_data(nx,ny,nz,0,phi,'phi')


end program test_SOR
