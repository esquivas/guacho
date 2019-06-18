module pic_module

  use parameters
  use constants
  implicit none
  integer, parameter :: N_MP = 100  !< # of macro particles
  real               :: posMP0(N_MP, 3) !< Particles positions
  real               :: PosMP1(N_MP, 3) !< Positions after predictor
  real               :: velMP (N_MP, 3) !< Particles velocities

contains

  !================================================================
  ! Se inicializa el modulo de particulas
  subroutine init_nbody()

    implicit none

  end subroutine init_nbody

  !================================================================
  subroutine predictor

    use globals, only : primit, dt_CFL
    use parameters, only : nx, ny, nz
    implicit none
    integer :: i_mp
    real    :: vp(3)   !< particle velocity (grid interpolated)

    do i=1, i_mp

      call interpBD()

    end do

  end subroutine predictor

  !================================================================
  subroutine corrector

    use globals, only : primit, dt_CFL
    use parameters, only : nx, ny, nz
    implicit none

  end subroutine corrector

  !================================================================
  !  @brief bilineal interpolator
  !  @ details: returns weights and indices corresponding to a bilineal
  !> interpolation at for the particle position. If the particle  lies outside
  !> the processor domain it returns a negative flag.
  !> @param real [in] pos(3) : Three dimensional position of particle,
  !> @param integer   ind(3) : Reference corner i0,j0,k0 indices of the cube of
  !> cells that are used in the interpolation
  !> @param real [out] weights(8) : weights associated with each corner of the
  !> interopolation in the following order
  !> 1-> i0,j0,k0,  3-> i0,j1,k0,  5-> i0,j0,k1,  7-> i0,j1,k1
  !> 2-> i1,j0,k0,  4-> i1,j1,k0,  6-> i1,j0,k1,  8-> i1,j1,k1
  !> @param logical [out] in_domain  : True if particle is inside processor domain
  !> False if outside
  subroutine interpBD(pos,ind,weights,in_domain)

    use parameters, only : nxtot, nytot, nztot, nx, ny, nz, dx, dy, dz
    use globals,    only : coords
    implicit none
    real,    intent(in)  :: pos(3)
    integer, intent(out) :: ind(3)
    logical, intent(out) :: in_domain
    real,    intent(out) :: weights(8)
    real                 :: x0, y0, z0

    !ind(1) = int(  pos(1) - real( coords(1)*nx*dx ) + 0.5  )
    !ind(2) = int(  pos(2) - real( coords(2)*ny*dy ) + 0.5  )
    !ind(3) = int(  pos(3) - real( coords(3)*nz*dz ) + 0.5  )

    ! get the index in the whole domain (integer part)
    x0 = pos(1)/nxtot
    y0 = pos(2)/nytot
    z0 = pos(3)/nztot

    ! shift to the local processor
    ind(1) = int(x0) - coords(0)*nx
    ind(2) = int(y0) - coords(1)*ny
    ind(3) = int(z0) - coords(2)*nz

    if ( ind(1)<0  .or. ind(2)<0  .or. ind(3)<0 .or. &
         ind(1)>nx .or. ind(2)>ny .or. ind(3)>nz ) then

      in_domain = .false.

      return

    else
      ! get the reminder
      x0 = x0 - int(x0)
      y0 = y0 - int(y0)
      z0 = z0 - int(z0)

      in_domain = .true.

      weights(1) = (1-x0)*(1-y0)*(1-z0)
      weights(2) =   x0  *(1-y0)*(1-z0)
      weights(3) = (1-x0)*  y0  *(1-z0)
      weights(4) =   x0  *  y0  *(1-z0)
      weights(5) = (1-x0)*(1-y0)*  z0
      weights(6) =   x0  *(1-y0)*  z0
      weights(7) = (1-x0)*  y0  *  z0
      weights(8) =   x0  *  y0  *  z0

      return

    endif

  end subroutine interpBD

  !================================================================


end module pic_module

!================================================================
