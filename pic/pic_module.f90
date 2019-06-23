module pic_module
  use parameters, only : N_MP
  use globals   , only : Q_MP0, Q_MP1, partOwner, n_activeMP
  implicit none

contains

  !================================================================
  ! @brief initialization of module
  subroutine init_pic()

    implicit none

    allocate( Q_MP0(N_MP,6) )
    allocate( Q_MP1(N_MP,3) )
    allocate( partOwner(N_MP))

  end subroutine init_pic

  !================================================================
  subroutine PICpredictor

    use globals,   only : primit, dt_CFL
    use utilities, only : isInDomain
    implicit none
    integer :: i_mp, i, j, k, l, ind(3)
    real    :: weights(8)

    do i_mp=1, N_MP

      ! exetute only if particle i_mp is in processor domain
      if ( isInDomain( Q_MP0(i_mp,1:3) ) ) then

        ! Calculate interpolation reference and weights
        call interpBD(Q_MP0(i_mp,1:3),ind,weights)

        !  Interpolate the velocity field to particle position
        l=1
        Q_MP0(i_mp,4:6) = 0
        do k= ind(3),ind(3)+1
          do j=ind(2),ind(2)+1
            do i=ind(1),ind(1)+1
              Q_MP0(i_mp,4) = Q_MP0(i_mp,4) + primit(2,i,j,k)*weights(l)
              Q_MP0(i_mp,5) = Q_MP0(i_mp,5) + primit(3,i,j,k)*weights(l)
              Q_MP0(i_mp,6) = Q_MP0(i_mp,6) + primit(4,i,j,k)*weights(l)
              l = l + 1
            end do
          end do
        end do

        !  predictor step
        Q_MP1(i_mp,1:3) = Q_MP0(i_mp,1:3)+ dt_CFL*Q_MP0(i_mp,4:6)

      end if

    end do

  end subroutine PICpredictor

  !================================================================
  subroutine PICcorrector

    use globals,   only : primit, dt_CFL
    use utilities, only : isInDomain
    implicit none
    integer :: i_mp, i, j, k, l, ind(3)
    real    :: weights(8), vel1(3)


    do i_mp=1, N_MP

      ! exetute only if particle i_mp is in processor domain
      if ( isInDomain( Q_MP1(i_mp,1:3) ) ) then

        ! Calculate interpolation reference and weights
        call interpBD(Q_MP1(i_mp,1:3),ind,weights)
        !  Interpolate the velocity field to particle position, and add
        !  to the velocity from the corrector step
        l=1
        vel1(:) = 0
        do k= ind(3),ind(3)+1
          do j=ind(2),ind(2)+1
            do i=ind(1),ind(1)+1
              vel1(1) = vel1(1) + primit(2,i,j,k)*weights(l)
              vel1(2) = vel1(2) + primit(3,i,j,k)*weights(l)
              vel1(3) = vel1(3) + primit(4,i,j,k)*weights(l)
              l = l + 1
            end do
          end do
        end do

        !  corrector step
        Q_MP0(i_mp,1:3) = Q_MP0(i_mp,1:3) &
                        + 0.5*dt_CFL*( Q_MP0(i_mp,4:6) + vel1(1:3) )

      end if

    end do

  end subroutine PICcorrector

  !================================================================
  !  @brief bilineal interpolator
  !  @ details: returns weights and indices corresponding to a bilineal
  !> interpolation at for the particle position.
  !> @param real [in] pos(3) : Three dimensional position of particle,
  !> @param integer   ind(3) : Reference corner i0,j0,k0 indices of the cube of
  !> cells that are used in the interpolation
  !> @param real [out] weights(8) : weights associated with each corner of the
  !> interopolation in the following order
  !> 1-> i0,j0,k0,  3-> i0,j1,k0,  5-> i0,j0,k1,  7-> i0,j1,k1
  !> 2-> i1,j0,k0,  4-> i1,j1,k0,  6-> i1,j0,k1,  8-> i1,j1,k1
    subroutine interpBD(pos,ind,weights)

    use parameters, only : nx, ny, nz
    use globals,    only : coords, dx, dy, dz
    implicit none
    real,    intent(in)  :: pos(3)
    integer, intent(out) :: ind(3)
    real,    intent(out) :: weights(8)
    real                 :: x0, y0, z0

    ! get the index in the whole domain (integer part)
    x0 = pos(1)/dx
    y0 = pos(2)/dy
    z0 = pos(3)/dz

    ! shift to the local processor
    ind(1) = int(x0) - coords(0)*nx
    ind(2) = int(y0) - coords(1)*ny
    ind(3) = int(z0) - coords(2)*nz

    ! get the reminder
    x0 = x0 - int(x0)
    y0 = y0 - int(y0)
    z0 = z0 - int(z0)

    weights(1) = (1-x0)*(1-y0)*(1-z0)
    weights(2) =   x0  *(1-y0)*(1-z0)
    weights(3) = (1-x0)*  y0  *(1-z0)
    weights(4) =   x0  *  y0  *(1-z0)
    weights(5) = (1-x0)*(1-y0)*  z0
    weights(6) =   x0  *(1-y0)*  z0
    weights(7) = (1-x0)*  y0  *  z0
    weights(8) =   x0  *  y0  *  z0

    return

  end subroutine interpBD

  !================================================================
  ! @brief Writes pic output
  ! @details Writes position and particles velocities in a single file
  !> format is binary, one integer with the number of points in domain, and then
  !> all the positions and velocities ( in the default real precision)
  subroutine write_pic(itprint)

    use parameters, only : outputpath
    use globals,    only : rank
    implicit none
    integer, intent(in) :: itprint
    character(len = 128) :: fileout
    integer              :: unitout, i_mp

#ifdef MPIP
    write(fileout,'(a,i3.3,a,i3.3,a)')  &
        trim(outputpath)//'BIN/pic',rank,'.',itprint,'.bin'
    unitout=rank+10
#else
    write(fileout,'(a,i3.3,a)')  trim(outputpath)//'BIN/pic',itprint,'.bin'
    unitout=10
#endif

  open(unit=unitout,file=fileout,status='unknown',access='stream')

  !  write how many particles are in rank domain
  write(unitout) n_activeMP

  !  loop over particles owned by processor and write the active ones
  do i_mp=1,N_MP
    if (partOwner(i_mp)==rank) then

      write(unitout) i_mp
      write(unitout) Q_MP0(i_mp,1:6)

    endif
  end do

end subroutine write_pic

end module pic_module

!================================================================
