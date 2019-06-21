module pic_module

  use parameters
  use constants
  implicit none
  integer, parameter   :: N_MP = 512    !< # of macro particles
  real, allocatable    :: posMP0(:,:)   !< Particles positions
  real, allocatable    :: posMP1(:,:)   !< Positions after predictor
  real, allocatable    :: velMP (:,:)   !< Particles velocities
  integer, allocatable :: partOwner(:)  !< Particle Owner (rank)
  integer              :: n_active
contains

  !================================================================
  ! @brief initialization of module
  subroutine init_pic()

    use constants,  only : pi
    use parameters, only : rsc
    use globals,    only : dz, rank
    implicit none
    integer :: ir , ith, i_mp, x, y, z
    real    :: pos(3)

    allocate( posMP0(N_MP,3) )
    allocate( posMP1(N_MP,3) )
    allocate( velMP (N_MP,3) )
    allocate( partOwner(N_MP))

    !  initialize Owners (-1 means no body has claimed the particle)
    partOwner(:) = -1
    n_active     =  0

    i_mp = 1
    !Stationary vortex setup
    do ir=1,8
      do ith=1,64

        pos(1)= 0.5*real(ir)*cos( real(ith)*(2.*pi)/64. ) + 5.
        pos(2)= 0.5*real(ir)*sin( real(ith)*(2.*pi)/64. ) + 5.
        pos(3)= 0.

        if(isInDomain(pos) ) then
          partOwner(i_mp) = rank
          posMP0(i_mp,:)  = pos(:)
          n_active        = n_active + 1
        endif

        i_mp = i_mp + 1

      end do
    end do

  end subroutine init_pic

  !================================================================
  ! @brief In domain function (logical)
  ! @details Determines if the position of a given point lies within the
  !> procesor domain, returns true if it is, false if it isn't
  ! @param real [in] : 3D position with respect to a cornet of the domain
  function isInDomain(pos)

    use parameters, only : nx, ny, nz
    use globals,    only : coords, dx, dy, dz
    implicit none
    logical  isInDomain
    real, intent(in) :: pos(3)
    integer          :: ind(3)

    ! shift to the local processor
    ind(1) = int(pos(1)/dx) - coords(0)*nx
    ind(2) = int(pos(2)/dy) - coords(1)*ny
    ind(3) = int(pos(3)/dz) - coords(2)*nz

    if ( ind(1)<0  .or. ind(2)<0  .or. ind(3)<0 .or. &
         ind(1)>=nx .or. ind(2)>=ny .or. ind(3)>=nz ) then

      isInDomain = .false.

    else

      isInDomain = .true.

    end if

  end function isInDomain

  !================================================================
  subroutine predictor

    use globals, only : primit, dt_CFL, rank
    implicit none
    integer :: i_mp, i, j, k, l, ind(3)
    real    :: weights(8)

    do i_mp=1, N_MP

      ! exetute only if particle i_mp is in processor domain
      if ( isInDomain( posMP0(i_mp,:) ) ) then

        ! Calculate interpolation reference and weights
        call interpBD(posMP0(i_mp,:),ind,weights)
        !print*,rank,posMP0(i_mp,:)
        !print*, rank,ind(:)
        !  Interpolate the velocity field to particle position
        l=1
        velMP(i_mp,:) = 0
        do k= ind(3),ind(3)+1
          do j=ind(2),ind(2)+1
            do i=ind(1),ind(1)+1
              velMP(i_mp,1) = velMP(i_mp,1) + primit(2,i,j,k)*weights(l)
              velMP(i_mp,2) = velMP(i_mp,2) + primit(3,i,j,k)*weights(l)
              velMP(i_mp,3) = velMP(i_mp,3) + primit(4,i,j,k)*weights(l)
              l = l + 1
            end do
          end do
        end do

        !  predictor step
        posMP1(i_mp,:) = posMP0(i_mp,:)+ dt_CFL*velMP(i_mp,:)

      end if

    end do

  end subroutine predictor

  !================================================================
  subroutine corrector

    use globals, only : primit, dt_CFL
    implicit none
    integer :: i_mp, i, j, k, l, ind(3)
    real    :: weights(8)


    do i_mp=1, N_MP

      ! exetute only if particle i_mp is in processor domain
      if ( isInDomain( posMP1(i_mp,:) ) ) then

        ! Calculate interpolation reference and weights
        call interpBD(posMP1(i_mp,:),ind,weights)
        !  Interpolate the velocity field to particle position, and add
        !  to the velocity from the corrector step
        l=1
        !velMP(i_mp,:) = 0
        do k= ind(3),ind(3)+1
          do j=ind(2),ind(2)+1
            do i=ind(1),ind(1)+1
              velMP(i_mp,1) = velMP(i_mp,1) + primit(2,i,j,k)*weights(l)
              velMP(i_mp,2) = velMP(i_mp,2) + primit(3,i,j,k)*weights(l)
              velMP(i_mp,3) = velMP(i_mp,3) + primit(4,i,j,k)*weights(l)
              l = l + 1
            end do
          end do
        end do

        !  predictor step
        posMP0(i_mp,:) = posMP0(i_mp,:)+ 0.5*dt_CFL*velMP(i_mp,:)

      end if

    end do

  end subroutine corrector

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
  write(unitout) n_active

  !  loop over particles owned by processor and write the active ones
  do i_mp=1,N_MP
    if (partOwner(i_mp)==rank) then

      write(unitout) posMP0(i_mp,:)
      write(unitout) velMP (i_mp,:)

    endif
  end do

end subroutine write_pic

end module pic_module

!================================================================
