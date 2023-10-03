!=======================================================================
!> @file stokes_lp.f90
!> @brief stokes parameters utilities
!> @author M. Schneiter, Alejandro Esquivel
!> @date 4/Oct/2019

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

!> @brief stokes lp utilities
!> @details Utilities to compute the Stokes parameters with the module of
!> Lagrangian particles

module stokes_lp_utilities

  implicit none
  real, allocatable :: stokesTab(:,:)
  integer           :: nTabLines

contains

!> @brief Initializes data
!> @details Initializes data, MPI and other stuff
subroutine init_stokes()

  !  Initializes MPI, data arrays, etc
  use parameters
  use globals, only : u, dx, dy, dz, coords, rank, left, right,                &
                      top, bottom, out, in, rank, comm3d,                      &
                      Q_MP0, MP_SED, P_DSA, partID, partOwner
  implicit none
  integer :: nps, err, i, inp
  integer, dimension(0:ndim-1) :: dims
  logical, dimension(0:ndim-1) :: period
  real :: xTab, fTab, gTab ! Tabulation of Modified Bessel Functions

  !initializes MPI
#ifdef MPIP
#ifdef PERIODX
  logical, parameter :: perx=.true.
#else
  logical, parameter :: perx=.false.
#endif
#ifdef PERIODY
  logical, parameter :: pery=.true.
#else
  logical, parameter :: pery=.false.
#endif
#ifdef PERIODZ
  logical, parameter :: perz=.true.
#else
  logical, parameter :: perz=.false.
#endif
  period(0)=perx
  period(1)=pery
  period(2)=perz
  dims(0)  =MPI_NBX
  dims(1)  =MPI_NBY
  dims(2)  =MPI_NBZ

  call mpi_init (err)
  call mpi_comm_rank (mpi_comm_world,rank,err)
  call mpi_comm_size (mpi_comm_world,nps,err)
  if (nps.ne.np) then
     print*,                                                                   &
     'processor number (',nps,') is not equal to pre-defined number (',np,')'
     call mpi_finalize(err)
     stop
  endif
#else
  rank=0
  coords(:)=0
#endif
  if(rank.eq.master) then
  print '(a)' ,"*******************************************"
     print '(a)' ,"                        _                 *"
     print '(a)' ,"  __   _   _  __ _  ___| |__   ___    3   *"
     print '(a)' ," / _ `| | | |/ _` |/ __| '_ \ / _ \    D  *"
     print '(a)' ,"| (_| | |_| | (_| | (__| | | | (_) |      *"
     print '(a)' ," \__, |\__,_|\__,_|\___|_| |_|\___/       *"
     print '(a)' ," |___/                                    *"
     print '(a)' ,"                                          *"
  endif
#ifdef MPIP
  if(rank.eq.master) then
     print '(a,i3,a)','*    running with mpi in', np , ' processors    *'
     print '(a)' ,'*******************************************'
     print '(a)', 'Calculating Stokes Parameters'
  end if
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, .true., comm3d, err)
  call mpi_comm_rank(comm3d, rank, err)
  call mpi_cart_coords(comm3d, rank, ndim, coords, err)
  print '(a,i3,a,3i4)', 'processor ', rank,                                    &
        ' ready w/coords',coords(0),coords(1),coords(2)
  call mpi_cart_shift(comm3d, 0, 1, left  , right, err)
  call mpi_cart_shift(comm3d, 1, 1, bottom, top  , err)
  call mpi_cart_shift(comm3d, 2, 1, out   , in   , err)
  call mpi_barrier(mpi_comm_world, err)
  !
#else
  print '(a)' ,'*******************************************'
  print '(a)' ,'*     running on a single processor       *'
  print '(a)' ,'*******************************************'
  print '(a)', 'Calculating Stokes Parameters'
#endif

  !   grid spacing
  dx=xmax/nxtot
  dy=ymax/nytot
  dz=zmax/nztot

  !   allocate big arrays in memory

  !  MHD
  allocate( u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )

  !  LP
  if (lmp_distf) then
    allocate( Q_MP0(N_MP,12) )
    !Q_MP0(i, eq) has the following info:
    ! eq = 1-3 : x, y, z
    ! eq = 4-6 : vx, vy, vz
    ! eq = 7   : b**2/2 = (bx**2+by**2+bz**2)/2
    ! eq = 8   : rho
    ! eq = 9   : P
    ! eq = 10  : shock flag (1 if shocked)
    ! eq = 11  : compression ratio (does not reset)
    ! eq = 12  : angle between the shock normal and the preshock field
    allocate( MP_SED(2,NBinsSEDMP,N_MP) )
    ! MP_SED(1,:,i) :  Energy (Lagrangian) bins
    ! MP_SED(1,:,i) :  Number of MP with Energy E_i +- Delta E
    allocate( P_DSA(N_MP,2,8))
    ! P_DSA(i, 1, :) : Pre  shock MHD info (U1 in Vaidya et al 2018)
    ! P_DSA(i, 2, :) : Post shock MHD info (U2 in Vaidya et al 2018)
    allocate( partID   (N_MP) )  ! Individual particle identifier
    allocate( partOwner(N_MP) )  ! Rank of the processor that owns said particle

    do inp=0,np-1  ! take turns to read tables
      if (rank==inp) then

        open(unit=10,file= trim(workdir)//'../src/LPlib/SynchroBessels.tab',   &
             status='old')
        read(10,*) nTabLines
        !print*, nlines
        allocate(stokesTab(3,nTabLines))

        do i=1,nTabLines
          read(10,*) xTab, fTab, gTab
          stokesTab(1,i)=xTab
          stokesTab(2,i)=fTab
          stokesTab(3,i)=gTab
        end do

        close(unit=10)
        print*, 'rank: ',rank,'  nTablines: ',nTabLines
      end if
#ifdef MPIP
      call mpi_barrier(mpi_comm_world, err)
#endif

    end do

  else
    print '(a)', "The SED of the Lagraingian particles is not enabled"
    call mpi_finalize(err)
    stop
  end if

end subroutine init_stokes

!=======================================================================
!> @brief reads data from file
!> @details reads data from file
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param integer [in] itprint : number of output
!> @param string [in] filepath : path where the output is
subroutine read_data(u,itprint,filepath)

  use parameters, only : np, neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax,    &
                         lmp_distf
  use globals, only : rank, comm3d, Q_MP0, MP_SED, P_DSA, partID, partOwner,   &
                      n_activeMP
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: itprint
  character (len=128), intent(in) :: filepath
  character                       :: byte_read
  integer :: nxp, nyp, nzp, x0p, y0p, z0p, &
             mpi_xp, mpi_yp, mpi_zp,neqp, neqdynp, nghostp
  integer :: unitin, ip, err, unitin2, npP, n_mpP, i_activeP, NBinsSEDMPP
  real    :: dxp, dyp, dzp, scal(3), cvp
  integer :: i_mp
  character (len=128) file1, file2

  take_turns : do ip=0,np-1
    if (rank == ip) then

#ifdef MPIP
        write(file1,'(a,i3.3,a,i3.3,a)')                                       &
             trim(filepath)//'BIN/points',rank,'.',itprint,'.bin'
        write(file2,'(a,i3.3,a,i3.3,a)')                                       &
              trim(filepath)//'BIN/lmp',rank,'.',itprint,'.bin'
        unitin =rank+10
        unitin2=rank+11
#else
         write(file1,'(a,i3.3,a)')                                             &
              trim(filepath)//'BIN/points',itprint,'.bin'
         write(file1,'(a,i3.3,a)')
              trim(filepath)//'BIN/lmp',itprint,'.bin'
         unitin =10
         unitin2=11
#endif
         open(unit=unitin,file=file1,status='unknown', access='stream')
         !   discard the ascii header
         do while (byte_read /= achar(255) )
           read(unitin) byte_read
           !print*, byte_read
         end do
         !  read bin header, sanity check to do
         read(unitin) byte_read
         read(unitin) byte_read
         read(unitin) nxp, nyp, nzp
         read(unitin) dxp, dyp, dzp
         read(unitin) x0p, y0p, z0p
         read(unitin) mpi_xp, mpi_yp, mpi_zp
         read(unitin) neqp, neqdynp
         read(unitin) nghostp
         read(unitin) scal(1:3)
         read(unitin) cvp
         read(unitin) u(:,:,:,:)
         close(unitin)

         print'(i3,a,a)',rank,' read file:',trim(file1)

         open(unit=unitin2,file=file2,status='unknown', access='stream')
         read(unitin2) npP, N_MPP, i_activeP, NBinsSEDMPP
         n_activeMP = i_activeP
         do i_mp=1,i_activeP
           read(unitin2) partID(i_mp)
           read(unitin2) Q_MP0(i_mp,1:3)
           read(unitin2) Q_MP0(i_mp,11:12)
           read(unitin2) MP_SED(1,:,i_mp)
           read(unitin2) MP_SED(2,:,i_mp)
           read(unitin2) P_DSA(i_mp,:,:)
           partOwner(i_mp) = rank
         end do

         close(unitin2)
         print'(i3,a,a,a,i0,a)',rank,' read file:',trim(file2),' (got ',       &
                                n_activeMP,' active LPs)'

    end if
    call mpi_barrier(comm3d,err)
  end do take_turns

end subroutine read_data

!=======================================================================
!> @brief gets position of a cell
!> @details Returns the position and spherical radius calculated with
!! respect to  the center of the grid
!> @param integer [in] i : cell index in the x direction
!> @param integer [in] j : cell index in the y direction
!> @param integer [in] k : cell index in the z direction
!> @param real    [in] x : x position in the grid
!> @param real    [in] y : y position in the grid
!> @param real    [in] z : z position in the grid
  subroutine getXYZ(i,j,k,x,y,z)
    use globals,    only : dx, dy, dz, coords
    use parameters, only : nx, ny, nz, nxtot, nytot, nztot
    implicit none
    integer, intent(in)  :: i, j, k
    real,    intent(out) :: x, y, z

    x=(real(i+coords(0)*nx-nxtot/2) - 0.5)*dx
    y=(real(j+coords(1)*ny-nytot/2) - 0.5)*dy
    z=(real(k+coords(2)*nz-nztot/2) - 0.5)*dz

  end subroutine getXYZ

!=======================================================================
!> @brief Rotation around the X axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid
subroutine rotation_x(theta,x,y,z,xn,yn,zn)
   implicit none
   real, intent(in ) :: theta, x, y, z
   reaL, intent(out) :: xn, yn, zn
   xn =   x
   yn =   y*cos(theta) - z*sin(theta)
   zn =   y*sin(theta) + z*cos(theta)
 end subroutine rotation_x

!=======================================================================
!> @brief Rotation around the Y axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid
 subroutine rotation_y(theta,x,y,z,xn,yn,zn)
   implicit none
   real, intent(in ) :: theta, x, y, z
   real, intent(out) :: xn, yn, zn
   xn =   x*cos(theta) + z*sin(theta)
   yn =   y
   zn = - x*sin(theta) + z*cos(theta)
 end subroutine rotation_y

!=======================================================================

!> @brief Rotation around the Z axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid
 subroutine rotation_z(theta,x,y,z,xn,yn,zn)

   implicit none
   real, intent(in ) :: theta, x, y, z
   real, intent(out) :: xn, yn, zn
   xn =   x*cos(theta) - y*sin(theta)
   yn =   x*sin(theta) + y*cos(theta)
   zn =   z
 end subroutine rotation_z

!=======================================================================
!> @brief Fill target map
!> @details Fills the target map of one MPI block
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param real [in] u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax) :
!! conserved variables
!> @param real [out] map(nxmap,mymap) : Target map
!> @param real [in] dxT : target pixel width
!> @param real [in] dyT : target pixel height
!> @param real [in] thetax : Rotation around X
!> @param real [in] thetay : Rotation around Y
!> @param real [in] thetaz : Rotation around Z
subroutine fill_map(nxmap, nymap, nmaps, map, freq_obs,dxT , dyT,                       &
                   theta_x, theta_y, theta_z)
  use globals,    only : u, Q_MP0, n_activeMP, dz
  use parameters, only : xmax, ymax, zmax, Bsc, rsc
  use lmp_module, only : interpBD
  implicit none
  integer, intent(in)  :: nxmap,nymap,nmaps
  real,    intent(in)  :: freq_obs,dxT, dYT, theta_x, theta_y, theta_z
  real,    intent(out) :: map(nxmap,nymap,nmaps)
  integer              :: i_mp, iobs, jobs, ind(3), i, j, k, l
  real                 :: x, xn, y, yn, z, zn
  real                 :: weights(8), Bx, By, Bz, Bxn, Byn, Bzn, SI, SQ, SU

  !  Clear target map
  map(:,:,:) = 0.0

  !  loop over al particles
  !  (there's no need to test if it's active, we just read active LPs)
  do i_mp=1, n_activeMP

    !  unpack the positions (just for clarity) and recenter
    x = Q_MP0(i_mp,1) - xmax/2.0
    y = Q_MP0(i_mp,2) - ymax/2.0
    z = Q_MP0(i_mp,3) - zmax/2.0

    !  Interpolate Bfield to each particle position (and scale it to Gauss)
    call interpBD(Q_MP0(i_mp,1:3),ind,weights)
    Bx = 0.0
    By = 0.0
    Bz = 0.0
    l = 1
    do k= ind(3),ind(3)+1
      do j=ind(2),ind(2)+1
        do i=ind(1),ind(1)+1
          Bx = Bx + u(6,i,j,k) * weights(l) * Bsc  !  cgs
          By = By + u(7,i,j,k) * weights(l) * Bsc  !  cgs
          Bz = Bz + u(8,i,j,k) * weights(l) * Bsc  !  cgs
          l  = l + 1
        end do
      end do
    end do

    !  rotate the coordinates
    call rotation_x(theta_x, x , y , z,  xn, yn, zn)
    call rotation_y(theta_y, xn, yn, zn, x,  y,  z )
    call rotation_z(theta_z, x , y,  z,  xn, yn, zn)
    !  rotate B
    call rotation_x(theta_x, Bx , By , Bz,  Bxn, Byn, Bzn)
    call rotation_y(theta_y, Bxn, Byn, Bzn, Bx,  By,  Bz )
    call rotation_z(theta_z, Bx , By,  Bz,  Bxn, Byn, Bzn)

    ! This is the position on the target (centered)
    ! Integration is along Z
    iobs=nint(xn/dxT) + nxmap/2
    jobs=nint(yn/dyT) + nymap/2

    !  Fill the maps
    if( (iobs >=1    ).and.(jobs >=1    ) .and.                                &
        (iobs <=nxmap).and.(jobs <=nymap) ) then

      !  obtain stokes parameters of a single LP in one cells
      !  the integrals in eqs 37 and 41 of Vaidya et al. is achieved by
      !  summing all the elements in a map.
      if (Q_MP0(i_mp, 11) /= 0) then

        !call get_stokes(i_mp,freq_obs,Bx,By,SI,SQ,SU)
        call get_stokes(i_mp,freq_obs,Bxn,Byn,SI,SQ,SU)

        map(iobs,jobs,1)= map(iobs,jobs,1) + SI*dz*rsc
        map(iobs,jobs,2)= map(iobs,jobs,2) + SQ*dz*rsc
        map(iobs,jobs,3)= map(iobs,jobs,3) + SU*dz*rsc

      end if

    end if

  end do

end subroutine fill_map

!=======================================================================
subroutine get_stokes(i_mp,freq_obs,Bx,By,I,Q,U)
  use parameters, only : NBinsSEDMP
  use globals,    only : MP_SED, rank
  implicit none
  integer, intent(in)  :: i_mp
  real,    intent(in)  :: freq_obs, Bx, By
  real,    intent(out) :: I, Q, U
  real, parameter :: Jconst = 1.8755e-23 ! sqrt(3)*e^3/(4pi me c^2)
  real, parameter :: xconst = 1.5754e-19 ! 4pi me^3 c^5 /(3 e)
  real            :: Bperp, x0, x1, Fsyn0, Fsyn1, Fpol0, Fpol1, Isyn, Ipol,    &
                     x, F, G, slopeF, slopeG
  integer         :: ibin
  real            :: Jpol

  Bperp = sqrt(Bx**2+By**2)

  Isyn  = 0.0
  Ipol  = 0.0

  do ibin = 1, NBinsSEDMP-1

    x0    = MP_SED(1,ibin  ,i_mp)
    x1    = MP_SED(1,ibin+1,i_mp)

    !x = xconst*freq_obs / (x0*x1*Bperp)
    x = xconst*freq_obs / (MP_SED(1,ibin  ,i_mp)**2*Bperp)

    if(isnan(x)) then

      print*, 'Invalid SED for particle', i_mp, ' ignoring it'
      I = 0.0
      Q = 0.0
      U = 0.0
      return

    endif

    call getBessels(x, F, G)
    Fsyn0 = MP_SED(2,ibin  ,i_mp)*F
    Fpol0 = MP_SED(2,ibin  ,i_mp)*G

    x = xconst*freq_obs / (MP_SED(1,ibin+1,i_mp)**2*Bperp)
    call getBessels(x, F, G)
    Fsyn1 = MP_SED(2,ibin+1,i_mp)*F
    Fpol1 = MP_SED(2,ibin+1,i_mp)*G

    if (x1 /= x0) then
        slopeF = ( log(Fsyn1/Fsyn0) )/( log(x1/x0) )
        slopeG = ( log(Fpol1/Fpol0) )/( log(x1/x0) )
      else
        Isyn = -1.
        Ipol = -1.
        stop
        return
      end if

      if (slopeF /= -1.) then
        Isyn = Isyn + Fsyn0/(slopeF+1.)*(x1*(x1/x0)**slopeF-x0)
      else
        Isyn = Isyn + Fsyn0 * x0 * log(x1/x0)
      end if

      if (slopeG /= -1.) then
        Ipol = Ipol + Fpol0/(slopeG+1.)*(x1*(x1/x0)**slopeG-x0)
      else
        Ipol = Ipol + Fpol0 * x0 * log(x1/x0)
      end if

  end do

  I    = Jconst*Bperp*Isyn      ! eq (48) Vaidya + 2018
  Jpol = Jconst*Bperp*Ipol      ! eq (41) ''

  ! Eqs (49-50) ''
  Q = Jpol * (Bx**2-By**2 ) / Bperp**2
  U = Jpol * ( -2.0*Bx*By ) / Bperp**2

  !I = 0.0
  !Q = 0.0
  !U = 0.0

end subroutine get_stokes

!=======================================================================
subroutine getBessels(x, F, G)
  implicit none
  real, intent(in ) :: x
  real, intent(out) :: F, G
  real              ::  x1, x2, f1, f2, g1, g2
  integer           :: i

  !  if x is smaller than the first element in the tables use assymptotic
  !  expression (see Mathematica notebook)
  if (x <= stokesTab(1,1) ) then
    !print*,'SMALL',x
    F = 2.14953*x**(1.0/3.0)
    G = 1.07476*x**(1.0/3.0)
    return
  else if (x >= stokesTab(1,nTabLines) ) then
    if (x < 30) then   ! needed to avoid underflows
      !print*,'LARGE',x
      F = 0.278823*x*exp(-4.0*x/3.0)                                           &
      + exp(-x)*(1.1147 + 0.938696*x + 0.092941*x**(1.5))/sqrt(x)
      G = 1.25331*exp(-x)/x**(1.5) * (-0.35108 + x*(0.0972222 + x) )
      return
    else
      !print*,'extra LARGE',x
      F = 7.61e-13
      G = 6.44e-13
      return
    end if
  end if

  !  Otherwise use the tables and interpolate them

  i = int( 1 + (nTabLines-1) *                                                 &
      log10(x/stokesTab(1,1))/log10(stokesTab(1,nTabLines)/stokesTab(1,1)) )
  !i = max(i,    1      )
  !i = min(i,nTabLines-1)
  !if (i < 1) print*, i, x, stokesTab(1,1), stokesTab(1,nTabLines)
  if (i==nTabLines-1) then
    F = stokesTab(2,i+1)
    G = stokesTab(3,i+1)
    !if (rank == 0) print'(i0,2es15.3)',i, x, f
    return
  else if (i==1) then
    F = stokesTab(2,i)
    G = stokesTab(3,i)
    !if (rank == 0) print'(i0,2es15.3)',i, x, f
    return
  else
    x1 = stokesTab(1,i  )
    x2 = stokesTab(1,i+1)
    f1 = stokesTab(2,i  )
    f2 = stokesTab(2,i+1)
    g1 = stokesTab(3,i  )
    g2 = stokesTab(3,i+1)

    F = f1 * (x/x1)**(log10(f2/f1)/log10(x2/x1))
    G = g1 * (x/x1)**(log10(g2/g1)/log10(x2/x1))
    !if (rank == 0) print'(i0,6es15.3)',i, x, x1,x2,f1,f2,f

  end if

end subroutine getBessels
!=======================================================================
!> @brief Writes projection to file
!> @details Writes projection to file
!> @param integer [in] itprint : number of output
!> @param string [in] filepath : path where to write
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param integer [in] nvmap : Number of maps (I,Q,U)
!> @param real [in] map(nxmap,mymap) : Target map
subroutine  write_stokes(itprint,filepath,nxmap,nymap,nmaps,map)
  implicit none
  integer, intent(in) :: nxmap, nymap,nmaps,itprint
  character (len=128), intent(in) :: filepath
  real, intent(in) :: map(nxmap,nymap,nmaps)
  character (len=128) file1
  integer ::  unitout

  write(file1,'(a,i3.3,a)')  trim(filepath)//'BIN/stokes-',itprint,'.bin'
  unitout=11
  open(unit=unitout,file=file1,status='unknown',access='stream')

  write (unitout) map(:,:,:)
  close(unitout)

  print'(a,a)'," wrote file:",trim(file1)

end subroutine write_stokes

!=======================================================================

end module stokes_lp_utilities

!=======================================================================
!> @brief Computes the Ly-alpha apbsorption
!> @details Computes the Ly-alpha apbsorption
!! @n It rotates the data along each of the coordinates axis
!! by an amount @f$ \theta_x, \theta_y, \theta_z @f$, and the LOS
!! is along the Z axis
program stokes_lp
  use constants, only : pi
  use parameters, only : xmax,master, mpi_real_kind, outputpath, nxtot, nytot
  use globals, only : u, rank, comm3d
  use stokes_lp_utilities
#ifdef MPIP
  use mpi
#endif
  implicit none

  character (len=128) :: filepath
  integer :: err
  integer :: itprint
  !
  real, parameter :: theta_x = 0.0 *pi/180.
  real, parameter :: theta_y = 0.0 *pi/180.
  real, parameter :: theta_z = 0.0 *pi/180.
  !   map and its dimensions
  integer, parameter :: nmaps= 3   !< (1=I, 2=Q, 3=U)
  integer            :: nxmap, nymap
  real :: dxT, dyT, freq_obs
  real, allocatable :: map(:,:,:), map1(:,:,:)
  !real :: map(nxmap, nymap,nvmap), map1(nxmap, nymap,nvmap)

  nxmap = nxtot
  nymap = nytot
  allocate( map(nxmap, nymap,nmaps))
  allocate(map1(nxmap, nymap,nmaps))

  ! initializes program
  call init_stokes()

  !  Target pixel size, relative to the simulation
  dxT= xmax/real(nxmap)
  dyT= dxT

  ! chose output (fix later to input form screen)
  filepath=trim(outputpath) !'/datos/esquivel/EXO-GUACHO/P1c/'

  freq_obs  =  1.40e9 !< frequency of observation (Hz)

  loop_over_outputs : do itprint=0,15

    !  read MHD and particles data
    call read_data(u,itprint,filepath)

    !  resets map
    map(:,:,:)  = 0.0
    map1(:,:,:) = 0.0

    if (rank == master) then
       print'(a)', 'Calculating projection with angles of rotaton'
       print'(f6.2,a,f6.2,a,f6.2,a)',theta_x*180./pi,'° around X, '            &
                                    ,theta_y*180./pi,'° around Y, '            &
                                    ,theta_z*180./pi,'° around Z, '
       print'(a,es10.3,a)', 'Stokes parameters for a frequency of ',freq_obs,' Hz'
    end if

    !  add info to the map
    call fill_map(nxmap,nymap,nmaps,map,freq_obs,dxT,dyT,theta_x, theta_y, theta_z)
    !  sum all the partial sums

    call mpi_reduce(map,map1,nxmap*nymap*nmaps, mpi_real_kind, mpi_sum, master,&
                    comm3d, err)

    !  write result
    if (rank == master) then
      call write_stokes(itprint,filepath,nxmap,nymap,nmaps,map1)
    end if

  end do loop_over_outputs

  if (rank == master) print*, 'my work here is done, have a  nice day'
#ifdef MPIP
  call mpi_finalize(err)
#endif

  stop

end program stokes_lp


!=======================================================================
