!=======================================================================
!  MOVING SNR
!=======================================================================
module msnr
  use parameters
  use constants
  implicit none

contains

!  subroutine init_ot()
!    implicit none
!    
!  end subroutine init_ot

  !--------------------------------------------------------------------
  ! Initial conditions for the orzag_tang
  subroutine impose_msnr(u)
    use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
         nxtot, nytot, nztot, neq, nx, ny, nz, &
         rsc, rhosc, vsc, Psc, cv
    use globals,    only : dx, dy, dz, coords
    use constants,  only : pi, amh, kb
    implicit none
    real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

    real :: xpared 
    real :: FETH, FEK, ek, e0
    real :: AMSNR
    real :: xm, n, d0, rho0, rhoc, dc
    real :: radioin, rradio, snrin
    real :: vstar ! Star velocity
    real :: exr 
    real ::  kk, xd
    real :: de, rrc, vmax
    real :: alfab, ts, temp
    real :: rad, x, y, z
    real, parameter :: mu_snr=0.6
    

    integer :: i, j, k


!      xpared=2.542d21 ! POSICION DE LA PARED
      AMSNR=1.4     ! Ejecta mass in solar masses
      E0=1.0        ! TOTAL ENERGY IN UNITS OF 10**51 ERG?
      
      FEK = 0.95    ! Kinetic energy fraction ! no funciona
      FETH=1.-FEK   ! Thermal energy fraction
      ek=fek*e0     ! Kinetic energy
      alfab=50.*pi/180.
      XM=3./7.      ! Mass fraction for outer region
      n=7.          ! Index for density decay in outer regoin
      d0=.05        ! Density of interstellar medium?
      rho0=4.*d0    ! Dont know
      RADIOIN=2.d18 ! Core radius
!      vstar=0.!1.d7    ! Star velocity (cm/s)
      vstar=3.d7 ! For the runaway SNR      
      exr=1./(3.-n)
      RRADIO=RADIOIN/3.1E18 ! Radius normalized to pc 
      SNRIN=RADIOIN*(1-((XM*(3-n)*AMSNR*3.2183)/(mu_snr* &
           rho0*RRADIO**3.)))**exr
      TS=2.864E9*FETH*E0*mu_snr/(CV*AMSNR)
      xd=100.*pc
      KK =1.e-6*xd**3
    
    ! FOR THE CORE
    
    rrc=snrin/3.1d18
    rhoc=6.7134e-23*3.*(1.-XM)*AMSNR/(4.*pi*rrc**3.)
    dc=rhoc/(mu_snr*amh)
    
    VMAX=((4.49679e9*ek**.5)/(2.*pi*mu_snr)**0.5)/ &
         (dc*rrc**5./(5.*rradio**2.)+&
         rho0*rradio**3.*(1.-(rradio/rrc)**(n-5.))/(5.-n))**0.5
    
    
    do k=nzmin,nzmax
       z=(float(k+coords(2)*nz-nztot/2) - 0.5)*dz
       do j=nymin,nymax
          y=(float(j+coords(1)*ny-nytot/2) - 0.5)*dy
          do i=nxmin,nxmax
             x=(float(i+coords(0)*nx-nxtot/2) - 0.5)*dx 

             !   measured from the centre of the computational mesh
             rad=sqrt(x**2+y**2+z**2)  

             if(rad.le.snrin/rsc) then                 
                ! CORE OF SRN
                u(1,i,j,k) = rhoc/rhosc!de*mh*mu
                u(2,i,j,k)=u(1,i,j,k)*(vmax*x/radioin+vstar)/vsc
                u(3,i,j,k)=u(1,i,j,k)*(vmax*y/radioin)/vsc
                u(4,i,j,k)=u(1,i,j,k)*(vmax*z/radioin)/vsc
                u(5,i,j,k) = cv*(dc*kb*ts)/psc + &
                     0.5*(u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)/u(1,i,j,k) + &
                     0.5*(u(6,i,j,k)**2+u(7,i,j,k)**2+u(8,i,j,k)**2)
                !                  
             elseif(rad.le.RADIOIN/rsc) then
                ! SHELL OF SRN
                u(1,i,j,k) = rhoc*(snrin/rsc/rad)**n/rhosc
                u(2,i,j,k)=u(1,i,j,k)*(vmax*x/radioin+vstar)/vsc
                u(3,i,j,k)=u(1,i,j,k)*(vmax*y/radioin)/vsc
                u(4,i,j,k)=u(1,i,j,k)*(vmax*z/radioin)/vsc
                u(5,i,j,k) = cv*dc*kb*ts/psc + &
                     0.5*(u(2,i,j,k)**2+ u(3,i,j,k)**2 + u(4,i,j,k)**2)/ u(1,i,j,k) + &
                     0.5*(u(6,i,j,k)**2+u(7,i,j,k)**2+u(8,i,j,k)**2) 
             end if





        end do
       end do
    end do

  end subroutine impose_msnr
  !--------------------------------------------------------------------
  
end module msnr
