module nbody

  use parameters
  use constants
  implicit none
  integer, parameter         :: N_star = 101  !< Numero de estrellas + 1
  real, parameter            :: MassBH = 4.e6*msun
  real, dimension(N_star)    :: x_star, y_star, z_star,  vx_star ,vy_star ,vz_star ,massp
  !real, dimension(N_star)    :: Fx,Fy,Fz
  real, dimension(N_star)    :: xp,yp,zp,vxp,vyp,vzp
  !real                       :: FFx,FFy,FFz
  CHARACTER (40)             :: Fileout,fileoutt
  integer          :: im,i,j
  !
contains
  !
  !================================================================
  ! Se inicializa en N-body
  subroutine init_nbody()
    use constants

    call openfile()

  end subroutine init_nbody
  !================================================================
  !
  subroutine do_nbody(t,dt)!,x_star,y_star,z_star)

    use constants
    implicit none
    real, INTENT(IN)        :: dt, t
    !real,  INTENT(OUT)  :: x_star(N_star), y_star(N_star), z_star(N_star)

    !
    !DO WHILE(TMIN.LE.TFIN)
    ! -------------------------
    !    Paso en el tiempo  !
    !--------------------------
    call nbody_step(T,DT)  ! Modifica Xp, vxp
    !                                  ^   ^
    do im=1,N_star
      x_star(im)=xp(im)
      y_star(im)=yp(im)
      z_star(im)=zp(im)
      vx_star(im)=vxp(im)
      vy_star(im)=vyp(im)
      vz_star(im)=vzp(im)
   enddo

   do i=2,N_star
      write(20,*)  T,x_star(I),y_star(I),z_star(I)
   enddo

end subroutine do_nbody
  !
  !================================================================
 subroutine openfile()

   implicit none
   character (40) :: File1
   Integer        :: i
   !
   x_star(1) =0.
   y_star(1) =0.
   z_star(1) =0.
   vx_star(1) = 0.
   vy_star(1) = 0.
   vz_star(1) = 0.
   massp(1)   = MassBH
   !
!   write(File1,'("./IC/stars_ran_Vy_20.dat") ')
   write(File1,'("./IC/IC_p1e.dat") ')
   open(unit=10,file=File1,status='old')
   do i=2,N_star
      read(10,*,end=2) x_star(i),y_star(i),z_star(i),vx_star(i),vy_star(i),vz_star(i),massp(i)
   enddo
2  continue
   close(unit=10)

   !
 end subroutine openfile
  !
  !================================================================
 subroutine Force(im,x,y,z,Fx,Fy,Fz)
   implicit none
   integer, intent(in)  :: im
   real, dimension(N_star), intent(in)  :: x, y, z
   real,   intent(out)  :: Fx, Fy, Fz
   integer  :: j
   real     :: FFx,FFy,FFz
   !
   Fx = 0.
   Fy = 0.
   Fz = 0.
   !
   do j=1,N_star
      if(j.ne.im) then
         call Gforce(im,j,x(im),y(im),z(im),x(j),y(j),z(j),FFx,FFy,FFz)
         Fx =Fx + FFx
         Fy =Fy + FFy
         Fz =Fz + FFz
      endif
   enddo
 end subroutine Force

  !================================================================

 subroutine Gforce(i,j,x1,y1,z1,x2,y2,z2,FFx,FFy,FFz)
   implicit none
   integer, intent(in) :: i,j
   real, intent(in)        :: x1, y1, z1, x2, y2, z2
   real, intent(out)     :: FFx, FFy, FFz
   real     :: R3,GR
   real     :: SOFT

   SOFT=1.e10
   !
   R3=sqrt(((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)+SOFT**2)
   GR=Ggrav*massp(i)*massp(j)/(R3)**3
   FFx=GR*(x2-x1)
   FFy=GR*(y2-y1)
   FFz=GR*(z2-z1)
   !
 end subroutine Gforce
  !================================================================
 subroutine nbody_step(t,dt)
   implicit none
   real, INTENT(IN)         ::  t, dt
   real, dimension(N_star)  :: Fxp,Fyp,Fzp
   real        :: dt2, Fx,Fy,Fz
   integer :: im,i
   !
   dt2=0.5*dt
   !
   do im=1,N_star
      call Force(im,x_star,y_star,z_star,Fx,Fy,Fz)
      xp(im)=x_star(im)+dt2*vx_star(im)
      yp(im)=y_star(im)+dt2*vy_star(im)
      zp(im)=z_star(im)+dt2*vz_star(im)
      !
      vxp(im)=vx_star(im)+dt2*Fx/massp(im)
      vyp(im)=vy_star(im)+dt2*Fy/massp(im)
      vzp(im)=vz_star(im)+dt2*Fz/massp(im)
      !
   enddo
   !
   do im=1,N_star
      call Force(im, xp, yp, zp, Fx, Fy, Fz)
      !
      x_star(im)=x_star(im)+dt*vxp(im)
      y_star(im)=y_star(im)+dt*vyp(im)
      z_star(im)=z_star(im)+dt*vzp(im)
      !
      vxp(im)=vx_star(im)+dt*Fx/massp(im)
      vyp(im)=vy_star(im)+dt*Fy/massp(im)
      vzp(im)=vz_star(im)+dt*Fz/massp(im)
      !
   enddo
   !
   do im=1,N_star
      xp(im)=x_star(im)
      yp(im)=y_star(im)
      zp(im)=z_star(im)
   end do
   !
   !T=T+DT
   !
 end subroutine nbody_step
  !
  !================================================================

end module nbody

!================================================================
