program test
  use nbody
  implicit none
  real :: t, dt

  dt= 1e6
  t = 0.

  call init_nbody()
  !write(6,*),x_star(2),y_star(2),z_star(2)

  !print*,x_star(:),vx_star(:)
  !print*,''
  !print*,y_star(:),vy_star(:)
  !print*,''
  !print*,z_star(:),vz_star(:)
  !print''
  !print*, massp(:)
  !print*,'---'

  do while (t <= 1e10)
    call do_nbody(t, dt)
    t=t +dt
    !print*,x_star(:)
    !print*,''
    !print*,y_star(:)
    !print*,''
    !print*,z_star(:)
    !write(6,*),x_star(2),y_star(2),z_star(2)
    write(6,'(10(xe25.12))') t,x_star(1),y_star(1),z_star(1),x_star(2),y_star(2),z_star(2),x_star(3),y_star(3),z_star(3)
  enddo

end program
