program test
  implicit none
  integer :: nx, nxtot,i0,j0
  real    :: x, y, posx, posy, dx, remx, remy, lx0,lx1, ly0,ly1,distx,disty

  nxtot = 8
  nx   =  4
  dx = 1./nxtot

  x = 0.7
  y = 0.88

  posx = x/dx
  posy = y/dx

  remx = posx - int(posx)
  remy = posy - int(posy)

  print*, posx, posy
  !print*, int(posx), int(posy)
  print*, remx, remy

  print*, 'nearest cell is', nint(posx)+1, nint(posy)+1

  if (remx < 0.5) then
    i0= int(posx)
  else
    i0= int(posx)+1
  end if

  if (remy < 0.5) then
    j0= int(posy)
  else
    j0 = int(posy)+1
  end if

  distx = posx - ( real(i0)-0.5 )
  disty = posy - ( real(j0)-0.5 )

  print*, "bounds are: ",i0,i0+1,j0,j0+1
  print*, distx, disty

end program
