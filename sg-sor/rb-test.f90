program test
	implicit none
	integer :: i,j,k, ipass, isw,jsw,ksw, nx, ny, nz,count,jpass, kpass
	nx=4
	ny=4
  nz=4

  count=0

  print*,'total',count

  isw=1
  do ipass=1,2
    jsw=isw
      do i=2, nx-1
        do j=jsw+1,ny-1,2

          !ksw=isw
          !do k=ksw+1,nz-1,2
                
            print'(4i4)',i,j,k, ipass
            count=count+1
          !end do
              
          !ksw=3-ksw
         end do
         
        jsw=3-jsw
        end do
     
        print*,'++++++'
      isw=3-isw
    !end do
  end do

print*,count

end program test