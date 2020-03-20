program test
	implicit none
	integer :: i,j,k, ipass, isw,jsw,ksw,count,jpass, kpass,i_rb
	integer, parameter :: nx=5, ny=5, nz=5
	integer :: map(nx,ny,nz)

	map(:,:,:) = 0
	ksw=1
	black_red: do kpass=1,2

		do ipass=1,2
			jsw=ksw
			do k=1,nz
				do j=jsw,ny,2
					do i=ipass,nx,2

						map(i,j,k) = map(i,j,k) + kpass

					end do
				end do
				jsw=3-jsw
			end do
			ksw=3-ksw
		end do
	end do black_red


	PRINT*,'=============================='
	do k=1,nz
		do j=ny,1,-1
			print'(5i3)', map(:,j,k)
		end do
		print*,'-------'
	end do
	PRINT*,'=============================='

	map(:,:,:) = 0
	do i_rb =1,2
		jsw = i_rb
		do k=1,nz
			isw =jsw
			do j=1,ny
				do i=isw,nx,2

					map(i,j,k) = map(i,j,k) + i_rb

				end do
				isw = 3-isw
			end do
			jsw = 3-jsw
		end do
	end do

	PRINT*,'=============================='
	do k=1,nz
		do j=ny,1,-1
			print'(5i3)', map(:,j,k)
		end do
		print*,'-------'
	end do
	PRINT*,'=============================='

end program test
!  this shouldn't work but the code commented above should, and it
!  doesn't hav to check this later
! ksw=1
! black_red: do kpass=1,2
!
!   do ipass=1,2
!     jsw=ksw
!     do k=1,nz
!       do j=jsw,ny,2
!         do i=ipass,nx,2
!
!           call get_residue(i,j,k,xi)
!
!           ph0             = phi_grav(i,j,k)
!           phi_grav(i,j,k) = ph0 - omega*xi/e_ijk
!
!           relative_error  = abs(phi_grav(i,j,k)-ph0)/abs(ph0 + 1e-30)
!           max_error_local = max(max_error_local,relative_error)
!
!         end do
!       end do
!       jsw=3-jsw
!     end do
!     ksw=3-ksw
!   end do
!
! end do black_red
else
