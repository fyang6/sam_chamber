subroutine test_temperature()

use vars, only:t
use params, only: tabs_s, tabs_t
use grid, only: nx, ny, nz, nzm, dimx1_s, dimx2_s, dimy1_s, dimy2_s
integer i,j,k

! print *,minval(t), minloc(t)
! print *,maxval(t), maxloc(t)

! Checks for the temperature and ensures that the values exceeding the bounds
! are diffused to the surrounding cells
! Uppermost Layer St
 k=nz
 do i=1,nx
  do j=1,ny
  t(i,j,k) = (t(i+1,j,k)+t(i-1,j,k)+t(i,j+1,k)+t(i,j-1,k))/4.
  end do
 end do
 ! Lowermost Layer St
 k=1
 do i=1,nx
  do j=1,ny
  t(i,j,k) = (t(i+1,j,k)+t(i-1,j,k)+t(i,j+1,k)+t(i,j-1,k))/4.
  end do
 end do

! do k=1,nzm     
!  do j=1,ny
!   do i=1,nx
!    if(t(i,j,k).gt.(tabs_s).or.t(i,j,k).lt.(tabs_t)) then
!    print*,i,j,k,t(i,j,k)
!    t(i,j,k) = (t(i+1,j,k)+t(i-1,j,k)+t(i,j+1,k)+t(i,j-1,k))/4.
!    end if
!   end do
!  end do
! end do

end subroutine


