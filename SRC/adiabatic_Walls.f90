subroutine adiabaticWalls
use vars, only:t
use params, only:tabs_s, tabs_t, tabs_w
use grid
use domain
integer i,j,k
! Left Wall St 
if(mod(rank,nsubdomains_x).eq.0) then
!print *,'l'
do k=1,nzm
! do i=dimx1_s,dimx1_s
 do i=dimx1_s,0
     do j=dimy1_s,dimy2_s
     t(i,j,k) = tabs_w
!     t(i,j,k) = tabs(1,j,k)
     end do
 end do
end do
end if
! Right Wall St 
if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
!print *,'r'
do k=1,nzm
! do i=dimx2_s,dimx2_s
 do i=nx+1,dimx2_s
     do j=dimy1_s,dimy2_s
     t(i,j,k) = tabs_w
!     t(i,j,k) = tabs(nx,j,k)
     end do
 end do
end do
end if
! Front Wall St 
if(rank.lt.nsubdomains_x) then
!print *,'f'
do k=1,nzm
 do i=dimx1_s,dimx2_s
!     do j=dimy1_s,dimy1_s
     do j=dimy1_s,0
     t(i,j,k) = tabs_w
!     t(i,j,k) = tabs(i,1,k) 
     end do
 end do
end do
end if
 
! Back Wall St 
if(rank.gt.nsubdomains-nsubdomains_x-1) then
!print *,'b'
do k=1,nzm
do i=dimx1_s,dimx2_s
!     do j=dimy2_s,dimy2_s
     do j=ny+1,dimy2_s
     t(i,j,k) = tabs_w
!     t(i,j,k) = tabs(i,ny,k) 
     end do
 end do
end do
end if

!Uppermost Layer St
!k=nz
!do i=dimx1_s,dimx2_s
!do j=dimx1_s, dimx2_s
!t(i,j,k) = tabs_t 
!end do
!end do
!
!Lowermost Layer St
!k=1
!do i=dimx1_s,dimx2_s
!do j=dimx1_s, dimx2_s
!t(i,j,k) = tabs_s 
!end do
!end do

end subroutine
