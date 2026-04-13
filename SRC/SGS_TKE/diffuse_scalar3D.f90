!modified based on origianl subroutine
!subroutine diffuse_scalar3D (field,fluxb,fluxt,tkh,rho,rhow,flux)
! FY @ MTU 2017
subroutine diffuse_scalar3D (field,fluxb,fluxt,fluxl,fluxr, &
                         fluxq,fluxh,tkh,rho,rhow,flux)

use grid
use params, only: docolumn,dowallx,dowally,dosgs,dowallsflux
!use vars
!use params
use sgs, only: grdf_x,grdf_y,grdf_z
implicit none
! input	
real field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real tkh(0:nxp1,1-YES3D:nyp1,nzm)	! eddy conductivity
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real fluxl(ny,nzm)
real fluxr(ny,nzm)
real fluxq(nx,nzm)
real fluxh(nx,nzm)
real rho(nzm)
real rhow(nz)
real flux(nz)
! local        
real flx(0:nx,0:ny,0:nzm)
real dfdt(nx,ny,nz)
real rdx2,rdy2,rdz2,rdx,rdy,rdz,rdx5,rdy5,rdz5,tmp
real dxy,dxz,dyx,dyz,dzx,dzy,tkx,tky,tkz,rhoi
real taux0, tauy0, xlmo
real coef, coef1
real t_s, q_s
integer i,j,k,ib,ic,jb,jc,kc,kb


if(.not.dosgs) return

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)
rdz2=1./(dz*dz)
rdx=1./dx
rdy=1./dy
rdz=1./dz
dxy=dx/dy
dxz=dx/dz
dyx=dy/dx
dyz=dy/dz
dzx=dz/dx
dzy=dz/dy

dfdt(:,:,:)=0.

!  Horizontal diffusion:
if (dowallx.or.dowally) then
! scalar flux for pi chamber
if (dowallsflux) then

 do k=1,nzm
 rdx5=0.5*rdx2*grdf_x(k)
 if (dowallx) then
  do j=1,ny
  if(mod(rank,nsubdomains_x).eq.0) then
!  left wall
   flx(0,j,k)=fluxl(j,k)*rdx
  else
   tkx=rdx5*(tkh(1,j,k)+tkh(0,j,k))
   flx(0,j,k)=-tkx*(field(1,j,k)-field(0,j,k))
  end if

  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
   flx(nx,j,k)=fluxr(j,k)*rdx
!  right wall
  else
   tkx=rdx5*(tkh(nx+1,j,k)+tkh(nx,j,k))
   flx(nx,j,k)=-tkx*(field(nx+1,j,k)-field(nx,j,k))
  end if
  end do ! end of j
!  normal condition
  do j=1,ny
   do i=1,nx-1
    ic=i+1
    tkx=rdx5*(tkh(i,j,k)+tkh(ic,j,k))
    flx(i,j,k)=-tkx*(field(ic,j,k)-field(i,j,k))
   end do
   do i=1,nx
    ib=i-1
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(ib,j,k))
   end do
  end do

 end if ! end of dowallx

 if (dowally) then
  rdy5=0.5*rdy2*grdf_y(k)
  do i=1,nx
  if(rank.lt.nsubdomains_x) then
!  front wall
   flx(i,0,k)=fluxq(i,k)*rdy
  else
   tky=rdy5*(tkh(i,0,k)+tkh(i+1,0,k))
   flx(i,0,k)=-tky*(field(i,1,k)-field(i,0,k))
  end if

  if(rank.gt.nsubdomains-nsubdomains_x-1) then
!  back wall
   flx(i,ny,k)=fluxh(i,k)*rdy
  else
   tky=rdy5*(tkh(i,ny+1,k)+tkh(i,ny,k))
   flx(i,ny,k)=-tky*(field(i,ny+1,k)-field(i,ny,k))
  end if
  end do
!  normal condition
  do j=1,ny-1
   do i=1,nx
    jc=j+1
    tky=rdy5*(tkh(i,jc,k)+tkh(i,j,k))
    flx(i,j,k)=-tky*(field(i,jc,k)-field(i,j,k))
   end do
  end do
  do j=1,ny
   jb=j-1
   do i=1,nx
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(i,jb,k))
   end do
  end do
 end if !end of dowally
end do ! end of do k loop

else ! else of dowallsflux
! no wall flux
 if (dowallx) then
  if (mod(rank,nsubdomains_x).eq.0) then
   do k=1,nzm
    do j=1,ny
     field(0,j,k) = field(1,j,k)
    end do
   end do
  end if

  if (mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
   do k=1,nzm
    do j=1,ny
     field(nx+1,j,k)=field(nx,j,k)
    end do
   end do
  end if
 end if ! end of dowallx

 if (dowally) then
  if (rank.lt.nsubdomains_x) then
   do k=1,nzm
    do i=1,nx
     field(i,1-YES3D,k) = field(i,1,k)
    end do
   end do
  end if

  if (rank.gt.nsubdomains-nsubdomains_x-1) then
   do k=1,nzm
    do i=1,nx
     field(i,ny+YES3D,k) = field(i,ny,k)
    end do
   end do
  end if
 end if ! end of dowally
end if ! end of dowallsflux

else ! periodic boundary-original code

do k=1,nzm
	
 rdx5=0.5*rdx2  * grdf_x(k)
 rdy5=0.5*rdy2  * grdf_y(k)

 do j=1,ny
  do i=0,nx
    ic=i+1
    tkx=rdx5*(tkh(i,j,k)+tkh(ic,j,k)) 	
    flx(i,j,k)=-tkx*(field(ic,j,k)-field(i,j,k))
  end do 
  do i=1,nx
    ib=i-1
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(ib,j,k))
  end do 
 end do 

 do j=0,ny
  jc=j+1
  do i=1,nx
   tky=rdy5*(tkh(i,j,k)+tkh(i,jc,k)) 	
   flx(i,j,k)=-tky*(field(i,jc,k)-field(i,j,k))
  end do 
 end do
 do j=1,ny
  jb=j-1
  do i=1,nx	    
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(i,jb,k))
  end do 
 end do 
 
end do ! k

end if ! end of dowallx/dowally

!  Vertical diffusion:

flux(1) = 0.
flux(nzm)=0.
tmp=1./adzw(nzm)
do j=1,ny
 do i=1,nx	
   flx(i,j,0)=fluxb(i,j)*rdz*rhow(1)
! add for pi chamber
!   flx(i,j,nzm)=-fluxt(i,j)*rdz*tmp*rhow(nz)
   flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz)
   flux(1) = flux(1) + flx(i,j,0)
 end do
end do

do k=1,nzm-1
 kc=k+1
 flux(kc)=0. 
 rhoi = rhow(kc)/adzw(kc)
 rdz5=0.5*rdz2 * grdf_z(k)
 do j=1,ny
  do i=1,nx
    tkz=rdz5*(tkh(i,j,k)+tkh(i,j,kc))
    flx(i,j,k)=-tkz*(field(i,j,kc)-field(i,j,k))*rhoi
    flux(kc) = flux(kc) + flx(i,j,k)
  end do 
 end do
end do

do k=1,nzm
 kb=k-1
 rhoi = 1./(adz(k)*rho(k))
 do j=1,ny
  do i=1,nx		 
   !dfdt(i,j,k)=dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)
   !field(i,j,k)=field(i,j,k)+dfdt(i,j,k)
   field(i,j,k)=field(i,j,k)+dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)
  end do 
 end do 	 
end do 

end subroutine diffuse_scalar3D
