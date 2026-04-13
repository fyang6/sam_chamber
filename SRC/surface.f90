subroutine surface()
! code modified for Pi chamber simulation when suing CL scheme. Fan Yang @ BNL
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor
use micro_prm, only: ncloud, naerosol
implicit none
	
real qvs(0:nx,1-YES3D:ny),t_s, q_s, u_h0, t_f, q_f
real qv_b(0:nx,1-YES3D:ny), qv_t(0:nx,1-YES3D:ny)
real qv_l(1-YES3D:ny,nzm), qv_r(1-YES3D:ny,nzm)
real qv_q(0:nx,nzm), qv_h(0:nx,nzm)
real qv_point, t_point,ql_point

real taux0, tauy0, xlmo
real diag_ustar, coef, coef1
integer i,j,k, j_gl,i_gl
integer ib, it, jb, jt, kb, kt
real(8) buffer(2), buffer1(2)

! LES mode: 

call t_startf ('surface')


if(.not.SFC_FLX_FXD) then

  if(LAND) then

    if(LES) then    
!   surface flux
      coef = (1000./pres0)**(rgas/cp)
      coef1 = (1000./pres(1))**(rgas/cp)
      t_s = ts_b
      q_s = soil_b*qsatw(t_s,pres(1))
! non-constant surface flux
      qv_b(0:nx,1-YES3D:ny) = micro_field(0:nx,1-YES3D:ny,1,index_water_vapor) ! index_water_vapor = 1, which is total water instead of qv.  AW@PNNL

      do j=1,ny
        do i=1,nx
          if (dos_b) then
            t_point = t(i,j,1) - gamaz(1) + fac_cond * qcl(i,j,1)   ! correct tabs calculation for cloud chamber condition
            call landflx(t_point*coef1, t_s*coef,  &
                         qv(i,j,1), q_s, 0.5*(u(i+1,j,1)+u(i,j,1)),       &
                         0.5*(v(i,j+YES3D,1)+v(i,j,1)), z(1), z0,         &
                         fluxt0, fluxq0, taux0, tauy0)
            fluxbt(i,j) = fluxt0
            fluxbq(i,j) = fluxq0
          else
            fluxbt(i,j) = 0e0
            fluxbq(i,j) = 0e0
          end if

          if (i==1) then
            ql_point = sum(micro_field(i-1,j,1,naerosol+2:naerosol+ncloud+1)) 
            qv_point = 0.5*(qv(i,j,1) + qv_b(i-1,j) - ql_point)
            t_point = 0.5*(t(i,j,1)+t(i-1,j,1)+fac_cond*(qcl(i,j,1)+ql_point))-gamaz(1)
          else
            qv_point = 0.5*(qv(i,j,1) + qv(i-1,j,1))
            t_point = 0.5*(t(i,j,1)+t(i-1,j,1)+fac_cond*(qcl(i,j,1)+qcl(i-1,j,1)))-gamaz(1)
          end if
          call landflx(t_point*coef1, & 
                      t_s*coef, qv_point, q_s, u(i,j,1)+ug,     & ! qv instead of qvs should be used here.  AW@PNNL
                      0.25*(v(i-1,j+YES3D,1)+v(i-1,j,1)+v(i,j+YES3D,1)+v(i,j,1))+vg, &
                      z(1), z0,fluxt0, fluxq0, taux0, tauy0)
          fluxbu(i,j)=taux0
              
          if (j==1) then
            ql_point = sum(micro_field(i,j-YES3D,1,naerosol+2:naerosol+ncloud+1))
            qv_point = 0.5*(qv(i,j,1) + qv_b(i,j-YES3D) - ql_point)
            t_point = 0.5*(t(i,j,1)+t(i,j-YES3D,1)+fac_cond*(qcl(i,j,1)+ql_point))-gamaz(1)
          else
            qv_point = 0.5*(qv(i,j,1)+qv(i,j-YES3D,1))
            t_point = 0.5*(t(i,j,1)+t(i,j-YES3D,1)+fac_cond*(qcl(i,j,1)+qcl(i,j-1,1)))-gamaz(1)
          end if
          call landflx(t_point*coef1, & 
                      t_s*coef, qv_point, q_s, & ! qv instead of qvs should be used here.  AW@PNNL
                      0.25*(u(i,j,1)+u(i+1,j,1)+u(i,j-YES3D,1)+u(i+1,j-YES3D,1))+ug, &
                      v(i,j,1)+vg, &
                      z(1), z0,fluxt0, fluxq0, taux0, tauy0)
          fluxbv(i,j)=tauy0
        end do
      end do
!---------------------------------------------------------------
! top flux add for pi chamber
      coef = (1000./(pres(nzm)+0.5*(pres(nzm)-pres(nzm-1))))**(rgas/cp)
      coef1 = (1000./pres(nzm))**(rgas/cp)
      t_s = ts_t
      q_s = soil_t*qsatw(t_s,pres(nzm))
! non-constant top surface flux
      qv_t(0:nx,1-YES3D:ny) = micro_field(0:nx,1-YES3D:ny,nzm,index_water_vapor) ! index_water_vapor = 1, which is total water instead of qv.  AW@PNNL
      do j=1,ny
        jb = j-YES3D
        jt = j+YES3D
        do i=1,nx
            ib = i-1
            it = i+1
            if (dos_t) then
            t_point = t(i,j,nzm) - gamaz(nzm) + fac_cond*qcl(i,j,nzm)
            call landflx(t_s*coef,t_point*coef1,  &
                       q_s,qv(i,j,nzm), 0.5*(u(it,j,nzm)+u(i,j,nzm)),       &
                       0.5*(v(i,jt,nzm)+v(i,j,nzm)), z(1), z0,          &
                       fluxt0, fluxq0, taux0, tauy0)
             fluxtt(i,j) = fluxt0
             fluxtq(i,j) = fluxq0
           else
             fluxtt(i,j) = 0e0
             fluxtq(i,j) = 0e0
           end if

           if (i==1) then
             ql_point = sum(micro_field(ib,j,nzm,naerosol+2:naerosol+ncloud+1))
             qv_point = 0.5*(qv(i,j,nzm) + qv_t(ib,j) - ql_point)
             t_point = 0.5*(t(i,j,nzm)+t(ib,j,nzm)+fac_cond*(qcl(i,j,nzm)+ql_point))-gamaz(nzm)
           else
             qv_point = 0.5*(qv(i,j,nzm) + qv(ib,j,nzm))
             t_point = 0.5*(t(i,j,nzm)+t(ib,j,nzm)+fac_cond*(qcl(i,j,nzm)+qcl(ib,j,nzm)))-gamaz(nzm)
           end if

           call landflx(t_s*coef,t_point*coef1, & 
                       q_s, qv_point, u(i,j,nzm),     & ! qv instead of qvs should be used here.  AW@PNNL
                       0.25*(v(ib,jt,nzm)+v(ib,j,nzm)+v(i,jt,nzm)+v(i,j,nzm)), &
                       z(1), z0,fluxt0, fluxq0, taux0, tauy0)
           fluxtu(i,j)= -taux0
               
           if (j==1) then
             ql_point = sum(micro_field(i,jb,nzm,naerosol+2:naerosol+ncloud+1))
             qv_point = 0.5*(qv(i,j,nzm) + qv_t(i,jb) - ql_point)
             t_point = 0.5*(t(i,j,nzm)+t(i,jb,nzm)+fac_cond*(qcl(i,j,nzm)+ql_point))-gamaz(nzm)
           else
             qv_point = 0.5*(qv(i,j,nzm) + qv(i,jb,nzm))
             t_point = 0.5*(t(i,j,nzm)+t(i,jb,nzm)+fac_cond*(qcl(i,j,nzm)+qcl(i,jb,nzm)))-gamaz(nzm)
           end if
           call landflx(t_s*coef,t_point*coef1, & 
                       q_s, qv_point, & ! qv instead of qvs should be used here.  AW@PNNL
                       0.25*(u(i,j,nzm)+u(it,j,nzm)+u(i,jb,nzm)+u(it,jb,nzm)), &
                       v(i,j,nzm), &
                       z(1), z0,fluxt0, fluxq0, taux0, tauy0)
                       fluxtv(i,j)= -tauy0
         end do
       end do
!---------------------------------------------------------------
! wall flux add for pi chamber
if(dowallx) then
! left wall flux
   qv_l(1-YES3D:ny,:) = micro_field(1,1-YES3D:ny,:,index_water_vapor) ! index_water_vapor = 1, which is total water instead of qv.  AW@PNNL
! set soil_wetness to 0.0 for dry side walls
  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
     kb = max(1,k-1)   
     do j=1,ny
       jb = j-YES3D
       jt = j+YES3D
       coef = (1000./pres(k))**(rgas/cp)
       coef1 = (1000./pres(k))**(rgas/cp)
! Side wall with checker-board pattern of temp
       j_gl = (int(rank/nsubdomains_y) * ny) +j
       if (j_gl .le. ny_gl/2) then 
          if (k .le. nzm/6) then
             t_s = ts_t
          else if ((k .gt. nzm/6) .and. (k .le. 2*nzm/6))  then
             t_s = ts_b
          else if ((k .gt. 2*nzm/6) .and. (k .le. 3*nzm/6))  then
             t_s = ts_t
          else if ((k .gt. 3*nzm/6) .and. (k .le. 4*nzm/6)) then
             t_s = ts_b
          else if ((k .gt. 4*nzm/6) .and. (k .le. 5*nzm/6)) then
             t_s = ts_t
          else if ((k .gt. 5*nzm/6) .and. (k .le. 6*nzm/6)) then
             t_s = ts_b
          end if

       else if (j_gl .gt. ny_gl/2) then
          if (k .le. nzm/6) then
             t_s = ts_b
          else if ((k .gt. nzm/6) .and. (k .le. 2*nzm/6))  then
             t_s = ts_t
          else if ((k .gt. 2*nzm/6) .and. (k .le. 3*nzm/6))  then
             t_s = ts_b
          else if ((k .gt. 3*nzm/6) .and. (k .le. 4*nzm/6)) then
             t_s = ts_t
          else if ((k .gt. 4*nzm/6) .and. (k .le. 5*nzm/6)) then
             t_s = ts_b
          else if ((k .gt. 5*nzm/6) .and. (k .le. 6*nzm/6)) then
             t_s = ts_t
          end if
       end if

       q_s = soil_w*qsatw(t_s,pres(k))
       if (dos_w) then
         t_point = t(1,j,k) - gamaz(k) + fac_cond*qcl(1,j,k)
         call landflxSW(t_point*coef1, t_s*coef, &
                     qv(1,j,k), q_s, 0.5*(v(1,j,k)+v(1,jt,k)), &
                     0.5*(w(1,j,k)+w(1,j,k+1)),z(1),z0, &
                     fluxt0, fluxq0, taux0, tauy0)
         fluxlt(j,k) = fluxt0
         fluxlq(j,k) = fluxq0
       else
         fluxlt(j,k) = 0e0
         fluxlq(j,k) = 0e0
       end if

       if (j==1) then
         ql_point = sum(micro_field(1,jb,k,naerosol+2:naerosol+ncloud+1))
         qv_point = 0.5*(qv(1,j,k) + qv_l(jb,k) - ql_point)
         t_point = 0.5*(t(1,j,k)+t(1,jb,k)+fac_cond*(qcl(1,j,k)+ql_point))-gamaz(k)
       else
         qv_point = 0.5*(qv(1,j,k) + qv(1,jb,k))
         t_point = 0.5*(t(1,j,k)+t(1,jb,k)+fac_cond*(qcl(1,j,k)+qcl(1,jb,k)))-gamaz(k)
       end if
       call landflxSW(t_point*coef1, &
            t_s*coef, qv_point, q_s, v(1,j,k),     & ! qv instead of qvy should be used here.  AW@PNNL
            0.25*(w(1,jb,k)+w(1,jb,k+1)+w(1,j,k+1)+w(1,j,k)), &
            z(1), z0,fluxt0, fluxq0, taux0, tauy0)
       fluxlv(j,k)=taux0

       qv_point = 0.5*(qv(1,j,k) + qv(1,j,kb))
       t_point = 0.5*(t(1,j,k)+t(1,j,kb)+fac_cond*(qcl(1,j,k)+qcl(1,j,kb)))-gamaz(k)

       call landflxSW(t_point*coef1, &
            t_s*coef, qv_point, q_s, & ! qv instead of qvy should be used here.  AW@PNNL
            0.25*(v(1,j,k)+v(1,j,kb)+v(1,jt,k)+v(1,jt,kb))+vg, &
            w(1,j,k), &
            z(1), z0,fluxt0, fluxq0, taux0, tauy0)
       fluxlw(j,k)=tauy0
     end do
    end do
  end if

! right wall flux
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
   qv_r(1-YES3D:ny,:) = micro_field(nx,1-YES3D:ny,:,index_water_vapor) ! index_water_vapor = 1, which is total water instead of qv.  AW@PNNL
   do k=1,nzm
     kb = max(1,k-1) 
     do j=1,ny
       jb = j-YES3D
       jt = j+YES3D
       coef = (1000./pres(k))**(rgas/cp)
       coef1 = (1000./pres(k))**(rgas/cp)
! Side wall with checker-board pattern of temp
    j_gl = (int(rank/nsubdomains_y)*ny)+j
       if (j_gl .le. ny_gl/2) then
          if (k .le. nzm/6) then
             t_s = ts_b
          else if ((k .gt. nzm/6) .and. (k .le. 2*nzm/6))  then
             t_s = ts_t
          else if ((k .gt. 2*nzm/6) .and. (k .le. 3*nzm/6))  then
             t_s = ts_b
          else if ((k .gt. 3*nzm/6) .and. (k .le. 4*nzm/6)) then
             t_s = ts_t
          else if ((k .gt. 4*nzm/6) .and. (k .le. 5*nzm/6)) then
             t_s = ts_b
          else if ((k .gt. 5*nzm/6) .and. (k .le. 6*nzm/6)) then
             t_s = ts_t
          end if

       else if (j_gl .gt. ny_gl/2) then
          if (k .le. nzm/6) then
             t_s = ts_t
          else if ((k .gt. nzm/6) .and. (k .le. 2*nzm/6))  then
             t_s = ts_b
          else if ((k .gt. 2*nzm/6) .and. (k .le. 3*nzm/6))  then
             t_s = ts_t
          else if ((k .gt. 3*nzm/6) .and. (k .le. 4*nzm/6)) then
             t_s = ts_b
          else if ((k .gt. 4*nzm/6) .and. (k .le. 5*nzm/6)) then
             t_s = ts_t
          else if ((k .gt. 5*nzm/6) .and. (k .le. 6*nzm/6)) then
             t_s = ts_b
          end if
       end if
       
       q_s = soil_e*qsatw(t_s,pres(k))
       t_f = max(t_s,(t(nx,j,k)-gamaz(k))*coef1)
       q_f = max(q_s,qv(nx,j,k))
       if (dos_e) then
         t_point = t(nx,j,k) - gamaz(k) + fac_cond*qcl(nx,j,k)
         call landflxSW(t_point*coef1, t_s*coef, &
                   qv(nx,j,k), q_s, 0.5*(v(nx,j,k)+v(nx,jt,k)), &
                   0.5*(w(nx,j,k)+w(nx,j,k+1)),z(1),z0, &
                   fluxt0, fluxq0, taux0, tauy0)
         fluxrt(j,k) = - fluxt0
         fluxrq(j,k) = - fluxq0
       else
         fluxrt(j,k) = 0e0
         fluxrq(j,k) = 0e0
       end if

       if (j==1) then
         ql_point = sum(micro_field(nx,jb,k,naerosol+2:naerosol+ncloud+1))
         qv_point = 0.5*(qv(nx,j,k) + qv_r(jb,k) - ql_point)
         t_point = 0.5*(t(nx,j,k)+t(nx,jb,k)+fac_cond*(qcl(nx,j,k)+ql_point))-gamaz(k)
       else
         qv_point = 0.5*(qv(nx,j,k) + qv(nx,jb,k))
         t_point = 0.5*(t(nx,j,k)+t(nx,jb,k)+fac_cond*(qcl(nx,j,k)+qcl(nx,jb,k)))-gamaz(k)
       end if
       call landflxSW(t_point*coef1, & 
            t_s*coef, qv_point, q_s, v(nx,j,k),     & ! qv instead of qvy should be used here.  AW@PNNL
            0.25*(w(nx,jb,k)+w(nx,jb,k+1)+w(nx,j,k+1)+w(nx,j,k)), &
            z(1), z0,fluxt0, fluxq0, taux0, tauy0)
       fluxrv(j,k)=-taux0

       qv_point = 0.5*(qv(nx,j,k) + qv(nx,j,kb))
       t_point = 0.5*(t(nx,j,k)+t(nx,j,kb)+fac_cond*(qcl(nx,j,k)+qcl(nx,j,kb)))-gamaz(k)

       call landflxSW(t_point*coef1, & 
            t_s*coef, qv_point, q_s, & ! qv instead of qvy should be used here.  AW@PNNL
            0.25*(v(nx,j,k)+v(nx,j,kb)+v(nx,jt,k)+v(nx,jt,kb))+vg, &
            w(nx,j,k), &
            z(1), z0,fluxt0, fluxq0, taux0, tauy0)
       fluxrw(j,k)=-tauy0
     end do
    end do
  end if

end if

if(dowally) then
! front wall flux
  if(rank.lt.nsubdomains_x) then   
    qv_q(0:nx,:) = micro_field(0:nx,1,:,index_water_vapor) ! index_water_vapor = 1, which is total water instead of qv.  AW@PNNL
    do k=1,nzm
     kb = max(1,k-1)
     do i=1,nx
       ib = i-1
       it = i+1
       coef = (1000./pres(k))**(rgas/cp)
       coef1 = (1000./pres(k))**(rgas/cp)
! Side wall with checker-board pattern of temp
       i_gl = (rank * nx)+i
       if (i_gl .le. nx_gl/2) then
          if (k .le. nzm/6) then
             t_s = ts_b
          else if ((k .gt. nzm/6) .and. (k .le. 2*nzm/6))  then
             t_s = ts_t
          else if ((k .gt. 2*nzm/6) .and. (k .le. 3*nzm/6))  then
             t_s = ts_b
          else if ((k .gt. 3*nzm/6) .and. (k .le. 4*nzm/6)) then
             t_s = ts_t
          else if ((k .gt. 4*nzm/6) .and. (k .le. 5*nzm/6)) then
             t_s = ts_b
          else if ((k .gt. 5*nzm/6) .and. (k .le. 6*nzm/6)) then
             t_s = ts_t
          end if
       
       else if (i_gl .gt. nx_gl/2) then
          if (k .le. nzm/6) then
             t_s = ts_t
          else if ((k .gt. nzm/6) .and. (k .le. 2*nzm/6))  then
             t_s = ts_b
          else if ((k .gt. 2*nzm/6) .and. (k .le. 3*nzm/6))  then
             t_s = ts_t
          else if ((k .gt. 3*nzm/6) .and. (k .le. 4*nzm/6)) then
             t_s = ts_b
          else if ((k .gt. 4*nzm/6) .and. (k .le. 5*nzm/6)) then
             t_s = ts_t
          else if ((k .gt. 5*nzm/6) .and. (k .le. 6*nzm/6)) then
             t_s = ts_b
          end if 
       end if       
       q_s = soil_s*qsatw(t_s,pres(k))
       if (dos_s) then
         t_point = t(i,1,k) - gamaz(k) + fac_cond*qcl(i,1,k)
         call landflxSW(t_point*coef1, t_s*coef, &
                   qv(i,1,k), q_s, 0.5*(w(i,1,k)+w(i,1,k+1)), &
                   0.5*(u(i,1,k)+u(it,1,k)),z(1),z0, &
                   fluxt0, fluxq0, taux0, tauy0)
         fluxqt(i,k) = fluxt0
         fluxqq(i,k) = fluxq0
       else
         fluxqt(i,k) = 0e0
         fluxqq(i,k) = 0e0
       end if

       qv_point = 0.5*(qv(i,1,k) + qv(i,1,kb))
       t_point = 0.5*(t(i,1,k)+t(i,1,kb)+fac_cond*(qcl(i,1,k)+qcl(i,1,kb)))-gamaz(k)

       call landflxSW(t_point*coef1, & 
            t_s*coef, qv_point, q_s, w(i,1,k),     & ! qv instead of qvx should be used here.  AW@PNNL
            0.25*(u(i,1,k)+u(it,1,k)+u(it,1,kb)+u(i,1,kb))+ug, &
            z(1), z0,fluxt0, fluxq0, taux0, tauy0)
       fluxqw(i,k)=taux0

       if (i==1) then
         ql_point = sum(micro_field(ib,1,k,naerosol+2:naerosol+ncloud+1))
         qv_point = 0.5*(qv(i,1,k) + qv_q(ib,k) - ql_point)
         t_point = 0.5*(t(i,1,k)+t(ib,1,k)+fac_cond*(qcl(i,1,k)+ql_point))-gamaz(k)
       else
         qv_point = 0.5*(qv(i,1,k) + qv(ib,1,k))
         t_point = 0.5*(t(i,1,k)+t(ib,1,k)+fac_cond*(qcl(i,1,k)+qcl(ib,1,k)))-gamaz(k)
       end if
       call landflxSW(t_point*coef1, & 
            t_s*coef, qv_point, q_s, & ! qv has been correctly applied here...  AW@PNNL
            0.25*(w(i,1,k)+w(i,1,k+1)+w(ib,1,k)+w(ib,1,k+1)), &
            u(i,1,k), &
            z(1), z0,fluxt0, fluxq0, taux0, tauy0)
       fluxqu(i,k)=tauy0
     end do
    end do
  end if
! back wall flux
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    qv_h(0:nx,:) = micro_field(0:nx,ny,:,index_water_vapor) ! index_water_vapor = 1, which is total water instead of qv.  AW@PNNL
    do k=1,nzm
     kb = max(1,k-1)
     do i=1,nx
       ib = i-1
       it = i+1
       coef = (1000./pres(k))**(rgas/cp)
       coef1 = (1000./pres(k))**(rgas/cp)
! Side wall with checker-board pattern of temp
       i_gl = ((rank-nx)*nx) + i

       if (i_gl .le. nx_gl/2) then
          if (k .le. nzm/6) then
             t_s = ts_b
          else if ((k .gt. nzm/6) .and. (k .le. 2*nzm/6))  then
             t_s = ts_t
          else if ((k .gt. 2*nzm/6) .and. (k .le. 3*nzm/6))  then
             t_s = ts_b
          else if ((k .gt. 3*nzm/6) .and. (k .le. 4*nzm/6)) then
             t_s = ts_t
          else if ((k .gt. 4*nzm/6) .and. (k .le. 5*nzm/6)) then
             t_s = ts_b
          else if ((k .gt. 5*nzm/6) .and. (k .le. 6*nzm/6)) then
             t_s = ts_t
          end if

       else if (i_gl .gt. nx_gl/2) then
          if (k .le. nzm/6) then
             t_s = ts_t
          else if ((k .gt. nzm/6) .and. (k .le. 2*nzm/6))  then
             t_s = ts_b
          else if ((k .gt. 2*nzm/6) .and. (k .le. 3*nzm/6))  then
             t_s = ts_t
          else if ((k .gt. 3*nzm/6) .and. (k .le. 4*nzm/6)) then
             t_s = ts_b
          else if ((k .gt. 4*nzm/6) .and. (k .le. 5*nzm/6)) then
             t_s = ts_t
          else if ((k .gt. 5*nzm/6) .and. (k .le. 6*nzm/6)) then
             t_s = ts_b
          end if
       end if
       q_s = soil_n*qsatw(t_s,pres(k))
       if (dos_n) then
         t_point = t(i,ny,k) - gamaz(k) + fac_cond*qcl(i,ny,k)
         call landflxSW(t_point*coef1, t_s*coef, &
                   qv(i,ny,k), q_s, 0.5*(w(i,ny,k)+w(i,ny,k+1)), &
                   0.5*(u(i,ny,k)+u(it,ny,k)),z(1),z0, &
                   fluxt0, fluxq0, taux0, tauy0)
         fluxht(i,k) = - fluxt0
         fluxhq(i,k) = - fluxq0
       else
         fluxht(i,k) = 0e0
         fluxhq(i,k) = 0e0
       end if

       qv_point = 0.5*(qv(i,ny,k) + qv(i,ny,kb))
       t_point = 0.5*(t(i,ny,k)+t(i,ny,kb)+fac_cond*(qcl(i,ny,k)+qcl(i,ny,kb)))-gamaz(k)

       call landflxSW(t_point*coef1, & 
            t_s*coef, qv_point, q_s, w(i,ny,k),     & ! qv instead of qvx should be used here.  AW@PNNL
            0.25*(u(i,ny,k)+u(it,ny,k)+u(it,ny,kb)+u(i,ny,kb))+ug, &
            z(1), z0,fluxt0, fluxq0, taux0, tauy0)
       fluxhw(i,k)= - taux0
     
       if (i==1) then
         ql_point = sum(micro_field(ib,ny,k,naerosol+2:naerosol+ncloud+1))
         qv_point = 0.5*(qv(i,ny,k) + qv_h(ib,k) - ql_point)
         t_point = 0.5*(t(i,ny,k)+t(ib,ny,k)+fac_cond*(qcl(i,ny,k)+ql_point))-gamaz(k)
       else
         qv_point = 0.5*(qv(ib,ny,k) + qv(i,ny,k))
         t_point = 0.5*(t(i,ny,k)+t(ib,ny,k)+fac_cond*(qcl(i,ny,k)+qcl(ib,ny,k)))-gamaz(k)
       end if
       call landflxSW(t_point*coef1, & 
            t_s*coef, qv_point, q_s, & ! qv instead of qvx should be used here.  AW@PNNL
            0.25*(w(i,ny,k)+w(i,ny,k+1)+w(ib,ny,k)+w(ib,ny,k+1)), &
            u(i,ny,k)+ug, &
            z(1), z0,fluxt0, fluxq0, taux0, tauy0)
       fluxhu(i,k)= - tauy0
     end do
    end do
  end if

end if

!-----------------------------------------
            end if ! LES

            if(CEM) then

              coef = (1000./pres0)**(rgas/cp)
              coef1 = (1000./pres(1))**(rgas/cp)
              qvs(0:nx,1-YES3D:ny) = micro_field(0:nx,1-YES3D:ny,1,index_water_vapor)

              do j=1,ny  
               do i=1,nx

               t_s = (sstxy(i,j)+t00)*coef
               q_s = soil_wetness*qsatw(sstxy(i,j)+t00,pres(1))
               call landflx((t(i,j,1)-gamaz(1))*coef1, t_s,   &
                      qv(i,j,1), q_s, 0.5*(u(i+1,j,1)+u(i,j,1))+ug,     &
                        0.5*(v(i,j+YES3D,1)+v(i,j,1))+vg, z(1), z0,        &
                      fluxt0, fluxq0, taux0, tauy0)
               fluxbt(i,j) = fluxt0
               fluxbq(i,j) = fluxq0

               t_s = (0.5*(sstxy(i-1,j)+sstxy(i,j))+t00)*coef
               q_s = soil_wetness*qsatw(0.5*(sstxy(i-1,j)+sstxy(i,j))+t00,pres(1))
               call landflx((0.5*(t(i-1,j,1)+t(i,j,1))-gamaz(1))*coef1, t_s,   &
                      0.5*(qvs(i-1,j)+qvs(i,j)), q_s, u(i,j,1)+ug,     &
                        0.25*(v(i-1,j+YES3D,1)+v(i-1,j,1)+v(i,j+YES3D,1)+v(i,j,1))+vg, &
                       z(1), z0, fluxt0, fluxq0, taux0, tauy0)
               if(SFC_TAU_FXD) then
                   u_h0 = max(1.,sqrt((u(i,j,1)+ug)**2+ &
                        (0.25*(v(i-1,j+YES3D,1)+v(i-1,j,1)+v(i,j+YES3D,1)+v(i,j,1))+vg)**2))
                   taux0 = -(u(i,j,1)+ug)/u_h0*tau0*rhow(1)
               end if
               fluxbu(i,j) = taux0

               t_s = (0.5*(sstxy(i,j-YES3D)+sstxy(i,j))+t00)*coef
               q_s = soil_wetness*qsatw(0.5*(sstxy(i,j-YES3D)+sstxy(i,j))+t00,pres(1))
               call landflx((0.5*(t(i,j-YES3D,1)+t(i,j,1))-gamaz(1))*coef1, t_s,   &
                      0.5*(qvs(i,j-YES3D)+qvs(i,j)), q_s,  &
                      0.25*(u(i+1,j-YES3D,1)+u(i,j-YES3D,1)+u(i+1,j,1)+u(i,j,1))+ug,     &
                      v(i,j,1)+vg, &
                      z(1), z0, fluxt0, fluxq0, taux0, tauy0)
               if(SFC_TAU_FXD) then
                  u_h0 = max(1.,sqrt( &
                       (0.25*(u(i+1,j-YES3D,1)+u(i,j-YES3D,1)+u(i+1,j,1)+u(i,j,1))+ug)**2+ &
                       (v(i,j,1)+vg)**2))
                  tauy0 = -(v(i,j,1)+vg)/u_h0*tau0*rhow(1)
               end if
               fluxbv(i,j) = tauy0

               end do
              end do

            end if ! CEM


  end if ! LAND

end if! .not.SFC_FLX_FXD



if(SFC_FLX_FXD) then

  u_h0 = max(1.,sqrt((u0(1)+ug)**2+(v0(1)+vg)**2))

  if(.not.SFC_TAU_FXD) then
    if(OCEAN) z0 = 0.0001  ! for LAND z0 should be set in namelist (default z0=0.035)

    tau0 = diag_ustar(z(1),  &
                bet(1)*(fluxt0+epsv*(t0(1)-gamaz(1))*fluxq0),u_h0,z0)**2  

  end if ! .not.SFC_TAU_FXD

  if(LES) then
    taux0 = -(u0(1)+ug)/u_h0*tau0
    tauy0 = -(v0(1)+vg)/u_h0*tau0
    fluxbu(:,:) = taux0
    fluxbv(:,:) = tauy0
  else
    fluxbu(:,:) = -(u(1:nx,1:ny,1)+ug)/u_h0*tau0
    fluxbv(:,:) = -(v(1:nx,1:ny,1)+vg)/u_h0*tau0
  end if

  fluxbt(:,:) = fluxt0
  fluxbq(:,:) = fluxq0

end if ! SFC_FLX_FXD

!
! Homogenize the surface scalar fluxes if needed for sensitivity studies
!
   if(dosfchomo) then

	fluxt0 = 0.
	fluxq0 = 0.
	do j=1,ny
         do i=1,nx
	   fluxt0 = fluxt0 + fluxbt(i,j)
	   fluxq0 = fluxq0 + fluxbq(i,j)
         end do
        end do
	fluxt0 = fluxt0 / float(nx*ny)
	fluxq0 = fluxq0 / float(nx*ny)
        if(dompi) then
            buffer(1) = fluxt0
            buffer(2) = fluxq0
            call task_sum_real8(buffer,buffer1,2)
	    fluxt0 = buffer1(1) /float(nsubdomains)
	    fluxq0 = buffer1(2) /float(nsubdomains)
        end if ! dompi
	fluxbt(:,:) = fluxt0
	fluxbq(:,:) = fluxq0

   end if

shf_xy(:,:) = shf_xy(:,:) + fluxbt(:,:) * dtfactor
lhf_xy(:,:) = lhf_xy(:,:) + fluxbq(:,:) * dtfactor
! Ssve top and side wall's flux. Aaron Wang @ PNNL
shft_xy(:,:) = shft_xy(:,:) + fluxtt(:,:) * dtfactor
lhft_xy(:,:) = lhft_xy(:,:) + fluxtq(:,:) * dtfactor
! Left wall
if(mod(rank,nsubdomains_x).eq.0) then
  shfl_yz(:,:) = shfl_yz(:,:) + fluxlt(:,:) * dtfactor
  lhfl_yz(:,:) = lhfl_yz(:,:) + fluxlq(:,:) * dtfactor
end if
! Right wall
if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
  shfr_yz(:,:) = shfr_yz(:,:) + fluxrt(:,:) * dtfactor
  lhfr_yz(:,:) = lhfr_yz(:,:) + fluxrq(:,:) * dtfactor
end if
! Front wall
if(rank.lt.nsubdomains_x) then
  shfq_xz(:,:) = shfq_xz(:,:) + fluxqt(:,:) * dtfactor
  lhfq_xz(:,:) = lhfq_xz(:,:) + fluxqq(:,:) * dtfactor
end if
! Back wall
if(rank.gt.nsubdomains-nsubdomains_x-1) then
  shfh_xz(:,:) = shfh_xz(:,:) + fluxht(:,:) * dtfactor
  lhfh_xz(:,:) = lhfh_xz(:,:) + fluxhq(:,:) * dtfactor
end if
! End of editing. Aaron Wang @ PNNL

call t_stopf ('surface')

end




! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below 
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!
! Code corrected 8th June 1999 (obukhov length was wrong way up,
! so now used as reciprocal of obukhov length)

      real function diag_ustar(z,bflx,wnd,z0)

      implicit none
      real, parameter      :: vonk =  0.4   ! von Karmans constant
      real, parameter      :: g    = 9.81   ! gravitational acceleration
      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "
      real, parameter      :: eps  = 1.e-10 ! non-zero, small number

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z
      real, intent (in)    :: z0            ! momentum roughness height

      integer :: iterate
      real    :: lnz, klnz, c1, x, psi1, zeta, rlmo, ustar

      lnz   = log(z/z0) 
      klnz  = vonk/lnz              
      c1    = 3.14159/2. - 3.*log(2.)

      ustar =  wnd*klnz
      if (bflx /= 0.0) then 
        do iterate=1,4
          rlmo   = -bflx * vonk/(ustar**3 + eps)   !reciprocal of
                                                   !obukhov length
          zeta  = z*rlmo
          if (zeta > 0.) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
          else
            x     = sqrt( sqrt( 1.0 - bm*zeta ) )
            psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
            ustar = wnd*vonk/(lnz - psi1)
          end if
        end do
      end if

      diag_ustar = ustar

      return
      end function diag_ustar
! ----------------------------------------------------------------------

