!------------------------------------------    
!  SUBROUTINE calc_statistics() 
!------------------------------------------    



!------------------------------------------    
!------------------------------------------    
   SUBROUTINE calc_statistics(stat_flag)

     USE global_parameters
     USE flow_parameters
     USE grid_arrays,ONLY : xc,yc,zc
     USE flow_arrays,ONLY : u,v,w
     USE pressure_arrays
     USE stat_arrays

     IMPLICIT NONE

!... Parameters 
     INTEGER, INTENT(IN):: stat_flag

!... Loop Variables
     INTEGER           :: i,j,k


!... Local Variables
     REAL(KIND=CGREAL) :: rStatCtr, rStatCtrSqr
     REAL(KIND=CGREAL) :: rStatCtrInv, rStatCtrSqrInv


     SELECT CASE(stat_flag)
       CASE (0) 
         DO k = 1,nzc
         DO j = 1,nyc
         DO i = 1,nxc
           uAv(i,j,k)   =  uAv(i,j,k) + u(i,j,k)
           vAv(i,j,k)   =  vAv(i,j,k) + v(i,j,k)
           wAv(i,j,k)   =  wAv(i,j,k) + w(i,j,k)
           pAv(i,j,k)   =  pAv(i,j,k) + p(i,j,k)
           uvAv(i,j,k)  =  uvAv(i,j,k)+ u(i,j,k)*v(i,j,k)
           vwAv(i,j,k)  =  vwAv(i,j,k)+ v(i,j,k)*w(i,j,k)
           uwAv(i,j,k)  =  uwAv(i,j,k)+ u(i,j,k)*w(i,j,k)
           uuAv(i,j,k)  =  uuAv(i,j,k)+ u(i,j,k)*u(i,j,k)
           vvAv(i,j,k)  =  vvAv(i,j,k)+ v(i,j,k)*v(i,j,k)
           wwAv(i,j,k)  =  wwAv(i,j,k)+ w(i,j,k)*w(i,j,k)
         ENDDO ! i
         ENDDO ! j
         ENDDO ! k

       CASE (1) 

! Open Files
!  Use 2 stat files to avoid clobbering problems.
 
         SELECT CASE(indexStat)
         CASE(1)
           OPEN(ifuStatOut,FILE='stat_out1.dat',FORM='UNFORMATTED')
           indexStat = 2
         CASE(2)
           OPEN(ifuStatOut,FILE='stat_out2.dat',FORM='UNFORMATTED')
           indexStat = 1
         END SELECT ! indexStat

         rStatCtr    = REAL(statCtr,KIND=CGREAL)
         rStatCtrSqr = rStatCtr**2
         rStatCtrInv = oned/rstatCtr
         rStatCtrSqrInv = oned/rstatCtrSqr

         WRITE(ifuStatOut)ntime-ntime_start+1,statCtr,dt,         &
                          uAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv,  &
                          vAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv,  &
                          wAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv,  &
                          pAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv,  &
                          uvAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv, &
                          vwAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv, &
                          uwAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv, &
                          uuAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv, &
                          vvAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv, &
                          wwAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInv

         close (ifuStatOut)

         OPEN(ifuStatPlot,FILE='stat_plot.dat')

         !WRITE(ifuStatPlot,*)'VARIABLES="X","Y","Z","UAV","VAV","WAV","PAV","UpUpAV","VpVpAV","WpWpAV","UpVpAV","VpWpAV","UpWpAV"'
         WRITE(ifuStatPlot,*)'VARIABLES="X","Y","Z","UAV","VAV","WAV","PAV"'
         WRITE(ifuStatPlot,*)'ZONE F=POINT, I=',nxc,',J=',nyc,'K=',nzc
         DO k = 1,nzc
         DO j = 1,nyc
         DO i = 1,nxc
!           WRITE(ifuStatPlot,101) xc(i),yc(j),zc(k),                                  &
!                        uAv(i,j,k) *rStatCtrInv,vAv(i,j,k)*rStatCtrInv                &
!                       ,wAv(i,j,k) *rStatCtrInv,pAv(i,j,k)*rStatCtrInv                &
!                       ,uuAv(i,j,k)*rStatCtrInv-uAv(i,j,k)*uAv(i,j,k)*rStatCtrSqrInv  &
!                       ,vvAv(i,j,k)*rStatCtrInv-vAv(i,j,k)*vAv(i,j,k)*rStatCtrSqrInv  &
!                       ,wwAv(i,j,k)*rStatCtrInv-wAv(i,j,k)*wAv(i,j,k)*rStatCtrSqrInv  &
!                       ,uvAv(i,j,k)*rStatCtrInv-uAv(i,j,k)*vAv(i,j,k)*rStatCtrSqrInv  &
!                       ,vwAv(i,j,k)*rStatCtrInv-vAv(i,j,k)*wAv(i,j,k)*rStatCtrSqrInv  &
!                       ,uwAv(i,j,k)*rStatCtrInv-uAv(i,j,k)*wAv(i,j,k)*rStatCtrSqrInv
           WRITE(ifuStatPlot,102) xc(i),yc(j),zc(k)                                   &
                                  ,uAv(i,j,k) *rStatCtrInv                             &
                                  ,vAv(i,j,k)*rStatCtrInv                              &
                                  ,wAv(i,j,k) *rStatCtrInv                             &
                                  ,pAv(i,j,k)*rStatCtrInv    
                                   
         ENDDO ! i
         ENDDO ! i
         ENDDO ! k

         CLOSE(ifuStatPlot)
     END SELECT ! stat_flag

101  FORMAT(13(2x,e14.7))
102  FORMAT(7(2x,e14.7))
   END SUBROUTINE calc_statistics
!------------------------------------------   

!------------------------------------------
!  SUBROUTINE calc_statistics_vorticity()
!
!  written by Reni
!------------------------------------------
   SUBROUTINE calc_statistics_vorticity(stat_flagv)
 
     USE global_parameters
     USE flow_parameters
     USE grid_arrays,ONLY : xc,yc,zc
     USE flow_arrays
     USE stat_vort_arrays
     USE multiuse_arrays
 
     IMPLICIT NONE
 
!... Parameters
     INTEGER, INTENT(IN):: stat_flagv
 
!... Loop Variables
     INTEGER           :: i,j,k
 
!... Local Variables
     REAL(KIND=CGREAL) :: rStatCtrv, rStatCtrSqrv
     REAL(KIND=CGREAL) :: rStatCtrInvv, rStatCtrSqrInvv
 
     call vorticity()
 
     SELECT CASE(stat_flagv)
       CASE (0)
         DO k = 1,nzc
         DO j = 1,nyc
         DO i = 1,nxc
           oxAv(i,j,k)   =  oxAv(i,j,k) + nlu(i,j,k)
           oyAv(i,j,k)   =  oyAv(i,j,k) + nlv(i,j,k)
           ozAv(i,j,k)   =  ozAv(i,j,k) + nlw(i,j,k)
           oxoxAv(i,j,k)  =  oxoxAv(i,j,k)+ nlu(i,j,k)*nlu(i,j,k)
           oyoyAv(i,j,k)  =  oyoyAv(i,j,k)+ nlv(i,j,k)*nlv(i,j,k)
           ozozAv(i,j,k)  =  ozozAv(i,j,k)+ nlw(i,j,k)*nlw(i,j,k)
         ENDDO ! i 
         ENDDO ! j 
         ENDDO ! k
      
       CASE (1)
         SELECT CASE(indexStatVort)
         CASE(1)
           OPEN(ifuStatOut,FILE='stat_vort_out1.dat',FORM='UNFORMATTED')
           indexStatVort = 2 
         CASE(2)
           OPEN(ifuStatOut,FILE='stat_vort_out2.dat',FORM='UNFORMATTED')
           indexStatVort = 1
         END SELECT ! indexStatVort

         rStatCtrv    = REAL(statCtrv,KIND=CGREAL)
         rStatCtrSqrv = rStatCtrv**2
         rStatCtrInvv = oned/rstatCtrv
         rStatCtrSqrInvv = oned/rstatCtrSqrv
          
         WRITE(ifuStatOut)ntime-ntime_start+1,statCtrv,dt,         &
                         oxAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInvv,  &
                         oyAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInvv,  &
                         ozAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInvv,  &
                        oxoxAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInvv, &
                        oyoyAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInvv, &
                        ozozAv(0:nx+1,0:ny+1,0:nz+1)*rStatCtrInvv

         close (ifuStatOut)

         OPEN(ifuStatPlot,FILE='stat_vort_plot.dat')
 
         WRITE(ifuStatPlot,*)'VARIABLES="X","Y","Z","OXAV","OYAV","OZAV","OXpOXpAV","OYpOYpAV","OZpOZpAV"'
         WRITE(ifuStatPlot,*)'ZONE F=POINT, I=',nxc,',J=',nyc,'K=',nzc
         DO k = 1,nzc
         DO j = 1,nyc
         DO i = 1,nxc
           WRITE(ifuStatPlot,102) xc(i),yc(j),zc(k),                                        &
                        oxAv(i,j,k) *rStatCtrInvv,oyAv(i,j,k)*rStatCtrInvv                  &
                       ,ozAv(i,j,k) *rStatCtrInvv                                           &
                       ,oxoxAv(i,j,k)*rStatCtrInvv-oxAv(i,j,k)*oxAv(i,j,k)*rStatCtrSqrInvv  &
                       ,oyoyAv(i,j,k)*rStatCtrInvv-oyAv(i,j,k)*oyAv(i,j,k)*rStatCtrSqrInvv  &
                       ,ozozAv(i,j,k)*rStatCtrInvv-ozAv(i,j,k)*ozAv(i,j,k)*rStatCtrSqrInvv
         ENDDO ! i       
         ENDDO ! i       
         ENDDO ! k
          
         CLOSE(ifuStatPlot)
     END SELECT ! stat_flagv 
          
102  FORMAT(9(2x,e14.7))
          
          
   END SUBROUTINE calc_statistics_vorticity 
!------------------------------------------    
 
