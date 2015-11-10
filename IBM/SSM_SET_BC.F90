SUBROUTINE SSM_set_bc_internal()

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE boundary_arrays
  USE grid_arrays
  USE gcm_arrays 
  USE unstructured_surface_arrays
  USE body_dynamics

  IMPLICIT NONE

    INTEGER               :: i,j,k,iBody,m,closest_marker,body_dist_min,count,cou,bcBody(nBody)
    INTEGER, DIMENSION(1) :: iclose,mVelLog
    REAL(KIND=CGREAL) :: dist_min,dist_min_body
    REAL(KIND=CGREAL) :: dist_sqr(nPtsMax)

    INTEGER :: iblankCellBound, iCellBound, nCellBound, anyWallPorous
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: iMapCellBound

    REAL(KIND=CGREAL) :: distX, distY, distZ, uT(nBody),vT(nBody),wT(nBody),velAmp(nBody)
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE   :: xyzCellBound
    REAL(KIND=CGREAL), DIMENSION(3)  :: Vec
    
    
    
    do k=1,nzc
    do j=1,nyc
    do i=1,nxc
        if(iblank(i,j,k)==1)then
            ium(i,j,k)=0
            iup(i,j,k)=0
            jum(i,j,k)=0
            jup(i,j,k)=0
            kum(i,j,k)=0
            kup(i,j,k)=0
        end if
    end do
    end do
    end do

    anyWallPorous = 0

    DO iBody = 1,nBody
      anyWallPorous = anyWallPorous + wall_type(iBody)
    ENDDO
 
    IF (boundary_motion /= MOVING_BOUNDARY .AND. &
        anyWallPorous == 0 ) THEN
       
      DO k=1,nzc    
      DO j=1,nyc   
      DO i=1,nxc    
        IF ( ium(i,j,k) == 1.or.iup(i,j,k)==1 ) THEN
          bcxu(i,j,k) = zero
          bcxv(i,j,k) = zero
          bcxw(i,j,k) = zero
        ENDIF
        IF ( jum(i,j,k) == 1.or.jup(i,j,k)==1 ) THEN
          bcyu(i,j,k) = zero
          bcyv(i,j,k) = zero
          bcyw(i,j,k) = zero
        ENDIF
        IF ( kum(i,j,k) == 1.or.kup(i,j,k)==1 ) THEN
          bczu(i,j,k) = zero
          bczv(i,j,k) = zero
          bczw(i,j,k) = zero
        ENDIF
      ENDDO ! i
      ENDDO ! j	
      ENDDO ! k

      RETURN

        ELSE


                
      bcBlank=0
      DO k=1,nzc
		  DO j=1,nyc
		  DO i=2,nxc-1  ! exclude outer boundary
              
            if(ium(i,j,k)==1.and.iblank(i-1,j,k)==0) iMarkm(i,j,k)=1
            if(iup(i,j,k)==1.and.iblank(i+1,j,k)==0) iMarkp(i,j,k)=1
            
            if(bcBlank(i,j,k)==1) cycle
            count=0
            do iBody=1,nBody
                if(btest(conflictCell(i,j,k),iBody-1))then
                    count=count+1
                    bcBody(count)=iBody
                end if
            end do
        
 
		    IF(ium(i,j,k) == 1.or.iup(i,j,k)==1) THEN
                
                if(ium(i,j,k)+iup(i,j,k)==2)then
                    Vec =(/xc(i), yc(j), zc(k)/)
                    !if(iblank(i-1,j,k)==0.and.iblank(i+1,j,k)==0)then
                    !    iMarkm(i,j,k)=1
                    !    iMarkp(i,j,k)=1
                    !    iMarkp(i-1,j,k)=1
                    !    iMarkm(i+1,j,k)=1
                    !end if
                else
                    if(ium(i,j,k)==1) Vec =(/x(i), yc(j), zc(k)/)
                    if(iup(i,j,k)==1) Vec =(/x(i+1), yc(j), zc(k)/)
                end if
            
                !if(ium(i,j,k)==1.and.iup(i,j,k)==1) write(5660,*) nTime,'i',i,j,k,count

                if(count==1)then
                    if(nBody_solid/=0.and.bcBody(1)==1)then
		                call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                    else
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(1))
                    end if
                    
                    bcxu(i,j,k) = uBodyMarker(body_dist_min,closest_marker)
		            bcxv(i,j,k) = vBodyMarker(body_dist_min,closest_marker)
		            bcxw(i,j,k) = wBodyMarker(body_dist_min,closest_marker)
                else
                    bodyNum(i,j,k)=0
                    !write(*,*) ium(i,j,k),iup(i,j,k)
                    
                    if(nBody_solid/=0.and.bcBody(1)==1.and.iblank(i-1,j,k)+iblank(i+1,j,k)>0)then
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                    else if(nBody_solid/=0.and.bcBody(1)==1)then
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(2),bcBody(count))
                    else
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(count))
                    end if
                    
                    !if(nBody_solid/=0.and.bcBody(1)==1.and.iblank(i-1,j,k)+iblank(i+1,j,k)>0)then
                    !    call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                    !else
                    !    call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(count))
                    !end if
                    
                    bcxu(i,j,k) = uBodyMarker(body_dist_min,closest_marker)
		            bcxv(i,j,k) = vBodyMarker(body_dist_min,closest_marker)
		            bcxw(i,j,k) = wBodyMarker(body_dist_min,closest_marker)
                    
                    !if(nBody_solid/=0.and.bcBody(1)==1)then
                    !    if(ium(i,j,k)==1.and.iblank(i-1,j,k)==0)then
                    !        bcBlank(i-1,j,k)=1
                    !        bcxu(i-1,j,k) = bcxu(i,j,k)
		                  !  bcxv(i-1,j,k) = bcxv(i,j,k)
		                  !  bcxw(i-1,j,k) = bcxw(i,j,k)
                    !    end if
                    !    if(iup(i,j,k)==1.and.iblank(i+1,j,k)==0)then
                    !        bcBlank(i+1,j,k)=1
                    !        bcxu(i+1,j,k) = bcxu(i,j,k)
		                  !  bcxv(i+1,j,k) = bcxv(i,j,k)
		                  !  bcxw(i+1,j,k) = bcxw(i,j,k)
                    !    end if
                    !    
                    !end if
                    
                    !if(ium(i,j,k)+iup(i,j,k)==2)then
                    !    if(iblank(i-1,j,k)==1.and.iblank(i+1,j,k)==0)then
                    !        bcBlank(i+1,j,k)=1
                    !        bcxu(i+1,j,k) = bcxu(i,j,k)
		                  !  bcxv(i+1,j,k) = bcxv(i,j,k)
		                  !  bcxw(i+1,j,k) = bcxw(i,j,k)
                    !    else if(iblank(i-1,j,k)==0.and.iblank(i+1,j,k)==1)then
                    !        bcBlank(i-1,j,k)=1
                    !        bcxu(i-1,j,k) = bcxu(i,j,k)
		                  !  bcxv(i-1,j,k) = bcxv(i,j,k)
		                  !  bcxw(i-1,j,k) = bcxw(i,j,k)
                    !    end if
                    !    
                    !end if
                    
              !      do cou=1,count
              !          if(nBody_solid/=0.and.bcBody(cou)==1)then
              !              call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
              !          else
              !              call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(cou),bcBody(cou))
              !          end if
              !          
              !          uT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          vT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          wT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          velAmp(cou)=uT(cou)*uT(cou)+vT(cou)*vT(cou)+wT(cou)*wT(cou)
              !      end do
              !      mVelLog=minloc(velAmp(1:count))
              !      bcxu(i,j,k) = uT(mVelLog(1))
		            !bcxv(i,j,k) = vT(mVelLog(1))
		            !bcxw(i,j,k) = wT(mVelLog(1))
                end if
		    END IF
		  END DO
		  END DO
		  END DO
 
      bcBlank=0
      DO k=1,nzc
		  DO j=2,nyc-1  ! exclude outer boundary
		  DO i=1,nxc  
              
            if(jum(i,j,k)==1.and.iblank(i,j-1,k)==0) jMarkm(i,j,k)=1
            if(jup(i,j,k)==1.and.iblank(i,j+1,k)==0) jMarkp(i,j,k)=1
            
            if(bcBlank(i,j,k)==1) cycle
            count=0
            do iBody=1,nBody
                if(btest(conflictCell(i,j,k),iBody-1))then
                    count=count+1
                    bcBody(count)=iBody
                end if
            end do
            
            
 
		    IF(jum(i,j,k) == 1.or.jup(i,j,k)==1) THEN
    

                if(jum(i,j,k)+jup(i,j,k)==2)then
                    Vec =(/xc(i), yc(j), zc(k)/)
                    !if(iblank(i,j-1,k)==0.and.iblank(i,j+1,k)==0)then
                    !    jMarkm(i,j,k)=1
                    !    jMarkp(i,j,k)=1
                    !    jMarkp(i,j-1,k)=1
                    !    jMarkm(i,j+1,k)=1
                    !end if
                else
                    if(jum(i,j,k)==1) Vec =(/xc(i), y(j), zc(k)/)
                    if(jup(i,j,k)==1) Vec =(/xc(i), y(j+1), zc(k)/)
                end if
                !if(jum(i,j,k)==1.and.jup(i,j,k)==1) write(5660,*) nTime,'j'

		      if(count==1)then
		            if(nBody_solid/=0.and.bcBody(1)==1)then
		                call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                    else
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(1))
                    end if
                    bcyu(i,j,k) = uBodyMarker(body_dist_min,closest_marker)
		            bcyv(i,j,k) = vBodyMarker(body_dist_min,closest_marker)
		            bcyw(i,j,k) = wBodyMarker(body_dist_min,closest_marker)
              else
                  bodyNum(i,j,k)=0
                  
                  if(nBody_solid/=0.and.bcBody(1)==1.and.iblank(i,j-1,k)+iblank(i,j+1,k)>0)then
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                    else if(nBody_solid/=0.and.bcBody(1)==1)then
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(2),bcBody(count))
                    else
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(count))
                    end if
                  
                  !if(nBody_solid/=0.and.bcBody(1)==1)then
                  !      call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                  !  else
                  !      call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(count))
                  !  end if
                    
                    bcyu(i,j,k) = uBodyMarker(body_dist_min,closest_marker)
		            bcyv(i,j,k) = vBodyMarker(body_dist_min,closest_marker)
		            bcyw(i,j,k) = wBodyMarker(body_dist_min,closest_marker)
                    
                    !if(nBody_solid/=0.and.bcBody(1)==1)then
                    !    if(jum(i,j,k)==1.and.iblank(i,j-1,k)==0)then
                    !        bcBlank(i,j-1,k)=1
                    !        bcyu(i,j-1,k) = bcyu(i,j,k)
		                  !  bcyv(i,j-1,k) = bcyv(i,j,k)
		                  !  bcyw(i,j-1,k) = bcyw(i,j,k)
                    !    end if
                    !    if(jup(i,j,k)==1.and.iblank(i,j+1,k)==0)then
                    !        bcBlank(i,j+1,k)=1
                    !        bcyu(i,j+1,k) = bcyu(i,j,k)
		                  !  bcyv(i,j+1,k) = bcyv(i,j,k)
		                  !  bcyw(i,j+1,k) = bcyw(i,j,k)
                    !    end if
                    !    
                    !end if
                    
                    !if(jum(i,j,k)+jup(i,j,k)==2)then
                    !    if(iblank(i,j-1,k)==1.and.iblank(i,j+1,k)==0)then
                    !        bcBlank(i,j+1,k)=1
                    !        bcyu(i,j+1,k) = bcyu(i,j,k)
		                  !  bcyv(i,j+1,k) = bcyv(i,j,k)
		                  !  bcyw(i,j+1,k) = bcyw(i,j,k)
                    !    else if(iblank(i,j-1,k)==0.and.iblank(i,j+1,k)==1)then
                    !        bcBlank(i,j-1,k)=1
                    !        bcyu(i,j-1,k) = bcyu(i,j,k)
		                  !  bcyv(i,j-1,k) = bcyv(i,j,k)
		                  !  bcyw(i,j-1,k) = bcyw(i,j,k)
                    !    end if
                    !    
                    !end if
              !      do cou=1,count
              !          if(nBody_solid/=0.and.bcBody(cou)==1)then
              !              call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
              !          else
              !              call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(cou),bcBody(cou))
              !          end if
              !          uT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          vT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          wT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          velAmp(cou)=uT(cou)*uT(cou)+vT(cou)*vT(cou)+wT(cou)*wT(cou)
              !      end do
              !      mVelLog=minloc(velAmp(1:count))
              !      bcyu(i,j,k) = uT(mVelLog(1))
		            !bcyv(i,j,k) = vT(mVelLog(1))
		            !bcyw(i,j,k) = wT(mVelLog(1))
                end if
   
		    END IF
 
		  END DO
		  END DO
		  END DO

      IF (nDim==DIM_3D) THEN
      bcBlank=0
      DO k=2,nzc-1  ! exclude outer boundary
		  DO j=1,nyc
		  DO i=1,nxc
              
            if(kum(i,j,k)==1.and.iblank(i,j,k-1)==0) kMarkm(i,j,k)=1
            if(kup(i,j,k)==1.and.iblank(i,j,k+1)==0) kMarkp(i,j,k)=1
              
            if(bcBlank(i,j,k)==1) cycle
            count=0
            do iBody=1,nBody
                if(btest(conflictCell(i,j,k),iBody-1))then
                    count=count+1
                    bcBody(count)=iBody
                end if
            end do
            
            
 
		    IF(kum(i,j,k) == 1.or.kup(i,j,k)==1) THEN

                if(kum(i,j,k)+kup(i,j,k)==2)then
                    Vec =(/xc(i), yc(j), zc(k)/)
                    !if(iblank(i,j,k-1)==0.and.iblank(i,j,k+1)==0)then
                    !    kMarkm(i,j,k)=1
                    !    kMarkp(i,j,k)=1
                    !    kMarkp(i,j,k-1)=1
                    !    kMarkm(i,j,k+1)=1
                    !end if
                else
                    if(kum(i,j,k)==1) Vec =(/xc(i), yc(j), z(k)/)
                    if(kup(i,j,k)==1) Vec =(/xc(i), yc(j), z(k+1)/)
                end if

                !if(kum(i,j,k)==1.and.kup(i,j,k)==1) write(5660,*) nTime,'k'

		      if(count==1)then
		            if(nBody_solid/=0.and.bcBody(1)==1)then
		                call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                    else
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(1))
                    end if
                    bczu(i,j,k) = uBodyMarker(body_dist_min,closest_marker)
		            bczv(i,j,k) = vBodyMarker(body_dist_min,closest_marker)
		            bczw(i,j,k) = wBodyMarker(body_dist_min,closest_marker)
              else
                  bodyNum(i,j,k)=0
                  
                  if(nBody_solid/=0.and.bcBody(1)==1.and.iblank(i,j,k-1)+iblank(i,j,k+1)>0)then
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                    else if(nBody_solid/=0.and.bcBody(1)==1)then
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(2),bcBody(count))
                    else
                        call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(count))
                    end if
                  
                  !if(nBody_solid/=0.and.bcBody(1)==1)then
                  !      call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
                  !  else
                  !      call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(1),bcBody(count))
                  !  end if
                    
                    bczu(i,j,k) = uBodyMarker(body_dist_min,closest_marker)
		            bczv(i,j,k) = vBodyMarker(body_dist_min,closest_marker)
		            bczw(i,j,k) = wBodyMarker(body_dist_min,closest_marker)
                    
                    !if(nBody_solid/=0.and.bcBody(1)==1)then
                    !    if(kum(i,j,k)==1.and.iblank(i,j,k-1)==0)then
                    !        bcBlank(i,j,k-1)=1
                    !        bczu(i,j,k-1) = bczu(i,j,k)
		                  !  bczv(i,j,k-1) = bczv(i,j,k)
		                  !  bczw(i,j,k-1) = bczw(i,j,k)
                    !    end if
                    !    if(kup(i,j,k)==1.and.iblank(i,j,k+1)==0)then
                    !        bcBlank(i,j,k+1)=1
                    !        bczu(i,j,k+1) = bczu(i,j,k)
		                  !  bczv(i,j,k+1) = bczv(i,j,k)
		                  !  bczw(i,j,k+1) = bczw(i,j,k)
                    !    end if
                    !    
                    !end if
                        
                    
                    !if(kum(i,j,k)+kup(i,j,k)==2)then
                    !    if(iblank(i,j,k-1)==1.and.iblank(i,j,k+1)==0)then
                    !        bcBlank(i,j,k+1)=1
                    !        bczu(i,j,k+1) = bczu(i,j,k)
		                  !  bczv(i,j,k+1) = bczv(i,j,k)
		                  !  bczw(i,j,k+1) = bczw(i,j,k)
                    !    else if(iblank(i,j,k-1)==0.and.iblank(i,j,k+1)==1)then
                    !        bcBlank(i,j,k-1)=1
                    !        bczu(i,j,k-1) = bczu(i,j,k)
		                  !  bczv(i,j,k-1) = bczv(i,j,k)
		                  !  bczw(i,j,k-1) = bczw(i,j,k)
                    !    end if
                    !    
                    !end if
              !      do cou=1,count
              !          if(nBody_solid/=0.and.bcBody(cou)==1)then
              !              call find_closest_element_modified(vec, body_dist_min, closest_marker,1,nBody_solid)
              !          else
              !              call find_closest_element_modified(vec, body_dist_min, closest_marker,bcBody(cou),bcBody(cou))
              !          end if
              !          uT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          vT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          wT(cou)=uBodyMarker(body_dist_min,closest_marker)
              !          velAmp(cou)=uT(cou)*uT(cou)+vT(cou)*vT(cou)+wT(cou)*wT(cou)
              !      end do
              !      mVelLog=minloc(velAmp(1:count))
              !      bczu(i,j,k) = uT(mVelLog(1))
		            !bczv(i,j,k) = vT(mVelLog(1))
		            !bczw(i,j,k) = wT(mVelLog(1))
                end if

		    END IF
 
		  END DO
		  END DO
		  END DO

      ENDIF

	  ENDIF ! boundary_motion
      
!        CALL write_dump_debug('bcyu',11,bcyu)
!        CALL write_dump_debug('bcyv',11,bcyv)

!        CALL write_dump_debug_body('uBodyMarker',niterFS,1,uBodyMarker)
!        CALL write_dump_debug_body('vBodyMarker',niterFS,1,vBodyMarker)

!          print *, 'min,maxval(uBodyMarker) =', minval(uBodyMarker), maxval(uBodyMarker)
!          print *, 'min,maxval(vBodyMarker) =', minval(vBodyMarker), maxval(vBodyMarker)
!          print *, 'min,maxval(bcxu) =', minval(bcxu), maxval(bcxu)
!          print *, 'min,maxval(bcyv) =', minval(bcyv), maxval(bcyv)
!          print *, 'min,maxval(bczw) =', minval(bczw), maxval(bczw)

        print *, 'max/min bcyv=',maxval(bcyv),minval(bcyv)

   END SUBROUTINE SSM_set_bc_internal

    
    
    
!SUBROUTINE SSM_set_bc_internal()
!
!USE global_parameters
!USE flow_parameters
!USE flow_arrays
!USE boundary_arrays
!USE grid_arrays
!USE gcm_arrays 
!USE unstructured_surface_arrays
!USE body_dynamics
!use operation
!
!IMPLICIT NONE
!
!INTEGER               :: i,j,k,ii,jj,kk,iBody,m,closest_marker(nBody),body_dist_min(1),count,conBody(nBody),iCon,bcBody(nBody),cBody,cMarker,t1,t2,t3
!INTEGER, DIMENSION(1) :: iclose
!REAL(KIND=CGREAL) :: dist_min(nBody),dist_min_body,uVelx,vVelx,wVelx,uVely,vVely,wVely,uVelz,vVelz,wVelz
!REAL(KIND=CGREAL) :: dist_sqr(nPtsMax),sumDist
!
!INTEGER :: iblankCellBound, iCellBound, nCellBound, anyWallPorous
!INTEGER, DIMENSION(:,:), ALLOCATABLE :: iMapCellBound
!
!REAL(KIND=CGREAL) :: distX, distY, distZ
!REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE   :: xyzCellBound
!
!REAL(KIND=CGREAL) :: velClosestMarker(nBody),dist(nBody),dist1(nBody),vec(3)
!logical :: solidBC
!
!anyWallPorous = 0
!
!uVelx=zero
!vVelx=zero
!wVelx=zero
!uVely=zero
!vVely=zero
!wVely=zero
!uVelz=zero
!vVelz=zero
!wVelz=zero
!
!DO iBody = 1,nBody
!    anyWallPorous = anyWallPorous + wall_type(iBody)
!ENDDO
! 
!IF (boundary_motion /= MOVING_BOUNDARY .AND. &
!    anyWallPorous == 0 ) THEN
!
!    DO k=1,nz-1    
!    DO j=1,ny-1    
!    DO i=1,nx-1    
!    IF ( (1-ium(i,j,k))*(1-iup(i,j,k)) == 0 ) THEN
!        bcxu(i,j,k) = zero
!        bcxv(i,j,k) = zero
!        bcxw(i,j,k) = zero
!    ENDIF
!    IF ( (1-jum(i,j,k))*(1-jup(i,j,k)) == 0 ) THEN
!        bcyu(i,j,k) = zero
!        bcyv(i,j,k) = zero
!        bcyw(i,j,k) = zero
!    ENDIF
!    IF ( (1-kum(i,j,k))*(1-kup(i,j,k)) == 0 ) THEN
!        bczu(i,j,k) = zero
!        bczv(i,j,k) = zero
!        bczw(i,j,k) = zero
!    ENDIF
!    ENDDO ! i
!    ENDDO ! j	
!    ENDDO ! k
!
!    if(channel_flow.and.nGate/=0)then
!    call identify_gates()
!    end if
!
!    RETURN
!
!ELSE
!
!    print *, 'SSM_SET_BC, max/minvBodyMarker:', maxval(vBodyMarker), minval(vBodyMarker)
!    print *, 'max/min(ium+iup)=', maxval(iup+ium),minval(iup+ium)
!
!    DO k=1,nz-1 
!	DO j=1,ny-1 
!	DO i=1,nx-1 
! 
!	!IF(ghostCellMemb(i,j,k) /= 0 .OR. ghostCellSolid(I,J,K)==1) THEN
!    IF(conflictCell(i,j,k)/=0) THEN
!        
!        !if(nTime==10.and.conflictCell(i,j,k)==3)then
!        !    write(*,*) nTime
!        !    write(*,*) u(i,j,k),v(i,j,k),w(i,j,k)
!        !end if
!        
!        velClosestMarker=zero
!        dist=zero
!        conBody=0
!        do iBody=1,nBody
!            dist_min(iBody)      = 1.0E10_CGREAL
!		    body_dist_min(1) = 0
!   
!            SELECT CASE(ndim)
!		        CASE(DIM_2D)
!  			        DO m=1,nPtsBodyMarker(iBody)
!			        distX = (xc(i)-xBodyMarker(iBody,m))
!			        distY = (yc(j)-yBodyMarker(iBody,m))
!			        dist_sqr(m) = distX**2 + distY**2
!			        ENDDO ! m
!		        CASE(DIM_3D)
!			        DO m=1,nPtsBodyMarker(iBody)
!			        distX = (xc(i)-xBodyMarker(iBody,m)) 
!			        distY = (yc(j)-yBodyMarker(iBody,m))
!			        distZ = (zc(k)-zBodyMarker(iBody,m))
!			        dist_sqr(m) = distX**2 + distY**2 + distZ**2
!			        ENDDO ! m
!		    END SELECT ! canonical_body_type
!
!		    dist_min_body  = MINVAL(dist_sqr(1:nPtsBodyMarker(iBody)))
!		    IF ( dist_min_body <= dist_min(iBody) ) THEN
!			    dist_min(iBody)      = dist_min_body
!			    iclose         = MINLOC(dist_sqr(1:nPtsBodyMarker(iBody)))
!			    closest_marker(iBody) = iclose(1)
!			    body_dist_min(1)  = iBody
!            ENDIF ! dist_min_body
!        end do
!        
!        iclose         = MINLOC(dist_min(1:nBody))
!!****************************************************************************************
!!****************************************************************************************
!!****************************************************************************************
!        
!        !count=0
!        !do iBody=1,nBody
!        !    if(btest(conflictCell(i,j,k),iBody-1))then
!        !        count=count+1
!        !        bcBody(count)=iBody
!        !    end if
!        !end do
!        !if(count==1)then
!        !    if(ium(i,j,k)==1.or.iup(i,j,k)==1)then
!        !        bcxu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        bcxv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        bcxw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        if(unstruc_surface_type(bcBody(1))==MEMBRANE.and.ium(i,j,k)==1) iMarkM(i,j,k)=-1
!        !        if(unstruc_surface_type(bcBody(1))==MEMBRANE.and.iup(i,j,k)==1) iMarkP(i,j,k)=-1
!        !    end if
!        !    if(jum(i,j,k)==1.or.jup(i,j,k)==1)then
!        !        bcyu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        bcyv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        bcyw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        if(unstruc_surface_type(bcBody(1))==MEMBRANE.and.jum(i,j,k)==1) jMarkM(i,j,k)=-1
!        !        if(unstruc_surface_type(bcBody(1))==MEMBRANE.and.jup(i,j,k)==1) jMarkP(i,j,k)=-1
!        !    end if
!        !    if(nDim == DIM_3D)then
!        !    if(kum(i,j,k)==1.or.kup(i,j,k)==1)then
!        !        bczu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        bczv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        bczw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!        !        if(unstruc_surface_type(bcBody(1))==MEMBRANE.and.kum(i,j,k)==1) kMarkM(i,j,k)=-1
!        !        if(unstruc_surface_type(bcBody(1))==MEMBRANE.and.kup(i,j,k)==1) kMarkP(i,j,k)=-1
!        !    end if
!        !    end if
!        !else
!        !    
!        !    if(unstruc_surface_type(bcBody(1))/=MEMBRANE)then
!        !        
!        !        do kk=k,k
!        !        do jj=j,j
!        !        do ii=i,i
!        !        
!        !        if(iblank(ii,jj,kk)==0) bodyNum(ii,jj,kk)=0
!        !        iblank_memb(ii,jj,kk)=0
!        !        ghostCellMemb(ii,jj,kk)=0
!        !        ghostCellMark(ii,jj,kk)=0
!        !        
!        !        if(iblank(ii,jj,kk)==0)then
!        !        !i-direction
!        !        if(ium(ii,jj,kk)==1)then
!        !            solidBC=.false.
!        !            do iBody=1,nBody
!        !                if(btest(conflictBCi(ii,jj,kk),iBody-1))then
!        !                    if(unstruc_surface_type(iBody)/=MEMBRANE)then
!        !                        solidBC=.true.
!        !                    end if
!        !                end if
!        !            end do
!        !            if(solidBC)then
!        !                iMarkM(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii-1,jj,kk)=0
!        !                vec=(/xc(ii),yc(jj),zc(kk)/)
!        !                call find_closest_element_modified(vec,cBody,cMarker,1,nBody_solid)
!        !                bcxu(ii,jj,kk)=uBodyMarker(cBody,cMarker)
!        !                bcxv(ii,jj,kk)=vBodyMarker(cBody,cMarker)
!        !                bcxw(ii,jj,kk)=wBodyMarker(cBody,cMarker)
!        !            else
!        !                iMarkM(ii,jj,kk)=0
!        !                iMarkP(ii-1,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii-1,jj,kk)=0
!        !                ium(ii,jj,kk)=0
!        !                iup(ii-1,jj,kk)=0
!        !                if(iup(ii,jj,kk)==0)then
!        !                    bcxu(ii,jj,kk)=zero
!        !                    bcxv(ii,jj,kk)=zero
!        !                    bcxw(ii,jj,kk)=zero
!        !                end if
!        !                if(ium(ii-1,jj,kk)==0)then
!        !                    bcxu(ii-1,jj,kk)=zero
!        !                    bcxv(ii-1,jj,kk)=zero
!        !                    bcxw(ii-1,jj,kk)=zero
!        !                end if
!        !            end if
!        !        end if
!        !        
!        !        if(iup(ii,jj,kk)==1)then
!        !            solidBC=.false.
!        !            do iBody=1,nBody
!        !                if(btest(conflictBCi(ii+1,jj,kk),iBody-1))then
!        !                    if(unstruc_surface_type(iBody)/=MEMBRANE)then
!        !                        solidBC=.true.
!        !                    end if
!        !                end if
!        !            end do
!        !            if(solidBC)then
!        !                iMarkP(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii+1,jj,kk)=0
!        !                vec=(/xc(ii),yc(jj),zc(kk)/)
!        !                call find_closest_element_modified(vec,cBody,cMarker,1,nBody_solid)
!        !                bcxu(ii,jj,kk)=uBodyMarker(cBody,cMarker)
!        !                bcxv(ii,jj,kk)=vBodyMarker(cBody,cMarker)
!        !                bcxw(ii,jj,kk)=wBodyMarker(cBody,cMarker)
!        !            else
!        !                iMarkP(ii,jj,kk)=0
!        !                iMarkM(ii+1,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii+1,jj,kk)=0
!        !                iup(ii,jj,kk)=0
!        !                ium(ii+1,jj,kk)=0
!        !                if(ium(ii,jj,kk)==0)then
!        !                    bcxu(ii,jj,kk)=zero
!        !                    bcxv(ii,jj,kk)=zero
!        !                    bcxw(ii,jj,kk)=zero
!        !                end if
!        !                if(iup(ii+1,jj,kk)==0)then
!        !                    bcxu(ii+1,jj,kk)=zero
!        !                    bcxv(ii+1,jj,kk)=zero
!        !                    bcxw(ii+1,jj,kk)=zero
!        !                end if
!        !            end if
!        !        end if
!        !        
!        !        !j-direction
!        !        if(jum(ii,jj,kk)==1)then
!        !            solidBC=.false.
!        !            do iBody=1,nBody
!        !                if(btest(conflictBCj(ii,jj,kk),iBody-1))then
!        !                    if(unstruc_surface_type(iBody)/=MEMBRANE)then
!        !                        solidBC=.true.
!        !                    end if
!        !                end if
!        !            end do
!        !            if(solidBC)then
!        !                jMarkM(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj-1,kk)=0
!        !                vec=(/xc(ii),yc(jj),zc(kk)/)
!        !                call find_closest_element_modified(vec,cBody,cMarker,1,nBody_solid)
!        !                bcyu(ii,jj,kk)=uBodyMarker(cBody,cMarker)
!        !                bcyv(ii,jj,kk)=vBodyMarker(cBody,cMarker)
!        !                bcyw(ii,jj,kk)=wBodyMarker(cBody,cMarker)
!        !            else
!        !                jMarkM(ii,jj,kk)=0
!        !                jMarkP(ii,jj-1,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj-1,kk)=0
!        !                jum(ii,jj,kk)=0
!        !                jup(ii,jj-1,kk)=0
!        !                if(jup(ii,jj,kk)==0)then
!        !                    bcyu(ii,jj,kk)=zero
!        !                    bcyv(ii,jj,kk)=zero
!        !                    bcyw(ii,jj,kk)=zero
!        !                end if
!        !                if(jum(ii,jj-1,kk)==0)then
!        !                    bcyu(ii,jj-1,kk)=zero
!        !                    bcyv(ii,jj-1,kk)=zero
!        !                    bcyw(ii,jj-1,kk)=zero
!        !                end if
!        !            end if
!        !        end if
!        !        
!        !        if(jup(ii,jj,kk)==1)then
!        !            solidBC=.false.
!        !            do iBody=1,nBody
!        !                if(btest(conflictBCj(ii,jj+1,kk),iBody-1))then
!        !                    if(unstruc_surface_type(iBody)/=MEMBRANE)then
!        !                        solidBC=.true.
!        !                    end if
!        !                end if
!        !            end do
!        !            if(solidBC)then
!        !                jMarkP(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj+1,kk)=0
!        !                vec=(/xc(ii),yc(jj),zc(kk)/)
!        !                call find_closest_element_modified(vec,cBody,cMarker,1,nBody_solid)
!        !                bcyu(ii,jj,kk)=uBodyMarker(cBody,cMarker)
!        !                bcyv(ii,jj,kk)=vBodyMarker(cBody,cMarker)
!        !                bcyw(ii,jj,kk)=wBodyMarker(cBody,cMarker)
!        !            else
!        !                jMarkP(ii,jj,kk)=0
!        !                jMarkM(ii,jj+1,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj+1,kk)=0
!        !                jup(ii,jj,kk)=0
!        !                jum(ii,jj+1,kk)=0
!        !                if(jum(ii,jj,kk)==0)then
!        !                    bcyu(ii,jj,kk)=zero
!        !                    bcyv(ii,jj,kk)=zero
!        !                    bcyw(ii,jj,kk)=zero
!        !                end if
!        !                if(jup(ii,jj+1,kk)==0)then
!        !                    bcyu(ii,jj+1,kk)=zero
!        !                    bcyv(ii,jj+1,kk)=zero
!        !                    bcyw(ii,jj+1,kk)=zero
!        !                end if
!        !            end if
!        !        end if
!        !        
!        !        
!        !        !k-direction
!        !        if(nDim == DIM_3D)then
!        !        if(kum(ii,jj,kk)==1)then
!        !            solidBC=.false.
!        !            do iBody=1,nBody
!        !                if(btest(conflictBCk(ii,jj,kk),iBody-1))then
!        !                    if(unstruc_surface_type(iBody)/=MEMBRANE)then
!        !                        solidBC=.true.
!        !                    end if
!        !                end if
!        !            end do
!        !            if(solidBC)then
!        !                kMarkM(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk-1)=0
!        !                vec=(/xc(ii),yc(jj),zc(kk)/)
!        !                call find_closest_element_modified(vec,cBody,cMarker,1,nBody_solid)
!        !                bczu(ii,jj,kk)=uBodyMarker(cBody,cMarker)
!        !                bczv(ii,jj,kk)=vBodyMarker(cBody,cMarker)
!        !                bczw(ii,jj,kk)=wBodyMarker(cBody,cMarker)
!        !            else
!        !                kMarkM(ii,jj,kk)=0
!        !                kMarkP(ii,jj,kk-1)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk-1)=0
!        !                kum(ii,jj,kk)=0
!        !                kup(ii,jj,kk-1)=0
!        !                if(kup(ii,jj,kk)==0)then
!        !                    bczu(ii,jj,kk)=zero
!        !                    bczv(ii,jj,kk)=zero
!        !                    bczw(ii,jj,kk)=zero
!        !                end if
!        !                if(kum(ii,jj,kk-1)==0)then
!        !                    bczu(ii,jj,kk-1)=zero
!        !                    bczv(ii,jj,kk-1)=zero
!        !                    bczw(ii,jj,kk-1)=zero
!        !                end if
!        !            end if
!        !        end if
!        !        
!        !        if(kup(ii,jj,kk)==1)then
!        !            solidBC=.false.
!        !            do iBody=1,nBody
!        !                if(btest(conflictBCk(ii,jj,kk+1),iBody-1))then
!        !                    if(unstruc_surface_type(iBody)/=MEMBRANE)then
!        !                        solidBC=.true.
!        !                    end if
!        !                end if
!        !            end do
!        !            if(solidBC)then
!        !                kMarkP(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk+1)=0
!        !                vec=(/xc(ii),yc(jj),zc(kk)/)
!        !                call find_closest_element_modified(vec,cBody,cMarker,1,nBody_solid)
!        !                bczu(ii,jj,kk)=uBodyMarker(cBody,cMarker)
!        !                bczv(ii,jj,kk)=vBodyMarker(cBody,cMarker)
!        !                bczw(ii,jj,kk)=wBodyMarker(cBody,cMarker)
!        !            else
!        !                kMarkP(ii,jj,kk)=0
!        !                kMarkM(ii,jj,kk+1)=0
!        !                !ghostCellMemb(ii,jj,kk)=0
!        !                !ghostCellMemb(ii,jj,kk+1)=0
!        !                kup(ii,jj,kk)=0
!        !                kum(ii,jj,kk+1)=0
!        !                if(kum(ii,jj,kk)==0)then
!        !                    bczu(ii,jj,kk)=zero
!        !                    bczv(ii,jj,kk)=zero
!        !                    bczw(ii,jj,kk)=zero
!        !                end if
!        !                if(kup(ii,jj,kk+1)==0)then
!        !                    bczu(ii,jj,kk+1)=zero
!        !                    bczv(ii,jj,kk+1)=zero
!        !                    bczw(ii,jj,kk+1)=zero
!        !                end if
!        !            end if
!        !        end if
!        !        end if
!        !        end if
!        !        end do
!        !        end do
!        !        end do
!        !    
!        !    
!        !    
!        !    end if
!        !    
!        !end if
!        
!        
!        !i-direction
!        if(ium(i,j,k)==1)then
!            count=0
!            do iBody=1,nBody
!                if(btest(conflictBCi(i,j,k),iBody-1))then
!                    count=count+1
!                    bcBody(count)=iBody
!                end if
!            end do
!            
!            if(count==1)then
!                if(iup(i,j,k)==1)then
!                    bcxu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                    bcxv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                    bcxw(i,j,k)=wBodyMarker(1,closest_marker(1))
!                else
!                    bcxu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bcxv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bcxw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                end if
!                if(unstruc_surface_type(bcBody(1))==MEMBRANE) iMarkM(i,j,k)=-1
!            else if(count>1)then
!                bcxu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                bcxv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                bcxw(i,j,k)=wBodyMarker(1,closest_marker(1))
!            end if
!        end if
!        
!        if(iup(i,j,k)==1)then
!            count=0
!            do iBody=1,nBody
!                if(btest(conflictBCi(i+1,j,k),iBody-1))then
!                    count=count+1
!                    bcBody(count)=iBody
!                end if
!            end do
!            
!            if(count==1)then
!                if(ium(i,j,k)==1)then
!                    bcxu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                    bcxv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                    bcxw(i,j,k)=wBodyMarker(1,closest_marker(1))
!                else
!                    bcxu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bcxv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bcxw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                end if
!                if(unstruc_surface_type(bcBody(1))==MEMBRANE) iMarkP(i,j,k)=-1
!            else if(count>1)then              
!                bcxu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                bcxv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                bcxw(i,j,k)=wBodyMarker(1,closest_marker(1))
!            end if
!        end if
!        
!        !j-direction
!        if(jum(i,j,k)==1)then
!            count=0
!            do iBody=1,nBody
!                if(btest(conflictBCj(i,j,k),iBody-1))then
!                    count=count+1
!                    bcBody(count)=iBody
!                end if
!            end do
!            
!            if(count==1)then
!                if(jup(i,j,k)==1)then
!                    bcyu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                    bcyv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                    bcyw(i,j,k)=wBodyMarker(1,closest_marker(1))
!                else
!                    bcyu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bcyv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bcyw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                end if
!                if(unstruc_surface_type(bcBody(1))==MEMBRANE) jMarkM(i,j,k)=-1
!            else if(count>1)then
!                bcyu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                bcyv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                bcyw(i,j,k)=wBodyMarker(1,closest_marker(1))
!            end if
!        end if
!        
!        if(jup(i,j,k)==1)then
!            count=0
!            do iBody=1,nBody
!                if(btest(conflictBCj(i,j+1,k),iBody-1))then
!                    count=count+1
!                    bcBody(count)=iBody
!                end if
!            end do
!            
!            if(count==1)then
!                if(jum(i,j,k)==1)then
!                    bcyu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                    bcyv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                    bcyw(i,j,k)=wBodyMarker(1,closest_marker(1))
!                else
!                    bcyu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bcyv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bcyw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                end if
!                if(unstruc_surface_type(bcBody(1))==MEMBRANE) jMarkP(i,j,k)=-1
!            else if(count>1)then
!                bcyu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                bcyv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                bcyw(i,j,k)=wBodyMarker(1,closest_marker(1))
!            end if
!        end if
!        
!        !k-direction
!        if(kum(i,j,k)==1)then
!            count=0
!            do iBody=1,nBody
!                if(btest(conflictBCk(i,j,k),iBody-1))then
!                    count=count+1
!                    bcBody(count)=iBody
!                end if
!            end do
!            
!            if(count==1)then
!                if(kup(i,j,k)==1)then
!                    bczu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                    bczv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                    bczw(i,j,k)=wBodyMarker(1,closest_marker(1))
!                else
!                    bczu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bczv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bczw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                end if
!                if(unstruc_surface_type(bcBody(1))==MEMBRANE) kMarkM(i,j,k)=-1
!            else if(count>1)then
!                bczu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                bczv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                bczw(i,j,k)=wBodyMarker(1,closest_marker(1))
!            end if
!        end if
!        
!        if(kup(i,j,k)==1)then
!            count=0
!            do iBody=1,nBody
!                if(btest(conflictBCk(i,j,k+1),iBody-1))then
!                    count=count+1
!                    bcBody(count)=iBody
!                end if
!            end do
!            
!            if(count==1)then
!                if(kum(i,j,k)==1)then
!                    bczu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                    bczv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                    bczw(i,j,k)=wBodyMarker(1,closest_marker(1))
!                else
!                    bczu(i,j,k)=uBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bczv(i,j,k)=vBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                    bczw(i,j,k)=wBodyMarker(bcBody(1),closest_marker(bcBody(1)))
!                end if
!                if(unstruc_surface_type(bcBody(1))==MEMBRANE) kMarkP(i,j,k)=-1
!            else if(count>1)then
!                bczu(i,j,k)=uBodyMarker(1,closest_marker(1))
!                bczv(i,j,k)=vBodyMarker(1,closest_marker(1))
!                bczw(i,j,k)=wBodyMarker(1,closest_marker(1))
!            end if
!        end if
!
!        !if(ium(i,j,k)==1.and.iup(i,j,k)==1)then    
!        !    iblank(i,j,k)=1
!        !    bcxu(i,j,k)=zero
!        !    bcxv(i,j,k)=zero
!        !    bcxw(i,j,k)=zero
!        !end if
!        !
!        !if(jum(i,j,k)==1.and.jup(i,j,k)==1)then    
!        !    iblank(i,j,k)=1
!        !    bcyu(i,j,k)=zero
!        !    bcyv(i,j,k)=zero
!        !    bcyw(i,j,k)=zero
!        !end if
!        !
!        !if(kum(i,j,k)==1.and.kup(i,j,k)==1)then    
!        !    iblank(i,j,k)=1
!        !    bczu(i,j,k)=zero
!        !    bczv(i,j,k)=zero
!        !    bczw(i,j,k)=zero
!        !end if
!        !
!        if(iblank(i,j,k)==1)then
!            ium(i,j,k)=0
!            iup(i,j,k)=0
!            jum(i,j,k)=0
!            jup(i,j,k)=0
!            kum(i,j,k)=0
!            kup(i,j,k)=0
!            bcxu(i,j,k)=zero
!            bcxv(i,j,k)=zero
!            bcxw(i,j,k)=zero
!            bcyu(i,j,k)=zero
!            bcyv(i,j,k)=zero
!            bcyw(i,j,k)=zero
!            bczu(i,j,k)=zero
!            bczv(i,j,k)=zero
!            bczw(i,j,k)=zero
!            iblank_memb(i,j,k)=0
!            ghostCellMemb(i,j,k)=0
!            ghostCellMark(i,j,k)=0
!        end if
!        
!        
!        !count=0
!        !do iBody=1,nBody
!        !    if(btest(conflictCell(i,j,k),iBody-1))then
!        !        count=count+1
!        !        dist(count)=sqrt(dist_min(iBody))
!        !        conBody(count)=ibody
!        !    end if
!        !end do
!        
!        !body_dist_min=MAXLOC(velClosestMarker)
!        !sumDist=sum(dist)
!        !uVel=zero
!        !vVel=zero
!        !wVel=zero
!        !do iCon=1,count
!        !    uVel(1)=uVel(1)+uBodyMarker(1,closest_marker(1))/real(count,CGREAL)
!        !    vVel(1)=vVel(1)+vBodyMarker(1,closest_marker(1))/real(count,CGREAL)
!        !    wVel(1)=wVel(1)+wBodyMarker(1,closest_marker(1))/real(count,CGREAL)
!        !    
!        !    uVel(2)=uVel(2)+uBodyMarker(2,closest_marker(2))/real(count,CGREAL)
!        !    vVel(2)=vVel(2)+vBodyMarker(2,closest_marker(2))/real(count,CGREAL)
!        !    wVel(2)=wVel(2)+wBodyMarker(2,closest_marker(2))/real(count,CGREAL)
!        !    
!        !end do
!        
!       
!		!bcxu(i,j,k) = uVelx*(ium(i,j,k)+iup(i,j,k))
!		!bcxv(i,j,k) = vVelx*(ium(i,j,k)+iup(i,j,k))
!		!bcxw(i,j,k) = wVelx*(ium(i,j,k)+iup(i,j,k))
!  ! 
!		!bcyu(i,j,k) = uVely*(jum(i,j,k)+jup(i,j,k))
!		!bcyv(i,j,k) = vVely*(jum(i,j,k)+jup(i,j,k))
!		!bcyw(i,j,k) = wVely*(jum(i,j,k)+jup(i,j,k))
!  !
!		!bczu(i,j,k) = uVelz*(kum(i,j,k)+kup(i,j,k))
!		!bczv(i,j,k) = vVelz*(kum(i,j,k)+kup(i,j,k))
!		!bczw(i,j,k) = wVelz*(kum(i,j,k)+kup(i,j,k))
!
!!print *,bcxu(i,j,k)
!	END IF
! 
!	END DO
!	END DO
!	END DO
!
!    ENDIF ! boundary_motion
!
!      
!!        CALL write_dump_debug('bcyu',11,bcyu)
!!        CALL write_dump_debug('bcyv',11,bcyv)
!
!!        CALL write_dump_debug_body('uBodyMarker',niterFS,1,uBodyMarker)
!!        CALL write_dump_debug_body('vBodyMarker',niterFS,1,vBodyMarker)
!
!!          print *, 'min,maxval(uBodyMarker) =', minval(uBodyMarker), maxval(uBodyMarker)
!!          print *, 'min,maxval(vBodyMarker) =', minval(vBodyMarker), maxval(vBodyMarker)
!!          print *, 'min,maxval(bcxu) =', minval(bcxu), maxval(bcxu)
!!          print *, 'min,maxval(bcyv) =', minval(bcyv), maxval(bcyv)
!!          print *, 'min,maxval(bczw) =', minval(bczw), maxval(bczw)
!
!    print *, 'max/min bcyv=',maxval(bcyv),minval(bcyv)
!    
!    !CALL write_dump_debug_i_2D('iMaP',1,iMarkP)            !yan_dbg
!    !CALL write_dump_debug_i_2D('iMaM',1,iMarkM)            !yan_dbg
!    !CALL write_dump_debug_i_2D('jMaP',1,jMarkP)            !yan_dbg
!    !CALL write_dump_debug_i_2D('jMaM',1,jMarkM)            !yan_dbg
!    !CALL write_dump_debug_i_2D('iiup',1,iup)            !yan_dbg
!    !CALL write_dump_debug_i_2D('iium',1,ium)            !yan_dbg
!    !CALL write_dump_debug_i_2D('jjup',1,jup)            !yan_dbg 
!    !CALL write_dump_debug_i_2D('jjum',1,jum)            !yan_dbg
!    !
!    !CALL write_dump_debug_f_2D('bcxu',1,iup)            !yan_dbg
!    !CALL write_dump_debug_f_2D('bcxv',1,ium)            !yan_dbg
!    !CALL write_dump_debug_f_2D('bcyu',1,jup)            !yan_dbg
!    !CALL write_dump_debug_f_2D('bcyv',1,jum)            !yan_dbg
!    
!    !write(555,*) nTime,bcxu(75,67,1),bcxv(75,67,1)
!    !write(556,*) nTime,bcxu(76,66,1),bcxv(76,66,1)
!    !CALL enforce_global_mass_consv()
!    
!    !if(nTime==11)then
!    !    t1=76
!    !    t2=67
!    !    t3=1
!    !    write(*,*) 'cell:',t1,t2,t3
!    !    write(*,*) ium(t1,t2,t3),iup(t1-1,t2,t3)
!    !    write(*,*) iup(t1,t2,t3),ium(t1+1,t2,t3)
!    !    
!    !    write(*,*) jum(t1,t2,t3),jup(t1,t2-1,t3)
!    !    write(*,*) jup(t1,t2,t3),jum(t1,t2+1,t3)
!    !    
!    !    write(*,*) kum(t1,t2,t3),kup(t1,t2,t3-1)
!    !    write(*,*) kup(t1,t2,t3),kum(t1,t2,t3+1)
!    !    
!    !    write(*,*) bcxu(t1,t2,t3),bcxv(t1,t2,t3),bcxw(t1,t2,t3)
!    !    write(*,*) bcyu(t1,t2,t3),bcyv(t1,t2,t3),bcyw(t1,t2,t3)
!    !    write(*,*) bczu(t1,t2,t3),bczv(t1,t2,t3),bczw(t1,t2,t3)
!    !    
!    !    write(*,*) ghostCellMemb(t1,t2,t3)
!    !    
!    !    
!    !    
!    !    
!    !    pause
!    !    write(*,*) jup(t1,t2,t3)
!    !end if
!
!    !bodyNum=0
!    
!END SUBROUTINE SSM_set_bc_internal

!----------------------------------------------------------------------------------

   subroutine identify_gates
   USE global_parameters
   USE flow_parameters
   USE flow_arrays
   USE boundary_arrays
   USE grid_arrays
   USE gcm_arrays 
   USE unstructured_surface_arrays
   USE body_dynamics

   IMPLICIT NONE

   integer :: i,j,k,m
   REAL(KIND=CGREAL) :: dist_min,distX,distY,distZ,dist
   integer :: mLabel

   gateTest=0

   do k=1,nzc
   do j=1,nyc
   do i=1,nxc
        if(iblank(i,j,k)==0.and.(iup(i,j,k)==1.or.ium(i,j,k)==1.or.&
                                 jup(i,j,k)==1.or.jum(i,j,k)==1.or.&
                                 kup(i,j,k)==1.or.kum(i,j,k)==1))then
            dist_min=1.0E10_CGREAL
            do m=1,nPtsBodymarker(1)
                distX=xc(i)-xbodymarker(1,m)
                distY=yc(j)-ybodymarker(1,m)
                distZ=zc(k)-zbodymarker(1,m)
                dist=sqrt(distX**twod+distY**twod+distZ**twod)
                if(dist<=dist_min)then
                    dist_min=dist
                    mLabel=m
                end if
            end do

            do m=1,nGate
                if(gateLabel(1,mLabel)==m)then

                gateTest(i,j,k)=m
                


                    select case(bcTypeGate(m))
                    case(3)
                        if(iup(i,j,k)==1.or.ium(i,j,k)==1)then
                            bcxu(i,j,k)=bcGateU(m)*sin(twod*pi*freqGate(m)*time)
                            bcxv(i,j,k)=bcGateV(m)*sin(twod*pi*freqGate(m)*time)
                            bcxw(i,j,k)=bcGateW(m)*sin(twod*pi*freqGate(m)*time)
                        end if
                        if(jup(i,j,k)==1.or.jum(i,j,k)==1)then
                            bcyu(i,j,k)=bcGateU(m)*sin(twod*pi*freqGate(m)*time)
                            bcyv(i,j,k)=bcGateV(m)*sin(twod*pi*freqGate(m)*time)
                            bcyw(i,j,k)=bcGateW(m)*sin(twod*pi*freqGate(m)*time)
                        end if
                        if(kup(i,j,k)==1.or.kum(i,j,k)==1)then
                            bczu(i,j,k)=bcGateU(m)*sin(twod*pi*freqGate(m)*time)
                            bczv(i,j,k)=bcGateV(m)*sin(twod*pi*freqGate(m)*time)
                            bczw(i,j,k)=bcGateW(m)*sin(twod*pi*freqGate(m)*time)
                        end if
                    case(2)
                        if(iup(i,j,k)==1.or.ium(i,j,k)==1)then
                            bcxu(i,j,k)=u(i,j,k)
                            bcxv(i,j,k)=v(i,j,k)
                            bcxw(i,j,k)=w(i,j,k)
                        end if
                        if(jup(i,j,k)==1.or.jum(i,j,k)==1)then
                            bcyu(i,j,k)=u(i,j,k)
                            bcyv(i,j,k)=v(i,j,k)
                            bcyw(i,j,k)=w(i,j,k)
                        end if
                        if(kup(i,j,k)==1.or.kum(i,j,k)==1)then
                            bczu(i,j,k)=u(i,j,k)
                            bczv(i,j,k)=v(i,j,k)
                            bczw(i,j,k)=w(i,j,k)
                        end if
                    case(1)
                        if(iup(i,j,k)==1.or.ium(i,j,k)==1)then
                            bcxu(i,j,k)=bcGateU(m)
                            bcxv(i,j,k)=bcGateV(m)
                            bcxw(i,j,k)=bcGateW(m)
                        end if
                        if(jup(i,j,k)==1.or.jum(i,j,k)==1)then
                            bcyu(i,j,k)=bcGateU(m)
                            bcyv(i,j,k)=bcGateV(m)
                            bcyw(i,j,k)=bcGateW(m)
                        end if
                        if(kup(i,j,k)==1.or.kum(i,j,k)==1)then
                            bczu(i,j,k)=bcGateU(m)
                            bczv(i,j,k)=bcGateV(m)
                            bczw(i,j,k)=bcGateW(m)
                        end if
                    end select
                end if

            end do
        end if

   



   end do
   end do
   end do

   !CALL enforce_global_mass_consv_channel()

    write(6678,*)  'VARIABLES="X","Y","Z","gateTest"'
    write(6678,*)  'ZONE F=POINT, I=32 , J=184 , K=240'
    do k=1,nzc
    do j=1,nyc
    do i=1,nxc
        write(6678,6678) xc(i),yc(j),zc(k),gateTest(i,j,k)
    end do
    end do
    end do
    close(6678)

    6678 format(3f16.8,i4)





   end subroutine identify_gates

   !--------------------------------------------------------------

      SUBROUTINE enforce_global_mass_consv_channel

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER             :: i,j,k,n
    REAL(KIND=CGREAL)   :: massflux,correction_vel,fluxBody

       outflow_area=zero
       do k=0,nz
       do j=0,ny
       do i=0,nx
            if(gateTest(i,j,k)==3)then 
                if(iup(i,j,k)==1) outflow_area=outflow_area+dy(j)*dz(k)
                if(ium(i,j,k)==1) outflow_area=outflow_area+dy(j)*dz(k)
                if(jup(i,j,k)==1) outflow_area=outflow_area+dx(i)*dz(k)
                if(jum(i,j,k)==1) outflow_area=outflow_area+dx(i)*dz(k)
                if(kup(i,j,k)==1) outflow_area=outflow_area+dx(i)*dy(j)
                if(kum(i,j,k)==1) outflow_area=outflow_area+dx(i)*dy(j)
            end if
       end do
       end do
       end do

    massflux = zero

       DO k=1,nzc
       DO j=1,nyc
       DO i=1,nxc
        IF(ghostCellMemb(i,j,k)==0) Then
         massflux = massflux +                            &
                  ( -bcxu(i,j,k)*ium(i,j,k)*dy(j)*dz(k)   &
                    +bcxu(i,j,k)*iup(i,j,k)*dy(j)*dz(k)   &
                    -bcyv(i,j,k)*jum(i,j,k)*dx(i)*dz(k)   &
                    +bcyv(i,j,k)*jup(i,j,k)*dx(i)*dz(k)   &
                    -bczw(i,j,k)*kum(i,j,k)*dx(i)*dy(j)   &
                    +bczw(i,j,k)*kup(i,j,k)*dx(i)*dy(j) ) &
                 *REAL(1-iblank(i,j,k),KIND=CGREAL)

         END IF
       ENDDO ! i
       ENDDO ! j
       ENDDO ! k

        correction_vel =-massflux/outflow_area  


!       DO k=0,nz
!       DO j=0,ny
!         IF (bcx1 == BC_TYPE_ZERO_GRADIENT) bcxu(1,j,k)   = bcxu(1,j,k)   - correction_vel
!         IF (bcx2 == BC_TYPE_ZERO_GRADIENT) bcxu(nxc,j,k) = bcxu(nxc,j,k) + correction_vel
!       ENDDO ! j
!       ENDDO ! k
!
!       DO k=0,nz
!       DO i=0,nx
!         IF (bcy1 == BC_TYPE_ZERO_GRADIENT) bcyv(i,1,k)   = bcyv(i,1,k)   - correction_vel
!         IF (bcy2 == BC_TYPE_ZERO_GRADIENT) bcyv(i,nyc,k) = bcyv(i,nyc,k) + correction_vel
!       ENDDO ! i
!       ENDDO ! k
!
!       DO j=0,ny
!       DO i=0,nx
!         IF (bcz1 == BC_TYPE_ZERO_GRADIENT) bczw(i,j,1)   = bczw(i,j,1)   - correction_vel
!         IF (bcz2 == BC_TYPE_ZERO_GRADIENT) bczw(i,j,nzc) = bczw(i,j,nzc) + correction_vel
!       ENDDO ! i
!       ENDDO ! j


       do k=0,nz
       do j=0,ny
       do i=0,nx
            if(gateTest(i,j,k)==3)then 
                if(iup(i,j,k)==1) bcxu(i,j,k)   = bcxu(i,j,k)   + correction_vel
                if(ium(i,j,k)==1) bcxu(i,j,k)   = bcxu(i,j,k)   - correction_vel
                if(jup(i,j,k)==1) bcxu(i,j,k)   = bcyv(i,j,k)   + correction_vel
                if(jum(i,j,k)==1) bcxu(i,j,k)   = bcyv(i,j,k)   - correction_vel
                if(kup(i,j,k)==1) bcxu(i,j,k)   = bczw(i,j,k)   + correction_vel
                if(kum(i,j,k)==1) bcxu(i,j,k)   = bczw(i,j,k)   - correction_vel
            end if
       end do
       end do
       end do


END SUBROUTINE enforce_global_mass_consv_channel
!---------------------------------------------------------------------
SUBROUTINE fill_cell

USE global_parameters
USE flow_parameters
USE flow_arrays
USE boundary_arrays
USE grid_arrays
USE gcm_arrays 
USE unstructured_surface_arrays
USE body_dynamics
use operation

IMPLICIT NONE

real(CGREAL) :: vec(3)
integer :: body_dist_min,closest_marker
integer :: i,j,k



do k=1,nzc
do j=1,nyc
do i=1,nxc
    if(iblank(i,j,k)==0)then
        if((iblank(i-1,j,k)==1.and.iup(i,j,k)==1).or.&
           (iblank(i+1,j,k)==1.and.ium(i,j,k)==1))then
            iblank(i,j,k)=1
            IF (fresh_cell(i,j,k) == 1) THEN
               fresh_cell(i,j,k) = 0
            end if
        
        end if
        
    end if
    
end do
end do
end do

do k=1,nzc
do j=1,nyc
do i=1,nxc
    if(iblank(i,j,k)==0)then
        if((iblank(i,j-1,k)==1.and.jup(i,j,k)==1).or.&
           (iblank(i,j+1,k)==1.and.jum(i,j,k)==1))then
            iblank(i,j,k)=1
            IF (fresh_cell(i,j,k) == 1) THEN
               fresh_cell(i,j,k) = 0
            end if
        
        end if
        
    end if
    
end do
end do
end do

if(ndim==dim_3D)then
do k=1,nzc
do j=1,nyc
do i=1,nxc
    if(iblank(i,j,k)==0)then
        if((iblank(i,j,k-1)==1.and.kup(i,j,k)==1).or.&
           (iblank(i,j,k+1)==1.and.kum(i,j,k)==1))then
            iblank(i,j,k)=1
            IF (fresh_cell(i,j,k) == 1) THEN
               fresh_cell(i,j,k) = 0
            end if
        
        end if
        
    end if
    
end do
end do
end do
    
    
end if






end subroutine fill_cell
