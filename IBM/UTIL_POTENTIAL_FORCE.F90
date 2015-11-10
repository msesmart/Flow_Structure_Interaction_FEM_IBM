

   SUBROUTINE drag_lift_potential() 
!
!   Compute the Lift and Drag coefficients for n-Bodies in the flow
!
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE usr_module ,ONLY : scx,scy,scz,scmx,scmy,scmz !VEERA

    IMPLICIT NONE

    INTEGER              :: i,j,k
    INTEGER              :: ibody

    REAL(KIND=CGREAL)    :: cxp,cyp,czp
    REAL(KIND=CGREAL)    :: cmxp,cmyp,cmzp
    REAL(KIND=CGREAL)    :: cx,cy,cz
    REAL(KIND=CGREAL)    :: cmx,cmy,cmz
    REAL(KIND=CGREAL)    :: xp,yp,zp
    REAL(KIND=CGREAL)    :: dPhidX,dPhidY,dPhidZ
    REAL(KIND=CGREAL)    :: slipVelTang1,slipVelTang2
   
    
    CHARACTER*9          :: dragfile
    CHARACTER*25         :: indragfile
  
    DO iBody = 1, nBody
 
!---- Pressure Term
      
      cxp  = zero
      cyp  = zero
      czp  = zero
      cmxp = zero
      cmyp = zero
      cmzp = zero

      DO k = 1, nz-1
      DO j = 1, ny-1
      DO i = 1, nx-1
        IF ( bodyNum(i-1,j  ,k  ) == iBody .OR. &
             bodyNum(i+1,j  ,k  ) == iBody .OR. &
             bodyNum(i  ,j+1,k  ) == iBody .OR. & 
             bodyNum(i  ,j-1,k  ) == iBody .OR. &
             bodyNum(i  ,j  ,k+1) == iBody .OR. &
             bodyNum(i  ,j  ,k-1) == iBody      ) THEN      

          xp = -p(i,j,k)*ium(i,j,k)*iblank(i-1,j,k)*dy(j)*dz(k) &
	       +p(i,j,k)*iup(i,j,k)*iblank(i+1,j,k)*dy(j)*dz(k)
		  
	  yp = -p(i,j,k)*jum(i,j,k)*iblank(i,j-1,k)*dx(i)*dz(k) &
	       +p(i,j,k)*jup(i,j,k)*iblank(i,j+1,k)*dx(i)*dz(k)
		  
	  zp = -p(i,j,k)*kum(i,j,k)*iblank(i,j,k-1)*dx(i)*dy(j) &
	       +p(i,j,k)*kup(i,j,k)*iblank(i,j,k+1)*dx(i)*dy(j)

          cxp = cxp + xp
	  cyp = cyp + yp
	  czp = czp + zp

          cmxp = cmxp + ( yc(j) - ycent(iBody) )*zp 
          cmyp = cmyp - ( xc(i) - xcent(iBody) )*zp
          cmzp = cmzp + ( xc(i) - xcent(iBody) )*yp &
                      - ( yc(j) - ycent(iBody) )*xp

          IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
            cmxp = cmxp - ( zc(k) - zcent(iBody) )*yp
            cmyp = cmyp + ( zc(k) - zcent(iBody) )*xp 
          ENDIF

        ENDIF ! bodyNum	

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k
      
      scx(ibody)  = cxp          !VEERA ... Fed to Flow Induced Motion subroutine
      scy(ibody)  = cyp          ! 
      scz(ibody)  = czp          !
      scmx(ibody) = cmxp         !
      scmy(ibody) = cmyp         !
      scmz(ibody) = cmzp         !VEERA ... Fed to Flow Induced Motion subroutine
 
!--- Construct Components

      cxp = 2.0_CGREAL*cxp
      cyp = 2.0_CGREAL*cyp
      czp = 2.0_CGREAL*czp

      cmxp= 2.0_CGREAL*cmxp
      cmyp= 2.0_CGREAL*cmyp
      cmzp= 2.0_CGREAL*cmzp

      IF ( ndim == DIM_2D ) THEN
        cxp  = cxp/zout   
        cyp  = cyp/zout  
        cmzp = cmzp/zout   
      ENDIF

      cx  = cxp          ! for non-cylinders these are
      cy  = cyp          ! raw
      cz  = czp          ! forces
      cmx = cmxp         ! and moments
      cmy = cmyp         ! need to non-dmensionalize them
      cmz = cmzp         ! during post-processing
      
!--- Write File
       
      WRITE(ifuDragOut+iBody-1,121) time,cxp,cyp,czp,cmxp,cmyp,cmzp
      
   END DO ! ibody 

121 FORMAT(20(3x,1PE12.5))
          
   END SUBROUTINE  drag_lift_potential
!-------------------------------------------------------------------------------

   SUBROUTINE GCM_drag_lift_potential() 
!
!   Compute the Lift and Drag coefficients for n-Bodies in the flow
!
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    USE usr_module ,ONLY : scx,scy,scz,scmx,scmy,scmz !VEERA

    IMPLICIT NONE

    INTEGER           :: i,j,k
    INTEGER           :: ibody,ifudrag
    INTEGER           :: iG, jG, kG, indZ
    INTEGER           :: nG,m ,m1, m2, mPointsOrig, nTrielemMax
    INTEGER           :: ii,jj,kk,n,iRow
    INTEGER           :: iClosest,jClosest,kClosest

    REAL(KIND=CGREAL) :: uIP, vIP, wIP,        &
                         uBI, vBI, wBI,        &
                         dUt1Dn,dUt2Dn,        &
                         alphaX,alphaY,alphaZ, &
                         rectElemCentX,rectElemCentY,rectElemCentZ,rectElemArea, &
                         dist, distMin, xp, yp, zp
    REAL(KIND=CGREAL) :: dist_int, p_inv_dist_int
                          


    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:)   :: cxp,cyp,czp,     &
                                                      cmxp,cmyp,cmzp

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:) :: pTriElemCent,slipVelTang1,slipVelTang2
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:)   :: InterceptElemCentX,InterceptElemCentY    						

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:)   :: cxpSurf, cypSurf
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:) :: cxpTemp, cypTemp
        
    
    CHARACTER*9          :: dragfile
    CHARACTER*25         :: indragfile
    
    CHARACTER*16         :: Ptname1		!output file names  ..abel
    CHARACTER*16         :: Pbname1
   
    ifudrag = 155

! allocate local arrays 

    ALLOCATE(cxp(nBody))
    ALLOCATE(cyp(nBody))
    ALLOCATE(czp(nBody))
    
    ALLOCATE(cxpSurf(nPtsBodyMarkerOrig(nBody)))	!... abel
    ALLOCATE(cypSurf(nPtsBodyMarkerOrig(nBody)))
    
    ALLOCATE(cxpTemp(nPtsBodyMarkerOrig(nBody),nz))	!... abel
    ALLOCATE(cypTemp(nPtsBodyMarkerOrig(nBody),nz))
    
    ALLOCATE(InterceptElemCentX(nPtsBodyMarkerOrig(nBody)))
    ALLOCATE(InterceptElemCentY(nPtsBodyMarkerOrig(nBody)))
 
    ALLOCATE(cmxp(nBody))
    ALLOCATE(cmyp(nBody))
    ALLOCATE(cmzp(nBody))
    
    nTriElemMax = MAXVAL(totNumTriElem(:))
    ALLOCATE(pTriElemCent(nBody,nTriElemMax))
    ALLOCATE(slipVelTang1(nBody,nTriElemMax))
    ALLOCATE(slipVelTang2(nBody,nTriElemMax))
               
! initialize variables

    cxp = zero
    cyp = zero
    czp = zero
    
    cxpSurf = zero		!...abel
    cypSurf = zero
    
    cxpTemp = zero		!...abel
    cypTemp = zero
       
    InterceptElemCentX = zero	!...abel
    InterceptElemCentY = zero
    
    cmxp = zero
    cmyp = zero
    cmzp = zero
   
    DO iBody = 1, nBody

      SELECT CASE( ndim )

        CASE(DIM_2D)

          mPointsOrig = totNumTriElem(iBody)/2/(nz-1)
	
          DO indZ = 1,nz-1

           DO m=1,mPointsOrig   ! looping over original marker points

             m1 = 2*mPointsOrig*(indZ-1)+(2*m-1)      !counter results in odd values
			m2 = m1+1				      !counter results in even values
	     	
  ! compute pressure at centroid of rectangle formed by these two triangular elements
             rectElemCentX = half*(triElemCentX(iBody,m1) + triElemCentX(iBody,m2))
             rectElemCentY = half*(triElemCentY(iBody,m1) + triElemCentY(iBody,m2))
             rectElemCentZ = half*(triElemCentZ(iBody,m1) + triElemCentZ(iBody,m2))
	
			rectElemArea  = triElemArea(iBody,m1) + triElemArea(iBody,m2)
	     
  ! Storing body intercept coordinates 			...abel
			 InterceptElemCentX(m) = rectElemCentX		
			 InterceptElemCentY(m) = rectElemCentY		

	   
  ! find closest ghost node
             distMin = 1.0E8_CGREAL
             DO n=1,nGhost
               i = iGhost(n)
               j = jGhost(n)
               k = kGhost(n)

               dist = (xc(i)-rectElemCentX)**2 &
                     +(yc(j)-rectElemCentY)**2 &
                     +(zc(k)-rectElemCentZ)**2 

               IF ( dist <= distMin ) THEN
                  distMin = dist
                  nG  = n
                  iG  = i
                  jG  = j
                  kG  = k
	       ENDIF
             ENDDO

             IF (kG /= indZ) THEN
                PRINT*,' WARNING #################something wrong in lift-drag routine'
                PRINT*,nG,iG,jG,kG
             ENDIF

  ! search vicinity of closest ghost node to find nodes surrounding element
             DO i = iG-3, iG+3
               IF ( ( xc(i) <= rectElemCentX ) .AND. ( xc(i+1) > rectElemCentX ) ) ii = i
             ENDDO
             DO j = jG-3, jG+3
               IF ( ( yc(j) <= rectElemCentY ) .AND. ( yc(j+1) > rectElemCentY ) ) jj = j
             ENDDO
             kk = indz

! inverse distance sqr weighted interpolation
             p_inv_dist_int = zero
             dist_int       = zero

             DO j=0,1
             DO i=0,1
                dist = (xc(ii+i)-rectElemCentX)**2 &
                      +(yc(jj+j)-rectElemCentY)**2 &
                      +(zc(kk)  -rectElemCentZ)**2
                IF ( iblank(ii+i,jj+j,kk) == 0 ) THEN
                   dist_int       = dist_int + (oned/dist)
                   p_inv_dist_int = p_inv_dist_int + (oned/dist)*p(ii+i,jj+j,kk)
                ENDIF
             ENDDO
             ENDDO

             pTriElemCent(iBody,m) = p_inv_dist_int/dist_int

!---- Pressure Force Term
!   C_p = Integral (p n.ds)  (no negative sign since n acts into the body)

             xp =  pTriElemCent(iBody,m)*rectElemArea*triElemNormx(iBody,m1)
             yp =  pTriElemCent(iBody,m)*rectElemArea*triElemNormy(iBody,m1)
             zp =  pTriElemCent(iBody,m)*rectElemArea*triElemNormz(iBody,m1)
  	  	     
	     cxp(iBody) = cxp(iBody) + xp
             cyp(iBody) = cyp(iBody) + yp 
             czp(iBody) = czp(iBody) + zp
	     
	     cmxp(iBody) = cmxp(iBody) + ( yc(j) - ycent(iBody) )*zp
             cmyp(iBody) = cmyp(iBody) - ( xc(i) - xcent(iBody) )*zp
             cmzp(iBody) = cmzp(iBody) + ( xc(i) - xcent(iBody) )*yp - ( yc(j) - ycent(iBody) )*xp
  
  !~~~~~~~~~~~~~~ Storing Surface Components ~~~~~~~~~~~~~~~~~~~~~~~~~~	
	     cxpTemp(m,indZ) = 2.0_CGREAL*xp/rectElemArea / (nz-1)	!remove Area, want comp. only
    	     cypTemp(m,indZ) = 2.0_CGREAL*yp/rectElemArea / (nz-1) 	! values normalized
  !~~~~~~~~~~~~~~ Storing Surface Components ~~~~~~~~~~~~~~~~~~~~~~~~~~	

           ENDDO ! m - markerPts loop

         ENDDO  ! indZ
	 

  !~~~~~~~~~~~~~~ Summing up components along the span ~~~~~~~~~~~~~~~~~~~~~~~~~ .. abel	
	DO m = 1, mPointsOrig	
	  cxpSurf(m) = cxpTemp(m,nz-2) + cxpTemp(m,nz-1)
	  cypSurf(m) = cypTemp(m,nz-2) + cypTemp(m,nz-1)
	ENDDO
  !~~~~~~~~~~~~~~ Summing up components along the span ~~~~~~~~~~~~~~~~~~~~~~~~~ .. abel	


!******************************************************************************************
!*********** Printing terms along the surface --- modified by abel **************************
 
!IF ( MOD(ntime,nrestart) == 0 .OR. ntime==ntime_start+no_tsteps) THEN 
 !IF ( MOD(ntime,ndump) == 0 )    THEN			

    IF ( ntime >= 0 .AND. ntime .le. 9  )            &
        write(Ptname1,341)ntime
    IF ( ntime >= 10 .AND. ntime .le.99 )            &
        write(Ptname1,342)ntime
    IF ( ntime >= 100 .AND. ntime .le. 999 )         &
        write(Ptname1,343)ntime
    IF ( ntime >= 1000 .AND. ntime .le. 9999 )       &
        write(Ptname1,344)ntime
    IF ( ntime >= 10000 .AND. ntime .le. 99999 )     &
        write(Ptname1,345)ntime
    IF ( ntime >= 100000 .AND. ntime .le. 999999 )   &
        write(Ptname1,346)ntime
    IF ( ntime >= 1000000 .AND. ntime .le. 9999999 ) &
        write(Ptname1,347)ntime
341     format('Ptop.000000',i1)
342     format('Ptop.00000',i2)
343     format('Ptop.0000',i3)
344     format('Ptop.000',i4)
345     format('Ptop.00',i5)
346     format('Ptop.0',i6)
347     format('Ptop.',i7)


   IF ( ntime >= 0 .AND. ntime .le. 9  )            &
        write(Pbname1,441)ntime
   IF ( ntime >= 10 .AND. ntime .le.99 )            &
        write(Pbname1,442)ntime
   IF ( ntime >= 100 .AND. ntime .le. 999 )         &
        write(Pbname1,443)ntime
   IF ( ntime >= 1000 .AND. ntime .le. 9999 )       &
        write(Pbname1,444)ntime
   IF ( ntime >= 10000 .AND. ntime .le. 99999 )     &
        write(Pbname1,445)ntime
   IF ( ntime >= 100000 .AND. ntime .le. 999999 )   &
        write(Pbname1,446)ntime
   IF ( ntime >= 1000000 .AND. ntime .le. 9999999 ) &
        write(Pbname1,447)ntime
441     format('Pbot.000000',i1)
442     format('Pbot.00000',i2)
443     format('Pbot.0000',i3)
444     format('Pbot.000',i4)
445     format('Pbot.00',i5)
446     format('Pbot.0',i6)
447     format('Pbot.',i7)


    OPEN(UNIT=1,FILE=Ptname1,STATUS='UNKNOWN')		!top surface
    OPEN(UNIT=2,FILE=Pbname1,STATUS='UNKNOWN')		!bott surface

      WRITE(1,*)'VARIABLES = x/c y Pressure C_x_p C_y_p C_x_s C_y_s'  
      WRITE(2,*)'VARIABLES = x/c y Pressure C_x_p C_y_p C_x_s C_y_s'

      WRITE(*,*)'  '
      WRITE(*,*)'****** Saving pressure distribution  *******'   

   !-------Determine top & bottom surface by the sign of the triElemNormy --------------

      DO m=1,mPointsOrig     

	   WRITE(1,122) InterceptElemCentX(m), InterceptElemCentY(m), pTriElemCent(iBody,m)
	   WRITE(981,122) 180.0_CGREAL*ATAN2(InterceptElemCentY(m), InterceptElemCentX(m))/PI &
                           , pTriElemCent(iBody,m)

      ENDDO

    CLOSE(1)
    CLOSE(2)

!ENDIF

!******************************* END OF 2D CASE  ****************************************
!****************************************************************************************

        CASE(DIM_3D)

          IF(canonical_body_type(iBody) ==  ELLIPSOID) THEN
            PRINT*,'CANNOT COMPUTE SURFACE QUANTITIES FOR "ELLIPSOID"'
            PRINT*,'ABORTING RUN'
            STOP
	  ENDIF

          DO m=1,totNumTriElem(iBody)


  ! find closest ghost node
             distMin = 1.0E8_CGREAL
             DO n=1,nGhost
               i = iGhost(n)
               j = jGhost(n)
               k = kGhost(n)

               dist = (xc(i)-triElemCentX(iBody,m))**2 &
                     +(yc(j)-triElemCentY(iBody,m))**2 &
                     +(zc(k)-triElemCentZ(iBody,m))**2 

               IF ( dist <= distMin ) THEN
                  distMin = dist
                  nG  = n
                  iG  = i
                  jG  = j
                  kG  = k
	       ENDIF
             ENDDO

  ! search vicinity of closest ghost node to find nodes surrounding element
             DO i = iG-3, iG+3
               IF ( ( xc(i)-triElemCentX(iBody,m) )*(xc(i+1)-triElemCentX(iBody,m) ) <= zero ) ii = i
             ENDDO
             DO j = jG-3, jG+3
               IF ( ( yc(j)-triElemCentY(iBody,m) )*(yc(j+1)-triElemCentY(iBody,m) ) <= zero ) jj = j
             ENDDO
             DO k = kG-3, kG+3
               IF ( ( zc(k)-triElemCentZ(iBody,m) )*(zc(k+1)-triElemCentZ(iBody,m) ) <= zero ) kk = k
             ENDDO

! inverse distance sqr weighted interpolation
             p_inv_dist_int = zero
             dist_int       = zero                
             
             DO k=0,1
             DO j=0,1
             DO i=0,1
                dist = (xc(ii+i)-rectElemCentX)**2 &
                      +(yc(jj+j)-rectElemCentY)**2 &
                      +(zc(kk+k)-rectElemCentZ)**2
                IF ( iblank(ii+i,jj+j,kk+k) == 0 .OR. ghostCellMark(ii+i,jj+j,kk+k) == 1 ) THEN
                   dist_int       = dist_int + (oned/dist)
                   p_inv_dist_int = p_inv_dist_int + (oned/dist)*p(ii+i,jj+j,kk+k)
                ENDIF
             ENDDO
             ENDDO
             ENDDO

             pTriElemCent(iBody,m) = p_inv_dist_int/dist_int


!---- Pressure Force Term
!   C_p = Integral (p n.ds)  (no negative sign since n acts into the body)

             xp = pTriElemCent(iBody,m)*triElemArea(iBody,m)*triElemNormx(iBody,m)
             yp = pTriElemCent(iBody,m)*triElemArea(iBody,m)*triElemNormy(iBody,m)
             zp = pTriElemCent(iBody,m)*triElemArea(iBody,m)*triElemNormz(iBody,m)

             cxp(iBody) = cxp(iBody) + xp
             cyp(iBody) = cyp(iBody) + yp
             czp(iBody) = czp(iBody) + zp

             cmxp = cmxp + ( yc(j) - ycent(iBody) )*zp - ( zc(k) - zcent(iBody) )*yp
             cmyp = cmyp - ( xc(i) - xcent(iBody) )*zp + ( zc(k) - zcent(iBody) )*xp
             cmzp = cmzp + ( xc(i) - xcent(iBody) )*yp - ( yc(j) - ycent(iBody) )*xp

          ENDDO ! m

      END SELECT 
  
    ENDDO ! iBody

    DO iBody = 1, nBody

!--- Construct  Coefficients  
!    Note    Cf = F / (0.5 rho U^2 Af )      where Af is frontal area
!    Note    Cm = M / (0.5 rho U^2 Af L )    where L is lenght scale
!    Here we assume that  U = L = 1, 
!    Af = 1 x zout for 2D ;  Af = 1 for 3D
!    If this is not the case, then divide by appropriate values
	
	scx(iBody) = cxp(iBody)    !VEERA ... Fed to Flow Induced Motion subroutine
        scy(iBody) = cyp(iBody)    !VEERA ... Fed to Flow Induced Motion subroutine
        scz(iBody) = czp(iBody)    !VEERA ... Fed to Flow Induced Motion subroutine

        scmx(iBody) = cmxp(iBody) !VEERA ... Fed to Flow Induced Motion subroutine
        scmy(iBody) = cmyp(iBody) !VEERA ... Fed to Flow Induced Motion subroutine
        scmz(iBody) = cmzp(iBody) !VEERA ... Fed to Flow Induced Motion subroutine
      
        cxp(iBody) = 2.0_CGREAL*cxp(iBody)
        cyp(iBody) = 2.0_CGREAL*cyp(iBody)
        czp(iBody) = 2.0_CGREAL*czp(iBody)

        cmxp(iBody) = 2.0_CGREAL*cmxp(iBody)
        cmyp(iBody) = 2.0_CGREAL*cmyp(iBody)
        cmzp(iBody) = 2.0_CGREAL*cmzp(iBody)
                                            
        IF ( ndim == DIM_2D ) THEN
          cxp(iBody)  = cxp(iBody)/zout               
          cyp(iBody)  = cyp(iBody)/zout
          cmzp(iBody) = cmzp(iBody)/zout
        ENDIF

!--- Write File
      
        WRITE(ifuDragOut+iBody-1,121) time,                                    &
                                    cxp(iBody), &
                                    cyp(iBody), &
                                    czp(iBody), &
                                    cmxp(iBody),&
                                    cmyp(iBody),&
                                    cmzp(iBody),&
                                    surfArea(iBody)

    ENDDO ! iBody
  
! deallocate local arrays

    DEALLOCATE(cxp)
    DEALLOCATE(cyp)
    DEALLOCATE(czp)
    
    DEALLOCATE(cxpSurf)			!... abel
    DEALLOCATE(cypSurf)
    
    DEALLOCATE(cxpTemp)			!... abel
    DEALLOCATE(cypTemp)
   
    DEALLOCATE(InterceptElemCentX)	!... abel
    DEALLOCATE(InterceptElemCentY)
    
    DEALLOCATE(cmxp)
    DEALLOCATE(cmyp)
    DEALLOCATE(cmzp)
    
    DEALLOCATE(pTriElemCent)
    
121 FORMAT(20(3x,1PE12.5))
122 FORMAT(20(3x,1PE14.7))		!...abel
          
   END SUBROUTINE  GCM_drag_lift_potential
!--------------------------------------------------------------------------------
   SUBROUTINE gradient(qVar)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,0:nz+1),INTENT(IN)  :: qVar

    INTEGER              :: i,j,k
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb,pgx,pgy,pgz

! gradient =  d( )/dx i + d( )/dy j + d( )/dz k

      DO k = 1,nz-1
      DO j = 1,ny-1
      DO i = 1,nx-1

        pe = ( fx(i+1)*qVar(i+1,j,k) + (oned-fx(i+1))*qVar(i,j,k)   )*(1-iup(i,j,k)) &
             + qVar(i,j,k)*iup(i,j,k)

        pw = ( fx(i)  *qVar(i,j,k)   + (oned-fx(i))  *qVar(i-1,j,k) )*(1-ium(i,j,k)) &
             + qVar(i,j,k)*ium(i,j,k)

        pn = ( fy(j+1)*qVar(i,j+1,k) + (oned-fy(j+1))*qVar(i,j,k)   )*(1-jup(i,j,k)) &
             + qVar(i,j,k)*jup(i,j,k)

        ps = ( fy(j)  *qVar(i,j,k)   + (oned-fy(j))  *qVar(i,j-1,k) )*(1-jum(i,j,k)) &
             + qVar(i,j,k)*jum(i,j,k)

        pf = ( fz(k+1)*qVar(i,j,k+1) + (oned-fz(k+1))*qVar(i,j,k)   )*(1-kup(i,j,k)) &
             + qVar(i,j,k)*kup(i,j,k)

        pb = ( fz(k)  *qVar(i,j,k)   + (oned-fz(k))  *qVar(i,j,k-1) )*(1-kum(i,j,k)) &
             + qVar(i,j,k)*kum(i,j,k)

        nlu(i,j,k)= (pe-pw)*dxinv(i)*REAL(1-iblank(i,j,k),KIND=CGREAL)
        nlv(i,j,k)= (pn-ps)*dyinv(j)*REAL(1-iblank(i,j,k),KIND=CGREAL)
        nlw(i,j,k)= (pf-pb)*dzinv(k)*REAL(1-iblank(i,j,k),KIND=CGREAL)

      ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE gradient

