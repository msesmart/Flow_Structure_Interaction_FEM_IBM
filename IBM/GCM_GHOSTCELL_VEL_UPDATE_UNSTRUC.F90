!------------------------------------------------------------------------------
   SUBROUTINE GCM_GhostCell_Vel() 
!
!  update u*, v* and w* values for ghost cell
!
!
!  Ghost cell velocity satisfies the following equations
!
!  Integral [ DEL.(u*) dv ] =- Uo.nDS  with coeff of "dead" faces = 0
!
! [                                                               ]
! [ U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)   ] . tau
! [  gp                       ip      b                           ]   ---
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i


    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE flow_arrays
    USE grid_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... Loop variables

    INTEGER :: iBody, iRow, k, n, iterGC
    
!... Local variables

    INTEGER           :: i,j, ii,jj,kk,iP(3),nCol,cElement
    REAL(KIND=CGREAL) :: Wk(3),rC,rVelGhost(3),resVel,resVelMax,uPrev,vPrev,wPrev
    REAL(KIND=CGREAL) :: uIP,vIP,wIP,volCell
    REAL(KIND=CGREAL) :: bmx,bpx,bcx,bmy,bpy,bcy,bmz,bcz,bpz,ghostVelMatrix(3,3),fluxBody
    
! Iterate to correct interdependent ghost points

! Initialize values
       iterGC    = 0
       resVelMax = 1.0E10_CGREAL

       DO WHILE ((iterGC .LT. iterMax_ad) .AND. (resVelMax .GT. restol_ad))

         resVelMax = zero

         DO n = 1, nGhost

           i=iGhost(n)
           j=jGhost(n)
           k=kGhost(n)

           uIP = zero
           vIP = zero
           wIP = zero

           DO iRow = 1, iRowMax

            ii = iCellIndex(n) + incI(iRow)
            jj = jCellIndex(n) + incJ(iRow)
            kk = kCellIndex(n) + incK(iRow)


            IF ( ii /= i .OR. jj /= j .OR. kk /= k) THEN
              uIP = uIP + coeffGCMD(iRow,n)* u(ii,jj,kk)
              vIP = vIP + coeffGCMD(iRow,n)* v(ii,jj,kk)
              wIP = wIP + coeffGCMD(iRow,n)* w(ii,jj,kk)
            ELSE
              uIP = uIP + coeffGCMD(iRow,n)* uBodyIntercept(n)
              vIP = vIP + coeffGCMD(iRow,n)* vBodyIntercept(n)
              wIP = wIP + coeffGCMD(iRow,n)* wBodyIntercept(n)
            ENDIF ! ii

           ENDDO ! iRow

           uPrev = u(i,j,k)
           vPrev = v(i,j,k)
           wPrev = w(i,j,k)

           u(i,j,k) = zero
           v(i,j,k) = zero
           w(i,j,k) = zero

           u(i,j,k) =  uBodyIntercept(n)*bodyInterceptWeight - uIP*imagePointWeight
           v(i,j,k) =  vBodyIntercept(n)*bodyInterceptWeight - vIP*imagePointWeight
           w(i,j,k) =  wBodyIntercept(n)*bodyInterceptWeight - wIP*imagePointWeight

! Compute residual

           resVel = ABS( u(i,j,k)-uPrev ) + ABS( v(i,j,k)-vPrev )  &
                  + ABS( w(i,j,k)-wPrev )

           IF (resVel > resVelMax ) resVelMax = resVel

         ENDDO ! n

         IF ( MOD(ntime,nmonitor) == 0 ) &
         PRINT*,' Ghostcell Velocity Convergence: ',iterGC,resVelMax

         iterGC = iterGC + 1

       ENDDO ! iterGC

       IF ( iterGC .EQ. iterMax_ad .AND. resVelMax .GT. restol_ad ) THEN
        PRINT*,'GCM_vel_set_bc_internal for iBody :', iBody
        PRINT*,'   GhostCell u* did not converge in ',iterMax_ad,' iterations'
        PRINT*,'   Final residual = ',resVelMax
        STOP
       ELSE
        IF (MOD(ntime,nmonitor) == 0) THEN
!         PRINT*,'GhostCell u* convergence : k=',k,iterGC,resVelMax
        END IF ! ntime
       ENDIF

      CALL GCM_enforce_global_mass_consv()

   END SUBROUTINE  GCM_GhostCell_Vel

!----------------------------------------------------------------------
  SUBROUTINE GCM_vel_set_bc_internal()

!----------------------------------------------------------------------
! Compute Ghost point values at internal boundary for the velocity field
!----------------------------------------------------------------------
!
!     
!   U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)
!    gp                       ip      b
!
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i


    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE

!... Loop variables

    INTEGER :: iBody,iRow,k,n

!... Local variables

    INTEGER :: iG,ii,ind,jG,jj,kG,kk,iterGC
    
    REAL(KIND=CGREAL) :: cMR, uIP, vIP, wIP, &
                         resVel, resVelMax,  &
                         uPrev, vPrev, wPrev  

!*****************************************************************************************

! Compute body intercept velocity field

    CALL GCM_SetBodyInterceptValues()
    
! Loop over all immersed bodies

      iterGC    = 0
      resVelMax = 1.0E10_CGREAL
      
! Iterate to correct interdependent ghost points

      DO WHILE ((iterGC .LT. iterMax_ad) .AND. (resVelMax .GT. restol_ad))

        resVelMax = zero

        DO n = 1,nGhost

          iG = iGhost(n)
          jG = jGhost(n)
          kG = kGhost(n)

! Initialize values

            uIP = zero
            vIP = zero
            wIP = zero

! Compute velocity field at Image points

            DO iRow = 1, iRowMax
              ii = iCellIndex(n) + incI(iRow) 
              jj = jCellIndex(n) + incJ(iRow)
              kk = kCellIndex(n) + inck(iRow)

              IF ( ii /= iG .OR. jj /= jG .OR. kk /=kG ) THEN
                uIP = uIP + coeffGCMD(iRow,n)* u(ii,jj,kk)
                vIP = vIP + coeffGCMD(iRow,n)* v(ii,jj,kk) 
                wIP = wIP + coeffGCMD(iRow,n)* w(ii,jj,kk)
              ELSE
                uIP = uIP + coeffGCMD(iRow,n)* uBodyIntercept(n)
                vIP = vIP + coeffGCMD(iRow,n)* vBodyIntercept(n) 
                wIP = wIP + coeffGCMD(iRow,n)* wBodyIntercept(n)
              ENDIF ! ii          

            ENDDO ! iRow

! Load temporary values
          
            uPrev = u(iG,jG,kG)
            vPrev = v(iG,jG,kG)
            wPrev = w(iG,jG,kG)
          
! Apply Dirichlet conditions on Ghost Nodes 

            u(iG,jG,kG) = uBodyIntercept(n) *bodyInterceptWeight - uIP*imagePointWeight
            v(iG,jG,kG) = vBodyIntercept(n) *bodyInterceptWeight - vIP*imagePointWeight
            w(iG,jG,kG) = wBodyIntercept(n) *bodyInterceptWeight - wIP*imagePointWeight

! Compute residual

            resVel = ABS( u(iG,jG,kG)-uPrev ) &
                   + ABS( v(iG,jG,kG)-vPrev ) &
                   + ABS( w(iG,jG,kG)-wPrev )

            IF (resVel > resVelMax ) resVelMax = resVel
 
        ENDDO ! n    

        iterGC = iterGC + 1

    ENDDO ! iterGC

    IF ( iterGC .EQ. iterMax_ad .AND. resVelMax .GT. restol_ad ) THEN
      PRINT*,'GCM_vel_set_bc_internal for iBody :', iBody
      PRINT*,'   GhostCell Velocity did not converge in ',iterMax_ad,' iterations'
      PRINT*,'   Final residual = ',resVelMax
    ELSE
      IF (MOD(ntime,nmonitor) == 0) THEN
        PRINT*,'GhostCell Velocity convergence : ',iterGC,resVelMax
      END IF ! ntime
    ENDIF
    
  END SUBROUTINE GCM_vel_set_bc_internal      
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  SUBROUTINE GCM_p_set_bc_internal(pres,mx,my,mz)

!----------------------------------------------------------------------
! Compute Ghost point values at internal boundary for the velocity field
!----------------------------------------------------------------------
!
!     
!   P  =   P   
!    gp     ip 
!
!
! and
!        (          ) 
!   P  = (coeffGCMN ) P
!    ip  (         i)  i


    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER :: mx, my, mz
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT(INOUT) :: pres   ! upper bound changed to nx from nx+1

!... Loop variables

    INTEGER :: iBody,iRow,k,n

!... Local variables

    INTEGER :: iG,ii,ind,jG,jj,kG,kk,iterGC
    
    REAL(KIND=CGREAL) :: cMR, pIP, &
                         resPres, resPresMax,  &
                         pPrev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(KIND=CGREAL) :: uNormal,xBIN,yBIN,zBIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*****************************************************************************************

! Loop over all immersed bodies

    iterGC    = 0
    resPresMax = 1.0E10_CGREAL
      
! Iterate to correct interdependent ghost points

    DO WHILE ((iterGC .LT. iterMax_ad) .AND. (resPresMax .GT. restol_ad))

      resPresMax = zero

      DO n = 1,nGhost

        iG = iGhost(n)
        jG = jGhost(n)
        kG = kGhost(n)

        IF (flow_type == POTENTIAL_FLOW) THEN
          iBody = bodyNum(iG,jG,kG)
          CALL GCM_calc_BIVelocity_Unstruc( iG, jG, kG, xC(iG), yC(jG), zC(kG),  &
                                            closestElementGC(n),                &
                                            uBodyIntercept(n),                  &
                                            vBodyIntercept(n),                  &
                                            wBodyIntercept(n)                     )
            
          xBIN = triElemNormx(iBody,closestElementGC(n))
          yBIN = triElemNormy(iBody,closestElementGC(n))
          zBIN = triElemNormz(iBody,closestElementGC(n))
          uNormal = uBodyIntercept(n)*XBIN &
                  +vBodyIntercept(n)*YBIN &
                  +wBodyIntercept(n)*ZBIN 
        ELSE
	        uNormal = 0.0_CGREAL
        ENDIF

! Initialize values

        pIP = zero

! Compute velocity field at Image points

        DO iRow = 1, iRowMax
          ii = iCellIndex(n) + incI(iRow) 
          jj = jCellIndex(n) + incJ(iRow)
          kk = kCellIndex(n) + incK(iRow)

          IF ( ii == iG .AND. jj == jG .AND. kk == kG ) THEN
            pIP = pIP + coeffGCMN(iRow,n)*uNormal
          ELSE
            pIP = pIP + coeffGCMN(iRow,n)* pres(ii,jj,kk)
          ENDIF ! ii          
        ENDDO ! iRow

! Load temporary values
          
        pPrev = pres(iG,jG,kG)
          
! Apply Dirichlet conditions on Ghost Nodes 

        pres(iG,jG,kG) = pIP + probeLength(n)*uNormal
!          pres(iG,jG,kG) = pPrev + omega*(pIP-pPrev)  ! using overrelaxation

! Compute residual

        resPres = ABS( pres(iG,jG,kG) - pPrev )

        IF (resPres > resPresMax ) resPresMax = resPres
      
 
      ENDDO ! n    

      iterGC = iterGC + 1

    ENDDO ! iterGC

    IF ( iterGC .EQ. iterMax_ad .AND. resPresMax .GT. restol_ad ) THEN
      PRINT*,'GCM_p_set_bc_internal for iBody :', iBody
      PRINT*,'GhostCell Pressure did not converge in ',iterMax_ad,' iterations'
      PRINT*,'Final residual = ',resPresMax
    ELSE
!     IF (MOD(ntime,nmonitor) == 0) THEN
!       PRINT*,'GhostCell pressure convergence : ',iterGC,resPresMax
!     END IF ! ntime
    ENDIF
    
  END SUBROUTINE GCM_p_set_bc_internal      
!----------------------------------------------------------------------

  SUBROUTINE GCM_set_face_vel_body()

!----------------------------------------------------------------------
! Zeroing out face velocities for solid cells.
!  Also for "dead" faces of Ghost Cells
!----------------------------------------------------------------------
!
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE

!... Loop variables

    INTEGER :: i,j,k

    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc

      IF ( iblank(i,j,k) == 1 .AND. ghostCellMark(i,j,k) == 0 ) THEN
        face_u(i+1,j,k) = zero
        face_u(i,j,k)   = zero 
        face_v(i,j+1,k) = zero
        face_v(i,j,k)   = zero
        face_w(i,j,k+1) = zero
        face_w(i,j,k)   = zero
      ENDIF ! iblank

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
    
  END SUBROUTINE GCM_set_face_vel_body      
!----------------------------------------------------------------------
