!-------------------------------------------------------------
  SUBROUTINE GCM_set_freshcell()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    
!   USE GCM_ModInterfaces, ONLY : GCM_calc_bodyIntercept2D
    
    IMPLICIT NONE

!... loop variables

    INTEGER :: i,iBody,iRow,j,k,m,n

!... local variables

    INTEGER :: iFC,jFC,kFC, mM,mP,nbdr,nMarker
    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax

    INTEGER :: info

    REAL(KIND=CGREAL) :: cosTheta,dsIntercept,dsInterceptInv, &
                         probeLengthFresh,rCond,sinTheta,     &
                         xBI,yBI,zBI,xC1,xC2,xC3,xFC,yFC,zFC,xFIP,yFIP,zFIP,   &
                         slopeX, slopeY, slopeZ,              &
                         xBITangFresh, yBITangFresh,          &
                         xBINormFresh, yBINormFresh,sumcoeff

!*****************************************************************************************

! get the number of fresh cells

      nbdr = 0
      DO k = 1, nz-1
      DO j = 1, ny-1
      DO i = 1, nx-1
        IF ( fresh_Cell(i,j,k) == 1 )   THEN
           nbdr = nbdr + 1
        ENDIF ! fresh_Cell
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k
      
      nFresh = nbdr         ! total number of fresh cells
      PRINT*,'nFresh = ',nFresh

! Deallocate arrays pertinent to Fresh Cells 

      IF ( ntime >= ntime_start+1 ) THEN
        CALL GCM_DeallocateFreshCellArrays
      ENDIF ! ntime

! Allocate Arrays pertinent to Fresh Cells

      CALL GCM_AllocateFreshCellArrays()

! Set appropriate values for iGhost and jGhost by doing search

      nbdr = 0
      DO k = 1, nz-1
      DO j = 1, ny-1
      DO i = 1, nx-1
        IF ( fresh_Cell(i,j,k) == 1 )   THEN
           nbdr         = nbdr + 1
           iFresh(nbdr) = i
           jFresh(nbdr) = j
           kFresh(nbdr) = k
        ENDIF ! fresh_Cell
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

! Find marker point closest to each fresh  node
! and compute normal at closest marker point

      DO n = 1, nFresh

        iFC = iFresh(n)
        jFC = jFresh(n)
        kFC = kFresh(n)

        xFC = xc(iFC)
        yFC = yc(jFC)
        zFC = zc(kFC)

        CALL GCM_Calc_BodyIntercept_Unstruc( iFC, jFC, kFC, xFC, yFC, zFC,  xBI, yBI, zBI,&
                                       closestElementFresh(n)                       )

!-- Extract coordinates of Body Intercept

         xBodyInterceptFresh(n) = xBI
         yBodyInterceptFresh(n) = yBI
         zBodyInterceptFresh(n) = zBI

!-- Get length of normal

         dsIntercept    = SQRT( (yFC-yBI)**2 + (xFC-xBI)**2 + (zFC-zBI)**2 )
         dsInterceptInv = oned/dsIntercept

!-- Now compute location of Image Point
!-- for fresh cells, the image point only needed to determine surrounding nodes.

         slopeX = (xFC - xBI)*dsInterceptInv
         slopeY = (yFC - yBI)*dsInterceptInv
         slopeZ = (zFC - zBI)*dsInterceptInv

         probeLengthFresh = dsIntercept*probeLengthNormalized

         xFIP = xBI + slopeX*probeLengthFresh
         yFIP = yBI + slopeY*probeLengthFresh
         zFIP = zBI + slopeZ*probeLengthFresh

! Find the lower left grid point to Image Point in Physical domain
         iMin = iFresh(n)-1
         iMax = iFresh(n)+1
         iMin = MAX(iMin,1)
         iMax = MIN(iMax,nx-1)

         DO i = iMin,iMax
           IF ( xc(i) <= xFIP .AND. xc(i+1) > xFIP ) iFreshCellIndex(n) = i
         ENDDO ! i

         jMin = jFresh(n)-1
         jMax = jFresh(n)+1
         jMin = MAX(jMin,1)
         jMax = MIN(jMax,ny-1)

         DO j = jMin,jMax
           IF ( yc(j) <= yFIP .AND. yc(j+1) > yFIP ) jFreshCellIndex(n) = j
         ENDDO ! j

         kMin = kFresh(n)-1
         kMax = kFresh(n)+1
         kMin = MAX(kMin,1)
         kMax = MIN(kMax,nz-1)

         DO k = kMin,kMax
           IF ( zc(k) <= zFIP .AND. zc(k+1) > zFIP ) kFreshCellIndex(n) = k
         ENDDO ! k

! Test Closest Marker with 9-point Stencil

         IF ( (iFreshCellIndex(n)-iFresh(n)) > +0 .OR.  &
              (iFreshCellIndex(n)-iFresh(n)) < -1 .OR.  &
              (jFreshCellIndex(n)-jFresh(n)) > +0 .OR.  &
              (jFreshCellIndex(n)-jFresh(n)) < -1 .OR.  &
              (kFreshCellIndex(n)-kFresh(n)) > +0 .OR.  &
              (kFreshCellIndex(n)-kFresh(n)) < -1       )  THEN

            PRINT*,'IMAGE POINT FOR FRESH CELL IS NOT INSIDE 9-POINT STENCIL'
            PRINT*,'iBody,n,iFresh,jFresh,kFresh,iFreshCellIndex,jFreshCellIndex,kFreshCellIndex'
            PRINT*,iBody,n,iFresh(n),jFresh(n),kFresh(n),iFreshCellIndex(n),jFreshCellIndex(n),kFreshCellIndex(n)
            PRINT*,'dsIntercept,probeLengthFresh'
            PRINT*,dsIntercept,probeLengthFresh
            PRINT*,'slopeX, slopeY,slopeZ'
            PRINT*,slopeX, slopeY,slopeZ
            PRINT*,'xFIP,yFIP,zFIP'
            PRINT*,xFIP,yFIP,zFIP
            STOP
         ENDIF ! iFreshCellIndex

! Perform bilinear interpolation

!   |-------|---------|----/--|-------|--         N : Nth fresh node
!   |   ii  |  iii    |   *   |       |           * : markers
!   |   0...|....O    |  /.   |   .   |           O : other nodes used in bilinear interpolation
!   |   .   |    .    | *     |       |           + : probe tip (Image Point)
!   |---.---|-+---.---|/------|-------|--         @ : body intercept
!   |   .   |  \  .   *       |       |
!   |   0.  |   \N . /|       |   .   |
!   |   i ....   \ .* |       |       |
!   |-------| .... @--|-------|-------|--
!                 /
! interpolant      U = a X X  + b X  + c X + d
!                         1 2      1      2
!
!   Van Matrix For Dirichlet conditions at Intersection Point (@)
!
!         [  X X     X     X   1  ]  [   ]     [     ]
!      i  [   1 2     1     2     ]  [ a ]     [ U   ]
!         [                       ]  [   ]     [  i  ]
!         [  X X     X     X   1  ]  [   ]     [     ]
!      ii [   1 2     1     2     ]  [ b ]     [ U   ]
!         [                       ]  [   ]  =  [  ii ]
!     iii [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ c ]     [ U   ]
!         [                       ]  [   ]     [  iii]
!      N  [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ d ]     [ U   ]
!         [                       ]  [   ]     [  @  ]
!
        DO iRow= 1, iRowMax
          i  = iFreshCellIndex(n) + incI(iRow)
          j  = jFreshCellIndex(n) + incJ(iRow)
          k  = kFreshCellIndex(n) + incK(iRow)

          xC1 = xc(i)
          xC2 = yc(j)
          xC3 = zc(k)

!-- Construct Vandermonde Matrices

!--- Dirichlet conditions for velocity field

          vanMatrixD(iRow,1) = xC1*xC2*xC3
          vanMatrixD(iRow,2) = xC1*xC2
          vanMatrixD(iRow,3) = xC1*xC3
          vanMatrixD(iRow,4) = xC2*xC3
          vanMatrixD(iRow,5) = xC1
          vanMatrixD(iRow,6) = xC2
          vanMatrixD(iRow,7) = xC3
          vanMatrixD(iRow,8) = oned

!-- Correct For Fresh node part of cell formation, switch to Body Intercept point

          IF ( i==iFC .AND. j == jFC .AND. k == kFC) THEN
            xC1 = xBI
            xC2 = yBI
            xC3 = zBI

            vanMatrixD(iRow,1) = xC1*xC2*xC3
            vanMatrixD(iRow,2) = xC1*xC2
            vanMatrixD(iRow,3) = xC1*xC3
            vanMatrixD(iRow,4) = xC2*xC3
            vanMatrixD(iRow,5) = xC1
            vanMatrixD(iRow,6) = xC2
            vanMatrixD(iRow,7) = xC3
            vanMatrixD(iRow,8) = oned

          ENDIF ! i

        ENDDO ! iRow

! Compute inverse of Vandermonde Matrices

        CALL DGETRF(8, 8, vanMatrixD,8,iPvt, info)
        CALL DGETRI(8, vanMatrixD,8,iPvt,work, 8, info)

! Load Coeff-Matrices

        xFC = xc(iFC)
        yFC = yc(jFC)
        zFC = zc(kFC)
!RRRRRRRRRRRRRRRRRRRRRRR
        sumcoeff = zero
!RRRRRRRRRRRRRRRRRRRRRRR
        DO iRow = 1, iRowMax
          coeffGCMFreshD(iRow,n) =   vanMatrixD(1,iRow)*xFC*yFC*zFC  &
                                   + vanMatrixD(2,iRow)*xFC*yFC      &
                                   + vanMatrixD(3,iRow)*xFC*zFC      &
                                   + vanMatrixD(4,iRow)*yFC*zFC      &
                                   + vanMatrixD(5,iRow)*xFC          &
                                   + vanMatrixD(6,iRow)*yFC          &
                                   + vanMatrixD(7,iRow)*zFC          &
                                   + vanMatrixD(8,iRow)
         ENDDO ! iRow

      ENDDO ! n
      

  END SUBROUTINE GCM_set_freshcell 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  SUBROUTINE GCM_SetBodyInterceptValuesFresh()

! Compute velocity components at
! body intercept points corresponding to fresh cells.
!----------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE pressure_arrays
    USE grid_arrays
    USE GCM_arrays

!   USE GCM_ModInterfaces, ONLY : GCM_calc_BIVelocity2D

    IMPLICIT NONE

!... Loop variables

    INTEGER :: n

!... Local variables

    INTEGER           :: iFr,jFr,kFr
    REAL(KIND=CGREAL) :: XBIFr,YBIFr,zBIFr

!*****************************************************************************************
! Loop over all fresh cell points

! velocity field

    DO n=1,nFresh

      iFr     = iFresh(n)
      jFr     = jFresh(n)
      kFr     = kFresh(n)

      xBIFr =  xBodyInterceptFresh(n)
      yBIFr =  yBodyInterceptFresh(n)
      zBIFr =  zBodyInterceptFresh(n)

      CALL GCM_calc_BIVelocity_Unstruc( iFr, jFr, kFr, xBIFr, yBIFr, zBIFr, closestElementFresh(n), &
                                        uBodyInterceptFresh(n),                                     &
                                        vBodyInterceptFresh(n),                                     &
                                        wBodyInterceptFresh(n)           )
    ENDDO ! n

  END SUBROUTINE GCM_SetBodyInterceptValuesFresh

