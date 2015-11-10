!------------------------------------------------------------
  SUBROUTINE GCM_set_internal_boundary()
!
! Give a set of marker points, this subroutine computes
!   1) normal intercepts and associated geometrical info. from ghost nodes to body
!   2) Image point location
!   3) Weights in stencil for computing values at the image points
!

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

!... loop variables

    INTEGER :: i,iBody,iRow,j,k,m,n

!... local variables

    INTEGER :: iG, jG, kG, nbdr, iCIndx, jCIndx, kCIndx , iCIndxS, jCIndxS
    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
    INTEGER :: iM, iP, jM, jP, kM, kP
    INTEGER :: iRange, jRange, kRange
    INTEGER :: node1,node2,node3

    REAL(KIND=CGREAL) :: cosTheta,dsIntercept,sinTheta,               &
                         xBI,xBIN,xGC,xIP,xIPS,yBI,yBIN,yGC,yIP,yIPS,zBI,zBIN,zGC,zIP, &
                         minProbeLengthShear,slopeX, slopeY, slopeZ, maxDelta

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: coeffGCMDirc, coeffGCMNeum

!*****************************************************************************************

! Allocate arrays

    ALLOCATE( coeffGCMDirc(iRowMax) )
    ALLOCATE( coeffGCMNeum(iRowMax) )

!
! Mark boundary cells (ie. cells inside body which have 
! at least one neighbor in fluid)

    DO k = 1, nz-1
    DO j = 1, ny-1
    DO i = 1, nx-1

      iM = MAX(i-1,1)
      iP = MIN(i+1,nx-1)
      jM = MAX(j-1,1)
      jP = MIN(j+1,ny-1)
      kM = MAX(k-1,1)
      kP = MIN(k+1,nz-1)

      IF ( iblank(i,j,k)*iblank(iM,j,k) == 0  .OR.  &
           iblank(i,j,k)*iblank(iP,j,k) == 0  .OR.  &
           iblank(i,j,k)*iblank(i,jM,k) == 0  .OR.  &
           iblank(i,j,k)*iblank(i,jP,k) == 0  .OR.  &
           iblank(i,j,k)*iblank(i,j,kM) == 0  .OR.  &
           iblank(i,j,k)*iblank(i,j,kP) == 0        )  THEN
        IF (iblank(i,j,k) == 1) THEN
          ghostCellMark(i,j,k) = 1
        ENDIF
      ENDIF

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

    nbdr = 0
    DO k = 1, nz-1
    DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( ghostCellMark(i,j,k) == 1 )   THEN
        nbdr = nbdr + 1
      ENDIF ! ghostCellMark    
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k
      
    nGhost = nbdr         ! total number of ghost points

    IF ( MOD(ntime,nmonitor) == 0 .OR. ntime == 1) &
      PRINT*,'GCM_set_internal_boundary: nGhost = ',nGhost 

! Deallocate arrays pertinent to Ghost Cells 

    IF ( ntime >= ntime_start+1 ) THEN
      CALL GCM_DeallocateGhostCellArrays
!      DEALLOCATE(iGhostP)
!      DEALLOCATE(jGhostP)
!      DEALLOCATE(kGhostP)
    ENDIF ! ntime

! Allocate Arrays pertinent to Ghost Cells

    CALL GCM_AllocateGhostCellArrays()

! Set appropriate values for iGhost and jGhost by doing search
      
    nbdr = 0
    DO k = 1, nz-1
    DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( ghostCellMark(i,j,k) == 1 )   THEN
        nbdr         = nbdr + 1
        iGhost(nbdr) = i
        jGhost(nbdr) = j
        kGhost(nbdr) = k
      ENDIF ! ghostCellMark
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k

	
! Find marker point closest to each boundary node 
! and compute normal at closest marker point

    DO n = 1, nGhost

      iG = iGhost(n) 
      jG = jGhost(n)
      kG = kGhost(n)
        
      iCellIndex(n) = -1
      jCellIndex(n) = -1
      kCellIndex(n) = -1

      iBody = bodyNum(iG,jG,kG)
 
      xGC = xc(iG)
      yGC = yc(jG)
      zGC = zc(kG)

      CALL GCM_Calc_BodyIntercept_Unstruc( iG, jG, kG, xGC, yGC, zGC,    &
                                          xBI, yBI, zBI, closestElementGC(n) )

!-- Extract coordinates of Body Intercept

      xBodyIntercept(n) = xBI
      yBodyIntercept(n) = yBI
      zBodyIntercept(n) = zBI

!-- Get length of normal

	    dsIntercept = SQRT( (xGC-xBI)**2 + (yGC-yBI)**2 + (zGC-zBI)**2 )

!-- check intercept length against cell size
!-- if intercept is longer than cel diagonal then potential
! --sign of  oblique intecept.
      IF ( dsIntercept > SQRT(dxc(iG)**2 + dyc(jG)**2 + dzc(kG)**2) ) THEN
        PRINT*,' Normal intercept might not be correct'
        PRINT*,n
        PRINT*,iG,jG,kG
        PRINT*,dsIntercept,SQRT(dxc(iG)**2 + dyc(jG)**2 + dzc(kG)**2)
        PRINT*,dsIntercept/SQRT(dxc(iG)**2 + dyc(jG)**2 + dzc(kG)**2)
        PRINT*,'Check fort.198 for more info.'

        node1   = triElemNeig(iBody,1,closestElementGC(n))
        node2   = triElemNeig(iBody,2,closestElementGC(n))
        node3   = triElemNeig(iBody,3,closestElementGC(n))

        write(198,*)'ZONE'
        write(198,124)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
        write(198,124)xBodyMarker(iBody,node2),yBodyMarker(iBody,node2),zBodyMarker(iBody,node2)
        write(198,124)xBodyMarker(iBody,node3),yBodyMarker(iBody,node3),zBodyMarker(iBody,node3)
        write(198,124)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
        write(198,*)'ZONE'
        write(198,124)xGC,yGC,zGC
        write(198,124)xBI,yBI,zBI
        write(198,124)xIP,yIP,zIP
124       format(3(2x,e14.7))

        IF (dsIntercept/SQRT(dxc(iG)**2 + dyc(jG)**2 + dzc(kG)**2) > 2.0_CGREAL) THEN
          PRINT*,'Intercept is too long'
          CALL write_dump()
          STOP
        ENDIF

      ENDIF

         
!-- Now compute location of probe-tip (Image Point)
!-- Equation of 3D line  (parametric form)
!-- x = xo + slopeX . s
!-- y = yo + slopeY . s
!-- z = zo + slopeZ . s

      slopeX = (xBI-xGC)/dsIntercept
      slopeY = (yBI-yGC)/dsIntercept
      slopeZ = (zBI-zGC)/dsIntercept

      probeLength(n) = dsIntercept*probeLengthNormalized

      xImagePoint(n) = xGC + slopeX*probeLength(n)
      yImagePoint(n) = yGC + slopeY*probeLength(n)
      zImagePoint(n) = zGC + slopeZ*probeLength(n)

      xIP = xImagePoint(n)
      yIP = yImagePoint(n)
      zIP = zImagePoint(n)

 
! Find the lower left grid point to Image Point in Physical domain

      maxDelta = MAX(dx(iG),dy(jG),dz(kG))

!  Base range on probeLength

      iRange  = NINT(probeLength(n)/dx(iG)) +1
      jRange  = NINT(probeLength(n)/dy(jG)) +1
      kRange  = NINT(probeLength(n)/dz(kG)) +1

      iMin = iGhost(n)-iRange
      iMax = iGhost(n)+iRange
      iMin = MAX(iMin,0)    ! note image point is allowed to be between xc(0) and x(1)
      iMax = MIN(iMax,nx)   ! note image point is allowed to be between x(nx) and xc(nx)

      DO i = iMin,iMax 

        IF ( ( xc(i) <= xIP .AND. xc(i+1) > xIP ) ) THEN
          iCellIndex(n) = i
        ENDIF ! xc

      ENDDO ! i

      jMin = jGhost(n)-jRange
      jMax = jGhost(n)+jRange
      jMin = MAX(jMin,0)
      jMax = MIN(jMax,ny)

      DO j = jMin,jMax

        IF ( ( yc(j) <= yIP .AND. yc(j+1) > yIP ) ) THEN
          jCellIndex(n) = j
        ENDIF ! xc

      ENDDO ! j

      kMin = kGhost(n)-kRange
      kMax = kGhost(n)+kRange
      kMin = MAX(kMin,0)
      kMax = MIN(kMax,nz)

      DO k = kMin,kMax

        IF ( ( zc(k) <= zIP .AND. zc(k+1) > zIP ) ) THEN
          kCellIndex(n) = k
        ENDIF ! xc

      ENDDO ! k

      IF ( iCellIndex(n) ==-1 .OR. jCellIndex(n) ==-1 .OR. kCellIndex(n) ==-1 ) THEN
	      PRINT*,'Failed to Find four nodes surrounding an image point'
	      PRINT*,n
        PRINT*,iG, jG, kG
        PRINT*,xgc,ygc,zgc
        PRINT*,xBI,yBI,zBI
        PRINT*,xIP,yIP,zIP
        PRINT*,'Aborting Run'
        STOP
      ENDIF

! Perform bilinear interpolation
 
      iCIndx = iCellIndex(n)
      jCIndx = jCellIndex(n)
      kCIndx = kCellIndex(n)
 
      xBIN = triElemNormx(iBody,closestElementGC(n))
      yBIN = triElemNormy(iBody,closestElementGC(n))
      zBIN = triElemNormz(iBody,closestElementGC(n))
      CALL GCM_Calc_vanMatrixDN( iG, jG, kG, iCIndx, jCIndx, kCIndx,             &
                                xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                                coeffGCMDirc, coeffGCMNeum      )
 
      coeffGCMD(1:iRowMax,n) = coeffGCMDirc(1:iRowMax)
      coeffGCMN(1:iRowMax,n) = coeffGCMNeum(1:iRowMax)

! Test Closest Marker with 25-point Stencil

!        IF ( (iCellIndex(n)-iGhost(n)) > +2 .OR.  &
!             (iCellIndex(n)-iGhost(n)) < -2 .OR.  &
!             (jCellIndex(n)-jGhost(n)) > +2 .OR.  &
!             (jCellIndex(n)-jGhost(n)) < -2 .OR.  &
!             (kCellIndex(n)-kGhost(n)) > +2 .OR.  &
!             (kCellIndex(n)-kGhost(n)) < -2       )  THEN
!
!           PRINT*,'CLOSEST MARKER IS NOT INSIDE 25-POINT STENCIL'
!        ENDIF ! iCellIndex

    ENDDO ! n 

! Deallocate arrays

    DEALLOCATE( coeffGCMDirc )
    DEALLOCATE( coeffGCMNeum )

  END SUBROUTINE GCM_set_internal_boundary 
!----------------------------------------------------------------------

  SUBROUTINE GCM_Calc_VanMatrixDN( iG, jG, kG, iCIndex,jCIndex, kCIndex,      &
                                  xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                                  coeffGCMDirc, coeffGCMNeum      )
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE gcm_arrays

    IMPLICIT NONE

!... parameters variables

    INTEGER, INTENT(IN)           :: iG, jG, kG, iCIndex,jCIndex, kCIndex
    REAL(KIND=CGREAL), INTENT(IN) :: xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN
    REAL(KIND=CGREAL), DIMENSION(iRowMax), INTENT(OUT) :: coeffGCMDirc, &
                                                          coeffGCMNeum

!... loop variables

    INTEGER :: i,j,k,iRow
    INTEGER :: info

!... local variables
    
    REAL(KIND=CGREAL) :: rCond, xC1,xC2,xC3, xN1,xN2,xN3

!*****************************************************************************************
  
!   |-------|-------|---/---|-------|--         N : Nth ghost point
!   |   ii  |  iii  |  *    |       |           * : markers
!   |   0...|...O   | / .   |   .   |           O : other nodes used in bilinear interpolation
!   |   .   |   .   |*      |       |           + : probe tip (Image Point) 
!   |---.---|--+.---/-------|-------|--
!   |   .   |   .  *|       |       |
!   |   0...|. .O / |   N   |   .   |
!   |   i   |  iv*  |       |       |
!   |-------| --/ --|-------|-------|--

! interpolant      U = a X X X  + b X X  + c X X  + d X X
!                         1 2 3      1 2      1 3      2 3
!
!                    + e X  + f X + g X  + h
!                         1      2     3
!
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
!     iv  [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ d ]     [ U   ]
!         [                       ]  [   ]     [  iv ]
!
!
!   Van Matrix For Dirichlet conditions at Intersection Point (N)
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
!         [                       ]  [   ]     [  N  ]
!
!   Van Matrix For Neumann conditions at Intersection point (N)
!    B1 = n_x, B2 = n_y (components of normal vectors)
!    F_m = value of normal derivative 
!
!         [  X X           X     X   1  ]  [   ]     [     ] 
!      i  [   1 2           1     2     ]  [ a ]     [ U   ]
!         [                             ]  [   ]     [  i  ]
!         [  X X           X     X   1  ]  [   ]     [     ]
!      ii [   1 2           1     2     ]  [ b ]     [ U   ]
!         [                             ]  [   ]  =  [  ii ]
!     iii [  X X           X     X   1  ]  [   ]     [     ]
!         [   1 2           1     2     ]  [ c ]     [ U   ]
!         [                             ]  [   ]     [  iii]
!      N  [  B X  + B X    B     B   0  ]  [   ]     [     ]
!         [   1 2    2  1   1     2     ]  [ d ]     [ F   ]
!         [                             ]  [   ]     [  N  ]
!

    DO iRow= 1, iRowMax
      i  = iCIndex + incI(iRow)
      j  = jCIndex + incJ(iRow)
      k  = kCIndex + incK(iRow)

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

!--- Neumann conditions for pressure field


      vanMatrixN(iRow,1) = xC1*xC2*xC3
      vanMatrixN(iRow,2) = xC1*xC2
      vanMatrixN(iRow,3) = xC1*xC3
      vanMatrixN(iRow,4) = xC2*xC3
      vanMatrixN(iRow,5) = xC1
      vanMatrixN(iRow,6) = xC2
      vanMatrixN(iRow,7) = xC3
      vanMatrixN(iRow,8) = oned

!-- Correct For Ghost node part of cell formation, switch to Body Intercept point

      IF ( i==iG .AND. j == jG  .AND. k== kG) THEN
        xC1 = xBI
        xC2 = yBI
        xC3 = zBI
        xN1 = xBIN
        xN2 = yBIN
        xN3 = zBIN

        vanMatrixD(iRow,1) = xC1*xC2*xC3
        vanMatrixD(iRow,2) = xC1*xC2
        vanMatrixD(iRow,3) = xC1*xC3
        vanMatrixD(iRow,4) = xC2*xC3
        vanMatrixD(iRow,5) = xC1
        vanMatrixD(iRow,6) = xC2
        vanMatrixD(iRow,7) = xC3
        vanMatrixD(iRow,8) = oned

        vanMatrixN(iRow,1) = xN1*xC2*XC3 + xN2*xC1*XC3 + xN3*XC1*XC2
        vanMatrixN(iRow,2) = xN1*xC2 + xN2*xC1
        vanMatrixN(iRow,3) = xN1*xC3 + xN3*xC1
        vanMatrixN(iRow,4) = xN2*xC3 + xN3*xC2
        vanMatrixN(iRow,5) = xN1
        vanMatrixN(iRow,6) = xN2
        vanMatrixN(iRow,7) = xN3
        vanMatrixN(iRow,8) = zero

      ENDIF ! i
    ENDDO ! iRow		

! Compute inverse of Vandermonde Matrices

    CALL DGETRF(8, 8, vanMatrixD,8,iPvt, info) 
    CALL DGETRI(8, vanMatrixD,8,iPvt,work, 8, info) 

    CALL DGETRF(8, 8, vanMatrixN,8,iPvt, info)
    CALL DGETRI(8, vanMatrixN,8,iPvt,work, 8, info)

! Load Coeff-Matrices

    DO iRow = 1, iRowMax
      coeffGCMDirc(iRow) = vanMatrixD(1,iRow)*xIP*yIP*zIP  &
                         + vanMatrixD(2,iRow)*xIP*yIP      &
                         + vanMatrixD(3,iRow)*xIP*zIP      &
                         + vanMatrixD(4,iRow)*yIP*zIP      &
                         + vanMatrixD(5,iRow)*xIP          &
                         + vanMatrixD(6,iRow)*yIP          &
                         + vanMatrixD(7,iRow)*zIP          &
                         + vanMatrixD(8,iRow)

      coeffGCMNeum(iRow) = vanMatrixN(1,iRow)*xIP*yIP*zIP  &
                         + vanMatrixN(2,iRow)*xIP*yIP      &
                         + vanMatrixN(3,iRow)*xIP*zIP      &
                         + vanMatrixN(4,iRow)*yIP*zIP      &
                         + vanMatrixN(5,iRow)*xIP          &
                         + vanMatrixN(6,iRow)*yIP          &
                         + vanMatrixN(7,iRow)*zIP          &
                         + vanMatrixN(8,iRow)
    ENDDO ! iRow 

  END SUBROUTINE GCM_Calc_VanMatrixDN
!----------------------------------------------------------------------
SUBROUTINE identify_ghostcells_solid()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

!... loop variables
    
    INTEGER :: i,j,k

!... local variables

    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
    INTEGER :: iM, iP, jM, jP, kM, kP



!      iblank         = iblank_solid
!      ghostCellMark  = 0
    ghostCellSolid = 0
!
! Mark boundary cells (ie. cells inside body which have 
! at least one neighbor in fluid)

    DO k = 1, nz-1
    DO j = 1, ny-1
    DO i = 1, nx-1

      iM = MAX(i-1,1)
      iP = MIN(i+1,nx-1)
      jM = MAX(j-1,1)
      jP = MIN(j+1,ny-1)
      kM = MAX(k-1,1)
      kP = MIN(k+1,nz-1)

      IF ( iblank(i,j,k)/=iblank(iM,j,k)   .OR.  &
           iblank(i,j,k)/=iblank(iP,j,k)   .OR.  &
           iblank(i,j,k)/=iblank(i,jM,k)   .OR.  &
           iblank(i,j,k)/=iblank(i,jP,k)   .OR.  &
           iblank(i,j,k)/=iblank(i,j,kM)   .OR.  &
           iblank(i,j,k)/=iblank(i,j,kP)         )  THEN
        IF (iblank(i,j,k)==1) THEN
          ghostCellSolid(i,j,k) = -1
!           IF (iblank(i,j,k) == 1) &
          ghostCellMark(i,j,k)  = -1
        ELSE
          ghostCellSolid(i,j,k) = 1
!           IF (iblank(i,j,k) == 1) &
          ghostCellMark(i,j,k)  = 1
		    END IF          
      ENDIF

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

END SUBROUTINE identify_ghostcells_solid 
!---------------------------------------------------------------------



SUBROUTINE identify_ghostcells_membrane(iBody)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

    INTEGER, INTENT(IN):: iBody
!... loop variables
    INTEGER :: i,iRow,j,k,m,n,e

!... local variables

    INTEGER(1) :: boundCell(0:nx+1,0:ny+1,0:nz+1), ghostCell(0:nx+1,0:ny+1,0:nz+1) !ghostCell is for current membrane
    INTEGER :: num_fresh
    INTEGER :: iG, jG, kG, nbdr, nbdrG,ic,jc,kc
    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
    INTEGER :: iM, iP, jM, jP, kM, kP
    INTEGER :: iRange, jRange, kRange
    INTEGER :: cElemG

    REAL(KIND=CGREAL) :: xBIT,yBIT,zBIT
    REAL(KIND=CGREAL) :: xGC,yGC,zGC
    
    REAL(KIND=CGREAL) :: vert1(3),vert2(3),vert3(3)
    REAL(KIND=CGREAL) :: xMin,xMax,yMin,yMax,zMin,zMax


    num_fresh  = 0

! Find all cells that contain a IB node

    boundCell = 0

    hybridMarkMemb=0

    !DO m = 1,nPtsBodyMarker(iBody)
    DO e = 1,totNumTriElem(iBody)
        
      vert1(1)=xBodyMarker(iBody,triElemNeig(iBody,1,e))
      vert1(2)=yBodyMarker(iBody,triElemNeig(iBody,1,e))
      vert1(3)=zBodyMarker(iBody,triElemNeig(iBody,1,e))
      
      vert2(1)=xBodyMarker(iBody,triElemNeig(iBody,2,e))
      vert2(2)=yBodyMarker(iBody,triElemNeig(iBody,2,e))
      vert2(3)=zBodyMarker(iBody,triElemNeig(iBody,2,e))
      
      vert3(1)=xBodyMarker(iBody,triElemNeig(iBody,3,e))
      vert3(2)=yBodyMarker(iBody,triElemNeig(iBody,3,e))
      vert3(3)=zBodyMarker(iBody,triElemNeig(iBody,3,e))
      
      xMin=minval((/vert1(1),vert2(1),vert3(1)/))
      xMax=maxval((/vert1(1),vert2(1),vert3(1)/))
      yMin=minval((/vert1(2),vert2(2),vert3(2)/))
      yMax=maxval((/vert1(2),vert2(2),vert3(2)/))
      zMin=minval((/vert1(3),vert2(3),vert3(3)/))
      zMax=maxval((/vert1(3),vert2(3),vert3(3)/))
      
      do i=1,nx
          if(xc(i)<=xMin.and.xc(i+1)>=xMin) iMin=i
          if(xc(i)<=xMax.and.xc(i+1)>=xMax) iMax=i+1
      end do
      
      do j=1,ny
          if(yc(j)<=yMin.and.yc(j+1)>=yMin) jMin=j
          if(yc(j)<=yMax.and.yc(j+1)>=yMax) jMax=j+1
      end do
      
      !IF (ndim /= DIM_3D ) THEN
      !    kMin=1
      !    kMax=2
      !else
          do k=1,nz
              if(zc(k)<=zMin.and.zc(k+1)>=zMin) kMin=k
              if(zc(k)<=zMax.and.zc(k+1)>=zMax) kMax=k+1
          end do
      !end if
      
      do k=kMin,kMax
      do j=jMin,jMax
      do i=iMin,iMax
          boundCell(i,j,k)=1
      end do
      end do
      end do
    

    !  iC = -1
    !  jC = -1
    !  kC = -1
    !  
    !  DO i = 1,nx
    !    IF ( ( xc(i) <= xBodyMarker(iBody,m)  .AND. xc(i+1) > xBodyMarker(iBody,m) ) ) THEN
    !      iC = i
    !    ENDIF ! xc
    !  ENDDO ! i
    !  DO j = 1,ny
    !    IF ( ( yc(j) <= yBodyMarker(iBody,m)  .AND. yc(j+1) > yBodyMarker(iBody,m) ) ) THEN
    !      jC = j
    !    ENDIF ! yc
    !  ENDDO ! j
    !  IF (ndim /= DIM_3D ) THEN
    !    kC = 1 
    !  ELSE 
    !    DO k = 1,nz
    !      IF ( ( zc(k) <= zBodyMarker(iBody,m)  .AND. zc(k+1) > zBodyMarker(iBody,m) ) ) THEN
    !        kC = k
    !      ENDIF ! zc
    !    ENDDO ! k
    !  ENDIF
    !  
    !  IF ( iC == -1 .OR. jC == -1 .OR. kC == -1 ) CALL abort_vicar3d(30)
    !      
    !  boundCell(iC,jC,kC)       = 1
    !  boundCell(iC+1,jC,kC)     = 1
    !  boundCell(iC,jC+1,kC)     = 1
    !  boundCell(iC,jC,kC+1)     = 1
    !  boundCell(iC+1,jC+1,kC)   = 1
    !  boundCell(iC,jC+1,kC+1)   = 1
    !  boundCell(iC+1,jC,kC+1)   = 1
    !  boundCell(iC+1,jC+1,kC+1) = 1
    !
    ENDDO
!
! Mark boundary cells (ie. cells inside body which have 
! at least one neighbor in fluid)

    ghostCell = 0
        
    DO k = 1, nz-1
    DO j = 1, ny-1
    DO i = 1, nx-1

      iM = MAX(i-1,1)
      iP = MIN(i+1,nx-1)
      jM = MAX(j-1,1)
      jP = MIN(j+1,ny-1)
      kM = MAX(k-1,1)
      kP = MIN(k+1,nz-1)

      IF (boundCell(i,j,k) == 1) THEN

        IF ( ( iblank(i,j,k)+iblank(iM,j,k) == 1 .AND. boundCell(iM,j,k) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(iP,j,k) == 1 .AND. boundCell(iP,j,k) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(i,jM,k) == 1 .AND. boundCell(i,jM,k) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(i,jP,k) == 1 .AND. boundCell(i,jP,k) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(i,j,kM) == 1 .AND. boundCell(i,j,kM) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(i,j,kP) == 1 .AND. boundCell(i,j,kP) == 1 )      )  THEN
            bodyNum(i,j,k)=iBody
            !conflictCell(i,j,k)=conflictCell(i,j,k)+2**(iBody-1)
            if (iblank(i,j,k)==1) then 
              ghostCell(i,j,k) = 1
            else
              ghostCell(i,j,k) = -1
            end if
! check fresh cell generated when membrane moves a little bit
! commented out for consistency with v11.0 since the it is hard to interpolate the value of fresh cell for membrane
!             IF (ghostCell(i,j,k)*ghostCellMemb(i,j,k)==-1 .AND. ntime > ntime_start+1) THEN
!              fresh_cell(i,j,k)=  1
!              num_fresh        = num_fresh+1
!              WRITE(ifuFreshCellOut,*)ntime,i,j,k,'   --- fresh cell'
!             END IF
        ENDIF

      ENDIF

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!    if(ntime==32)then
!    write(568,*) 'VARIABLES="X","Y","ghostCell"'
!    write(568,*) 'ZONE F=POINT, I=',nxc,' , J=',nyc
!    do j=1,nyc
!    do i=1,nxc
!        write(568,*) xc(i),yc(j),ghostCell(i,j,1)
!    end do
!    end do
!    end if

    hybridMarkMemb=ghostCell

    IF (boundary_formulation == GCM_METHOD) THEN

!  Determine which ghost cells are real and calculate body intercept 
!  and closest element for these ghost cells.

      nbdr  = 0
      nbdrG = 0

      DO k = 1, nz-1
      DO j = 1, ny-1
      DO i = 1, nx-1

        IF ( ghostCell(i,j,k) == 1 .or. ghostCell(i,j,k) == -1)   THEN
          iG    = i
          jG    = j
          kG    = k
          xGC   = xc(iG)
          yGC   = yc(jG)
          zGC   = zc(kG)

          nbdr = nbdr + 1
     
          CALL GCM_Calc_BodyIntercept_Unstruc( iG, jG, kG, xGC, yGC, zGC, xBIT, yBIT, zBIT, cElemG )

          IF ( cElemG == -1 )  THEN    ! did not find an intercept
            ghostCell(i,j,k) = -1
          ENDIF
        ENDIF ! ghostCellMemb

      ENDDO ! i 
      ENDDO ! j
      ENDDO ! k

      IF ( MOD(ntime,nmonitor) == 0 .OR. ntime == 1) &
        PRINT*,'GCM_set_internal_boundary: nGhost = ',nbdr

    ENDIF ! boundary_formulation

    
    DO k=1,nz-1
    DO j=1,ny-1
    DO i=1,nx-1
      ghostCellMark(i,j,k) = ghostCellMark(i,j,k) + ghostCell(i,j,k)
      IF ( ghostCellMark(i,j,k) > 1 .or. ghostCellMark(i,j,k) < -1) THEN
        PRINT*,'Warning: absolute value of ghostCellMark is > 1; ', i,j,k,ghostCellMark(i,j,k)
!        STOP
      ENDIF 
      IF (iblank_memb(i,j,k)>0 .OR. ghostCellMark(i,j,k)/=0) &
        iblank_memb(i,j,k) = 1
    ENDDO
    ENDDO
    ENDDO

    IF ( ntime==1           .OR. &
         ntime==ntime_start .OR. &
         MOD(ntime,nmonitor) == 0 )  THEN
      num_iblank(iBody)=sum(INT(iblank_memb(1:nx-1,1:ny-1,1:nz-1)))-num_iblank(iBody-1)
      WRITE(*,*) 'SUM of IBLANK = ', num_iblank(iBody)
      PRINT*,'Number of Fresh Cells = ',num_fresh
    END IF ! ntime

END SUBROUTINE identify_ghostcells_membrane
!---------------------------------------------------------------------


