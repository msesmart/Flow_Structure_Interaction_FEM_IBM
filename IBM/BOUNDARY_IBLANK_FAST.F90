!------------------------------------------------------------------------------
   SUBROUTINE set_iblank_canonical_body_fast(nBody_Begin,nBody_End)
!------------------------------------------------------------------------------
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nbody_Begin, nbody_End
    INTEGER, PARAMETER :: UNDECIDED_VALUE = 1

    INTEGER           :: i,j,k,m,n,m_theta,m_phi,np,num_fresh,num_dead,inside
    INTEGER           :: cMarker,cElement,nBodySelect
    INTEGER           :: ii,jj,kk
    INTEGER           :: iMin,iMax,jMin,jMax,kMin,kMax
    INTEGER           :: iCellMin,iCellMax,jCellMin,jCellMax,kCellMin,kCellMax
    INTEGER           :: scanRange,sumIblank,tmpIblank
    INTEGER           :: iBeg,jBeg,kBeg
    INTEGER           :: iErr,nLevelRefine,nSubDiv
    INTEGER           :: iclock1,iclock2,iclock3,iclock4,iclock5,iclock6,clock_rate
    INTEGER           :: iVert,mVert1,mVert2,mVert3,nVertTot
    INTEGER           :: iTri,jTri,kTri
    INTEGER           :: cElement2,nBodySelect2,tmpIblank2
    INTEGER(1) :: iblankUndecided(0:nx+1,0:ny+1,0:nz+1)
    INTEGER(1) :: iblankTemp(0:nx+1,0:ny+1,0:nz+1)
    REAL(KIND=CGREAL) :: rad,theta,phi,x_prime,y_prime,z_prime
    REAL(KIND=CGREAL) :: distMarker,distElement,dMin,dotNorm,dMinUnstruc
    REAL(KIND=CGREAL) :: xBM,yBM,zBM
    REAL(KIND=CGREAL) :: xTemp,yTemp,zTemp
    REAL(KIND=CGREAL) :: lenGridMin,lenElemMax,ratioElemGrid,checkDist
    REAL(KIND=CGREAL) :: rSubDivInv,uTri,vTri,wTri
    REAL(KIND=CGREAL) :: dMinUnstruc0,dMinUnstruc1,dMinUnstruc2,dMinUnstruc3
    REAL(KIND=CGREAL) :: dotNorm1,dotNorm2,dotNorm3
    REAL(KIND=CGREAL) :: xBoundMax,xBoundMin,yBoundMax,yBoundMin,&
                         zBoundMax,zBoundMin
    REAL(KIND=CGREAL) :: xElem,yElem,zElem
    REAL(KIND=CGREAL) :: distMin1,distMin2
    REAL(KIND=CGREAL) :: distxBoundMin,distxBoundMax
    REAL(KIND=CGREAL) :: distyBoundMin,distyBoundMax
    REAL(KIND=CGREAL) :: distzBoundMin,distzBoundMax
    REAL(KIND=CGREAL), DIMENSION(3) :: delta,lenElem,lenGrid 
    REAL(KIND=CGREAL), DIMENSION(3) :: xVert,yVert,zVert
 
! ============================================================================
!   Initialize values
! ============================================================================

    num_fresh  = 0
    num_dead   = 0
 
! ============================================================================
!   Allocate local iblank array
! ============================================================================

    iblankUndecided = 0
    iblankTemp      = 0
    
! ============================================================================
!   Initialize iblank to undecided value over all bodies in domain
! ============================================================================

!    iblank = 0
! DEBUG
!   FMN: To test is all iblanks are being filled.
!     iblank = -20*UNDECIDED_VALUE
! END DEBUG

!    CALL write_dump_debug_i4('ByM ',nBody_Begin,Bodynum)

    SELECT CASE (nDim)
      CASE (DIM_2D)
        k=1
        DO j=0,ny
        DO i=0,nx
          iblankUndecided(i,j,k) = UNDECIDED_VALUE
        END DO ! i
        END DO ! j
      CASE (DIM_3D)
        DO k=0,nz
        DO j=0,ny
        DO i=0,nx
          iblankUndecided(i,j,k) = UNDECIDED_VALUE
        END DO ! i
        END DO ! j
        END DO ! k
    END SELECT ! nDim

! ============================================================================
!   Loop over all bodies in domain and body surface elements
! ============================================================================

    CALL system_clock(iclock1)
    DO n=nbody_Begin, nbody_End
      DO m=1,totNumTriElem(n)

! ******************************************************************************
!       Extract cell indices and vertices of triangular element
! ******************************************************************************

         mVert1    = triElemNeig(n,1,m)
         mVert2    = triElemNeig(n,2,m)
         mVert3    = triElemNeig(n,3,m)

         xElem     = triElemCentx(n,m)
         yElem     = triElemCenty(n,m)
         zElem     = triElemCentz(n,m)

         xVert(1)  = xBodyMarker(n,mVert1)
         yVert(1)  = yBodyMarker(n,mVert1)
         zVert(1)  = zBodyMarker(n,mVert1)

         xVert(2)  = xBodyMarker(n,mVert2)
         yVert(2)  = yBodyMarker(n,mVert2)
         zVert(2)  = zBodyMarker(n,mVert2)

         xVert(3)  = xBodyMarker(n,mVert3)
         yVert(3)  = yBodyMarker(n,mVert3)
         zVert(3)  = zBodyMarker(n,mVert3)

! ******************************************************************************
!        cycle if element is outside computational domain
! ******************************************************************************

         IF ( ( xVert(1) <= x(1)  .AND. xVert(2) <= x(1)          &
                                  .AND. xVert(3) <= x(1)  ) .OR.  &
              ( xVert(1) >= x(nx) .AND. xVert(2) >= x(nx)         &
                                  .AND. xVert(3) >= x(nx) ) .OR.  &
              ( yVert(1) <= y(1)  .AND. yVert(2) <= y(1)          &
                                  .AND. yVert(3) <= y(1)  ) .OR.  &
              ( yVert(1) >= y(ny) .AND. yVert(2) >= y(ny)         &
                                  .AND. yVert(3) >= y(ny) ) .OR.  &
              ( zVert(1) <= z(1)  .AND. zVert(2) <= z(1)          &
                                  .AND. zVert(3) <= z(1)  ) .OR.  &
              ( zVert(1) >= z(nz) .AND. zVert(2) >= z(nz)         &
                                  .AND. zVert(3) >= z(nz) )     ) CYCLE

! ******************************************************************************
!        Find all cells within the bounding box of the vertices  
!        for each element
! ******************************************************************************


         xBoundMin = MINVAL(xVert(1:3))
         xBoundMax = MAXVAL(xVert(1:3))

         yBoundMin = MINVAL(yVert(1:3))
         yBoundMax = MAXVAL(yVert(1:3))

         zBoundMin = MINVAL(zVert(1:3))
         zBoundMax = MAXVAL(zVert(1:3))

!DEBUG
 !        WRITE(*,*)'n,m = ',n,m
 !        WRITE(*,*)'xBoundMin-Max =',xBoundMin,xBoundMax
 !        WRITE(*,*)'yBoundMin-Max =',yBoundMin,yBoundMax
 !        WRITE(*,*)'zBoundMin-Max =',zBoundMin,zBoundMax
 !        WRITE(*,*)'xMin-Max     =',MINVAL(x(1:nx)),MAXVAL(x(1:nx))
 !        WRITE(*,*)'yMin-Max     =',MINVAL(y(1:ny)),MAXVAL(y(1:ny))
 !        WRITE(*,*)'zMin-Max     =',MINVAL(z(1:nz)),MAXVAL(z(1:nz))
!END DEBUG

         iCellMin = 0
         iCellMax = 0
         jCellMin = 0
         jCellMax = 0
         kCellMin = 0
         kCellMax = 0

! ******************************************************************************
!        i-direction 
! ******************************************************************************

         iMin = 0
         iMax = nx
         DO i = iMin, iMax-1
           IF ( (xc(i)-xBoundMin)*(xc(i+1)-xBoundMin) <= zero ) &
              iCellMin = i
           IF ( (xc(i)-xBoundMax)*(xc(i+1)-xBoundMax) <= zero ) &
              iCellMax = i+1
         ENDDO
         IF (iCellMin == 0 )                    iCellMin = 1
         IF (iCellMax == 0 .OR. iCellMax == nx) iCellMax = nxc

! ******************************************************************************
!        j-direction 
! ******************************************************************************

         jMin = 0
         jMax = ny
         DO j = jMin, jMax-1
           IF ( (yc(j)-yBoundMin)*(yc(j+1)-yBoundMin) <= zero ) &
              jCellMin = j
           IF ( (yc(j)-yBoundMax)*(yc(j+1)-yBoundMax) <= zero ) &
              jCellMax = j+1
         ENDDO

         IF (jCellMin == 0 )                    jCellMin = 1
         IF (jCellMax == 0 .OR. jCellMax == ny) jCellMax = nyc


! ******************************************************************************
!        k-direction 
! ******************************************************************************

         SELECT CASE (nDim)
           CASE (DIM_2D)
             kCellMin = 1
             kCellMax = nzc

           CASE (DIM_3D)

            kMin = 0
            kMax = nz
            DO k = kMin, kMax-1
              IF ( (zc(k)-zBoundMin)*(zc(k+1)-zBoundMin) <= zero ) &
              kCellMin = k
              IF ( (zc(k)-zBoundMax)*(zc(k+1)-zBoundMax) <= zero ) &
              kCellMax = k+1
            ENDDO
            IF (kCellMin == 0 ) kCellMin = 1
            IF (kCellMax == 0 .OR. kCellMax == nz) kCellMax = nzc

         END SELECT ! nDim

! DEBUG
!      WRITE(*,*) ' m   = ',m
!      WRITE(*,*) 'iCellMin-Max = ',iCellMin,iCellMax
!      WRITE(*,*) 'jCellMin-Max = ',jCellMin,jCellMax
!      WRITE(*,*) 'kCellMin-Max = ',kCellMin,kCellMax
!END DEBUG

! ******************************************************************************
!        set iblankUndecided to NEGATIVE undecided value for the ring 
!            of cells extracted. Also save bodyNum
! ******************************************************************************

         DO k = kCellMin, kCellMax
         DO j = jCellMin, jCellMax
         DO i = iCellMin, iCellMax
              iblankUndecided(i,j,k)  = -UNDECIDED_VALUE
         END DO ! i 
         END DO ! j 
         END DO ! k 

      END DO ! m
    END DO ! n
    CALL system_clock(iclock3)

! ============================================================================
!   Loop over cells whose iblank is -UNDECIDED_VALUE for unstructured surfaces
!    invoking dot normal algorithm 
! ============================================================================

    CALL system_clock(iclock4)
    SELECT CASE(nDim)
      CASE (DIM_2D)
        kMin=1
        kMax=1

      CASE (DIM_3D)
        kMin=1
        kMax=nzc
    END SELECT ! nDim

!    CALL write_dump_debug_i('unde',nBody_Begin,iblankUndecided)

    DO k=kMin,kMax
    DO j=1,nyc
    DO i=1,nxc

      tmpiblank    = 0

      IF ( iblankUndecided(i,j,k) == -UNDECIDED_VALUE ) THEN
! ******************************************************************************
!       Search for vertex closest to cells in band
!       Extract the elements that share that vertex
!       Drop the normal and find BI point
!       Check if normal intercept lies inside closest triangular element
!       Authors: Fady and Rajat Oct 18,2005
! ******************************************************************************
        DO n = nbody_Begin, nbody_End
          
          CALL search_vertex_dotNorm(n,i,j,k,cElement,dotNorm)

          IF (dotNorm >= zero) THEN
            iblankTemp(i,j,k)  = 1
            bodyNum(i,j,k)     = n
            GOTO 888
          ENDIF
 
        ENDDO

        iblankTemp(i,j,k)  = 0

888 CONTINUE

      END IF ! iblankUndecided

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
    CALL system_clock(iclock5)

!    CALL write_dump_debug_i('blkt',nBody_Begin,iblankTemp)

!DEBUG
!    write(921,*)'VARIABLES="X","Y","Z","IBUNDEC","IBLANK", "BodyNum"'
!    write(921,*)'ZONE F=POINT, I=',nxc,', J=',nyc,' ,K=',nzc
!    do k=1,nzc
!    do j=1,nyc
!    do i=1,nxc
!       write(921,'(3(3X,1PE12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblankUndecided(i,j,k),iblankTemp(i,j,k),bodyNum(i,j,k)
!    enddo
!    enddo
!    enddo
!DEBUG

     IF (unstruc_surface_type(nBody_Begin) == MEMBRANE) GOTO 1000

! ============================================================================
!   Set undecided iblank values outside the ring of cells for body
!    by searching horizontal, similar to a ray tracing routine. 
!    Set iblank value at grid cell by searching for first DECIDED VALUE
!    Move along i-direction, j-direction, then k-direction
! ============================================================================

! ******************************************************************************
!   i-direction 
! ******************************************************************************

    CALL system_clock(iclock6)
    iBeg = 0
    DO k=kMin,kMax
    DO j=1,nyc
      initLoopI: DO ii=1,nx
          IF ( iblankUndecided(ii,j,k) == UNDECIDED_VALUE ) CYCLE
            iblankTemp(iBeg,j,k) = iblankTemp(ii,j,k)
!           bodyNum(iBeg,j,k)    = bodyNum(ii,j,k)
            EXIT initLoopI
      END DO  initLoopI
    END DO ! j
    END DO ! k

    DO i=1,nx
      DO k=kMin,kMax
      DO j=1,nyc
        IF ( iblankUndecided(i,j,k) /= UNDECIDED_VALUE ) CYCLE
            iblankTemp(i,j,k)  = iblankTemp(i-1,j,k)
            IF (iblankTemp(i,j,k)==1) bodyNum(i,j,k)     = bodyNum(i-1,j,k)
      END DO ! j
      END DO ! k
    END DO ! i


! ******************************************************************************
!   j-direction 
! ******************************************************************************

    jBeg = 0
    DO k=kMin,kMax
    DO i=1,nxc
      initLoopJ: DO jj=1,ny
          IF ( iblankUndecided(i,jj,k) == UNDECIDED_VALUE ) CYCLE
          iblankTemp(i,jBeg,k) = iblankTemp(i,jj,k)
!         bodyNum(i,jBeg,k)    = bodyNum(i,jj,k)
          EXIT initLoopJ
      END DO  initLoopJ
    END DO ! i
    END DO ! k

    DO j=1,ny
      DO k=kMin,kMax
      DO i=1,nxc
        IF ( iblankUndecided(i,j,k) /= UNDECIDED_VALUE ) CYCLE
        iblankTemp(i,j,k) = iblankTemp(i,j-1,k)
        IF (iblankTemp(i,j,k)==1) bodyNum(i,j,k)    = bodyNum(i,j-1,k)
      END DO ! i
      END DO ! k
    END DO ! j

! ******************************************************************************
!   k-direction 
! ******************************************************************************

    IF (nDim == DIM_3D) THEN
      kBeg = 0
      DO j=1,nyc
      DO i=1,nxc
        initLoopK: DO kk=kMin,kMax+1
          IF ( iblankUndecided(i,j,kk) == UNDECIDED_VALUE ) CYCLE
            iblankTemp(i,j,kBeg)  = iblankTemp(i,j,kk)
!           bodyNum(i,j,kBeg)     = bodyNum(i,j,kk)
            EXIT initLoopK
        END DO  initLoopK
      END DO ! i
      END DO ! j

      DO k=kMin,kMax+1
        DO j=1,nyc
        DO i=1,nxc
          IF ( iblankUndecided(i,j,k) /= UNDECIDED_VALUE ) CYCLE
          iblankTemp(i,j,k)  = iblankTemp(i,j,k-1)
          IF (iblankTemp(i,j,k)==1) bodyNum(i,j,k)     = bodyNum(i,j,k-1)
        END DO ! i
        END DO ! j
      END DO   ! k

    ENDIF !nDim

1000 CONTINUE
! ============================================================================
!   Extend iblank for 2D simulations 
! ============================================================================

    IF (nDim == DIM_2D) THEN
      DO k=kMin+1,nz
      DO j=1,nyc
      DO i=1,nxc
        iblankTemp(i,j,k)  = iblankTemp(i,j,1)
        IF ( iblankTemp(i,j,k) == 1 ) bodyNum(i,j,k) = bodyNum(i,j,1)
      END DO !i
      END DO !j
      END DO !k
    END IF ! nDim
    
!!! The fresh cell for membrance is set in identify_ghostcells_membrane
    IF (unstruc_surface_type(nBody_Begin) == MEMBRANE) THEN
    
      DO k=1,nzc
      DO j=1,nyc
      DO i=1,nxc
        iblank(i,j,k)  = iblankTemp(i,j,k)
      END DO ! i
      END DO ! j
      END DO ! k

    ELSE
    DO k=1,nzc
    DO j=1,nyc
    DO i=1,nxc
          IF  ( iblankTemp(i,j,k) == 1 ) THEN
            IF ( iblank(i,j,k) == 0  .AND. ntime > ntime_start+1 ) THEN
               WRITE(ifuFreshCellOut,*)ntime,i,j,k,'   --- dead cell'
               num_dead     = num_dead+1
            ENDIF
          ENDIF
          IF  ( iblankTemp(i,j,k) == 0 ) THEN
            IF ( iblank(i,j,k) == 1 .AND. ntime > 1 ) THEN
              fresh_cell(i,j,k)=  1
              bodyNum(i,j,k)   = bodyNum(i,j,k)  ! this is done to show that this array has correct value
              num_fresh        = num_fresh+1
              WRITE(ifuFreshCellOut,*)ntime,i,j,k,'   --- fresh cell'
            ENDIF
          ENDIF
        iblank(i,j,k)  = iblankTemp(i,j,k)
    END DO ! i
    END DO ! j
    END DO   ! k
    
    IF (MOD(ntime,nmonitor)==0) THEN
      PRINT*,'Number of Dead  Cells = ',num_dead
      PRINT*,'Number of Fresh Cells = ',num_fresh
    ENDIF ! ntime

    CALL system_clock(iclock2,clock_rate)
    IF ( ntime==1           .OR. &
         ntime==ntime_start .OR. &
         MOD(ntime,nmonitor) == 0 )  THEN
      WRITE(*,*) 'CPU Time for Fast Initial Iblank = ',&
      REAL(iclock2-iclock1,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
      WRITE(*,*) '    CPU Time for Initial Setup = ',&
      REAL(iclock3-iclock1,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
      WRITE(*,*) '    CPU Time for dotNorm = ',&
      REAL(iclock5-iclock4,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
      WRITE(*,*) '    CPU Time for Filling Iblank = ',&
      REAL(iclock2-iclock6,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)

      num_iblank(nBody_Begin)=sum(INT(IBLANK(1:nxc,1:nyc,1:nzc)))
      WRITE(*,*) 'SUM of IBLANK = ', num_iblank(nBody_Begin)
      write(5670,*) nTime,num_iblank(nBody_Begin)
    END IF ! ntime

    ENDIF
    
    !CALL write_dump_debug_i('iblk',nBody_Begin,iblank)

!!    DO k=1,nzc
!    DO j=1,nyc
!    DO i=1,nxc
!	if (BodyNum(i,j,k)==1) then
!	print *, i,j,k
!end if
!enddo
!enddo
!enddo

   END SUBROUTINE set_iblank_canonical_body_fast

!------------------------------------------------------------------------------
   SUBROUTINE search_vertex_dotNorm(iBody,iCell,jCell,kCell,&
                                    closestElement,dotNorm)
!------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    use operation

    IMPLICIT NONE

!... parameters

    INTEGER, INTENT(IN)             :: iBody,iCell,jCell,kCell
    INTEGER, INTENT(OUT)            :: closestElement
    REAL(KIND=CGREAL) , INTENT(OUT) :: dotNorm

!... loop variables

    INTEGER :: m,n,ns,i
    
!... local variables
 
    INTEGER, PARAMETER       :: NSIZE = 1000
    INTEGER, PARAMETER       :: MSIZE = 20
    INTEGER,           DIMENSION(:),ALLOCATABLE   :: NeighElemInd
    REAL(KIND=CGREAL), DIMENSION(:),ALLOCATABLE   :: distMarkerS

    INTEGER                  :: iErr,nBodySelect,numNeighElement
    INTEGER                  :: elemInd,node1,node2,node3,nMarker
    INTEGER                  :: nSave
    INTEGER,DIMENSION(1)     :: iDummy(1)
    INTEGER                  :: shortestProbe
    INTEGER,DIMENSION(1:MSIZE) :: cElement,closestVert


    REAL(KIND=CGREAL)        :: xCell,yCell,zCell
    REAL(KIND=CGREAL)        :: dMin,dsIntercept,xM,yM,zM
    REAL(KIND=CGREAL)        :: areaDiffMin
    REAL(KIND=CGREAL)        :: distBIElem, distBIElemMin
    REAL(KIND=CGREAL)        :: planeConst,distanceToPlane,distPointToPlane
    REAL(KIND=CGREAL)        :: side12,side23,side31,side14,side24,side34
    REAL(KIND=CGREAL)        :: area123,area124,area234,area314
    REAL(KIND=CGREAL)        :: semiPerimeter123,semiPerimeter124,semiPerimeter234,semiPerimeter314
    REAL(KIND=CGREAL)        :: epsiArea,areaDiff
    REAL(KIND=CGREAL)        :: xBI, yBI, zBI
    REAL(KIND=CGREAL)        :: xBITemp, yBITemp, zBITemp
    REAL(KIND=CGREAL),DIMENSION(3) :: xVert, yVert, zVert
    REAL(KIND=CGREAL),DIMENSION(1:MSIZE) :: xBIT,yBIT,zBIT,dist,dist1
    REAL(KIND=CGREAL)        :: distInside,distIn
    logical :: inside,insideJudge
    
    REAL(KIND=CGREAL),DIMENSION(3) :: pCell,pMarker,eNorm,eRef,pOut,n1Coord,n2Coord,n3Coord
    REAL(KIND=CGREAL) :: distCell,distMarker,area1,area2,area3,areaSum,eArea,moMarker
    integer :: choice
    logical :: changeMarker,edgeMarker,ifPointIn

!******************************************************************************

    xCell = xc(iCell)
    yCell = yc(jCell)
    zCell = zc(kCell)
    
    !if(icell==60.and.jcell==103.and.kcell==63)then
    !    write(*,*) xcell,ycell,zcell
    !end if
    
    pCell(1)=xCell
    pCell(2)=yCell
    pCell(3)=zCell

    dotNorm = -oned

    nMarker = nPtsBodyMarker(iBody)

! nSave:  Number of closesest nodes to save for later use

    nSave  = 7

    IF (nSave > MSIZE) THEN
       PRINT*,'nSave in GCM_calc_bodyIntercept_Unstruc is limited to', MSIZE
       PRINT*,'Increase array size'
       STOP
    ENDIF

! ============================================================================
!   Allocate local array 
! ============================================================================

    ALLOCATE(distMarkerS(nMarker),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Allocation Error for distMarker'
      STOP
    ENDIF ! ierr

    ALLOCATE(NeighElemInd(NSIZE),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Allocation Error for NeighElemInd'
      STOP
    ENDIF ! ierr 

! ============================================================================
!   Get closestMarker for cell 
! ============================================================================

    dMin = 1.0E+16_CGREAL

    DO m = 1, nMarker
      xM = xBodyMarker(iBody,m)
      yM = yBodyMarker(iBody,m)
      zM = zBodyMarker(iBody,m)
      distMarkerS(m) = (xM-xCell)**2 + (yM-yCell)**2 + (zM-zCell)**2
    ENDDO ! m

    DO ns = 1,nSave
        iDummy                        = MINLOC(distMarkerS(1:nMarker))
        closestVert(ns)               = iDummy(1)
        distMarkerS(closestVert(ns))   = 1.0E20_CGREAL
    ENDDO


! ============================================================================
!   Find elements that share the three closest nodes/markers
! ============================================================================
    dist        = 1.0E+16_CGREAL
    dist1       = 1.0E+16_CGREAL
    insideJudge=.false.

    numNeighElement = 0
    DO m=1,totNumTriElem(iBody)
       IF ( triElemNeig(iBody,1,m) == closestVert(1) .OR. &
            triElemNeig(iBody,2,m) == closestVert(1) .OR. &
            triElemNeig(iBody,3,m) == closestVert(1) .OR. &
            triElemNeig(iBody,1,m) == closestVert(2) .OR. &
            triElemNeig(iBody,2,m) == closestVert(2) .OR. &
            triElemNeig(iBody,3,m) == closestVert(2) .OR. &
            triElemNeig(iBody,1,m) == closestVert(3) .OR. &
            triElemNeig(iBody,2,m) == closestVert(3) .OR. &
            triElemNeig(iBody,3,m) == closestVert(3) .OR. &
            triElemNeig(iBody,1,m) == closestVert(4) .OR. &
            triElemNeig(iBody,2,m) == closestVert(4) .OR. &
            triElemNeig(iBody,3,m) == closestVert(4) .OR. &
            triElemNeig(iBody,1,m) == closestVert(5) .OR. &
            triElemNeig(iBody,2,m) == closestVert(5) .OR. &
            triElemNeig(iBody,3,m) == closestVert(5) .OR. &
            triElemNeig(iBody,1,m) == closestVert(6) .OR. &
            triElemNeig(iBody,2,m) == closestVert(6) .OR. &
            triElemNeig(iBody,3,m) == closestVert(6) .OR. &
            triElemNeig(iBody,1,m) == closestVert(7) .OR. &
            triElemNeig(iBody,2,m) == closestVert(7) .OR. &
            triElemNeig(iBody,3,m) == closestVert(7)) THEN
          numNeighElement               = numNeighElement + 1
          NeighElemInd(numNeighElement) = m
       ENDIF
    ENDDO ! m

! ============================================================================
!   Trap error if array NeighElemenInd overflows 
! ============================================================================

    IF ( numNeighElement > NSIZE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Overflow Error for NeighElemInd'
      WRITE(STDOUT,*) ' Allocated size = ',NSIZE
      WRITE(STDOUT,*) ' Current size   = ',numNeighElement
      WRITE(STDOUT,*) ' Aborting Run'

      STOP
    ENDIF ! NeighElemInd
    
! ============================================================================
!   change the ref marker (closest marker) if the line segment "marker>cell" 
!   intersect with body surface
! ============================================================================
    choice=0
    do i=1,nSave
        changeMarker=.false.
        
        pMarker(1)=xBodyMarker(iBody,closestVert(i))
        pMarker(2)=yBodyMarker(iBody,closestVert(i))
        pMarker(3)=zBodyMarker(iBody,closestVert(i))
        
        do n=1,numNeighElement
            
            elemInd=NeighElemInd(n)
            
            node1=triElemNeig(iBody,1,elemInd)
            node2=triElemNeig(iBody,2,elemInd)
            node3=triElemNeig(iBody,3,elemInd)
            
            if(node1==closestVert(i).or.&
               node2==closestVert(i).or.&
               node3==closestVert(i))then
                cycle
            end if
            
            eNorm(1)=triElemNormx(iBody,elemInd)
            eNorm(2)=triElemNormy(iBody,elemInd)
            eNorm(3)=triElemNormz(iBody,elemInd)       
        
            eRef(1)=triElemCentx(iBody,elemInd)
            eRef(2)=triElemCenty(iBody,elemInd)
            eRef(3)=triElemCentz(iBody,elemInd)
        
            n1Coord(1)=xBodyMarker(iBody,node1)
            n1Coord(2)=yBodyMarker(iBody,node1)
            n1Coord(3)=zBodyMarker(iBody,node1)
            
            n2Coord(1)=xBodyMarker(iBody,node2)
            n2Coord(2)=yBodyMarker(iBody,node2)
            n2Coord(3)=zBodyMarker(iBody,node2)
            
            n3Coord(1)=xBodyMarker(iBody,node3)
            n3Coord(2)=yBodyMarker(iBody,node3)
            n3Coord(3)=zBodyMarker(iBody,node3)
                      
            call get_distance_face(eNorm,eRef,pCell,distCell)
            call get_distance_face(eNorm,eRef,pMarker,distMarker)
            
            if(distCell*distMarker>zero)then
                cycle
            else
                call get_inter(eNorm,eRef,pCell,pMarker,pOut)
            end if
            
            call pointInTriangle(n1Coord,n2Coord,n3Coord,pOut,ifPointIn)
            
            if(.not.ifPointIn)then
                cycle
            else
                changeMarker=.true.
                exit
            end if
            
            !call triangle_area(n1Coord,n2Coord,n3Coord,eArea)
            !
            !call triangle_area(n1Coord,n2Coord,pOut,area1)
            !call triangle_area(n2Coord,n3Coord,pOut,area2)
            !call triangle_area(n3Coord,n1Coord,pOut,area3)
            !
            !areaSum=area1+area2+area3
            !moMarker=mo(pMarker-pOut,3)
            !
            !if(abs(areaSum-eArea)>1e-30)then
            !    cycle
            !else if(moMarker<1e-30)then
            !    cycle
            !else
            !    changeMarker=.true.
            !    
            !    !if(icell==149.and.jcell==120.and.kcell==64)then
            !    !    write(*,*) changeMarker
            !    !pause
            !    !end if
            !    
            !    exit
            !end if
            
        end do
        
        if(.not.changeMarker)then
            choice=i
            exit
        else
            if(i==nSave)then
                write(*,*) 'Closest marker point is not found, you may need to incease "nSave" in "search_vertex_dotNorm"!'
                write(*,*) icell,jcell,kcell
                
                stop
                !dotNorm=-1.0d0
                !return
            end if
        end if
        
    end do
    
    call check_if_edge_marker(iBody,closestVert(choice),edgeMarker,iCell,jCell,kCell,dotNorm)
    if(edgeMarker)then
        dotNorm=-1.0d0
        return
    end if
    
    numNeighElement = 0
    DO m=1,totNumTriElem(iBody)
        IF ( triElemNeig(iBody,1,m) == closestVert(choice) .OR. &
            triElemNeig(iBody,2,m) == closestVert(choice) .OR. &
            triElemNeig(iBody,3,m) == closestVert(choice)) THEN
            numNeighElement               = numNeighElement + 1
            NeighElemInd(numNeighElement) = m
        ENDIF
    ENDDO ! m
    
    !if(choice/=1)then
    !    write(*,*) choice,icell,jcell,kcell
    !    write(*,*) xcell,ycell,zcell
    !    write(*,*) closestVert(1)
    !    write(*,*) '================'
    !end if
    
    
! ============================================================================
!   Check which element contains normal intercept
! ============================================================================

    distBIElemMin = 1.0E+16_CGREAL
    areaDiffMin   = 1.0E+16_CGREAL
    epsiArea      = 1.0E-4_CGREAL
    distIn        = 1.0E+16_CGREAL

    closestElement = 0
    
    inside=.false.
    DO n = 1,numNeighElement

     elemInd = NeighElemInd(n)

     node1   = triElemNeig(iBody,1,elemInd)
     node2   = triElemNeig(iBody,2,elemInd)
     node3   = triElemNeig(iBody,3,elemInd)

! ******************************************************************************
!    Check if BI inside the triangle of the surface element
!     through area differences
! ******************************************************************************
    !if(icell==149.and.jcell==120.and.kcell==64)then
    !    write(*,*) 'debug',xcell,ycell,zcell,choice
    !    write(*,*) closestVert(1),closestVert(choice)
    !    write(*,*) 'element:',elemInd
    !    write(*,*) triElemNormx(nBody,elemInd),triElemNormy(nBody,elemInd),triElemNormz(nBody,elemInd)
    !end if
     CALL check_BIInsideTriangle(iBody,elemInd,node1,node2,node3,xCell,yCell,zCell,&
                                 xBITemp,yBITemp,zBITemp,area123,areaDiff,distInside)

!!!!!!!!!!!!!!!!!!!!
!    IF (iCell == 1 .and. jCell == 1 .and. kCell == 1) then
!   IF (iCell == 42 .and. jCell == 47 .and. kCell == 91) then
!      WRITE(355,*)'ZONE'
!      WRITE(355,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!      WRITE(355,*)xBodyMarker(iBody,node2),yBodyMarker(iBody,node2),zBodyMarker(iBody,node2)
!      WRITE(355,*)xBodyMarker(iBody,node3),yBodyMarker(iBody,node3),zBodyMarker(iBody,node3)
!      WRITE(355,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!      WRITE(355,*)'ZONE'
!      WRITE(355,*)xCell,yCell,zCell
!      WRITE(355,*)xBItemp,yBItemp,zBItemp
!      WRITE(356,*) 'n,elemInd,areaDiff,area123 = ',n,elemInd,areadiff,area123
!   ENDIF
!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------------
!    Select closest Elem and BI coordinates:
!     If BI falls inside the element use that
!     Else Base the selection on the minimum distance 
!       between BI and either the norm to closest side or vertices of side
!--------------------------------------------------------------------------

     IF ( ABS(areaDiff) < epsiArea*area123) THEN
         
        inside=.true.
        insideJudge=insideJudge.or.inside
        if(distInside <= distIn)then
            distIn=distInside
        
            xBI = xBITemp
            yBI = yBITemp
            zBI = zBITemp
            closestElement = elemInd
            !dist(nc) = zero
            !GOTO 999
        end if
     ELSE if(.not.inside)then

        CALL calc_BIOutsideTriangle(iBody,elemInd,node1,node2,node3,xcell,ycell,zcell,closestVert(choice),  &
                                    xBITemp,yBITemp,zBITemp,distBIElem)
         
        
        !xBI = xBITemp
        !yBI = yBITemp
        !zBI = zBITemp
        !closestElement = elemInd
        !dist(nc) = zero
        !GOTO 999

!!!!!!!!!!!!!!!!!!!!
!    IF (iCell == 1 .and. jCell == 1 .and. kCell == 1) then
!   IF (iCell == 42 .and. jCell == 47 .and. kCell == 91) then
!       WRITE(356,*) 'n,distBIElem = ',n,distBIElem,distBIElemMin
!   ENDIF
!!!!!!!!!!!!!!!!!!!!!

        IF (distBIElem <= distBIElemMin) THEN
          distBIElemMin = distBIElem
          closestElement = elemInd
          xBI = xBITemp
          yBI = yBITemp
          zBI = zBITemp

!!!!!!!!!!!!!!!!!!!!
!    IF (iCell == 1 .and. jCell == 1 .and. kCell == 1) then
!   IF (iCell == 42 .and. jCell == 47 .and. kCell == 91) then
!       WRITE(356,*) 'n,closestElem-distBI = ',n,elemInd
!       WRITE(356,*) 'n,xyzBI-distBI = ',n,xBI,yBI,zBI
!   ENDIF
!!!!!!!!!!!!!!!!!!!!!

        ENDIF ! distBIElem

        !dist(nc) = distBIElem          !corrected by yan

     ENDIF ! areaDiff


    ENDDO ! n

999  CONTINUE

    xBIT(1) = xBI
    yBIT(1) = yBI
    zBIT(1) = zBI
    if(inside)then
        dist(1) = distIn         !corrected by yan
    else
        dist1(1) = distBIElemMin
    end if
    cElement(1) = closestElement



!DEBUG!!!!!!!!!!!!!
!   IF (iCell == 73 .and. jCell == 40 .and. kCell == 18) then
!    node1   = triElemNeig(iBody,1,closestElement)
!    node2   = triElemNeig(iBody,2,closestElement)
!    node3   = triElemNeig(iBody,3,closestElement)

!      WRITE(359,*)'ZONE'
!      WRITE(359,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!      WRITE(359,*)xBodyMarker(iBody,node2),yBodyMarker(iBody,node2),zBodyMarker(iBody,node2)
!      WRITE(359,*)xBodyMarker(iBody,node3),yBodyMarker(iBody,node3),zBodyMarker(iBody,node3)
!      WRITE(359,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!      WRITE(359,*)'ZONE'
!      WRITE(359,*)xCell,yCell,zCell
!      WRITE(359,*)xBI,yBI,zBI
!   ENDIF
!!!!!!!!!!!!!!!!!!!!!

!DEBUG!!!!!!!!!!!!


    if(insideJudge)then
        iDummy           = MINLOC(dist(1:1))
        shortestProbe    = iDummy(1)
        closestElement   = cElement(shortestProbe)
    else
        iDummy           = MINLOC(dist1(1:1))
        shortestProbe    = iDummy(1)
        closestElement   = cElement(shortestProbe)
    end if

!DEBUG!!!!!!!!!!!!!
!   IF (iCell == 73 .and. jCell == 40 .and. kCell == 18) then
!    node1   = triElemNeig(iBody,1,closestElement)
!    node2   = triElemNeig(iBody,2,closestElement)
!    node3   = triElemNeig(iBody,3,closestElement)

!      WRITE(359,*)'ZONE'
!      WRITE(359,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!      WRITE(359,*)xBodyMarker(iBody,node2),yBodyMarker(iBody,node2),zBodyMarker(iBody,node2)
!      WRITE(359,*)xBodyMarker(iBody,node3),yBodyMarker(iBody,node3),zBodyMarker(iBody,node3)
!      WRITE(359,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!      WRITE(359,*)'ZONE'
!      WRITE(359,*)xCell,yCell,zCell
!      WRITE(359,*)xBIT(shortestProbe),yBIT(shortestProbe),zBIT(shortestProbe)
!   ENDIF
!!!!!!!!!!!!!!!!!!!!!



! ============================================================================
!   Perform the dot product with element that has shortest distance or area.
! ============================================================================

    dotNorm = (xCell - triElemCentx(iBody,closestElement)) &
                      *triElemNormx(iBody,closestElement)  &
             +(yCell - triElemCenty(iBody,closestElement)) &
                      *triElemNormy(iBody,closestElement)  &
             +(zCell - triElemCentz(iBody,closestElement)) &
                      *triElemNormz(iBody,closestElement)

!!!!!!!!!!!!!!!!!!!!!
!    IF (iCell == 1 .and. jCell == 1 .and. kCell == 1) then
!   IF (iCell == 40 .and. jCell == 46 .and. kCell == 89) then
!WRITE(356,*) ' DotNorm = ' ,dotnorm
!       WRITE(356,*) ' xyzBI   = ' ,xBI,yBI,zBI
!WRITE(356,*) ' ClosestElement = ' ,ClosestElement
!WRITE(356,*) ' triElemNormxyz = ' ,triElemNormx(iBody,closestElement),&
!                                          triElemNormy(iBody,closestElement),&
!                                          triElemNormz(iBody,closestElement)
!       WRITE(355,*)'ZONE'
!       WRITE(355,*)xCell,yCell,zCell
!       WRITE(355,*)xBI,yBI,zBI
!   ENDIF
!!!!!!!!!!!!!!!!!!!!!

! ============================================================================
!   Deallocate local array 
! ============================================================================

    DEALLOCATE(NeighElemInd,STAT=iErr)
    DEALLOCATE(distMarkerS, STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Deallocation Error for NeighElemInd'
      STOP
    ENDIF ! ierr 

   END SUBROUTINE search_vertex_dotNorm
!------------------------------------------------------------------------------
                            
!------------------------------------------------------------------------------
   SUBROUTINE check_BIInsideTriangle(iBody,elemInd,node1,node2,node3,xCell,yCell,zCell,&
                                     xBITemp,yBITemp,zBITemp,area123,areaDiff,distInside)
!------------------------------------------------------------------------------
           
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... parameters

    INTEGER, INTENT(IN) :: elemInd,iBody,node1,node2,node3
    REAL(KIND=CGREAL), INTENT(IN)  :: xCell,yCell,zCell
    REAL(KIND=CGREAL), INTENT(OUT) :: xBITemp,yBITemp,zBITemp,area123,areaDiff,distInside

!... loop variables

    INTEGER :: iside

!... local variables

    REAL(KIND=CGREAL)        :: planeConst,distanceToPlane,distPointToPlane
    REAL(KIND=CGREAL)        :: side12,side23,side31,side14,side24,side34
    REAL(KIND=CGREAL)        :: area124,area234,area314
    REAL(KIND=CGREAL)        :: semiPerimeter123,semiPerimeter124
    REAL(KIND=CGREAL)        :: semiPerimeter234,semiPerimeter314

! ******************************************************************************
! equation of plane (note our normals are unit normals)
     
!  n  x + n  y + n  z + planeConst = 0
!   x         y      z
! ******************************************************************************

    planeConst =- triElemNormx(iBody,elemInd)*xBodyMarker(iBody,node1) &
                - triElemNormy(iBody,elemInd)*yBodyMarker(iBody,node1) &
                - triElemNormz(iBody,elemInd)*zBodyMarker(iBody,node1)
       
! ******************************************************************************
! Compute coordinates of normal intercept
!      
! Consider point Po = (xo,yo,zo)  not on plane
!   and  point   P1 = (x1,y1,z1)  on the plane
!                                                               ^
! equation of line through Po normal to plane is  P(s) = Po + s n
!      
! normal distance from Po to Plane is given by
!            ^                               ^         ^ ^
!            n. ( P(s) - P1 ) = 0   => so = -n.(Po-P1)/n.n                     ^ ^
!                                         = -(n xo + n yo + n zo + planeConst)/n.n
!                                              x      y      z
!
!                                                 ^
!   subsequently normal intersection point = Po + so n  
!                   
! ******************************************************************************

     distanceToPlane = -(  triElemNormx(iBody,elemInd)*xCell  &
                         + triElemNormy(iBody,elemInd)*yCell  &
                         + triElemNormz(iBody,elemInd)*zCell  &
                         + planeConst )
     distInside=abs(distanceToPlane)

     xBITemp = xCell + triElemNormx(iBody,elemInd)*distanceToPlane
     yBITemp = yCell + triElemNormy(iBody,elemInd)*distanceToPlane
     zBITemp = zCell + triElemNormz(iBody,elemInd)*distanceToPlane
! ******************************************************************************
! Check to see if normal intercept lies inside the closest trianglular element
!               3 
!               *  .
!              /  \   .
!             /    \    .
!            /      \    * 4=BI
!           /        \  .
!         1*__________*2
!         
! Basic Idea :  IF [ AREA(124) + AREA(234) + AREA(314) ] > AREA(123) THEN  POINT(4) is
! outside triangle (123)
!
! using Heron formula for area of triangle
! AREA(123) = SQRT[ S * ( S - S12) * (S - S23) * (S - S31) ]
! S = 0.5*(S12 + S23 + S31) 
! ******************************************************************************

     side12 =  SQRT( (xBodyMarker(iBody,node2)-xBodyMarker(iBody,node1))**2  &
                    +(yBodyMarker(iBody,node2)-yBodyMarker(iBody,node1))**2  &
                    +(zBodyMarker(iBody,node2)-zBodyMarker(iBody,node1))**2  )
     side23 =  SQRT( (xBodyMarker(iBody,node3)-xBodyMarker(iBody,node2))**2  &
                    +(yBodyMarker(iBody,node3)-yBodyMarker(iBody,node2))**2  &
                    +(zBodyMarker(iBody,node3)-zBodyMarker(iBody,node2))**2  )
     side31 =  SQRT( (xBodyMarker(iBody,node1)-xBodyMarker(iBody,node3))**2  &
                    +(yBodyMarker(iBody,node1)-yBodyMarker(iBody,node3))**2  &
                    +(zBodyMarker(iBody,node1)-zBodyMarker(iBody,node3))**2  )
     side14 =  SQRT( (xBITemp-xBodyMarker(iBody,node1))**2  &
                    +(yBITemp-yBodyMarker(iBody,node1))**2  &
                    +(zBITemp-zBodyMarker(iBody,node1))**2  )
     side24 =  SQRT( (xBITemp-xBodyMarker(iBody,node2))**2  &
                    +(yBITemp-yBodyMarker(iBody,node2))**2  &
                    +(zBITemp-zBodyMarker(iBody,node2))**2  )
     side34 =  SQRT( (xBITemp-xBodyMarker(iBody,node3))**2  &
                    +(yBITemp-yBodyMarker(iBody,node3))**2  &
                    +(zBITemp-zBodyMarker(iBody,node3))**2  )

     semiPerimeter123 = half*(side12 + side23 + side31)
     semiPerimeter124 = half*(side12 + side24 + side14)
     semiPerimeter234 = half*(side23 + side24 + side34)
     semiPerimeter314 = half*(side31 + side34 + side14)

     area123       = SQRT( semiPerimeter123*(semiPerimeter123-side12) &
                                           *(semiPerimeter123-side23) &
                                           *(semiPerimeter123-side31)   )
   
     area124       = SQRT( semiPerimeter124*(semiPerimeter124-side12) &
                                           *(semiPerimeter124-side24) &
                                           *(semiPerimeter124-side14)   )
    
     area234       = SQRT( semiPerimeter234*(semiPerimeter234-side23) &
                                           *(semiPerimeter234-side24) &
                                           *(semiPerimeter234-side34)   )

     area314       = SQRT( semiPerimeter314*(semiPerimeter314-side31) &
                                           *(semiPerimeter314-side34) &
                                           *(semiPerimeter314-side14)   )

     areaDiff  = area124 + area234 + area314 - area123

   END SUBROUTINE check_BIInsideTriangle 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE calc_BIOutsideTriangle(iBody,elemInd,node1,node2,node3,xcell,ycell,zcell,cMarker,  &
                                     xBITemp,yBITemp,zBITemp,distBIElem)
!------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    use operation

    IMPLICIT NONE

!... parameters

    INTEGER, INTENT(IN)  :: elemInd,iBody,node1,node2,node3,cMarker
    REAL(KIND=CGREAL), INTENT(IN)  :: xBITemp,yBITemp,zBITemp
    REAL(KIND=CGREAL), INTENT(IN)  :: xCell,yCell,zCell
    REAL(KIND=CGREAL), INTENT(OUT) :: distBIElem

!... loop variables

    INTEGER :: iside

!... local variables

    INTEGER :: isideSelect,mside
    INTEGER :: nodeSelect1,nodeSelect2

    REAL(KIND=CGREAL) :: aCrossbVectMagn,dotVal
    REAL(KIND=CGREAL) :: distIntBI,distIntCG,distIntMin
    REAL(KIND=CGREAL) :: distNorm,distNode1BINorm,distNode2BINorm
    REAL(KIND=CGREAL) :: distVert1BI,distVert2BI
    REAL(KIND=CGREAL) :: magnitude12,magnitudeBICG,magnitude,projectedLength
    REAL(KIND=CGREAL) :: node12x,node12y,node12z
    REAL(KIND=CGREAL) :: vec01x,vec01y,vec01z
    REAL(KIND=CGREAL) :: vec12x,vec12y,vec12z
    REAL(KIND=CGREAL) :: xBINorm,yBINorm,zBINorm
    REAL(KIND=CGREAL) :: xCG,yCG,zCG
    REAL(KIND=CGREAL), DIMENSION(3) :: aVect,bVect,cVect,aCrossbVect,cCrossbVect
    REAL(KIND=CGREAL), DIMENSION(3) :: xInt,yInt,zInt
    REAL(KIND=CGREAL), DIMENSION(3) :: xVert,yVert,zVert
    REAL(KIND=CGREAL), DIMENSION(3) :: vect1,vect2,vect3,vect4,norm
    REAL(KIND=CGREAL), DIMENSION(3,3) :: vectInt
    
    REAL(KIND=CGREAL) :: xClose,yClose,zClose,angleTest
    REAL(KIND=CGREAL), DIMENSION(3) :: vecMarker

!******************************************************************************
          
    xVert(1)  = xBodyMarker(iBody,node1)
    yVert(1)  = yBodyMarker(iBody,node1)
    zVert(1)  = zBodyMarker(iBody,node1)

    xVert(2)  = xBodyMarker(iBody,node2)
    yVert(2)  = yBodyMarker(iBody,node2)
    zVert(2)  = zBodyMarker(iBody,node2)

    xVert(3)  = xBodyMarker(iBody,node3)
    yVert(3)  = yBodyMarker(iBody,node3)
    zVert(3)  = zBodyMarker(iBody,node3)
    
    xClose=xBodyMarker(iBody,cMarker)
    yClose=yBodyMarker(iBody,cMarker)
    zClose=zBodyMarker(iBody,cMarker)

    xCG       = triElemCentx(iBody,elemInd)
    yCG       = triElemCenty(iBody,elemInd)
    zCG       = triElemCentz(iBody,elemInd)
    norm      = (/triElemNormx(iBody,elemInd),triElemNormy(iBody,elemInd),triElemNormz(iBody,elemInd)/)
!========================================================================================
!Modified algorithm for finding the closest element when BI outside triangle: add by Yan 
!========================================================================================
    vecMarker=(/xcell,ycell,zcell/)-(/xClose,yClose,zClose/)
    !vecMarker=(/xcell,ycell,zcell/)-(/xCG,yCG,zCG/)
    norm=-norm
    
    angleTest=vector_angle(vecMarker,norm)
    
    distBIElem=angleTest
    
    
    
    
!! ============================================================================
!!   Construct Intersection points between 
!!     line linking BI and Centroid of Surface Element and the triangle sides
!!     L1: BI-CG, L2:Sides of Vertices
!!
!!   use formula for intersection point between 2 co-planar lines from Mathworld
!!   http://mathworld.wolfram.com/Line-LineIntersect.html
!!
!!             x4
!!             *
!!             |
!!             |
!!    x1       | Int      x2
!!    *--------*----------*
!!             |
!!             |  
!!             | 
!!             * x3
!!
!! ============================================================================
!   
!    vect1(1:3) = (/xBITemp,yBITemp,zBITemp/)
!    vect2(1:3) = (/xCG,yCG,zCG/)
!
!    aVect(1:3) = vect2(1:3)-vect1(1:3)
!
!    DO iside = 1,3
!      mside = iside +1
!      IF(iside == 3) mside = 1
!
!      vect3(1:3) =(/xVert(iside),yVert(iside),zVert(iside)/)
!      vect4(1:3) =(/xVert(mside),yVert(mside),zVert(mside)/)
!
!      bVect(1:3)  = vect4(1:3) -vect3(1:3)
!      cVect(1:3)  = vect3(1:3) -vect1(1:3)
!
!      call calc_crossProduct(aVect,bVect,aCrossbVect)
!      call calc_crossProduct(cVect,bVect,cCrossbVect)
!
!      aCrossbVectMagn = aCrossbVect(1)**2 + aCrossbVect(2)**2 +aCrossbVect(3)**2
!      dotVal = DOT_PRODUCT(cCrossbVect,aCrossbVect)
!
!      vectInt(1:3,iside) = vect1(1:3) + aVect(1:3) *dotVal/aCrossbVectMagn
!    END DO ! iside 
!      
!! ============================================================================
!!   Choose closest intersection point lying between BI and CG
!!     Normalsize value with L1
!! ============================================================================
!
!    magnitudeBICG = SQRT( (vect1(1)-vect2(1))**2 &
!                        + (vect1(2)-vect2(2))**2 &
!                        + (vect1(3)-vect2(3))**2 )
!
!    distIntMin = 1.0E+16_CGREAL
!    isideSelect = -1000
!
!    DO iside = 1,3
!      distIntBI = SQRT( (vect1(1)-vectInt(1,iside))**2 &
!                      + (vect1(2)-vectInt(2,iside))**2 &
!                      + (vect1(3)-vectInt(3,iside))**2 )/magnitudeBICG
!
!      distIntCG = SQRT( (vect2(1)-vectInt(1,iside))**2 &
!                      + (vect2(2)-vectInt(2,iside))**2 &
!                      + (vect2(3)-vectInt(3,iside))**2 )/magnitudeBICG
!    !write(*,*) 'IntBi,IntCG',iside,xBItemp,yBItemp,zBItemp
!      !IF ( distIntBI <= oned .AND. distIntCG <= oned) THEN
!      IF ( abs(distIntBI+distIntCG-oned) <= 1e-3) THEN
!        isideSelect = iside
!      END IF ! distIntBI
!    END DO ! iside 
!
!! ============================================================================
!!   Trap error for isideSelect 
!! ============================================================================
!   
!     IF ( isideSelect < 0 ) THEN
!      WRITE(STDOUT,*) &
!       'calc_BIdistMin: Incorrect selection of iside (Should be either 1, 2 or 3'
!      WRITE(STDOUT,*) &
!       '                default value selected = ',isideSelect
!      STOP
!     END IF ! isideSelect 
!
!! ============================================================================
!!   Select appropriate vertices from isideSelect 
!! ============================================================================
!   
!     SELECT CASE(isideSelect)
!       CASE(1)
!         nodeSelect1 = 1
!         nodeSelect2 = 2
!       CASE(2)
!         nodeSelect1 = 2
!         nodeSelect2 = 3
!       CASE(3)
!         nodeSelect1 = 3
!         nodeSelect2 = 1
!     END SELECT ! isideSelect
!
!! ============================================================================
!!   Drop normals from BI to selected side 
!!    and find coordinates of intersection point
!!
!!   unit vector from node 1 to node 2
!!   use formula for distance between point and line from Mathworld
!!   http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
!! 
!!    x1                  x2
!!    *-------------------*
!!             |                              |(x2-x1) x (x1-x0)|
!!             | d                       d =  --------------------
!!             |                                   |x2-x1|
!!             * x0
!!   x0: BI
!! ============================================================================
!
!     vec12x = xVert(nodeSelect2) - xVert(nodeSelect1)
!     vec12y = yVert(nodeSelect2) - yVert(nodeSelect1)
!     vec12z = zVert(nodeSelect2) - zVert(nodeSelect1)
!
!     magnitude12 = SQRT(vec12x**2 + vec12y**2 + vec12z**2)
!
!     vec01x = xVert(nodeSelect1) - xBITemp
!     vec01y = yVert(nodeSelect1) - yBITemp
!     vec01z = zVert(nodeSelect1) - zBITemp
!     
!     distNorm = SQRT(  (vec12y*vec01z - vec12z*vec01y)**2  &
!                     + (vec12z*vec01x - vec01z*vec12x)**2  &
!                     + (vec12x*vec01y - vec01x*vec12y)**2  )/magnitude12
!
!! ============================================================================
!!    Project vector BI-node1 onto node12 to find body intercept point
!! ============================================================================
!
!     node12x = xVert(nodeSelect2) - xVert(nodeSelect1)
!     node12y = yVert(nodeSelect2) - yVert(nodeSelect1)
!     node12z = zVert(nodeSelect2) - zVert(nodeSelect1)
!
!     magnitude = SQRT(node12x**2 + node12y**2 + node12z**2)
!
!     node12x = node12x/magnitude
!     node12y = node12y/magnitude
!     node12z = node12z/magnitude
!
!     projectedLength = (xBITemp - xVert(nodeSelect1))*node12x  &
!                      +(yBITemp - yVert(nodeSelect1))*node12y  &
!                      +(zBITemp - zVert(nodeSelect1))*node12z
!
!     xBINorm = xVert(nodeSelect1) + projectedLength*node12x
!     yBINorm = yVert(nodeSelect1) + projectedLength*node12y
!     zBINorm = zVert(nodeSelect1) + projectedLength*node12z
!
!! ============================================================================
!!    Determine distance between BINorm and vertices of selected side.
!!     If normal point lies inside the side, select that distance.
!!     If it lies outside, find the minimum distance with vertices
!!     Use normalized length
!! ============================================================================
!
!        !distBIElem  =  distNorm
!     
!      distNode1BINorm = SQRT( (xVert(nodeSelect1)-xBINorm)**2 &
!                            + (yVert(nodeSelect1)-yBINorm)**2 &
!                            + (zVert(nodeSelect1)-zBINorm)**2 )/magnitude
!      
!      !distNode1BINorm = SQRT( (xVert(nodeSelect2)-xBINorm)**2 &    !previous
!      distNode2BINorm = SQRT( (xVert(nodeSelect2)-xBINorm)**2 &     !changed by Yan
!                            + (yVert(nodeSelect2)-yBINorm)**2 &
!                            + (zVert(nodeSelect2)-zBINorm)**2 )/magnitude
!      
!      !IF ( distNode1BINorm <= oned) THEN                           !previous
!      IF ( distNode1BINorm + distNode2BINorm - oned <= 1e-3) THEN   !changed by Yan
!        distBIElem  =  distNorm
!      ELSE
!        distVert1BI = SQRT( (xVert(nodeSelect1)-xBITemp)**2 &
!                          + (yVert(nodeSelect1)-yBITemp)**2 &
!                          + (zVert(nodeSelect1)-zBITemp)**2 )
!        
!        distVert2BI = SQRT( (xVert(nodeSelect2)-xBITemp)**2 &
!                          + (yVert(nodeSelect2)-yBITemp)**2 &
!                          + (zVert(nodeSelect2)-zBITemp)**2 )
!        
!        distBIElem = DMIN1(distVert1BI,distVert2BI)
!      END IF ! distNode1BINorm

   END SUBROUTINE calc_BIOutsideTriangle 
!------------------------------------------------------------------------------

   SUBROUTINE calc_crossProduct(r,s,cross_product)
   
    USE global_parameters
    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(3), INTENT(IN)  :: r,s
    REAL(KIND=CGREAL), DIMENSION(3), INTENT(OUT) :: cross_product 

    INTEGER :: component,i,j

    DO component = 1,3
      i = MODULO(component,3) + 1
      j = MODULO(i,3) + 1
      cross_product(component) = r(i)*s(j) - s(i)*r(j)
    END DO ! component 

   END SUBROUTINE calc_crossProduct
!------------------------------------------------------------------------------

SUBROUTINE calc_area_n_volumn

  USE global_parameters
  USE flow_parameters
  USE grid_arrays
  USE flow_arrays
  USE boundary_arrays
  USE GCM_arrays
  USE unstructured_surface_arrays

  implicit none

  integer :: i, j, k, cutting
  REAL(KIND=CGREAL) :: CALC_CELL_VOLUMN, TMP, cut_area
  integer :: im, ip, jm, jp, km, kp
  INTEGER :: CALC_CUT_FACE
  
  cutting=sum(abs(int(ghostCellSolid)))
  IF (ASSOCIATED(CELL)) THEN
	DEALLOCATE(CELL)
    deallocate(face)
    deallocate(edge_CUT)
  ENDIF
    
  allocate(CELL(cutting*8))
!  allocate(point_pair(cutting*8))
  allocate(face(-1:cutting*7*8))
  allocate(edge_CUT(cutting*12*8))
  
  allocate(iblank_in(0:nx+1,0:ny+1,0:nz+1)) ! iblank for cell-centered, iblank_in for vectex_centered
  allocate(iex(0:nx+1,0:ny+1,0:nz+1),iey(0:nx+1,0:ny+1,0:nz+1),iez(0:nx+1,0:ny+1,0:nz+1))
  allocate(ifx(0:nx+1,0:ny+1,0:nz+1),ify(0:nx+1,0:ny+1,0:nz+1),ifz(0:nx+1,0:ny+1,0:nz+1))
  
  CUTTING_EDGE_N=0
  CUTTING_FACE_N=0
  CUTTING_CELL_N=0

  FACE(-1)%A=zero
  FACE( 0)%A=oned
  
  ivc=0
  iblank_in=0
  
  iex=-1 ! -1: not tested; 0: tested, no interception; >0: tested, has interception.
  iey=-1
  iez=-1
  
  ifx=-1
  ify=-1
  ifz=-1

  CALL FIND_EDGES_INTERCEPTING_TRIANGLE() ! now can deal with the case where one and only one intercepted point on each line segment.
  
  DO k = 1, nzc
  DO j = 1, nyc
  DO i = 1, nxc
    iM = MAX(i-1,1)
    iP = MIN(i+1,nxc)
    jM = MAX(j-1,1)
    jP = MIN(j+1,nyc)
    kM = MAX(k-1,1)
    kP = MIN(k+1,nzc)
    IF (ANY(iblank_solid(im:ip,jm:jp,km:kp)>0) .AND. ANY(iblank_solid(im:ip,jm:jp,km:kp)==0)) THEN
!  IF ( ghostCellSolid(i,j,k)/=0 .OR. ghostCellSolid(i+1,j,k)/=0 .OR. ghostCellSolid(i,j,k+1)/=0 .OR. ghostCellSolid(i,j+1,k)/=0 )  THEN

    IF (ifz(I,J,K+1)<0)  CALL check_front_face(i,j,k)
    IF (ifx(I+1,J,K)<0)  CALL check_right_face(i,j,k)
    IF (ify(I,J+1,K)<0)  CALL check_up_face(i,j,k)

    IF (ifz(I,J,K)<0)  CALL check_back_face(i,j,k)
    IF (ifx(I,J,K)<0)  CALL check_left_face(i,j,k)
    IF (ify(I,J,K)<0)  CALL check_low_face(i,j,k)
  ENDIF

  ENDDO ! i
  ENDDO ! j
  ENDDO ! k
  
  DO k = 1, nzc
  DO j = 1, nyc
  DO i = 1, nxc
  IF ( ifx(i,j,k)>0 .or. ifx(i+1,j,k)>0 .or. &
       ify(i,j,k)>0 .or. ify(i,j+1,k)>0 .or. &
       ifz(i,j,k)>0 .or. ifz(i,j,k+1)>0)  THEN
    
	CUTTING_CELL_N=CUTTING_CELL_N+1
	ivc(i,j,k)=CUTTING_CELL_N
	CELL(CUTTING_CELL_N)%I=I
	CELL(CUTTING_CELL_N)%J=J
	CELL(CUTTING_CELL_N)%k=k

	CELL(CUTTING_CELL_N)%F_im=ifx(i,j,k)
	CELL(CUTTING_CELL_N)%F_ip=ifx(i+1,j,k)
	CELL(CUTTING_CELL_N)%F_jm=ify(i,j,k)
	CELL(CUTTING_CELL_N)%F_jp=ify(i,j+1,k)
	CELL(CUTTING_CELL_N)%F_km=ifz(i,j,k)
	CELL(CUTTING_CELL_N)%F_kp=ifz(i,j,k+1)

    CELL(CUTTING_CELL_N)%F_slice=CALC_CUT_FACE(I,J,K,CELL(CUTTING_CELL_N)%slice_normal)
    
	CELL(CUTTING_CELL_N)%VOLUMN=CALC_CELL_VOLUMN(I,J,K)
!	CELL(CUTTING_CELL_N)%volumn_fraction=CELL(CUTTING_CELL_N)%VOLUMN*dxinv(i)*dyinv(j)*dzinv(k)
  ENDIF

  ENDDO ! i
  ENDDO ! j
  ENDDO ! k

  deallocate(iblank_in)
  deallocate(iex,iey,iez)
  deallocate(ifx,ify,ifz)

END SUBROUTINE


SUBROUTINE check_front_face(I,J,K)

  USE global_parameters
  USE flow_parameters
  USE grid_arrays
  USE flow_arrays
  USE boundary_arrays
  USE GCM_arrays
  USE unstructured_surface_arrays

  implicit none

  integer :: i, j, k
  INTEGER(4) :: FIND_INTERSECTION
  REAL(KIND=CGREAL) :: CALC_CUTTING_AREA, TMP, CALC_CUTTING_AREA_NORMALIZED
  REAL(KIND=CGREAL) :: cx, cy
  
!  IF (iex(i,  j,k+1)<0) iex(i,  j,k+1) = FIND_INTERSECTION(I,  J,K+1,5)
!  IF (iex(i,j+1,k+1)<0) iex(i,j+1,k+1) = FIND_INTERSECTION(I,J+1,K+1,6)
!  IF (iey(i,  j,k+1)<0) iey(i,  j,k+1) = FIND_INTERSECTION(I,  J,K+1,7)
!  IF (iey(i+1,j,K+1)<0) iey(i+1,j,K+1) = FIND_INTERSECTION(I+1,J,K+1,8)
  INTEGER :: TM(4)

  tm(1)=iex(i  ,j  ,k+1) ! 5
  tm(2)=iey(i+1,j  ,k+1) ! 8
  tm(3)=iex(i  ,j+1,k+1) ! 6
  tm(4)=iey(i  ,j  ,k+1) ! 7
!  TMP=CALC_CUTTING_AREA(cx, cy, tm, &
!    xc(i), yc(j), x(i), x(i+1), y(j), y(j+1), dx(i), dy(j), &
!    iblank_in(i+1,j,k+1), iblank_in(i+1,j+1,k+1), iblank_in(i,j+1,k+1), iblank_in(i,j,k+1))
  TMP=CALC_CUTTING_AREA_NORMALIZED(cx, cy, tm, &
    xc(i), yc(j), x(i), x(i+1), y(j), y(j+1), dxinv(i), dyinv(j), &
    iblank_in(i+1,j,k+1), iblank_in(i+1,j+1,k+1), iblank_in(i,j+1,k+1), iblank_in(i,j,k+1))

!  TMP=CALC_CUTTING_AREA_Z(I,J,K+1,centroid)
  IF (TMP>ZERO) THEN
	CUTTING_FACE_N=CUTTING_FACE_N+1
	ifZ(I,J,K+1)=CUTTING_FACE_N
!	face(CUTTING_FACE_N)%area=TMP
	face(CUTTING_FACE_N)%a=TMP !/(dx(i)*dy(j))
	face(CUTTING_FACE_N)%centroid=(/cx,cy,z(k+1)/)
  ELSE
	ifZ(I,J,K+1)=0
  ENDIF
END SUBROUTINE


SUBROUTINE check_right_face(I,J,K)

  USE global_parameters
  USE flow_parameters
  USE grid_arrays
  USE flow_arrays
  USE boundary_arrays
  USE GCM_arrays
  USE unstructured_surface_arrays

  implicit none

  integer :: i, j, k
  INTEGER(4) :: FIND_INTERSECTION
  REAL(KIND=CGREAL) :: CALC_CUTTING_AREA, TMP, CALC_CUTTING_AREA_NORMALIZED
  REAL(KIND=CGREAL) :: cx, cy

!  IF (iez(i+1,  j,  k)<0) iez(i+1,  j,  k) = FIND_INTERSECTION(I+1,  J,K  ,11)
!  IF (iez(i+1,j+1,  K)<0) iez(i+1,j+1,  K) = FIND_INTERSECTION(I+1,J+1,K  ,12)
!  IF (iey(i+1,  j,  k)<0) iey(i+1,  j,  k) = FIND_INTERSECTION(I+1,  J,K  ,4)
!  IF (iey(i+1,  j,k+1)<0) iey(i+1,  j,k+1) = FIND_INTERSECTION(I+1,  J,K+1,8)
  INTEGER :: TM(4)

  tm(1)=iez(i+1,j  ,k)   ! 11
  tm(2)=iey(i+1,j  ,k+1) ! 8
  tm(3)=iez(i+1,j+1,k)   ! 12
  tm(4)=iey(i+1,j  ,k)   ! 4
!  TMP=CALC_CUTTING_AREA(cx, cy, tm, &
!    zc(k), yc(j), z(k), z(k+1), y(j), y(j+1), dz(k), dy(j), &
!    iblank_in(i+1,j,k+1), iblank_in(i+1,j+1,k+1), iblank_in(i+1,j+1,k), iblank_in(i+1,j,k))
  TMP=CALC_CUTTING_AREA_NORMALIZED(cx, cy, tm, &
    zc(k), yc(j), z(k), z(k+1), y(j), y(j+1), dzinv(k), dyinv(j), &
    iblank_in(i+1,j,k+1), iblank_in(i+1,j+1,k+1), iblank_in(i+1,j+1,k), iblank_in(i+1,j,k))
  
!  TMP=CALC_CUTTING_AREA_X(I+1,J,K,centroid)
  IF (TMP>ZERO) THEN
	CUTTING_FACE_N=CUTTING_FACE_N+1
	ifx(I+1,J,K)=CUTTING_FACE_N
!	face(CUTTING_FACE_N)%area=TMP
	face(CUTTING_FACE_N)%a   =TMP !/(dy(j)*dz(k))
	face(CUTTING_FACE_N)%centroid=(/x(i+1),cy,cx/)
  ELSE
    ifx(I+1,J,K)=0
  ENDIF
END SUBROUTINE

SUBROUTINE check_up_face(I,J,K)

  USE global_parameters
  USE flow_parameters
  USE grid_arrays
  USE flow_arrays
  USE boundary_arrays
  USE GCM_arrays
  USE unstructured_surface_arrays

  implicit none

  integer :: i, j, k
  INTEGER(4) :: FIND_INTERSECTION
  REAL(KIND=CGREAL) :: CALC_CUTTING_AREA, TMP, CALC_CUTTING_AREA_NORMALIZED
  REAL(KIND=CGREAL) :: cx, cy

!  IF (iex(i  ,j+1,  k)<0) iex(i  ,j+1,  k) = FIND_INTERSECTION(I  ,J+1,  K,2)
!  IF (iex(i  ,j+1,K+1)<0) iex(i  ,j+1,K+1) = FIND_INTERSECTION(I  ,J+1,K+1,6)
!  IF (iez(i  ,j+1,  k)<0) iez(i  ,j+1,  k) = FIND_INTERSECTION(I  ,J+1,  K,10)
!  IF (iez(i+1,j+1,  k)<0) iez(i+1,j+1,  k) = FIND_INTERSECTION(I+1,J+1,  K,12)
  INTEGER :: TM(4)

  tm(1)=iex(i  ,j+1,k)   ! 2
  tm(2)=iez(i+1,j+1,k)   ! 12
  tm(3)=iex(i  ,j+1,k+1) ! 6
  tm(4)=iez(i  ,j+1,k)   ! 10

!  TMP=CALC_CUTTING_AREA(cx, cy, tm, &
!    xc(i), zc(k), x(i), x(i+1), z(k), z(k+1), dx(i), dz(k), &
!    iblank_in(i+1,j+1,k), iblank_in(i+1,j+1,k+1), iblank_in(i,j+1,k+1), iblank_in(i,j+1,k))
!  if (ntime==83 .and. i==802 .and. j==92) then
!    write(*,*) i,j,k
!    write(*,*) tm
!    write(*,*) xc(i), zc(k)
!    write(*,*) x(i), x(i+1) 
!    write(*,*) z(k), z(k+1) 
!    write(*,*) dxinv(i), dzinv(k)
!    write(*,*) iblank_in(i+1,j+1,k), iblank_in(i+1,j+1,k+1), iblank_in(i,j+1,k+1), iblank_in(i,j+1,k)
!    write(*,*)
!  endif  
  TMP=CALC_CUTTING_AREA_NORMALIZED(cx, cy, tm, &
    xc(i), zc(k), x(i), x(i+1), z(k), z(k+1), dxinv(i), dzinv(k), &
    iblank_in(i+1,j+1,k), iblank_in(i+1,j+1,k+1), iblank_in(i,j+1,k+1), iblank_in(i,j+1,k))

!  TMP=CALC_CUTTING_AREA_Y(I,J+1,K,centroid)
  IF (TMP>ZERO) THEN
	CUTTING_FACE_N=CUTTING_FACE_N+1
	ifY(I,J+1,K)=CUTTING_FACE_N
!	face(CUTTING_FACE_N)%area=TMP
	face(CUTTING_FACE_N)%a   =TMP !/(dx(i)*dz(k))
	face(CUTTING_FACE_N)%centroid=(/cx,y(j+1),cy/)
  ELSE
    ify(I,J+1,K)=0
  ENDIF
END SUBROUTINE

SUBROUTINE check_back_face(I,J,K)

  USE global_parameters
  USE flow_parameters
  USE grid_arrays
  USE flow_arrays
  USE boundary_arrays
  USE GCM_arrays
  USE unstructured_surface_arrays

  implicit none

  integer :: i, j, k
  INTEGER(4) :: FIND_INTERSECTION
  REAL(KIND=CGREAL) :: CALC_CUTTING_AREA, TMP, CALC_CUTTING_AREA_NORMALIZED
  REAL(KIND=CGREAL) :: cx, cy

!  IF (iex(i,  j,k)<0) iex(i,  j,k) = FIND_INTERSECTION(I,  J,K,1)
!  IF (iex(i,j+1,k)<0) iex(i,j+1,k) = FIND_INTERSECTION(I,J+1,K,2)
!  IF (iey(i,  j,k)<0) iey(i,  j,k) = FIND_INTERSECTION(I,  J,K,3)
!  IF (iey(i+1,j,K)<0) iey(i+1,j,K) = FIND_INTERSECTION(I+1,J,K,4)
  INTEGER :: TM(4)

  tm(1)=iex(i  ,j  ,k) ! 1
  tm(2)=iey(i+1,j  ,k) ! 4
  tm(3)=iex(i  ,j+1,k) ! 2
  tm(4)=iey(i  ,  j,k) ! 3
!  TMP=CALC_CUTTING_AREA(cx, cy, tm, &
!    xc(i), yc(j), x(i), x(i+1), y(j), y(j+1), dx(i), dy(j), &
!    iblank_in(i+1,j,k), iblank_in(i+1,j+1,k), iblank_in(i,j+1,k), iblank_in(i,j,k))
  TMP=CALC_CUTTING_AREA_NORMALIZED(cx, cy, tm, &
    xc(i), yc(j), x(i), x(i+1), y(j), y(j+1), dxinv(i), dyinv(j), &
    iblank_in(i+1,j,k), iblank_in(i+1,j+1,k), iblank_in(i,j+1,k), iblank_in(i,j,k))

!  TMP=CALC_CUTTING_AREA_Z(I,J,K,centroid)
  IF (TMP>ZERO) THEN
	CUTTING_FACE_N=CUTTING_FACE_N+1
	ifZ(I,J,K)=CUTTING_FACE_N
!	face(CUTTING_FACE_N)%area=TMP
	face(CUTTING_FACE_N)%a   =TMP !/(dx(i)*dy(j))
	face(CUTTING_FACE_N)%centroid=(/cx,cy,z(k)/)
  ELSE
	ifZ(I,J,K)=0
  ENDIF
END SUBROUTINE


SUBROUTINE check_left_face(I,J,K)

  USE global_parameters
  USE flow_parameters
  USE grid_arrays
  USE flow_arrays
  USE boundary_arrays
  USE GCM_arrays
  USE unstructured_surface_arrays

  implicit none

  integer :: i, j, k
  INTEGER(4) :: FIND_INTERSECTION
  REAL(KIND=CGREAL) :: CALC_CUTTING_AREA, TMP, CALC_CUTTING_AREA_NORMALIZED
  REAL(KIND=CGREAL) :: cx, cy

!  IF (iez(i,  j,  k)<0) iez(i,  j,  k) = FIND_INTERSECTION(I,  J,K  ,9)
!  IF (iez(i,j+1,  K)<0) iez(i,j+1,  K) = FIND_INTERSECTION(I,J+1,K  ,10)
!  IF (iey(i,  j,  k)<0) iey(i,  j,  k) = FIND_INTERSECTION(I,  J,K  ,3)
!  IF (iey(i,  j,k+1)<0) iey(i,  j,k+1) = FIND_INTERSECTION(I,  J,K+1,7)
  INTEGER :: TM(4)

  tm(1)=iez(i,j  ,k)   ! 9
  tm(2)=iey(i,j  ,k+1) ! 7
  tm(3)=iez(i,j+1,k)   ! 10
  tm(4)=iey(i,j  ,k)   ! 3
!  TMP=CALC_CUTTING_AREA(cx, cy, tm, &
!    zc(k), yc(j), z(k), z(k+1), y(j), y(j+1), dz(k), dy(j), &
!    iblank_in(i,j,k+1), iblank_in(i,j+1,k+1), iblank_in(i,j+1,k), iblank_in(i,j,k))
  TMP=CALC_CUTTING_AREA_NORMALIZED(cx, cy, tm, &
    zc(k), yc(j), z(k), z(k+1), y(j), y(j+1), dzinv(k), dyinv(j), &
    iblank_in(i,j,k+1), iblank_in(i,j+1,k+1), iblank_in(i,j+1,k), iblank_in(i,j,k))
  
!  TMP=CALC_CUTTING_AREA_X(I,J,K,centroid)
  IF (TMP>ZERO) THEN
	CUTTING_FACE_N=CUTTING_FACE_N+1
	ifx(I,J,K)=CUTTING_FACE_N
!	face(CUTTING_FACE_N)%area=TMP
	face(CUTTING_FACE_N)%a   =TMP !/(dy(j)*dz(k))
	face(CUTTING_FACE_N)%centroid=(/x(i),cy,cx/)
  ELSE
    ifx(I,J,K)=0
  ENDIF
END SUBROUTINE

SUBROUTINE check_low_face(I,J,K)

  USE global_parameters
  USE flow_parameters
  USE grid_arrays
  USE flow_arrays
  USE boundary_arrays
  USE GCM_arrays
  USE unstructured_surface_arrays

  implicit none

  integer :: i, j, k
  INTEGER(4) :: FIND_INTERSECTION
  REAL(KIND=CGREAL) :: CALC_CUTTING_AREA, TMP, CALC_CUTTING_AREA_NORMALIZED
  REAL(KIND=CGREAL) :: cx, cy

!  IF (iex(i  ,j,  k)<0) iex(i  ,j,  k) = FIND_INTERSECTION(I  ,J,  K,1)
!  IF (iex(i  ,j,K+1)<0) iex(i  ,j,K+1) = FIND_INTERSECTION(I  ,J,K+1,5)
!  IF (iez(i  ,j,  k)<0) iez(i  ,j,  k) = FIND_INTERSECTION(I  ,J,  K,9)
!  IF (iez(i+1,j,  k)<0) iez(i+1,j,  k) = FIND_INTERSECTION(I+1,J,  K,11)
  INTEGER :: TM(4)

  tm(1)=iex(i  ,j,k)   ! 1
  tm(2)=iez(i+1,j,k)   ! 11
  tm(3)=iex(i  ,j,k+1) ! 5
  tm(4)=iez(i  ,j,k)   ! 9
!  TMP=CALC_CUTTING_AREA(cx, cy, tm, &
!    xc(i), zc(k), x(i), x(i+1), z(k), z(k+1), dx(i), dz(k), &
!    iblank_in(i+1,j,k), iblank_in(i+1,j,k+1), iblank_in(i,j,k+1), iblank_in(i,j,k))
  TMP=CALC_CUTTING_AREA_NORMALIZED(cx, cy, tm, &
    xc(i), zc(k), x(i), x(i+1), z(k), z(k+1), dxinv(i), dzinv(k), &
    iblank_in(i+1,j,k), iblank_in(i+1,j,k+1), iblank_in(i,j,k+1), iblank_in(i,j,k))

!  TMP=CALC_CUTTING_AREA_Y(I,J,K,centroid)
  IF (TMP>ZERO) THEN
	CUTTING_FACE_N=CUTTING_FACE_N+1
	ifY(I,J,K)=CUTTING_FACE_N
!	face(CUTTING_FACE_N)%area=TMP
	face(CUTTING_FACE_N)%a   =TMP !/(dx(i)*dz(k))
	face(CUTTING_FACE_N)%centroid=(/cx,y(j),cy/)
  ELSE
    ify(I,J,K)=0
  ENDIF
END SUBROUTINE

!FUNCTION CALC_CUTTING_AREA_X(I,J,K, centroid)
!  USE global_parameters
!!  USE unstructured_surface_arrays
!  USE flow_arrays
!  USE grid_arrays
!  USE boundary_arrays
!
!  implicit none
!
!  integer :: i, j, k, ii, mn
!  TYPE (VECTOR) :: centroid
!  REAL(KIND=CGREAL) :: z11, y4, CALC_CUTTING_AREA_X, a, a1, h1, h2
!
!  INTEGER :: TMP(4) ! last two elements are to avoid bound check
!
!  centroid%x(1)=x(i)
!  tmp=0
!  tmp(1)=iez(i,j,k)
!  tmp(2)=iey(i,j,k)
!  tmp(3)=iez(i,j+1,k)
!  tmp(4)=iey(i,j,k+1)
!  mn=0
!  do ii=1, 4
!  if (tmp(ii)<0) mn=mn+1
!  end do
!  if (mn/=2 .and. mn/=4) then
!    write(*,*)  mn, 'intercepted points in face x.', i,j,k
!    stop
!  endif
!
!  a1=dy(j)*dz(k)
!  if (tmp(1)>0) then
!    z11=edge_cut(tmp(1))
!    if (tmp(2)>0) then
!      h1=z11-z(k)
!      h2=edge_cut(tmp(2))-y(j)
!      a=half*(h1)*(h2)
!      centroid%x(2)=oned/3*(h2)+y(j)
!      centroid%x(3)=oned/3*(h1)+z(k)
!      if (iblank_in(i,j,k)>0) a=a1-a
!    else if (tmp(3)>0) then
!      h1=z11-z(k)
!      h2=edge_cut(tmp(3))-z(k)
!      a=half*(h1+h2)*dy(j)
!      centroid%x(2)=dy(j)*(h1+2*h2)/(3*(h1+h2))+y(j)
!      centroid%x(3)=(h1*h1/(h1+h2)+h2)/3+z(k)
!      if (iblank_in(i,j,k)>0) a=a1-a
!    else
!      h1=z(k+1)-z11
!      h2=edge_cut(tmp(4))-y(j)
!      a=half*(h1)*(h2)
!      centroid%x(2)=oned/3*(h2)+y(j)
!      centroid%x(3)=z(k+1)-oned/3*(h1)
!      if (iblank_in(i,j,k+1)>0) a=a1-a
!    endif
!  else if (tmp(2)>0) then
!    y4=edge_cut(tmp(2))
!    h1=y(j+1)-y4
!    if (tmp(3)>0) then
!      h2=edge_cut(tmp(3))-z(k)
!      a=half*(h1)*(h2)
!      centroid%x(2)=y(j+1)-oned/3*(h1)
!      centroid%x(3)=oned/3*(h2)+z(k)
!      if (iblank_in(i,j+1,k)>0) a=a1-a
!    else
!      h2=y(j+1)-edge_cut(tmp(4))
!      a=half*(h1+h2)*dz(k)
!      centroid%x(2)=y(j+1)-(h1*h1/(h1+h2)+h2)/3
!      centroid%x(3)=dz(k)*(h1+2*h2)/(3*(h1+h2))+z(k)
!      if (iblank_in(i,j+1,k)>0) a=a1-a
!    endif
!  else if (tmp(3)>0) then
!      h1=z(k+1)-edge_cut(tmp(3))
!      h2=y(j+1)-edge_cut(tmp(4))
!      a=half*(h1)*(h2)
!      centroid%x(2)=y(j+1)-oned/3*(h2)
!      centroid%x(3)=z(k+1)-oned/3*(h1)
!      if (iblank_in(i,j+1,k+1)>0) a=a1-a
!  else
!    a=zero
!  endif
!  CALC_CUTTING_AREA_X=a
!END FUNCTION
!
!FUNCTION CALC_CUTTING_AREA_Y(I,J,K, centroid)
!  USE global_parameters
!!  USE unstructured_surface_arrays
!  USE flow_arrays
!  USE grid_arrays
!  USE boundary_arrays
!
!  implicit none
!
!  integer :: i, j, k, ii, mn
!  TYPE (VECTOR) :: centroid
!  REAL(KIND=CGREAL) :: x6, z12, CALC_CUTTING_AREA_Y, a, a1, h1, h2
!  INTEGER :: TMP(4) ! last two elements are to avoid bound check
!
!  centroid%x(2)=y(j)
!  tmp=0
!  tmp(1)=iex(i,j,k+1)
!  tmp(2)=iez(i+1,j,k)
!  tmp(3)=iex(i,j,k)
!  tmp(4)=iez(i,j,k)
!  mn=0
!  do ii=1, 4
!  if (tmp(ii)<0) mn=mn+1
!  end do
!  if (mn/=2 .and. mn/=4) then
!    write(*,*)  mn, 'intercepted points in face y.', i,j,k
!    stop
!  endif
!
!  a1=dx(i)*dz(k)
!  if (tmp(1)>0) then
!    x6=edge_cut(tmp(1))
!    if (tmp(2)>0) then
!      h1=x(i+1)-x6
!      h2=z(k+1)-edge_cut(tmp(2))
!      a=half*(h1)*(h2)
!      centroid%x(1)=x(i+1)-oned/3*(h2)
!      centroid%x(3)=z(k+1)-oned/3*(h1)
!      if (iblank_in(i+1,j,k+1)>0) a=a1-a
!    else if (tmp(3)>0) then
!      h1=x(i+1)-x6
!      h2=x(i+1)-edge_cut(tmp(3))
!      a=half*(h1+h2)*dz(k)
!      centroid%x(1)=dy(j)*(h1+2*h2)/(3*(h1+h2))+y(j)
!      centroid%x(3)=(h1*h1/(h1+h2)+h2)/3+z(k)
!      if (iblank_in(i+1,j,k+1)>0) a=a1-a
!    else
!      a=half*(edge_cut(tmp(1))-x(i))*(z(k+1)-edge_cut(tmp(4)))
!      if (iblank_in(i,j,k+1)>0) a=a1-a
!    endif
!  else if (tmp(2)>0) then
!    if (tmp(3)>0) then
!      a=half*(edge_cut(tmp(2))-z(k))*(x(i+1)-edge_cut(tmp(3)))
!      if (iblank_in(i+1,j,k)>0) a=a1-a
!    else
!      a=half*(edge_cut(tmp(2))-z(k)+edge_cut(tmp(4))-z(k))*dx(i)
!      if (iblank_in(i+1,j,k)>0) a=a1-a
!    endif
!  else if (tmp(3)>0) then
!      a=half*(edge_cut(tmp(3))-x(i))*(edge_cut(tmp(4))-z(k))
!      if (iblank_in(i,j,k)>0) a=a1-a
!  else
!    a=zero
!  endif
!  CALC_CUTTING_AREA_Y=a
!END FUNCTION
!
!FUNCTION CALC_CUTTING_AREA_Z(I,J,K, centroid)
!  USE global_parameters
!!  USE unstructured_surface_arrays
!  USE flow_arrays
!  USE grid_arrays
!  USE boundary_arrays
!
!  implicit none
!
!  integer :: i, j, k, mn, ii
!  TYPE (VECTOR) :: centroid
!  REAL(KIND=CGREAL) :: x5, y8, x6, y7, CALC_CUTTING_AREA_Z,a, a1, h1, h2
!  INTEGER :: TMP(4) ! last two elements are to avoid bound check
!
!  tmp=0
!  tmp(1)=iex(i,j,k)
!  tmp(2)=iey(i+1,j,k)
!  tmp(3)=iex(i,j+1,k)
!  tmp(4)=iey(i,j,k)
!  mn=0
!  do ii=1, 4
!  if (tmp(ii)<0) mn=mn+1
!  end do
!  if (mn/=2 .and. mn/=4) then
!    write(*,*)  mn, 'intercepted points in face z.', i,j,k
!    stop
!  endif
!  a1=dx(i)*dy(j)
!  if (tmp(1)>0) then
!    x5=edge_cut(tmp(1))
!	if (tmp(2)>0) then
!      a=half*(x(i+1)-x5)*(edge_cut(tmp(2))-y(j))
!      if (iblank_in(i+1,j,k)>0) a=a1-a
!	else if (tmp(3)>0) then
!	  a=half*(x(i+1)-x5+x(i+1)-edge_cut(tmp(3)))*dy(j)
!      if (iblank_in(i+1,j,k)>0) a=a1-a
!	else
!	  a=half*(x5-x(i))*(edge_cut(tmp(4))-y(j))
!	  if (iblank_in(i,j,k)>0) a=a1-a
!	endif
!  else if (tmp(2)>0) then
!    y8=edge_cut(tmp(2))
!    if (tmp(3)>0) then
!      a=half*(y(j+1)-y8)*(x(i+1)-edge_cut(tmp(3)))
!      if (iblank_in(i+1,j+1,k)>0) a=a1-a
!    else
!      a=half*(y(j+1)-y8+y(j+1)-edge_cut(tmp(4)))*dx(i)
!      if (iblank_in(i+1,j+1,k)>0) a=a1-a
!    endif
!  else if (tmp(3)>0) then
!    x6=edge_cut(tmp(3))
!    y7=edge_cut(tmp(4))
!    a=half*(x6-x(i))*(y(j+1)-y7)
!    if (iblank_in(i,j+1,k)>0) a=a1-a
!  else
!    a=zero
!  endif
!  CALC_CUTTING_AREA_Z=a
!END FUNCTION

SUBROUTINE invert_centroid(xc,yc,cx,cy,a,a1)
  USE global_parameters

  REAL(KIND=CGREAL) :: xc, yc
  REAL(KIND=CGREAL) :: cx, cy
  REAL(KIND=CGREAL) :: a, a1

  cx=xc*a1-cx*a
  cy=yc*a1-cy*a
  a=a1-a
  cx=cx/a
  cy=cy/a

END SUBROUTINE

FUNCTION CALC_CUTTING_AREA(cx, cy, tmp, xc, yc, x0, x1, y0, y1, dx, dy, in1, in2, in3, in4)
  USE global_parameters
  USE flow_arrays
  USE boundary_arrays

  implicit none

  integer :: mn, ii
  REAL(KIND=CGREAL) :: cx, cy
  INTEGER :: TMP(4)
  REAL(KIND=CGREAL) :: xc, yc, x0, x1, y0, y1, dx, dy
  INTEGER(1) :: in1, in2, in3, in4
  REAL(KIND=CGREAL) :: x5, y8, x6, y7, CALC_CUTTING_AREA,a, a1, h1, h2

  mn=0
  do ii=1, 4
  if (tmp(ii)>0) mn=mn+1
  end do
  if (mn==1) then
!        CALL write_dump_debug_i('ibln',0,iblank)
    write(*,*)  'Only one intercepted point in face.', x0, x1, y0, y1
    stop
  else if (mn==0) then
    CALC_CUTTING_AREA=zero
    return
  endif

  a=zero
  a1=dx*dy
  if (tmp(1)>0) then
    x5=edge_cut(tmp(1))
    IF (X5/=x1) THEN
	if (tmp(2)>0) then
      h1=x1-x5
	  h2=edge_cut(tmp(2))-y0
      a=half*(h1)*(h2)
      cx=x1 - oned/3*(h1)
      cy=y0 + oned/3*(h2)
      if (in1>0) call invert_centroid(xc,yc,cx,cy,a,a1)
	else if (tmp(3)>0) then
      h1=x1-x5
	  h2=x1-edge_cut(tmp(3))
	  a=half*(h1+h2)*dy
      cx=x1 - (h1*h1/(h1+h2)+h2)/3
      cy=y0 + dy*(h1+2*h2)/(3*(h1+h2))
      if (in1>0) call invert_centroid(xc,yc,cx,cy,a,a1)
	else
	  h1=x5-x0
	  h2=edge_cut(tmp(4))-y0
	  a=half*(h1)*(h2)
      cx=x0 + oned/3*(h1)
      cy=y0 + oned/3*(h2)
      if (in4>0) call invert_centroid(xc,yc,cx,cy,a,a1)
	endif
	ENDIF
  else if (tmp(2)>0) then
    y8=edge_cut(tmp(2))
    IF (Y8/=Y1) THEN
    h1=y1-y8
    if (tmp(3)>0) then
      h2=x1-edge_cut(tmp(3))
      a=half*(h1)*(h2)
      cx=x1-oned/3*(h2)
      cy=y1-oned/3*(h1)
      if (in2>0) call invert_centroid(xc,yc,cx,cy,a,a1)
    else
      h2=y1-edge_cut(tmp(4))
      a=half*(h1+h2)*dx
      cx=x1-dx*(h1+2*h2)/(3*(h1+h2))
      cy=y1-(h1*h1/(h1+h2)+h2)/3
      if (in2>0) call invert_centroid(xc,yc,cx,cy,a,a1)
    endif
    ENDIF
  else if (tmp(3)>0) then
    x6=edge_cut(tmp(3))
    IF (X6/=X0) THEN
    y7=edge_cut(tmp(4))
    h1=x6-x0
    h2=y1-y7
    a=half*(h1)*(h2)
    cx=x0 + oned/3*(h1)
    cy=y1 - oned/3*(h2)
    if (in3>0) call invert_centroid(xc,yc,cx,cy,a,a1)
    ENDIF
  endif
  CALC_CUTTING_AREA=a
END FUNCTION

FUNCTION CALC_CUTTING_AREA_NORMALIZED (cx, cy, tmp, xc, yc, x0, x1, y0, y1, dxinv, dyinv, in1, in2, in3, in4)

  USE global_parameters
  USE flow_arrays
  USE boundary_arrays

  implicit none

  integer :: mn, ii
  REAL(KIND=CGREAL) :: cx, cy
  INTEGER :: TMP(4)
  REAL(KIND=CGREAL) :: xc, yc, x0, x1, y0, y1, dxinv, dyinv
  INTEGER(1) :: in1, in2, in3, in4
  REAL(KIND=CGREAL) :: x5, y8, x6, y7, CALC_CUTTING_AREA_NORMALIZED,a, a1, h1, h2

  mn=0
  do ii=1, 4
  if (tmp(ii)>0) mn=mn+1
  end do
  if (mn==1) then
!        CALL write_dump_debug_i('ibln',0,iblank)
    write(*,*)  'Only one intercepted point in face.', x0, x1, y0, y1
    stop
  else if (mn==0) then
    CALC_CUTTING_AREA_NORMALIZED=zero
    return
  endif
  
  a=zero
  a1=oned
  if (tmp(1)>0) then
    x5=edge_cut(tmp(1))
    IF (X5/=x1) THEN
	if (tmp(2)>0) then
      h1=x1-x5
	  h2=edge_cut(tmp(2))-y0
      a=half*dxinv*(h1)*dyinv*(h2)
      cx=x1 - oned/3*(h1)
      cy=y0 + oned/3*(h2)
      if (in1>0) call invert_centroid(xc,yc,cx,cy,a,a1)
	else if (tmp(3)>0) then
      h1=x1-x5
	  h2=x1-edge_cut(tmp(3))
	  a=half*dxinv*(h1+h2)
      cx=x1 - (h1*h1/(h1+h2)+h2)/3
      cy=y0 + (h1+2*h2)/(3*(h1+h2))/dyinv
      if (in1>0) call invert_centroid(xc,yc,cx,cy,a,a1)
	else
	  h1=x5-x0
	  h2=edge_cut(tmp(4))-y0
	  a=half*dxinv*(h1)*dyinv*(h2)
      cx=x0 + oned/3*(h1)
      cy=y0 + oned/3*(h2)
      if (in4>0) call invert_centroid(xc,yc,cx,cy,a,a1)
	endif
	ENDIF
  else if (tmp(2)>0) then
    y8=edge_cut(tmp(2))
    IF (Y8/=Y1) THEN
    h1=y1-y8
    if (tmp(3)>0) then
      h2=x1-edge_cut(tmp(3))
      a=half*dyinv*(h1)*dxinv*(h2)
      cx=x1-oned/3*(h2)
      cy=y1-oned/3*(h1)
      if (in2>0) call invert_centroid(xc,yc,cx,cy,a,a1)
    else
      h2=y1-edge_cut(tmp(4))
      a=half*dyinv*(h1+h2)
      cx=x1-(h1+2*h2)/(3*(h1+h2))/dxinv
      cy=y1-(h1*h1/(h1+h2)+h2)/3
      if (in2>0) call invert_centroid(xc,yc,cx,cy,a,a1)
    endif
    ENDIF
  else if (tmp(3)>0) then
    x6=edge_cut(tmp(3))
    IF (X6/=X0) THEN
    y7=edge_cut(tmp(4))
    h1=x6-x0
    h2=y1-y7
    a=half*dxinv*(h1)*dyinv*(h2)
    cx=x0 + oned/3*(h1)
    cy=y1 - oned/3*(h2)
    if (in3>0) call invert_centroid(xc,yc,cx,cy,a,a1)
    ENDIF
  endif
  if (a==oned) &
    a=zero
  CALC_CUTTING_AREA_NORMALIZED=a
END FUNCTION

!FUNCTION FIND_INTERSECTION(I,J,K,DIR)
!  USE global_parameters
!  USE unstructured_surface_arrays
!  USE flow_parameters
!  USE boundary_arrays
!  USE grid_arrays
!  USE flow_arrays
!  
!  implicit none
!
!  integer :: i, j, k,DIR
!  LOGICAL :: check_POINT_IN_TRIANGLE
!  TYPE (VECTOR) :: NORM, V0, V1
!  integer :: iBody,m, mVert1
!  INTEGER(4) :: FIND_INTERSECTION
!
!  DO iBody=1, nBody_Solid
!    DO m=1,totNumTriElem(iBody)
!      NORM%X(1)=triElemNormX(iBody,m)
!      NORM%X(2)=triElemNormY(iBody,m)
!      NORM%X(3)=triElemNormZ(iBody,m)
!
!      mVert1   = triElemNeig(iBody,1,m)
!
!      V0%X(1)=xBodyMarker(iBody,mVert1)
!      V0%X(2)=yBodyMarker(iBody,mVert1)
!      V0%X(3)=zBodyMarker(iBody,mVert1)
!
!      SELECT CASE (DIR)
!        CASE (1,2,5,6)
!          IF (ABS(NORM%X(1))<1E-8) CYCLE ! parallel to the line
!          V1%X(2)=y(j)
!          V1%X(3)=z(k)
!          V1%X(1)=(norm%x(2)*(V0%X(2)-V1%X(2))+norm%x(3)*(V0%X(3)-V1%X(3)))/norm%x(1)+V0%X(1)
!          if (V1%X(1)>=x(i) .and. V1%X(1)<=x(i+1)) then
!            IF (check_POINT_IN_TRIANGLE(iBody,m,V1)) THEN
!              CUTTING_EDGE_N=CUTTING_EDGE_N+1
!              edge_CUT(CUTTING_EDGE_N)=V1%X(1)
!              if (norm%x(1)<zero)  iblank_in(i,j,k)=1 !positive normal vector points into body
!              FIND_INTERSECTION=CUTTING_EDGE_N
!              return
!            ENDIF
!          ENDIF
!
!        CASE (9,10,11,12)
!          IF (ABS(NORM%X(3))<1E-8) CYCLE
!          V1%X(1)=x(i)
!          V1%X(2)=y(j)
!          V1%X(3)=(norm%x(1)*(V0%X(1)-V1%X(1))+norm%x(2)*(V0%X(2)-V1%X(2)))/norm%x(3)+V0%X(3)
!          if (V1%X(3)>=Z(K) .and. V1%X(3)<=z(k+1)) then
!            IF (check_POINT_IN_TRIANGLE(iBody,m,V1)) THEN
!              CUTTING_EDGE_N=CUTTING_EDGE_N+1
!              edge_CUT(CUTTING_EDGE_N)=V1%X(3)
!              if (norm%x(3)<zero)  iblank_in(i,j,k)=1
!              FIND_INTERSECTION=CUTTING_EDGE_N
!              return
!            ENDIF
!          ENDIF
!
!        CASE (3,4,7,8)
!          IF (ABS(NORM%X(2))<1E-8) CYCLE
!          V1%X(1)=x(i)
!          V1%X(3)=z(k)
!          V1%X(2)=(norm%x(1)*(V0%X(1)-V1%X(1))+norm%x(3)*(V0%X(3)-V1%X(3)))/norm%x(2)+V0%X(2)
!          if (V1%X(2)>=y(j) .and. V1%X(2)<=y(j+1)) then
!            IF (check_POINT_IN_TRIANGLE(iBody,m,V1)) THEN
!              CUTTING_EDGE_N=CUTTING_EDGE_N+1
!              edge_CUT(CUTTING_EDGE_N)=V1%X(2)
!              if (NORM%X(2)<zero)  iblank_in(i,j,k)=1
!              FIND_INTERSECTION=CUTTING_EDGE_N
!              return
!            ENDIF
!          ENDIF
!
!      END SELECT
!    END DO ! m
!  END DO ! iBody
!  FIND_INTERSECTION=0
!END FUNCTION



FUNCTION check_POINT_IN_TRIANGLE2(A,B,C,P)
  USE global_parameters
  USE unstructured_surface_arrays
  USE boundary_arrays

  implicit none

  LOGICAL :: check_POINT_IN_TRIANGLE2
  REAL(KIND=CGREAL), DIMENSION(3) :: P, A, B, C, V0, V1, V2
  REAL(KIND=CGREAL) :: dot00, dot01, dot02, dot11, dot12
  REAL(KIND=CGREAL) :: invDenom, U, V
  
! Compute vectors
  v0 = C - A
  v1 = B - A
  v2 = P - A

! Compute dot products
  dot00 = DOT_PRODUCT(v0, v0)
  dot01 = DOT_PRODUCT(v0, v1)
  dot02 = DOT_PRODUCT(v0, v2)
  dot11 = DOT_PRODUCT(v1, v1)
  dot12 = DOT_PRODUCT(v1, v2)

! Compute barycentric coordinates
  invDenom = 1.D0 / (dot00 * dot11 - dot01 * dot01)
  u = (dot11 * dot02 - dot01 * dot12) * invDenom
  v = (dot00 * dot12 - dot01 * dot02) * invDenom

! Check if point is in triangle
  IF (u >= -1E-10 .AND. v >= -1E-10 .AND. u + v <oned + 1E-10) THEN
    check_POINT_IN_TRIANGLE2=.TRUE.
  ELSE
    check_POINT_IN_TRIANGLE2=.FALSE.
  ENDIF

!  TYPE (VECTOR) :: cross01, cross12, cross20
!  REAL(KIND=CGREAL) :: dot12, dot01, dot02
!!  REAL(KIND=CGREAL) :: invDenom, U, V
!! Compute vectors
!  v0%X = A%X - P%X
!  v1%X = B%X - P%X
!  v2%X = C%X - P%X
!
!! Compute dot products
!  cross01 = cross(v0, v1)
!  cross12 = cross(v1, v2)
!  cross20 = cross(v2, v0)
!
!  dot01 = dot(cross01,cross12)
!  dot02 = dot(cross01,cross20)
!!  dot12 = dot(cross12,cross20)
!  if ((dot01>=zero .and. dot02>=zero) .or. &
!     (dot01<=zero .and. dot02<=zero)) then 
!    check_POINT_IN_TRIANGLE2=.TRUE.
!  ELSE
!    check_POINT_IN_TRIANGLE2=.FALSE.
!  ENDIF
END FUNCTION



SUBROUTINE FIND_INTERSECTION_X(Norm, A, B, C, J,K)
  USE global_parameters
  USE unstructured_surface_arrays
  USE flow_parameters
  USE boundary_arrays
  USE grid_arrays
  USE flow_arrays
  
  implicit none

  integer :: i, j, k
  LOGICAL :: check_POINT_IN_TRIANGLE2
  REAL(KIND=CGREAL), DIMENSION(3) :: NORM, A, B, C, V1

  IF (ABS(NORM(1))<1E-10) return ! parallel to the line
  V1(2)=y(j)
  V1(3)=z(k)
  V1(1)=(norm(2)*(A(2)-V1(2))+norm(3)*(A(3)-V1(3)))/norm(1)+A(1)
  IF (V1(1)>=x(1) .and. V1(1)<=x(nx)) THEN
    IF (check_POINT_IN_TRIANGLE2(A,B,C,V1)) THEN
      DO i=1, nx
        IF (x(i)<=V1(1) .AND. V1(1)<x(i+1)) EXIT
      ENDDO
      if (iex(i,j,k) >0) then
        if (edge_CUT(iex(i,j,k))/=V1(1)) then
!          write(*,*) 'Edge in X dir was cut again!', i,j,k
!          CALL write_dump()
!          CALL write_dump_debug_i('ibln',0,iblank)
!          WRITE(*,*) y(j), z(k), dy(j), dz(k)
!          WRITE(*,*) A
!          WRITE(*,*) B
!          WRITE(*,*) C
!          WRITE(*,*) Norm
!          WRITE(*,*) V1
!          WRITE(*,*) edge_CUT(iex(i,j,k))
          if (iblank(i,j,k)==1) then
            iblank_in(i,j,k)=1
            iblank_in(i+1,j,k)=1
          else
            iblank_in(i,j,k)=0
            iblank_in(i+1,j,k)=0
          endif
!          stop
          iex(i,j,k)=-1
        endif
      else
        CUTTING_EDGE_N=CUTTING_EDGE_N+1
        edge_CUT(CUTTING_EDGE_N)=V1(1)
        iex(i,j,k) = CUTTING_EDGE_N
        IF (x(i)==V1(1)) iex(i-1,j,k) = CUTTING_EDGE_N
        if (norm(1)<zero) then
          iblank_in(i,j,k)=1 !positive normal vector points into body
        else
          iblank_in(i+1,j,k)=1 !positive normal vector points into body
        endif
      endif
    ENDIF
  ENDIF

END SUBROUTINE

SUBROUTINE FIND_INTERSECTION_Y(Norm, A, B, C, I,K)
  USE global_parameters
  USE unstructured_surface_arrays
  USE flow_parameters
  USE boundary_arrays
  USE grid_arrays
  USE flow_arrays
  
  implicit none

  integer :: i, j, k
  LOGICAL :: check_POINT_IN_TRIANGLE2
  REAL(KIND=CGREAL), DIMENSION(3) :: NORM, A, B, C, V1

  IF (ABS(NORM(2))<1E-10) return  ! parallel to the line
  V1(1)=x(i)
  V1(3)=z(k)
  V1(2)=(norm(1)*(A(1)-V1(1))+norm(3)*(A(3)-V1(3)))/norm(2)+A(2)
  IF (V1(2)>=y(1) .and. V1(2)<=y(ny)) THEN
    IF (check_POINT_IN_TRIANGLE2(A,B,C,V1)) THEN
      DO j=1, ny
        IF (y(j)<=V1(2) .AND. V1(2)<y(j+1)) EXIT
      ENDDO
      if (iey(i,j,k) >0) then
        if (edge_CUT(iey(i,j,k))/=V1(2)) then
!          write(*,*) 'Edge in Y dir was cut again!', i,j,k
!          CALL write_dump()
!          CALL write_dump_debug_i('ibln',0,iblank)
!          WRITE(*,*) x(i), z(k), dx(i), dz(k)
!          WRITE(*,*) A
!          WRITE(*,*) B
!          WRITE(*,*) C
!          WRITE(*,*) Norm
!          WRITE(*,*) V1
!          WRITE(*,*) edge_CUT(iey(i,j,k))
          if (iblank(i,j,k)==1) then
            iblank_in(i,j,k)=1
            iblank_in(i,j+1,k)=1
          else
            iblank_in(i,j,k)=0
            iblank_in(i,j+1,k)=0
          endif
!          stop
          iey(i,j,k)=-1
        endif
      else
        CUTTING_EDGE_N=CUTTING_EDGE_N+1
        edge_CUT(CUTTING_EDGE_N)=V1(2)
        iey(i,j,k) = CUTTING_EDGE_N
        if (y(j)==V1(2)) iey(i,j-1,k) = CUTTING_EDGE_N
        if (NORM(2)<zero) then
          iblank_in(i,j,k)=1
        else
          iblank_in(i,j+1,k)=1
        endif
      endif
    ENDIF
  ENDIF

END SUBROUTINE

SUBROUTINE FIND_INTERSECTION_Z(Norm, A, B, C, i, j)
  USE global_parameters
  USE unstructured_surface_arrays
  USE flow_parameters
  USE boundary_arrays
  USE grid_arrays
  USE flow_arrays
  
  implicit none

  integer :: i, j, k
  LOGICAL :: check_POINT_IN_TRIANGLE2
  REAL(KIND=CGREAL), DIMENSION(3) :: NORM, A, B, C, V1

  IF (ABS(NORM(3))<1E-10) return ! parallel to the line
  V1(1)=x(i)
  V1(2)=y(j)
  V1(3)=(norm(1)*(A(1)-V1(1))+norm(2)*(A(2)-V1(2)))/norm(3)+A(3)
  IF (V1(3)>=z(1) .and. V1(3)<=z(nz)) THEN
    IF (check_POINT_IN_TRIANGLE2(A,B,C,V1)) THEN
      DO k=1, nz
        IF (z(k)<=V1(3) .AND. V1(3)<z(k+1)) EXIT
      ENDDO
      if (iez(i,j,k) >0) then
        if (edge_CUT(iez(i,j,k))/=V1(3)) then
!          write(*,*) 'Edge in Z dir was cut again!', i,j,k
!          CALL write_dump()
!          CALL write_dump_debug_i('ibln',0,iblank)
!          WRITE(*,*) x(i), y(j), dx(i), dy(j)
!          WRITE(*,*) A
!          WRITE(*,*) B
!          WRITE(*,*) C
!          WRITE(*,*) Norm
!          WRITE(*,*) V1
!          WRITE(*,*) edge_CUT(iez(i,j,k))
          if (iblank(i,j,k)==1) then
            iblank_in(i,j,k)=1
            iblank_in(i,j,k+1)=1
          else
            iblank_in(i,j,k)=0
            iblank_in(i,j,k+1)=0
          endif
!          stop
          iez(i,j,k)=-1
        endif
      else
        CUTTING_EDGE_N=CUTTING_EDGE_N+1
        edge_CUT(CUTTING_EDGE_N)=V1(3)
        iez(i,j,k) = CUTTING_EDGE_N
        if (z(k)==V1(3)) iez(i,j,k-1) = CUTTING_EDGE_N
        if (norm(3)<zero) then
          iblank_in(i,j,k)=1 !positive normal vector points into body
        else
          iblank_in(i,j,k+1)=1 !positive normal vector points into body
        endif
      endif
    ENDIF
  ENDIF

END SUBROUTINE



SUBROUTINE FIND_EDGES_INTERCEPTING_TRIANGLE()
  USE global_parameters
  USE unstructured_surface_arrays
  USE flow_parameters
  USE boundary_arrays
  USE grid_arrays
  USE flow_arrays
  
  implicit none

  integer :: i, j, k,DIR
  REAL(KIND=CGREAL), DIMENSION(3) :: NORM, A, B, C
  integer :: iBody, m, mVert1, mVert2, mVert3
  REAL(KIND=CGREAL) :: xBoundMax,xBoundMin,yBoundMax,yBoundMin,&
                       zBoundMax,zBoundMin
  INTEGER          :: iCellMin,iCellMax,jCellMin,jCellMax,kCellMin,kCellMax

  DO iBody=1, nBody_Solid
    DO m=1,totNumTriElem(iBody)
	  NORM(1)=triElemNormX(iBody,m)
	  NORM(2)=triElemNormY(iBody,m)
	  NORM(3)=triElemNormZ(iBody,m)

	  mVert1  = triElemNeig(iBody,1,m)
	  mVert2  = triElemNeig(iBody,2,m)
	  mVert3  = triElemNeig(iBody,3,m)

	  A(1)  = xBodyMarker(iBody,mVert1)
	  A(2)  = yBodyMarker(iBody,mVert1)
	  A(3)  = zBodyMarker(iBody,mVert1)

	  B(1)  = xBodyMarker(iBody,mVert2)
	  B(2)  = yBodyMarker(iBody,mVert2)
	  B(3)  = zBodyMarker(iBody,mVert2)

	  C(1)  = xBodyMarker(iBody,mVert3)
	  C(2)  = yBodyMarker(iBody,mVert3)
	  C(3)  = zBodyMarker(iBody,mVert3)
	  
	  IF ( ( A(1) < x(1)  .AND. B(1) < x(1)  .AND. C(1) < x(1)  ) .OR.  &
		   ( A(1) > x(nx) .AND. B(1) > x(nx) .AND. C(1) > x(nx) ) .OR.  &
		   ( A(2) < y(1)  .AND. B(2) < y(1)  .AND. C(2) < y(1)  ) .OR.  &
		   ( A(2) > y(ny) .AND. B(2) > y(ny) .AND. C(2) > y(ny) ) .OR.  &
		   ( A(3) < z(1)  .AND. B(3) < z(1)  .AND. C(3) < z(1)  ) .OR.  &
		   ( A(3) > z(nz) .AND. B(3) > z(nz) .AND. C(3) > z(nz) )     ) CYCLE
	  
	  xBoundMin = MIN(A(1),B(1),C(1))
	  xBoundMax = MAX(A(1),B(1),C(1))

	  yBoundMin = MIN(A(2),B(2),C(2))
	  yBoundMax = MAX(A(2),B(2),C(2))

	  zBoundMin = MIN(A(3),B(3),C(3))
	  zBoundMax = MAX(A(3),B(3),C(3))

	  iCellMin = 0
	  iCellMax = -1
	  jCellMin = 0
	  jCellMax = -1
	  kCellMin = 0
	  kCellMax = -1

! ******************************************************************************
!        i-direction 
! ******************************************************************************

	  DO i = 1, nx
	    IF ( (x(i)-xBoundMin)*(x(i)-xBoundMax) <= zero ) THEN
		  iCellMin = i
		  EXIT
		ENDIF
      ENDDO
      DO i = nx, 1, -1
	    IF ( (x(i)-xBoundMin)*(x(i)-xBoundMax) <= zero ) THEN
		  iCellMax = i
		  EXIT
		ENDIF
	  ENDDO

! ******************************************************************************
!        j-direction 
! ******************************************************************************

	  DO j = 1, ny
	    IF ( (y(j)-yBoundMin)*(y(j)-yBoundMax) <= zero ) THEN
		  jCellMin = j
		  EXIT
		ENDIF
	  ENDDO
	  DO j = ny, 1, -1
	    IF ( (y(j)-yBoundMin)*(y(j)-yBoundMax) <= zero ) THEN
		  jCellMax = j
		  EXIT
		ENDIF
	  ENDDO

!	  IF (jCellMin == 0 ) CYCLE
!	  IF (jCellMax == 0 ) CYCLE


! ******************************************************************************
!        k-direction 
! ******************************************************************************

	  DO k = 1, nz
		IF ( (z(k)-zBoundMin)*(z(k)-zBoundMax) <= zero ) THEN
		  kCellMin = k
		  EXIT
		ENDIF
	  ENDDO
	  DO k = nz, 1, -1
		IF ( (z(k)-zBoundMin)*(z(k)-zBoundMax) <= zero ) THEN
		  kCellMax = k
		  EXIT
		ENDIF
	  ENDDO
	  
!	  IF (kCellMin == 0 ) CYCLE
!	  IF (kCellMax == 0 ) CYCLE
!
!	  IF (iCellMin == 0 .OR. iCellMax == 0 ) THEN
!	    IF (
!	  ENDIF

      DO K=kCellMin, kCellMax
      DO J=jCellMin, jCellMax
        CALL FIND_INTERSECTION_X(Norm, A, B, C, j, k)
      ENDDO
      ENDDO

      DO K=kCellMin, kCellMax
      DO I=iCellMin, iCellMax
        CALL FIND_INTERSECTION_Y(Norm, A, B, C, i, k)
      ENDDO
      ENDDO

      DO J=jCellMin, jCellMax
      DO I=iCellMin, iCellMax
        CALL FIND_INTERSECTION_Z(Norm, A, B, C, i, j)
      ENDDO
      ENDDO
  ENDDO !m
  ENDDO !iBody
END SUBROUTINE

!FUNCTION check_POINT_IN_TRIANGLE(iBody,m,P)
!  USE global_parameters
!  USE unstructured_surface_arrays
!  USE boundary_arrays
!
!  implicit none
!
!  integer :: iBody,m
!  LOGICAL :: check_POINT_IN_TRIANGLE
!  integer :: mVert1, mVert2, mVert3
!  TYPE (VECTOR) :: P, A, B, C, V0, V1, V2
!  REAL(KIND=CGREAL) :: dot00, dot01, dot02, dot11, dot12
!  REAL(KIND=CGREAL) :: invDenom, U, V
!  
!  mVert1  = triElemNeig(iBody,1,m)
!  mVert2  = triElemNeig(iBody,2,m)
!  mVert3  = triElemNeig(iBody,3,m)
!
!  A%X(1)  = xBodyMarker(iBody,mVert1)
!  A%X(2)  = yBodyMarker(iBody,mVert1)
!  A%X(3)  = zBodyMarker(iBody,mVert1)
!
!  B%X(1)  = xBodyMarker(iBody,mVert2)
!  B%X(2)  = yBodyMarker(iBody,mVert2)
!  B%X(3)  = zBodyMarker(iBody,mVert2)
!
!  C%X(1)  = xBodyMarker(iBody,mVert3)
!  C%X(2)  = yBodyMarker(iBody,mVert3)
!  C%X(3)  = zBodyMarker(iBody,mVert3)
!! Compute vectors
!  v0%X = C%X - A%X
!  v1%X = B%X - A%X
!  v2%X = P%X - A%X
!
!! Compute dot products
!  dot00 = dot(v0, v0)
!  dot01 = dot(v0, v1)
!  dot02 = dot(v0, v2)
!  dot11 = dot(v1, v1)
!  dot12 = dot(v1, v2)
!
!! Compute barycentric coordinates
!  invDenom = 1.D0 / (dot00 * dot11 - dot01 * dot01)
!  u = (dot11 * dot02 - dot01 * dot12) * invDenom
!  v = (dot00 * dot12 - dot01 * dot02) * invDenom
!
!! Check if point is in triangle
!  IF (u > 0 .AND. v > 0 .AND. u + v < 1.D0) THEN
!    check_POINT_IN_TRIANGLE=.TRUE.
!  ELSE
!    check_POINT_IN_TRIANGLE=.FALSE.
!  ENDIF
!END FUNCTION
SUBROUTINE add_point(n, A, x, y, z)
  USE global_parameters

  implicit none

  REAL(KIND=CGREAL)  :: A(3,7) 
  integer :: i, n
  REAL(KIND=CGREAL) :: x,y,z
  
  do i=1, n
    IF (A(1,i)==x .and. A(2,i)==Y .and. A(3,i)==Z) return
  end do

  n=n+1
  A(:,n)=(/x ,y, z/)

END SUBROUTINE
SUBROUTINE insert_point(n, A, dxinv, dyinv, dzinv)
  USE global_parameters

  implicit none

  REAL(KIND=CGREAL)  :: A(3,7) 
  REAL(KIND=CGREAL), DIMENSION(3) :: Norm, V1, V2, V3, tmp, dd
  integer :: i, j, k, n
  REAL(KIND=CGREAL) :: dot1, dxinv, dyinv, dzinv
  
!  V1%x=A(2)%x-A(1)%x
!  V2%x=A(3)%x-A(1)%x
!  Norm=Cross(V1,V2)

  dd=(/dxinv, dyinv, dzinv/)
  
  V1=(A(:,2)-A(:,1)) * dd
  V2=(A(:,3)-A(:,1)) * dd

  Norm=Cross(V1,V2) * dd

  do i=4, n
    do j=2, i-1
	    V1=(A(:,j)-A(:,i))*dd
	    V2=(A(:,j-1)-A(:,i))*dd

      V3=Cross(V1,V2)*dd
      dot1=DOT_PRODUCT(V3,NORM)
      IF (dot1>1E-10) THEN
        tmp=A(:,i)
        DO K=i, j+1, -1
          A(:,k)=A(:,k-1)
        enddo 
        A(:,j)=tmp
      ENDIF
    enddo
  end do
  
END SUBROUTINE

FUNCTION construct_slice_2(I,J,K,A)
  USE global_parameters
  USE flow_arrays
  USE grid_arrays
  USE boundary_arrays

  implicit none

  INTEGER :: construct_slice_2
  REAL(KIND=CGREAL)  :: A(3,7) 
  integer :: i, j, k, n, m, N1, N2, N3, N4, N5, N6
  LOGICAL CC(5), uc(6)
  
  n=0
  if (iex(i,j  ,k)>0)   call add_point(n,A,edge_cut(iex(i,j,k))    ,Y(J)  ,Z(K))
  if (iex(i,j+1,k)>0)   call add_point(n,A,edge_cut(iex(i,j+1,k))  ,Y(J+1),Z(K))
  if (iex(i,j  ,k+1)>0) call add_point(n,A,edge_cut(iex(i,j,k+1))  ,Y(J)  ,Z(K+1))
  if (iex(i,j+1,k+1)>0) call add_point(n,A,edge_cut(iex(i,j+1,k+1)),Y(J+1),Z(K+1))
  
  if (iey(i  ,j,k)>0)   call add_point(n,A,X(i)  ,edge_cut(iey(i,j,k))  ,Z(K))
  if (iey(i+1,j,k)>0)   call add_point(n,A,X(i+1),edge_cut(iey(i+1,j,k)),Z(K))
  if (iey(i  ,j,k+1)>0) call add_point(n,A,X(i)  ,edge_cut(iey(i,j,k+1)),Z(K+1))
  if (iey(i+1,j,k+1)>0) call add_point(n,A,X(i+1),edge_cut(iey(i+1,j,k+1)),Z(K+1))
  
  if (iez(i  ,j  ,k)>0) call add_point(n,A,X(i)  ,y(j)  ,edge_cut(iez(i,j,k)))
  if (iez(i  ,j+1,k)>0) call add_point(n,A,X(i)  ,y(j+1),edge_cut(iez(i,j+1,k)))
  if (iez(i+1,j  ,k)>0) call add_point(n,A,X(i+1),y(j)  ,edge_cut(iez(i+1,j,k)))
  if (iez(i+1,j+1,k)>0) call add_point(n,A,X(i+1),y(j+1),edge_cut(iez(i+1,j+1,k)))

  IF (n>5) then
    write(*,*) 'Find n>5in cell', i,j,k,n
    WRITE(*,*) A(:,1:n)
    stop
  endif
  
  if (N>3) call insert_point(n, A, dxinv(i), dyinv(j), dzinv(k))
  
  construct_slice_2=N
END FUNCTION


FUNCTION construct_slice(I,J,K,A)
  USE global_parameters
  USE flow_arrays
  USE grid_arrays
  USE boundary_arrays

  implicit none

  TYPE :: point_pair
    INTEGER :: index
    REAL(KIND=CGREAL), DIMENSION(3) :: pv
  END TYPE

  INTEGER :: construct_slice
  TYPE (point_pair) :: B(5), A(5)
  integer :: i, j, k, n, m, N1, N2, N3, N4, N5, N6
!  TYPE (VECTOR) :: Norm
  INTEGER :: face_line_index(4,6)
  DATA face_line_index /5, 8, 6, 7, 11, 8, 12, 4, 2, 12, 6, 10, 1, 4, 2, 3, 9, 7, 10, 3, 1, 11, 5, 9/
  LOGICAL CC(5), uc(6)
  
  n=0
  if (iex(i,j,k)>0) then 
    n=n+1
    B(n)%index=1
    B(n)%pv(1)=edge_cut(iex(i,j,k))
    B(n)%pv(2)=Y(J)
    B(n)%pv(3)=Z(K)
  endif
  if (iex(i,j+1,k)>0) then 
    n=n+1
    B(n)%index=2
    B(n)%pv(1)=edge_cut(iex(i,j+1,k))
    B(n)%pv(2)=Y(J+1)
    B(n)%pv(3)=Z(K)
  endif
  if (iex(i,j,k+1)>0) then 
    n=n+1
    B(n)%index=5
    B(n)%pv(1)=edge_cut(iex(i,j,k+1))
    B(n)%pv(2)=Y(J)
    B(n)%pv(3)=Z(K+1)
  endif
  if (iex(i,j+1,k+1)>0) then 
    n=n+1
    B(n)%index=6
    B(n)%pv(1)=edge_cut(iex(i,j+1,k+1))
    B(n)%pv(2)=Y(J+1)
    B(n)%pv(3)=Z(K+1)
  endif
  
  if (iey(i,j,k)>0) then 
    n=n+1
    B(n)%index=3
    B(n)%pv(1)=X(i)
    B(n)%pv(2)=edge_cut(iey(i,j,k))
    B(n)%pv(3)=Z(K)
  endif
  if (iey(i+1,j,k)>0) then 
    n=n+1
    B(n)%index=4
    B(n)%pv(1)=X(i+1)
    B(n)%pv(2)=edge_cut(iey(i+1,j,k))
    B(n)%pv(3)=Z(K)
  endif
  if (iey(i,j,k+1)>0) then 
    n=n+1
    B(n)%index=7
    B(n)%pv(1)=X(i)
    B(n)%pv(2)=edge_cut(iey(i,j,k+1))
    B(n)%pv(3)=Z(K+1)
  endif
  if (iey(i+1,j,k+1)>0) then 
    n=n+1
    B(n)%index=8
    B(n)%pv(1)=X(i+1)
    B(n)%pv(2)=edge_cut(iey(i+1,j,k+1))
    B(n)%pv(3)=Z(K+1)
  endif
  
  if (iez(i,j,k)>0) then 
    n=n+1
    B(n)%index=9
    B(n)%pv(1)=X(i)
    B(n)%pv(2)=y(j)
    B(n)%pv(3)=edge_cut(iez(i,j,k))
  endif
  if (iez(i,j+1,k)>0) then 
    n=n+1
    B(n)%index=10
    B(n)%pv(1)=X(i)
    B(n)%pv(2)=y(j+1)
    B(n)%pv(3)=edge_cut(iez(i,j+1,k))
  endif
  if (iez(i+1,j,k)>0) then 
    n=n+1
    B(n)%index=11
    B(n)%pv(1)=X(i+1)
    B(n)%pv(2)=y(j)
    B(n)%pv(3)=edge_cut(iez(i+1,j,k))
  endif
  if (iez(i+1,j+1,k)>0) then 
    n=n+1
    B(n)%index=12
    B(n)%pv(1)=X(i+1)
    B(n)%pv(2)=y(j+1)
    B(n)%pv(3)=edge_cut(iez(i+1,j+1,k))
  endif
  IF (n>5) then
    write(*,*) 'Find n>5 in cell', i,j,k,n
    stop
  endif
  
  CC(1:N)=.TRUE.
  do n1=1, 6
  do n2=1, 4
	if (face_line_index(n2,n1)==B(1)%index) goto 10
  enddo
  enddo
10  continue
  A(1)=B(1)
  CC(1)=.FALSE.
  
  do m=2, n  
    do n3=1, 4
	do n4=1, n
	  if (face_line_index(n3,n1)==B(n4)%index .and. cc(N4)) goto 20
	enddo
	enddo
20  continue
    A(m)=B(n4)
    CC(N4)=.FALSE.
    
	do n5=1, 6
    do n6=1, 4
	  if (face_line_index(n6,n5)==face_line_index(n3,n1) .AND. N5/=N1) goto 30
	enddo
	enddo
30  continue
    
    n1=n5
  enddo
  
  construct_slice=N
!! find uncut face, set its a to zero
!  uc(1:6)=.TRUE.
!  do n1=1, 6
!  do n2=1, 4
!	if (A(n1)%index==face_line_index(n2,n1)) uc(n1)=.false.
!  enddo
!  enddo
!  
!  do n1=1, 6
!    if (uc(n1)) then
!	  select case (n1)
!	  case (1)
!		CELL(CUTTING_CELL_N)%F_kp=-1
!	  case (2)
!		CELL(CUTTING_CELL_N)%F_ip=-1
!	  case (3)
!		CELL(CUTTING_CELL_N)%F_jp=-1
!	  case (4)
!		CELL(CUTTING_CELL_N)%F_km=-1
!	  case (5)
!		CELL(CUTTING_CELL_N)%F_im=-1
!	  case (6)
!		CELL(CUTTING_CELL_N)%F_jm=-1
!	  end select
!    endif
!  enddo
END FUNCTION

! area3D_Polygon(): computes the area of a 3D planar polygon
!    Input:  int n = the number of vertices in the polygon
!            Point* V = an array of n+2 vertices in a plane
!                       with V[n]=V[0] and V[n+1]=V[1]
!            Point N = unit normal vector of the polygon's plane
!    Return: the (float) area of the polygon
FUNCTION area3D_Polygon( n, A, Norm )
  USE global_parameters
  USE unstructured_surface_arrays
  
  REAL(KIND=CGREAL) ::  area, area3D_Polygon
  REAL(KIND=CGREAL) ::  an, ax, ay, az  ! abs value of normal and its coords
  REAL(KIND=CGREAL)  :: A(3, n+2)
  REAL(KIND=CGREAL), DIMENSION(3)  :: Norm
  integer   coord              ! coord to ignore: 1=x, 2=y, 3=z
  integer   i, j, k, n            ! loop indices

  ! select largest abs coordinate to ignore for projection
  IF (Norm(1)>0) THEN ! abs x-coord
    ax = Norm(1)
  else
    ax = -Norm(1)
  endif
  IF (Norm(2)>0) THEN ! abs y-coord
    ay = Norm(2)
  else
    ay = -Norm(2)
  endif
  IF (Norm(3)>0) THEN ! abs z-coord
    az = Norm(3)
  else
    az = -Norm(3)
  endif

  coord = 3                          ! ignore z-coord
  if (ax > ay) then
     if (ax > az) coord = 1    ! ignore x-coord
  else if (ay > az) then
     coord = 2   ! ignore y-coord
  endif
  area = zero
  ! compute area of the 2D projection
  j=3
  k=1
  do i=2, n+1
      select case (coord)
      case (1)
          area = area + (A(2,i) * (A(3,j) - A(3,k)))
      case (2)
          area = area + (A(1,i) * (A(3,j) - A(3,k)))
      case (3)
          area = area + (A(1,i) * (A(2,j) - A(2,k)))
      end select
      j=j+1
      k=k+1
  enddo

  select case (coord)
  case (1)
      area = area / (2*ax)
  case (2)
      area = area / (2*ay)
  case (3)
      area = area / (2*az)
  end select

  area3D_Polygon = area
END FUNCTION

FUNCTION area3D_Polygon_Normalized( n, A, Norm, dxinv, dyinv, dzinv )
  USE global_parameters
  USE unstructured_surface_arrays

  REAL(KIND=CGREAL) ::  area, area3D_Polygon_Normalized
  REAL(KIND=CGREAL) ::  an, ax, ay, az  ! abs value of normal and its coords
  REAL(KIND=CGREAL) ::  dxinv, dyinv, dzinv
  REAL(KIND=CGREAL)  :: A(3, n+2)
  REAL(KIND=CGREAL), DIMENSION(3)  :: Norm
  integer   coord              ! coord to ignore: 1=x, 2=y, 3=z
  integer   i, j, k, n            ! loop indices

  ! select largest abs coordinate to ignore for projection
  IF (Norm(1)>0) THEN ! abs x-coord
    ax = Norm(1)
  else
    ax = -Norm(1)
  endif
  IF (Norm(2)>0) THEN ! abs y-coord
    ay = Norm(2)
  else
    ay = -Norm(2)
  endif
  IF (Norm(3)>0) THEN ! abs z-coord
    az = Norm(3)
  else
    az = -Norm(3)
  endif

  coord = 3                          ! ignore z-coord
  if (ax > ay) then
     if (ax > az) coord = 1    ! ignore x-coord
  else if (ay > az) then
     coord = 2   ! ignore y-coord
  endif
  area = zero
  ! compute area of the 2D projection
  j=3
  k=1
  do i=2, n+1
      select case (coord)
      case (1)
          area = area + (A(2,i) * dyinv * (A(3,j) - A(3,k))) * dzinv
      case (2)
          area = area + (A(1,i) * dxinv * (A(3,j) - A(3,k))) * dzinv
      case (3)
          area = area + (A(1,i) * dxinv * (A(2,j) - A(2,k))) * dyinv
      end select
      j=j+1
      k=k+1
  enddo

  select case (coord)
  case (1)
      area = area / (2*ax) * dxinv
  case (2)
      area = area / (2*ay) * dyinv
  case (3)
      area = area / (2*az) * dzinv
  end select

  area3D_Polygon_Normalized = area
END FUNCTION

SUBROUTINE POLYGON_CENTROID_3D(n, A, CENTROID)

  USE global_parameters
  USE unstructured_surface_arrays

  IMPLICIT NONE

  REAL(KIND=CGREAL) ::  area, areatri
  REAL(KIND=CGREAL) ::  centroid_x,centroid_y,centroid_z

  REAL(KIND=CGREAL)  :: A(3, n)
  REAL(KIND=CGREAL), DIMENSION(3)  :: CENTROID, V1, V2, NORM

  integer   i, n         ! loop indices

  centroid_x = zero
  centroid_y = zero
  centroid_z = zero
  area = zero

  DO i = 1, n-2
     V1=A(:,i+1)-A(:,i)
     V2=A(:,n)-A(:,i)
     
     Norm=Cross(V1,V2)
     areatri=half*SQRT(DOT_PRODUCT(NORM,NORM))
     
     area = area + areatri
     centroid_x = centroid_x + areatri*(A(1,i)+A(1,i+1)+A(1,n) )/3.0_CGREAL
     centroid_y = centroid_y + areatri*(A(2,i)+A(2,i+1)+A(2,n) )/3.0_CGREAL
     centroid_z = centroid_z + areatri*(A(3,i)+A(3,i+1)+A(3,n) )/3.0_CGREAL
  ENDDO

  CENTROID(1) = centroid_x/area
  CENTROID(2) = centroid_y/area
  CENTROID(3) = centroid_z/area

END SUBROUTINE POLYGON_CENTROID_3D

SUBROUTINE POLYGON_CENTROID_3D_Normalized(n, A, dd, CENTROID)

  USE global_parameters
  USE unstructured_surface_arrays

  IMPLICIT NONE

  REAL(KIND=CGREAL) ::  area, areatri
  REAL(KIND=CGREAL)  :: A(3, n)
  REAL(KIND=CGREAL), DIMENSION(3)  :: CENTROID, V1, V2, NORM, dd, centroid_x

  integer   i, n         ! loop indices

  centroid_x = zero
  area = zero

  DO i = 1, n-2
     V1=(A(:,i+1)-A(:,i)) * dd
     V2=(A(:,n)-A(:,i)) * dd

     Norm=Cross(V1,V2) * dd
     areatri=half*SQRT(DOT_PRODUCT(NORM,NORM))

     area = area + areatri
     centroid_x = centroid_x + areatri*(A(:,i)+A(:,i+1)+A(:,n) )/3.0_CGREAL
  ENDDO

  CENTROID = centroid_x/area

END SUBROUTINE POLYGON_CENTROID_3D_Normalized
!===================================================================
FUNCTION CALC_CUT_FACE(I,J,K,norm)
  USE global_parameters
  USE unstructured_surface_arrays
  USE flow_arrays
  USE boundary_arrays
  USE grid_arrays
  
  implicit none

  integer :: i, j, k, n, ii, i1, j1, k1
  REAL(KIND=CGREAL)  :: B(3, 7)
  REAL(KIND=CGREAL), DIMENSION(3)  :: Norm, V1, V2, dd
  REAL(KIND=CGREAL), DIMENSION(3)  :: POLYGON_CENTROID_3D
!  TYPE (point_pair) :: A(5)
  REAL(KIND=CGREAL) :: area3D_Polygon, area3D_Polygon_Normalized, S
  INTEGER :: construct_slice, CALC_CUT_FACE, construct_slice_2
  REAL(KIND=CGREAL) ::  an, ax, ay, az  ! abs value of normal and its coords
  integer   coord              ! coord to ignore: 1=x, 2=y, 3=z
  
!  n = construct_slice(i, j, k, A)
  n = construct_slice_2(i, j, k, B)
  if (n<3) then 
  endif
!  do ii=1, n
!    B(ii)=A(ii)%pv
!  enddo
  B(:,n+1)=B(:,1)
  B(:,n+2)=B(:,2)
  
!  V1%x=B(1)%x-B(2)%x
!  V2%x=B(3)%x-B(2)%x
!  
!  NORM=CROSS(V1,V2)
  dd=(/dxinv(i), dyinv(j), dzinv(k)/)

  V1=(B(:,2)-B(:,1))*dd
  V2=(B(:,3)-B(:,1))*dd

  Norm=Cross(V1,V2)*dd
  
  S=SQRT(DOT_PRODUCT(NORM,NORM))
  NORM=NORM/S
  
  i1=Nint(i+Norm(1))
  j1=Nint(j+Norm(2))
  k1=Nint(k+Norm(3))
  IF (iblank(i1, j1, k1)==0) norm=-norm

  CUTTING_FACE_N=CUTTING_FACE_N+1
!  face(CUTTING_FACE_N)%area=abs(area3D_Polygon( n, B, Norm ))
!  face(CUTTING_FACE_N)%a=abs(area3D_Polygon( n, B, Norm ))
  face(CUTTING_FACE_N)%a=abs(area3D_Polygon_Normalized( n, B, Norm,dxinv(i),dyinv(j),dzinv(k) ))
!  face(CUTTING_FACE_N)%centroid=POLYGON_CENTROID_3D(n, B)
  CALL POLYGON_CENTROID_3D_Normalized(n, B, dd, face(CUTTING_FACE_N)%centroid)
  
  CALC_CUT_FACE=CUTTING_FACE_N
  
  IF ( CELL(CUTTING_CELL_N)%F_im==0 ) THEN
    IF (Norm(1)<ZERO) THEN
	  CELL(CUTTING_CELL_N)%F_im=-1
	ELSE
	  CUTTING_FACE_N=CUTTING_FACE_N+1
	  CELL(CUTTING_CELL_N)%F_im=CUTTING_FACE_N
	  face(CUTTING_FACE_N)%a=oned
	  face(CUTTING_FACE_N)%centroid=(/x(i),yc(j),zc(k)/)
	ENDIF
  ENDIF
  IF ( CELL(CUTTING_CELL_N)%F_ip==0 ) THEN
    IF (Norm(1)>ZERO) THEN
      CELL(CUTTING_CELL_N)%F_ip=-1
	ELSE
	  CUTTING_FACE_N=CUTTING_FACE_N+1
	  CELL(CUTTING_CELL_N)%F_ip=CUTTING_FACE_N
	  face(CUTTING_FACE_N)%a=oned
	  face(CUTTING_FACE_N)%centroid=(/x(i+1),yc(j),zc(k)/)
	ENDIF
  ENDIF
  IF ( CELL(CUTTING_CELL_N)%F_jm==0 ) THEN
    IF (Norm(2)<ZERO) THEN
      CELL(CUTTING_CELL_N)%F_jm=-1
	ELSE
	  CUTTING_FACE_N=CUTTING_FACE_N+1
	  CELL(CUTTING_CELL_N)%F_jm=CUTTING_FACE_N
	  face(CUTTING_FACE_N)%a=oned
	  face(CUTTING_FACE_N)%centroid=(/xc(i),y(j),zc(k)/)
	ENDIF
  ENDIF
  IF ( CELL(CUTTING_CELL_N)%F_jp==0 ) THEN
    IF (Norm(2)>ZERO) THEN
      CELL(CUTTING_CELL_N)%F_jp=-1
	ELSE
	  CUTTING_FACE_N=CUTTING_FACE_N+1
	  CELL(CUTTING_CELL_N)%F_jp=CUTTING_FACE_N
	  face(CUTTING_FACE_N)%a=oned
	  face(CUTTING_FACE_N)%centroid=(/xc(i),y(j+1),zc(k)/)
	ENDIF
  ENDIF
  IF ( CELL(CUTTING_CELL_N)%F_km==0 ) THEN
    IF (Norm(3)<ZERO) THEN
      CELL(CUTTING_CELL_N)%F_km=-1
	ELSE
	  CUTTING_FACE_N=CUTTING_FACE_N+1
	  CELL(CUTTING_CELL_N)%F_km=CUTTING_FACE_N
	  face(CUTTING_FACE_N)%a=oned
	  face(CUTTING_FACE_N)%centroid=(/xc(i),yc(j),z(k)/)
	ENDIF
  ENDIF
  IF ( CELL(CUTTING_CELL_N)%F_kp==0 ) THEN
    IF (Norm(3)>ZERO) THEN
      CELL(CUTTING_CELL_N)%F_kp=-1
	ELSE
	  CUTTING_FACE_N=CUTTING_FACE_N+1
	  CELL(CUTTING_CELL_N)%F_kp=CUTTING_FACE_N
	  face(CUTTING_FACE_N)%a=oned
	  face(CUTTING_FACE_N)%centroid=(/xc(i),yc(j),z(k+1)/)
	ENDIF
  ENDIF

END FUNCTION

FUNCTION CALC_CELL_VOLUMN(I,J,K)
  USE global_parameters
  USE unstructured_surface_arrays
  USE boundary_arrays
  USE grid_arrays
  USE flow_arrays
  
  implicit none

  integer :: i, j, k, n, fc
  REAL(KIND=CGREAL) :: CALC_CELL_VOLUMN, Volumn
  REAL(KIND=CGREAL) :: Ae, Aw, An, As, Af, Ab, Ag
  REAL(KIND=CGREAL), DIMENSION(3)  :: Xe, Xw, Xn, Xs, Xf, Xb, Xg
  REAL(KIND=CGREAL), DIMENSION(3)  :: Ne, Nw, Nn, Ns, Nf, Nb, Ng
  REAL(KIND=CGREAL) :: dote, dotw, dotn, dots, dotf, dotb, dotg
  
  DATA Ne /oned, zero, zero/
  DATA Nw /-oned, zero, zero/
  DATA Nn /zero, oned, zero/
  DATA Ns /zero, -oned, zero/
  DATA Nf /zero, zero, oned/
  DATA Nb /zero, zero, -oned/

  fc=cell(cutting_cell_N)%F_ip
  IF (fc==-1) THEN
    Ae=ZERO
    Xe=ZERO
  ELSE
    Ae=face(fc)%a*dxinv(i)
    Xe=face(fc)%centroid
  ENDIF
  
  fc=cell(cutting_cell_N)%F_im
  IF (fc==-1) THEN
	  Aw=ZERO
    Xw=ZERO
  ELSE
	Aw=face(fc)%a*dxinv(i)
	Xw=face(fc)%centroid
  ENDIF
  
  fc=cell(cutting_cell_N)%F_jp
  IF (fc==-1) THEN
	  An=ZERO
    Xn=ZERO
  ELSE
    An=face(fc)%a*dyinv(j)
    Xn=face(fc)%centroid
  ENDIF
  
  fc=cell(cutting_cell_N)%F_jm
  IF (fc==-1) THEN
    As=zero
    Xs=zero
  ELSE
    As=face(fc)%a*dyinv(j)
    Xs=face(fc)%centroid
  ENDIF
  
  fc=cell(cutting_cell_N)%F_kp
  IF (fc==-1) THEN
    Af=zero
    Xf=zero
  ELSE
    Af=face(fc)%a*dzinv(k)
    Xf=face(fc)%centroid
  ENDIF
  
  fc=cell(cutting_cell_N)%F_km
  IF (fc==-1) THEN
	  Ab=zero
    Xb=zero
  ELSE
	  Ab=face(fc)%a*dzinv(k)
	  Xb=face(fc)%centroid
  ENDIF

  Ag=face(cell(cutting_cell_N)%F_slice)%a
  Xg=face(cell(cutting_cell_N)%F_slice)%centroid
  
  Ng=cell(cutting_cell_N)%slice_normal
  
  dote=DOT_PRODUCT(Xe, Ne)
  dotw=DOT_PRODUCT(Xw, Nw)
  dotn=DOT_PRODUCT(Xn, Nn)
  dots=DOT_PRODUCT(Xs, Ns)
  dotf=DOT_PRODUCT(Xf, Nf)
  dotb=DOT_PRODUCT(Xb, Nb)
  dotg=DOT_PRODUCT(Xg, Ng)
  Volumn=oned/3*(dote*Ae + dotw*Aw +  &
                 dotn*An + dots*As +  &
                 dotf*Af + dotb*Ab +  &
                 dotg*Ag)
  IF (Volumn<-1E-10 .or. Volumn>(oned+1E-10)) THEN
    WRITE(*,*) 'Negtive cut volumn in cell', i,j,k
!    CALL write_dump()
    CALL write_dump_debug_i('ibln',0,iblank)
    WRITE(*,*) x(i), y(j), z(k)
    WRITE(*,*) x(i+1), y(j+1), z(k+1)
    WRITE(*,*) xc(i), yc(j), zc(k)
    WRITE(*,*) dx(i), dy(j), dz(k)
    WRITE(*,*) dxinv(i), dyinv(j), dzinv(k)
    WRITE(*,*) cell(cutting_cell_N)
    write(*,*) edge_cut(iex(i,j,k)), edge_cut(iex(i,j,k+1))
    write(*,*) edge_cut(iex(i,j+1,k)), edge_cut(iex(i,j+1,k+1))
    write(*,*) x(i+1)-edge_cut(iex(i,j,k))
    write(*,*) x(i+1)-edge_cut(iex(i,j+1,k))
    write(*,*) (x(i+1)-edge_cut(iex(i,j+1,k))+x(i+1)-edge_cut(iex(i,j,k)))/twod
    WRITE(*,*) dote, dotw
    WRITE(*,*) dotn, dots
    WRITE(*,*) dotf, dotb
    WRITE(*,*) dotg
    WRITE(*,*) Xe, Ne, Ae
    WRITE(*,*) Xw, Nw, Aw
    WRITE(*,*) Xn, Nn, An
    WRITE(*,*) Xs, Ns, As
    WRITE(*,*) Xf, Nf, Af
    WRITE(*,*) Xb, Nb, Ab
    WRITE(*,*) Xg, Ng, Ag
    WRITE(*,*) Volumn
    stop
  ENDIF
  CALC_CELL_VOLUMN=Volumn
  
END FUNCTION

!-------------------------------------------------------------------------
subroutine get_inter(norm,ref,p1,p2,pout)
implicit none

real(8),intent(in) :: norm(3),ref(3),p1(3),p2(3)
real(8),intent(out) :: pout(3)

real*8 :: lhs1,lhs2

lhs1=norm(1)*(p1(1)-ref(1))+norm(2)*(p1(2)-ref(2))+norm(3)*(p1(3)-ref(3))
lhs2=norm(1)*(p2(1)-ref(1))+norm(2)*(p2(2)-ref(2))+norm(3)*(p2(3)-ref(3))
if(abs(lhs1)<1e-10)then
    pout=p1
else if(abs(lhs2)<1e-10)then
    pout=p2
else
    pout(1)=(norm(1)*ref(1)*p1(1) - norm(1)*ref(1)*p2(1) + norm(2)*p1(1)*ref(2) - norm(2)*p2(1)*ref(2) - norm(2)*p1(1)*p2(2) + norm(2)*p2(1)*p1(2) + norm(3)*p1(1)*ref(3) - norm(3)*p2(1)*ref(3) - norm(3)*p1(1)*p2(3) + norm(3)*p2(1)*p1(3))/(norm(1)*p1(1) - norm(1)*p2(1) + norm(2)*p1(2) - norm(2)*p2(2) + norm(3)*p1(3) - norm(3)*p2(3))
    pout(2)=(norm(1)*ref(1)*p1(2) - norm(1)*ref(1)*p2(2) + norm(1)*p1(1)*p2(2) - norm(1)*p2(1)*p1(2) + norm(2)*ref(2)*p1(2) - norm(2)*ref(2)*p2(2) + norm(3)*p1(2)*ref(3) - norm(3)*p2(2)*ref(3) - norm(3)*p1(2)*p2(3) + norm(3)*p2(2)*p1(3))/(norm(1)*p1(1) - norm(1)*p2(1) + norm(2)*p1(2) - norm(2)*p2(2) + norm(3)*p1(3) - norm(3)*p2(3))
    pout(3)=(norm(1)*ref(1)*p1(3) - norm(1)*ref(1)*p2(3) + norm(1)*p1(1)*p2(3) - norm(1)*p2(1)*p1(3) + norm(2)*ref(2)*p1(3) - norm(2)*ref(2)*p2(3) + norm(2)*p1(2)*p2(3) - norm(2)*p2(2)*p1(3) + norm(3)*ref(3)*p1(3) - norm(3)*ref(3)*p2(3))/(norm(1)*p1(1) - norm(1)*p2(1) + norm(2)*p1(2) - norm(2)*p2(2) + norm(3)*p1(3) - norm(3)*p2(3))
end if

end subroutine get_inter
    
subroutine get_distance_face(norm,ref,p,dist)
use operation
implicit none

real(8),intent(in) :: norm(3),ref(3),p(3)
real(8),intent(out) :: dist

real(8) :: vector_p(3),angle

if(p(1)==ref(1).and.p(2)==ref(2).and.p(3)==ref(3))then
    dist=0.0d0
else

    vector_p=p-ref
    angle=vector_angle(norm,vector_p)
    
    dist=mo(vector_p,3)*cos(angle)

end if

end subroutine get_distance_face
    
subroutine triangle_area(n1Coord,n2Coord,n3Coord,areaOut)
use operation
implicit none

real*8,intent(in) :: n1Coord(3),n2Coord(3),n3Coord(3)
real*8,intent(out) :: areaOut

real*8 :: vec12(3),vec13(3)
real*8 :: angle1,length12,length13

vec12=n2Coord-n1Coord
vec13=n3Coord-n1Coord

angle1=vector_angle(vec12,vec13)
length12=mo(vec12,3)
length13=mo(vec13,3)

areaOut=0.5d0*length12*length13*sin(angle1)

end subroutine triangle_area


!--------------------------------
subroutine check_if_edge_marker(numBd,numMarker,edgeM,iCell,jCell,kCell,dotNorm)
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    use operation

IMPLICIT NONE

integer,intent(in) :: numBd,numMarker,iCell,jCell,kCell
logical,intent(out) :: edgeM
REAL(KIND=CGREAL) , INTENT(OUT) :: dotNorm

integer :: m,n
integer :: neighElemInd(1000)
integer :: numNeighElement
logical :: sharp

real(CGREAL) :: eNorm(3),eNormRef(3),eCentRef(3),cCent(3),eCent(3),vecAngle,dist

numNeighElement = 0
DO m=1,totNumTriElem(numBd)
    IF (triElemNeig(numBd,1,m) == numMarker .OR. &
        triElemNeig(numBd,2,m) == numMarker .OR. &
        triElemNeig(numBd,3,m) == numMarker) THEN
        numNeighElement = numNeighElement + 1
        NeighElemInd(numNeighElement) = m
    ENDIF
ENDDO ! m

do m=1,numNeighElement-1
do n=m+1,numNeighElement
    eNormRef(1)=triElemNormx(numBd,NeighElemInd(m))
    eNormRef(2)=triElemNormy(numBd,NeighElemInd(m))
    eNormRef(3)=triElemNormz(numBd,NeighElemInd(m))
    
    eNorm(1)=triElemNormx(numBd,NeighElemInd(n))
    eNorm(2)=triElemNormy(numBd,NeighElemInd(n))
    eNorm(3)=triElemNormz(numBd,NeighElemInd(n))
    eNormRef=-eNormRef
    eNorm=-eNorm

    vecAngle=vector_angle(eNormRef,eNorm)
    if(vecAngle>=pi/2.0d0)then
        edgeM=.true.
        return
        !goto 888
    else
        edgeM=.false.
    end if
    
end do
end do

!888 continue
!if(edgeM)then
!    
!    eCentRef(1)=triElemCentx(numBd,NeighElemInd(m))
!    eCentRef(2)=triElemCenty(numBd,NeighElemInd(m))
!    eCentRef(3)=triElemCentz(numBd,NeighElemInd(m))
!
!    eNormRef(1)=triElemNormx(numBd,NeighElemInd(m))
!    eNormRef(2)=triElemNormy(numBd,NeighElemInd(m))
!    eNormRef(3)=triElemNormz(numBd,NeighElemInd(m))
!
!    eCent(1)=triElemCentx(numBd,NeighElemInd(n))
!    eCent(2)=triElemCenty(numBd,NeighElemInd(n))
!    eCent(3)=triElemCentz(numBd,NeighElemInd(n))
!
!    call get_distance_face(eNormRef,eCentRef,eCent,dist)
!    if(dist>0.0d0)then
!        sharp=.true.
!        dotNorm=-1.0d0
!        return
!    else
!        sharp=.false.
!        edgeM=.false.
!        return
!    end if
!end if

!if(edgeM)then
!    do m=1,numNeighElement
!    
!        eCentRef(1)=triElemCentx(numBd,NeighElemInd(m))
!        eCentRef(2)=triElemCenty(numBd,NeighElemInd(m))
!        eCentRef(3)=triElemCentz(numBd,NeighElemInd(m))
!    
!        eNormRef(1)=triElemNormx(numBd,NeighElemInd(m))
!        eNormRef(2)=triElemNormy(numBd,NeighElemInd(m))
!        eNormRef(3)=triElemNormz(numBd,NeighElemInd(m))
!        
!        cCent(1)=xc(iCell)
!        cCent(2)=yc(jCell)
!        cCent(3)=zc(kCell)
!            
!        call get_distance_face(eNormRef,eCentRef,cCent,dist)
!        
!        if(sharp)then
!            if(dist>0.0d0)then
!                dotNorm=1.0d0
!            else
!                dotNorm=-1.0d0
!                return
!            end if
!        else
!            if(dist<0.0d0)then
!                dotNorm=-1.0d0
!            else
!                dotNorm=1.0d0
!                return
!            end if
!            
!        end if
!        
!            
!    
!    end do
!end if

end subroutine check_if_edge_marker

!---------------------------------------------------
subroutine identify_conflict_cell
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE GCM_arrays
USE unstructured_surface_arrays

IMPLICIT NONE

integer :: i,j,k,e,iBody
REAL(KIND=CGREAL) :: vert1(3),vert2(3),vert3(3)
REAL(KIND=CGREAL) :: xMin,xMax,yMin,yMax,zMin,zMax,dotNorm
INTEGER :: iMin,iMax,jMin,jMax,kMin,kMax
INTEGER :: iM,iP,jM,jP,kM,kP
integer :: cElement


INTEGER(1) :: boundCell(0:nx+1,0:ny+1,0:nz+1)

if(nBody_solid/=0)then
    do k=1,nzc
    do j=1,nyc
    do i=1,nxc
        if(iblank(i,j,k)==0)then
            if(iblank(i,j,k)+iblank(i-1,j,k)==1.or.&
               iblank(i,j,k)+iblank(i+1,j,k)==1.or.&
               iblank(i,j,k)+iblank(i,j-1,k)==1.or.&
               iblank(i,j,k)+iblank(i,j+1,k)==1.or.&
               iblank(i,j,k)+iblank(i,j,k-1)==1.or.&
               iblank(i,j,k)+iblank(i,j,k+1)==1)then
                
                conflictCell(i,j,k)=1

                if(iblank(i,j,k)+iblank(i-1,j,k)==1) conflictBCi(i,j,k)=1
                if(iblank(i,j,k)+iblank(i,j-1,k)==1) conflictBCj(i,j,k)=1
                if(iblank(i,j,k)+iblank(i,j,k-1)==1) conflictBCk(i,j,k)=1
                
                if(iblank(i,j,k)+iblank(i+1,j,k)==1) conflictBCi(i+1,j,k)=1
                if(iblank(i,j,k)+iblank(i,j+1,k)==1) conflictBCj(i,j+1,k)=1
                if(iblank(i,j,k)+iblank(i,j,k+1)==1) conflictBCk(i,j,k+1)=1
        
            end if
        end if
    end do
    end do
    end do
end if


if(nBody_membrane/=0)then
do iBody=nBody_solid+1,nBody
    boundCell=0
    do e=1,totNumTriElem(iBody)
        
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
    end do
    
    do k=1,nzc
    do j=1,nyc
    do i=1,nxc

        if(boundCell(i,j,k)==1)then
            call search_vertex_dotNorm(iBody,i,j,k,cElement,dotNorm)

            if(dotNorm>=zero)then
                boundCell(i,j,k)=-1
            else
                boundCell(i,j,k)=1
            end if
        end if
        
    end do
    end do
    end do
    
    do k=1,nzc
    do j=1,nyc
    do i=1,nxc
        
        if(abs(boundCell(i,j,k))==1)then
        
        if((boundCell(i,j,k)*boundCell(i-1,j,k)==-1.and.abs(boundCell(i-1,j,k))==1).or.&
        (boundCell(i,j,k)*boundCell(i+1,j,k)==-1.and.abs(boundCell(i+1,j,k))==1).or.&
        (boundCell(i,j,k)*boundCell(i,j-1,k)==-1.and.abs(boundCell(i,j-1,k))==1).or.&
        (boundCell(i,j,k)*boundCell(i,j+1,k)==-1.and.abs(boundCell(i,j+1,k))==1).or.&
        (boundCell(i,j,k)*boundCell(i,j,k-1)==-1.and.abs(boundCell(i,j,k-1))==1).or.&
        (boundCell(i,j,k)*boundCell(i,j,k+1)==-1.and.abs(boundCell(i,j,k+1))==1))then

            conflictCell(i,j,k)=conflictCell(i,j,k)+2**(iBody-1)

            if(boundCell(i,j,k)*boundCell(i-1,j,k)==-1) conflictBCi(i,j,k)=conflictBCi(i,j,k)+2**(iBody-1)
            if(boundCell(i,j,k)*boundCell(i,j-1,k)==-1) conflictBCj(i,j,k)=conflictBCj(i,j,k)+2**(iBody-1)
            if(boundCell(i,j,k)*boundCell(i,j,k-1)==-1) conflictBCk(i,j,k)=conflictBCk(i,j,k)+2**(iBody-1)
        end if
        end if

        if(iblank(i,j,k)==1) conflictCell(i,j,k)=0
        
    end do
    end do
    end do

end do
end if


end subroutine identify_conflict_cell
    
    
!------------------------------------------------
subroutine pointInTriangle(n1Coord,n2Coord,n3Coord,pOut,ifPointIn)
use global_parameters
use operation
implicit none

REAL(KIND=CGREAL),intent(in) :: n1Coord(3),n2Coord(3),n3Coord(3),pOut(3)
logical,intent(out) :: ifPointIn 
REAL(KIND=CGREAL) :: pi=4.0_CGREAL*atan(1.0_CGREAL)

REAL(KIND=CGREAL) :: vec1(3),vec2(3),vec3(3),vec11(3),vec22(3),vec33(3)
REAL(KIND=CGREAL) :: c1(3),c2(3),c3(3),moC1,moC2,moC3,a12,a23

vec1=n2Coord-n1Coord
vec2=n3Coord-n2Coord
vec3=n1Coord-n3Coord

vec11=pOut-n1Coord
vec22=pOut-n2Coord
vec33=pOut-n3Coord

c1=cross(vec1,vec11)
c2=cross(vec2,vec22)
c3=cross(vec3,vec33)

moC1=mo(c1,3)
moC2=mo(c2,3)
moC3=mo(c3,3)

if(moC1<1e-20.or.moC2<1e-20.or.moC3<1e-20)then
    ifPointIn=.false.
    return
else

    c1=c1/moC1
    c2=c2/moC1
    c3=c3/moC1

    a12=vector_angle(c1,c2)
    a23=vector_angle(c2,c3)

    if(a12<pi/10.0_CGREAL.and.a23<pi/10.0_CGREAL)then
        ifPointIn=.true.
    else
        ifPointIn=.false.
    end if
    return
end if

end subroutine pointInTriangle