!******************************************************************************
!
! Purpose: generalized kernel to compute the test filtering operation
!
! Description: none.
!
! Input: field variables
!
! Output: test-filtered field variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: Exp $
!
! Copyright: (c) 2004 by the George Washington University
!
!******************************************************************************

!------------------------------------------------------------------------------
   SUBROUTINE TURB_ApplyTestFiltering( nVars, q, qF)

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE turb_parameters
    USE grid_arrays
    USE boundary_arrays
    USE turb_arrays
!    USE TURB_ModInterfaces, ONLY : TURB_TestFiltering

    IMPLICIT NONE

!... Parameters

    INTEGER, INTENT(IN) :: nVars
    REAL(KIND=CGREAL), DIMENSION(nVars,0:nx+1,0:ny+1,0:nz+1), INTENT(IN)  :: q
    REAL(KIND=CGREAL), DIMENSION(nVars,0:nx+1,0:ny+1,0:nz+1), INTENT(OUT) :: qF

!... Loop variables

    INTEGER :: iVar,i,j,k
     
!... Local variables

    INTEGER :: iErr
    REAL(KIND=CGREAL), DIMENSION(:,:,:)  , ALLOCATABLE :: qLocal3D,qFLocal3D
    REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: qFTemp

!******************************************************************************

!------------------------------------------------------------------------------
! Allocate local arrays
!------------------------------------------------------------------------------

   ALLOCATE(qLocal3D(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
   IF ( iErr /= ERR_NONE ) THEN
     WRITE(STDOUT,*) &
      'TURB_ApplyTestFiltering: Memory Allocation Error for qLocal3D'
     STOP
   ENDIF ! iErr 
   
   ALLOCATE(qFLocal3D(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
   IF ( iErr /= ERR_NONE ) THEN
     WRITE(STDOUT,*) &
       'TURB_ApplyTestFiltering: Memory Allocation Error for qFLocal3D'
     STOP
   ENDIF ! iErr    

   IF ( SUM(testFilterDir) > 1 ) THEN
     ALLOCATE(qFTemp(nVars,0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
     IF ( iErr /= ERR_NONE ) THEN
       WRITE(STDOUT,*) &
      'TURB_ApplyTestFiltering: Memory Allocation Error for qFTemp'
       STOP
     ENDIF ! iErr 
   ENDIF ! SUM(testFilterDir)    

!------------------------------------------------------------------------------
! Compute test filter at cell centers
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! i-direction only active
!------------------------------------------------------------------------------

    IF ( testFilterDir(DIRX) == ACTIVE .AND. &
         testFilterDir(DIRY) /= ACTIVE .AND. & 
         testFilterDir(DIRZ) /= ACTIVE       ) THEN

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = q(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRX, q(iVar,0:,0:,0:), qF(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRX, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qF(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

    ENDIF ! testFilterDir

!------------------------------------------------------------------------------    
! j-direction only active
!------------------------------------------------------------------------------

    IF ( testFilterDir(DIRX) /= ACTIVE .AND. &
         testFilterDir(DIRY) == ACTIVE .AND. & 
         testFilterDir(DIRZ) /= ACTIVE       ) THEN

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = q(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRY, q(iVar,0:,0:,0:), qF(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRY, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qF(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

    ENDIF ! testFilterDir

!------------------------------------------------------------------------------    
! k-direction only active
!------------------------------------------------------------------------------

    IF ( testFilterDir(DIRX) /= ACTIVE .AND. &
         testFilterDir(DIRY) /= ACTIVE .AND. & 
         testFilterDir(DIRZ) == ACTIVE       ) THEN

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = q(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRZ, q(iVar,0:,0:,0:), qF(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRZ, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qF(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

    ENDIF ! testFilterDir

!------------------------------------------------------------------------------
! ij-directions active
!------------------------------------------------------------------------------

    IF ( testFilterDir(DIRX) == ACTIVE .AND. &
         testFilterDir(DIRY) == ACTIVE .AND. & 
         testFilterDir(DIRZ) /= ACTIVE       ) THEN
 

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = q(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRX, q(iVar,0:,0:,0:), qFTemp(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRX, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qFTemp(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = qFTemp(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRY, qFTemp(iVar,0:,0:,0:), qF(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRY, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qF(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k
      ENDDO ! iVar
   
    ENDIF ! testFilterDir

!------------------------------------------------------------------------------
! ik-directions active
!------------------------------------------------------------------------------

    IF ( testFilterDir(DIRX) == ACTIVE .AND. &
         testFilterDir(DIRY) /= ACTIVE .AND. & 
         testFilterDir(DIRZ) == ACTIVE       ) THEN

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = q(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRX, q(iVar,0:,0:,0:), qFTemp(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRX, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qFTemp(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = qFTemp(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRZ, qFTemp(iVar,0:,0:,0:), qF(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRZ, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qF(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

    ENDIF ! testFilterDir

!------------------------------------------------------------------------------
! jk-directions active
!------------------------------------------------------------------------------

    IF ( testFilterDir(DIRX) /= ACTIVE .AND. &
         testFilterDir(DIRY) == ACTIVE .AND. & 
         testFilterDir(DIRZ) == ACTIVE       ) THEN

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = q(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRY, q(iVar,0:,0:,0:), qFTemp(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRY, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qFTemp(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = qFTemp(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRZ, qFTemp(iVar,0:,0:,0:), qF(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRZ, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qF(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

    ENDIF ! testFilterDir

!------------------------------------------------------------------------------
! all directions active
!------------------------------------------------------------------------------

    IF ( testFilterDir(DIRX) == ACTIVE .AND. &
         testFilterDir(DIRY) == ACTIVE .AND. & 
         testFilterDir(DIRZ) == ACTIVE       ) THEN
   
      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = q(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRX, q(iVar,0:,0:,0:), qF(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRX, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qF(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      ENDDO ! iVar

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = qF(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRY, qF(iVar,0:,0:,0:), qFTemp(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRY, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qFTemp(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k
      ENDDO ! iVar

      DO iVar =1, nVars         
        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qLocal3D(i,j,k) = qFTemp(iVar,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

!        CALL TURB_TestFiltering( DIRZ, qFTemp(iVar,0:,0:,0:), qF(iVar,0:,0:,0:) )
        CALL TURB_TestFiltering( DIRZ, qLocal3D, qFLocal3D )

        DO k = 0, nz+1
        DO j = 0, ny+1
        DO i = 0, nx+1
          qF(iVar,i,j,k) = qFLocal3D(i,j,k) 
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k
      ENDDO ! iVar

    ENDIF ! testFilterDir

!------------------------------------------------------------------------------
! Allocate local arrays
!------------------------------------------------------------------------------
 
   DEALLOCATE(qLocal3D,STAT=iErr)
   IF ( iErr /= ERR_NONE ) THEN
     WRITE(STDOUT,*) &
      'TURB_ApplyTestFiltering: Memory DeAllocation Error for qLocal3D'
     STOP
   ENDIF ! iErr 
   
   DEALLOCATE(qFLocal3D,STAT=iErr)
   IF ( iErr /= ERR_NONE ) THEN
     WRITE(STDOUT,*) &
       'TURB_ApplyTestFiltering: Memory DeAllocation Error for qFLocal3D'
     STOP
    ENDIF ! iErr     
   
   IF ( SUM(testFilterDir) > 1 ) THEN
     DEALLOCATE(qFTemp,STAT=iErr)
     IF ( iErr /= ERR_NONE ) THEN
       WRITE(STDOUT,*) &
      'TURB_ApplyTestFiltering: Memory DeAllocation Error for qFTemp'
       STOP
     ENDIF ! iErr 
   ENDIF ! SUM(testFilterDir)    

   END SUBROUTINE TURB_ApplyTestFiltering
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    SUBROUTINE TURB_TestFiltering(iDir, q, qT)

!==============================================================================
!  Purpose: Apply test filtering based on direction 
!           using trapezoidal integration rule
!==============================================================================
      
      USE global_parameters
      USE turb_global_parameters      
      USE flow_parameters
      USE turb_parameters
      USE grid_arrays
      USE boundary_arrays
      
      IMPLICIT NONE

!... Parameters

      INTEGER :: iDir
      REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(IN)  :: q
      REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(OUT) :: qT

!... Loop Variables
      
      INTEGER :: i,j,k

!... Local variables

      REAL(KIND=CGREAL) :: wc, wm, wp, wT
      REAL(KIND=CGREAL) :: qe, qw, qn, qs, qf, qb

!******************************************************************************
    
      SELECT CASE (iDir)

!------------------------------------------------------------------------------
! x-direction 
!------------------------------------------------------------------------------
        CASE(DIRX)
          DO k =1, nz-1
          DO j =1, ny-1
          DO i =1, nx-1
            wT = ( dx(i+1) +dx(i) ) *(1-iup(i,j,k)) +dx(i) *iup(i,j,k) &
               + ( dx(i-1) +dx(i) ) *(1-ium(i,j,k)) +dx(i) *ium(i,j,k)

	    wm = dx(i-1) *(1-ium(i,j,k))/wT
	    wp = dx(i+1) *(1-iup(i,j,k))/wT
	    wc = oned -wm -wp    

            qT(i,j,k) = ( wm *q(i-1,j,k) +wc *q(i,j,k) +wp *q(i+1,j,k) ) &
	              *REAL(1-iblank(i,j,k),KIND=CGREAL)
          ENDDO ! i
          ENDDO ! j
          ENDDO ! k

!------------------------------------------------------------------------------
! y-direction 
!------------------------------------------------------------------------------
      
        CASE(DIRY)
          DO k =1, nz-1
          DO j =1, ny-1
          DO i =1, nx-1
            wT = ( dy(j+1) +dy(j) ) *(1-jup(i,j,k)) +dy(j) *jup(i,j,k) &
               + ( dy(j-1) +dy(j) ) *(1-jum(i,j,k)) +dy(j) *jum(i,j,k)

            wm = dy(j-1) *(1-jum(i,j,k))/wT
            wp = dy(j+1) *(1-jup(i,j,k))/wT 
	    wc = oned -wm -wp 

            qT(i,j,k) = ( wm *q(i,j-1,k) +wc *q(i,j,k) +wp *q(i,j+1,k) ) &
	              *REAL(1-iblank(i,j,k),KIND=CGREAL)
          ENDDO ! 
          ENDDO ! j
          ENDDO ! k  

!------------------------------------------------------------------------------      
! z-direction
!------------------------------------------------------------------------------
      
        CASE(DIRZ) 
          DO k =1, nz-1
          DO j =1, ny-1
          DO i =1, nx-1
            wT = ( dz(k+1) +dz(k) ) *(1-kup(i,j,k)) +dz(k) *kup(i,j,k) &
               + ( dz(k-1) +dz(k) ) *(1-kum(i,j,k)) +dz(k) *kum(i,j,k)

            wm = dz(k-1) *(1-kum(i,j,k))/wT
            wp = dz(k+1) *(1-kup(i,j,k))/wT
            wc = oned -wm -wp

            qT(i,j,k) = ( wm *q(i,j,k-1) +wc *q(i,j,k) +wp *q(i,j,k+1) ) &
	              *REAL(1-iblank(i,j,k),KIND=CGREAL)
          ENDDO ! i
          ENDDO ! j
          ENDDO ! k  
                  
      END SELECT ! iDiR
       
    END SUBROUTINE TURB_TestFiltering
!------------------------------------------------------------------------------
