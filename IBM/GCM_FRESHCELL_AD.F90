!--------------------------------------------------------------------------
   SUBROUTINE GCM_correct_rhs_ad(mDirection,iFr,jFr,kFr,n1,rhs,var)

    USE global_parameters
    USE   flow_parameters
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

!... Parameters

    INTEGER, INTENT(IN) :: mDirection,iFr,jFr,kFr,n1
    REAL(KIND=CGREAL), DIMENSION(0:n1),                 INTENT(INOUT) :: rhs
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(IN)    :: var

!... Loop Variables

    INTEGER :: n,iRow

!... Local Variables

    INTEGER :: i,j,k

!**************************************************************
! Off diagonal elements contribution
!  First select appropriate fresh cell from list

    
    DO n = 1, nFresh
      IF ( iFr == iFresh(n) .AND. jFr == jFresh(n) .AND. kFr == kFresh(n) ) THEN

        DO iRow = 1,iRowMax
          i  = iFreshCellIndex(n) + incI(iRow)
          j  = jFreshCellIndex(n) + incJ(iRow)
          k  = kFreshCellIndex(n) + incK(iRow)

          IF ( (i==iFr-1 .AND. j==jFr-1) .OR. &     ! these are the stencil elements
               (i==iFr-1 .AND. j==jFr+1) .OR. &     !
               (i==iFr+1 .AND. j==jFr-1) .OR. &     !
               (i==iFr+1 .AND. j==jFr+1) .OR. &     ! that do not fit
               (i==iFr-1 .AND. k==kFr-1) .OR. &     !
               (i==iFr-1 .AND. k==kFr+1) .OR. &     !
               (i==iFr+1 .AND. k==kFr-1) .OR. &     ! within the normal 7-point stencil
               (i==iFr+1 .AND. k==kFr+1) .OR. &     !
               (j==jFr-1 .AND. k==kFr-1) .OR. &     !
               (j==jFr-1 .AND. k==kFr+1) .OR. &     !
               (j==jFr+1 .AND. k==kFr-1) .OR. &     !
               (j==jFr+1 .AND. k==kFr+1)      ) THEN

            SELECT CASE(mDirection)
              CASE(ICOORD)
                rhs(iFr) = rhs(iFr)+coeffGCMFreshD(iRow,n)*var(i,j,k)

              CASE(JCOORD)
                rhs(jFr) = rhs(jFr)+coeffGCMFreshD(iRow,n)*var(i,j,k)

              CASE(KCOORD)
                rhs(kFr) = rhs(kFr)+coeffGCMFreshD(iRow,n)*var(i,j,k)

            END SELECT ! mDirection

          ENDIF ! i
        ENDDO ! iRow

      ENDIF ! iFr
    ENDDO ! n


   END SUBROUTINE GCM_correct_rhs_ad
!--------------------------------------------------------------------------

   SUBROUTINE GCM_correct_res_ad(iFr,jFr,kFr,var,res)

    USE global_parameters
    USE   flow_parameters
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

!... Parameters

    INTEGER,           INTENT(IN)    :: iFr, jFr, kFr
    REAL(KIND=CGREAL), INTENT(INOUT) :: res
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(IN) :: var

!... Loop Variables

    INTEGER :: n,iRow

!... Local Variables

    INTEGER :: i,j,k

!**************************************************************
! Off diagonal elements contribution
!  First select appropriate fresh cell from list

    IF (nDim == DIM_2D) THEN
      k = kFr
    ENDIF ! nDim

    DO n = 1, nFresh

      IF ( iFr == iFresh(n) .AND. jFr == jFresh(n) .AND. kFr == kFresh(n) ) THEN

        DO iRow = 1,iRowMax

          i  = iFreshCellIndex(n) + incI(iRow)
          j  = jFreshCellIndex(n) + incJ(iRow)
          k  = kFreshCellIndex(n) + incK(iRow)

          IF ( (i==iFr-1 .AND. j==jFr-1) .OR. &     ! these are the stencil elements
               (i==iFr-1 .AND. j==jFr+1) .OR. &     !
               (i==iFr+1 .AND. j==jFr-1) .OR. &     !
               (i==iFr+1 .AND. j==jFr+1) .OR. &     ! that do not fit
               (i==iFr-1 .AND. k==kFr-1) .OR. &     !
               (i==iFr-1 .AND. k==kFr+1) .OR. &     !
               (i==iFr+1 .AND. k==kFr-1) .OR. &     ! within the normal 7-point stencil
               (i==iFr+1 .AND. k==kFr+1) .OR. &     !
               (j==jFr-1 .AND. k==kFr-1) .OR. &     !
               (j==jFr-1 .AND. k==kFr+1) .OR. &     !
               (j==jFr+1 .AND. k==kFr-1) .OR. &     !
               (j==jFr+1 .AND. k==kFr+1)      ) THEN

            res = res + coeffGCMFreshD(iRow,n)*var(i,j,k)

          ENDIF ! i

        ENDDO ! iRow

      ENDIF ! iFr

    ENDDO ! n

   END SUBROUTINE GCM_correct_res_ad

