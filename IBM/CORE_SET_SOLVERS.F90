!******************************************************************************
!
! Purpose: generalized kernel to compute the coefficients 
!          for the diffusion term
!
! Description: none.
!
! Input: grid size, total viscosity
!
! Output: amx_ad, amy_ad, amz_ad, apx_ad, apy_ad, apz_ad. 
!
!******************************************************************************
!
! $Id: Exp $
!
! Copyright: (c) 2004 by the George Washington University
!
!******************************************************************************
!------------------------------------------------------------------------------
   SUBROUTINE set_solve_ad()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays

    IMPLICIT NONE

!... Loop variables 
    
    INTEGER :: i,j,k,n,iRow
    
!... Local variables

    INTEGER :: iErr    
    INTEGER :: iFr,jFr,kFr
    REAL(KIND=CGREAL) :: rFreshCell, rnDim
    REAL(KIND=CGREAL) :: nuE,nuW,nuS,nuN,nuF,nuB
!******************************************************************************

!------------------------------------------------------------------------------ 
! Initialize coefficients 
!------------------------------------------------------------------------------ 

    amx_ad = zero
    apx_ad = zero
    
    amy_ad = zero
    apy_ad = zero
    
    amz_ad = zero
    apz_ad = zero
    
    rnDim  = REAL((ndim - DIM_2D),KIND=CGREAL)

    DO k=1,nz-1    
    DO j=1,ny-1   
    DO i=1,nx-1
      nuE = (       fx(i+1)   *viscTot(i+1,j,k)                            &
          +  ( oned-fx(i+1) ) *viscTot(i,j,k)   )*(1-iup(i,j,k))  &
          + bcxvisc(i,j,k)*iup(i,j,k)

      nuW = (       fx(i)     *viscTot(i,j,k)                              &
          +  ( oned-fx(i)   ) *viscTot(i-1,j,k) )*(1-ium(i,j,k))  &
          + bcxvisc(i,j,k)*ium(i,j,k)

      nuN = (       fy(j+1)   *viscTot(i,j+1,k)                            &
          +  ( oned-fy(j+1) ) *viscTot(i,j,k)   )*(1-jup(i,j,k))  &
          + bcyvisc(i,j,k)*jup(i,j,k)

      nuS = (       fy(j)     *viscTot(i,j,k)                              &
          +  ( oned-fy(j)   ) *viscTot(i,j-1,k) )*(1-jum(i,j,k))  &
          + bcyvisc(i,j,k)*jum(i,j,k)
        
      nuF = (       fz(k+1)   *viscTot(i,j,k+1)                            &
          +  ( oned-fz(k+1) ) *viscTot(i,j,k)   )*(1-kup(i,j,k))  &
          + bczvisc(i,j,k)*kup(i,j,k)
 
      nuB = (       fz(k)     *viscTot(i,j,k)                              &
          +  ( oned-fz(k)   ) *viscTot(i,j,k-1) )*(1-kum(i,j,k))  &
          + bczvisc(i,j,k)*kum(i,j,k)

      amx_ad(i,j,k) = ( dxcinv(i)*(1-ium(i,j,k))            &
                       +dxinv(i)*ium(i,j,k)*twod )*dxinv(i)
      apx_ad(i,j,k) = ( dxcinv(i+1)*(1-iup(i,j,k))          &
                       +dxinv(i)*iup(i,j,k)*twod )*dxinv(i)

      amx_ad(i,j,k) =- (half *dt*nuW)*amx_ad(i,j,k)
      apx_ad(i,j,k) =- (half *dt*nuE)*apx_ad(i,j,k)
       
      amy_ad(i,j,k) = ( dycinv(j)*(1-jum(i,j,k))            &
                       +dyinv(j)*jum(i,j,k)*twod )*dyinv(j)
      apy_ad(i,j,k) = ( dycinv(j+1)*(1-jup(i,j,k))          &
                       +dyinv(j)*jup(i,j,k)*twod )*dyinv(j)

      amy_ad(i,j,k) =- (half *dt*nuS)*amy_ad(i,j,k)
      apy_ad(i,j,k) =- (half *dt*nuN)*apy_ad(i,j,k)
      
      amz_ad(i,j,k) = ( dzcinv(k)*(1-kum(i,j,k))            &
                       +dzinv(k)*kum(i,j,k)*twod )*dzinv(k)
      apz_ad(i,j,k) = ( dzcinv(k+1)*(1-kup(i,j,k))          &
                       +dzinv(k)*kup(i,j,k)*twod )*dzinv(k)

      amz_ad(i,j,k) =- (half *dt*nuB)*amz_ad(i,j,k)*rnDim
      apz_ad(i,j,k) =- (half *dt*nuF)*apz_ad(i,j,k)*rnDim 
    END DO ! i
    END DO ! j
    END DO !k

!------------------------------------------------------------------------------ 
!...TAKE CARE OF FRESH CELLS 
!------------------------------------------------------------------------------ 
      
    IF ( boundary_motion == MOVING_BOUNDARY ) THEN

      SELECT CASE(boundary_formulation)
        CASE(SSM_METHOD)
!
!        Value of fresh cell is computed through interpolation from six neighbors
!        use Inversed Distance Weighted Intepolation (Shepards Method)
!
!                          1   6  [  -p    ]
!         u              = -  SUM [ h   u  ]   
!          j               H  i=1 [  ij  i ]
!                         j 
!        where 
!                              6  [  -p  ]
!         H              =    SUM [ h    ]   
!          j                  i=1 [  ij  ]
!
!         h  : distance between location i and j
!          ij
!         p  : free parameter (usually taken as 2)
!

         DO k=1,nz-1    
         DO j=1,ny-1   
         DO i=1,nx-1
           rFreshCell = REAL(fresh_cell(i,j,k),KIND=CGREAL)
           amx_ad(i,j,k) = amx_ad(i,j,k)*(oned - rFreshCell)     &
                         - (  (      dxcinv(i)**sidw  )   *(1-ium(i,j,k))       &
                             +( twod*dxcinv(i)   )**sidw  *   ium(i,j,k)    )*  &
                           rFreshCell
           apx_ad(i,j,k) = apx_ad(i,j,k)*(oned - rFreshCell)     &
                         - (  (      dxcinv(i+1)**sidw )  *(1-iup(i,j,k))       &
                             +( twod*dxcinv(i+1) )**sidw  *   iup(i,j,k)    )*  &
                           rFreshCell

           amy_ad(i,j,k) = amy_ad(i,j,k)*(oned - rFreshCell)     &
                         - (  (      dycinv(j)**sidw  )   *(1-jum(i,j,k))       &
                             +( twod*dycinv(j)   )**sidw  *   jum(i,j,k)    )*  &
                           rFreshCell
           apy_ad(i,j,k) = apy_ad(i,j,k)*(oned - rFreshCell)     &
                         - (  (      dycinv(j+1)**sidw )  *(1-jup(i,j,k))       &
                             +( twod*dycinv(j+1) )**sidw  *   jup(i,j,k)    )*  &
                           rFreshCell

           amz_ad(i,j,k) = amz_ad(i,j,k)*(oned - rFreshCell)     &
                         - (  (      dzcinv(k)**sidw  )   *(1-kum(i,j,k))       &
                             +( twod*dzcinv(k)   )**sidw  *   kum(i,j,k)    )*  &
                          rFreshCell

           apz_ad(i,j,k) = apz_ad(i,j,k)*(oned - rFreshCell)     &
                         - ( (       dzcinv(k+1)**sidw )  *(1-kup(i,j,k))       &
                             +( twod*dzcinv(k+1) )**sidw  *   kup(i,j,k)    )*  &
                          rFreshCell

          amz_ad(i,j,k) = amz_ad(i,j,k)*rnDim
          apz_ad(i,j,k) = apz_ad(i,j,k)*rnDim

         END DO ! i
         END DO ! j
         END DO ! k

!------------------------------------------------------------------------------ 
!       GCM Method
!------------------------------------------------------------------------------ 
 
        CASE(GCM_METHOD)

         PRINT*,'CALL GCM_set_freshcell()i in set_solve_ad'
         CALL GCM_set_freshcell()

         DO n=1,nFresh

            iFr = iFresh(n)
            jFr = jFresh(n)
            kFr = kFresh(n)

!------------------------------------------------------------------------------ 
!           Initialize coefficient to zero for fresh cells
!------------------------------------------------------------------------------

            amx_ad (iFr,jFr,kFr) = zero
            apx_ad (iFr,jFr,kFr) = zero
            amy_ad (iFr,jFr,kFr) = zero
            apy_ad (iFr,jFr,kFr) = zero
            amz_ad (iFr,jFr,kFr) = zero
            apz_ad (iFr,jFr,kFr) = zero

            DO iRow = 1,iRowMax
               i  = iFreshCellIndex(n) + incI(iRow)
               j  = jFreshCellIndex(n) + incJ(iRow)
               k  = kFreshCellIndex(n) + incK(iRow)

               IF ( i==iFr-1 .AND. j==jFr   .AND. k==kFr   ) amx_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr+1 .AND. j==jFr   .AND. k==kFr   ) apx_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr   .AND. j==jFr-1 .AND. k==kFr   ) amy_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr   .AND. j==jFr+1 .AND. k==kFr   ) apy_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr   .AND. j==jFr   .AND. k==kFr-1 ) amz_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr   .AND. j==jFr   .AND. k==kFr+1 ) apz_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
            ENDDO ! iRow

         ENDDO ! n

      END SELECT 

    ENDIF !  boundary_motion 

  END SUBROUTINE set_solve_ad
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

   SUBROUTINE set_solve_ad_Old()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_ad_arrays
    USE GCM_arrays

    IMPLICIT NONE
    
    INTEGER :: i,j,k,iFr,jFr,kFr,n,iRow
    REAL(KIND=CGREAL) :: rnDim

    amx_ad = zero
    apx_ad = zero
    
    amy_ad = zero
    apy_ad = zero
    
    amz_ad = zero
    apz_ad = zero
    
    rnDim  = REAL((ndim - DIM_2D),KIND=CGREAL)

    DO k=1,nz-1    
    DO j=1,ny-1   
    DO i=1,nx-1
       amx_ad(i,j,k) = ( dxcinv(i)*(1-ium(i,j,k))            &
                        +dxinv(i)*ium(i,j,k)*twod )*dxinv(i)
       apx_ad(i,j,k) = ( dxcinv(i+1)*(1-iup(i,j,k))          &
                        +dxinv(i)*iup(i,j,k)*twod )*dxinv(i)
       amx_ad(i,j,k) =- (half *dt*reinv)*amx_ad(i,j,k)
       apx_ad(i,j,k) =- (half *dt*reinv)*apx_ad(i,j,k)
       
       amy_ad(i,j,k) = ( dycinv(j)*(1-jum(i,j,k))            &
                        +dyinv(j)*jum(i,j,k)*twod )*dyinv(j)
       apy_ad(i,j,k) = ( dycinv(j+1)*(1-jup(i,j,k))          &
                        +dyinv(j)*jup(i,j,k)*twod )*dyinv(j)
       amy_ad(i,j,k) =- (half *dt*reinv)*amy_ad(i,j,k)
       apy_ad(i,j,k) =- (half *dt*reinv)*apy_ad(i,j,k)
      
       amz_ad(i,j,k) = ( dzcinv(k)*(1-kum(i,j,k))            &
                        +dzinv(k)*kum(i,j,k)*twod )*dzinv(k)
       apz_ad(i,j,k) = ( dzcinv(k+1)*(1-kup(i,j,k))          &
                        +dzinv(k)*kup(i,j,k)*twod )*dzinv(k)
       amz_ad(i,j,k) =- (half *dt*reinv)*amz_ad(i,j,k)*rnDim
       apz_ad(i,j,k) =- (half *dt*reinv)*apz_ad(i,j,k)*rnDim 
    END DO ! i
    END DO ! j
    END DO !k

!...TAKE CARE OF FRESH CELLS       
    IF ( boundary_motion == MOVING_BOUNDARY ) THEN

      SELECT CASE(boundary_formulation)
        CASE(SSM_METHOD)
!
!        Value of fresh cell is computed through interpolation from six neighbors
!        use Inversed Distance Weighted Intepolation (Shepards Method)
!
!                          1   6  [  -p    ]
!         u              = -  SUM [ h   u  ]   
!          j               H  i=1 [  ij  i ]
!                           j 
!        where 
!                              6  [  -p  ]
!         H              =    SUM [ h    ]   
!          j                  i=1 [  ij  ]
!
!         h  : distance between location i and j
!          ij
!         p  : free parameter (usually taken as 2)
!

         DO k=1,nz-1    
         DO j=1,ny-1   
         DO i=1,nx-1

          amx_ad(i,j,k) = amx_ad(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)     &
                         - (  (      dxcinv(i)**sidw  )   *(1-ium(i,j,k))       &
                             +( twod*dxcinv(i)   )**sidw  *   ium(i,j,k)    )*  &
                          REAL(fresh_cell(i,j,k),KIND=CGREAL)
          apx_ad(i,j,k) = apx_ad(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)     &
                         - (  (      dxcinv(i+1)**sidw )  *(1-iup(i,j,k))       &
                             +( twod*dxcinv(i+1) )**sidw  *   iup(i,j,k)    )*  &
                          REAL(fresh_cell(i,j,k),KIND=CGREAL)

          amy_ad(i,j,k) = amy_ad(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)     &
                         - (  (      dycinv(j)**sidw  )   *(1-jum(i,j,k))       &
                             +( twod*dycinv(j)   )**sidw  *   jum(i,j,k)    )*  &
                          REAL(fresh_cell(i,j,k),KIND=CGREAL)
          apy_ad(i,j,k) = apy_ad(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)     &
                         - (  (      dycinv(j+1)**sidw )  *(1-jup(i,j,k))       &
                             +( twod*dycinv(j+1) )**sidw  *   jup(i,j,k)    )*  &
                          REAL(fresh_cell(i,j,k),KIND=CGREAL)

          amz_ad(i,j,k) = amz_ad(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)     &
                         - (  (      dzcinv(k)**sidw  )   *(1-kum(i,j,k))       &
                             +( twod*dzcinv(k)   )**sidw  *   kum(i,j,k)    )*  &
                          REAL(fresh_cell(i,j,k),KIND=CGREAL)

          apz_ad(i,j,k) = apz_ad(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)     &
                         - ( (       dzcinv(k+1)**sidw )  *(1-kup(i,j,k))       &
                             +( twod*dzcinv(k+1) )**sidw  *   kup(i,j,k)    )*  &
                          REAL(fresh_cell(i,j,k),KIND=CGREAL)

          amz_ad(i,j,k) = amz_ad(i,j,k)*rnDim
          apz_ad(i,j,k) = apz_ad(i,j,k)*rnDim

         END DO ! i
         END DO ! j
         END DO ! k

        CASE(GCM_METHOD)

         PRINT*,'CALL GCM_set_freshcell()i in set_solve_ad'
         CALL GCM_set_freshcell()

         DO n=1,nFresh

            iFr = iFresh(n)
            jFr = jFresh(n)
            kFr = kFresh(n)

! initialize coefficient to zero for fresh cells
            amx_ad (iFr,jFr,kFr) = zero
            apx_ad (iFr,jFr,kFr) = zero
            amy_ad (iFr,jFr,kFr) = zero
            apy_ad (iFr,jFr,kFr) = zero
            amz_ad (iFr,jFr,kFr) = zero
            apz_ad (iFr,jFr,kFr) = zero

            DO iRow = 1,iRowMax
               i  = iFreshCellIndex(n) + incI(iRow)
               j  = jFreshCellIndex(n) + incJ(iRow)
               k  = kFreshCellIndex(n) + incK(iRow)

               IF ( i==iFr-1 .AND. j==jFr   .AND. k==kFr   ) amx_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr+1 .AND. j==jFr   .AND. k==kFr   ) apx_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr   .AND. j==jFr-1 .AND. k==kFr   ) amy_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr   .AND. j==jFr+1 .AND. k==kFr   ) apy_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr   .AND. j==jFr   .AND. k==kFr-1 ) amz_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
               IF ( i==iFr   .AND. j==jFr   .AND. k==kFr+1 ) apz_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
            ENDDO

         ENDDO

      END SELECT 

    ENDIF !  boundary_motion 
       
  END SUBROUTINE set_solve_ad_Old
!------------------------------------------------------------------------------
