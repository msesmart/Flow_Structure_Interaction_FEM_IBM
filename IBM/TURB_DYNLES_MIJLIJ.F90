!******************************************************************************
!
! Purpose: generalized kernel to compute L_{ij} and M_{ij} and 
!          contracting M_{ij}*M_{ij} and M_{ij}*L_{ij}
!
! Description: none.
!
! Input: u,v,w  = velocity Field
!
! Output: L_{ij}, M_{ij}, M_{ij}*M_{ij} and M_{ij}*L_{ij} 
!
! Notes: Mij is trace-free since it is based on trace-free Sij
!        When contracting Lij with Mij, LijMij becomes trace-free
!
!******************************************************************************
!
! $Id: Exp $
!
! Copyright: (c) 2004 by the George Washington University
!
!******************************************************************************

!------------------------------------------------------------------------------
   SUBROUTINE TURB_CalcMij()

!==============================================================================
!  Purpose: Comput Mij term
! Note    : Definition of M_{ij} has a negative (-) sign difference with 
!          literature as to align with Dynamic Lagrangian LES formalism 
!==============================================================================
      
    USE global_parameters
    USE turb_global_parameters
    USE turb_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays, ONLY    : iblank
    USE turb_arrays
!    USE TURB_ModInterfaces, ONLY : TURB_ApplyTestFiltering, &
!                                   TURB_CalcStrainRate,     &
!                                   TURB_CalcFilterWidth
    
    IMPLICIT NONE

!... Loop variables

    INTEGER :: i,j,k,iVar
    
!... Local variables
    
    INTEGER :: iErr, nVarsSij, nVarsVel
    INTEGER :: iFlagFilter, iFlagTestFilter
    INTEGER, DIMENSION(3) :: indexNode
    
    REAL(KIND=CGREAL) :: deltaSqr, deltaTestSqr, fTerm, fTerm2, &
                         strainMagnTest,strainMagnGrid

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:,:,:) :: q3Field, q3Face

!------------------------------------------------------------------------------
! Set dimensions 
!------------------------------------------------------------------------------

    nVarsSij = 6
    nVarsVel = 3
  
    iFlagFilter     = ACTIVE-1
    iFlagTestFilter = ACTIVE

!------------------------------------------------------------------------------
! Allocate local dynamic arrays 
!------------------------------------------------------------------------------

    ALLOCATE(q3Field(nVarsVel,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'TURB_CalcMij: Memory Allocation Error for q3Field'
       STOP
    ENDIF ! iErr 

    ALLOCATE(q3Face(nVarsVel,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'TURB_CalcMij: Memory Allocation Error for q3Face'
       STOP
    ENDIF ! iErr 

!------------------------------------------------------------------------------
! Initialize arrays 
!------------------------------------------------------------------------------

!    DO k=0, nz+1
!    DO j=0, ny+1
!    DO i=0, nx+1
!    DO iVar=1, nVarsSij
!      mij(iVar,i,j,k)        = zero
!      sij(iVar,i,j,k)        = zero
!      qTestField(iVar,i,j,k) = zero          
!    ENDDO ! iVar
!    ENDDO ! i
!    ENDDO ! j
!    ENDDO ! k
!
!    DO k=0, nz+1
!    DO j=0, ny+1
!    DO i=0, nx+1
!    DO iVar=1, nVarsVel
!      q3Field(iVar,i,j,k)    = zero
!      q3Face (iVar,i,j,k)    = zero
!    ENDDO ! iVar
!    ENDDO ! i
!    ENDDO ! j
!    ENDDO ! k

!------------------------------------------------------------------------------    
! Load grid-based values 
!------------------------------------------------------------------------------

    DO k=0, nz+1
    DO j=0, ny+1
    DO i=0, nx+1
      q3Field(1,i,j,k) = u(i,j,k)
      q3Field(2,i,j,k) = v(i,j,k)
      q3Field(3,i,j,k) = w(i,j,k)

      q3Face(1,i,j,k)  = face_u(i,j,k)
      q3Face(2,i,j,k)  = face_v(i,j,k)
      q3Face(3,i,j,k)  = face_w(i,j,k)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------  
! Evaluate strain rate on basic field 
!------------------------------------------------------------------------------  

!    WRITE(*,*) ' Calculating Strain Rate'
    CALL TURB_CalcStrainRate( q3Field(:,0:,0:,0:),q3Face(:,0:,0:,0:), &
                              sij(:,0:,0:,0:) )

!------------------------------------------------------------------------------
! Correct strain rate for ghost cells near immersed body 
!------------------------------------------------------------------------------ 
    
    IF ( boundary_formulation == GCM_METHOD ) THEN
      CALL TURB_GCM_CalcStrainRate
    ENDIF ! boundary_formulation

!------------------------------------------------------------------------------  
! Apply test filter on Strain Rate to maintain trace-free condition 
!------------------------------------------------------------------------------  

!    WRITE(*,*) ' Applying Test Filter'
    CALL TURB_ApplyTestFiltering( nVarsSij,         sij(:,0:,0:,0:), &
                                           qTestField(:,0:,0:,0:)  )

!------------------------------------------------------------------------------
! Evaluate the second term in M_{ij} 
!------------------------------------------------------------------------------

!    WRITE(*,*) ' Computing Second Term of Mij'
    DO k=1, nz-1
    DO j=1, ny-1
    DO i=1, nx-1

!------------------------------------------------------------------------------
!     Compute Strain Magnitude of Sij(u): |S(u)| = SQRT( 2 S_{ij}*S_{ij} )
!------------------------------------------------------------------------------

      strainMagnTest = twoSqrt * SQRT( qTestField(S11,i,j,k)*qTestField(S11,i,j,k) &
                                     + qTestField(S22,i,j,k)*qTestField(S22,i,j,k) &
                                     + qTestField(S33,i,j,k)*qTestField(S33,i,j,k) &    
                         + 2.0_CGREAL* qTestField(S12,i,j,k)*qTestField(S12,i,j,k) &
                         + 2.0_CGREAL* qTestField(S13,i,j,k)*qTestField(S13,i,j,k) &
                         + 2.0_CGREAL* qTestField(S23,i,j,k)*qTestField(S23,i,j,k) )

!------------------------------------------------------------------------------
!     Compute test filter width 
!------------------------------------------------------------------------------

      indexNode(DIRX) = i
      indexNode(DIRY) = j
      indexNode(DIRZ) = k

      CALL TURB_CalcFilterWidth(iFlagTestFilter,indexNode,deltaTestSqr)

!------------------------------------------------------------------------------
!     Assemble second term 
!------------------------------------------------------------------------------
      fTerm2 = -deltaTestSqr*strainMagnTest

      mij(S11,i,j,k) = fTerm2 *qTestField(S11,i,j,k)
      mij(S12,i,j,k) = fTerm2 *qTestField(S12,i,j,k)
      mij(S13,i,j,k) = fTerm2 *qTestField(S13,i,j,k)
      mij(S22,i,j,k) = fTerm2 *qTestField(S22,i,j,k)
      mij(S23,i,j,k) = fTerm2 *qTestField(S23,i,j,k)
      mij(S33,i,j,k) = fTerm2 *qTestField(S33,i,j,k)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k  

!------------------------------------------------------------------------------
!  Monitor output
!------------------------------------------------------------------------------

    IF ( MOD(ntime,nmonitor) == 0 ) THEN
      PRINT*,'-----------------------------------------------------------------'
      WRITE(*,*) 'Min-Max Vals of Sii-TestFiltered = ',&
      MINVAL(qTestField(S11,1:nx-1,1:ny-1,1:nz-1)&
            +qTestField(S22,1:nx-1,1:ny-1,1:nz-1)&
            +qTestField(S33,1:nx-1,1:ny-1,1:nz-1)),&
      MAXVAL(qTestField(S11,1:nx-1,1:ny-1,1:nz-1)&
            +qTestField(S22,1:nx-1,1:ny-1,1:nz-1)&
            +qTestField(S33,1:nx-1,1:ny-1,1:nz-1))

      PRINT*,'-----------------------------------------------------------------'
    ENDIF ! ntime

111 FORMAT(2(2X,1PE15.7))

!------------------------------------------------------------------------------
! Evaluate the first term in M_{ij} 
!------------------------------------------------------------------------------

!    DO k=0, nz+1
!    DO j=0, ny+1
!    DO i=0, nx+1
!    DO iVar=1, nVarsSij
!      qField(iVar,i,j,k)     = zero
!      qTestField(iVar,i,j,k) = zero
!    ENDDO ! iVar
!    ENDDO ! i
!    ENDDO ! j
!    ENDDO ! k

!    WRITE(*,*) ' Computing First Term of Mij'

    DO k=1, nz-1
    DO j=1, ny-1
    DO i=1, nx-1

!------------------------------------------------------------------------------
!     Compute filter width
!------------------------------------------------------------------------------

      indexNode(DIRX) = i
      indexNode(DIRY) = j
      indexNode(DIRZ) = k

      CALL TURB_CalcFilterWidth(iFlagFilter,indexNode,deltaSqr)

! Compute Strain Magnitude of Sij(u): |S(u)| = SQRT( 2 S_{ij}*S_{ij} ) --------
 
      strainMagnGrid = twoSqrt * SQRT( sij(S11,i,j,k)*sij(S11,i,j,k) &
                                     + sij(S22,i,j,k)*sij(S22,i,j,k) &
                                     + sij(S33,i,j,k)*sij(S33,i,j,k) &    
                         + 2.0_CGREAL* sij(S12,i,j,k)*sij(S12,i,j,k) &
                         + 2.0_CGREAL* sij(S13,i,j,k)*sij(S13,i,j,k) &
                         + 2.0_CGREAL* sij(S23,i,j,k)*sij(S23,i,j,k) )      
      
      fTerm = deltaSqr*strainMagnGrid

!------------------------------------------------------------------------------
!     Assemble D^2 |S| S_ij
!------------------------------------------------------------------------------
     
      qField(S11,i,j,k) =  fTerm *sij(S11,i,j,k)
      qField(S12,i,j,k) =  fTerm *sij(S12,i,j,k)
      qField(S13,i,j,k) =  fTerm *sij(S13,i,j,k)
      qField(S22,i,j,k) =  fTerm *sij(S22,i,j,k)
      qField(S23,i,j,k) =  fTerm *sij(S23,i,j,k)
      qField(S33,i,j,k) =  fTerm *sij(S33,i,j,k)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------ 
! Perform test filtering on first term of M_{ij} 
!------------------------------------------------------------------------------

    CALL TURB_ApplyTestFiltering( nVarsSij,    qField(:,0:,0:,0:), &
                                           qTestField(:,0:,0:,0:)  )

!------------------------------------------------------------------------------
! Assemble M_{ij}
!------------------------------------------------------------------------------

    DO k=1, nz-1
    DO j=1, ny-1
    DO i=1, nx-1
      mij(S11,i,j,k) = mij(S11,i,j,k) +qTestField(S11,i,j,k)
      mij(S12,i,j,k) = mij(S12,i,j,k) +qTestField(S12,i,j,k)
      mij(S13,i,j,k) = mij(S13,i,j,k) +qTestField(S13,i,j,k)
      mij(S22,i,j,k) = mij(S22,i,j,k) +qTestField(S22,i,j,k)
      mij(S23,i,j,k) = mij(S23,i,j,k) +qTestField(S23,i,j,k)
      mij(S33,i,j,k) = mij(S33,i,j,k) +qTestField(S33,i,j,k)   
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
! Apply immersed boundary conditions
!------------------------------------------------------------------------------
    
    DO k=1, nz-1
    DO j=1, ny-1
    DO i=1, nx-1 
      mij(S11,i,j,k) = 2.0_CGREAL*mij(S11,i,j,k) 
      mij(S12,i,j,k) = 2.0_CGREAL*mij(S12,i,j,k)
      mij(S13,i,j,k) = 2.0_CGREAL*mij(S13,i,j,k)
      mij(S22,i,j,k) = 2.0_CGREAL*mij(S22,i,j,k) 
      mij(S23,i,j,k) = 2.0_CGREAL*mij(S23,i,j,k) 
      mij(S33,i,j,k) = 2.0_CGREAL*mij(S33,i,j,k) 
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
!  Monitor output
!------------------------------------------------------------------------------

    IF ( MOD(ntime,nmonitor) == 0 ) THEN
      PRINT*,'-----------------------------------------------------------------'
      WRITE(STDOUT, *)' Min-Max Vals of M11 = ',&
        MINVAL(mij(S11,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(mij(S11,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of M12 = ',&
        MINVAL(mij(S12,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(mij(S12,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of M13 = ',&
        MINVAL(mij(S13,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(mij(S13,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of M22 = ',&
        MINVAL(mij(S22,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(mij(S22,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of M23 = ',&
        MINVAL(mij(S23,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(mij(S23,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of M33 = ',&
        MINVAL(mij(S33,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(mij(S33,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of MII = ',&
      MINVAL(mij(S11,1:nx-1,1:ny-1,1:nz-1)&
            +mij(S22,1:nx-1,1:ny-1,1:nz-1)&
            +mij(S33,1:nx-1,1:ny-1,1:nz-1)),&
      MAXVAL(mij(S11,1:nx-1,1:ny-1,1:nz-1)&
            +mij(S22,1:nx-1,1:ny-1,1:nz-1)&
            +mij(S33,1:nx-1,1:ny-1,1:nz-1))
  
      PRINT*,'-----------------------------------------------------------------'
    ENDIF ! ntime

!------------------------------------------------------------------------------
! Deallocate local dynamic arrays
!------------------------------------------------------------------------------

    DEALLOCATE(q3Field,STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'TURB_CalcMij: Memory Deallocation Error for q3Field'
       STOP
    ENDIF ! iErr

    DEALLOCATE(q3Face,STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'TURB_CalcMij: Memory Deallocation Error for q3Face'
       STOP
    ENDIF ! iErr

   END SUBROUTINE TURB_CalcMij
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_CalcLij()

!==============================================================================
!  Purpose: Compute Lij term
!  Note   : Definition of L_{ij} does not need to be trace free since 
!           when contracting Lij with Mij, LijMij becomes trace-free 
!==============================================================================
      
    USE global_parameters
    USE turb_global_parameters
    USE turb_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE turb_arrays
!    USE TURB_ModInterfaces, ONLY : TURB_ApplyTestFiltering
    
    IMPLICIT NONE

!... Loop variables

    INTEGER :: i,j,k,iVar

!... Local variables
    
    INTEGER :: iErr, nVarsConv, nVarsVel
    
    REAL(KIND=CGREAL) :: fTerm1, fTerm2

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:,:,:) :: q3Field, q3TestField

!------------------------------------------------------------------------------    
! Set dimensions
!------------------------------------------------------------------------------

    nVarsConv = 6
    nVarsVel  = 3

!------------------------------------------------------------------------------
! Allocate local dynamic arrays
!------------------------------------------------------------------------------

    ALLOCATE(q3Field(nVarsVel,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'TURB_CalcLij: Memory Allocation Error for q3Field'
       STOP
    ENDIF ! iErr
    
    ALLOCATE(q3TestField(nVarsVel,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'TURB_CalcLij: Memory Allocation Error for q3TestField'
      STOP
    ENDIF ! iErr

!------------------------------------------------------------------------------
! Initialize arrays
!------------------------------------------------------------------------------

!    DO k=0, nz+1
!    DO j=0, ny+1
!    DO i=0, nx+1
!    DO iVar=1, nVarsConv
!          qField(iVar,i,j,k) = zero
!      qTestField(iVar,i,j,k) = zero
!    ENDDO ! iVar
!    ENDDO ! i
!    ENDDO ! j
!    ENDDO ! k
!
!    DO k=0, nz+1
!    DO j=0, ny+1
!    DO i=0, nx+1
!    DO iVar=1, nVarsVel
!        q3Field(iVar,i,j,k) = zero
!    q3TestField(iVar,i,j,k) = zero
!    ENDDO ! iVar
!    ENDDO ! i
!    ENDDO ! j
!    ENDDO ! k

!------------------------------------------------------------------------------
! Evaluate the six components of u_{ij}*u_{ij}
!------------------------------------------------------------------------------

    DO k=1, nz-1
    DO j=1, ny-1
    DO i=1, nx-1
      qField(S11,i,j,k) = u(i,j,k)*u(i,j,k)
      qField(S12,i,j,k) = u(i,j,k)*v(i,j,k)
      qField(S13,i,j,k) = u(i,j,k)*w(i,j,k)
      qField(S22,i,j,k) = v(i,j,k)*v(i,j,k)
      qField(S23,i,j,k) = v(i,j,k)*w(i,j,k)
      qField(S33,i,j,k) = w(i,j,k)*w(i,j,k)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
! Perform test filtering on velocity field
!------------------------------------------------------------------------------

    CALL TURB_ApplyTestFiltering( nVarsConv,    qField(:,0:,0:,0:), &
                                            qTestField(:,0:,0:,0:)  )

!------------------------------------------------------------------------------
! Evaluate the first term in L_{ij} 
!------------------------------------------------------------------------------

    DO k=1, nz-1
    DO j=1, ny-1
    DO i=1, nx-1
      lij(S11,i,j,k) = qTestField(S11,i,j,k)
      lij(S12,i,j,k) = qTestField(S12,i,j,k)
      lij(S13,i,j,k) = qTestField(S13,i,j,k) 
      lij(S22,i,j,k) = qTestField(S22,i,j,k)
      lij(S23,i,j,k) = qTestField(S23,i,j,k) 
      lij(S33,i,j,k) = qTestField(S33,i,j,k) 
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
! Load grid-based values 
!------------------------------------------------------------------------------

!    DO k=0, nz+1
!    DO j=0, ny+1
!    DO i=0, nx+1
!    DO iVar=1, nVarsConv
!          qField(iVar,i,j,k) = zero
!      qTestField(iVar,i,j,k) = zero
!    ENDDO ! iVar
!    ENDDO ! i
!    ENDDO ! j
!    ENDDO ! k
      
     DO k=0, nz+1
     DO j=0, ny+1
     DO i=0, nx+1
       q3Field(1,i,j,k) = u(i,j,k)
       q3Field(2,i,j,k) = v(i,j,k)
       q3Field(3,i,j,k) = w(i,j,k)
     ENDDO ! i
     ENDDO ! j
     ENDDO ! k

!------------------------------------------------------------------------------
! Perform test filtering on velocity field 
!------------------------------------------------------------------------------

     CALL TURB_ApplyTestFiltering( nVarsVel,    q3Field(:,0:,0:,0:), &
                                            q3TestField(:,0:,0:,0:)  )

!------------------------------------------------------------------------------
! Assemble L_{ij} 
!------------------------------------------------------------------------------

      DO k=1, nz-1
      DO j=1, ny-1
      DO i=1, nx-1
        lij(S11,i,j,k) = lij(S11,i,j,k) -q3TestField(1,i,j,k)*q3TestField(1,i,j,k)
        lij(S12,i,j,k) = lij(S12,i,j,k) -q3TestField(1,i,j,k)*q3TestField(2,i,j,k)
        lij(S13,i,j,k) = lij(S13,i,j,k) -q3TestField(1,i,j,k)*q3TestField(3,i,j,k)
        lij(S22,i,j,k) = lij(S22,i,j,k) -q3TestField(2,i,j,k)*q3TestField(2,i,j,k)
        lij(S23,i,j,k) = lij(S23,i,j,k) -q3TestField(2,i,j,k)*q3TestField(3,i,j,k)
        lij(S33,i,j,k) = lij(S33,i,j,k) -q3TestField(3,i,j,k)*q3TestField(3,i,j,k)     
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

!------------------------------------------------------------------------------
!  Monitor output
!------------------------------------------------------------------------------

    IF ( MOD(ntime,nmonitor) == 0 ) THEN
      PRINT*,'-----------------------------------------------------------------'
      WRITE(STDOUT, *)' Min-Max Vals of L11 = ',&
        MINVAL(lij(S11,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(lij(S11,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of L12 = ',&
        MINVAL(lij(S12,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(lij(S12,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of L13 = ',&
        MINVAL(lij(S13,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(lij(S13,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of L22 = ',&
        MINVAL(lij(S22,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(lij(S22,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of L23 = ',&
        MINVAL(lij(S23,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(lij(S23,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of L33 = ',&
        MINVAL(lij(S33,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(lij(S33,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of LII = ',&
      MINVAL(lij(S11,1:nx-1,1:ny-1,1:nz-1)&
            +lij(S22,1:nx-1,1:ny-1,1:nz-1)&
            +lij(S33,1:nx-1,1:ny-1,1:nz-1)),&
      MAXVAL(lij(S11,1:nx-1,1:ny-1,1:nz-1)&
            +lij(S22,1:nx-1,1:ny-1,1:nz-1)&
            +lij(S33,1:nx-1,1:ny-1,1:nz-1))
  
      PRINT*,'-----------------------------------------------------------------'
    ENDIF ! ntime

111 FORMAT(2(2X,1PE15.7))
!------------------------------------------------------------------------------
! Deallocate local dynamic arrays
!------------------------------------------------------------------------------

    DEALLOCATE(q3Field,STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'TURB_CalcLij: Memory Deallocation Error for q3Field'
       STOP
    ENDIF ! iErr

    DEALLOCATE(q3TestField,STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'TURB_CalcLij: Memory Deallocation Error for q3TestField'
      STOP
    ENDIF ! iErr

   END SUBROUTINE TURB_CalcLij
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_ContractLijMij()

!==============================================================================
!  Purpose: Contract Lij with Mij applying appropriate homogenization
!  Note   : Definition of L_{ij} does not need to be trace free since 
!           when contracting Lij with Mij, LijMij becomes trace-free 
!==============================================================================
      
    USE global_parameters
    USE turb_global_parameters
    USE turb_parameters
    USE flow_parameters
    USE flow_arrays
    USE turb_arrays
    USE boundary_arrays
    USE grid_arrays
    
    IMPLICIT NONE

!... Loop variables

    INTEGER :: i,j,k

!... Local variables

    REAL(KIND=CGREAL) :: mijMij, mijLij 

!------------------------------------------------------------------------------
! initialize smagorinsky Cs
!------------------------------------------------------------------------------

!    DO k=0, nz+1
!    DO j=0, ny+1
!    DO i=0, nx+1
!      viscTot(i,j,k) = zero
!    ENDDO ! i
!    ENDDO ! j
!    ENDDO ! k

!------------------------------------------------------------------------------
! first contract the numerator and denominator
!------------------------------------------------------------------------------

    DO k=1, nz-1
    DO j=1, ny-1
    DO i=1, nx-1

!------------------------------------------------------------------------------
!      M_{ij}*M_{ij}
!------------------------------------------------------------------------------

       mijMij =     mij(S11,i,j,k)*mij(S11,i,j,k) &
                   +mij(S22,i,j,k)*mij(S22,i,j,k) &
                   +mij(S33,i,j,k)*mij(S33,i,j,k) &
         +2._CGREAL*mij(S12,i,j,k)*mij(S12,i,j,k) &
         +2._CGREAL*mij(S13,i,j,k)*mij(S13,i,j,k) &
         +2._CGREAL*mij(S23,i,j,k)*mij(S23,i,j,k) 

!------------------------------------------------------------------------------
!      M_{ij}*L_{ij} 
!------------------------------------------------------------------------------

       mijLij =     mij(S11,i,j,k)*lij(S11,i,j,k) &
                   +mij(S22,i,j,k)*lij(S22,i,j,k) &
                   +mij(S33,i,j,k)*lij(S33,i,j,k) &
         +2._CGREAL*mij(S12,i,j,k)*lij(S12,i,j,k) &
         +2._CGREAL*mij(S13,i,j,k)*lij(S13,i,j,k) &
         +2._CGREAL*mij(S23,i,j,k)*lij(S23,i,j,k) 

!------------------------------------------------------------------------------
! Overwrite arrays to save storage 
!------------------------------------------------------------------------------

       mij(S11,i,j,k) = mijMij*REAL(1-iblank(i,j,k),KIND=CGREAL)
       lij(S11,i,j,k) = mijLij*REAL(1-iblank(i,j,k),KIND=CGREAL)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
!  Monitor output
!------------------------------------------------------------------------------

    IF ( MOD(ntime,nmonitor) == 0 ) THEN
      PRINT*,'-----------------------------------------------------------------'
      WRITE(STDOUT, *)' Min-Max Vals of LijMij = ',&
        MINVAL(lij(S11,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(lij(S11,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of MijMij = ',&
        MINVAL(mij(S11,1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(mij(S11,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *)' Min-Max Vals of CSmag = ',&
        MINVAL(viscTot(1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(viscTot(1:nx-1,1:ny-1,1:nz-1))
      PRINT*,'-----------------------------------------------------------------'
    ENDIF ! ntime

111 FORMAT(2(2X,1PE15.7))

   END SUBROUTINE TURB_ContractLijMij
!------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE TURB_CalcCsDynSmag()

!==============================================================================
!  Purpose: Compute Dynamic Smagorinsky Constant  
!           applying appropriate homogenization
!==============================================================================
      
    USE global_parameters
    USE turb_global_parameters
    USE turb_parameters
    USE flow_parameters
    USE flow_arrays
    USE turb_arrays
    USE boundary_arrays
    USE grid_arrays
    
    IMPLICIT NONE

!... Loop variables
      
    INTEGER :: i,j,k
            
!... Local variables
      
    INTEGER :: homogDirSum
    REAL(KIND=CGREAL) :: denomSum, eps, mijMij, mijLij, numSum

!******************************************************************************

!------------------------------------------------------------------------------       
! Invoke homogenization if needed
!------------------------------------------------------------------------------

    homogDirSum = SUM(homogDir)
 
!------------------------------------------------------------------------------
! set eps based on Tom Lund formalism 
!  this takes care of the case where MijMij is much smaller than LijMij
!  making the ratio substantially large.
!------------------------------------------------------------------------------

    IF ( homogDirSum == ZERO_DIR ) &
        eps = oned

    IF ( homogDir(DIRX) == ACTIVE ) &
      eps = (nx-1)
    IF ( homogDir(DIRY) == ACTIVE ) &
      eps = (ny-1)
    IF ( homogDir(DIRZ) == ACTIVE ) &
      eps = (nz-1)

    IF ( homogDir(DIRX) == ACTIVE .AND. homogDir(DIRY) == ACTIVE ) &
      eps = (nx-1)*(ny-1)
    IF ( homogDir(DIRX) == ACTIVE .AND. homogDir(DIRZ) == ACTIVE ) &
      eps = (nx-1)*(nz-1)
    IF ( homogDir(DIRY) == ACTIVE .AND. homogDir(DIRZ) == ACTIVE ) &
      eps = (ny-1)*(nz-1)

    IF (  homogDir(DIRX) == ACTIVE .AND. homogDir(DIRY) == ACTIVE .AND. &
          homogDir(DIRZ) == ACTIVE ) & 
      eps = (nx-1)*(ny-1)*(nz-1)

!------------------------------------------------------------------------------       
! Apply contraction operation based on  homogenization 
!  Note ratio is taken for (1-iblank)/=0, i.e. fluids cells
!------------------------------------------------------------------------------

    SELECT CASE (homogDirSum )

      CASE( ZERO_DIR )

!------------------------------------------------------------------------------
!       no-homogeneous direction
!------------------------------------------------------------------------------

        numSum   = zero
        denomSum = zero
        DO k=1, nz-1 
        DO j=1, ny-1
        DO i=1, nx-1
          numSum   = lij(S11,i,j,k)
          denomSum = mij(S11,i,j,k) +eps
     
          IF ( (1-iblank(i,j,k)) /= 0 ) THEN
            viscTot(i,j,k) = numSum/denomSum
          ENDIF
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      CASE( ONE_DIR )

!------------------------------------------------------------------------------
!       x-direction homogeneous only
!------------------------------------------------------------------------------

        IF ( homogDir(DIRX) == ACTIVE ) THEN
          DO k=1, nz-1 
          DO j=1, ny-1
            numSum   = zero
            denomSum = zero
            DO i=1, nx-1
              numSum   = numSum   +lij(S11,i,j,k)
              denomSum = denomSum +mij(S11,i,j,k) 
            ENDDO ! i 
            
            denomSum = denomSum +eps

            DO i=1, nx-1    
              IF ( (1-iblank(i,j,k)) /= 0 ) THEN
                viscTot(i,j,k) = numSum/denomSum
              ENDIF
            ENDDO ! i             
          ENDDO ! j
          ENDDO ! k
        ENDIF ! homogDir(DIRX)

!------------------------------------------------------------------------------
!       y-direction homogeneous only
!------------------------------------------------------------------------------
        
        IF ( homogDir(DIRY) == ACTIVE ) THEN
          DO k=1, nz-1 
          DO i=1, nx-1
            numSum   = zero
            denomSum = zero
            DO j=1, ny-1
              numSum   = numSum   +lij(S11,i,j,k)
              denomSum = denomSum +mij(S11,i,j,k) 
            ENDDO ! j 
            
            denomSum = denomSum +eps

            DO j=1, ny-1    
              IF ( (1-iblank(i,j,k)) /= 0 ) THEN
                viscTot(i,j,k) = numSum/denomSum
              ENDIF
            ENDDO ! j             
          ENDDO ! i
          ENDDO ! k
        ENDIF ! homogDir(DIRY)

!------------------------------------------------------------------------------
!       z-direction homogeneous only
!------------------------------------------------------------------------------

        IF ( homogDir(DIRZ) == ACTIVE ) THEN
          DO j=1, ny-1
          DO i=1, nx-1
            numSum   = zero
            denomSum = zero
            DO k=1, nz-1 
              numSum   = numSum   +lij(S11,i,j,k)
              denomSum = denomSum +mij(S11,i,j,k) 
            ENDDO ! k 
            
            denomSum = denomSum + eps

            DO k=1, nz-1    
              IF ( (1-iblank(i,j,k)) /= 0 ) THEN
                viscTot(i,j,k) = numSum/denomSum
              ENDIF
            ENDDO ! k             
          ENDDO ! i
          ENDDO ! j
        ENDIF ! homogDir(DIRZ)
                                
      CASE( TWO_DIR )

!------------------------------------------------------------------------------
!       xy-direction homogeneous 
!------------------------------------------------------------------------------

        IF ( homogDir(DIRX) == ACTIVE .AND. homogDir(DIRY) == ACTIVE ) THEN
          DO k=1, nz-1 
            numSum   = zero
            denomSum = zero
            DO j=1, ny-1
            DO i=1, nx-1
              numSum   = numSum   +lij(S11,i,j,k)
              denomSum = denomSum +mij(S11,i,j,k) 
            ENDDO ! i 
            ENDDO ! j
            
            denomSum = denomSum +eps

            DO j=1, ny-1
            DO i=1, nx-1    
              IF ( (1-iblank(i,j,k)) /= 0 ) THEN
                viscTot(i,j,k) = numSum/denomSum
              ENDIF
            ENDDO ! i             
            ENDDO ! j         
          ENDDO ! k
        ENDIF ! homogDir(DIRX)

!------------------------------------------------------------------------------
!       xz-direction homogeneous 
!------------------------------------------------------------------------------
        
        IF ( homogDir(DIRX) == ACTIVE .AND. homogDir(DIRZ) == ACTIVE ) THEN
          DO j=1, ny-1 
            numSum   = zero
            denomSum = zero
            DO k=1, nz-1
            DO i=1, nx-1
              numSum   = numSum   +lij(S11,i,j,k)
              denomSum = denomSum +mij(S11,i,j,k)
            ENDDO ! i 
            ENDDO ! k
            
            denomSum = denomSum +eps    

            DO k=1, nz-1
            DO i=1, nx-1    
              IF ( (1-iblank(i,j,k)) /= 0 ) THEN
                viscTot(i,j,k) = numSum/denomSum
              ENDIF
            ENDDO ! i             
            ENDDO ! k         
          ENDDO ! j
        ENDIF ! homogDir(DIRX)

!------------------------------------------------------------------------------
!       yz-direction homogeneous 
!------------------------------------------------------------------------------

        IF ( homogDir(DIRY) == ACTIVE .AND. homogDir(DIRZ) == ACTIVE ) THEN
          DO i=1, nx-1
            numSum   = zero
            denomSum = zero
            DO k=1, nz-1 
            DO j=1, ny-1
              numSum   = numSum   +lij(S11,i,j,k)
              denomSum = denomSum +mij(S11,i,j,k)
            ENDDO ! j 
            ENDDO ! k
            
            denomSum = denomSum +eps

            DO k=1, nz-1 
            DO j=1, ny-1    
              IF ( (1-iblank(i,j,k)) /= 0 ) THEN
                viscTot(i,j,k) = numSum/denomSum
              ENDIF
            ENDDO ! j             
            ENDDO ! k         
          ENDDO ! i
        ENDIF ! homogDir(DIRY)

!------------------------------------------------------------------------------
!     xyz-direction homogeneous
!------------------------------------------------------------------------------
                  
      CASE( THREE_DIR )
        numSum   = zero
        denomSum = zero
        DO k=1, nz-1 
        DO j=1, ny-1
        DO i=1, nx-1
          numSum   = numSum   + lij(S11,i,j,k)
          denomSum = denomSum + mij(S11,i,j,k)
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

        denomSum = denomSum +eps

        DO k=1, nz-1 
        DO j=1, ny-1
        DO i=1, nx-1
          IF ( (1-iblank(i,j,k)) /= 0 ) THEN
            viscTot(i,j,k) = numSum/denomSum
          ENDIF
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k
      
    END SELECT ! homogDirSum  

!------------------------------------------------------------------------------
!  Monitor output
!------------------------------------------------------------------------------

    IF ( MOD(ntime,nmonitor) == 0 ) THEN
      PRINT*,'-----------------------------------------------------------------'
      WRITE(STDOUT, *)' Min-Max Vals of CSmag-DynLes = ',&
        MINVAL(viscTot(1:nx-1,1:ny-1,1:nz-1)),&
        MAXVAL(viscTot(1:nx-1,1:ny-1,1:nz-1))
      PRINT*,'-----------------------------------------------------------------'
    ENDIF ! ntime

111 FORMAT(2(2X,1PE15.7))

999 CONTINUE
       
   END SUBROUTINE TURB_CalcCsDynSmag
!------------------------------------------------------------------------
