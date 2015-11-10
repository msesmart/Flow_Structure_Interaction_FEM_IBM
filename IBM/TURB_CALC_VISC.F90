!******************************************************************************
!
! Purpose: generalized kernel to compute the turbulence model based on LES
!
! Description: none.
!
! Input: u,v,w      = velocity Field
!
! Output: viscTurb = Turbulent viscosity
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
   SUBROUTINE TURB_print_inputs

!==============================================================================
!  Purpose: Set restart file and print input data
!==============================================================================
   
    USE global_parameters
    USE turb_global_parameters
    USE turb_parameters
    USE turb_arrays
    
    IMPLICIT NONE

!------------------------------------------------------------------------------
!   Print to output file information pertinent to turbulence model
!------------------------------------------------------------------------------
    
    WRITE(STDOUT,'(A,2X,I2)') &
      ' TurbModel (0: NoModel, 1: Dynamic Smagorinsky, 2: Dynamic Lagrangian)', &
        turbModel
    WRITE(STDOUT,'(A,2X,1PE12.5)') ' cSmagFix  = ',cSmagFix
    WRITE(STDOUT,'(A,2X,I2)')&
      ' filterWidthModel (1: Cube Root, 2: Geometric)  = ',filterWidthModel
    WRITE(STDOUT,'(A,3(2X,I2))')&
      ' testFilterDir    (0: None, 1: Active)          = ',testFilterDir(DIRX:DIRZ)
    WRITE(STDOUT,'(A,3(2X,I2))')&
      ' homogDir         (0: None, 1: Active)          = ',homogDIr(DIRX:DIRZ)
    WRITE(STDOUT,'(A,3(2X,I2))')&
      ' timeScale (1: Tom Lund Strain Magnitude, 2: JFM Formalism)  = ',turbLagrTimeScale

!------------------------------------------------------------------------------
!   Set restart file for Dynamic Lagrangian model
!------------------------------------------------------------------------------
    
    IF ( turbModel == TURB_MODEL_DYNLAGR ) THEN
       homogDir(1:3) = 0
       OPEN(ifuRstrtTurbIn,FILE='restart_turb_in.dat',FORM='UNFORMATTED') 
    ENDIF ! turbModel 

!------------------------------------------------------------------------------
!   Set global values
!------------------------------------------------------------------------------

    oneThird = oned/3.0_CGREAL
    twoThird = 2.0_CGREAL/3.0_CGREAL
    twoSqrt  = SQRT(2.0_CGREAL)

   END SUBROUTINE TURB_print_inputs
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_check_inputs

!==============================================================================
!  Purpose: Check input data
!==============================================================================
   
    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE turb_parameters
    USE turb_arrays
    
    IMPLICIT NONE
    INTEGER ::iDir
    REAL(KIND=CGREAL) :: sumTestFilter

    IF ( turbActive /= ACTIVE ) GOTO 999

    IF ( nDim == DIM_2D .AND. turbActive == ACTIVE ) THEN
      WRITE(STDOUT,'(A)')'Unable to run 2D Simulations with Turbulence Active'
      WRITE(STDOUT,'(A)')'Use turbActive: 0 with nDim: 2'
      WRITE(STDOUT,'(A)')' or turbActive: 1 with nDim: 2 or 3'
      STOP 
    ENDIF ! nDim_2D
                      
    IF ( turbModel < TURB_NOMODEL .OR. turbModel > TURB_MODEL_DYNLAGR ) THEN
      WRITE(STDOUT,'(A,2X,I2)')'Incorrect Value for turbModel in Turb Input',turbModel 
      WRITE(STDOUT,'(A)') 'Use turbModel: 0 (NoModel), 2 (Dynamic Smagorinsky)'
      WRITE(STDOUT,'(A)') '               2 (Dynamic Lagrangian) '
      STOP      
    ENDIF ! turbModel
 
    IF ( filterWidthModel < TURB_FWIDTH_CUBEROOT .OR. &
         filterWidthModel > TURB_FWIDTH_GEOM          ) THEN
      WRITE(STDOUT,'(A,2X,I2)') 'Incorrect Value for filterWidthModel in Turb Input',filterWidthModel 
      WRITE(STDOUT,'(A)') 'Use filterWidthModel: 1 (CubeRoot), 1 (Geometric)'
      STOP      
    ENDIF ! filterWidthModel

    DO iDir = 1, 3
      IF ( testFilterDir(iDir) < 0 .OR. testFilterDir(iDir) > 1 ) THEN
        WRITE(STDOUT,'(2(A,2X,I2))') &
         'Incorrect Value for testFilterDir in Turb Input',testFilterDir(iDir),'  in iDir = ',iDir
        WRITE(STDOUT,'(A)') 'testFilterDir is 0 (None), 1 (Active)'
        STOP      
      ENDIF ! homogDir

      IF ( homogDir(iDir) < 0 .OR. homogDir(iDir) > 1 ) THEN
        WRITE(STDOUT,'(2(A,2X,I2))') &
         'Incorrect Value for homogDir in Turb Input',homogDir(iDir),'  in iDir = ',iDir
        WRITE(STDOUT,'(A)') 'homogDir is 0 (Non-Homogeneous), or 1 (Homogeneous)'
        STOP      
      ENDIF ! homogDir
      
    ENDDO ! iDir
    
    sumTestFilter = SUM(testFilterDir)
    IF ( turbModel >= TURB_MODEL_DYNSMAG .AND. &
         turbModel <= TURB_MODEL_DYNLAGR .AND. &
         sumTestFilter == 0                    ) THEN
      WRITE(STDOUT,'(A,2X,I2)') 'Incorrect Value for testFilterDir in Turb Input',testFilterDir
      WRITE(STDOUT,'(A)') 'At Least one direction should have test filtering' 
      STOP
    ENDIF ! turbModel
    
    IF ( turbLagrTimeScale < TURB_LAGR_TSCALE_TSL .OR. &
         turbLagrTimeScale > TURB_LAGR_TSCALE_JFM      ) THEN
      WRITE(STDOUT,'(A,2X,I2)')'Incorrect Value for turbLagrTimeScale in Turb Input',turbLagrTimeScale 
      WRITE(STDOUT,'(A)') 'Use turbLagrTimeScale: 1 (Tom Lunds Magnitude Time Scalee)'
      WRITE(STDOUT,'(A)') '                       2 (JFM Time Scale) '
      STOP      
    ENDIF ! turbLagrTimeScale         

999 CONTINUE

   END SUBROUTINE TURB_check_inputs
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_allocate_memory()

!==============================================================================
!  Purpose: Allocate and initialize arrays pertinent to turbulence model
!==============================================================================

    USE global_parameters
    USE turb_global_parameters
    USE turb_parameters
    USE flow_parameters
    USE turb_arrays
    
    IMPLICIT NONE
    INTEGER :: iErr,i,j,k,iVar,nVarSij

!------------------------------------------------------------------------------
! Data pertinent to LES models
!------------------------------------------------------------------------------

    nVarSij = TENSOR_SYM_NELM
      
    PRINT*,'Allocate Lij-Mij-Sij Arrays'

    ALLOCATE(lij(nVarSij,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'TURB_allocate_memory: Memory Allocation Error for Lij'
      STOP
    ENDIF ! ierr 

    ALLOCATE(mij(nVarSij,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'TURB_allocate_memory: Memory Allocation Error for Mij'
      STOP
    ENDIF ! iErr 

    ALLOCATE(sij(nVarSij,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'TURB_allocate_memory: Memory Allocation Error for sij'
      STOP
    ENDIF ! iErr 

    ALLOCATE(qField(nVarSij,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'TURB_allocate_memory: Memory Allocation Error for qField'
      STOP
    ENDIF ! iErr 

    ALLOCATE(qTestField(nVarSij,0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'TURB_allocate_memory: Memory Allocation Error for qTestField'
      STOP
    ENDIF ! iErr 

!------------------------------------------------------------------------------
! Data pertinent to dynamic Lagrangian LES model
!------------------------------------------------------------------------------
           
    IF ( turbModel == TURB_MODEL_DYNLAGR ) THEN
      ALLOCATE(phiLM(0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'TURB_allocate_memory: Memory Allocation Error for PhiLM'
        STOP
      ENDIF ! iErr 

      ALLOCATE(phiMM(0:nx+1,0:ny+1,0:nz+1),STAT=ierr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'TURB_allocate_memory: Memory Allocation Error for PhiMM'
        STOP
      ENDIF ! iErr 

    ENDIF ! turbModel 

!------------------------------------------------------------------------------
! Initialize arrays
!------------------------------------------------------------------------------

    DO k = 0,nz+1
    DO j = 0,ny+1
    DO i = 0,nx+1
      DO iVar = 1, nVarSij
        lij(iVar,i,j,k) = zero
        mij(iVar,i,j,k) = zero      
        sij(iVar,i,j,k) = zero
     
       qField(iVar,i,j,k)      = zero
       qTestField(iVar,i,j,k)  = zero
     END DO ! iVar
     
     IF ( turbModel == TURB_MODEL_DYNLAGR ) THEN
      phiLM(i,j,k) = zero
      phiMM(i,j,k) = zero
     ENDIF ! turbModel 
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
  
   END SUBROUTINE TURB_allocate_memory
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_CalcVisc()

!==============================================================================
!  Purpose: Main driver to compute turbulent viscosity based on LES models
!==============================================================================

    USE turb_global_parameters
    USE flow_parameters
    USE flow_arrays
    USE turb_parameters
    USE turb_arrays
    USE boundary_arrays

    IMPLICIT NONE

!... Local variables
    
    INTEGER :: i,j,k
    INTEGER, DIMENSION(3) :: indexNode
    INTEGER  :: iFlagFilter
   
    REAL(KIND=CGREAL)     :: deltaSqr, strainMagnitude, turbViscLocal

!------------------------------------------------------------------------------
!   Set filter flag
!------------------------------------------------------------------------------

    iFlagFilter = ACTIVE-1

!------------------------------------------------------------------------------
!   Initialize arrays 
!------------------------------------------------------------------------------

    DO k = 0,nz+1
    DO j = 0,ny+1
    DO i = 0,nx+1
      viscTot(i,j,k) = zero
      bcxvisc(i,j,k) = zero
      bcyvisc(i,j,k) = zero
      bczvisc(i,j,k) = zero
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------    
!   Calculate Lij and Mij terms 
!------------------------------------------------------------------------------

    IF(MOD(ntime,nmonitor)==0) WRITE(STDOUT,*) '   Entering TURB_CalcMij '     
    CALL TURB_CalcMij()
    
    IF(MOD(ntime,nmonitor)==0) WRITE(STDOUT,*) '   Entering TURB_CalcLij ' 
    CALL TURB_CalcLij()

!------------------------------------------------------------------------------
!   Contract terms and apply appropriate averaging
!------------------------------------------------------------------------------

    IF(MOD(ntime,nmonitor)==0) WRITE(STDOUT,*) '   Entering TURB_ContractLijMij ' 
    CALL TURB_ContractLijMij()

!------------------------------------------------------------------------------
!   Compute Smagorinsky Constant 
!------------------------------------------------------------------------------
 
    SELECT CASE ( turbModel )
      CASE ( TURB_MODEL_DYNSMAG )
        CALL TURB_CalcCsDynSmag()

      CASE ( TURB_MODEL_DYNLAGR ) 
        CALL TURB_CalcCsDynLagr()

    END SELECT ! turbModel 

!------------------------------------------------------------------------------
!   Compute turbulent viscosity
!------------------------------------------------------------------------------

    DO k=1, nz-1
    DO j=1, ny-1
    DO i=1, nx-1

!------------------------------------------------------------------------------
!      Compute filter width
!------------------------------------------------------------------------------

      indexNode(DIRX) = i
      indexNode(DIRY) = j
      indexNode(DIRZ) = k

      CALL TURB_CalcFilterWidth(iFlagFilter,indexNode,deltaSqr)

!------------------------------------------------------------------------------
!     Evaluate strain magnitude 
!     |S(u)| = SQRT( 2 S_{ij}*S_{ij} )
!------------------------------------------------------------------------------

      strainMagnitude = twoSqrt * SQRT( sij(S11,i,j,k)*sij(S11,i,j,k) &
                                      + sij(S22,i,j,k)*sij(S22,i,j,k) &
                                      + sij(S33,i,j,k)*sij(S33,i,j,k) &    
                          + 2.0_CGREAL* sij(S12,i,j,k)*sij(S12,i,j,k) &
                          + 2.0_CGREAL* sij(S13,i,j,k)*sij(S13,i,j,k) &
                          + 2.0_CGREAL* sij(S23,i,j,k)*sij(S23,i,j,k) ) 

!------------------------------------------------------------------------------
!     Compute turbulent viscosity and clip off negative values
!------------------------------------------------------------------------------

      turbViscLocal   = deltaSqr *strainMagnitude *viscTot(i,j,k)

!!    viscTot(i,j,k) = DMAX1(turbViscLocal,-reInv) &
!!                    * ( oned -REAL(iblank(i,j,k),KIND=CGREAL) )

      viscTot(i,j,k) = half* (ABS(turbViscLocal)+turbViscLocal) &
                     * REAL(1-iblank(i,j,k),KIND=CGREAL)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------   
!  Monitor output
!------------------------------------------------------------------------------

    IF (MOD(ntime,nmonitor)==0) THEN
      WRITE(STDOUT, *) ' Min-Max Vals of viscTurb-Clipped = ',&
         MINVAL(viscTot(1:nx-1,1:ny-1,1:nz-1)),&
         MAXVAL(viscTot(1:nx-1,1:ny-1,1:nz-1))
    ENDIF ! ntime 

111 FORMAT(2(2X,1PE15.7))
    
   END SUBROUTINE TURB_CalcVisc
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_Visc_SetBoundCells()

!==============================================================================
!  Purpose: Set viscosity in fluids boundary cells near an immersed body
!==============================================================================

    USE turb_global_parameters
    USE flow_parameters
    USE flow_arrays
    USE turb_parameters
    USE turb_arrays
    USE boundary_arrays

    IMPLICIT NONE

!... Local variables
    
    INTEGER :: i,j,k,iblankCellBound

!------------------------------------------------------------------------------
! Set total viscosity in boundary cells in Fluid to molecular viscosity
! (ie. cells outside body which have 
! at least one neighbor in body)
!------------------------------------------------------------------------------

      DO k = 1, nz-1
      DO j = 1, ny-1
      DO i = 1, nx-1

!------------------------------------------------------------------------------
!       Mark boundary fluid cells near the immersed body
!------------------------------------------------------------------------------

        iblankCellBound = iblank(i-1,j  ,k)   &
                        + iblank(i+1,j  ,k)   &
                        + iblank(i  ,j-1,k)   &
                        + iblank(i  ,j+1,k)   &
                        + iblank(i  ,j  ,k-1) &
                        + iblank(i  ,j  ,k+1) 

        IF ( iblank(i,j,k) == 0 .AND. iblankCellBound > 0 ) THEN

          IF ( iblank(i+1,j,k) == 1 ) viscTot(i,j,k) = reinv
          IF ( iblank(i-1,j,k) == 1 ) viscTot(i,j,k) = reinv
          IF ( iblank(i,j+1,k) == 1 ) viscTot(i,j,k) = reinv
          IF ( iblank(i,j-1,k) == 1 ) viscTot(i,j,k) = reinv
          IF ( iblank(i,j,k+1) == 1 ) viscTot(i,j,k) = reinv
          IF ( iblank(i,j,k-1) == 1 ) viscTot(i,j,k) = reinv
        ENDIF ! iblank(i,j,k) == 0
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
          
   END SUBROUTINE TURB_Visc_SetBoundCells
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_Visc_SetFreshCell

!  Update Advection-Diffusion RHS for fresh cell

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE nlold_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k

    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1
      viscTot(i,j,k) = viscTot(i,j,k)*(oned - REAL(fresh_cell(i,j,k),KIND=CGREAL)) &
                     + reInv *REAL(fresh_cell(i,j,k),KIND=CGREAL)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

   END SUBROUTINE  TURB_Visc_SetFreshCell
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_Visc_set_bc()

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE turb_parameters
    USE boundary_arrays
    USE flow_arrays

    IMPLICIT NONE

!... loop variables

    INTEGER :: i,j,k

!******************************************************************************

! Set internal Boundaries Conditions

    DO k=1,nz-1    
    DO j=1,ny-1    
    DO i=1,nx-1    
      IF ( (1-ium(i,j,k))*(1-iup(i,j,k)) == 0 ) THEN
        bcxvisc(i,j,k) = zero
      ENDIF
      IF ( (1-jum(i,j,k))*(1-jup(i,j,k)) == 0 ) THEN
        bcyvisc(i,j,k) = zero
      ENDIF
      IF ( (1-kum(i,j,k))*(1-kup(i,j,k)) == 0 ) THEN
        bczvisc(i,j,k) = zero
      ENDIF
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k 
    
! Set outer boundary conditions
!  Note: Setting turbulent viscosity at cell center to zero for 
!        non-outflow boundary conditions

    DO k=0,nz
    DO j=0,ny

! left boundary ---------------------------------------------------------------
      i = 1
      SELECT CASE (bcx1)
        CASE (BC_TYPE_DIRICHLET)            ! dirichlet bc
          bcxvisc(i,j,k) = zero 
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_ZERO_GRADIENT)        ! outflow bc ( zero gradient ; explicit)
          bcxvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bcxvisc(i,j,k) = zero            
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SYMMETRY)             ! symmetry bc (explicit)
          bcxvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc  (explicit & dirty implementation)
          bcxvisc(i,j,k) = half*( viscTot(i,j,k) + viscTot(nx-1,j,k) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          bcxvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SHEAR)                ! shear bc 
          bcxvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
       END SELECT 

! right boundary---------------------------------------------------------------
      i = nx-1
      SELECT CASE (bcx2)
        CASE (BC_TYPE_DIRICHLET)            ! dirichlet bc
          bcxvisc(i,j,k) = zero 
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_ZERO_GRADIENT)        ! outflow bc ( zero gradient ; explicit)
          bcxvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bcxvisc(i,j,k) = zero            
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SYMMETRY)             ! symmetry bc (explicit)
          bcxvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc
          bcxvisc(i,j,k) = half*( viscTot(i,j,k) + viscTot(1,j,k) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          bcxvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SHEAR)                ! shear bc
          bcxvisc(i,j,k) = zero 
          viscTot(i,j,k) = zero 
      END SELECT 

    ENDDO ! j
    ENDDO ! k

    DO k=0,nz
    DO i=0,nx

! bottom boundary -------------------------------------------------------------
      j = 1
      SELECT CASE (bcy1)
        CASE (BC_TYPE_DIRICHLET)            ! dirichlet bc
          bcyvisc(i,j,k) = zero 
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_ZERO_GRADIENT)        ! outflow bc ( zero gradient ; explicit)
          bcyvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bcyvisc(i,j,k) = zero            
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SYMMETRY)             ! symmetry bc (explicit)
          bcyvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc 
          bcyvisc(i,j,k) = half*( viscTot(i,j,k) + viscTot(i,ny-1,k) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          bcyvisc(i,j,k) = zero 
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SHEAR)                ! shear bc
          bcyvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
      END SELECT 

! top boundary ----------------------------------------------------------------
      j = ny-1
      SELECT CASE (bcy2)
        CASE (BC_TYPE_DIRICHLET)            ! dirichlet bc
          bcyvisc(i,j,k) = zero 
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_ZERO_GRADIENT)        ! outflow bc ( zero gradient ; explicit)
          bcyvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bcyvisc(i,j,k) = zero            
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SYMMETRY)             ! symmetry bc (explicit)
          bcyvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc 
          bcyvisc(i,j,k) = half*( viscTot(i,j,k) + viscTot(i,1,k) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          bcyvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SHEAR)                ! shear bc
          bcyvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
      END SELECT 

    ENDDO ! i
    ENDDO ! k

    DO j=0,ny
    DO i=0,nx

! front boundary --------------------------------------------------------------
      k = 1
      SELECT CASE (bcz1)
        CASE (BC_TYPE_DIRICHLET)            ! diriclet bc
          bczvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_ZERO_GRADIENT)        ! outflow bc ( zero gradient ; explicit)
          bczvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bczvisc(i,j,k) = zero         
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SYMMETRY)             ! symmetry bc (explicit)
          bczvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc 
          bczvisc(i,j,k) = half*( viscTot(i,j,k) + viscTot(i,j,nz-1) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          bczvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SHEAR)          ! shear bc
          bczvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
      END SELECT 

! back boundary ---------------------------------------------------------------
      k = nz-1
      SELECT CASE (bcz2)
        CASE (BC_TYPE_DIRICHLET)            ! diriclet bc
          bczvisc(i,j,k) = zero 
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_ZERO_GRADIENT)        ! outflow bc ( zero gradient ; explicit)
          bczvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bczvisc(i,j,k) = zero            
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SYMMETRY)             ! symmetry bc (explicit)
          bczvisc(i,j,k) = viscTot(i,j,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc 
          bczvisc(i,j,k) = half*( viscTot(i,j,k) + viscTot(i,j,1) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          bczvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
        CASE (BC_TYPE_SHEAR)          ! shear bc
          bczvisc(i,j,k) = zero
          viscTot(i,j,k) = zero 
      END SELECT 

    ENDDO ! i
    ENDDO ! j

   END SUBROUTINE TURB_Visc_set_bc
!------------------------------------------------------------------------------

!******************************************************************************
!
! Purpose: Read, write restart files for LES turbulence model
!
!******************************************************************************
!
! $Id: Exp $
!
! Copyright: (c) 2004 by the George Washington University
!
!******************************************************************************
!------------------------------------------------------------------------------
   SUBROUTINE TURB_read_restart()

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE turb_parameters
    USE turb_arrays  
    
    IMPLICIT NONE
    
    INTEGER :: nTimeT, nxT, nyT, nzT

! Read restart file for dynamic lagrangian model 

    IF ( turbModel == TURB_MODEL_DYNLAGR ) THEN
      READ(ifuRstrtTurbIn) phiLM(0:nx+1,0:ny+1,0:nz+1), &
                           phiMM(0:nx+1,0:ny+1,0:nz+1)    
    END IF ! turbModel
    CLOSE(ifuRstrtFlowIn)
    
   END SUBROUTINE TURB_read_restart
!------------------------------------------------------------------------------

   SUBROUTINE TURB_write_restart()

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE turb_parameters
    USE turb_arrays  
    
    IMPLICIT NONE

! Write restart file for dynamic lagrangian model 

    IF ( turbModel == TURB_MODEL_DYNLAGR ) THEN
      WRITE(ifuRstrtTurbOut) phiLM(0:nx+1,0:ny+1,0:nz+1), &
                             phiMM(0:nx+1,0:ny+1,0:nz+1)    
    END IF ! turbModel

   END SUBROUTINE TURB_write_restart
!------------------------------------------------------------------------------
