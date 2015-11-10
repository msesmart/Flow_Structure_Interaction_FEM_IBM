!-------------------------------------------------------------------
!
!     Geometric MultiGrid Method PPE Solver (GMG-PPE)
!
!--------------------------------------------------------------------
! This files contains all subroutines for multi-grid method.
! It includes:
!  SUBROUTINE MG_initial()
!  SUBROUTINE MG_itsolv()
!  SUBROUTINE MG_residual()
!  SUBROUTINE MG_restrict()
!  SUBROUTINE MG_prolong()
!  SUBROUTINE MG_solver()
!  SUBROUTINE MG_solv_x, MG_solv_y, MG_solv_z
!  SUBROUTINE MG_prepare()
!---------------------------------------------------------------------
! Version 2.0 done by Z. Liang and Dr. H. Dong
!                           May 26th, 2008
!
! Copyright: (c) 2008 by Wright State University
!
!----------------------------------------------------------------------
!
! Original work done by Dr. H. Dong, instructed by Dr. R. Mittal
!
! Original Copyright: (c) 2006 by the George Washington University
!
!---------------------------------------------------------------------


SUBROUTINE MG_initial() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    INTEGER, PARAMETER :: MG_LEVEL_MAX=10
    INTEGER :: i, j, k, nLevel, ii, iim, jj, jjm, kk, kkm, nCoarseGrids


    nCoarseGrids = 4

    IF (mgLevels_X==0) THEN
      ii = nxc
      DO nLevel = 1, MG_LEVEL_MAX
        ii = ii/2 
        IF (ii>=nCoarseGrids) CYCLE
        EXIT 
      END DO ! nLevel
      mgLevels_X = nLevel
    ENDIF ! mgLevels_X

    IF (mgLevels_Y==0) THEN
      jj = nyc
      DO nLevel = 1, MG_LEVEL_MAX
        jj = jj/2
        IF (jj>=nCoarseGrids) CYCLE
        EXIT  
      END DO ! nLevel 
      mgLevels_Y = nLevel
    ENDIF ! mgLevels_Y

    IF (mgLevels_Z==0) THEN
      kk = NZ-1
      DO nLevel = 1, MG_LEVEL_MAX
        kk = kk/2
        IF (kk>=nCoarseGrids) CYCLE
        EXIT  
      END DO ! nLevel
      mgLevels_Z = nLevel
    ENDIF ! mgLevels_Z

    WRITE(STDOUT,'(5X,A,1X,I5)') 'mgLevels_X = ',mgLevels_X
    WRITE(STDOUT,'(5X,A,1X,I5)') 'mgLevels_Y = ',mgLevels_Y
    WRITE(STDOUT,'(5X,A,1X,I5)') 'mgLevels_Z = ',mgLevels_Z

    CALL MG_Memory_Allocation(mgLevels_X,ICOORD)

    IF (Full_Coarsening) THEN


    ELSE

	  CALL MG_Memory_Allocation(mgLevels_Y,JCOORD)
	  IF (ndim == DIM_3D) &
		CALL MG_Memory_Allocation(mgLevels_Z,KCOORD)
    ENDIF
      
	CALL MG_Allocate_Memory_Grid

    mgrid_I(1) = nxc
    mgrid_J(1) = nyc
    mgrid_K(1) = nzc

    DO nLevel = 2, mgLevels_X
      IF ( MOD(mgrid_I(nLevel-1),2) .EQ. 0 ) THEN
        mgrid_I(nLevel) = mgrid_I(nLevel-1)/2
      ELSE
        mgrid_I(nLevel) = (mgrid_I(nLevel-1)-1)/2
      ENDIF ! mgrid_I
    END DO ! nLevel

    DO nLevel = 2, mgLevels_Y
      IF ( MOD(mgrid_J(nLevel-1),2) .EQ. 0 ) THEN
        mgrid_J(nLevel) = mgrid_J(nLevel-1)/2
      ELSE
        mgrid_J(nLevel) = (mgrid_J(nLevel-1)-1)/2 
      ENDIF ! mgrid_J
    END DO ! n

    IF (ndim .EQ. DIM_3D) THEN
      DO nLevel = 2, mgLevels_Z
        IF ( MOD(mgrid_K(nLevel-1),2) .EQ. 0 ) THEN
          mgrid_K(nLevel) = mgrid_K(nLevel-1)/2
        ELSE
          mgrid_K(nLevel) = (mgrid_K(nLevel-1)-1)/2 
        ENDIF ! mgrid_K
      END DO ! nLevel
    ENDIF  ! ndim

    WRITE(STDOUT,'(5X,A)') 'mgrid_I '
    DO nLevel = 2, mgLevels_X
      WRITE(STDOUT,'(5X,A,2(3X,I5))')'nLevel, mgrid_I ',nLevel,mgrid_I(nLevel)
    END DO ! nLevel
    WRITE(STDOUT,*) '    '
    
    WRITE(STDOUT,'(5X,A)') 'mgrid_J '
    DO nLevel = 2, mgLevels_Y
      WRITE(STDOUT,'(5X,A,2(3X,I5))')'nLevel, mgrid_J ',nLevel,mgrid_J(nLevel)
    END DO ! nLevel
    WRITE(STDOUT,*) '    '

    IF (ndim .EQ. DIM_3D) THEN
      WRITE(STDOUT,'(5X,A)') 'mgrid_K '
      DO nLevel = 2, mgLevels_Z
        WRITE(STDOUT,'(5X,A,2(3X,I5))')'nLevel, mgrid_K ',nLevel,mgrid_K(nLevel)
      END DO ! nLevel
      WRITE(STDOUT,*) '    '
    ENDIF  ! ndim

    DO i = 0, nx
      dxc_MG(i,1) = dxc(i)
      dx_MG(i,1)  = dx(i)
      dxcinv_MG(i,1) = dxcinv(i)
      dxinv_MG(i,1)  = dxinv(i)
      xc_MG(i,1)     = xc(i)
    END DO ! i

    DO j = 0, ny 
      dyc_MG(j,1) = dyc(j)
      dy_MG(j,1)  = dy(j)
      dycinv_MG(j,1) = dycinv(j)
      dyinv_MG(j,1)  = dyinv(j)
      yc_MG(j,1)     = yc(j) 
    END DO  !j

    DO k = 0, nz 
      dzc_MG(k,1) = dzc(k)
      dz_MG(k,1)  = dz(k)
      dzcinv_MG(k,1) = dzcinv(k)
      dzinv_MG(k,1)  = dzinv(k)
      zc_MG(k,1)     = zc(k) 
    END DO 

    DO i = 0, nx+1
      x_MG(i,1) = x(i)
    END DO

    DO j = 0, ny+1
      y_MG(j,1) = y(j) 
    END DO

    DO k = 0, nz+1
      z_MG(k,1) = z(k) 
    END DO

    DO nLevel = 2, mgLevels_X

      DO i=1,mgrid_I(nLevel)+1
        iim = 2*i-1
        x_MG(i,nLevel) = x_MG(iim, nLevel-1)
      ENDDO ! i
      x_MG(0,nLevel)=-x_MG(2,nLevel)
      x_MG(mgrid_I(nLevel)+2,nLevel)=twod*x_MG(mgrid_I(nLevel)+1,nLevel)-x_MG(mgrid_I(nLevel),nLevel)

      DO i=0,mgrid_I(nLevel)+1
        dx_MG(i,nLevel) = x_MG(i+1,nLevel)-x_MG(i,nLevel)
        dxinv_MG(i,nLevel) = oned / dx_MG(i,nLevel)
      ENDDO ! i 

      DO i=0,mgrid_I(nLevel)+1
        xc_MG(i,nLevel) = half*(x_MG(i,nLevel)+x_MG(i+1,nLevel))
      ENDDO ! i

      DO i= 1,mgrid_I(nLevel)+1 
        dxc_MG(i,nLevel) = xc_MG(i,nLevel)-xc_MG(i-1,nLevel)
        dxcinv_MG(i,nLevel) = oned/dxc_MG(i,nLevel)
      ENDDO ! i

!      dxc_MG(1,nLevel)                 = dxc_MG(1,nLevel-1)*2.0_CGREAL
!      dxc_MG(mgrid_I(nLevel)+1,nLevel) = dxc_MG(mgrid_I(nLevel),nLevel)
!      dxcinv_MG(1,nLevel)                 = dxcinv_MG(1,nLevel-1)/2.0_CGREAL
!      dxcinv_MG(mgrid_I(nLevel)+1,nLevel) = dxcinv_MG(mgrid_I(nLevel),nLevel)

    END DO ! nLevel

    DO nLevel = 2, mgLevels_Y

      DO j=1,mgrid_J(nLevel)+1
        jjm = 2*j -1
        y_MG(j,nLevel) = y_MG(jjm, nLevel-1)
      ENDDO ! j
      y_MG(0,nLevel)=-y_MG(2,nLevel)
      y_MG(mgrid_J(nLevel)+2,nLevel)=twod*y_MG(mgrid_J(nLevel)+1,nLevel)-y_MG(mgrid_J(nLevel),nLevel)

      DO j=0,mgrid_J(nLevel)+1
        dy_MG(j,nLevel) = y_MG(j+1,nLevel)-y_MG(j,nLevel)
        dyinv_MG(j,nLevel) = oned/dy_MG(j,nLevel)
      ENDDO ! j

      DO j=0,mgrid_J(nLevel)+1
        yc_MG(j,nLevel) = half*(y_MG(j,nLevel)+y_MG(j+1,nLevel))
      ENDDO ! j

      DO j= 1,mgrid_J(nLevel)+1
        dyc_MG(j,nLevel) = yc_MG(j,nLevel)-yc_MG(j-1,nLevel)
        dycinv_MG(j,nLevel) = oned/dyc_MG(j,nLevel)
      ENDDO ! j

!      dyc_MG(1,nLevel)                 = dyc_MG(1,nLevel-1)*2.0_CGREAL
!      dyc_MG(mgrid_J(nLevel)+1,nLevel) = dyc_MG(mgrid_J(nLevel),nLevel)
!      dycinv_MG(1,nLevel)                 = dycinv_MG(1,nLevel-1)/2.0_CGREAL
!      dycinv_MG(mgrid_J(nLevel)+1,nLevel) = dycinv_MG(mgrid_J(nLevel),nLevel)

    END DO ! nLevel

    IF (ndim .EQ. DIM_3D) THEN
      DO nLevel = 2, mgLevels_Z

        DO k=1,mgrid_K(nLevel)+1
          kkm = 2*k-1
          z_MG(k,nLevel) = z_MG(kkm, nLevel-1)
        ENDDO ! k
		z_MG(0,nLevel)=-z_MG(2,nLevel)
		z_MG(mgrid_K(nLevel)+2,nLevel)=twod*z_MG(mgrid_K(nLevel)+1,nLevel)-z_MG(mgrid_K(nLevel),nLevel)

        DO k=1,mgrid_K(nLevel)
          dz_MG(k,nLevel) = z_MG(k+1,nLevel)-z_mg(k,nLevel)
          dzinv_MG(k,nLevel) = oned/dz_MG(k,nLevel)
        ENDDO

        DO k=1,mgrid_K(nLevel)
          zc_MG(k,nLevel) = half*(z_MG(k,nLevel)+z_MG(k+1,nLevel))
        ENDDO

        DO k= 1,mgrid_K(nLevel)+1
          dzc_MG(k,nLevel) = zc_MG(k,nLevel)-zc_MG(k-1,nLevel)
          dzcinv_MG(k,nLevel) = oned/dzc_MG(k,nLevel)
        ENDDO

!        dzc_MG(1,nLevel)                 = dzc_MG(1,nLevel-1)*2.0_CGREAL
!        dzc_MG(mgrid_K(nLevel)+1,nLevel) = dzc_MG(mgrid_K(nLevel),nLevel)
!        dzcinv_MG(1,nLevel)                 = dzcinv_MG(1,nLevel-1)/2.0_CGREAL
!        dzcinv_MG(mgrid_K(nLevel)+1,nLevel) = dzcinv_MG(mgrid_K(nLevel),nLevel)

      END DO ! nLevel
    ENDIF ! ndim

    WRITE(STDOUT,'(5X,A)') '             '

END SUBROUTINE MG_initial
!---------------------------------------------------------------------

SUBROUTINE MG_Solver(var,rr, MG_itSolver, MG_itResidual, compTime)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual
    REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,0:nz+1),INTENT (IN)    :: rr
    REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,0:nz+1),INTENT (INOUT) :: var
    REAL(KIND=CGREAL),INTENT (OUT) :: compTime

    INTEGER :: clock1, clock2, clock_rate
    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz) :: rr1
    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz) :: var1

!******************************************************************************
 
    CALL system_clock(clock1)
    
    rr1(0:nx,0:ny,0:nz)  = rr(0:nx,0:ny,0:nz)
    var1(0:nx,0:ny,0:nz) = var(0:nx,0:ny,0:nz)
    
    CALL MG_solv_x(var1, rr1, MG_itSolver, MG_itResidual)

    var(0:nx,0:ny,0:nz) = var1(0:nx,0:ny,0:nz)
!    CALL write_dump_debug('var ',0,var)

    CALL MG_solv_y(var1, rr1, MG_itSolver, MG_itResidual)

    var(0:nx,0:ny,0:nz) = var1(0:nx,0:ny,0:nz)
!    CALL write_dump_debug('var ',1,var)

    IF ( nDim  .EQ.  DIM_3D ) THEN
       CALL MG_solv_z(var1, rr1, MG_itSolver, MG_itResidual)
    ENDIF ! nDim
    
    var(0:nx,0:ny,0:nz) = var1(0:nx,0:ny,0:nz)

    CALL system_clock(clock2, clock_rate)
 
!    IF ( MOD(ntime,nmonitor) == 0 ) &
!      WRITE(STDOUT,*) '*****Total time for this cycle is:', &
!                       REAL(clock2-clock1)/REAL(clock_rate)
    compTime=REAL(clock2-clock1)/REAL(clock_rate)
!    WRITE(STDOUT,*) compTime

END SUBROUTINE MG_Solver
!------------------------------------------------------------------------------


SUBROUTINE MG_Solver_FMG(var,rr, MG_itSolver, MG_itResidual, compTime)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual
    REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,0:nz+1),INTENT (IN)    :: rr
    REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,0:nz+1),INTENT (INOUT) :: var
    REAL(KIND=CGREAL),INTENT (OUT) :: compTime

    INTEGER :: clock1, clock2, clock_rate
    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz) :: rr1
    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz) :: var1
    
    INTEGER  :: GMM     ! NEW ZX

!******************************************************************************
 
    CALL system_clock(clock1)
    
    rr1(0:nx,0:ny,0:nz)  = rr(0:nx,0:ny,0:nz)
    var1(0:nx,0:ny,0:nz) = var(0:nx,0:ny,0:nz)
    
    IF (mgCycleX == 1) THEN
              GMM = 1
    ELSE
              GMM = 2
    ENDIF
 
    CALL MG_FMG(var1, rr1, 1, mgLevels_X, nx, ny, nz, mgrid_I(2)+1, 1, GMM, mgCycleX, MG_itSolver, MG_itResidual)

    var(0:nx,0:ny,0:nz) = var1(0:nx,0:ny,0:nz)

    CALL system_clock(clock2, clock_rate)
 
!    IF ( MOD(ntime,nmonitor) == 0 ) &
!      WRITE(STDOUT,*) '*****Total time for this cycle is:', &
!                       REAL(clock2-clock1)/REAL(clock_rate)
    compTime=REAL(clock2-clock1)/REAL(clock_rate)
!    WRITE(STDOUT,*) compTime

END SUBROUTINE MG_Solver_FMG
!------------------------------------------------------------------------------


SUBROUTINE MG_Solv_x(var, rr, MG_itSolver, MG_itResidual)

! ---------------------------------------------------------------
! An Instroduction to Multigrid Methods by Pieter Wesseling
!----------------------------------------------------------------
!            V cycle      F cycle                    W cycle
!*                *            O                       @
! \              /            /                       / 
!  *            *            O                       @
!   \          /            / \                     /
!    *        *    O       O   @       @           @
!     \      /    / \     /     \     / \         /
!      *    *    O   O   O       @   @   @   @   @
!       \  / \  /     \ /         \ /     \ / \ /
!        *    O        O           @       @   @
! ---------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual

    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz),INTENT (INOUT) :: rr 
    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz),INTENT (INOUT) :: var

    INTEGER  :: GMM     ! NEW ZX

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
    CALL MG_Memory_Allocation_iblank(mgLevels_X,ICOORD)
    CALL MG_Prepare_IBlank(ICOORD)
    ENDIF

    IF (mgCycleX == 1) THEN
              GMM = 1
    ELSE
              GMM = 2
    ENDIF
 
    CALL MG_X(var, rr, 1, mgLevels_X, nx, ny, nz, mgrid_I(2)+1, 1, GMM, mgCycleX, MG_itSolver, MG_itResidual)  ! NEW ZX

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
    CALL MG_Memory_Deallocation_iblank(mgLevels_X,ICOORD)
    ENDIF

END SUBROUTINE MG_Solv_x
!---------------------------------------------------------------------



SUBROUTINE MG_Solv_y(var, rr, MG_itSolver, MG_itResidual)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual

    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz),INTENT (INOUT) :: rr 
    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz),INTENT (INOUT) :: var

    INTEGER  :: GMM     ! NEW ZX

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
    CALL MG_Memory_Allocation_iblank(mgLevels_Y,JCOORD)
    CALL MG_Prepare_IBLANK(JCOORD)
    ENDIF

    IF (mgCycleY == 1) THEN
              GMM = 1
    ELSE
              GMM = 2
    ENDIF
    
    CALL MG_Y(var, rr, 1, mgLevels_Y, nx, ny, nz, mgrid_J(2)+1, 1, GMM, mgCycleY, MG_itSolver, MG_itResidual)  ! NEW ZX

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
    CALL MG_Memory_Deallocation_iblank(mgLevels_Y,JCOORD)
    ENDIF

END SUBROUTINE MG_Solv_y
!---------------------------------------------------------------------


SUBROUTINE MG_Solv_z(var, rr, MG_itSolver, MG_itResidual)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual

    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz),INTENT (INOUT) :: rr 
    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz),INTENT (INOUT) :: var

!   Local variables
    INTEGER  :: GMM     ! NEW ZX

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
    CALL MG_Memory_Allocation_iblank(mgLevels_Z,KCOORD)
    CALL MG_Prepare_IBLANK(KCOORD)
    ENDIF

    IF (mgCycleX == 1) THEN
              GMM = 1
    ELSE
              GMM = 2
    ENDIF

    CALL MG_Z(var, rr, 1, mgLevels_Z, nx, ny, nz, mgrid_K(2)+1, 1, GMM, mgCycleZ, MG_itSolver, MG_itResidual)  ! NEW ZX

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
    CALL MG_Memory_Deallocation_iblank(mgLevels_Z,KCOORD)
    ENDIF

END SUBROUTINE MG_Solv_z
!---------------------------------------------------------------------

SUBROUTINE MG_Memory_Allocation(nlev, NC_dim)

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER,INTENT (IN) :: nlev, NC_dim
    INTEGER :: ierr, i

!   Allocate MG arrays
!   ------------------
    SELECT CASE (NC_dim)
    CASE(ICOORD)
        ALLOCATE( MGX(nlev),STAT=ierr )
    CASE(JCOORD)
        ALLOCATE( MGY(2:nlev),STAT=ierr )
    CASE(KCOORD)
        ALLOCATE( MGZ(2:nlev),STAT=ierr )
    END SELECT

    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Memory_Allocation: Memory Allocation Error for MG', NC_dim
      STOP
    ENDIF ! ierr

END SUBROUTINE MG_Memory_Allocation
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE MG_Memory_Deallocation(nlev, NC_dim)

    USE MG_parameters
    USE MG_arrays

    INTEGER,INTENT (IN) :: nlev, NC_dim

    INTEGER :: ierr, i

    SELECT CASE (NC_dim)
    CASE(ICOORD)
        DEALLOCATE(MGX,STAT=ierr)
    CASE(JCOORD)
        DEALLOCATE(MGY,STAT=ierr)
    CASE(KCOORD)
        DEALLOCATE(MGZ,STAT=ierr)
    END SELECT
        
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Memory_DeAllocation: Memory DeAllocation Error for MG', NC_dim
      STOP
    ENDIF ! ierr 
    
END SUBROUTINE MG_Memory_deallocation
!---------------------------------------------------------------------


SUBROUTINE MG_Memory_Allocation_iblank(nlev, NC_dim)

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER,INTENT (IN) :: nlev, NC_dim
    INTEGER :: ierr, i
    
!   Allocate MG arrays
!   ------------------
    IF (Full_Coarsening) THEN
    
    do i=2, nlev
        ALLOCATE( MGX(i)%iblank(0:mgrid_I(i)+1,0:mgrid_I(i)+1,0:nz), STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for iblank_MGX'
          STOP
        ENDIF ! ierr
        MGX(i)%iblank = 0
    end do
    
    ELSE
    
    SELECT CASE (NC_dim)
    CASE(ICOORD)
    do i=2, nlev
        ALLOCATE( MGX(i)%iblank(0:mgrid_I(i)+1,0:ny,0:nz), STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for iblank_MGX'
          STOP
        ENDIF ! ierr
        MGX(i)%iblank = 0
    end do
    CASE(JCOORD)
    do i=2, nlev
        ALLOCATE( MGY(i)%iblank(0:nx,0:mgrid_J(i)+1,0:nz), STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for iblank_MGY'
          STOP
        ENDIF ! ierr
        MGY(i)%iblank = 0
    end do
    CASE(KCOORD)
    do i=2, nlev
        ALLOCATE( MGZ(i)%iblank(0:nx,0:ny,0:mgrid_K(i)+1), STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for iblank_MGZ'
          STOP
        ENDIF ! ierr
        MGZ(i)%iblank = 0
    end do
    END SELECT
    ENDIF
    
END SUBROUTINE MG_Memory_Allocation_iblank
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE MG_Memory_Deallocation_iblank(nlev, NC_dim)

    USE MG_parameters
    USE MG_arrays

    INTEGER,INTENT (IN) :: nlev, NC_dim

    INTEGER :: ierr, i

    SELECT CASE (NC_dim)
    CASE(ICOORD)
        do i=2, nlev
            DEALLOCATE(MGX(i)%iblank,STAT=ierr)
            IF ( iErr /= ERR_NONE ) THEN
              WRITE(STDOUT,*) &
              'MG_Memory_DeAllocation: Memory DeAllocation Error for MGX%iblank', I
              STOP
            ENDIF ! ierr 
        end do
    CASE(JCOORD)
        do i=2, nlev
            DEALLOCATE(MGY(i)%iblank,STAT=ierr)
            IF ( iErr /= ERR_NONE ) THEN
              WRITE(STDOUT,*) &
              'MG_Memory_DeAllocation: Memory DeAllocation Error for MGY%iblank', I
              STOP
            ENDIF ! ierr 
        end do
    CASE(KCOORD)
        do i=2, nlev
            DEALLOCATE(MGZ(i)%iblank,STAT=ierr)
            IF ( iErr /= ERR_NONE ) THEN
              WRITE(STDOUT,*) &
              'MG_Memory_DeAllocation: Memory DeAllocation Error for MGZ%iblank', I
              STOP
            ENDIF ! ierr 
        end do
    END SELECT
        
END SUBROUTINE MG_Memory_Deallocation_iblank
!---------------------------------------------------------------------

SUBROUTINE MG_Prepare_Iblank(NC_dim)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    
    INTEGER, INTENT(IN)   :: NC_dim

    INTEGER           :: i, j, k, n, ii, iim, jj, jjm, kk, kkm
    REAL(KIND=CGREAL) :: ttemp, area1, area2, area


    IF ( internal_boundary_present == 1 ) THEN

!      DO k = 1, nzc
!      DO j = 1, nyc
!      DO i = 1, nxc
!        iblk_MG(1)%iblank_MG(i,j,k) = iblank(i,j,k)
!      END DO ! i
!      END DO ! j
!      END DO ! k      

      SELECT CASE (NC_dim)

      CASE(ICOORD)
        DO n = 2, mgLevels_X  
          DO k = 1, nzc
          DO j = 1, nyc
          DO i = 1, mgrid_I(n)
            ii  = 2*i 
            iim = 2*i -1

            IF (N>2) THEN
            MGX(n)%iblank(i,j,k) = MGX(n-1)%iblank(iim,j,k) &
                                 * MGX(n-1)%iblank(ii,j,k) 
            ELSE
            MGX(n)%iblank(i,j,k) = iblank(iim,j,k) * iblank(ii,j,k) 
            ENDIF
          END DO ! i
          END DO ! j
          END DO ! k
        END DO ! n

      CASE(JCOORD)
        DO n = 2, mgLevels_Y  
          DO k = 1, nzc
          DO i = 1, nxc
          DO j = 1, mgrid_J(n)
            jj  = 2*j
            jjm = 2*j -1

            IF (N>2) THEN
            MGY(n)%iblank(i,j,k) = MGY(n-1)%iblank(i,jjm,k) &
                                 * MGY(n-1)%iblank(i,jj,k)
            ELSE
            MGY(n)%iblank(i,j,k) = iblank(i,jjm,k) * iblank(i,jj,k)
            ENDIF
          END DO ! j
          END DO ! i
          END DO ! k
        END DO ! n

      CASE(KCOORD)
        DO n = 2, mgLevels_Z
          DO j = 1, nyc
          DO i = 1, nxc
          DO k = 1, mgrid_K(n)
            kk  = 2*k
            kkm = 2*k -1

            IF (N>2) THEN
            MGZ(n)%iblank(i,j,k) = MGZ(n-1)%iblank(i,j,kkm) &
                                 * MGZ(n-1)%iblank(i,j,kk)
            ELSE
            MGZ(n)%iblank(i,j,k) = iblank(i,j,kkm) * iblank(i,j,kk)
            ENDIF
          END DO ! k
          END DO ! j
          END DO ! i
        END DO ! n
      END SELECT ! NC_dim 

    ELSE  ! internal_boundary
      
      SELECT CASE (NC_dim)

      CASE(ICOORD)

        DO n = 2, mgLevels_X
           MGX(n)%iblank(1:mgrid_I(n),1:nyc,1:nzc) = 0
        ghostcellMark_MG(1:mgrid_I(n),1:nyc,1:nzc) = 0        !!! Added by ZX Liang
        END DO ! n
          


      CASE(JCOORD)

        DO n = 2, mgLevels_Y
           MGY(n)%iblank(1:nxc,1:mgrid_J(n),1:nzc) = 0
        ghostcellMark_MG(1:nxc,1:mgrid_J(n),1:nzc) = 0        !!! Added by ZX Liang
        END DO ! n
          

      CASE(KCOORD)

        DO n = 2, mgLevels_Z
           MGZ(n)%iblank(1:nxc,1:nyc,1:mgrid_K(n)) = 0
        ghostcellMark_MG(1:nxc,1:nyc,1:mgrid_K(n)) = 0        !!! Added by ZX Liang
        END DO ! n
          

      END SELECT ! NC_dim

    ENDIF ! internal_boundary_present

END SUBROUTINE MG_Prepare_Iblank

SUBROUTINE MG_Prepare_Iblank_FMG

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER           :: i, j, k, n, ii, iim, jj, jjm, kk, kkm
    REAL(KIND=CGREAL) :: ttemp, area1, area2, area


    IF ( internal_boundary_present == 1 ) THEN

        DO n = 2, mgLevels_X  
          DO k = 1, nzc
          DO j = 1, mgrid_I(n)
            jj  = 2*j
            jjm = 2*j - 1
          DO i = 1, mgrid_I(n)
            ii  = 2*i 
            iim = 2*i - 1

            IF (N>2) THEN
            MGX(n)%iblank(i,j,k) = MGX(n-1)%iblank(iim,jjm,k) &
                                 * MGX(n-1)%iblank(ii ,jjm,k) &
                                 * MGX(n-1)%iblank(iim,jj,k) &
                                 * MGX(n-1)%iblank(ii ,jj,k) 
            ELSE
            MGX(n)%iblank(i,j,k) = iblank(iim,jjm,k) * iblank(ii,jjm,k) * iblank(iim,jj,k) * iblank(ii,jj,k) 
            ENDIF
          END DO ! i
          END DO ! j
          END DO ! k
        END DO ! n

    ELSE  ! internal_boundary
      

        DO n = 2, mgLevels_X
           MGX(n)%iblank(1:mgrid_I(n),1:mgrid_I(n),1:nzc) = 0
        ghostcellMark_MG(1:mgrid_I(n),1:mgrid_I(n),1:nzc) = 0
        END DO ! n
          

    ENDIF ! internal_boundary_present

END SUBROUTINE MG_Prepare_Iblank_FMG
!---------------------------------------------------------------------

!!---------------------------------------------------------------------
!
!RECURSIVE SUBROUTINE MG(var, rhs, n, nlev, mx, my, mz, NC_dim)
!
!    USE global_parameters
!    USE flow_parameters
!    USE grid_arrays
!    USE boundary_arrays
!    USE multiuse_arrays
!    USE GCM_arrays
!    USE MG_parameters
!    USE MG_arrays
! 
!    IMPLICIT NONE
!
!    INTEGER, INTENT(IN) :: n, nlev, mx, my, mz, NC_dim  
!    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (IN)    :: rhs
!    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: var
!    
!    INTEGER  :: i, j, k, nIter, totIter, iErr
!    REAL(KIND=CGREAL),DIMENSION(:,:,:), ALLOCATABLE   :: r, r1, phi
!
!!******************************************************************************
!
!!------------------------------------------------------------------------------
!!   Initialize variable
!!------------------------------------------------------------------------------
!
!    iblank_MG => iblk_MG(n)%iblank_MG
!    
!!------------------------------------------------------------------------------
!!   Allocate local arrays
!!------------------------------------------------------------------------------
!
!    IF (n /= nlev ) THEN
!        ALLOCATE(r(0:mx,0:my,0:mz),STAT=iErr)
!        IF ( iErr /= ERR_NONE ) THEN
!          WRITE(STDOUT,*) &
!          'MG_X: Memory Allocation Error for r'
!          STOP
!        ENDIF ! iErr
!        
!        SELECT CASE (NC_dim)
!        CASE (ICOORD)
!            ALLOCATE(r1(0:mgrid_I(n+1)+1,0:my,0:mz),STAT=iErr)
!            IF ( iErr /= ERR_NONE ) THEN
!              WRITE(STDOUT,*) &
!              'MG_X: Memory Allocation Error for r1'
!              STOP
!            ENDIF ! iErr
!            
!            ALLOCATE(phi(0:mgrid_I(n+1)+1,0:my,0:mz),STAT=iErr)
!            IF ( iErr /= ERR_NONE ) THEN
!              WRITE(STDOUT,*) &
!              'MG_X: Memory Allocation Error for phi'
!              STOP
!            ENDIF ! iErr
!            
!        CASE (JCOORD)
!            ALLOCATE(r1(0:mx,0:mgrid_J(n+1)+1,0:mz),STAT=iErr)
!            IF ( iErr /= ERR_NONE ) THEN
!              WRITE(STDOUT,*) &
!              'MG_X: Memory Allocation Error for r1'
!              STOP
!            ENDIF ! iErr
!            
!            ALLOCATE(phi(0:mx,0:mgrid_J(n+1)+1,0:mz),STAT=iErr)
!            IF ( iErr /= ERR_NONE ) THEN
!              WRITE(STDOUT,*) &
!              'MG_X: Memory Allocation Error for phi'
!              STOP
!            ENDIF ! iErr
!        
!        CASE (KCOORD)
!
!            ALLOCATE(r1(0:mx,0:my,0:mgrid_K(n+1)+1),STAT=iErr)
!            IF ( iErr /= ERR_NONE ) THEN
!              WRITE(STDOUT,*) &
!              'MG_X: Memory Allocation Error for r1'
!              STOP
!            ENDIF ! iErr
!            
!            ALLOCATE(phi(0:mx,0:my,0:mgrid_K(n+1)+1),STAT=iErr)
!            IF ( iErr /= ERR_NONE ) THEN
!              WRITE(STDOUT,*) &
!              'MG_X: Memory Allocation Error for phi'
!              STOP
!            ENDIF ! iErr
!        END SELECT
!    ENDIF
!    
!    IF (n == 1) THEN
!        totIter = iterFinest
!
!    ELSE ! n > 1
!        
!        IF ((n == mgLevels_X .and. NC_dim == ICOORD) .or. &
!            (n == mgLevels_Y .and. NC_dim == JCOORD) .or. &
!            (n == mgLevels_Z .and. NC_dim == KCOORD)) THEN
!          totIter = iterCoarsest
!        ELSE
!          totIter = 1
!        ENDIF ! n==mgLevels_X
!
!        var = zero
!        
!    ENDIF ! n==1
!
!    CALL MG_Prepare_BC(n)
!
!    IF ( n /= 1 ) THEN
!        CALL MG_Prepare( n, mx-1, my-1, mz-1)
!    ENDIF ! n
!
!    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) THEN
!        CALL set_outer_ghost_pres(var)
!        CALL GCM_p_set_bc_internal(var)
!        CALL GCM_enforce_p_compatibility(var)
!
!        CALL MG_Prepare_BC(n)
!    ENDIF ! boundary_formulation 
!
!    DO nIter = 1, totIter
!        SELECT CASE (NC_dim)
!        CASE (ICOORD)
!            CALL MG_itsolv(var, rhs, n, 1, 1, mx, my, mz)
!        CASE (JCOORD)
!            CALL MG_itsolv(var, rhs, 1, n, 1, mx, my, mz) 
!        CASE (KCOORD)
!            CALL MG_itsolv(var, rhs, 1, 1, n, mx, my, mz) 
!        END SELECT
!    END DO ! nIter   
!
!    IF ( n == nlev ) return
!        
!    SELECT CASE (NC_dim)
!    CASE (ICOORD)
!        CALL MG_Residual(var, rhs, r, n, 1, 1, mx, my, mz)
!        CALL MG_Restrict_X(r, r1, n, mx, my, mz, mgrid_I(n+1)+1)
!        CALL MG(phi, r1, n+1, nlev, mgrid_I(n+1)+1, my, mz, ICOORD)
!        iblank_MG => iblk_MG(n+1)%iblank_MG
!        CALL MG_Prolong_X(r, phi, n, mgrid_I(n+1)+1, my, mz, mx)
!    CASE (JCOORD)
!        CALL MG_Residual(var, rhs, r, 1, n, 1, mx, my, mz)
!        CALL MG_Restrict_Y(r, r1, n, mx, my, mz, mgrid_J(n+1)+1)
!        CALL MG(phi, r1, n+1, nlev, mx, mgrid_J(n+1)+1, mz, JCOORD)
!        iblank_MG => iblk_MG(n+1)%iblank_MG
!        CALL MG_Prolong_Y(r, phi, n, mx, mgrid_J(n+1)+1, mz, my)
!    CASE (KCOORD)
!        CALL MG_Residual(var, rhs, r, 1, 1, n, mx, my, mz)
!        CALL MG_Restrict_Z(r, r1, n, mx, my, mz, mgrid_K(n+1)+1)
!        CALL MG(phi, r1, n+1, nlev, mx, my, mgrid_K(n+1)+1, KCOORD)
!        iblank_MG => iblk_MG(n+1)%iblank_MG
!        CALL MG_Prolong_Z(r, phi, n, mx, my, mgrid_K(n+1)+1, mz)
!    END SELECT
!    
!    iblank_MG => iblk_MG(n)%iblank_MG
!    
!    var(1:mx-1,1:my-1,1:mz-1) = var(1:mx-1,1:my-1,1:mz-1) +r(1:mx-1,1:my-1,1:mz-1) &
!                                *(oned - REAL(iblank_MG(1:mx-1,1:my-1,1:mz-1),KIND=CGREAL))
!
!    CALL MG_Prepare_BC(n)
!
!    IF ( n /= 1 ) THEN
!        CALL MG_Prepare( n, mx-1, my-1, mz-1)
!    ENDIF ! n
!
!    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) THEN
!        CALL set_outer_ghost_pres(var)
!        CALL GCM_p_set_bc_internal(var)
!        CALL GCM_enforce_p_compatibility(var)
!
!        CALL MG_Prepare_BC(n)
!    ENDIF ! boundary_formulation 
!
!    IF (n == 1) THEN
!        IF ( ndim == 2 .and. NC_dim /= ICOORD) THEN
!            totIter = iterFinest
!        ELSE
!            totIter = 0
!        ENDIF ! ndim
!    ELSE
!        totIter = 1 
!    ENDIF ! n
!
!    DO nIter = 1, totIter
!        SELECT CASE (NC_dim)
!        CASE (ICOORD)
!            CALL MG_itsolv(var, rhs, n, 1, 1, mx, my, mz)
!        CASE (JCOORD)
!            CALL MG_itsolv(var, rhs, 1, n, 1, mx, my, mz) 
!        CASE (KCOORD)
!            CALL MG_itsolv(var, rhs, 1, 1, n, mx, my, mz) 
!        END SELECT
!    END DO ! nIter   
!       
!    IF ( n == 1 ) THEN
!
!        IF (boundary_formulation == GCM_METHOD .and. (NC_dim == JCOORD .or. NC_dim == KCOORD) ) THEN
!             CALL set_outer_ghost_pres(var)
!             CALL GCM_p_set_bc_internal(var)
!             CALL GCM_enforce_p_compatibility(var)
!             CALL MG_Prepare_BC(n)
!        ENDIF ! boundary_formulation
!
!    ELSE
!    ENDIF ! n==1
!
!    IF (infoconv == 1) THEN
!        CALL MG_Residual(var, rhs, r, n, 1, 1, mx-1, my-1, mz-1)
!    ENDIF ! infoConv
!
!    DEALLOCATE(r,STAT=iErr)
!    IF ( iErr /= ERR_NONE ) THEN
!          WRITE(STDOUT,*) &
!          'MG: Memory Deallocation Error for r'
!          STOP
!    ENDIF ! iErr
!    
!    DEALLOCATE(r1,STAT=iErr)
!    IF ( iErr /= ERR_NONE ) THEN
!          WRITE(STDOUT,*) &
!          'MG: Memory Deallocation Error for r1'
!          STOP
!    ENDIF ! iErr
!    
!    DEALLOCATE(phi,STAT=iErr)
!    IF ( iErr /= ERR_NONE ) THEN
!          WRITE(STDOUT,*) &
!          'MG: Memory Deallocation Error for phi'
!          STOP
!    ENDIF ! iErr
!
!END SUBROUTINE MG

RECURSIVE SUBROUTINE MG_FMG(var, rhs, n, nlev, mx, my, mz, mx1, WF, GMM, mgCycle, MG_itSolver, MG_itResidual)
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual

    INTEGER, INTENT(IN) :: n, nlev, mx, my, mz, mx1, WF, mgCycle   ! NEW ZX 
    INTEGER, INTENT(INOUT) :: GMM
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: rhs
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: var
    
    INTEGER  :: i, j, k, nIter, totIter, iErr
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz) :: r
    REAL(KIND=CGREAL),DIMENSION(0:mx1,0:mx1,0:mz) :: r1, phi

!------------------------------------------------------------------------------
!   Initialize variable
!------------------------------------------------------------------------------

    IF (N > 1) THEN
      iblank_MG => MGX(n)%iblank
    ELSE
      iblank_MG => iblank
    ENDIF

    IF (n == 1) THEN
      totIter = iterFinest

    ELSE ! n > 1
        
      IF (n == mgLevels_X) THEN
        totIter = iterCoarsest
      ELSE
        totIter = 1
      ENDIF ! n==mgLevels_X

    ENDIF ! n==1

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
      CALL MG_Prepare_BC(n)

      IF ( n /= 1 ) THEN
        CALL MG_Prepare( n, mx-1, my-1, mz-1)
      ENDIF ! n
    ENDIF

    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

!	CALL write_dump_debug_MG('Bvar',10+n,var,n, 1, 1, 0,mx,0,my,0,mz)
!	CALL write_dump_debug_i_MG('ium ',10+n,ium_mg, n, 1, 1, 0,nx+1,0,ny+1,0,nz+1)
!	CALL write_dump_debug_i_MG('iup ',10+n,iup_mg, n, 1, 1, 0,nx+1,0,ny+1,0,nz+1)
!	CALL write_dump_debug_i_MG('jum ',10+n,jum_mg, n, 1, 1, 0,nx+1,0,ny+1,0,nz+1)
!	CALL write_dump_debug_i_MG('jup ',10+n,jup_mg, n, 1, 1, 0,nx+1,0,ny+1,0,nz+1)

    IF (WF == 1) THEN   ! NEW ZX
        DO nIter = 1, totIter
          CALL MG_itSolver(var, rhs, n, 1, 1, mx, my, mz)
        END DO ! nIter   
    ENDIF

    IF ( n == nlev ) THEN   ! NEW ZX
        IF (mgCycle == 3) GMM=1
        return
    ENDIF
        
    CALL MG_itResidual(var, rhs, r, n, 1, 1, mx, my, mz)

    CALL MG_Restrict_FMG(r, r1, n, mx, my, mz, mgrid_I(n+1)+1)

!	CALL write_dump_debug_MG('RtrX',10+n,r1, n+1, 1, 1, 0,mx1,0,my,0,mz)

    PHI = zero    ! NEW ZX

    DO I=1, MERGE(GMM, 1, N /= nlev - 1)     ! NEW ZX
        IF ( n == nlev - 1 ) THEN
            CALL MG_FMG(phi, r1, n+1, nlev, mgrid_I(n+1)+1, mgrid_I(n+1)+1, mz, 0, I, GMM, mgCycle, MG_itSolver, MG_itResidual)
        ELSE
            CALL MG_FMG(phi, r1, n+1, nlev, mgrid_I(n+1)+1, mgrid_I(n+1)+1, mz, mgrid_I(n+2)+1, I, GMM, mgCycle, MG_itSolver, MG_itResidual)
        ENDIF
    END DO
    
    iblank_MG => MGX(n+1)%iblank

    CALL MG_Prolong_FMG(r, phi, n, mgrid_I(n+1)+1, my, mz, mx)
        
    IF (N > 1) THEN
        iblank_MG => MGX(n)%iblank
    ELSE
        iblank_MG => iblank
    ENDIF

    var(1:mx-1,1:my-1,1:mz-1) = var(1:mx-1,1:my-1,1:mz-1) +r(1:mx-1,1:my-1,1:mz-1) &
                         *(oned - REAL(iblank_MG(1:mx-1,1:my-1,1:mz-1),KIND=CGREAL))

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
      CALL MG_Prepare_BC(n)

      IF ( n /= 1 ) THEN
        CALL MG_Prepare( n, mx-1, my-1, mz-1)
      ENDIF ! n
    ENDIF

    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

    DO nIter = 1, totIter
      CALL MG_itSolver(var, rhs, n, 1, 1, mx, my, mz)
    END DO ! nIter
       
    IF (N == 1 .AND. mgCycle == 3) GMM=2  ! NEW ZX

    IF (infoconv == 1) THEN
        CALL MG_itResidual(var, rhs, r, n, 1, 1, mx-1, my-1, mz-1)
    ENDIF ! infoConv

END SUBROUTINE MG_FMG


RECURSIVE SUBROUTINE MG_X(var, rhs, n, nlev, mx, my, mz, mx1, WF, GMM, mgCycle, MG_itSolver, MG_itResidual) ! NEW ZX

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual

    INTEGER, INTENT(IN) :: n, nlev, mx, my, mz, mx1, WF, mgCycle   ! NEW ZX 
    INTEGER, INTENT(INOUT) :: GMM
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: rhs
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: var
    
    INTEGER  :: i, j, k, nIter, totIter, iErr
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz) :: r
    REAL(KIND=CGREAL),DIMENSION(0:mx1,0:my,0:mz) :: r1, phi

!------------------------------------------------------------------------------
!   Initialize variable
!------------------------------------------------------------------------------

    IF (N > 1) THEN
        iblank_MG => MGX(n)%iblank
    ELSE
        iblank_MG => iblank
    ENDIF

    IF (n == 1) THEN
        totIter = iterFinest

    ELSE ! n > 1
        
        IF (n == mgLevels_X) THEN
          totIter = iterCoarsest
        ELSE
          totIter = 1
        ENDIF ! n==mgLevels_X

    ENDIF ! n==1

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
      CALL MG_Prepare_BC(n)

      IF ( n /= 1 ) THEN
        CALL MG_Prepare( n, mx-1, my-1, mz-1)
      ENDIF ! n
    ENDIF

    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

!	CALL write_dump_debug_MG('Bvar',10+n,var,n, 1, 1, 0,mx,0,my,0,mz)
!	CALL write_dump_debug_i_MG('ium ',10+n,ium_mg, n, 1, 1, 0,nx+1,0,ny+1,0,nz+1)
!	CALL write_dump_debug_i_MG('iup ',10+n,iup_mg, n, 1, 1, 0,nx+1,0,ny+1,0,nz+1)
!	CALL write_dump_debug_i_MG('jum ',10+n,jum_mg, n, 1, 1, 0,nx+1,0,ny+1,0,nz+1)
!	CALL write_dump_debug_i_MG('jup ',10+n,jup_mg, n, 1, 1, 0,nx+1,0,ny+1,0,nz+1)

!	  CALL write_dump_debug_MG('var ',0,var, n, 1, 1, 0,mx,0,my,0,mz)

    IF (WF == 1) THEN   ! NEW ZX
        DO nIter = 1, totIter
          CALL MG_itSolver(var, rhs, n, 1, 1, mx, my, mz)
        END DO ! nIter   
    ENDIF

    IF ( n == nlev ) THEN   ! NEW ZX
        IF (mgCycle == 3) GMM=1
        return
    ENDIF
        
    CALL MG_itResidual(var, rhs, r, n, 1, 1, mx, my, mz)

!	  CALL write_dump_debug_MG('var ',1,var, n, 1, 1, mx-1,my-1,0,mx,0,my,0,mz)
!	  CALL write_dump_debug_MG('R   ',1,r, n, 1, 1, mx-1,my-1, 0,mx,0,my,0,mz)
	
!	IF (Conformal_Mapping) THEN
!    CALL MG_Restrict_X_Conformal_Mapping(r, r1, n, mx, my, mz, mgrid_I(n+1)+1)
!	ELSE
    CALL MG_Restrict_X(r, r1, n, mx, my, mz, mgrid_I(n+1)+1)
!	ENDIF

!	  CALL write_dump_debug_MG('RtrX',1,r1, n+1, 1, 1, mx1-1,my-1, 0,mx1,0,my,0,mz)
!    stop
    PHI = zero    ! NEW ZX

    DO I=1, MERGE(GMM, 1, N /= nlev - 1)     ! NEW ZX
        IF ( n == nlev - 1 ) THEN
            CALL MG_X(phi, r1, n+1, nlev, mgrid_I(n+1)+1, my, mz, 0, I, GMM, mgCycle, MG_itSolver, MG_itResidual)
        ELSE
            CALL MG_X(phi, r1, n+1, nlev, mgrid_I(n+1)+1, my, mz, mgrid_I(n+2)+1, I, GMM, mgCycle, MG_itSolver, MG_itResidual)
        ENDIF
    END DO
    
    iblank_MG => MGX(n+1)%iblank

    IF (Conformal_Mapping) THEN
    CALL MG_Prolong_X_Conformal_Mapping(r, phi, n, mgrid_I(n+1)+1, my, mz, mx)
    ELSE
    CALL MG_Prolong_X(r, phi, n, mgrid_I(n+1)+1, my, mz, mx)
    ENDIF

!	  CALL write_dump_debug_MG('ProX',1,r, n, 1, 1, mx-1,my-1, 0,mx,0,my,0,mz)
        
    IF (N > 1) THEN
        iblank_MG => MGX(n)%iblank
    ELSE
        iblank_MG => iblank
    ENDIF

    var(1:mx-1,1:my-1,1:mz-1) = var(1:mx-1,1:my-1,1:mz-1) +r(1:mx-1,1:my-1,1:mz-1) &
                         *(oned - REAL(iblank_MG(1:mx-1,1:my-1,1:mz-1),KIND=CGREAL))

!	  CALL write_dump_debug_MG('var ',2,var, n, 1, 1, 0,mx,0,my,0,mz)

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
      CALL MG_Prepare_BC(n)

      IF ( n /= 1 ) THEN
        CALL MG_Prepare( n, mx-1, my-1, mz-1)
      ENDIF ! n
    ENDIF

    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

    IF (n == 1) THEN
        totIter = 0
    ELSE
        totIter = 1
    ENDIF ! n

    DO nIter = 1, totIter
      CALL MG_itSolver(var, rhs, n, 1, 1, mx, my, mz)
    END DO ! nIter

!	  CALL write_dump_debug_MG('var ',3,var, n, 1, 1, 0,mx,0,my,0,mz)
      
    IF (N == 1 .AND. mgCycle == 3) GMM=2  ! NEW ZX

    IF (infoconv == 1) THEN
        CALL MG_itResidual(var, rhs, r, n, 1, 1, mx-1, my-1, mz-1)
    ENDIF ! infoConv

END SUBROUTINE MG_X


RECURSIVE SUBROUTINE MG_Y(var, rhs, n, nlev, mx, my, mz, my1, WF, GMM, mgCycle, MG_itSolver, MG_itResidual) ! NEW ZX

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual

    INTEGER, INTENT(IN) :: n, nlev, mx, my, mz, my1, WF, mgCycle   ! NEW ZX
    INTEGER, INTENT(INOUT) :: GMM
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: rhs
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: var
    
    INTEGER  :: i, j, k, nIter, totIter, iErr
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz) :: r
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my1,0:mz) :: r1, phi

!------------------------------------------------------------------------------
!   Initialize variable
!------------------------------------------------------------------------------
    IF (N > 1) THEN
        iblank_MG => MGY(n)%iblank
    ELSE
        iblank_MG => iblank
    ENDIF

    IF (n == 1) THEN
        totIter = iterFinest

    ELSE ! n > 1
        
        IF (n == mgLevels_Y) THEN
          totIter = iterCoarsest
        ELSE
          totIter = 1
        ENDIF ! n==mgLevels_Y

    ENDIF ! n==1

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
      CALL MG_Prepare_BC(n)

      IF ( n /= 1 ) THEN
        CALL MG_Prepare( n, mx-1, my-1, mz-1)
      ENDIF ! n
    ENDIF

    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

    IF (WF == 1) THEN   ! NEW ZX
        DO nIter = 1, totIter
          CALL MG_itSolver(var, rhs, 1, n, 1, mx, my, mz) 
        END DO ! nIter
    ENDIF
    
    IF ( n == nlev ) THEN   ! NEW ZX
        IF (mgCycle == 3) GMM=1
        return
    ENDIF
        
    CALL MG_itResidual(var, rhs, r, 1, n, 1, mx, my, mz)

!	CALL write_dump_debug_MG('var ',20+n,var,1, n, 1, 0,mx,0,my,0,mz)

!	IF (Conformal_Mapping) THEN
!    CALL MG_Restrict_Y_Conformal_Mapping(r, r1, n, mx, my, mz, mgrid_J(n+1)+1)
!	ELSE
    CALL MG_Restrict_Y(r, r1, n, mx, my, mz, mgrid_J(n+1)+1)
!    ENDIF
    
    PHI = zero    ! NEW ZX

    DO I=1, MERGE(GMM, 1, N /= nlev - 1)     ! NEW ZX
        IF ( n == nlev - 1 ) THEN
            CALL MG_Y(phi, r1, n+1, nlev, mx, mgrid_J(n+1)+1, mz, 0, I, GMM, mgCycle, MG_itSolver, MG_itResidual)
        ELSE
            CALL MG_Y(phi, r1, n+1, nlev, mx, mgrid_J(n+1)+1, mz, mgrid_J(n+2)+1, I, GMM, mgCycle, MG_itSolver, MG_itResidual)
        ENDIF
    END DO
    
    iblank_MG => MGY(n+1)%iblank
    
    IF (Conformal_Mapping) THEN
    CALL MG_Prolong_Y_Conformal_Mapping(r, phi, n, mx, mgrid_J(n+1)+1, mz, my)
    ELSE
    CALL MG_Prolong_Y(r, phi, n, mx, mgrid_J(n+1)+1, mz, my)
    ENDIF
    
    IF (N > 1) THEN
        iblank_MG => MGY(n)%iblank
    ELSE
        iblank_MG => iblank
    ENDIF

    var(1:mx-1,1:my-1,1:mz-1) = var(1:mx-1,1:my-1,1:mz-1) +r(1:mx-1,1:my-1,1:mz-1) &
                                *(oned - REAL(iblank_MG(1:mx-1,1:my-1,1:mz-1),KIND=CGREAL))

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
      CALL MG_Prepare_BC(n) 

      IF ( n /= 1 ) THEN
        CALL MG_Prepare( n, mx-1, my-1, mz-1)
      ENDIF ! n
    ENDIF

    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

    IF (n == 1) THEN
        IF ( ndim == 2) THEN
            totIter = iterFinest
        ELSE
            totIter = 0
        ENDIF ! ndim
    ELSE
        totIter = 1 
    ENDIF ! n

    DO nIter = 1, totIter
        CALL MG_itSolver(var, rhs, 1, n, 1, mx, my, mz) 
    END DO ! nIter   
       
    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

    IF (N == 1 .AND. mgCycle == 3) GMM=2  ! NEW ZX

    IF (infoconv == 1) THEN
        CALL MG_itResidual(var, rhs, r, n, 1, 1, mx-1, my-1, mz-1)
    ENDIF ! infoConv

END SUBROUTINE MG_Y
!------------------------------------------------------------

RECURSIVE SUBROUTINE MG_Z(var, rhs, n, nlev, mx, my, mz, mz1, WF, GMM, mgCycle, MG_itSolver, MG_itResidual) ! NEW ZX

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    EXTERNAL MG_itSolver, MG_itResidual

    INTEGER, INTENT(IN) :: n, nlev, mx, my, mz, mz1, WF, mgCycle   ! NEW ZX
    INTEGER, INTENT(INOUT) :: GMM
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: rhs
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz),INTENT (INOUT) :: var
    
    INTEGER  :: i, j, k, nIter, totIter, iErr
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz) :: r
    REAL(KIND=CGREAL),DIMENSION(0:mx,0:my,0:mz1) :: r1, phi

!------------------------------------------------------------------------------
!   Initialize variable
!------------------------------------------------------------------------------

    IF (N > 1) THEN
        iblank_MG => MGZ(n)%iblank
    ELSE
        iblank_MG => iblank
    ENDIF
        
    IF (n == 1) THEN
        totIter = iterFinest

    ELSE ! n > 1
        
        IF (n == mgLevels_Z) THEN
          totIter = iterCoarsest
        ELSE
          totIter = 1
        ENDIF ! n==mgLevels_Z

    ENDIF ! n==1

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
      CALL MG_Prepare_BC(n) 

      IF ( n /= 1 ) THEN
        CALL MG_Prepare( n, mx-1, my-1, mz-1)
      ENDIF ! n
    ENDIF
    
    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

    IF (WF == 1) THEN   ! NEW ZX
        DO nIter = 1, totIter
          CALL MG_itSolver(var, rhs, 1, 1, n, mx, my, mz) 
        END DO ! nIter   
    ENDIF

    IF ( n == nlev ) THEN   ! NEW ZX
        IF (mgCycle == 3) GMM=1
        return
    ENDIF
        
    CALL MG_itResidual(var, rhs, r, 1, 1, n, mx, my, mz)

!   IF (Conformal_Mapping) THEN
!    CALL MG_Restrict_Z_Conformal_Mapping(r, r1, n, mx, my, mz, mgrid_K(n+1)+1)
!   ELSE
    CALL MG_Restrict_Z(r, r1, n, mx, my, mz, mgrid_K(n+1)+1)
!    ENDIF
    
    PHI = zero    ! NEW ZX

    DO I=1, MERGE(GMM, 1, N /= nlev - 1)     ! NEW ZX
        IF ( n == nlev - 1 ) THEN
            CALL MG_Z(phi, r1, n+1, nlev, mx, my, mgrid_K(n+1)+1, 0, I, GMM, mgCycle, MG_itSolver, MG_itResidual)
        ELSE
            CALL MG_Z(phi, r1, n+1, nlev, mx, my, mgrid_K(n+1)+1, mgrid_K(n+2)+1, I, GMM, mgCycle, MG_itSolver, MG_itResidual)
        ENDIF
    END DO
    
    iblank_MG => MGZ(n+1)%iblank
    
    IF (Conformal_Mapping) THEN
    CALL MG_Prolong_Z_Conformal_Mapping(r, phi, n, mx, my, mgrid_K(n+1)+1, mz)
    ELSE
    CALL MG_Prolong_Z(r, phi, n, mx, my, mgrid_K(n+1)+1, mz)
    ENDIF
    
    IF (N > 1) THEN
        iblank_MG => MGZ(n)%iblank
    ELSE
        iblank_MG => iblank
    ENDIF

    var(1:mx-1,1:my-1,1:mz-1) = var(1:mx-1,1:my-1,1:mz-1) +r(1:mx-1,1:my-1,1:mz-1) &
                            *(oned - REAL(iblank_MG(1:mx-1,1:my-1,1:mz-1),KIND=CGREAL))

    IF (n == 1) THEN
        totIter = iterFinest
    ELSE
        totIter = 1 
    ENDIF ! n

    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) THEN
      CALL MG_Prepare_BC(n)
      IF ( n /= 1 ) THEN
        CALL MG_Prepare( n, mx-1, my-1, mz-1)
      ENDIF ! n
    ENDIF
 
!    DO nIter = 1, totIter
!        CALL MG_itsolv(var, rhs, 1, 1, n, mx, my, mz) 
!    END DO ! nIter   

    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

    DO nIter = 1, totIter
      CALL MG_itSolver(var, rhs, 1, 1, n, mx, my, mz) 
    END DO ! nIter   
       
    IF ( n == 1 .and. boundary_formulation == GCM_METHOD ) CALL GCM_Pressure(var, rhs)

    IF (N == 1 .AND. mgCycle == 3)  GMM=2  ! NEW ZX
    
    IF (infoconv == 1) THEN
        CALL MG_itResidual(var, rhs, r, n, 1, 1, mx-1, my-1, mz-1)
    ENDIF ! infoConv

END SUBROUTINE MG_Z
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE MG_Prepare_BC(n)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays 
    USE MG_parameters
    USE MG_arrays
   
    IMPLICIT NONE

    INTEGER,INTENT (IN) 	:: N
    INTEGER :: i, j, k

    IF (N .EQ. 1) THEN

        ium_MG = ium
        iup_MG = iup
        jum_MG = jum
        jup_MG = jup
        kum_MG = kum
        kup_MG = kup
        
      IF (internal_boundary_present.EQ.1) THEN
          ghostcellMark_MG = ghostcellMark      !!! Added by ZX Liang
      ENDIF

      IF ( boundary_formulation == GCM_METHOD ) THEN
        CALL GCM_correct_pressure
      ELSE
        nlw=zero
      ENDIF
    ELSE

      IF (internal_boundary_present .EQ. 1) THEN

        ghostcellMark_MG = 0       !!! Added by ZX Liang

        nlw=zero
 
      ENDIF
       
    ENDIF

END SUBROUTINE MG_Prepare_BC
!---------------------------------------------------------------------

SUBROUTINE MG_Prepare(nlev,mx,my,mz) 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE
 
    INTEGER, INTENT(IN) :: nlev, mx, my, mz

    INTEGER           :: i, j, k
    REAL(KIND=CGREAL) :: ttemp1, ttemp2, ttemp3, ttemp4

    IF ( internal_boundary_present == 1 ) THEN
      DO k = 1, mz  

        DO j = 1, my 
        DO i = 1, mx-1
            ttemp1 = REAL((iblank_MG(i+1,j,k) - iblank_MG(i,j,k)),KIND=CGREAL)/2.0_CGREAL
            iup_MG(i,j,k) = REAL(INT(ttemp1 + half), KIND=CGREAL) 
        END DO ! i
               
        DO i = 2, mx
            ttemp2 = REAL((iblank_MG(i-1,j,k) - iblank_MG(i,j,k)),KIND=CGREAL)/2.0_CGREAL
            ium_mg(i,j,k) = REAL(INT(ttemp2 + half), KIND=CGREAL)
        END DO ! i
        END DO ! j

        DO j = 1, my-1
        DO i = 1, mx
            ttemp1 = REAL((iblank_MG(i,j+1,k) - iblank_MG(i,j,k)),KIND=CGREAL)/2.0_CGREAL
            jup_MG(i,j,k) = REAL(INT(ttemp1 + half), KIND=CGREAL)
        END DO ! i
        END DO ! j

        DO j = 2, my 
        DO i = 1, mx
            ttemp2 = REAL((iblank_MG(i,j-1,k) - iblank_MG(i,j,k)),KIND=CGREAL)/2.0_CGREAL
            jum_MG(i,j,k) = REAL(INT(ttemp2 + half), KIND=CGREAL)
        END DO ! i
        END DO ! j

      END DO !  K

      DO j = 1, my
      DO i = 1, mx
        DO k = 1, mz-1
            ttemp1 = REAL((iblank_MG(i,j,k+1) - iblank_MG(i,j,k)),KIND=CGREAL)/2.0_CGREAL
            kup_MG(i,j,k) = REAL(INT(ttemp1 + half), KIND=CGREAL)
        END DO ! k

        DO k = 2, mz
            ttemp2 = REAL((iblank_MG(i,j,k-1) - iblank_MG(i,j,k)),KIND=CGREAL)/2.0_CGREAL
            kum_MG(i,j,k) = REAL(INT(ttemp2 + half), KIND=CGREAL)
        END DO ! k
      END DO ! i
      END DO ! j

    ENDIF ! internal_boundary_present

    IF ( bcx1 == BC_TYPE_PERIODIC .AND. bcx2 == BC_TYPE_PERIODIC ) THEN
      DO k=1,mz
      DO j=1,my
        ium_MG(1,j,k)  = 0
        iup_MG(mx,j,k) = 0
      ENDDO ! j
      ENDDO ! k
    ELSE
      DO k=1,mz
      DO j=1,my
        ium_MG(1,j,k)  = 1
        iup_MG(mx,j,k) = 1
      ENDDO ! j
      ENDDO ! k
    ENDIF ! bcx1

    IF ( pbcx1 == PBC_DIRICHLET ) THEN
      DO k=1,mz
      DO j=1,my
        ium_MG(1,j,k)  = - 1
      ENDDO ! j
      ENDDO ! k
    ENDIF ! pbcx1

    IF ( pbcx2 == PBC_DIRICHLET ) THEN
      DO k=1,mz
      DO j=1,my
        iup_MG(mx,j,k) = - 1
      ENDDO ! j
      ENDDO ! k
    ENDIF ! pbcx2

    IF ( bcy1 == BC_TYPE_PERIODIC .AND. bcy2 == BC_TYPE_PERIODIC ) THEN
      DO k=1,mz
      DO i=1,mx
        jum_MG(i,1,k)  = 0
        jup_MG(i,my,k) = 0
      ENDDO ! i
      ENDDO ! k
    ELSE
      DO k=1,mz
      DO i=1,mx
        jum_MG(i,1,k)  = 1
        jup_MG(i,my,k) = 1
      ENDDO ! i
      ENDDO ! k
    ENDIF ! bcy1

    IF ( pbcy1 == PBC_DIRICHLET ) THEN
      DO k=1,mz
      DO i=1,mx
        jum_MG(i,1,k)  = - 1
      ENDDO ! i
      ENDDO ! k
    ENDIF ! pbcy1

    IF ( pbcy2 == PBC_DIRICHLET ) THEN
      DO k=1,mz
      DO i=1,mx
        jup_MG(i,my,k) = - 1
      ENDDO ! j
      ENDDO ! k
    ENDIF ! pbcy2
    
    IF ( bcz1 == BC_TYPE_PERIODIC .AND. bcz2 == BC_TYPE_PERIODIC ) THEN
      DO j=1,my
      DO i=1,mx
        kum_MG(i,j,1)  = 0
        kup_MG(i,j,mz) = 0
      ENDDO ! i
      ENDDO ! j
    ELSE
      DO j=1,my
      DO i=1,mx
        kum_MG(i,j,1)  = 1
        kup_MG(i,j,mz) = 1
      ENDDO ! i
      ENDDO ! j
    ENDIF ! bcz1

    IF ( pbcz1 == PBC_DIRICHLET ) THEN
      DO j=1,my
      DO i=1,mx
        kum_MG(i,j,1)  = - 1
      ENDDO ! i
      ENDDO ! j
    ENDIF ! pbcz1

    IF ( pbcz2 == PBC_DIRICHLET ) THEN
      DO j=1,my
      DO i=1,mx
        kup_MG(i,j,mz) = - 1
      ENDDO ! i
      ENDDO ! j
    ENDIF ! pbcy2

END SUBROUTINE MG_Prepare
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE MG_itsolv(var,R, nLevX,nLevY,nLevZ,mx,my,mz) 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nLevX, nLevY, nLevZ, mx, my, mz
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)     :: R
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (INOUT)  :: var

    INTEGER :: i,j,k, nLev
    INTEGER :: iBody, iRow, iG,jG, iNode, jNode, n
    INTEGER :: iinit, jinit, kinit, Ncolor   
    REAL(KIND=CGREAL)  :: AAe, AAw, AAn, AAs, AAf, AAb


!-----  Line solve in the x-direction -------------------
 
! IF statement for pressure periodic boundary conditions
!    if (nlevx==1 .and. nlevy==3) then
!      call write_dump_debug_MG('var ',0,var,nlevx,nlevy,nlevz,mx,my,mz,0,mx,0,my,0,mz)
!    endif
 
    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
       CALL enforce_p_periodic(var)
    ENDIF
 
  DO Ncolor = 1, TNcolorX
    kinit = (Ncolor-1)*(ndim-DIM_2D)+ 1
    jinit = Ncolor
 
    DO k = kinit, mz-1, kStep
    DO j = jinit, my-1, jStep
      DO i=1, mx-1 
       IF (cure_pressure_oscillations .AND. iblank_MG(i,j,k)==0 .AND. ivc(i,j,k)>0 .AND. nLevX==1 .AND. nLevY==1 .AND. nLevZ==1) THEN
          AAe=face(cell(ivc(i,j,k))%F_ip)%a
          AAw=face(cell(ivc(i,j,k))%F_im)%a
          AAn=face(cell(ivc(i,j,k))%F_jp)%a
          AAs=face(cell(ivc(i,j,k))%F_jm)%a
          AAf=face(cell(ivc(i,j,k))%F_kp)%a
          AAb=face(cell(ivc(i,j,k))%F_km)%a
       amx(i) =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) ) * AAw
       apx(i) =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) ) * AAe
       acx(i) = - ( amx(i) + apx(i) )
 
       amy(j) =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) ) * AAs
       apy(j) =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) ) * AAn
       acy(j) = - ( amy(j) + apy(j) )
 
       amz(k) =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL) * AAb
       apz(k) =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL) * AAf
       acz(k) = - ( amz(k) + apz(k) )
       ELSE
       amx(i) =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) )
       apx(i) =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) )
       acx(i) = - ( amx(i) + apx(i) )
 
       amy(j) =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) )
       apy(j) =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) )
       acy(j) = - ( amy(j) + apy(j) )
 
       amz(k) =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL)
       apz(k) =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL)
       acz(k) = - ( amz(k) + apz(k) )
       ENDIF

       rhs(i) = r(i,j,k) - var(i,j-1,k)*amy(j)  &
                         - var(i,j+1,k)*apy(j)  &
                         - var(i,j,k-1)*amz(k)  &
                         - var(i,j,k+1)*apz(k)  &
                         + nlw(i,j,k)
 
       amx(i) = amx(i)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
       apx(i) = apx(i)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
       acx(i) = (acx(i)+acy(j)+acz(k))*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL)  ) &
               + REAL(iblank_MG(i,j,k),KIND=CGREAL)
       rhs(i) = rhs(i)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) ) &
                + REAL(ghostcellMark_mg(i,j,k),KIND=CGREAL)*var(i,j,k)
        
       IF(abs(acx(i)).LT.1.0e-12) THEN
             amx(i) = zero 
             apx(i)  = zero      
             acx(i) = oned 
             rhs(i) = var(i,j,k)
       ENDIF
        
      ENDDO ! i           
                          
! Pressure periodic boundary condition in x direction
! Remove two terms to the right hand side. 
                          
      IF (bcx1 == BC_TYPE_PERIODIC .AND. &
          bcx2 == BC_TYPE_PERIODIC) THEN
         rhs(1) = rhs(1) - var(nxc,j,k)*amx(1)
         rhs(nxc) = rhs(nxc) - var(1,j,k)*apx(nxc)
      ENDIF
        
      IF (pbcx1 == PBC_DIRICHLET .and.      &           !Added by H. Luo
          nlevx*nlevy*nlevz == 1) THEN                  !
         rhs(1 ) = rhs(1 ) - pppx1 * amx(1 )            !The non-zero B.C. is only for
      ENDIF 

      IF (pbcx2 == PBC_DIRICHLET .and.      &           !
          nlevx*nlevy*nlevz == 1) THEN                  !
         rhs(mx-1) = rhs(mx-1) - pppx2 * apx(mx-1)            !the finest level.
      ENDIF
 
      CALL tdma(amx,acx,apx,rhs,dummy,1,mx-1)

      DO i=1,mx-1
        var(i,j,k) = var(i,j,k)*(oned-omega) + omega*dummy(i)
      ENDDO
 
!      if (nlevx==1 .and. nlevy==2) then
!        do i=1, mx-1
!          write(225,'(3I5,6(2X,1PE22.15))') i, j, k, amx(i), acx(i), apx(i), rhs(i), dummy(i), var(i,j,k)
!        enddo
!!        stop
!      endif

      IF (pbcx1 == PBC_DIRICHLET .and.      &    !
          nlevx*nlevy*nlevz == 1) THEN           !
         var(0,   j,k) = pppx1                   !Added by H. Luo
      ENDIF
 
      IF (pbcx2 == PBC_DIRICHLET .and.      &    !
          nlevx*nlevy*nlevz == 1) THEN           !
         var(mx,j,k) = pppx2                   !
      ENDIF
 
    ENDDO ! j
    ENDDO ! k

    ENDDO ! NcolorX 
           
    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
       CALL enforce_p_periodic(var) 
    ENDIF
       
!    if (nlevx==1 .and. nlevy==2) then
!      call write_dump_debug_MG('var ',1,var,nlevx,nlevy,nlevz,mx,my,mz,0,mx,0,my,0,mz)
!    endif

  DO Ncolor = 1, TNcolorY
 
    kinit = (Ncolor-1)*(ndim-DIM_2D)+ 1      
    iinit = Ncolor   
          
!------- Line solve in the y-direction ------------------
    DO k= kinit, mz-1,kstep
    DO i= iinit, mx-1,istep 
       DO j=1, my-1 
          
       IF (cure_pressure_oscillations .AND. iblank_MG(i,j,k)==0 .AND. ivc(i,j,k)>0 .AND. nLevX==1 .AND. nLevY==1 .AND. nLevZ==1) THEN
          AAe=face(cell(ivc(i,j,k))%F_ip)%a
          AAw=face(cell(ivc(i,j,k))%F_im)%a
          AAn=face(cell(ivc(i,j,k))%F_jp)%a
          AAs=face(cell(ivc(i,j,k))%F_jm)%a
          AAf=face(cell(ivc(i,j,k))%F_kp)%a
          AAb=face(cell(ivc(i,j,k))%F_km)%a
       amx(i) =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) ) * AAw
       apx(i) =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) ) * AAe
       acx(i) = - ( amx(i) + apx(i) )
 
       amy(j) =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) ) * AAs
       apy(j) =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) ) * AAn
       acy(j) = - ( amy(j) + apy(j) )
 
       amz(k) =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL) * AAb
       apz(k) =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL) * AAf
       acz(k) = - ( amz(k) + apz(k) )
       ELSE
       amx(i) =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) )
       apx(i) =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) )
       acx(i) = - ( amx(i) + apx(i) )
     
       amy(j) =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) )
       apy(j) =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) )
       acy(j) = - ( amy(j) + apy(j) )
 
       amz(k) =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL)
       apz(k) =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL)
       acz(k) = - ( amz(k) + apz(k) )
       ENDIF
        
       rhs(j) = r(i,j,k) - var(i-1,j,k)*amx(i)  &
                         - var(i+1,j,k)*apx(i)  &
                         - var(i,j,k-1)*amz(k)  &
                         - var(i,j,k+1)*apz(k)  &
                         + nlw(i,j,k)
 
!-- Modify rhs and coefficients for Ghost cells in GCM
 
       amy(j) = amy(j)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
       apy(j) = apy(j)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
       acy(j) = (acx(i)+acy(j)+acz(k))*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) ) &
               + REAL(iblank_MG(i,j,k),KIND=CGREAL)
       rhs(j) = rhs(j)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )&
                + REAL(ghostcellMark_mg(i,j,k),KIND=CGREAL)*var(i,j,k)  
                                       
       IF(abs(acy(J)).LT.1.0e-12) THEN
             amy(J) = zero       
             apy(J)  = zero  
             acy(J) = oned
             rhs(J) = var(i,j,k)  
       ENDIF
                          
      ENDDO               
                          
      IF (bcy1 == BC_TYPE_PERIODIC .AND. &
          bcy2 == BC_TYPE_PERIODIC) THEN
         rhs(1) = rhs(1) - var(i,nyc,k)*amy(1)
         rhs(nyc) = rhs(nyc) - var(i,1,k)*apy(nyc)
      ENDIF
 
      IF (pbcy1 == PBC_DIRICHLET .and.      &           !
          nlevx*nlevy*nlevz == 1) THEN                  !
         rhs(1 ) = rhs(1 ) - pppy1 * amy(1 )            !Added by H. Luo
      ENDIF 
 
      IF (pbcy2 == PBC_DIRICHLET .and.      &           !
          nlevx*nlevy*nlevz == 1) THEN                  !
         rhs(my-1) = rhs(my-1) - pppy2 * apy(my-1)            !
      ENDIF
 
      CALL tdma(amy,acy,apy,rhs,dummy,1,my-1)

!      if (nlevx==1 .and. nlevy==2) then
!        do j=1, my-1
!          write(225,'(I3,5(2X,1PE22.15))') j, amy(j), acy(j), apy(j), rhs(j), dummy(j)
!        enddo
!      endif
 
      DO j=1,my-1
       var(i,j,k) = var(i,j,k)*(oned-omega) + omega*dummy(j)
      ENDDO
 
      IF (pbcy1 == PBC_DIRICHLET .and.      &    !
          nlevx*nlevy*nlevz == 1) THEN           !
         var(i,   0,k) = pppy1                   !Added by H. Luo
      ENDIF
 
      IF (pbcy2 == PBC_DIRICHLET .and.      &    !
          nlevx*nlevy*nlevz == 1) THEN           !
         var(i,my,k) = pppy2                   !
      ENDIF
 
     ENDDO
     ENDDO
 
    ENDDO ! Ncolor  
           
!-------- Line solver in the z-direction ----------------------
    IF (ndim == DIM_3D) THEN
 
    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
       CALL enforce_p_periodic(var) 
    ENDIF
 
  DO Ncolor = 1, TNcolorZ  
           
    jinit = Ncolor 
    iinit = Ncolor
 
    DO j = jinit, my-1, jstep  
    DO i = iinit, mx-1, istep  
      DO k = 1, mz-1 
        IF (cure_pressure_oscillations .AND. iblank_MG(i,j,k)==0 .AND. ivc(i,j,k)>0 .AND. nLevX==1 .AND. nLevY==1 .AND. nLevZ==1) THEN
          AAe=face(cell(ivc(i,j,k))%F_ip)%a
          AAw=face(cell(ivc(i,j,k))%F_im)%a
          AAn=face(cell(ivc(i,j,k))%F_jp)%a
          AAs=face(cell(ivc(i,j,k))%F_jm)%a
          AAf=face(cell(ivc(i,j,k))%F_kp)%a
          AAb=face(cell(ivc(i,j,k))%F_km)%a
          amx(i) =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) ) * AAw
          apx(i) =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) ) * AAe
          acx(i) = - ( amx(i) + apx(i) )
 
          amy(j) =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) ) * AAs
          apy(j) =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) ) * AAn
          acy(j) = - ( amy(j) + apy(j) )
 
          amz(k) =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                        *REAL((ndim - DIM_2D),KIND=CGREAL) * AAb
          apz(k) =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                        *REAL((ndim - DIM_2D),KIND=CGREAL) * AAf
          acz(k) = - ( amz(k) + apz(k) )
        ELSE
       
          amx(i) =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) )
          apx(i) =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) )
          acx(i) = - ( amx(i) + apx(i) )
 
          amy(j) =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) )
          apy(j) =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) )
          acy(j) = - ( amy(j) + apy(j) )
 
          amz(k) =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                        *REAL((ndim - DIM_2D),KIND=CGREAL)
          apz(k) =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                        *REAL((ndim - DIM_2D),KIND=CGREAL)
          acz(k) = - ( amz(k) + apz(k) )
        ENDIF 
         rhs(k) = r(i,j,k) - var(i,j-1,k)*amy(j)  &
                           - var(i,j+1,k)*apy(j)  &
                           - var(i-1,j,k)*amx(i)  &
                           - var(i+1,j,k)*apx(i)  &
                           + nlw(i,j,k)
 
         amz(k) = amz(k)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
         apz(k) = apz(k)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
         acz(k) = (acx(i)+acy(j)+acz(k))*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) ) &
                + REAL(iblank_MG(i,j,k),KIND=CGREAL)        
         rhs(k) = rhs(k)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) ) & 
                + REAL(ghostcellMark_mg(i,j,k),KIND=CGREAL)*var(i,j,k)
 
         IF(abs(acz(k)).LT.1.0e-12) THEN
             amz(k) = zero           
             apz(k)  = zero   
             acz(k) = oned           
             rhs(k) = var(i,j,k)   
         ENDIF
          
        ENDDO ! k           
                            
        IF (bcz1 == BC_TYPE_PERIODIC .AND. &
            bcz2 == BC_TYPE_PERIODIC) THEN
           rhs(1) = rhs(1) - var(i,j,nzc)*amz(1)
           rhs(nzc) = rhs(nzc) - var(i,j,1)*apz(nzc)
        ENDIF
                            
      IF (pbcz1 == PBC_DIRICHLET .and.      &           !
          nlevx*nlevy*nlevz == 1) THEN                  !
         rhs(1 ) = rhs(1 ) - pppz1 * amz(1 )            !Added by H. Luo 
      ENDIF 
 
      IF (pbcz2 == PBC_DIRICHLET .and.      &           !
          nlevx*nlevy*nlevz == 1) THEN                  !
         rhs(mz-1) = rhs(mz-1) - pppz2 * apz(mz-1)            !
      ENDIF
 
        CALL tdma(amz,acz,apz,rhs,dummy,1,mz-1)
 
        DO k=1,mz-1
         var(i,j,k) =  var(i,j,k)*(oned-omega) + omega*dummy(k)
        ENDDO ! k
 
      IF (pbcz1 == PBC_DIRICHLET .and.      &    !
          nlevx*nlevy*nlevz == 1) THEN           !
         var(i,j,0   ) = pppz1                   !Added by H. Luo
      ENDIF
 
      IF (pbcz2 == PBC_DIRICHLET .and.      &    !
          nlevx*nlevy*nlevz == 1) THEN           !
         var(i,j,mz) = pppz2                   !
      ENDIF
 
      ENDDO ! i
      ENDDO ! j
 
     END DO ! Ncolor 
           
    ENDIF ! ndim  

END SUBROUTINE MG_itsolv
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE MG_residual(var,R,resL,nLevX,nLevY,nLevZ,mx, my, mz) 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nLevX, nLevY, nLevZ, mx, my, mz

    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (INOUT)  :: resL
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)     :: var, R

    INTEGER            :: i,j,k,iErr
    INTEGER            :: iG,jG,iBody,iRow,n
    REAL(KIND=CGREAL)  :: res, resMax
    REAL(KIND=CGREAL)  :: bmx,bpx,bcx,bc
    REAL(KIND=CGREAL)  :: bmy,bpy,bcy
    REAL(KIND=CGREAL)  :: bmz,bpz,bcz
    REAL(KIND=CGREAL)  :: AAe, AAw, AAn, AAs, AAf, AAb

    resMax = zero

    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
       CALL enforce_p_periodic(var) 
    ENDIF ! bcx1

!   ------------------------------------
    DO k = 1, mz-1
    DO j = 1, my-1
    DO i = 1, mx-1
       IF (cure_pressure_oscillations .AND. iblank_MG(i,j,k)==0 .AND. ivc(i,j,k)>0 .AND. nLevX==1 .AND. nLevY==1 .AND. nLevZ==1) THEN
          AAe=face(cell(ivc(i,j,k))%F_ip)%a
          AAw=face(cell(ivc(i,j,k))%F_im)%a
          AAn=face(cell(ivc(i,j,k))%F_jp)%a
          AAs=face(cell(ivc(i,j,k))%F_jm)%a
          AAf=face(cell(ivc(i,j,k))%F_kp)%a
          AAb=face(cell(ivc(i,j,k))%F_km)%a

            bmx =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) ) * AAw
            bpx =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) ) * AAe
            bcx = - ( bmx + bpx ) 
     
            bmy =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) ) * AAs
            bpy =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) ) * AAn
            bcy = - ( bmy + bpy ) 
     
            bmz =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL) * AAb
            bpz =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL) * AAf
            bcz = - ( bmz + bpz ) 
     
            bc =   (bcx+bcy+bcz)

	  ELSE
            bmx =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) )
            bpx =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) )
            bcx = - ( bmx + bpx ) 
     
            bmy =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) )
            bpy =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) )
            bcy = - ( bmy + bpy ) 
     
            bmz =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL)
            bpz =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                      *REAL((ndim - DIM_2D),KIND=CGREAL)
            bcz = - ( bmz + bpz ) 
     
            bmx = bmx*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
            bpx = bpx*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
 
            bmy = bmy*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
            bpy = bpy*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
 
            bmz = bmz*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
            bpz = bpz*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
 
            bc =   (bcx+bcy+bcz)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )  &
                   + REAL(iblank_MG(i,j,k),KIND=CGREAL)
	  ENDIF
	   
            res    = R(i,j,k) - var(i,j,k)*bc  &
                         - var(i-1,j,k)*bmx                &
                         - var(i+1,j,k)*bpx                &
                         - var(i,j-1,k)*bmy                &
                         - var(i,j+1,k)*bpy                &
                         - var(i,j,k-1)*bmz                &
                         - var(i,j,k+1)*bpz                &
                        + nlw(i,j,k)

           resL(i,j,k)  = res*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
       
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

    IF (infoconv == 1) THEN
      resMax = MAXVAL(ABS(resL))
        WRITE(*,100) nLevX, nLevY, nLevZ, resmax
    ENDIF ! infoconv

100   FORMAT('MG: residual check : ',1x,3I6,2x,E19.11)
END SUBROUTINE MG_Residual
!-------------------------------------------------------------------------------

SUBROUTINE MG_Restrict_FMG(var,r,nlev, mx, my, mz, mx1)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, mx1

    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)     :: var
    REAL(KIND=CGREAL), DIMENSION(0:mx1,0:mx1,0:mz), INTENT (OUT)    :: r

    INTEGER  :: ii, jj, kk, iim, jjm, kkm, i, j, k, iip, jjp, kkp, iErr
    REAL(KIND=CGREAL)  :: sum1, temp, area1, area2, area3, area4, area, grid_wht1, grid_wht2, grid_wht3, grid_wht4
    
    REAL(KIND=CGREAL), DIMENSION(mx,my,mz) ::  diblank_MG  
     
!   Start Restriction
!   -----------------
    
!   Restrict in x-direction
!   -----------------------
    diblank_MG = oned
    diblank_MG(1:mx,1:my,1:mz) = REAL(1-iblank_MG(1:mx,1:my,1:mz),KIND=CGREAL) 

    DO k=1,mz-1
    DO j=1,mx1-1
        jj  = 2*j
        jjm = 2*j - 1
    DO i=1,mx1-1
        ii  = 2*i
        iim = 2*i - 1

        area1 = dx_MG(iim,nlev)*dy_MG(jjm,nlev)
        area2 = dx_MG(ii ,nlev)*dy_MG(jjm,nlev)
        area3 = dx_MG(iim,nlev)*dy_MG(jj ,nlev)
        area4 = dx_MG(ii ,nlev)*dy_MG(jj ,nlev)
        area = area1 + area2 + area3 + area4

        grid_wht1 = area1/area
        grid_wht2 = area2/area
        grid_wht3 = area3/area
        grid_wht4 = area4/area

        sum1 = diblank_MG(iim,jjm,k) *grid_wht1 *var(iim,jjm,k) &
             + diblank_MG(ii ,jjm,k) *grid_wht2 *var(ii ,jjm,k) &
             + diblank_MG(iim,jj ,k) *grid_wht3 *var(iim,jj ,k) &
             + diblank_MG(ii ,jj ,k) *grid_wht4 *var(ii ,jj ,k)

        temp = diblank_MG(iim,jjm,k) *grid_wht1  &
             + diblank_MG(ii ,jjm,k) *grid_wht2  &
             + diblank_MG(iim,jj ,k) *grid_wht3  &
             + diblank_MG(ii ,jj ,k) *grid_wht4
              
        IF ( temp > zero ) THEN
          r(i,j,k) = sum1/temp
        ELSE
          r(i,j,k) = sum1
        ENDIF ! temp
    END DO ! i
    END DO ! j
    END DO ! k
   
END SUBROUTINE MG_Restrict_FMG
!---------------------------------------------------------------------

SUBROUTINE MG_Restrict_X(var,r,nlev, mx, my, mz, mx1)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, mx1

    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)     :: var
    REAL(KIND=CGREAL), DIMENSION(0:mx1,0:my,0:mz), INTENT (OUT)    :: r

    INTEGER  :: ii, jj, kk, iim, jjm, kkm, i, j, k, iip, jjp, kkp, iErr
    REAL(KIND=CGREAL)  :: sum1, temp, area1, area2, area, grid_wht1, grid_wht2
    
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE ::  diblank_MG
     
!   Allocate local array
!   --------------------
    ALLOCATE(diblank_MG(0:mx,0:my,0:mz),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Restrict_X: Memory Allocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

!   Start Restriction
!   -----------------
    
!   Restrict in x-direction
!   -----------------------
    DO k=1, mz-1
    DO j=1, my-1
    DO i=1, mx-1
        diblank_MG(i,j,k) = oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) 
    END DO ! i
    END DO ! j
    END DO ! k

    DO k=1,mz-1
    DO j=1,my-1
    DO i=1,mx1-1
        ii  = 2*i
        iim = 2*i -1
        iip = 2*i +1

        area1 = x_MG(ii ,nlev) -x_MG(iim,nlev)
        area2 = x_MG(iip,nlev) -x_MG(ii ,nlev)
        area = area1 + area2

        grid_wht1 = area1/area
        grid_wht2 = area2/area

        sum1 = diblank_MG(ii ,j,k) *grid_wht2 *var(ii ,j,k) &
             + diblank_MG(iim,j,k) *grid_wht1 *var(iim,j,k) 

        temp = diblank_MG(ii ,j,k) *grid_wht2  &
             + diblank_MG(iim,j,k) *grid_wht1 

        IF ( temp > zero ) THEN
          r(i,j,k) = sum1/temp
        ELSE
          r(i,j,k) = sum1
        ENDIF ! temp
    END DO ! i
    END DO ! j
    END DO ! k

!   Deallocate local array
!   ----------------------
    DEALLOCATE(diblank_MG,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Restrict_X: Memory Deallocation Error for diblank_MG'
      STOP
    ENDIF ! iErr
   
END SUBROUTINE MG_Restrict_X
!---------------------------------------------------------------------

SUBROUTINE MG_Restrict_Y(var,r,nlev,mx, my, mz, my1)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, my1

    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)   :: var
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my1,0:mz), INTENT (OUT) :: r

    INTEGER  :: ii, jj, kk, iim, jjm, kkm, i, j, k, iip, jjp, kkp, iErr
    REAL(KIND=CGREAL)  :: sum1, temp, area1, area2, area, grid_wht1, grid_wht2
    
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE ::  diblank_MG

!   Allocate local array
!   --------------------
    ALLOCATE(diblank_MG(0:mx,0:my,0:mz),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Restrict: Memory Allocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

!   Start Restriction
!   -----------------
!   Restrict in y-direction
!   -----------------------
      DO k = 1, mz-1
      DO j = 1, my-1
      DO i = 1, mx-1 
        diblank_MG(i,j,k) = oned-REAL(iblank_MG(i,j,k),KIND=CGREAL)
      END DO ! i
      END DO ! j
      END DO ! k

      DO k=1, mz-1
      DO j=1, my1-1
        jj  = 2*j
        jjm = 2*j -1
        jjp = 2*j +1

        area1 = y_MG(jj ,nlev) -y_MG(jjm,nlev)
        area2 = y_MG(jjp,nlev) -y_MG(jj ,nlev)
        area = area1 + area2 

        grid_wht1 = area1/area
        grid_wht2 = area2/area

        DO i=1, mx-1 
          sum1 = diblank_MG(i,jj ,k) *grid_wht2 *var(i,jj ,k) &
               + diblank_MG(i,jjm,k) *grid_wht1 *var(i,jjm,k) 
            
          temp = diblank_MG(i,jjm,k) *grid_wht1  &
               + diblank_MG(i,jj ,k) *grid_wht2 
            
          IF ( temp > zero ) THEN
            r(i,j,k) = sum1/temp
          ELSE
            r(i,j,k) = sum1
          ENDIF ! temp
        END DO ! i
      END DO ! j
      END DO ! k

!   Deallocate local array
!   ----------------------
    DEALLOCATE(diblank_MG,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Restrict: Memory Deallocation Error for diblank_MG'
      STOP
    ENDIF ! iErr
   
END SUBROUTINE MG_Restrict_Y
!---------------------------------------------------------------------

SUBROUTINE MG_Restrict_Z(var,r,nlev,mx, my, mz, mz1)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, mz1

    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)    :: var
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz1),  INTENT (OUT) :: r

    INTEGER  :: ii, jj, kk, iim, jjm, kkm, i, j, k, iip, jjp, kkp, iErr
    REAL(KIND=CGREAL)  :: sum1, temp, area1, area2, area, grid_wht1, grid_wht2
    
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE ::  diblank_MG


!   Allocate local array
!   --------------------
    ALLOCATE(diblank_MG(0:mx,0:my,0:mz),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Restrict: Memory Allocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

!   Start Restriction
!   -----------------
!   Restrict in z-direction
!   -----------------------
    DO k = 1, mz-1
    DO j = 1, my-1
    DO i = 1, mx-1       
        diblank_MG(i,j,k) = oned-REAL(iblank_MG(i,j,k),KIND=CGREAL)
    END DO ! i
    END DO ! j
    END DO ! k

    DO k = 1, mz1-1
        kk  = 2*k
        kkm = 2*k -1
        kkp = 2*k +1

        area1 = z_MG(kk ,nlev) -z_MG(kkm,nlev)
        area2 = z_MG(kkp,nlev) -z_MG(kk ,nlev)
        area = area1 + area2  

        grid_wht1 = area1/area
        grid_wht2 = area2/area

        DO j = 1, my-1
        DO i = 1, mx-1
          sum1 = diblank_MG(i,j,kk ) *grid_wht2 *var(i,j,kk )  &
               + diblank_MG(i,j,kkm) *grid_wht1 *var(i,j,kkm) 

          temp = diblank_MG(i,j,kkm) *grid_wht1 &
               + diblank_MG(i,j,kk ) *grid_wht2 
          IF ( temp > zero ) THEN
            r(i,j,k) = sum1/temp
          ELSE
            r(i,j,k) = sum1
          ENDIF ! temp
        END DO ! i
        END DO ! j
    END DO ! k

!   Deallocate local array
!   ----------------------
    DEALLOCATE(diblank_MG,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Restrict: Memory Deallocation Error for diblank_MG'
      STOP
    ENDIF ! iErr
   
END SUBROUTINE MG_Restrict_Z
!---------------------------------------------------------------------

SUBROUTINE MG_Prolong_FMG(phi,correc,nlev, mx1, my, mz, mx)
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, mx1
 
    REAL(KIND=CGREAL), DIMENSION(0:mx1,0:mx1,0:mz), INTENT (IN)  :: correc
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (OUT) :: phi

    INTEGER           :: ii, iim, jj, jjm, kk, kkm, i, j, k, nxyzMax, iErr
    REAL(KIND=CGREAL) :: var1, var2, var3, var4
    REAL(KIND=CGREAL) :: delta1, delta2, delta3, delta4, delta, d_1, d_2
    REAL(KIND=CGREAL), DIMENSION(0:mx1+1) :: pro_wht1, pro_wht2, pro_wht3, pro_wht4
    REAL(KIND=CGREAL), DIMENSION(0:mx1+1) :: pro_wht5, pro_wht6, pro_wht7, pro_wht8
    REAL(KIND=CGREAL), DIMENSION(0:mx1,0:mx1,0:mz) :: diblank_MG

!   Start Prolongation
!   ------------------
!   Prolongate in x-direction
!   -------------------------
    diblank_MG = oned
    diblank_MG(1:mx1-1,1:mx1-1,1:mz-1) = REAL(1-iblank_MG(1:mx1-1,1:mx1-1,1:mz-1),KIND=CGREAL)

    DO i = 1, mx1-1 !mgrid_I(nlev+1)
        ii  = 2*i
        iim = 2*i -1

        d_1 = oned/(xc_MG(i  ,nlev+1) -xc_MG(iim,nlev)) !**2
        d_2 = oned/(xc_MG(iim,nlev) -xc_MG(i-1,nlev+1)) !**2
        delta = d_1 +d_2

        pro_wht1(i) = d_1/delta
        pro_wht2(i) = d_2/delta

        d_1 = oned/(xc_MG(ii,nlev) -xc_MG(i  ,nlev+1)) !**2
        d_2 = oned/(xc_MG(i+1,nlev+1) -xc_MG(ii,nlev)) !**2
        delta = d_1 +d_2

        pro_wht3(i) = d_1/delta
        pro_wht4(i) = d_2/delta
    END DO ! i
    
    DO j = 1, mx1-1 !mgrid_I(nlev+1)
        jj  = 2*j
        jjm = 2*j -1

        d_1 = oned/(yc_MG(j  ,nlev+1) -yc_MG(jjm,nlev)) !**2
        d_2 = oned/(yc_MG(jjm,nlev) -yc_MG(j-1,nlev+1)) !**2
        delta = d_1 +d_2

        pro_wht5(j) = d_1/delta
        pro_wht6(j) = d_2/delta

        d_1 = oned/(yc_MG(jj,nlev) -yc_MG(j  ,nlev+1)) !**2
        d_2 = oned/(yc_MG(j+1,nlev+1) -yc_MG(jj,nlev)) !**2
        delta = d_1 +d_2

        pro_wht7(j) = d_1/delta
        pro_wht8(j) = d_2/delta
    END DO ! i

    DO k=1,mz-1
    DO j=1,mx1-1

        jj  = 2*j
        jjm = 2*j - 1

      DO i = 1, mx1-1 !mgrid_I(nlev+1)-1
        ii  = 2*i
        iim = 2*i - 1

!         prolongation
!         ------------
        var1 = correc(i  ,j  ,k)
        var2 = correc(i-1,j  ,k)
        var3 = correc(i  ,j-1,k)
        var4 = correc(i-1,j-1,k)

        delta1 = pro_wht1(i) * pro_wht5(j) * diblank_MG(i  ,j  ,k)  
        delta2 = pro_wht2(i) * pro_wht5(j) * diblank_MG(i-1,j  ,k)
        delta3 = pro_wht1(i) * pro_wht6(j) * diblank_MG(i  ,j-1,k)  
        delta4 = pro_wht2(i) * pro_wht6(j) * diblank_MG(i-1,j-1,k)
        delta  = delta1 + delta2 + delta3 + delta4

        phi(iim,jjm,k) = delta1*var1 +delta2*var2 + delta3*var3 +delta4*var4 

        IF ( delta > zero ) THEN
          phi(iim,jjm,k)  = phi(iim,jjm,k)/delta
        ENDIF ! delta
!
        var1 = correc(i  ,j  ,k)
        var2 = correc(i+1,j  ,k)
        var3 = correc(i  ,j-1,k)
        var4 = correc(i+1,j-1,k)

        delta1 = pro_wht3(i) * pro_wht5(j) * diblank_MG(i  ,j  ,k)  
        delta2 = pro_wht4(i) * pro_wht5(j) * diblank_MG(i+1,j  ,k)
        delta3 = pro_wht3(i) * pro_wht6(j) * diblank_MG(i  ,j-1,k)  
        delta4 = pro_wht4(i) * pro_wht6(j) * diblank_MG(i+1,j-1,k)
        delta  = delta1 + delta2 + delta3 + delta4

        phi(ii,jjm,k) = delta1*var1 +delta2*var2 + delta3*var3 +delta4*var4 

        IF ( delta > zero ) THEN
          phi(ii,jjm,k) =  phi(ii,jjm,k)/delta
        ENDIF ! delta
!
        var1 = correc(i  ,j  ,k)
        var2 = correc(i-1,j  ,k)
        var3 = correc(i  ,j+1,k)
        var4 = correc(i-1,j+1,k)

        delta1 = pro_wht1(i) * pro_wht7(j) * diblank_MG(i  ,j  ,k)  
        delta2 = pro_wht2(i) * pro_wht7(j) * diblank_MG(i-1,j  ,k)
        delta3 = pro_wht1(i) * pro_wht8(j) * diblank_MG(i  ,j+1,k)  
        delta4 = pro_wht2(i) * pro_wht8(j) * diblank_MG(i-1,j+1,k)
        delta  = delta1 + delta2 + delta3 + delta4

        phi(iim,jj ,k) = delta1*var1 +delta2*var2 + delta3*var3 +delta4*var4 

        IF ( delta > zero ) THEN
          phi(iim,jj ,k)  = phi(iim,jj,k)/delta
        ENDIF ! delta
!
        var1 = correc(i  ,j  ,k)
        var2 = correc(i+1,j  ,k)
        var3 = correc(i  ,j+1,k)
        var4 = correc(i+1,j+1,k)

        delta1 = pro_wht3(i) * pro_wht7(j) * diblank_MG(i  ,j  ,k)  
        delta2 = pro_wht4(i) * pro_wht7(j) * diblank_MG(i+1,j  ,k)
        delta3 = pro_wht3(i) * pro_wht8(j) * diblank_MG(i  ,j+1,k)  
        delta4 = pro_wht4(i) * pro_wht8(j) * diblank_MG(i+1,j+1,k)
        delta  = delta1 + delta2 + delta3 + delta4

        phi(ii ,jj ,k) = delta1*var1 +delta2*var2 + delta3*var3 +delta4*var4 

        IF ( delta > zero ) THEN
          phi(ii ,jj ,k)  = phi(ii ,jj ,k)/delta
        ENDIF ! delta

      END DO ! i
 
    END DO ! j
    END DO ! k
    
END SUBROUTINE MG_Prolong_FMG


SUBROUTINE MG_Prolong_X(phi,correc,nlev, mx1, my, mz, mx)
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, mx1
 
    REAL(KIND=CGREAL), DIMENSION(0:mx1,0:my,0:mz), INTENT (IN)  :: correc
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (OUT) :: phi

    INTEGER           :: ii, iim, jj, jjm, kk, kkm, i, j, k, nxyzMax, iErr
    REAL(KIND=CGREAL) :: var1, var2, var3, delta1, delta2, delta, d_1, d_2
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: pro_wht1, pro_wht2, pro_wht3, pro_wht4
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: diblank_MG

!   Allocate local arrays
!   ---------------------
    ALLOCATE(diblank_MG(0:mx,0:my,0:mz),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_X: Memory Allocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

    nxyzMax = MAX(nx,ny,nz)
    
    ALLOCATE(pro_wht1(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_X: Memory Allocation Error for pro_wht1'
      STOP
    ENDIF ! iErr    
    
    ALLOCATE(pro_wht2(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_X: Memory Allocation Error for pro_wht2'
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(pro_wht3(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_X: Memory Allocation Error for pro_wht3'
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(pro_wht4(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_X: Memory Allocation Error for pro_wht4'
      STOP
    ENDIF ! iErr 

!   Start Prolongation
!   ------------------
!   Prolongate in x-direction
!   -------------------------
    DO k=1, mz-1
    DO j=1, my-1
    DO i=1, mx1-1 !mgrid_I(nlev+1) 
        diblank_MG(i,j,k) = oned-REAL(iblank_MG(i,j,k),KIND=CGREAL)
    END DO ! i
    END DO ! j
    END DO ! k

    DO i = 1, mx1-1 !mgrid_I(nlev+1)
        ii  = 2*i
        iim = 2*i -1

        d_1 = oned/(xc_MG(iim,nlev) -xc_MG(i  ,nlev+1))**2
        d_2 = oned/(xc_MG(iim,nlev) -xc_MG(i-1,nlev+1))**2
        delta = d_1 +d_2

        pro_wht1(i) = d_1/delta
        pro_wht2(i) = d_2/delta

        d_1 = oned/(xc_MG(ii,nlev) -xc_MG(i  ,nlev+1))**2
        d_2 = oned/(xc_MG(ii,nlev) -xc_MG(i+1,nlev+1))**2
        delta = d_1 +d_2

        pro_wht3(i) = d_1/delta
        pro_wht4(i) = d_2/delta
    END DO ! i

    DO k=1,mz-1
    DO j=1,my-1
        phi(1, j, k) = correc(1, j, k)

        delta1  = pro_wht3(1) *diblank_MG(1,j,k)
        delta2  = pro_wht4(1) *diblank_MG(2,j,k)
        delta   = delta1 + delta2

        phi(2,j,k) = delta1 *correc(1,j,k) &
                   + delta2 *correc(2,j,k)
         
        IF ( delta > zero ) THEN
          phi(2,j,k) = phi(2,j,k)/delta
        ENDIF ! delta

        DO i = 2, mx1-2 !mgrid_I(nlev+1)-1
          ii  = 2*i
          iim = 2*i - 1

          var1 = correc(i,j,k)
          var2 = correc(i-1,j,k)
          var3 = correc(i+1,j,k)

!         prolongation
!         ------------
          delta1 = pro_wht1(i) *diblank_MG(i  ,j,k)  
          delta2 = pro_wht2(i) *diblank_MG(i-1,j,k)
          delta  = delta1 +delta2

          phi(iim,j,k) = delta1*var1 +delta2*var2

          IF ( delta > zero ) THEN
            phi(iim,j,k)  = phi(iim,j,k)/delta
          ENDIF ! delta

          delta1 = pro_wht3(i) *diblank_MG(i  ,j,k)    
          delta2 = pro_wht4(i) *diblank_MG(i+1,j,k)
          delta  = delta1 +delta2

          phi(ii,j,k) = delta1*var1 +delta2*var3

          IF ( delta > zero ) THEN
            phi(ii,j,k) =  phi(ii,j,k)/delta
          ENDIF ! delta
!        IF (iim==85 .and. j==473 .and. nlev==1) THEN
!           PRINT *, I,J,K
!           PRINT *, IIm, II, NLEV
!           PRINT *, dx_MG(iim,nlev), dx_MG(ii,nlev)
!           PRINT *, iblank_MG(i-1,j,k), Rblank_MG(i-1,j,k)
!           PRINT *, iblank_MG(i  ,j,k), Rblank_MG(i  ,j,k)
!           PRINT *, iblank_MG(i+1,j,k), Rblank_MG(i+1,j,k)
!           PRINT *, dx_MG(i-1,nlev+1), dx_MG(i,nlev+1), dx_MG(i+1,nlev+1)
!           PRINT *, pro_wht1(i), pro_wht2(i), pro_wht3(i), pro_wht4(i)    
!           PRINT *, VAR1, var2, var3
!           PRINT *, phi(iim,j,k), phi(ii,j,k)
!           STOP
!        ENDIF
        END DO ! i
 
        i   = mx1-1 !mgrid_I(nlev+1)
        ii  = 2*i
        iim = 2*i -1

        delta1 = pro_wht1(i) *diblank_MG(i  ,j,k)
        delta2 = pro_wht2(i) *diblank_MG(i-1,j,k)
        delta  = delta1 +delta2

        phi(iim,j,k) = delta1 *correc(i  ,j,k) &
                     + delta2 *correc(i-1,j,k)

        IF ( delta > zero ) THEN
          phi(iim,j,k) = phi(iim,j,k)/delta 
        ENDIF ! delta

        phi(ii,j,k) = correc(i,j,k)

!       new modification -H. Dong
!       -------------------------
        IF ( nlev == 1 ) THEN
          IF ( mx-1 > 2*i ) THEN
            phi(mx-1,j,k) = correc(i,j,k)
          ENDIF ! nxc
        ENDIF ! nlev
!       ------------------------------
!       end new modification - H. Dong

    END DO ! j
    END DO ! k

!   Deallocate local array
!   ----------------------
    DEALLOCATE(diblank_MG,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong: Memory Deallocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

    DEALLOCATE(pro_wht1,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong: Memory Deallocation Error for pro_wht1'
      STOP
    ENDIF ! iErr    
    
    DEALLOCATE(pro_wht2,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong: Memory Deallocation Error for pro_wht2'
      STOP
    ENDIF ! iErr 
    
    DEALLOCATE(pro_wht3,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong: Memory Deallocation Error for pro_wht3'
      STOP
    ENDIF ! iErr 

    DEALLOCATE(pro_wht4,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong: Memory Dellocation Error for pro_wht4'
      STOP
    ENDIF ! iErr 
    
END SUBROUTINE MG_Prolong_X

SUBROUTINE MG_Prolong_Y(phi,correc,nlev, mx, my1, mz, my)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, my1
 
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my1,0:mz), INTENT (IN)  :: correc
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (OUT) :: phi

    INTEGER           :: ii, iim, jj, jjm, kk, kkm, i, j, k, nxyzMax, iErr
    REAL(KIND=CGREAL) :: var1, var2, var3, delta1, delta2, delta, d_1, d_2
    REAL(KIND=CGREAL), DIMENSION(:)    , ALLOCATABLE :: pro_wht1, pro_wht2, pro_wht3, pro_wht4
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: diblank_MG

!   Allocate local arrays
!   ---------------------
    ALLOCATE(diblank_MG(0:mx,0:my,0:mz),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Allocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

    nxyzMax = MAX(nx,ny,nz)
    
    ALLOCATE(pro_wht1(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Allocation Error for pro_wht1'
      STOP
    ENDIF ! iErr    
    
    ALLOCATE(pro_wht2(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Allocation Error for pro_wht2'
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(pro_wht3(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Allocation Error for pro_wht3'
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(pro_wht4(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Allocation Error for pro_wht4'
      STOP
    ENDIF ! iErr 

!   Start Prolongation
!   ------------------
!   Prolongate in y-direction
!   -------------------------
    DO k=1,mz-1
    DO j=1,my1-1 !mgrid_J(nlev+1) 
    DO i=1,mx-1
        diblank_MG(i,j,k) = oned-REAL(iblank_MG(i,j,k),KIND=CGREAL)
    END DO ! i
    END DO ! j
    END DO ! k

    DO j = 1, my1-1 !mgrid_J(nlev+1) 
        jj  = 2*j 
        jjm = 2*j -1 

        d_1   = oned/(yc_MG(jjm,nlev) -yc_MG(j  ,nlev+1))**2 
        d_2   = oned/(yc_MG(jjm,nlev) -yc_MG(j-1,nlev+1))**2
        delta = d_1 + d_2  

        pro_wht1(j) = d_1/delta 
        pro_wht2(j) = d_2/delta

        d_1   = oned/(yc_MG(jj,nlev) -yc_MG(j  ,nlev+1))**2
        d_2   = oned/(yc_MG(jj,nlev) -yc_MG(j+1,nlev+1))**2
        delta = d_1 + d_2

        pro_wht3(j) = d_1/delta
        pro_wht4(j) = d_2/delta 
    END DO ! j
 
    DO k = 1, mz-1
    DO i = 1, mx-1
        phi(i,1, k) = correc(i, 1, k)

        delta1 =  pro_wht3(1) *diblank_MG(i,1,k)
        delta2 =  pro_wht4(1) *diblank_MG(i,2,k)
        delta  = delta1 + delta2

        phi(i,2,k) = delta1 *correc(i,1,k) &
                   + delta2 *correc(i,2,k)

        IF ( delta > zero ) THEN
          phi(i,2,k) = phi(i,2,k)/delta
        ENDIF ! delta

        DO j = 2, my1-2 !mgrid_J(nlev+1)-1
          jj = 2*j
          jjm = 2*j - 1

          var1 = correc(i,j  ,k)
          var2 = correc(i,j-1,k)
          var3 = correc(i,j+1,k)

!         prolongation
!         ------------
          delta1 = pro_wht1(j) *diblank_MG(i,j  ,k)
          delta2 = pro_wht2(j) *diblank_MG(i,j-1,k)
          delta  = delta1 +delta2

          phi(i,jjm,k) = delta1*var1 +delta2*var2

          IF ( delta > zero ) THEN
            phi(i,jjm,k) = phi(i,jjm,k)/delta
          ENDIF ! delta

          delta1 = pro_wht3(j) *diblank_MG(i,j  ,k)
          delta2 = pro_wht4(j) *diblank_MG(i,j+1,k)
          delta  = delta1 +delta2

          phi(i,jj,k) = delta1*var1 +delta2*var3

          IF ( delta > zero ) THEN
            phi(i,jj,k) =  phi(i,jj,k)/delta
          ENDIF ! delta
        END DO ! j

        j   = my1-1 !mgrid_J(nlev+1)
        jj  = 2*j
        jjm = 2*j -1

        delta1 = pro_wht1(j) *diblank_MG(i,j  ,k)
        delta2 = pro_wht2(j) *diblank_MG(i,j-1,k)
        delta  = delta1 +delta2

        phi(i,jjm,k) = delta1 *correc(i,j  ,k) &
                     + delta2 *correc(i,j-1,k)

        IF ( delta > zero ) THEN
          phi(i,jjm,k) = phi(i,jjm,k)/delta  
        ENDIF

        phi(i,jj,k) = correc(i,j,k) 

!       new modification - H. Dong
!       --------------------------
        IF ( nlev == 1 .and. my-1 > 2*j ) THEN
          phi(i,my-1,k) = correc(i,j,k)
        ENDIF ! nLev
!       --------------------
!       end new modification

    END DO ! i  
    END DO ! k

!   Deallocate local array
!   ----------------------
    DEALLOCATE(diblank_MG,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Deallocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

    DEALLOCATE(pro_wht1,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Deallocation Error for pro_wht1'
      STOP
    ENDIF ! iErr    
    
    DEALLOCATE(pro_wht2,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Deallocation Error for pro_wht2'
      STOP
    ENDIF ! iErr 
    
    DEALLOCATE(pro_wht3,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Deallocation Error for pro_wht3'
      STOP
    ENDIF ! iErr 

    DEALLOCATE(pro_wht4,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Y: Memory Dellocation Error for pro_wht4'
      STOP
    ENDIF ! iErr 
    
END SUBROUTINE MG_Prolong_Y
!---------------------------------------------------------------------

SUBROUTINE MG_Prolong_Z(phi,correc,nlev, mx, my, mz1, mz)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev,  mx, my, mz, mz1
 
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz1), INTENT (IN)  :: correc
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (OUT) :: phi

    INTEGER           :: ii, iim, jj, jjm, kk, kkm, i, j, k, nxyzMax, iErr
    REAL(KIND=CGREAL) :: var1, var2, var3, delta1, delta2, delta, d_1, d_2
    REAL(KIND=CGREAL), DIMENSION(:)    , ALLOCATABLE :: pro_wht1, pro_wht2, pro_wht3, pro_wht4
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: diblank_MG
     
!   Allocate local arrays
!   ---------------------
    ALLOCATE(diblank_MG(0:mx,0:my,0:mz),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Allocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

    nxyzMax = MAX(nx,ny,nz)
    
    ALLOCATE(pro_wht1(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Allocation Error for pro_wht1'
      STOP
    ENDIF ! iErr    
    
    ALLOCATE(pro_wht2(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Allocation Error for pro_wht2'
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(pro_wht3(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Allocation Error for pro_wht3'
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(pro_wht4(0:nxyzMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Allocation Error for pro_wht4'
      STOP
    ENDIF ! iErr 

!   Start Prolongation
!   ------------------
!   Prolongate in z-direction
!   -------------------------
    DO k=1, mz1-1 !mgrid_K(nlev+1)
    DO j=1, my-1
    DO i=1, mx-1
        diblank_MG(i,j,k) = oned-REAL(iblank_MG(i,j,k),KIND=CGREAL)
    END DO
    END DO
    END DO

    DO k = 1, mz1-1 !mgrid_K(nlev+1) 
        kk  = 2*k
        kkm = 2*k -1 

        d_1   = oned/(zc_MG(kkm,nlev) -zc_MG(k  ,nlev+1))**2 
        d_2   = oned/(zc_MG(kkm,nlev) -zc_MG(k-1,nlev+1))**2
        delta = d_1 + d_2  

        pro_wht1(k) = d_1/delta 
        pro_wht2(k) = d_2/delta

        d_1 = oned/(zc_MG(kk,nlev) -zc_MG(k  ,nlev+1))**2
        d_2 = oned/(zc_MG(kk,nlev) -zc_MG(k+1,nlev+1))**2
        delta = d_1 +d_2

        pro_wht3(k) = d_1/delta
        pro_wht4(k) = d_2/delta 

    END DO ! k

    DO j = 1, my-1
    DO i = 1, mx-1
        phi(i,j,1) = correc(i, j, 1)

        delta1 = pro_wht3(1) *diblank_MG(i,j,1)
        delta2 = pro_wht4(1) *diblank_MG(i,j,2)
        delta  = delta1 +delta2

        phi(i,j,2) = delta1 *correc(i,j,1) &
                   + delta2 *correc(i,j,2)

        IF ( delta > zero ) THEN
          phi(i,j,2) = phi(i,j,2)/delta
        ENDIF ! delta

        DO k = 2, mz1-2 !mgrid_K(nlev+1)-1
          kk  = 2*k
          kkm = 2*k -1

          var1 = correc(i,j,k)
          var2 = correc(i,j,k-1)
          var3 = correc(i,j,k+1)

!         prolongation
!         ------------
          delta1 = pro_wht1(k) *diblank_MG(i,j,k  )
          delta2 = pro_wht2(k) *diblank_MG(i,j,k-1)
          delta  = delta1 +delta2
 
          phi(i,j,kkm) = delta1*var1 +delta2*var2

          IF ( delta > zero ) THEN
            phi(i,j,kkm) = phi(i,j,kkm)/delta
          ENDIF ! delta
             
          delta1 = pro_wht3(k) *diblank_MG(i,j,k  )
          delta2 = pro_wht4(k) *diblank_MG(i,j,k+1)
          delta  = delta1 +delta2

          phi(i,j,kk) = delta1*var1 +delta2*var3
 
          IF ( delta > zero ) THEN
            phi(i,j,kk) =  phi(i,j,kk)/delta
          ENDIF ! delta
        END DO ! k
 
        k  = mz1-1 !mgrid_K(nlev+1)
        kk  = 2*k
        kkm = 2*k -1

        delta1 = pro_wht1(k) *diblank_MG(i,j,k  )
        delta2 = pro_wht2(k) *diblank_MG(i,j,k-1)
        delta  = delta1 +delta2

        phi(i,j,kkm) = delta1*correc(i,j,k  ) &
                     + delta2*correc(i,j,k-1)

        IF ( delta > zero ) THEN
          phi(i,j,kkm) = phi(i,j,kkm)/delta  
        ENDIF ! delta

        phi(i,j,kk) = correc(i,j,k)  

!       new modification -H. Dong
!       -------------------------
        IF ( nlev == 1 .and. mz-1 > 2*k ) THEN
          phi(i,j,mz-1) = correc(i,j,k)
        ENDIF ! nLev
!       --------------------
!       end new modification
              
    END DO ! i  
    END DO ! j

!   Deallocate local array
!   ----------------------
    DEALLOCATE(diblank_MG,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Deallocation Error for diblank_MG'
      STOP
    ENDIF ! iErr

    DEALLOCATE(pro_wht1,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Deallocation Error for pro_wht1'
      STOP
    ENDIF ! iErr    
    
    DEALLOCATE(pro_wht2,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Deallocation Error for pro_wht2'
      STOP
    ENDIF ! iErr 
    
    DEALLOCATE(pro_wht3,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Deallocation Error for pro_wht3'
      STOP
    ENDIF ! iErr 

    DEALLOCATE(pro_wht4,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Prolong_Z: Memory Dellocation Error for pro_wht4'
      STOP
    ENDIF ! iErr 
    
END SUBROUTINE MG_Prolong_Z
!---------------------------------------------------------------------

SUBROUTINE MG_Allocate_Memory_Grid()
           
    USE global_parameters
    USE flow_parameters
    USE flow_arrays 
    USE grid_arrays 
    USE boundary_arrays
    USE multiuse_arrays
    USE MG_parameters 
    USE MG_arrays 
 
    IMPLICIT NONE 
     
    INTEGER :: iErr
    INTEGER :: I
 
!   Allocate MG arrays 
!   ------------------
    ALLOCATE(mgrid_I(mgLevels_X),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for mgrid_I'
      STOP 
    ENDIF ! iErr 
     
    ALLOCATE(mgrid_J(mgLevels_Y),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for mgrid_J'
      STOP 
    ENDIF ! iErr
     
    ALLOCATE(mgrid_K(mgLevels_Z),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for mgrid_K'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(dxcinv_MG(0:nx+1, mgLevels_X),dxc_MG(0:nx+1, mgLevels_X),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for dxcinv_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(dxinv_MG(0:nx+1, mgLevels_X),dx_MG(0:nx+1, mgLevels_X),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for dxinv_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(dycinv_MG(0:ny+1, mgLevels_Y),dyc_MG(0:ny+1, mgLevels_Y),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for dycinv_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(dyinv_MG(0:ny+1, mgLevels_Y),dy_MG(0:ny+1, mgLevels_Y),STAT=iErr) 
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for dyinv_MG' 
      STOP
    ENDIF ! iErr
 
    ALLOCATE(dzcinv_MG(0:nz+1, mgLevels_Z),dzc_MG(0:nz+1, mgLevels_Z),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for dzcinv_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(dzinv_MG(0:nz+1, mgLevels_Z),dz_MG(0:nz+1, mgLevels_Z),STAT=iErr) 
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for dzinv_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(x_MG(0:nx+1, mgLevels_X),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for x_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(xc_MG(0:nx+1, mgLevels_X),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for xc_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(y_MG(0:ny+1, mgLevels_Y),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for y_MG' 
      STOP
    ENDIF ! iErr
 
    ALLOCATE(yc_MG(0:ny+1, mgLevels_Y),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for yc_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(z_MG(0:nz+1, mgLevels_Z),STAT=iErr) 
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for z_MG' 
      STOP
    ENDIF ! iErr
 
    ALLOCATE(zc_MG(0:nz+1, mgLevels_Z),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for zc_MG'
      STOP
    ENDIF ! iErr
    
    IF ( pp_solver_type == PP_SOLVER_TYPE_MG .or. pp_solver_type == PP_SOLVER_TYPE_MG_Point_Jacobi) &
        CALL MG_Allocate_Memory_IUP

!     ALLOCATE( pgradx1_MG(0:ny+1,0:nz+1), pgradx2_MG(0:ny+1,0:nz+1))
!     ALLOCATE( pgrady1_MG(0:nx+1,0:nz+1), pgrady2_MG(0:nx+1,0:nz+1))
!     ALLOCATE( pgradz1_MG(0:nx+1,0:ny+1), pgradz2_MG(0:nx+1,0:ny+1))
    ALLOCATE(ghostcellMark_MG(0:nx+1,0:ny+1,0:nz+1))        !!! Added by ZX Liang
        
END SUBROUTINE MG_Allocate_Memory_Grid

!---------------------------------------------------------------------
SUBROUTINE MG_Allocate_Memory_IUP()
           
    USE global_parameters
    USE flow_parameters
    USE MG_arrays 
 
    IMPLICIT NONE 
     
    INTEGER :: iErr


    ALLOCATE(ium_MG(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for ium_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(iup_MG(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for iup_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(jum_MG(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for jum_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(jup_MG(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for jup_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(kum_MG(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for kum_MG'
      STOP
    ENDIF ! iErr
 
    ALLOCATE(kup_MG(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'MG_Allocate_Memory_Grid: Memory Allocation Error for kup_MG'
      STOP
    ENDIF ! iErr

END SUBROUTINE MG_Allocate_Memory_IUP
!---------------------------------------------------------------------

SUBROUTINE MG_itsolv_Point_Jacobi(var, r, nLevX, nLevY, nLevZ, mx,my,mz) 
 
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE
 
!... parameters
 
    INTEGER, INTENT(IN) :: nLevX, nLevY, nLevZ, mx, my, mz
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)     ::r
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (INOUT)  ::var
 
 
!... Local variables
 
    INTEGER :: i,j,k
    INTEGER :: iBody, iRow, iG,jG, iNode, jNode, n

    REAL(KIND=CGREAL) 	:: bmx,bpx,bcx,bc
    REAL(KIND=CGREAL) 	:: bmy,bpy,bcy
    REAL(KIND=CGREAL) 	:: bmz,bpz,bcz, rhss
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz) ::var1
 
    var1=var
	
    DO k=1, mz-1
    DO j=1, my-1
      DO i=1, mx-1

        bmx =   dxcinv_mg(i,nLevX)  *dxinv_mg(i,nLevX)*(oned - ium_mg(i,j,k) )
        bpx =   dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*(oned - iup_mg(i,j,k) )
        bcx = - ( bmx + bpx )

        bmy =   dycinv_mg(j,nLevY)  *dyinv_mg(j,nLevY)*(oned - jum_mg(i,j,k) )
        bpy =   dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*(oned - jup_mg(i,j,k) )
        bcy = - ( bmy + bpy )

        bmz =   dzcinv_mg(k,nLevZ)  *dzinv_mg(k,nLevZ)*(oned - kum_mg(i,j,k) ) &
                                  *REAL((ndim - DIM_2D),KIND=CGREAL)
        bpz =   dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*(oned - kup_mg(i,j,k) ) &
                                  *REAL((ndim - DIM_2D),KIND=CGREAL)
        bcz = - ( bmz + bpz )

        bmx = bmx*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
        bpx = bpx*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )

        bmy = bmy*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
        bpy = bpy*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )

        bmz = bmz*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
        bpz = bpz*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )

        bc =   (bcx+bcy+bcz)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )  &
               + REAL(iblank_MG(i,j,k),KIND=CGREAL)

        rhss     = r(i,j,k)           &
                 - var(i-1,j,k)*bmx                &
                 - var(i+1,j,k)*bpx                &
                 - var(i,j-1,k)*bmy                &
                 - var(i,j+1,k)*bpy                &
                 - var(i,j,k-1)*bmz                &
                 - var(i,j,k+1)*bpz                &
                + nlw(i,j,k)
             
       rhss = (rhss/bc)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) ) &
                + REAL(ghostcellMark_mg(i,j,k),KIND=CGREAL)*var(i,j,k)
 
       var1(i,j,k) = var(i,j,k) + omega*(rhss - var(i,j,k))
      ENDDO
 
    ENDDO ! j
    ENDDO ! k
    
    var=var1
	
END SUBROUTINE MG_itsolv_Point_Jacobi
!-------------------------------------------------------------------------------

SUBROUTINE write_dump_debug_MG(vname, dbg,var, nLevX, nLevY, nLevZ, MX, MY, MZ, ILO,IHI, JLO,JHI, KLO,KHI)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays
  USE MG_arrays

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nLevX, nLevY, nLevZ, MX, MY, MZ, ILO,IHI, JLO,JHI, KLO,KHI
  INTEGER :: i, j, k, iG, jG, dbg
  REAL(KIND=CGREAL) :: VAR(ILO:IHI, JLO:JHI, KLO:KHI)
!  REAL(KIND=CGREAL) :: VAR(:, :, :)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
!  WRITE(fname1,"(A,'.000.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
!  WRITE(70,*)'VARIABLES="X","Y","U"'
!  WRITE(70,*) 'ZONE F=POINT, I=',(ihi-ilo-1)/2,', J=',(jhi-jlo-1)/2
!  WRITE(fname1,"(A,'.001.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=71,FILE=fname1,FORM='FORMATTED')
!  WRITE(71,*)'VARIABLES="X","Y","U"'
!  WRITE(71,*) 'ZONE F=POINT, I=',(ihi-ilo-1)/2,', J=',(jhi-jlo-1)/2
!  WRITE(fname1,"(A,'.002.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=72,FILE=fname1,FORM='FORMATTED')
!  WRITE(72,*)'VARIABLES="X","Y","U"'
!  WRITE(72,*) 'ZONE F=POINT, I=',(ihi-ilo-1)/2,', J=',(jhi-jlo-1)/2
!  WRITE(fname1,"(A,'.003.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=73,FILE=fname1,FORM='FORMATTED')
!  WRITE(73,*)'VARIABLES="X","Y","U"'
!  WRITE(73,*) 'ZONE F=POINT, I=',(ihi-ilo-1)/2,', J=',(jhi-jlo-1)/2

  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'_',I2.2,'_',I2.2,'_',I2.2,'.dat')") trim(vname),ntime,dbg,nlevx,nlevy,nlevz
  OPEN(UNIT=177,FILE=fname1,FORM='FORMATTED')
  WRITE(177,*)'VARIABLES="X","Y","U"'

  do k=1,mz
  WRITE(177,*) 'ZONE F=POINT, I=',MX,', J=',MY
  DO j=1,MY
  DO i=1,MX
!        WRITE(70,'(2(3X,1PE12.5),(3X,1PE19.11))') xc_MG(i,nLevX),yc_MG(j,nLevY), var(i,j,k)
        WRITE(177,'(2(3X,I8),(3X,1PE22.15))') i,j, var(i,j,k)
!      IF (j<=(jhi-jlo-1)/2) THEN
!        IF (i<=(ihi-ilo-1)/2) THEN
!        WRITE(70,'(2(3X,1PE12.5),(3X,1PE19.11))') xc_MG(i,nLevX), yc_MG(j,nLevY), var(i,j,k)
!        ELSE
!        WRITE(72,'(2(3X,1PE12.5),(3X,1PE19.11))') xc_MG(i,nLevX), yc_MG(j,nLevY), var(i,j,k)
!        ENDIF
!      ELSE
!        IF (i<=(ihi-ilo-1)/2) THEN
!        WRITE(71,'(2(3X,1PE12.5),(3X,1PE19.11))') xc_MG(i,nLevX), yc_MG(j,nLevY), var(i,j,k)
!        ELSE
!        WRITE(73,'(2(3X,1PE12.5),(3X,1PE19.11))') xc_MG(i,nLevX), yc_MG(j,nLevY), var(i,j,k)
!        ENDIF
!      ENDIF
  END DO
  END DO
  END DO

  CLOSE(177)
!  CLOSE(71)
!  CLOSE(72)
!  CLOSE(73)

  RETURN
END SUBROUTINE

!---------------------------------------------------------------------
SUBROUTINE write_dump_debug_i_MG(vname, dbg,var, nLevX, nLevY, nLevZ, ILO, IHI, JLO, JHI, KLO, KHI)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays
  USE MG_arrays

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nLevX, nLevY, nLevZ, ILO, IHI, JLO, JHI, KLO, KHI
  INTEGER :: i, j, k, iG, jG, dbg
!  REAL(KIND=CGREAL) :: VAR(ILO:IHI, JLO:JHI, KLO:KHI)
  INTEGER :: VAR(ILO:IHI, JLO:JHI, KLO:KHI)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
!  WRITE(fname1,"(A,'.000.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
!  WRITE(70,*)'VARIABLES="X","Y","U"'
!  WRITE(70,*) 'ZONE F=POINT, I=',(ihi-ilo-1)/2,', J=',(jhi-jlo-1)/2
!  WRITE(fname1,"(A,'.001.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=71,FILE=fname1,FORM='FORMATTED')
!  WRITE(71,*)'VARIABLES="X","Y","U"'
!  WRITE(71,*) 'ZONE F=POINT, I=',(ihi-ilo-1)/2,', J=',(jhi-jlo-1)/2
!  WRITE(fname1,"(A,'.002.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=72,FILE=fname1,FORM='FORMATTED')
!  WRITE(72,*)'VARIABLES="X","Y","U"'
!  WRITE(72,*) 'ZONE F=POINT, I=',(ihi-ilo-1)/2,', J=',(jhi-jlo-1)/2
!  WRITE(fname1,"(A,'.003.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=73,FILE=fname1,FORM='FORMATTED')
!  WRITE(73,*)'VARIABLES="X","Y","U"'
!  WRITE(73,*) 'ZONE F=POINT, I=',(ihi-ilo-1)/2,', J=',(jhi-jlo-1)/2

  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'_',I2.2,'_',I2.2,'_',I2.2,'.dat')") trim(vname),ntime,dbg,nlevx,nlevy,nlevz
  OPEN(UNIT=177,FILE=fname1,FORM='FORMATTED')
  WRITE(177,*)'VARIABLES="X","Y","U"'
  WRITE(177,*) 'ZONE F=POINT, I=',(ihi-ilo-1),', J=',(jhi-jlo-1)

  k=1
  DO j=1,jhi-jlo-1
  DO i=1,ihi-ilo-1
        WRITE(177,'(2(3X,1PE12.5),(3X,I8))') xc_MG(i,nLevX), yc_MG(j,nLevY), (var(i,j,k))
!      IF (j<=(jhi-jlo-1)/2) THEN
!        IF (i<=(ihi-ilo-1)/2) THEN
!        WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc_MG(i,nLevX), yc_MG(j,nLevY), int(var(i,j,k))
!        ELSE
!        WRITE(71,'(2(3X,1PE12.5),(3X,I8))') xc_MG(i,nLevX), yc_MG(j,nLevY), int(var(i,j,k))
!        ENDIF
!      ELSE
!        IF (i<=(ihi-ilo-1)/2) THEN
!        WRITE(72,'(2(3X,1PE12.5),(3X,I8))') xc_MG(i,nLevX), yc_MG(j,nLevY), int(var(i,j,k))
!        ELSE
!        WRITE(73,'(2(3X,1PE12.5),(3X,I8))') xc_MG(i,nLevX), yc_MG(j,nLevY), int(var(i,j,k))
!        ENDIF
!      ENDIF
  END DO
  END DO

  CLOSE(177)
!  CLOSE(71)
!  CLOSE(72)
!  CLOSE(73)

  RETURN
END SUBROUTINE


SUBROUTINE MG_Prolong_X_Conformal_Mapping(phi,correc,nlev, mx1, my, mz, mx)
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, mx1
 
    REAL(KIND=CGREAL), DIMENSION(0:mx1,0:my,0:mz), INTENT (IN)  :: correc
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (OUT) :: phi

    INTEGER           :: ii, iim, jj, jjm, kk, kkm, i, j, k, nxyzMax, iErr
    REAL(KIND=CGREAL) :: var1, var2, var3, delta1, delta2, delta, d_1, d_2, BC2, BA1, BA2
    REAL(KIND=CGREAL), DIMENSION(0:nx+1) :: pro_wht1, pro_wht2, pro_wht3, pro_wht4
    INTEGER(1), DIMENSION(:,:,:), POINTER :: iblank_MG1

    IF (nlev > 1) THEN
        iblank_MG1 => MGX(nlev)%iblank
    ELSE
        iblank_MG1 => iblank
    ENDIF


!   Start Prolongation
!   ------------------
!   Prolongate in x-direction
!   -------------------------

    DO i = 1, mx1-1 !mgrid_I(nlev+1)
        ii  = 2*i
        iim = 2*i -1

        d_1 = oned/(xc_MG(iim,nlev) -xc_MG(i  ,nlev+1))**2
        d_2 = oned/(xc_MG(iim,nlev) -xc_MG(i-1,nlev+1))**2
        delta = d_1 +d_2

        pro_wht1(i) = d_1/delta
        pro_wht2(i) = d_2/delta

        d_1 = oned/(xc_MG(ii,nlev) -xc_MG(i  ,nlev+1))**2
        d_2 = oned/(xc_MG(ii,nlev) -xc_MG(i+1,nlev+1))**2
        delta = d_1 +d_2

        pro_wht3(i) = d_1/delta
        pro_wht4(i) = d_2/delta
    END DO ! i

    DO k=1,mz-1
    DO j=1,my-1
        phi(1, j, k) = correc(1, j, k)

        delta1  = pro_wht3(1) * (1-iblank_MG(1,j,k))
        delta2  = pro_wht4(1) * (1-iblank_MG(2,j,k))
        
        phi(2,j,k) = delta1 *correc(1,j,k) &
                   + delta2 *correc(2,j,k)
         
        delta   = delta1 + delta2
        IF ( delta > zero ) THEN
          phi(2,j,k) = phi(2,j,k)/delta
        ENDIF ! delta

        DO i = 2, mx1-2 !mgrid_I(nlev+1)-1
          ii  = 2*i
          iim = 2*i - 1

          var1 = correc(i,j,k)
          var2 = correc(i-1,j,k)
          var3 = correc(i+1,j,k)

!         prolongation
!         ------------
          delta1 = pro_wht1(i) * (1-iblank_MG(i  ,j,k))
          delta2 = pro_wht2(i) * (1-iblank_MG(i-1,j,k))
          phi(iim,j,k) = delta1*var1 +delta2*var2

          delta  = delta1 +delta2
          IF ( delta > zero ) THEN
            phi(iim,j,k)  = phi(iim,j,k)/delta
          ENDIF ! delta

          delta1 = pro_wht3(i) * (1-iblank_MG(i  ,j,k))
          delta2 = pro_wht4(i) * (1-iblank_MG(i+1,j,k))
          phi(ii,j,k) = delta1*var1 +delta2*var3

          delta  = delta1 +delta2
          IF ( delta > zero ) THEN
            phi(ii,j,k) =  phi(ii,j,k)/delta
          ENDIF ! delta
          
          IF (iblank_MG1(iim,j,k)==1 .and. iblank_MG1(ii,j,k)==0 ) THEN
            BC2=(half*dx_mg(i+1,nlev+1)+dx_mg(ii,nlev))**2
            BA1=(half*dx_mg(i  ,nlev+1)-dx_mg(ii,nlev))**2
            BA2=(half*dx_mg(ii,nlev))**2
            phi(ii,j,k) = var1+(BA2-BA1)/(BC2-BA1)*(var3-var1)
            ! validation: uniform grid
            ! BC2=(0.5*2L+L)**2=4L2
            ! BA1=(0.5*2L-L)**2=0
            ! BA2=(0.5*L)**2=0.25L2
            !phi(ii)=var1+(0.25L2-0)/(4L2-0)*(var3-var1)=15/16var1+1/16var3
          ENDIF
          IF (iblank_MG1(iim,j,k)==0 .and. iblank_MG1(ii,j,k)==1 ) THEN
            BC2=(half*dx_mg(i-1,nlev+1)+dx_mg(iim,nlev))**2
            BA1=(half*dx_mg(i  ,nlev+1)-dx_mg(iim,nlev))**2
            BA2=(half*dx_mg(iim,nlev))**2
            phi(iim,j,k) = var1+(BA2-BA1)/(BC2-BA1)*(var2-var1)
            ! validation: uniform grid
            ! BC2=(0.5*2L+L)**2=4L2
            ! BA1=(0.5*2L-L)**2=0
            ! BA2=(0.5*L)**2=0.25L2
            !phi(iim)=var1+(0.25L2-0)/(4L2-0)*(var2-var1)=15/16var1+1/16var2
          ENDIF
        END DO ! i
 
        i   = mx1-1 !mgrid_I(nlev+1)
        ii  = 2*i
        iim = 2*i -1

        delta1 = pro_wht1(i) * (1-iblank_MG(i  ,j,k))
        delta2 = pro_wht2(i) * (1-iblank_MG(i-1,j,k))

        phi(iim,j,k) = delta1 *correc(i  ,j,k) &
                     + delta2 *correc(i-1,j,k)
        
        delta  = delta1 +delta2
        IF ( delta > zero ) THEN
          phi(iim,j,k) = phi(iim,j,k)/delta 
        ENDIF ! delta

        phi(ii,j,k) = correc(i,j,k)

!       new modification -H. Dong
!       -------------------------
        IF ( nlev == 1 ) THEN
          IF ( mx-1 > 2*i ) THEN
            phi(mx-1,j,k) = correc(i,j,k)
          ENDIF ! nxc
        ENDIF ! nlev
!       ------------------------------
!       end new modification - H. Dong

    END DO ! j
    END DO ! k

    
END SUBROUTINE MG_Prolong_X_Conformal_Mapping

SUBROUTINE MG_Prolong_Y_Conformal_Mapping(phi,correc,nlev, mx, my1, mz, my)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev, mx, my, mz, my1
 
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my1,0:mz), INTENT (IN)  :: correc
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (OUT) :: phi

    INTEGER           :: ii, iim, jj, jjm, kk, kkm, i, j, k, nxyzMax, iErr
    REAL(KIND=CGREAL) :: var1, var2, var3, delta1, delta2, delta, d_1, d_2, BC2, BA1, BA2
    REAL(KIND=CGREAL), DIMENSION(0:ny+1) :: pro_wht1, pro_wht2, pro_wht3, pro_wht4
    INTEGER(1), DIMENSION(:,:,:), POINTER :: iblank_MG1

    IF (nlev > 1) THEN
        iblank_MG1 => MGY(nlev)%iblank
    ELSE
        iblank_MG1 => iblank
    ENDIF

!   Start Prolongation
!   ------------------
!   Prolongate in y-direction
!   -------------------------
    DO j = 1, my1-1 !mgrid_J(nlev+1) 
        jj  = 2*j 
        jjm = 2*j -1 

        d_1   = oned/(yc_MG(jjm,nlev) -yc_MG(j  ,nlev+1))**2 
        d_2   = oned/(yc_MG(jjm,nlev) -yc_MG(j-1,nlev+1))**2
        delta = d_1 + d_2  

        pro_wht1(j) = d_1/delta 
        pro_wht2(j) = d_2/delta

        d_1   = oned/(yc_MG(jj,nlev) -yc_MG(j  ,nlev+1))**2
        d_2   = oned/(yc_MG(jj,nlev) -yc_MG(j+1,nlev+1))**2
        delta = d_1 + d_2

        pro_wht3(j) = d_1/delta
        pro_wht4(j) = d_2/delta 
    END DO ! j
 
    DO k = 1, mz-1
    DO i = 1, mx-1
        phi(i,1, k) = correc(i, 1, k)

        delta1 =  pro_wht3(1) * (1-iblank_MG(i,1,k))
        delta2 =  pro_wht4(1) * (1-iblank_MG(i,2,k))
        delta  = delta1 + delta2

        phi(i,2,k) = delta1 *correc(i,1,k) &
                   + delta2 *correc(i,2,k)

        IF ( delta > zero ) THEN
          phi(i,2,k) = phi(i,2,k)/delta
        ENDIF ! delta

        DO j = 2, my1-2 !mgrid_J(nlev+1)-1
          jj = 2*j
          jjm = 2*j - 1

          var1 = correc(i,j  ,k)
          var2 = correc(i,j-1,k)
          var3 = correc(i,j+1,k)

!         prolongation
!         ------------
          delta1 = pro_wht1(j) * (1-iblank_MG(i,j  ,k))
          delta2 = pro_wht2(j) * (1-iblank_MG(i,j-1,k))
          delta  = delta1 +delta2

          phi(i,jjm,k) = delta1*var1 +delta2*var2

          IF ( delta > zero ) THEN
            phi(i,jjm,k) = phi(i,jjm,k)/delta
          ENDIF ! delta

          delta1 = pro_wht3(j) * (1-iblank_MG(i,j  ,k))
          delta2 = pro_wht4(j) * (1-iblank_MG(i,j+1,k))
          delta  = delta1 +delta2

          phi(i,jj,k) = delta1*var1 +delta2*var3

          IF ( delta > zero ) THEN
            phi(i,jj,k) =  phi(i,jj,k)/delta
          ENDIF ! delta

          IF (iblank_MG1(i,jjm,k)==1 .and. iblank_MG1(i,jj,k)==0 ) THEN
            BC2=(half*dy_mg(j+1,nlev+1)+dy_mg(jj,nlev))**2
            BA1=(half*dy_mg(j  ,nlev+1)-dy_mg(jj,nlev))**2
            BA2=(half*dy_mg(jj,nlev))**2
            phi(i,jj,k) = var1+(BA2-BA1)/(BC2-BA1)*(var3-var1)
          ENDIF
          IF (iblank_MG1(i,jjm,k)==0 .and. iblank_MG1(i,jj,k)==1 ) THEN
            BC2=(half*dy_mg(j-1,nlev+1)+dy_mg(jjm,nlev))**2
            BA1=(half*dy_mg(j  ,nlev+1)-dy_mg(jjm,nlev))**2
            BA2=(half*dy_mg(jjm,nlev))**2
            phi(i,jjm,k) = var1+(BA2-BA1)/(BC2-BA1)*(var2-var1)
          ENDIF

        END DO ! j

        j   = my1-1 !mgrid_J(nlev+1)
        jj  = 2*j
        jjm = 2*j -1

        delta1 = pro_wht1(j) * (1-iblank_MG(i,j  ,k))
        delta2 = pro_wht2(j) * (1-iblank_MG(i,j-1,k))
        delta  = delta1 +delta2

        phi(i,jjm,k) = delta1 *correc(i,j  ,k) &
                     + delta2 *correc(i,j-1,k)

        IF ( delta > zero ) THEN
          phi(i,jjm,k) = phi(i,jjm,k)/delta  
        ENDIF

        phi(i,jj,k) = correc(i,j,k) 

!       new modification - H. Dong
!       --------------------------
        IF ( nlev==1 .and. my-1 > 2*j ) THEN
          phi(i,my-1,k) = correc(i,j,k)
        ENDIF ! nLev
!       --------------------
!       end new modification

    END DO ! i  
    END DO ! k
    
END SUBROUTINE MG_Prolong_Y_Conformal_Mapping
!---------------------------------------------------------------------

SUBROUTINE MG_Prolong_Z_Conformal_Mapping(phi,correc,nlev, mx, my, mz1, mz)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlev,  mx, my, mz, mz1
 
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz1), INTENT (IN)  :: correc
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (OUT) :: phi

    INTEGER           :: ii, iim, jj, jjm, kk, kkm, i, j, k, nxyzMax, iErr
    REAL(KIND=CGREAL) :: var1, var2, var3, delta1, delta2, delta, d_1, d_2
    REAL(KIND=CGREAL), DIMENSION(0:nz+1) :: pro_wht1, pro_wht2, pro_wht3, pro_wht4
    INTEGER(1), DIMENSION(:,:,:), POINTER :: iblank_MG1
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: diblank_MG

    IF (nlev > 1) THEN
        iblank_MG1 => MGZ(nlev)%iblank
    ELSE
        iblank_MG1 => iblank
    ENDIF
     
!   Start Prolongation
!   ------------------
!   Prolongate in z-direction
!   -------------------------

    DO k = 1, mz1-1 !mgrid_K(nlev+1) 
        kk  = 2*k
        kkm = 2*k -1 

        d_1   = oned/(zc_MG(kkm,nlev) -zc_MG(k  ,nlev+1))**2 
        d_2   = oned/(zc_MG(kkm,nlev) -zc_MG(k-1,nlev+1))**2
        delta = d_1 + d_2  

        pro_wht1(k) = d_1/delta 
        pro_wht2(k) = d_2/delta

        d_1 = oned/(zc_MG(kk,nlev) -zc_MG(k  ,nlev+1))**2
        d_2 = oned/(zc_MG(kk,nlev) -zc_MG(k+1,nlev+1))**2
        delta = d_1 +d_2

        pro_wht3(k) = d_1/delta
        pro_wht4(k) = d_2/delta 

    END DO ! k

    DO j = 1, my-1
    DO i = 1, mx-1
        phi(i,j,1) = correc(i, j, 1)

        delta1 = pro_wht3(1) * (1-iblank_MG(i,j,1))
        delta2 = pro_wht4(1) * (1-iblank_MG(i,j,2))
        delta  = delta1 +delta2

        phi(i,j,2) = delta1 *correc(i,j,1) &
                   + delta2 *correc(i,j,2)

        IF ( delta > zero ) THEN
          phi(i,j,2) = phi(i,j,2)/delta
        ENDIF ! delta

        DO k = 2, mz1-2 !mgrid_K(nlev+1)-1
          kk  = 2*k
          kkm = 2*k -1

          var1 = correc(i,j,k)
          var2 = correc(i,j,k-1)
          var3 = correc(i,j,k+1)

!         prolongation
!         ------------
          delta1 = pro_wht1(k) * (1-iblank_MG(i,j,k  ))
          delta2 = pro_wht2(k) * (1-iblank_MG(i,j,k-1))
          delta  = delta1 +delta2
 
          phi(i,j,kkm) = delta1*var1 +delta2*var2

          IF ( delta > zero ) THEN
            phi(i,j,kkm) = phi(i,j,kkm)/delta
          ENDIF ! delta
             
          delta1 = pro_wht3(k) * (1-iblank_MG(i,j,k  ))
          delta2 = pro_wht4(k) * (1-iblank_MG(i,j,k+1))
          delta  = delta1 +delta2

          phi(i,j,kk) = delta1*var1 +delta2*var3
 
          IF ( delta > zero ) THEN
            phi(i,j,kk) =  phi(i,j,kk)/delta
          ENDIF ! delta
        END DO ! k
 
        k  = mz1-1 !mgrid_K(nlev+1)
        kk  = 2*k
        kkm = 2*k -1

        delta1 = pro_wht1(k) *(1-iblank_MG(i,j,k  ))
        delta2 = pro_wht2(k) *(1-iblank_MG(i,j,k-1))
        delta  = delta1 +delta2

        phi(i,j,kkm) = delta1*correc(i,j,k  ) &
                     + delta2*correc(i,j,k-1)

        IF ( delta > zero ) THEN
          phi(i,j,kkm) = phi(i,j,kkm)/delta  
        ENDIF ! delta

        phi(i,j,kk) = correc(i,j,k)  

!       new modification -H. Dong
!       -------------------------
        IF ( nlev == 1 .and. mz-1 > 2*k ) THEN
          phi(i,j,mz-1) = correc(i,j,k)
        ENDIF ! nLev
!       --------------------
!       end new modification
              
    END DO ! i  
    END DO ! j
    
END SUBROUTINE MG_Prolong_Z_Conformal_Mapping
!---------------------------------------------------------------------
!---------------------------------------------------------------------

SUBROUTINE MG_Prepare_Iblank_Conformal_Mapping(NC_dim)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

    
    INTEGER, INTENT(IN)   :: NC_dim

    INTEGER           :: i, j, k, n, ii, jj, kk
    REAL(KIND=CGREAL) :: ttemp, area1, area2, area


    IF ( internal_boundary_present == 1 ) THEN

      SELECT CASE (NC_dim)

      CASE(ICOORD)
        DO n = 2, mgLevels_X  
          DO k = 1, nzc
          DO j = 1, nyc
          DO i = 1, mgrid_I(n)
            area=zero
            DO ii=(i-1)*2**(n-1)+1, i*2**(n-1)
              area=area+dx(ii)*iblank(ii,j,k)
            ENDDO
            area1 = x_MG(i+1,n) - x_MG(i,n)
            area=area/area1
            IF (1-area<1e-6) THEN
            MGX(n)%iblank(i,j,k) = 1
            ELSE
            MGX(n)%iblank(i,j,k) = 0
            ENDIF
          END DO ! i
          END DO ! j
          END DO ! k
        END DO ! n

      CASE(JCOORD)
        DO n = 2, mgLevels_Y  
          DO k = 1, nzc
          DO i = 1, nxc
          DO j = 1, mgrid_J(n)
            area=zero
            DO jj=(j-1)*2**(n-1)+1, j*2**(n-1)
              area=area+dy(jj)*iblank(i,jj,k)
            ENDDO
            area1 = y_MG(j+1,n) - y_MG(j,n)
            area=area/area1
            IF ( 1-area<1e-6) THEN
            MGY(n)%iblank(i,j,k) = 1
            ELSE
            MGY(n)%iblank(i,j,k) = 0
            ENDIF
          END DO ! j
          END DO ! i
          END DO ! k
        END DO ! n

      CASE(KCOORD)
        DO n = 2, mgLevels_Z
          DO j = 1, nyc
          DO i = 1, nxc
          DO k = 1, mgrid_K(n)
            area=zero
            DO kk=(k-1)*2**(n-1)+1, k*2**(n-1)
              area=area+dz(kk)*iblank(i,j,kk)
            ENDDO
            area1 = z_MG(k+1,n) - z_MG(k,n)
            area=area/area1
            IF ( 1-area<1e-6) THEN
            MGZ(n)%iblank(i,j,k) = 1
            ELSE
            MGZ(n)%iblank(i,j,k) = 0
            ENDIF
          END DO ! k
          END DO ! j
          END DO ! i
        END DO ! n
      END SELECT ! NC_dim 

    ELSE  ! internal_boundary
      
      SELECT CASE (NC_dim)

      CASE(ICOORD)

        DO n = 2, mgLevels_X
           MGX(n)%iblank(1:mgrid_I(n),1:nyc,1:nzc) = 0
        ghostcellMark_MG(1:mgrid_I(n),1:nyc,1:nzc) = 0
        END DO ! n

      CASE(JCOORD)
      
        DO n = 2, mgLevels_Y
           MGY(n)%iblank(1:nxc,1:mgrid_J(n),1:nzc) = 0
        ghostcellMark_MG(1:nxc,1:mgrid_J(n),1:nzc) = 0
        END DO ! n

      CASE(KCOORD)
      
        DO n = 2, mgLevels_Z
           MGZ(n)%iblank(1:nxc,1:nyc,1:mgrid_K(n)) = 0
        ghostcellMark_MG(1:nxc,1:nyc,1:mgrid_K(n)) = 0
        END DO ! n

      END SELECT ! NC_dim

    ENDIF ! internal_boundary_present

END SUBROUTINE MG_Prepare_Iblank_Conformal_Mapping
!---------------------------------------------------------------------


!-------------------------------------------------------------------------------
  SUBROUTINE MG_Precondition_MSIP_Conformal_Mapping(nLevX, nLevY, nLevZ, MX, MY, MZ)
    USE global_parameters
    USE flow_parameters
!    USE flow_arrays
!    USE grid_arrays
!    USE boundary_arrays
!    USE multiuse_arrays
!    USE solver_arrays
    USE GCM_arrays
!    USE pressure_arrays
    USE mg_parameters
    USE mg_arrays    

    IMPLICIT NONE

    INTEGER :: nLevX, nLevY, nLevZ, MX, MY, MZ
    
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: A, B, C, D, E, F, G, H, P, R, S, U, V
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER :: I,J,K, II,JJ,KK
    REAL(KIND=CGREAL) :: PHI1, PHI2, PHI3, PHI4, BETA
    REAL(KIND=CGREAL) :: PHI5, PHI6, PHI7, PHI8
    REAL(KIND=CGREAL) :: PHI9, PHI10, PHI11, PHI12
    
    IF (nLevY==1 .AND. nLevZ==1) THEN
        B=>MGX(nLevX)%LU(:,:,:,1)
        C=>MGX(nLevX)%LU(:,:,:,2)
        D=>MGX(nLevX)%LU(:,:,:,3)
        E=>MGX(nLevX)%LU(:,:,:,4)
        F=>MGX(nLevX)%LU(:,:,:,5)
        G=>MGX(nLevX)%LU(:,:,:,6)
        H=>MGX(nLevX)%LU(:,:,:,7)
        
        AS=>MGX(nLevX)%CA(1,:,:,:)
        AW=>MGX(nLevX)%CA(2,:,:,:)
        AP=>MGX(nLevX)%CA(3,:,:,:)
        AE=>MGX(nLevX)%CA(4,:,:,:)
        AN=>MGX(nLevX)%CA(5,:,:,:)

        IF (ndim == DIM_3D) THEN
            A=>MGX(nLevX)%LU(:,:,:,8)
            P=>MGX(nLevX)%LU(:,:,:,9)
            R=>MGX(nLevX)%LU(:,:,:,10)
            S=>MGX(nLevX)%LU(:,:,:,11)
            U=>MGX(nLevX)%LU(:,:,:,12)
            V=>MGX(nLevX)%LU(:,:,:,13)

            AF=>MGX(nLevX)%CA(6,:,:,:)
            AB=>MGX(nLevX)%CA(7,:,:,:)
        ENDIF
    ELSE IF (NLEVX==1 .AND. nLevZ==1) THEN
        B=>MGY(nLevY)%LU(:,:,:,1)
        C=>MGY(nLevY)%LU(:,:,:,2)
        D=>MGY(nLevY)%LU(:,:,:,3)
        E=>MGY(nLevY)%LU(:,:,:,4)
        F=>MGY(nLevY)%LU(:,:,:,5)
        G=>MGY(nLevY)%LU(:,:,:,6)
        H=>MGY(nLevY)%LU(:,:,:,7)
        
        AS=>MGY(nLevY)%CA(1,:,:,:)
        AW=>MGY(nLevY)%CA(2,:,:,:)
        AP=>MGY(nLevY)%CA(3,:,:,:)
        AE=>MGY(nLevY)%CA(4,:,:,:)
        AN=>MGY(nLevY)%CA(5,:,:,:)

        IF (ndim == DIM_3D) THEN
            A=>MGY(nLevY)%LU(:,:,:,8)
            P=>MGY(nLevY)%LU(:,:,:,9)
            R=>MGY(nLevY)%LU(:,:,:,10)
            S=>MGY(nLevY)%LU(:,:,:,11)
            U=>MGY(nLevY)%LU(:,:,:,12)
            V=>MGY(nLevY)%LU(:,:,:,13)

            AF=>MGY(nLevY)%CA(6,:,:,:)
            AB=>MGY(nLevY)%CA(7,:,:,:)
        ENDIF
    ELSE
        B=>MGZ(nLevZ)%LU(:,:,:,1)
        C=>MGZ(nLevZ)%LU(:,:,:,2)
        D=>MGZ(nLevZ)%LU(:,:,:,3)
        E=>MGZ(nLevZ)%LU(:,:,:,4)
        F=>MGZ(nLevZ)%LU(:,:,:,5)
        G=>MGZ(nLevZ)%LU(:,:,:,6)
        H=>MGZ(nLevZ)%LU(:,:,:,7)
        
        AS=>MGZ(nLevZ)%CA(1,:,:,:)
        AW=>MGZ(nLevZ)%CA(2,:,:,:)
        AP=>MGZ(nLevZ)%CA(3,:,:,:)
        AE=>MGZ(nLevZ)%CA(4,:,:,:)
        AN=>MGZ(nLevZ)%CA(5,:,:,:)

        A=>MGZ(nLevZ)%LU(:,:,:,8)
        P=>MGZ(nLevZ)%LU(:,:,:,9)
        R=>MGZ(nLevZ)%LU(:,:,:,10)
        S=>MGZ(nLevZ)%LU(:,:,:,11)
        U=>MGZ(nLevZ)%LU(:,:,:,12)
        V=>MGZ(nLevZ)%LU(:,:,:,13)

        AF=>MGZ(nLevZ)%CA(6,:,:,:)
        AB=>MGZ(nLevZ)%CA(7,:,:,:)
    ENDIF

    IF (ndim /= DIM_3D) THEN
    K=1
!    IF (nLevX==1 .AND. nLevY==1 .AND. nLevZ==1) THEN
      DO J=1, MY
	  DO I=1, MX
	  IF (iblank_MG(I,J,K)==1) CYCLE

		  ae(i,j,K)=dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)
		  aw(i,j,K)=dxcinv_mg(i,  nLevX)*dxinv_mg(i,nLevX)
		  an(i,j,K)=dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)
		  as(i,j,K)=dycinv_mg(j,  nLevY)*dyinv_mg(j,nLevY)
		  ap(i,j,K)=-( ae(i,j,K) + aw(i,j,K) + an(i,j,K) +  as(i,j,K))
          
		  ap(i,j,K)=ap(i,j,K) + ae(i,j,K)*iup_MG(i,j,k) + aw(i,j,K)*ium_MG(i,j,k) &
							  + an(i,j,K)*jup_MG(i,j,k) + as(i,j,K)*jum_MG(i,j,k)
		  ae(i,j,K)=ae(i,j,K)*(oned - iup_MG(i,j,k) )
		  aw(i,j,K)=aw(i,j,K)*(oned - ium_MG(i,j,k) )
		  an(i,j,K)=an(i,j,K)*(oned - jup_MG(i,j,k) )
		  as(i,j,K)=as(i,j,K)*(oned - jum_MG(i,j,k) )

	  END DO
	  END DO
!    ELSE
!	  DO J=1, MY
!	  DO I=1, MX
!	  IF (iblank_MG(I,J,K)==1) CYCLE
!        IF (iblank_MG(i+1,j,k)/=1) THEN
!		  ae(i,j,K)=2*dxinv_mg(i,nLevX)/(dx_mg(i,nLevX)+rblank_MG(i+1,j,k)*dx_mg(i+1,nLevX))
!		ELSE
!		  IF (abs(rblank_MG(i,j,k)-oned)<1e-6) THEN ! Pure fluid
!		  ae(i,j,K)=dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)
!		  ELSE
!		  ae(i,j,K)=dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)
!		  ENDIF
!		ENDIF
!        IF (iblank_MG(i-1,j,k)/=1) THEN
!		  aw(i,j,K)=2*dxinv_mg(i,nLevX)/(dx_mg(i,nLevX)+rblank_MG(i-1,j,k)*dx_mg(i-1,nLevX))
!		ELSE
!		  aw(i,j,K)=dxcinv_mg(i,  nLevX)*dxinv_mg(i,nLevX)
!        ENDIF  
!        IF (iblank_MG(i,j+1,k)/=1) THEN
!		  an(i,j,K)=2*dyinv_mg(j,nLevY)/(dy_mg(j,nLevY)+rblank_MG(i,j+1,k)*dy_mg(j+1,nLevY))
!		ELSE
!		  an(i,j,K)=dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)
!        ENDIF  
!        IF (iblank_MG(i,j-1,k)/=1) THEN
!		  as(i,j,K)=2*dyinv_mg(j,nLevY)/(dy_mg(j,nLevY)+rblank_MG(i,j-1,k)*dy_mg(j-1,nLevY))
!		ELSE
!		  as(i,j,K)=dycinv_mg(j,  nLevY)*dyinv_mg(j,nLevY)
!        ENDIF  
!		ap(i,j,K)=-( ae(i,j,K) + aw(i,j,K) + an(i,j,K) +  as(i,j,K))
!            
!		ap(i,j,K)=ap(i,j,K) + ae(i,j,K)*iup_MG(i,j,k) + aw(i,j,K)*ium_MG(i,j,k) &
!  						    + an(i,j,K)*jup_MG(i,j,k) + as(i,j,K)*jum_MG(i,j,k)
!		ae(i,j,K)=ae(i,j,K)*(oned - iup_MG(i,j,k) )
!		aw(i,j,K)=aw(i,j,K)*(oned - ium_MG(i,j,k) )
!		an(i,j,K)=an(i,j,K)*(oned - jup_MG(i,j,k) )
!		as(i,j,K)=as(i,j,K)*(oned - jum_MG(i,j,k) )
!	  END DO
!	  END DO
!    ENDIF
    
    BETA=OMEGA
    F=zero
    G=zero
    H=zero
    
    KK=1
    DO J=1, MY
    JJ=J+1
    DO I=1, MX
    II=I+1
    IF (iblank_MG(I,J,K)==1) CYCLE

        B(II,JJ,KK)=AS(I,J,K)/(1.D0-BETA*F(II,JJ-1,KK)*F(II+1,JJ-1,KK))      ! b 
        C(II,JJ,KK)=-B(II,JJ,KK)*F(II,JJ-1,KK)                             ! c
        D(II,JJ,KK)=(AW(I,J,K)-B(II,JJ,KK)*G(II,JJ-1,KK))/(1+2.D0*BETA*G(II-1,JJ,KK))

        PHI1=C(II,JJ,KK)*F(II+1,JJ-1,KK)
        PHI4=D(II,JJ,KK)*G(II-1,JJ,KK)

        E(II,JJ,KK)=1.D0/(AP(I,J,K)-B(II,JJ,KK)*H(II,JJ-1,KK)-C(II,JJ,KK)*G(II+1,JJ-1,KK)- &
                                    D(II,JJ,KK)*F(II-1,JJ,KK)+2.D0*BETA*(PHI1+PHI4))            ! e
        F(II,JJ,KK)=(AE(I,J,K)-C(II,JJ,KK)*H(II+1,JJ-1,KK)-2.D0*BETA*PHI1)*E(II,JJ,KK)            ! f
        G(II,JJ,KK)=-D(II,JJ,KK)*H(II-1,JJ,KK)*E(II,JJ,KK)  
        H(II,JJ,KK)=(AN(I,J,K)-BETA*PHI4)*E(II,JJ,KK)  

    END DO ! I
    END DO ! J
    
    ELSE

    ENDIF
  END SUBROUTINE MG_Precondition_MSIP_Conformal_Mapping
!-------------------------------------------------------------------------------

  SUBROUTINE GCM_Pressure(var, rhs)

    USE global_parameters
    USE flow_parameters
!    USE grid_arrays
!    USE boundary_arrays
!    USE multiuse_arrays
!    USE GCM_arrays
!    USE MG_parameters
!    USE MG_arrays
 
    IMPLICIT NONE

    REAL(KIND=CGREAL),DIMENSION(0:nx,0:ny,0:nz),INTENT (INOUT) :: var, rhs
    
    CALL set_outer_ghost_pres(var,nx,ny,nz)
    CALL GCM_p_set_bc_internal(var,nx,ny,nz)
    CALL GCM_enforce_p_compatibility(var,nx,ny,nz)

    CALL MG_Prepare_BC(1)

!    SELECT CASE ( pp_solver_type )
!    CASE (PP_SOLVER_TYPE_MG, PP_SOLVER_TYPE_MG_Point_Jacobi)
!
!      CALL MG_Prepare_BC(1)
!
!    CASE (PP_SOLVER_TYPE_MG_MSIP)
!
!      IF (ndim == DIM_2D) THEN
!        CALL GCM_correct_pressure_MSIP_2D(rhs,nx,ny,nz)
!      ELSE
!        CALL GCM_MG_Precondition_MSIP_3D
!      ENDIF
!
!    END SELECT

  END SUBROUTINE GCM_Pressure


