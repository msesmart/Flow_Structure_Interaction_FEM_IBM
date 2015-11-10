!******************************************************************************
!
! Purpose: generalized kernel to evolve phi_{LM} and phi_{MM} equations
!
! Description: none.
!
! Input: field variables
!
! Output: phi_{LM} and phi_{MM} 
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
   SUBROUTINE TURB_InitPhiField( )

!==============================================================================
!  Purpose: Set initial field for phiLM and phiMM based on a smagorinsky model
!==============================================================================

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE turb_parameters
    USE turb_arrays
    USE grid_arrays

    IMPLICIT NONE

!... Parameters

!... Loop variables

    INTEGER :: i, j, k
 
!... Local variables

    REAL(KIND=CGREAL) :: cFixLagrSqr, mijLij, mijMij 

!******************************************************************************

!------------------------------------------------------------------------------
! Set dimensions
!------------------------------------------------------------------------------   

    cFixLagrSqr = cSmagFix**2

!------------------------------------------------------------------------------
! Calculate Lij and Mij terms 
!------------------------------------------------------------------------------

    WRITE(STDOUT,*) '   Entering TURB_CalcMij ' 
    CALL TURB_CalcMij()

    WRITE(STDOUT,*) '   Entering TURB_CalcLij ' 
    CALL TURB_CalcLij()

!------------------------------------------------------------------------------
! Contract M_{ij}*M_{ij} and M_{ij}*L_{ij}
!  No homogenization is imposed
!------------------------------------------------------------------------------

    WRITE(STDOUT,*) '   Entering TURB_ContractLijMij '
    CALL TURB_ContractLijMij()

!------------------------------------------------------------------------------
! Compute Phi_{LM} and Phi_{MM} at initial step 
!------------------------------------------------------------------------------

    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1

!------------------------------------------------------------------------------
!     Load values
!------------------------------------------------------------------------------

      mijLij = lij(S11,i,j,k)
      mijMij = mij(S11,i,j,k)

!------------------------------------------------------------------------------      
!     Construct initial values of phiLM and phiMM
!------------------------------------------------------------------------------

      phiLM(i,j,k) = cFixLagrSqr * mijLij  
      phiMM(i,j,k) = cFixLagrSqr * mijMij

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
!   Monitor output
!------------------------------------------------------------------------------

    PRINT*,' MIN/MAX Vals of phiLM-Init = ',MINVAL(phiLM(1:nx-1,1:ny-1,1:nz-1)),&
                                            MAXVAL(phiLM(1:nx-1,1:ny-1,1:nz-1))

    PRINT*,' MIN/MAX Vals of phiMM-Init = ',MINVAL(phiMM(1:nx-1,1:ny-1,1:nz-1)),&
                                            MAXVAL(phiMM(1:nx-1,1:ny-1,1:nz-1))
        
   END SUBROUTINE TURB_InitPhiField
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_CalcCsDynLagr( )

!==============================================================================
!  Purpose: Evolve phiLM and phiMM based on ODE
!==============================================================================

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE flow_arrays
    USE turb_parameters
    USE turb_arrays
    USE grid_arrays
    USE boundary_arrays
    
!    USE TURB_ModInterfaces, ONLY : TURB_CalcFilterWidth

    IMPLICIT NONE

!... Parameters

!... Loop variables

    INTEGER :: i, j, k
 
!... Local variables

    INTEGER :: iFlagFilter,iPrev, jPrev, kPrev, ii, jj, kk
    INTEGER, DIMENSION(3) :: indexCurr, indexPrev 
    
    REAL(KIND=CGREAL), PARAMETER :: oneEight = -oned/8.0_CGREAL
    REAL(KIND=CGREAL) :: deltaSqr, eps, mijMij, mijLij, tauLagr
    REAL(KIND=CGREAL) :: phiLMPrev, phiLMCurr, phiLMTerm, phiMMCurr, phiMMPrev
    REAL(KIND=CGREAL) :: diffPosX, diffPosY, diffPosZ
    REAL(KIND=CGREAL) :: posPrevX, posPrevY, posPrevZ, strainMagnBar
    REAL(KIND=CGREAL) :: wx,wy,wz
    REAL(KIND=CGREAL), DIMENSION(3) :: posPrev
    REAL(KIND=CGREAL), DIMENSION(3) :: posCurr

    REAL(KIND=CGREAL) :: epsLagr

!******************************************************************************

!------------------------------------------------------------------------------
! Set dimensions
!------------------------------------------------------------------------------   

    eps      = EPSILON(oneEight)
    iFlagFilter = ACTIVE-1

!------------------------------------------------------------------------------
! Initialize arrays
!------------------------------------------------------------------------------
 
!    DO k = 0,nz+1
!    DO j = 0,ny+1
!    DO i = 0,nx+1
!      viscTot(i,j,k) = zero
!    ENDDO ! i
!    ENDDO ! j
!    ENDDO ! k

!------------------------------------------------------------------------------
! Compute Phi_{LM} and Phi_{MM} for the particle coordinates 
!  at the previous time step 
!------------------------------------------------------------------------------

    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1

!------------------------------------------------------------------------------
!      Load current cell indices
!------------------------------------------------------------------------------

      indexCurr(DIRX) = i
      indexCurr(DIRY) = j
      indexCurr(DIRZ) = k

!------------------------------------------------------------------------------
!     Compute particle location at previous timestep
!------------------------------------------------------------------------------

      posPrevX = xc(i) -u(i,j,k)*dt
      posPrevY = yc(j) -v(i,j,k)*dt
      posPrevZ = zc(k) -w(i,j,k)*dt

!------------------------------------------------------------------------------
!     Determine the cell where the particle was at t-dt
!     Simple search kernel based on the sign of the velocity field.
!     Since CFL <= 1, the particle can not move more than one grid
!     cell in either direction.
!------------------------------------------------------------------------------

      iPrev = i
      jPrev = j
      kPrev = k
      
      IF ( u(i,j,k) > zero ) iPrev = MAX( i-1, 1 ) 
      IF ( v(i,j,k) > zero ) jPrev = MAX( j-1, 1 )
      IF ( w(i,j,k) > zero ) kPrev = MAX( k-1, 1 )

      ii=iPrev; jj=jPrev; kk=kPrev;

      indexPrev(DIRX) = iPrev
      indexPrev(DIRY) = jPrev
      indexPrev(DIRZ) = kPrev

      posPrev(DIRX) = posPrevX
      posPrev(DIRY) = posPrevY
      posPrev(DIRZ) = posPrevZ

!------------------------------------------------------------------------------
!     Interpolate solution onto previous particle cell location 
!      using tri-linear scheme
!------------------------------------------------------------------------------     

      diffPosX = posPrevX - xc(iPrev)
      diffPosY = posPrevY - yc(jPrev)
      diffPosZ = posPrevZ - zc(kPrev)

      wx = diffPosX / ( xc(iPrev+1) - xc(iPrev) )
      wy = diffPosY / ( yc(jPrev+1) - yc(jPrev) )
      wz = diffPosZ / ( zc(kPrev+1) - zc(kPrev) )

      phiLMPrev = wz *( (oned-wy)*( (oned-wx)*phiLM(ii  ,jj  ,kk+1)      &
                +                                     wx *phiLM(ii+1,jj  ,kk+1) )    &
                +                   wy *( (oned-wx)*phiLM(ii  ,jj+1,kk+1)      &
                +                                     wx *phiLM(ii+1,jj+1,kk+1) ) )  &
    + (oned-wz)*( (oned-wy)*( (oned-wx)*phiLM(ii  ,jj  ,kk  )      &
                +                                     wx *phiLM(ii+1,jj  ,kk  )   )  &
                +                   wy *( (oned-wx)*phiLM(ii  ,jj+1,kk  )      &
                +                                     wx *phiLM(ii+1,jj+1,kk  )   )  )

      phiMMPrev = wz *( (oned-wy)*( (oned-wx)*phiLM(ii  ,jj  ,kk+1)      &
                +                                     wx *phiLM(ii+1,jj  ,kk+1) )    &
                +                   wy *( (oned-wx)*phiLM(ii  ,jj+1,kk+1)      &
                +                                     wx *phiLM(ii+1,jj+1,kk+1) ) )  &
    + (oned-wz)*( (oned-wy)*( (oned-wx)*phiLM(ii  ,jj  ,kk  )      &
                +                                     wx *phiLM(ii+1,jj  ,kk  )   )  &
                +                   wy *( (oned-wx)*phiLM(ii  ,jj+1,kk  )      &
                +                                     wx *phiLM(ii+1,jj+1,kk  )   )  )

!------------------------------------------------------------------------------
!     Evaluate strain magnitude 
!     |S(u)| = SQRT( 2 S_{ij}*S_{ij} )
!------------------------------------------------------------------------------

      strainMagnBar = twoSqrt * SQRT( sij(S11,i,j,k)*sij(S11,i,j,k) &
                                    + sij(S22,i,j,k)*sij(S22,i,j,k) &
                                    + sij(S33,i,j,k)*sij(S33,i,j,k) &    
                        + 2.0_CGREAL* sij(S12,i,j,k)*sij(S12,i,j,k) &
                        + 2.0_CGREAL* sij(S13,i,j,k)*sij(S13,i,j,k) &
                        + 2.0_CGREAL* sij(S23,i,j,k)*sij(S23,i,j,k) ) 


!------------------------------------------------------------------------------
!     Compute time scale and epsLagr based on activated model
!------------------------------------------------------------------------------

      SELECT CASE(turbLagrTimeScale)
        CASE(TURB_LAGR_TSCALE_TSL)
          tauLagr = strainMagnBar
          epsLagr = tauLagr*dt/(oned+ tauLagr*dt)

        CASE(TURB_LAGR_TSCALE_JFM)
!        Setting formalism from JFM Paper (Meneuveau et al)
!        Do not use ratio to avoid division by zero

         CALL TURB_CalcFilterWidth(iFlagFilter,indexCurr,deltaSqr)

         IF ( iblank(i,j,k) /= 1 .OR. ABS( phiLMPrev*phiMMPrev ) > eps ) THEN
           tauLagr = 1.5_CGREAL*SQRT(deltaSqr) &
                   * ABS( phiLMPrev*phiMMPrev )**oneEight
         ENDIF ! iblank 

         epsLagr = dt/(tauLagr+dt)

      END SELECT ! turbLagrTimeScale

!------------------------------------------------------------------------------
!     Load values
!------------------------------------------------------------------------------

      mijLij = lij(S11,i,j,k)
      mijMij = mij(S11,i,j,k)

!------------------------------------------------------------------------------      
!    Construct phiLM and phiMMCurr w/o a ramp function 
!    to allow backscatter (see Donghyun comment)
!------------------------------------------------------------------------------

      phiLMCurr = epsLagr *mijLij + (oned-epsLagr) *phiLMPrev   

      phiMMCurr = epsLagr *mijMij + (oned-epsLagr) *phiMMPrev  

!------------------------------------------------------------------------------ 
!     Apply internal boundary conditions
!------------------------------------------------------------------------------ 

      phiLMCurr = phiLMCurr * REAL(1-iblank(i,j,k),KIND=CGREAL)
      phiMMCurr = phiMMCurr * REAL(1-iblank(i,j,k),KIND=CGREAL)

!------------------------------------------------------------------------------ 
!     Compute dynamic Lagrangian constant
!       add small value for degeneracy case of ratio of two small numbers
!------------------------------------------------------------------------------ 
      
      IF ( iblank(i,j,k) /= 1 ) THEN
        viscTot(i,j,k) = phiLMCurr/( oned +phiMMCurr )
      ENDIF ! iblank 

!------------------------------------------------------------------------------ 
!     Load values in intermediate arrays as not to clobber phi values 
!------------------------------------------------------------------------------ 

      mij(S12,i,j,k) = phiMMCurr
      lij(S12,i,j,k) = phiLMCurr

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------ 
!   Update phi values
!------------------------------------------------------------------------------ 

    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1
        phiLM(i,j,k) = mij(S12,i,j,k)
        phiMM(i,j,k) = lij(S12,i,j,k)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------ 
!   Apply boundary conditions
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 
!   x-direction
!------------------------------------------------------------------------------ 

    DO k=0,nz
    DO j=0,ny

! left boundary ---------------------------------------------------------------
      i = 0
      SELECT CASE (bcx1)
        CASE (BC_TYPE_DIRICHLET,BC_TYPE_PULSATILE_INFLOW, &
              BC_TYPE_USER_SPECIFIED,BC_TYPE_SHEAR)
          phiLM(i,j,k) = zero 
          phiMM(i,j,k) = oned 
        CASE (BC_TYPE_ZERO_GRADIENT,BC_TYPE_SYMMETRY)
          phiLM(i,j,k) = phiLM(i+1,j,k) 
          phiMM(i,j,k) = phiLM(i+1,j,k) 
        CASE (BC_TYPE_PERIODIC)             ! periodic bc
          phiLM(i,j,k) = phiLM(nx-1,j,k) 
          phiMM(i,j,k) = phiLM(nx-1,j,k) 
      END SELECT 

! right boundary---------------------------------------------------------------
      i = nx
      SELECT CASE (bcx2)
        CASE (BC_TYPE_DIRICHLET,BC_TYPE_PULSATILE_INFLOW, &
              BC_TYPE_USER_SPECIFIED,BC_TYPE_SHEAR)
          phiLM(i,j,k) = zero
          phiMM(i,j,k) = oned
        CASE (BC_TYPE_ZERO_GRADIENT,BC_TYPE_SYMMETRY)  
          phiLM(i,j,k) = phiLM(i-1,j,k)         
          phiMM(i,j,k) = phiLM(i-1,j,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc
          phiLM(i,j,k) = phiLM(2,j,k)
          phiMM(i,j,k) = phiLM(2,j,k) 
      END SELECT 
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------ 
!   y-direction
!------------------------------------------------------------------------------ 

    DO k=0,nz
    DO i=0,nx

! bottom boundary -------------------------------------------------------------
      j = 1
      SELECT CASE (bcy1)
        CASE (BC_TYPE_DIRICHLET,BC_TYPE_PULSATILE_INFLOW, &
              BC_TYPE_USER_SPECIFIED,BC_TYPE_SHEAR)
          phiLM(i,j,k) = zero
          phiMM(i,j,k) = oned
        CASE (BC_TYPE_ZERO_GRADIENT,BC_TYPE_SYMMETRY)  
          phiLM(i,j,k) = phiLM(i,j+1,k)         
          phiMM(i,j,k) = phiLM(i,j+1,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc
          phiLM(i,j,k) = phiLM(i,ny-1,k)
          phiMM(i,j,k) = phiLM(i,ny-1,k) 
      END SELECT 

! top boundary ----------------------------------------------------------------
      j = ny
      SELECT CASE (bcy2)
        CASE (BC_TYPE_DIRICHLET,BC_TYPE_PULSATILE_INFLOW, &
              BC_TYPE_USER_SPECIFIED,BC_TYPE_SHEAR)
          phiLM(i,j,k) = zero
          phiMM(i,j,k) = oned
        CASE (BC_TYPE_ZERO_GRADIENT,BC_TYPE_SYMMETRY)
          phiLM(i,j,k) = phiLM(i,j-1,k)
          phiMM(i,j,k) = phiLM(i,j-1,k)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc
          phiLM(i,j,k) = phiLM(i,2,k)
          phiMM(i,j,k) = phiLM(i,2,k)
      END SELECT 

    ENDDO ! i
    ENDDO ! k

    DO j=0,ny
    DO i=0,nx

!------------------------------------------------------------------------------ 
!   k-direction
!------------------------------------------------------------------------------ 

! front boundary --------------------------------------------------------------
      k = 1
      SELECT CASE (bcz1)
        CASE (BC_TYPE_DIRICHLET,BC_TYPE_PULSATILE_INFLOW, &
              BC_TYPE_USER_SPECIFIED,BC_TYPE_SHEAR)
          phiLM(i,j,k) = zero
          phiMM(i,j,k) = oned
        CASE (BC_TYPE_ZERO_GRADIENT,BC_TYPE_SYMMETRY)
          phiLM(i,j,k) = phiLM(i,j,k+1)
          phiMM(i,j,k) = phiLM(i,j,k+1)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc
          phiLM(i,j,k) = phiLM(i,j,nz-1)
          phiMM(i,j,k) = phiLM(i,j,nz-1)
      END SELECT 

! back boundary ---------------------------------------------------------------
      k = nz-1
      SELECT CASE (bcz2)
        CASE (BC_TYPE_DIRICHLET,BC_TYPE_PULSATILE_INFLOW, &
              BC_TYPE_USER_SPECIFIED,BC_TYPE_SHEAR)
          phiLM(i,j,k) = zero
          phiMM(i,j,k) = oned
        CASE (BC_TYPE_ZERO_GRADIENT,BC_TYPE_SYMMETRY)
          phiLM(i,j,k) = phiLM(i,j,k-1)
          phiMM(i,j,k) = phiLM(i,j,k-1)
        CASE (BC_TYPE_PERIODIC)             ! periodic bc
          phiLM(i,j,k) = phiLM(i,j,2)
          phiMM(i,j,k) = phiLM(i,j,2)
      END SELECT 

    ENDDO ! i
    ENDDO ! j

!------------------------------------------------------------------------------ 
! Monitor output
!------------------------------------------------------------------------------ 

    IF ( MOD(ntime,nmonitor) == 0 ) THEN 
      WRITE(STDOUT, *) ' MIN/MAX Vals of LijMij = ',&
                                         MINVAL(lij(S11,1:nx-1,1:ny-1,1:nz-1)),&
                                         MAXVAL(lij(S11,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *) ' MIN/MAX Vals of MijMij = ',&
                                         MINVAL(mij(S11,1:nx-1,1:ny-1,1:nz-1)),&
                                         MAXVAL(mij(S11,1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *) ' MIN/MAX Vals of phiLM = ',&
                                         MINVAL(phiLM(1:nx-1,1:ny-1,1:nz-1)),&
                                         MAXVAL(phiLM(1:nx-1,1:ny-1,1:nz-1))

      WRITE(STDOUT, *) ' MIN/MAX Vals of phiMM = ',&
                                         MINVAL(phiMM(1:nx-1,1:ny-1,1:nz-1)),&
                                         MAXVAL(phiMM(1:nx-1,1:ny-1,1:nz-1))
 
      WRITE(STDOUT, *) ' MIN/MAX Vals of cS-DynLagr = ',&
                                         MINVAL(viscTot(1:nx-1,1:ny-1,1:nz-1)),&
                                         MAXVAL(viscTot(1:nx-1,1:ny-1,1:nz-1))

111 FORMAT(2(2X,1PE15.7))

    ENDIF ! ntime

   END SUBROUTINE TURB_CalcCsDynLagr
!------------------------------------------------------------------------------
