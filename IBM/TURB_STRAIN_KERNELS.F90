!******************************************************************************
!
! Purpose: generalized kernel to compute the strain rate 
!          and filter widths
!
! Description: none.
!
! Input: field variables
!
! Output: strain rates, grid and test filter widths
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
   SUBROUTINE TURB_CalcStrainRate( q, qFace, strn )

!==============================================================================
!  Purpose: Compute strain-rate tensor based on collocated and face velocities
!  Note   : Face-velocities are used to insure divergence-free condition
!           on S_ii tensor
!==============================================================================

    USE global_parameters
    USE turb_global_parameters    
    USE flow_parameters
    USE turb_parameters
    USE grid_arrays
    USE boundary_arrays
    USE flow_arrays

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(3,0:nx+1,0:ny+1,0:nz+1), INTENT(IN)  :: q, qFace
    REAL(KIND=CGREAL), DIMENSION(6,0:nx+1,0:ny+1,0:nz+1), INTENT(OUT) :: strn

!... Loop variables

    INTEGER :: i, j, k
 
!... Local variables
    
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: ue, uw, un, us, uf, ub
    REAL(KIND=CGREAL) :: ve, vw, vn, vs, vf, vb
    REAL(KIND=CGREAL) :: we, ww, wn, ws, wf, wb
    REAL(KIND=CGREAL) :: qGrad(3,3)

!------------------------------------------------------------------------------
! Compute gradients at cell centers based on cell faces to maintain 
!  divergence-free condition
!------------------------------------------------------------------------------
    
    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1

!------------------------------------------------------------------------------
!     u-velocity
!------------------------------------------------------------------------------

      ue = qFace(DIRX,i+1,j,k)

      uw = qFace(DIRX,i  ,j,k)
      
      un = ( fy(j+1)*q(DIRX,i,j+1,k) + (oned-fy(j+1))*q(DIRX,i,j,k)   )*(1-jup(i,j,k)) &
         + q(DIRX,i,j,k)*jup(i,j,k)

      us = ( fy(j)  *q(DIRX,i,j,k)   + (oned-fy(j))  *q(DIRX,i,j-1,k) )*(1-jum(i,j,k)) &
         + q(DIRX,i,j,k)*jum(i,j,k)
      
      uf = ( fz(k+1)*q(DIRX,i,j,k+1) + (oned-fz(k+1))*q(DIRX,i,j,k)   )*(1-kup(i,j,k)) &
         + q(DIRX,i,j,k)*kup(i,j,k)

      ub = ( fz(k)  *q(DIRX,i,j,k)   + (oned-fz(k))  *q(DIRX,i,j,k-1) )*(1-kum(i,j,k)) &
         + q(DIRX,i,j,k)*kum(i,j,k)      

      qGrad(DIRX,GRADX)= (ue-uw)*dxinv(i)
      qGrad(DIRX,GRADY)= (un-us)*dyinv(j)
      qGrad(DIRX,GRADZ)= (uf-ub)*dzinv(k)

!------------------------------------------------------------------------------
!     v-velocity
!------------------------------------------------------------------------------

      ve = ( fx(i+1)*q(DIRY,i+1,j,k) + (oned-fx(i+1))*q(DIRY,i,j,k)   )*(1-iup(i,j,k)) &
         + q(DIRY,i,j,k)*iup(i,j,k)

      vw = ( fx(i)  *q(DIRY,i,j,k)   + (oned-fx(i))  *q(DIRY,i-1,j,k) )*(1-ium(i,j,k)) &
         + q(DIRY,i,j,k)*ium(i,j,k)

      vn = qFace(DIRY,i,j+1,k)

      vs = qFace(DIRY,i,j  ,k)

      vf = ( fz(k+1)*q(DIRY,i,j,k+1) + (oned-fz(k+1))*q(DIRY,i,j,k)   )*(1-kup(i,j,k)) &
         + q(DIRY,i,j,k)*kup(i,j,k)

      vb = ( fz(k)  *q(DIRY,i,j,k)   + (oned-fz(k))  *q(DIRY,i,j,k-1) )*(1-kum(i,j,k)) &
         + q(DIRY,i,j,k)*kum(i,j,k)

      qGrad(DIRY,GRADX)= (ve-vw)*dxinv(i)
      qGrad(DIRY,GRADY)= (vn-vs)*dyinv(j)
      qGrad(DIRY,GRADZ)= (vf-vb)*dzinv(k)

!------------------------------------------------------------------------------
!     w-velocity
!------------------------------------------------------------------------------

      we = ( fx(i+1)*q(DIRZ,i+1,j,k) + (oned-fx(i+1))*q(DIRZ,i,j,k)   )*(1-iup(i,j,k)) &
         + q(DIRZ,i,j,k)*iup(i,j,k)

      ww = ( fx(i)  *q(DIRZ,i,j,k)   + (oned-fx(i))  *q(DIRZ,i-1,j,k) )*(1-ium(i,j,k)) &
         + q(DIRZ,i,j,k)*ium(i,j,k)

      wn = ( fy(j+1)*q(DIRZ,i,j+1,k) + (oned-fy(j+1))*q(DIRZ,i,j,k)   )*(1-jup(i,j,k)) &
         + q(DIRZ,i,j,k)*jup(i,j,k)

      ws = ( fy(j)  *q(DIRZ,i,j,k)   + (oned-fy(j))  *q(DIRZ,i,j-1,k) )*(1-jum(i,j,k)) &
         + q(DIRZ,i,j,k)*jum(i,j,k)

      wf = qFace(DIRZ,i,j,k+1)

      wb = qFace(DIRZ,i,j,k  )

      qGrad(DIRZ,GRADX)= (we-ww)*dxinv(i)
      qGrad(DIRZ,GRADY)= (wn-ws)*dyinv(j)
      qGrad(DIRZ,GRADZ)= (wf-wb)*dzinv(k)

!------------------------------------------------------------------------------
!     Assemble individual components
!------------------------------------------------------------------------------

      strn(S11,i,j,k) = qGrad(DIRX,GRADX)                                   ! s11 = ux
      strn(S12,i,j,k) = half*(qGrad(DIRX,GRADY) + qGrad(DIRY,GRADX) ) ! s12 = 1/2(uy+vx)
      strn(S13,i,j,k) = half*(qGrad(DIRX,GRADZ) + qGrad(DIRZ,GRADX) ) ! s13 = 1/2(uz+wx)
      
      strn(S22,i,j,k) = qGrad(DIRY,GRADY)                                   ! s22 = vy
      strn(S23,i,j,k) = half*(qGrad(DIRY,GRADZ) + qGrad(DIRZ,GRADY) ) ! s13 = 1/2(vz+wy)
      
      strn(S33,i,j,k) = qGrad(DIRZ,GRADZ)                                   ! s33 = wz

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
! Apply immersed boundary conditions
!------------------------------------------------------------------------------

    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1
      riblank         = oned-REAL(iblank(i,j,k),KIND=CGREAL)
      strn(S11,i,j,k) = strn(S11,i,j,k) *riblank
      strn(S12,i,j,k) = strn(S12,i,j,k) *riblank
      strn(S13,i,j,k) = strn(S13,i,j,k) *riblank
      strn(S22,i,j,k) = strn(S22,i,j,k) *riblank
      strn(S23,i,j,k) = strn(S23,i,j,k) *riblank
      strn(S33,i,j,k) = strn(S33,i,j,k) *riblank
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
! Compute strain rate tensor for ghost cells based 
!  on weighted averaged formalism
!------------------------------------------------------------------------------

    IF ( boundary_formulation == GCM_METHOD ) CALL TURB_GCM_CalcStrainRate

   END SUBROUTINE TURB_CalcStrainRate
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE TURB_CalcFilterWidth(iFlag, index,deltaTestSqr)

!==============================================================================
!  Purpose: compute width for grid and test filters
!==============================================================================

    USE global_parameters
    USE turb_global_parameters
    USE turb_parameters
    USE flow_parameters 
    USE grid_arrays
    USE boundary_arrays
    USE turb_arrays
    
    IMPLICIT NONE

!... Parameters

    INTEGER,               INTENT(IN) :: iFlag
    INTEGER, DIMENSION(3), INTENT(IN) :: index
    
    REAL(KIND=CGREAL), INTENT(OUT) :: deltaTestSqr
    
!... Loop variables

!... Local variables

    INTEGER :: i,j,k
 
    REAL(KIND=CGREAL) :: alphaTwoThird, alphaX, alphaY, alphaZ

!------------------------------------------------------------------------------
! Set dimensions  
!------------------------------------------------------------------------------

    alphaX   = oned
    alphaY   = oned
    alphaZ   = oned

    i = index(DIRX)
    j = index(DIRY)
    k = index(DIRZ)

!------------------------------------------------------------------------------    
! Set appropriate alpha
!------------------------------------------------------------------------------

    IF ( iFlag == ACTIVE .AND. testFilterDir(DIRY) == ACTIVE ) alphaX = fWidthRatio(DIRX)    
    IF ( iFlag == ACTIVE .AND. testFilterDir(DIRY) == ACTIVE ) alphaY = fWidthRatio(DIRY)     
    IF ( iFlag == ACTIVE .AND. testFilterDir(DIRZ) == ACTIVE ) alphaZ = fWidthRatio(DIRZ) 

    alphaTwoThird = (alphaX* alphaY* alphaZ)**twoThird

!------------------------------------------------------------------------------  
! Select appropriate model
!------------------------------------------------------------------------------

    SELECT CASE( filterWidthModel )
      CASE( TURB_FWIDTH_CUBEROOT )
        deltaTestSqr =  alphaTwoThird * ( dx(i) *dy(j) *dz(k) )**twoThird

      CASE( TURB_FWIDTH_GEOM )
        deltaTestSqr = ( alphaX* dx(i) )**2 &
                     + ( alphaY* dy(j) )**2 &
                     + ( alphaZ* dz(k) )**2
    END SELECT ! filterWidthModel 

   END SUBROUTINE TURB_CalcFilterWidth
!------------------------------------------------------------------------------
