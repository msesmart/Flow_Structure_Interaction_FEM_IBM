!******************************************************************************
!
! Purpose: Compute velocity components and normal pressure gradient at internal 
!          body intercept points. These are needed as BC in solvers.
!
! Description: none.
!
! Input: [u,v,w]BodyMarker  = velocity components of BM points
!        closestMarker      = closest Marker for Ghost point,
!        closestMarkerRatio = closest Marker ratio for Ghost point
!
! Output: [u,v,w]BodyIntercept = velocity components of BI point,
!         dpdnBodyIntercept    = normal pressure gradient of BI point,
!         dpdtBodyIntercept    = tangential pressure gradient of BI point.
!
! Notes: Currently dpdnBodyIntercept and dpdtBodyIntercept are set to zero.
!
!******************************************************************************
!
! $Id: Exp $
!
! Copyright: (c) 2003 by the George Washington University
!
!******************************************************************************

  SUBROUTINE GCM_SetBodyInterceptValues()

!----------------------------------------------------------------------
! Compute velocity components and normal pressure gradient at internal 
! body intercept points. These are needed as BC in solvers
!----------------------------------------------------------------------
!

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE pressure_arrays
    USE grid_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

!... Loop variables

    INTEGER :: n, k

!... Local variables

    INTEGER           :: iG,jG,kG
    REAL(KIND=CGREAL) :: xGC, yGC, zGC

!******************************************************************************

! Loop over all ghost points --------------------------------------------------

!  compute velocity field

    DO n=1,nGhost
      iG      = iGhost(n)
      jG      = jGhost(n)
      kG      = kGhost(n)
      xGC     = xc(iG)
      yGC     = yc(jG)
      zGC     = zc(kG)
      
!     CALL GCM_calc_BIVelocity2D( iG, jG, closestMarker(n), &
!                                 closestMarkerRatio(n),    &
!                                 uBodyIntercept(n),        &
!                                 vBodyIntercept(n),        &
!                                 wBodyIntercept(n)         )

      CALL GCM_calc_BIVelocity_Unstruc( iG, jG, kG, xGC, yGC, zGC, closestElementGC(n), &
                                       uBodyIntercept(n),        &
                                       vBodyIntercept(n),        &
                                       wBodyIntercept(n)         )

    ENDDO ! n

    IF (MOD(ntime,nmonitor) == 0 ) THEN 
      PRINT*, 'Min-Max uBI-Actual= ',MINVAL(uBodyIntercept(1:nGhost)),&
                                     MAXVAL(uBodyIntercept(1:nGhost))
      PRINT*, 'Min-Max vBI-Actual= ',MINVAL(vBodyIntercept(1:nGhost)),&
                                     MAXVAL(vBodyIntercept(1:nGhost))
      PRINT*, 'Min-Max wBI-Actual= ',MINVAL(wBodyIntercept(1:nGhost)),&
                                     MAXVAL(wBodyIntercept(1:nGhost))

      PRINT*, 'Min-Max xBINorm-Actual= ',MINVAL(xBodyInterceptNorm(1:nGhost)),&
                                         MAXVAL(xBodyInterceptNorm(1:nGhost))
      PRINT*, 'Min-Max yBINorm-Actual= ',MINVAL(yBodyInterceptNorm(1:nGhost)),&
                                         MAXVAL(yBodyInterceptNorm(1:nGhost))
    ENDIF ! ntime
 
!
! compute body intercept values for pressure field
!    
!      p   - p   = (dP/dn)  *  (total probe length)/2
!       gp    bi          bi           
!

    DO n = 1,nGhost
      iG = iGhost(n)
      jG = jGhost(n)

      DO k = 1,nz-1
!
!       IF (pressure_bc == INVISCID) THEN
!         dpdnBodyIntercept(n) = probeLength(n) *    
!                            ( -axBodyIntercept(n)*xBodyInterceptNorm(n) &  
!                              -ayBodyIntercept(n)*yBodyInterceptNorm(n) &
!                              -azBodyIntercept(n)*zBodyInterceptNorm(n) )
!       ELSE

          dpdnBodyIntercept(n) = zero

!       ENDIF ! pressure_bc 

! compute dp/dtau (tangential component)
!           dpdtBodyIntercept(n) = zero

      ENDDO ! k 
    ENDDO ! n
    
  END SUBROUTINE GCM_SetBodyInterceptValues      
!------------------------------------------------------------------------------

!******************************************************************************
!
! Purpose: Post compute following quantities following at body intercept points 
!          p        - for computing pressure force
!          u_x, u_y - for checking interpolation scheme
!          u_n, u_t - for checking interpolation scheme
!          du_t/dn  - form computing shear stress
!
! Description: none.
!
! Input: [u,v,w]  = velocity components,
!        p        = pressure.
!
! Output: pBodyIntercept = pressure at BI points.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: Exp $
!
! Copyright: (c) 2003 by the George Washington University
!
!******************************************************************************

  SUBROUTINE GCM_PostComputeBIValues()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE pressure_arrays
    USE grid_arrays
    USE GCM_arrays
    
    IMPLICIT NONE

!... Loop variables

    INTEGER :: iBody,iRow,k,n, kG

!... Local variables

    INTEGER :: iCM,iCMs,ind,iG,ii,jG,jj,i,j,kk 
    
    REAL(KIND=CGREAL) :: cMR, coeffBI, coeffBI_Inv, ratio_IPBI,  &
                         uIP, vIP, wIP, uBI, vBI, wBI,           &
                         uBIRhs, vBIRhs, wBIRhs,                 &
                         u_tGP, u_tIP, bodyInterceptWeightS, pIP

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: u_nBIPost, u_tBIPost, &
                                                    wBIPost, du_tdnBI
    
!******************************************************************************

! allocate local arrays -------------------------------------------------------

    ALLOCATE( u_nBIPost(nGhost) )
    ALLOCATE( u_tBIPost(nGhost) )
    ALLOCATE(   wBIPost(nGhost) )
    ALLOCATE( du_tdnBI(nGhost)  )
    
! initialize local arrays

    u_nBIPost = zero
    u_tBIPost = zero
    wBIPost   = zero
    du_tdnBI  = zero

    ratio_IPBI = imagePointWeight/bodyInterceptWeight
    
! compute pressure at body intercept points
!    
!      p   - p   = (dP/dn)  *  (DS)
!       gp    bi          bi       bi-gp     
!

    DO n = 1,nGhost
      iG = iGhost(n)
      jG = jGhost(n)
      DO k = 1,nz-1
        pBodyIntercept(n) = p(iG,jG,k)  - probeLength(n)*ratio_IPBI * dpdnBodyIntercept(n)
      END DO ! k
    END DO ! n 

!-----------------------------new 

! Compute pressure based on coefficient matrix from shear stress

    DO n = 1,nGhost
      iG = iGhost(n)
      jG = jGhost(n)
      kG = kGhost(n)

! Initialize values

        pIP = zero

! Compute pressure field at Image points

        DO iRow = 1, iRowMax
          ii = iCellIndexS(n) + incI(iRow)
          jj = jCellIndexS(n) + incJ(iRow)
!         kk = kCellIndexS(n) + incK(iRow)
          kk = 1

          IF ( ii /= iG .OR. jj /= jG .OR. kk /= kG ) THEN
            pIP = pIP + coeffGCMNS(iRow,n)*p(ii,jj,kk)
          ELSE
            pIP = pIP + coeffGCMNS(iRow,n)*dpdnBodyIntercept(n)
          ENDIF ! ii
        ENDDO ! iRow

        pBodyIntercept(n) = pIP  + probeLength(n)*(oned-ratio_IPBI) &
                                 * dpdnBodyIntercept(n)
    ENDDO   ! n

!-----------------------------new 
 
! post compute velocity (u_n, u_t)  and du_t/dn  at body intercept points
    
!     
!   U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)
!    gp                       ip      b
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i
!
!
    DO n = 1,nGhost
      iG = iGhost(n)
      jG = jGhost(n)
      kG = kGhost(n)

! Initialize values

        uIP = zero
        vIP = zero
        wIP = zero

        coeffBI = zero

! Compute velocity field at Image points

        DO iRow = 1, iRowMax
          ii = iCellIndexS(n) + incI(iRow) 
          jj = jCellIndexS(n) + incJ(iRow)
!         kk = kCellIndexS(n) + incK(iRow) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          kk =  1

          IF ( ii == iG .AND. jj == jG .AND. kk == kG ) THEN
            coeffBI = coeffGCMDS(iRow,n)*imagePointWeightS(n)       
          ELSE
            uIP = uIP + coeffGCMDS(iRow,n)* u(ii,jj,kk)
            vIP = vIP + coeffGCMDS(iRow,n)* v(ii,jj,kk) 
            wIP = wIP + coeffGCMDS(iRow,n)* w(ii,jj,kk)
          ENDIF ! ii          

        ENDDO ! iRow

! Extract values for Body Intercept points

        bodyInterceptWeightS = probeLengthNormalizedS(n)/&
                              ( probeLengthNormalizedS(n)-oned) 
        coeffBI_Inv = oned/(-coeffBI+bodyInterceptWeightS)

        uBI = coeffBI_Inv * ( u(iG,jG,kG) + uIP*imagePointWeightS(n) )
        vBI = coeffBI_Inv * ( v(iG,jG,kG) + vIP*imagePointWeightS(n) )
        wBI = coeffBI_Inv * ( w(iG,jG,kG) + wIP*imagePointWeightS(n) )

! Note :  e_n = e_z X e_t
! Note :  e_t = e_n X e_z
! Note :  e_z = e_t X e_n
!  =>   tx = ny; ty = -nx

! u_t = uIP*n_y - vIP*n_x

        u_nBIPost(n) = uBI*xBodyInterceptNorm(n) &
                      +vBI*yBodyInterceptNorm(n)
        u_tBIPost(n) = uBI*yBodyInterceptNorm(n) &
                      -vBI*yBodyInterceptNorm(n)
        wBIPost(n)   = wBI

        u_tIP        = uIP*yBodyInterceptNorm(n) &
                      -vIP*xBodyInterceptNorm(n)
        u_tGP        = u(iG,jG,k)*yBodyInterceptNorm(n) &
                      -v(iG,jG,k)*xBodyInterceptNorm(n)

        du_tdnBI(n) = ( u_tIP - u_tGP )/probeLengthS(n)

    END DO ! n
      
    IF (MOD(ntime,nmonitor) == 0 ) THEN 
      PRINT*, 'Min-Max uNBI*-POST= ',MINVAL(u_nBIPost(1:nGhost)),&
                                     MAXVAL(u_nBIPost(1:nGhost))
      PRINT*, 'Min-Max uTau*-POST= ',MINVAL(u_tBIPost(1:nGhost)),&
                                     MAXVAL(u_tBIPost(1:nGhost))
      PRINT*, 'Min-Max wBI*-POST = ',MINVAL(wBIPost(1:nGhost)),&
                                     MAXVAL(wBIPost(1:nGhost))                                  

      PRINT*, 'Min-Max duTau_dnBI*-POST= ',MINVAL(du_tdnBI(1:nGhost)),&
                                           MAXVAL(du_tdnBI(1:nGhost))
    ENDIF ! ntime

! deallocate local arrays

    DEALLOCATE(u_nBIPost)
    DEALLOCATE(u_tBIPost) 
    DEALLOCATE(wBIPost)
    DEALLOCATE(du_tdnBI)

  END SUBROUTINE GCM_PostComputeBIValues
!----------------------------------------------------------------------
