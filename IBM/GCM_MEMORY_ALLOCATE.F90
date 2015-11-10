!------------------------------------------------------------------------
   SUBROUTINE BOUNDARY_allocate_memory()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays
    USE GCM_arrays

    IMPLICIT NONE
    INTEGER ::iBody
! arrays required by all internal boundaries

    ALLOCATE(       iblank(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(       bcBlank(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  iblank_memb(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  gateTest(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  iblank_solid(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  ghostCellMemb(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  ghostCellSolid(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  conflictCell(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  conflictBCi(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  conflictBCj(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(  conflictBCk(0:nx+1,0:ny+1,0:nz+1))

    ALLOCATE(   fresh_cell(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(   exp_weight(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(ghostCellMark(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(hybridMarkMemb(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(      bodyNum(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(      zoneCheck(0:nx+1,0:ny+1,0:nz+1))   !added by Chengyu
    ALLOCATE(          ivc(0:nx+1,0:ny+1,0:nz+1))

    ALLOCATE(num_iblank(0:nBody))

!    ALLOCATE(    boundCell(0:nx+1,0:ny+1,0:nz+1))
    exp_weight    = zero
    ghostCellMark = 0
    bodyNum       = 0
    zoneCheck     = 0   !addedy by Chengyu
!    boundCell     = 0
	ivc=0
    iblank_memb = 0
    iblank_solid = 0
    ghostCellMemb = 0

    ALLOCATE(iup(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(ium(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(jup(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(jum(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(kup(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(kum(0:nx+1,0:ny+1,0:nz+1))

    ALLOCATE(iMarkp(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(iMarkm(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(jMarkp(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(jMarkm(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(kMarkp(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(kMarkm(0:nx+1,0:ny+1,0:nz+1))
    iup = 0
    ium = 0
    jup = 0
    jum = 0
    kup = 0
    kum = 0

    iMarkp = 0
    iMarkm = 0
    jMarkp = 0
    jMarkm = 0
    kMarkp = 0
    kMarkm = 0

    ALLOCATE(pot_flag(0:nx+1,0:ny+1,0:nz+1))
    pot_flag = zero

    IF (Hybrid) THEN
      ALLOCATE(iupp(0:nx+1,0:ny+1,0:nz+1))  !
      ALLOCATE(iumm(0:nx+1,0:ny+1,0:nz+1))  !
      ALLOCATE(jupp(0:nx+1,0:ny+1,0:nz+1))  ! Added by Rupesh (used for 2nd Upwinding)
      ALLOCATE(jumm(0:nx+1,0:ny+1,0:nz+1))  !
      ALLOCATE(kupp(0:nx+1,0:ny+1,0:nz+1))  !
      ALLOCATE(kumm(0:nx+1,0:ny+1,0:nz+1))  !
    ENDIF

! since elliptic & general cylinder will be converted into
! 3D unstruc surfaces, we need to determine memory requirement for these
    DO iBody = 1, nBody

      SELECT CASE (canonical_body_type(iBody))

        CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)
           nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
           nPtsBodyMarker(iBody)    = nPtsBodyMarkerOrig(iBody)*nz
           totNumTriElem(iBody)     = 2*nPtsBodyMarkerOrig(iBody)*(nz-1)

        CASE(ELLIPSOID)

          IF (boundary_formulation == GCM_METHOD) THEN
           nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
           totNumTriElem(iBody)     = 2*nPtsBodyMarker(iBody) + 5

          ELSE
           nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
           totNumTriElem(iBody)     = 1
          END IF

        CASE(UNSTRUCTURED_SURFACE)
           nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)

      END SELECT ! canonical_body_type

    ENDDO ! iBody

    CALL MARKER_allocate_memory()
    CALL UNSTRUC_allocate_memory()

    ALLOCATE(iGhostP(1:nx*ny))   !Changed dimension of i,j,kGhost on 05/21/10
    ALLOCATE(jGhostP(1:nx*ny))
    ALLOCATE(kGhostP(1:nx*ny))

    IF (gcmFLAG == 1) CALL GCM_allocate_static_arrays
    if (pressure_osc_velocity.or.pressure_osc_pressure) call hybrid_mem_allocate_static

   END SUBROUTINE BOUNDARY_allocate_memory
!------------------------------------------------------------------------
!------------------------------------------------------------------------
   SUBROUTINE GCM_AllocateGhostCellArrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER :: nFaceMax

    iRowMax  = 8
    nFaceMax = 6

! Allocate Memory for various arrays pertinent to GCM Ghost Cell arrays

    ALLOCATE(iGhost(1:nGhost))
    ALLOCATE(jGhost(1:nGhost))
    ALLOCATE(kGhost(1:nGhost))

    ALLOCATE(coeffGCMD(iRowMax,nGhost))
    ALLOCATE(coeffGCMN(iRowMax,nGhost))

    ALLOCATE(xBodyInterceptTang(nGhost))
    ALLOCATE(yBodyInterceptTang(nGhost))
    ALLOCATE(zBodyInterceptTang(nGhost))

    ALLOCATE(xBodyInterceptNorm(nGhost))
    ALLOCATE(yBodyInterceptNorm(nGhost))
    ALLOCATE(zBodyInterceptNorm(nGhost))

    ALLOCATE(xBodyIntercept(nGhost))
    ALLOCATE(yBodyIntercept(nGhost))
    ALLOCATE(zBodyIntercept(nGhost))

    ALLOCATE(uBodyIntercept(nGhost))
    ALLOCATE(vBodyIntercept(nGhost))
    ALLOCATE(wBodyIntercept(nGhost))
    ALLOCATE(pBodyIntercept(nGhost))
    ALLOCATE(dpdnBodyIntercept(nGhost))
    ALLOCATE(dpdtBodyIntercept(nGhost))

    ALLOCATE(closestMarker(nGhost))
    ALLOCATE(closestMarkerRatio(nGhost))
    ALLOCATE(closestElementGC(nGhost))

    ALLOCATE(iCellIndex(nGhost))
    ALLOCATE(jCellIndex(nGhost))
    ALLOCATE(kCellIndex(nGhost))

    ALLOCATE(xImagePoint(nGhost))
    ALLOCATE(yImagePoint(nGhost))
    ALLOCATE(zImagePoint(nGhost))
    ALLOCATE(probeLength(nGhost))

! Allocate infrastructure for shear stress

    ALLOCATE(coeffGCMDS(iRowMax,nGhost))
    ALLOCATE(coeffGCMNS(iRowMax,nGhost))

    ALLOCATE(iCellIndexS(nGhost))
    ALLOCATE(jCellIndexS(nGhost))
    ALLOCATE(kCellIndexS(nGhost))

    ALLOCATE(xImagePointS(nGhost))
    ALLOCATE(yImagePointS(nGhost))
    ALLOCATE(zImagePointS(nGhost))
    ALLOCATE(probeLengthS(nGhost))
    ALLOCATE(probeLengthNormalizedS(nGhost))
    ALLOCATE(imagePointWeightS(nGhost))
!   Initialize arrays
    iGhost             = 0
    jGhost             = 0
    kGhost             = 0

    coeffGCMD          = zero
    coeffGCMN          = zero

    xBodyInterceptTang = zero
    yBodyInterceptTang = zero
    zBodyInterceptTang = zero

    xBodyInterceptNorm = zero
    yBodyInterceptNorm = zero
    zBodyInterceptNorm = zero

    xBodyIntercept     = zero
    yBodyIntercept     = zero
    zBodyIntercept     = zero

    uBodyIntercept     = zero
    vBodyIntercept     = zero
    wBodyIntercept     = zero

    pBodyIntercept     = zero
    dpdnBodyIntercept  = zero
    dpdtBodyIntercept  = zero

    xImagePoint        = zero
    yImagePoint        = zero
    zImagePoint        = zero
    probeLength        = zero

    closestMarker      = 0
    closestMarkerRatio = zero
    closestElementGC   = zero

    coeffGCMDS         = zero
    coeffGCMNS         = zero

    xImagePointS       = zero
    yImagePointS       = zero
    zImagePointS       = zero

    probeLengthS           = zero
    probeLengthNormalizedS = zero
    imagePointWeightS      = zero

   END SUBROUTINE GCM_AllocateGhostCellArrays
!------------------------------------------------------------------------
!------------------------------------------------------------------------
   SUBROUTINE GCM_AllocateFreshCellArrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

! Allocate Memory for various arrays pertinent to GCM Fresh Cells

    ALLOCATE(iFresh(nFresh))
    ALLOCATE(jFresh(nFresh))
    ALLOCATE(kFresh(nFresh))
    ALLOCATE(closestMarkerFresh(nFresh))
    ALLOCATE(closestElementFresh(nFresh))
    ALLOCATE(closestMarkerRatioFresh(nFresh))
    ALLOCATE(xBodyInterceptFresh(nFresh))
    ALLOCATE(yBodyInterceptFresh(nFresh))
    ALLOCATE(zBodyInterceptFresh(nFresh))
    ALLOCATE(iFreshCellIndex(nFresh))
    ALLOCATE(jFreshCellIndex(nFresh))
    ALLOCATE(kFreshCellIndex(nFresh))
    ALLOCATE(coeffGCMFreshD(iRowMax,nFresh))
    ALLOCATE(uBodyInterceptFresh(nFresh))
    ALLOCATE(vBodyInterceptFresh(nFresh))
    ALLOCATE(wBodyInterceptFresh(nFresh))

!   Initialize arrays
    iFresh                   = 0
    jFresh                   = 0
    kFresh                   = 0
    closestMarkerFresh       = 0
    closestElementFresh      = 0
    closestMarkerRatioFresh  = zero
    xBodyInterceptFresh      = zero
    yBodyInterceptFresh      = zero
    zBodyInterceptFresh      = zero
    iFreshCellIndex          = 0
    jFreshCellIndex          = 0
    kFreshCellIndex          = 0
    coeffGCMFreshD           = zero
    uBodyInterceptFresh      = zero
    vBodyInterceptFresh      = zero
    wBodyInterceptFresh      = zero
   END SUBROUTINE GCM_AllocateFreshCellArrays

!------------------------------------------------------------------------
   SUBROUTINE MARKER_allocate_memory()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays
    USE GCM_arrays
    USE usr_module, ONLY : uBodyMarker_iter,vBodyMarker_iter,wBodyMarker_iter, &         !Added by Wanh for FSI
                           uBodyMarker_worelax,vBodyMarker_worelax,wBodyMarker_worelax   !Added by Wanh for FSI
    USE fea_unstructure_surface
    USE body_dynamics

    IMPLICIT NONE

    INTEGER :: nBodyPtsMax,iBody
    INTEGER :: I
! Marker point arrays
    nPtsMax = 1
    IF ( internal_boundary_present == INTR_BOUND_PRESENT .AND. body_type == CANONICAL) THEN
      nPtsMax = MAXVAL(nPtsBodyMarker(:))
    ELSE
      nPtsMax = 1
    ENDIF
    PRINT*,' GCM_MEMORY_ALLOCATE: nPtsMax = ',nPtsMax

    ALLOCATE(xBodyMarker(nBody,nPtsMax))
    ALLOCATE(yBodyMarker(nBody,nPtsMax))
    ALLOCATE(zBodyMarker(nBody,nPtsMax))
    ALLOCATE(xBodyMarkerInit(nBody,nPtsMax))
    ALLOCATE(yBodyMarkerInit(nBody,nPtsMax))
    ALLOCATE(zBodyMarkerInit(nBody,nPtsMax))
    ALLOCATE(xBodyMarkerNonIner(nBody,nPtsMax))
    ALLOCATE(yBodyMarkerNonIner(nBody,nPtsMax))
    ALLOCATE(zBodyMarkerNonIner(nBody,nPtsMax))
    ALLOCATE(gateLabel(nBody,nPtsMax))
    ALLOCATE(zoneMarker(nBody,nPtsMax))  !added by Chengyu

    xBodyMarker     = zero
    yBodyMarker     = zero
    zBodyMarker     = zero
    gateLabel       =   0
    zoneMarker   = zero  !added by Chengyu

    ALLOCATE(sBodyMarker(nBody,nPtsMax+1))
    ALLOCATE(dsBodyMarker(nBody,nPtsMax+1))
    ALLOCATE(xNormBodyMarker(nBody,nPtsMax+1))
    ALLOCATE(yNormBodyMarker(nBody,nPtsMax+1))
    sBodyMarker     = zero
    dsBodyMarker    = zero
    xNormBodyMarker = zero
    yNormBodyMarker = zero

    ALLOCATE(uBodyMarker(nBody,nPtsMax))
    ALLOCATE(vBodyMarker(nBody,nPtsMax))
    ALLOCATE(wBodyMarker(nBody,nPtsMax))
    ALLOCATE(uBodyMarkerDyn(nBody,nPtsMax))
    ALLOCATE(vBodyMarkerDyn(nBody,nPtsMax))
    ALLOCATE(wBodyMarkerDyn(nBody,nPtsMax))

    ALLOCATE(uBodyMarkerDefor(nBody,nPtsMax)) ! Added by G. Liu
    ALLOCATE(vBodyMarkerDefor(nBody,nPtsMax))
    ALLOCATE(wBodyMarkerDefor(nBody,nPtsMax))
    ALLOCATE(uBodyMarkerDefor_prvs(nBody,nPtsMax)) ! Added by G. Liu
    ALLOCATE(vBodyMarkerDefor_prvs(nBody,nPtsMax))
    ALLOCATE(wBodyMarkerDefor_prvs(nBody,nPtsMax))


    uBodyMarker = zero
    vBodyMarker = zero
    wBodyMarker = zero
    uBodyMarkerDyn = zero
    vBodyMarkerDyn = zero
    wBodyMarkerDyn = zero

    uBodyMarkerDefor = zero ! Added by G. Liu
    vBodyMarkerDefor = zero
    wBodyMarkerDefor = zero
    uBodyMarkerDefor_prvs = zero ! Added by G. Liu
    vBodyMarkerDefor_prvs = zero
    wBodyMarkerDefor_prvs = zero

    ALLOCATE(axBodyMarker(nBody,nPtsMax))
    ALLOCATE(ayBodyMarker(nBody,nPtsMax))
    ALLOCATE(azBodyMarker(nBody,nPtsMax))
    axBodyMarker = zero
    ayBodyMarker = zero
    azBodyMarker = zero

    ! Added by Wanh for FSI
    IF (FSI_ON) THEN
       ALLOCATE(nPrescribedMarker(nBody))
       ALLOCATE(PrescribedMarker(nBody,nPtsMax/10))
       ALLOCATE(uBodyMarker_iter(nBody,nPtsMax))
       ALLOCATE(vBodyMarker_iter(nBody,nPtsMax))
       ALLOCATE(wBodyMarker_iter(nBody,nPtsMax))
       ALLOCATE(uBodyMarker_worelax(nBody,nPtsMax))
       ALLOCATE(vBodyMarker_worelax(nBody,nPtsMax))
       ALLOCATE(wBodyMarker_worelax(nBody,nPtsMax))
       ALLOCATE(deltau_Aitken(nPtsMax))
       ALLOCATE(deltav_Aitken(nPtsMax))
       ALLOCATE(deltaw_Aitken(nPtsMax))
       ALLOCATE(deltau_Aitken_prev(nPtsMax))    !Here _prev means previous iteration in a time-step
       ALLOCATE(deltav_Aitken_prev(nPtsMax))
       ALLOCATE(deltaw_Aitken_prev(nPtsMax))
    ENDIF

    ! Added by Wanh for Partial dynamic coupling
    DO i=1, nbody
    IF (boundary_motion_type(i)==PARTIAL_DYNAMICS_COUPLED) THEN
       ALLOCATE(DynamicMarker(nBody,nSection,nPtsMax))
       ALLOCATE(SectionMarker(nBody,nSection))
       ALLOCATE(nPrescribedMarker(nBody))
       ALLOCATE(PrescribedMarker(nBody,nPtsMax))
    ENDIF
    END DO

  END SUBROUTINE MARKER_allocate_memory

!-------------------------------------------------------------------
! Allocate Memory for various arrays pertinent to unstructured surface

   SUBROUTINE UNSTRUC_allocate_memory()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER :: nBodyPtsMax,nTriElemMax,iBody
    LOGICAL :: unstruc

    nTriElemMax = 0
    unstruc     = .FALSE.
    DO iBody = 1,nBody
      IF ( canonical_body_type(iBody) == ELLIPTIC_CYLINDER  .OR.   &
           canonical_body_type(iBody) == GENERAL_CYLINDER   .OR.   &
           canonical_body_type(iBody) == ELLIPSOID   .OR.   &
           canonical_body_type(iBody) == UNSTRUCTURED_SURFACE ) unstruc  = .TRUE.
    ENDDO
    IF ( unstruc .EQV. .TRUE. ) THEN
       nTriElemMax = MAXVAL(totNumTriElem(:))
       ALLOCATE(triElemNeig(nBody,3,nTriElemMax))
       ALLOCATE(triElemtang1X(nBody,nTriElemMax))
       ALLOCATE(triElemtang1Y(nBody,nTriElemMax))
       ALLOCATE(triElemtang1Z(nBody,nTriElemMax))
       ALLOCATE(triElemtang2X(nBody,nTriElemMax))
       ALLOCATE(triElemtang2Y(nBody,nTriElemMax))
       ALLOCATE(triElemtang2Z(nBody,nTriElemMax))
       ALLOCATE(triElemNormX(nBody,nTriElemMax))
       ALLOCATE(triElemNormY(nBody,nTriElemMax))
       ALLOCATE(triElemNormZ(nBody,nTriElemMax))
       ALLOCATE(triElemCentX(nBody,nTriElemMax))
       ALLOCATE(triElemCentY(nBody,nTriElemMax))
       ALLOCATE(triElemCentZ(nBody,nTriElemMax))
       ALLOCATE(triElemArea(nBody,nTriElemMax))
       ALLOCATE(pointOutsideBodyX(nBody))
       ALLOCATE(pointOutsideBodyY(nBody))
       ALLOCATE(pointOutsideBodyZ(nBody))
       ALLOCATE(surfArea(nBody))
    ENDIF
   END SUBROUTINE UNSTRUC_allocate_memory

!------------------------------------------------------------------------
! static arrays for GCM that need to be declared only once

   SUBROUTINE GCM_allocate_static_arrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays
    USE GCM_arrays

    IMPLICIT NONE

    iRowMax = 8
    ALLOCATE(incI(iRowMax))
    ALLOCATE(incJ(iRowMax))
    ALLOCATE(incK(iRowMax))
    ALLOCATE(iPvt(iRowMax))
    ALLOCATE(work(iRowMax))
    ALLOCATE(vanMatrixD(iRowMax,iRowMax))
    ALLOCATE(vanMatrixN(iRowMax,iRowMax))

! These allow us to define stencil image point
! Assumed clockwise from lower left corner.

      incI(1) = 0
      incJ(1) = 0
      incK(1) = 0
      incI(2) = 0
      incJ(2) = 1
      incK(2) = 0
      incI(3) = 1
      incJ(3) = 1
      incK(3) = 0
      incI(4) = 1
      incJ(4) = 0
      incK(4) = 0

      incI(5) = 0  !----?????????
      incJ(5) = 0
      incK(5) = 1
      incI(6) = 0
      incJ(6) = 1
      incK(6) = 1
      incI(7) = 1
      incJ(7) = 1
      incK(7) = 1
      incI(8) = 1
      incJ(8) = 0
      incK(8) = 1  ! ----????????

    END SUBROUTINE GCM_allocate_static_arrays
