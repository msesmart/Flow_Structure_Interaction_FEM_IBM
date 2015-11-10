!------------------------------------------------------------------------
   SUBROUTINE GCM_DeallocateGhostCellArrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE


! Deallocate Memory for various arrays pertinent to GCM Body Markers

    DEALLOCATE(iGhost)
    DEALLOCATE(jGhost) 
    DEALLOCATE(kGhost) 

!    DEALLOCATE(iGhostP)   !Added on 05/31/10
!    DEALLOCATE(jGhostP) 
!    DEALLOCATE(kGhostP) 

    DEALLOCATE(coeffGCMD)
    DEALLOCATE(coeffGCMN)
    
    DEALLOCATE(xBodyInterceptTang)
    DEALLOCATE(yBodyInterceptTang)
    DEALLOCATE(zBodyInterceptTang)

    DEALLOCATE(xBodyInterceptNorm)
    DEALLOCATE(yBodyInterceptNorm)
    DEALLOCATE(zBodyInterceptNorm)
    
    DEALLOCATE(xBodyIntercept)
    DEALLOCATE(yBodyIntercept)
    DEALLOCATE(zBodyIntercept)

    DEALLOCATE(uBodyIntercept)
    DEALLOCATE(vBodyIntercept) 
    DEALLOCATE(wBodyIntercept) 
    DEALLOCATE(pBodyIntercept) 
    DEALLOCATE(dpdnBodyIntercept) 
    DEALLOCATE(dpdtBodyIntercept) 
                        
    DEALLOCATE(closestMarker) 
    DEALLOCATE(closestMarkerRatio) 
    DEALLOCATE(closestElementGC)
          
    DEALLOCATE(iCellIndex)
    DEALLOCATE(jCellIndex)
    DEALLOCATE(kCellIndex)
 
    DEALLOCATE(xImagePoint)
    DEALLOCATE(yImagePoint) 
    DEALLOCATE(zImagePoint) 
    DEALLOCATE(probeLength)

! Deallocate infrastructure for shear stress

    DEALLOCATE(coeffGCMDS)
    DEALLOCATE(coeffGCMNS)

    DEALLOCATE(iCellIndexS)
    DEALLOCATE(jCellIndexS)    
    DEALLOCATE(kCellIndexS)    

    DEALLOCATE(xImagePointS)
    DEALLOCATE(yImagePointS) 
    DEALLOCATE(zImagePointS) 
    DEALLOCATE(probeLengthS)
    DEALLOCATE(probeLengthNormalizedS)
    DEALLOCATE(imagePointWeightS)
    
   END SUBROUTINE GCM_DeallocateGhostCellArrays
!------------------------------------------------------------------------
!------------------------------------------------------------------------
   SUBROUTINE GCM_DeallocateFreshCellArrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE

! Deallocate Memory for various arrays pertinent to GCM Fresh Cells

    DEALLOCATE(iFresh)
    DEALLOCATE(jFresh)
    DEALLOCATE(kFresh)
    DEALLOCATE(closestMarkerFresh)
    DEALLOCATE(closestElementFresh)
    DEALLOCATE(closestMarkerRatioFresh)
    DEALLOCATE(xBodyInterceptFresh)
    DEALLOCATE(yBodyInterceptFresh)
    DEALLOCATE(zBodyInterceptFresh)
    DEALLOCATE(iFreshCellIndex)
    DEALLOCATE(jFreshCellIndex)
    DEALLOCATE(kFreshCellIndex)
    DEALLOCATE(coeffGCMFreshD)
    DEALLOCATE(uBodyInterceptFresh)
    DEALLOCATE(vBodyInterceptFresh)
    DEALLOCATE(wBodyInterceptFresh)

   END SUBROUTINE GCM_DeallocateFreshCellArrays
!------------------------------------------------------------------------
