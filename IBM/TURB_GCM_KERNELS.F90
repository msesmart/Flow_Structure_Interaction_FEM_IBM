!******************************************************************************
!
! Purpose: compute that the grid-filtered strain rates for the ghost cells
!          near an immersed boundary.
!
! Description: none.
!
! Input: field variables
!
! Output: test-filtered field variables
!
! Notes: Routines pertinent to GCM formalism
!
!******************************************************************************
!
! $Id: Exp $
!
! Copyright: (c) 2004 by the George Washington University
!
!******************************************************************************

!------------------------------------------------------------------------------
   SUBROUTINE TURB_GCM_CalcStrainRate

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE turb_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE turb_arrays
 
    IMPLICIT NONE

!... Parameters

!... Loop variables

    INTEGER :: iVar, n
     
!... Local variables

    INTEGER :: iG,jG,kG,nVarSij
    REAL(KIND=CGREAL) :: weightSij

!******************************************************************************

    nVarSij = 6

!------------------------------------------------------------------------------
! Compute strain rate in ghost cells (GC) based on inverse weighted method
!------------------------------------------------------------------------------
      
    DO n = 1, nGhost
      iG=iGhost(n)
      jG=jGhost(n)
      kG=kGhost(n)
      
      weightSij = (oned-iblank(iG-1,jG  ,kG  ))/dx(iG-1) &
        + (oned-iblank(iG+1,jG  ,kG  ))/dx(iG+1) &
		+ (oned-iblank(iG  ,jG-1,kG  ))/dy(jG-1) &
        + (oned-iblank(iG  ,jG+1,kG  ))/dy(jG+1) &
		+ (oned-iblank(iG  ,jG  ,kG-1))/dz(kG-1) &
        + (oned-iblank(iG  ,jG  ,kG+1))/dz(kG+1)
      
      DO iVar = 1, nVarSij
        sij(iVar,iG,jG,kG) = (oned-iblank(iG-1,jG  ,kG  ))/dx(iG-1)*sij(iVar,iG-1,jG  ,kG  ) &
	                   + (oned-iblank(iG+1,jG  ,kG  ))/dx(iG+1)*sij(iVar,iG+1,jG  ,kG  ) &
			   + (oned-iblank(iG  ,jG-1,kG  ))/dy(jG-1)*sij(iVar,iG  ,jG-1,kG  ) &
	                   + (oned-iblank(iG  ,jG+1,kG  ))/dy(jG+1)*sij(iVar,iG  ,jG+1,kG  ) &
			   + (oned-iblank(iG  ,jG  ,kG-1))/dz(kG-1)*sij(iVar,iG  ,jG  ,kG-1) &
	                   + (oned-iblank(iG  ,jG  ,kG+1))/dz(kG+1)*sij(iVar,iG  ,jG  ,kG+1) 
                  
        sij(iVar,iG,jG,kG) = sij(iVar,iG,jG,kG) *weightSij 
      ENDDO ! iVar 
    ENDDO   ! n
       
    END SUBROUTINE TURB_GCM_CalcStrainRate
   
!------------------------------------------------------------------------------



