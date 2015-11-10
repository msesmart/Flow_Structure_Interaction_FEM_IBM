!-------------------------------------------     
   SUBROUTINE FreshCell_CalcExpWeight()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k


!    Fresh Cells are created in Move_Boundary.
!... Adjust weight for NL term evaluation
!    AB2 for normal cell
!    FE  for cell that was fresh in the previous time step.

      DO k=1,nzc
      DO j=1,nyc
      DO i=1,nxc
           exp_weight(i,j,k) = half* (REAL(3-fresh_cell(i,j,k),KIND=CGREAL))
      ENDDO
      ENDDO
      ENDDO
       
   END SUBROUTINE FreshCell_CalcExpWeight

!-------------------------------------------------------------------------------
   SUBROUTINE FreshCell_UpdateRhs() 

!  Update Advection-Diffusion RHS for fresh cell

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE nlold_arrays

    IMPLICIT NONE

    INTEGER              :: i,j,k

    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      nlu(i,j,k) = nlu(i,j,k)*(REAL(1-fresh_cell(i,j,k),KIND=CGREAL))
      nlv(i,j,k) = nlv(i,j,k)*(REAL(1-fresh_cell(i,j,k),KIND=CGREAL))
      nlw(i,j,k) = nlw(i,j,k)*(REAL(1-fresh_cell(i,j,k),KIND=CGREAL))
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

    IF (boundary_motion_type(1) == FEA_FLOW_STRUC_INTERACTION .OR.  &
        boundary_motion_type(1) == PARTIAL_DYNAMICS_COUPLED .OR.    &
        boundary_motion_type(1) == DYNAMICS_COUPLED .OR. &
        boundary_motion_type(1) == BIO_DYNAMICS_COUPLED) THEN
      nlu_FSI = nlu        !Added by Wanh
      nlv_FSI = nlv        !Added by Wanh
      nlw_FSI = nlw        !Added by Wanh
    ENDIF

   END SUBROUTINE  FreshCell_UpdateRhs
!-------------------------------------------------------------------------------
 
