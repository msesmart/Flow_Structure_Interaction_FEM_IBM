   SUBROUTINE fea_structure(iBody)

   USE fea_unstructure_surface
   USE unstructured_surface_arrays
   USE flow_parameters

   IMPLICIT NONE

   INTEGER :: iBody

   struc_disp_iter(:,:) = struc_disp(:,:)

   struc_disp(:,:) = 0

   CALL fea_stiff(iBody)

   IF (.not. TRNS) THEN
      CALL fea_static(iBody)

   ELSE IF (TRNS) THEN

!      IF (ntime == 1) THEN   ! Here restart is not taken care of yet
         DISP(:) = 0.0
         VEL(:) = 0.0
         ACC(:) = 0.0
         CALL FEA_MASS(iBody)
!      ENDIF


      CALL FEA_TRANSIENT(iBody)

   ENDIF

   END SUBROUTINE fea_structure

