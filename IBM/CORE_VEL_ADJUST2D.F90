!---------------------------------
!  SUBROUTINE vel_adjust2D() 
!---------------------------------

   SUBROUTINE vel_adjust2D() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays

    IMPLICIT NONE

    INTEGER              :: i,j,k

! Set Velocity field for 2D calculations

! copy k=1 plane to other planes
   
    DO k = 2,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1
      u(i,j,k) = u(i,j,1)
      v(i,j,k) = v(i,j,1)
      face_u(i,j,k) = face_u(i,j,1)
      face_v(i,j,k) = face_v(i,j,1)
      face_u(i+1,j,k) = face_u(i+1,j,1)
      face_v(i,j+1,k) = face_v(i,j+1,1)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
    
! zero w-component
   
    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1
      w(i,j,k) = zero
      face_w(i,j,k) = zero
      face_w(i,j,k+1) = zero
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

   END SUBROUTINE  vel_adjust2D
!-------------------------------------------------------------------------------
