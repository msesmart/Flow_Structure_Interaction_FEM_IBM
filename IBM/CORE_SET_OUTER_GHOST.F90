   SUBROUTINE set_outer_ghost_vel()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE blasius_profile


    IMPLICIT NONE

    INTEGER             :: i,j,k


! outer ghost boundary conditions

      DO k=0,nz
      DO j=0,ny
! left boundary
        i = 0
        u(i,j,k) = 2.0_CGREAL*bcxu(i+1,j,k) - u(i+1,j,k)
        v(i,j,k) = 2.0_CGREAL*bcxv(i+1,j,k) - v(i+1,j,k)
        w(i,j,k) = 2.0_CGREAL*bcxw(i+1,j,k) - w(i+1,j,k)
! right boundary
        i = nx
        u(i,j,k) = 2.0_CGREAL*bcxu(i-1,j,k) - u(i-1,j,k)
        v(i,j,k) = 2.0_CGREAL*bcxv(i-1,j,k) - v(i-1,j,k)
        w(i,j,k) = 2.0_CGREAL*bcxw(i-1,j,k) - w(i-1,j,k)
      ENDDO ! j
      ENDDO ! k

      DO k=0,nz
      DO i=0,nx
! bottom boundary
        j = 0
        u(i,j,k) = 2.0_CGREAL*bcyu(i,j+1,k) - u(i,j+1,k)
        v(i,j,k) = 2.0_CGREAL*bcyv(i,j+1,k) - v(i,j+1,k)
        w(i,j,k) = 2.0_CGREAL*bcyw(i,j+1,k) - w(i,j+1,k)
! top boundary
        j = ny
        u(i,j,k) = 2.0_CGREAL*bcyu(i,j-1,k) - u(i,j-1,k)
        v(i,j,k) = 2.0_CGREAL*bcyv(i,j-1,k) - v(i,j-1,k)
        w(i,j,k) = 2.0_CGREAL*bcyw(i,j-1,k) - w(i,j-1,k)
      ENDDO ! i
      ENDDO ! k

      DO j=0,ny
      DO i=0,nx
! back boundary
        k = 0
        u(i,j,k) = 2.0_CGREAL*bczu(i,j,k+1) - u(i,j,k+1)
        v(i,j,k) = 2.0_CGREAL*bczv(i,j,k+1) - v(i,j,k+1)
        w(i,j,k) = 2.0_CGREAL*bczw(i,j,k+1) - w(i,j,k+1)
! front boundary
        k = nz
        u(i,j,k) = 2.0_CGREAL*bczu(i,j,k-1) - u(i,j,k-1)
        v(i,j,k) = 2.0_CGREAL*bczv(i,j,k-1) - v(i,j,k-1)
        w(i,j,k) = 2.0_CGREAL*bczw(i,j,k-1) - w(i,j,k-1)
     ENDDO ! i
     ENDDO ! j

   END SUBROUTINE set_outer_ghost_vel
!-------------------------------------------------------------------------

   SUBROUTINE set_outer_ghost_pres(pres, mx, my, mz)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE blasius_profile


    IMPLICIT NONE

    INTEGER :: mx, my, mz
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT(INOUT) :: pres   ! upper bound changed to nx from nx+1


    INTEGER             :: i,j,k


! outer ghost boundary conditions

   IF (bcx1 == BC_TYPE_PERIODIC .AND. &
       bcx2 == BC_TYPE_PERIODIC) THEN
     DO k=0,nz
     DO j=0,ny
        pres(0,j,k)  = pres(nx-1,j,k)
        pres(nx,j,k) = pres(1,j,k)
     ENDDO ! j
     ENDDO ! k
   END IF

   IF (pbcx1 == PBC_NEUMANN) THEN
     DO k=0,nz
     DO j=0,ny
        pres(0,j,k)  = pres(1,j,k)
     ENDDO ! j
     ENDDO ! k
   ELSE
	DO j=0,ny
	DO k=0,nz
        pres(0,j,k)  = - pres(1,j,k)
	ENDDO
	ENDDO
   
   END IF
 
   IF (pbcx2 == PBC_NEUMANN) THEN
     DO k=0,nz
     DO j=0,ny
        pres(nx,j,k) = pres(nx-1,j,k)
     ENDDO ! j
     ENDDO ! k
   ELSE
     DO k=0,nz
     DO j=0,ny
        pres(nx,j,k) = - pres(nx-1,j,k)
     ENDDO ! j
     ENDDO ! k
   END IF


   IF (bcy1 .EQ. BC_TYPE_PERIODIC .AND. &
       bcy2 .EQ. BC_TYPE_PERIODIC) THEN
     DO k=0,nz
     DO i=0,nx
        pres(i,0,k)  = pres(i,ny-1,k)
        pres(i,ny,k) = pres(i,1,k)
     ENDDO ! i
     ENDDO ! k
   END IF

   IF (pbcy1 == PBC_NEUMANN ) THEN
     DO k=0,nz
     DO i=0,nx
        pres(i,0,k)  = pres(i,1,k)
     ENDDO ! i
     ENDDO ! k
   ELSE
     DO k=0,nz
     DO i=0,nx
        pres(i,0,k)  = - pres(i,1,k)
     ENDDO ! i
     ENDDO ! k
   END IF
 
   IF (pbcy2 == PBC_NEUMANN ) THEN
     DO k=0,nz
     DO i=0,nx
        pres(i,ny,k) = pres(i,ny-1,k)
     ENDDO ! i
     ENDDO ! k
   ELSE
     DO k=0,nz
     DO i=0,nx
        pres(i,ny,k) = - pres(i,ny-1,k)
     ENDDO ! i
     ENDDO ! k
   END IF


 IF (ndim .EQ. DIM_3D) THEN
   IF (bcz1 .EQ. BC_TYPE_PERIODIC .AND. &
       bcz2 .EQ. BC_TYPE_PERIODIC) THEN
     DO j=0,ny
     DO i=0,nx
        pres(i,j,0)  = pres(i,j,nz-1)
        pres(i,j,nz) = pres(i,j,1)
     ENDDO ! i
     ENDDO ! j
   ELSE
     DO j=0,ny
     DO i=0,nx
        pres(i,j,0)  = pres(i,j,1)
        pres(i,j,nz) = pres(i,j,nz-1)
     ENDDO ! i
     ENDDO ! j
   END IF
 END IF

   IF (pbcz1 == PBC_NEUMANN ) THEN
     DO j=0,ny
     DO i=0,nx
        pres(i,j,0)  = pres(i,j,1)
     ENDDO ! i
     ENDDO ! j
   ELSE
     DO j=0,ny
     DO i=0,nx
        pres(i,j,0)  = - pres(i,j,1)
     ENDDO ! i
     ENDDO ! j
   END IF
 
   IF (pbcz2 == PBC_NEUMANN ) THEN
     DO j=0,ny
     DO i=0,nx
        pres(i,j,nz) = pres(i,j,nz-1)
     ENDDO ! i
     ENDDO ! j
   ELSE
     DO j=0,ny
     DO i=0,nx
        pres(i,j,nz) = - pres(i,j,nz-1)
     ENDDO ! i
     ENDDO ! j
   END IF

   END SUBROUTINE set_outer_ghost_pres
!-------------------------------------------------------------------------
