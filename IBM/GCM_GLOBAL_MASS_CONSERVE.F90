!-------------------------------------------------------------------------

! compute mass flux at all boundaries and adjust outflow BC so as to satisfy
! global mass conservation.

   SUBROUTINE GCM_enforce_global_mass_consv()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER             :: i,j,k,n
    REAL(KIND=CGREAL)   :: massflux,correction_vel,fluxBody

    massflux = zero

    
    IF ( bcx1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcx2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcy1 == BC_TYPE_ZERO_GRADIENT .OR. & 
         bcy2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz2 == BC_TYPE_ZERO_GRADIENT       ) THEN

       CALL face_vel()

       DO k=1,nz-1
       DO j=1,ny-1
       DO i=1,nx-1
         massflux = massflux +                            &
                  ( -face_u(i,j,k)  *dy(j)*dz(k)   &
                    +face_u(i+1,j,k)*dy(j)*dz(k)   &
                    -face_v(i,j,k)  *dx(i)*dz(k)   &
                    +face_v(i,j+1,k)*dx(i)*dz(k)   &
                    -face_w(i,j,k)  *dx(i)*dy(j)   &
                    +face_w(i,j,k+1)*dx(i)*dy(j) ) &
                 *REAL(1-iblank(i,j,k),KIND=CGREAL)

       ENDDO ! i
       ENDDO ! j
       ENDDO ! k

! remove outer boundary condition correction before modifying BC 
       CALL remove_rhs_adjust_bc()

! adjust BC at outflow to satisfy global mass conservation
       correction_vel =-massflux/outflow_area  

!      IF ( MOD(ntime,nmonitor) == 0 ) THEN
!        PRINT*,' SET_BC:massflux       = ',massflux
!        PRINT*,' SET_BC:outflow_area   = ',outflow_area
!        PRINT*,' SET_BC:correction_vel = ',correction_vel
!      END IF ! ntime

       DO k=0,nz
       DO j=0,ny
         IF (bcx1 == BC_TYPE_ZERO_GRADIENT) bcxu(1,j,k)    = bcxu(1,j,k)    - correction_vel
         IF (bcx2 == BC_TYPE_ZERO_GRADIENT) bcxu(nx-1,j,k) = bcxu(nx-1,j,k) + correction_vel
       ENDDO ! j
       ENDDO ! k

       DO k=0,nz
       DO i=0,nx
         IF (bcy1 == BC_TYPE_ZERO_GRADIENT) bcyv(i,1,k)    = bcyv(i,1,k)    - correction_vel
         IF (bcy2 == BC_TYPE_ZERO_GRADIENT) bcyv(i,ny-1,k) = bcyv(i,ny-1,k) + correction_vel
       ENDDO ! i
       ENDDO ! k

       DO j=0,ny
       DO i=0,nx
         IF (bcz1 == BC_TYPE_ZERO_GRADIENT) bczw(i,j,1)    = bczw(i,j,1)    - correction_vel
         IF (bcz2 == BC_TYPE_ZERO_GRADIENT) bczw(i,j,nz-1) = bczw(i,j,nz-1) + correction_vel
       ENDDO ! i
       ENDDO ! j

! readjust outer boundary condition correction after modifying BC 
       CALL rhs_readjust_bc()

    ENDIF ! bcx1

   END SUBROUTINE GCM_enforce_global_mass_consv
!-------------------------------------------------------------------------

!-------------------------------------------------
   SUBROUTINE remove_rhs_adjust_bc()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER           :: i,j,k,n,iFr,jFr,iRow
    REAL(KIND=CGREAL) :: amxd,apxd,acxd
    REAL(KIND=CGREAL) :: amyd,apyd,acyd
    REAL(KIND=CGREAL) :: amzd,apzd,aczd
    REAL(KIND=CGREAL) :: rnDim

    rnDim = REAL((ndim - DIM_2D),KIND=CGREAL) 

    DO k=1,nz-1
    DO j=1,ny-1
    DO i=1,nx-1
       amxd = ( dxcinv(i)*(1-ium(i,j,k))            &
              +dxinv(i)*ium(i,j,k)*twod )*dxinv(i)
       apxd = ( dxcinv(i+1)*(1-iup(i,j,k))          &
              +dxinv(i)*iup(i,j,k)*twod )*dxinv(i)

       amyd = ( dycinv(j)*(1-jum(i,j,k))            &
              +dyinv(j)*jum(i,j,k)*twod )*dyinv(j)
       apyd = ( dycinv(j+1)*(1-jup(i,j,k))          &
              +dyinv(j)*jup(i,j,k)*twod )*dyinv(j)

       amzd = ( dzcinv(k)*(1-kum(i,j,k))            &
              +dzinv(k)*kum(i,j,k)*twod )*dzinv(k)
       apzd = ( dzcinv(k+1)*(1-kup(i,j,k))          &
              +dzinv(k)*kup(i,j,k)*twod )*dzinv(k)

       amxd =     - (0.50_CGREAL*dt*reinv)*amxd
       apxd =     - (0.50_CGREAL*dt*reinv)*apxd

       amyd =     - (0.50_CGREAL*dt*reinv)*amyd
       apyd =     - (0.50_CGREAL*dt*reinv)*apyd

       amzd =     - (0.50_CGREAL*dt*reinv)*amzd*rnDim 
       apzd =     - (0.50_CGREAL*dt*reinv)*apzd*rnDim 


       nlu(i,j,k) =  nlu(i,j,k)                 &
                     + amxd*bcxu(i,j,k)*ium(i,j,k) &
                     + apxd*bcxu(i,j,k)*iup(i,j,k) &
                     + amyd*bcyu(i,j,k)*jum(i,j,k) &
                     + apyd*bcyu(i,j,k)*jup(i,j,k) &
                     + amzd*bczu(i,j,k)*kum(i,j,k) &
                     + apzd*bczu(i,j,k)*kup(i,j,k)

        nlv(i,j,k) =  nlv(i,j,k)                 &
                     + amxd*bcxv(i,j,k)*ium(i,j,k) &
                     + apxd*bcxv(i,j,k)*iup(i,j,k) &
                     + amyd*bcyv(i,j,k)*jum(i,j,k) &
                     + apyd*bcyv(i,j,k)*jup(i,j,k) &
                     + amzd*bczv(i,j,k)*kum(i,j,k) &
                     + apzd*bczv(i,j,k)*kup(i,j,k)  

        nlw(i,j,k) =  nlw(i,j,k)                 &
                     + amxd*bcxw(i,j,k)*ium(i,j,k) &
                     + apxd*bcxw(i,j,k)*iup(i,j,k) &
                     + amyd*bcyw(i,j,k)*jum(i,j,k) &
                     + apyd*bcyw(i,j,k)*jup(i,j,k) &
                     + amzd*bczw(i,j,k)*kum(i,j,k) &
                     + apzd*bczw(i,j,k)*kup(i,j,k)



    ENDDO
    ENDDO
    ENDDO
 

  END SUBROUTINE remove_rhs_adjust_bc

 
!-------------------------------------------------
   SUBROUTINE rhs_readjust_bc()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER           :: i,j,k,n,iFr,jFr,iRow
    REAL(KIND=CGREAL) :: amxd,apxd,acxd
    REAL(KIND=CGREAL) :: amyd,apyd,acyd
    REAL(KIND=CGREAL) :: amzd,apzd,aczd
    REAL(KIND=CGREAL) :: rnDim

    rnDim = REAL((ndim - DIM_2D),KIND=CGREAL) 

    DO k=1,nz-1
    DO j=1,ny-1
    DO i=1,nx-1
       amxd = ( dxcinv(i)*(1-ium(i,j,k))            &
              +dxinv(i)*ium(i,j,k)*twod )*dxinv(i)
       apxd = ( dxcinv(i+1)*(1-iup(i,j,k))          &
              +dxinv(i)*iup(i,j,k)*twod )*dxinv(i)

       amyd = ( dycinv(j)*(1-jum(i,j,k))            &
              +dyinv(j)*jum(i,j,k)*twod )*dyinv(j)
       apyd = ( dycinv(j+1)*(1-jup(i,j,k))          &
              +dyinv(j)*jup(i,j,k)*twod )*dyinv(j)

       amzd = ( dzcinv(k)*(1-kum(i,j,k))            &
              +dzinv(k)*kum(i,j,k)*twod )*dzinv(k)
       apzd = ( dzcinv(k+1)*(1-kup(i,j,k))          &
              +dzinv(k)*kup(i,j,k)*twod )*dzinv(k)

       amxd =     - (0.50_CGREAL*dt*reinv)*amxd
       apxd =     - (0.50_CGREAL*dt*reinv)*apxd

       amyd =     - (0.50_CGREAL*dt*reinv)*amyd
       apyd =     - (0.50_CGREAL*dt*reinv)*apyd

       amzd =     - (0.50_CGREAL*dt*reinv)*amzd*rnDim 
       apzd =     - (0.50_CGREAL*dt*reinv)*apzd*rnDim 


       nlu(i,j,k) =  nlu(i,j,k)                 &
                     - amxd*bcxu(i,j,k)*ium(i,j,k) &
                     - apxd*bcxu(i,j,k)*iup(i,j,k) &
                     - amyd*bcyu(i,j,k)*jum(i,j,k) &
                     - apyd*bcyu(i,j,k)*jup(i,j,k) &
                     - amzd*bczu(i,j,k)*kum(i,j,k) &
                     - apzd*bczu(i,j,k)*kup(i,j,k)

        nlv(i,j,k) =  nlv(i,j,k)                 &
                     - amxd*bcxv(i,j,k)*ium(i,j,k) &
                     - apxd*bcxv(i,j,k)*iup(i,j,k) &
                     - amyd*bcyv(i,j,k)*jum(i,j,k) &
                     - apyd*bcyv(i,j,k)*jup(i,j,k) &
                     - amzd*bczv(i,j,k)*kum(i,j,k) &
                     - apzd*bczv(i,j,k)*kup(i,j,k)  

        nlw(i,j,k) =  nlw(i,j,k)                 &
                     - amxd*bcxw(i,j,k)*ium(i,j,k) &
                     - apxd*bcxw(i,j,k)*iup(i,j,k) &
                     - amyd*bcyw(i,j,k)*jum(i,j,k) &
                     - apyd*bcyw(i,j,k)*jup(i,j,k) &
                     - amzd*bczw(i,j,k)*kum(i,j,k) &
                     - apzd*bczw(i,j,k)*kup(i,j,k)



    ENDDO
    ENDDO
    ENDDO
 

  END SUBROUTINE rhs_readjust_bc

 
!-------------------------------------------------------------------------

! compute mass flux at all boundaries and adjust outflow BC so as to satisfy
! global mass conservation.

   SUBROUTINE GCM_enforce_p_compatibility(pres,mx,my,mz)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE pressure_arrays
    USE grid_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER :: mx, my, mz
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT(IN) :: pres   ! upper bound changed to nx from nx+1

    INTEGER             :: i,j,k,n
    REAL(KIND=CGREAL)   :: pgradflux,correction_pgrad
    REAL(KIND=CGREAL)   :: pgxw,pgxe,pgys,pgyn,pgzb,pgzf

    pgradflux = zero

    
    IF ( bcx1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcx2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcy1 == BC_TYPE_ZERO_GRADIENT .OR. & 
         bcy2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz2 == BC_TYPE_ZERO_GRADIENT       ) THEN

! add pressure gradient to face velocities
   DO k = 1,nz-1
    DO j = 1,ny-1
     DO i = 1,nx-1

     IF (bcx1 .EQ. BC_TYPE_PERIODIC .AND. &
           bcx2 .EQ. BC_TYPE_PERIODIC ) THEN
        pgxw     = (pres(i,j,k)  -pres(i-1,j,k))*dxcinv(i)
        pgxe     = (pres(i+1,j,k)-pres(i,j,k)  )*dxcinv(i+1)
     ELSE
        pgxw     = pgradx1(j,k)*ium(i,j,k) + (1 - ium(i,j,k))*(pres(i,j,k)  -pres(i-1,j,k))*dxcinv(i)
        pgxe     = pgradx2(j,k)*iup(i,j,k) + (1 - iup(i,j,k))*(pres(i+1,j,k)-pres(i,j,k)  )*dxcinv(i+1)
     END IF

     IF (bcy1 .EQ. BC_TYPE_PERIODIC .AND. &
          bcy2 .EQ. BC_TYPE_PERIODIC) THEN
        pgys     = (pres(i,j,k)  -pres(i,j-1,k))*dycinv(j)
        pgyn     = (pres(i,j+1,k)-pres(i,j,k)  )*dycinv(j+1)
     ELSE
        pgys     = pgrady1(i,k)*jum(i,j,k) + (1 - jum(i,j,k))*(pres(i,j,k)  -pres(i,j-1,k))*dycinv(j)
        pgyn     = pgrady2(i,k)*jup(i,j,k) + (1 - jup(i,j,k))*(pres(i,j+1,k)-pres(i,j,k)  )*dycinv(j+1)
     END IF

     IF (bcz1 .EQ. BC_TYPE_PERIODIC .AND. &
         bcz2 .EQ. BC_TYPE_PERIODIC ) THEN
        pgzb     = (pres(i,j,k)  -pres(i,j,k-1))*dzcinv(k)
        pgzf     = (pres(i,j,k+1)-pres(i,j,k)  )*dzcinv(k+1)
     ELSE
        pgzb     = pgradz1(i,j)*kum(i,j,k) + (1 - kum(i,j,k))*(pres(i,j,k)  -pres(i,j,k-1))*dzcinv(k)
        pgzf     = pgradz2(i,j)*kup(i,j,k) + (1 - kup(i,j,k))*(pres(i,j,k+1)-pres(i,j,k)  )*dzcinv(k+1)
     END IF

        pgradflux = pgradflux +                            &
                 ( -pgxw*dy(j)*dz(k)   &
                   +pgxe*dy(j)*dz(k)   &
                   -pgys*dx(i)*dz(k)   &
                   +pgyn*dx(i)*dz(k)   &
                   -pgzb*dx(i)*dy(j)   &
                   +pgzf*dx(i)*dy(j) ) &
                 *REAL(1-iblank(i,j,k),KIND=CGREAL)
      ENDDO
      ENDDO
      ENDDO

! adjust BC at outflow to satisfy global mass conservation
       correction_pgrad =-pgradflux/outflow_area  

!      IF ( MOD(ntime,nmonitor) == 0 ) THEN
!        PRINT*,' SET_BC:pgrad flux       = ',pgradflux
!        PRINT*,' SET_BC:outflow_area     = ',outflow_area
!        PRINT*,' SET_BC:correction_pgrad = ',correction_pgrad
!      END IF ! ntime

       DO k=0,nz
       DO j=0,ny
         IF (bcx1 == BC_TYPE_ZERO_GRADIENT) pgradx1(j,k) = pgradx1(j,k) - correction_pgrad
         IF (bcx2 == BC_TYPE_ZERO_GRADIENT) pgradx2(j,k) = pgradx2(j,k) + correction_pgrad
       ENDDO ! j
       ENDDO ! k

       DO k=0,nz
       DO i=0,nx
         IF (bcy1 == BC_TYPE_ZERO_GRADIENT) pgrady1(i,k) = pgrady1(i,k) - correction_pgrad 
         IF (bcy2 == BC_TYPE_ZERO_GRADIENT) pgrady2(i,k) = pgrady2(i,k) + correction_pgrad 
       ENDDO ! i
       ENDDO ! k

       DO j=0,ny
       DO i=0,nx
         IF (bcz1 == BC_TYPE_ZERO_GRADIENT) pgradz1(i,j) = pgradz1(i,j) - correction_pgrad 
         IF (bcz2 == BC_TYPE_ZERO_GRADIENT) pgradz2(i,j) = pgradz2(i,j) + correction_pgrad 
       ENDDO ! i
       ENDDO ! j

    ENDIF ! bcx1

!   pgradflux = zero

!     DO k = 1,nz-1
!     DO j = 1,ny-1
!     DO i = 1,nx-1
!       pgxw     = pgradx1(j,k)*ium(i,j,k) + (oned - ium(i,j,k))*(pres(i,j,k)  -pres(i-1,j,k))*dxcinv(i)
!       pgxe     = pgradx2(j,k)*iup(i,j,k) + (oned - iup(i,j,k))*(pres(i+1,j,k)-pres(i,j,k)  )*dxcinv(i+1)
!       pgys     = pgrady1(i,k)*jum(i,j,k) + (oned - jum(i,j,k))*(pres(i,j,k)  -pres(i,j-1,k))*dycinv(j)
!       pgyn     = pgrady2(i,k)*jup(i,j,k) + (oned - jup(i,j,k))*(pres(i,j+1,k)-pres(i,j,k)  )*dycinv(j+1)
!       pgzb     = pgradz1(i,j)*kum(i,j,k) + (oned - kum(i,j,k))*(pres(i,j,k)  -pres(i,j,k-1))*dzcinv(k)
!       pgzf     = pgradz2(i,j)*kup(i,j,k) + (oned - kup(i,j,k))*(pres(i,j,k+1)-pres(i,j,k)  )*dzcinv(k+1)

!       pgradflux = pgradflux +                            &
!                ( -pgxw*dy(j)*dz(k)   &
!                  +pgxe*dy(j)*dz(k)   &
!                  -pgys*dx(i)*dz(k)   &
!                  +pgyn*dx(i)*dz(k)   &
!                  -pgzb*dx(i)*dy(j)   &
!                  +pgzf*dx(i)*dy(j) ) &
!                *(oned-REAL(iblank(i,j,k),KIND=CGREAL))

!     ENDDO
!     ENDDO
!     ENDDO

!     IF ( MOD(ntime,nmonitor) == 0 ) THEN
!       PRINT*,' check:pgrad flux     = ',pgradflux
!     END IF ! ntime

   END SUBROUTINE GCM_enforce_p_compatibility
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

! compute mass flux at all boundaries and adjust outflow BC so as to satisfy
! global mass conservation.

   SUBROUTINE GCM_enforce_p_compatibility_outter(pres)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE pressure_arrays
    USE grid_arrays
    USE GCM_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(IN) :: pres   ! upper bound changed to nx from nx+1

    INTEGER             :: i,j,k,n
    REAL(KIND=CGREAL)   :: pgradflux,correction_pgrad
    REAL(KIND=CGREAL)   :: pgxw,pgxe,pgys,pgyn,pgzb,pgzf

    pgradflux = zero

    
    IF ( bcx1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcx2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcy1 == BC_TYPE_ZERO_GRADIENT .OR. & 
         bcy2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz2 == BC_TYPE_ZERO_GRADIENT       ) THEN

! add pressure gradient to face velocities
   DO k = 1,nz-1
    DO j = 1,ny-1
     DO i = 1,nx-1

     IF (bcx1 .EQ. BC_TYPE_PERIODIC .AND. &
           bcx2 .EQ. BC_TYPE_PERIODIC ) THEN
        pgxw     = (pres(i,j,k)  -pres(i-1,j,k))*dxcinv(i)
        pgxe     = (pres(i+1,j,k)-pres(i,j,k)  )*dxcinv(i+1)
     ELSE
        pgxw     = pgradx1(j,k)*ium(i,j,k) + (1 - ium(i,j,k))*(pres(i,j,k)  -pres(i-1,j,k))*dxcinv(i)
        pgxe     = pgradx2(j,k)*iup(i,j,k) + (1 - iup(i,j,k))*(pres(i+1,j,k)-pres(i,j,k)  )*dxcinv(i+1)
     END IF

     IF (bcy1 .EQ. BC_TYPE_PERIODIC .AND. &
          bcy2 .EQ. BC_TYPE_PERIODIC) THEN
        pgys     = (pres(i,j,k)  -pres(i,j-1,k))*dycinv(j)
        pgyn     = (pres(i,j+1,k)-pres(i,j,k)  )*dycinv(j+1)
     ELSE
        pgys     = pgrady1(i,k)*jum(i,j,k) + (1 - jum(i,j,k))*(pres(i,j,k)  -pres(i,j-1,k))*dycinv(j)
        pgyn     = pgrady2(i,k)*jup(i,j,k) + (1 - jup(i,j,k))*(pres(i,j+1,k)-pres(i,j,k)  )*dycinv(j+1)
     END IF

     IF (bcz1 .EQ. BC_TYPE_PERIODIC .AND. &
         bcz2 .EQ. BC_TYPE_PERIODIC ) THEN
        pgzb     = (pres(i,j,k)  -pres(i,j,k-1))*dzcinv(k)
        pgzf     = (pres(i,j,k+1)-pres(i,j,k)  )*dzcinv(k+1)
     ELSE
        pgzb     = pgradz1(i,j)*kum(i,j,k) + (1 - kum(i,j,k))*(pres(i,j,k)  -pres(i,j,k-1))*dzcinv(k)
        pgzf     = pgradz2(i,j)*kup(i,j,k) + (1 - kup(i,j,k))*(pres(i,j,k+1)-pres(i,j,k)  )*dzcinv(k+1)
     END IF

        pgradflux = pgradflux +                            &
                 ( -pgxw*dy(j)*dz(k)   &
                   +pgxe*dy(j)*dz(k)   &
                   -pgys*dx(i)*dz(k)   &
                   +pgyn*dx(i)*dz(k)   &
                   -pgzb*dx(i)*dy(j)   &
                   +pgzf*dx(i)*dy(j) ) &
                 *REAL(1-iblank(i,j,k),KIND=CGREAL)
      ENDDO
      ENDDO
      ENDDO

! adjust BC at outflow to satisfy global mass conservation
       correction_pgrad =-pgradflux/outflow_area  

!      IF ( MOD(ntime,nmonitor) == 0 ) THEN
!        PRINT*,' SET_BC:pgrad flux       = ',pgradflux
!        PRINT*,' SET_BC:outflow_area     = ',outflow_area
!        PRINT*,' SET_BC:correction_pgrad = ',correction_pgrad
!      END IF ! ntime

       DO k=0,nz
       DO j=0,ny
         IF (bcx1 == BC_TYPE_ZERO_GRADIENT) pgradx1(j,k) = pgradx1(j,k) - correction_pgrad
         IF (bcx2 == BC_TYPE_ZERO_GRADIENT) pgradx2(j,k) = pgradx2(j,k) + correction_pgrad
       ENDDO ! j
       ENDDO ! k

       DO k=0,nz
       DO i=0,nx
         IF (bcy1 == BC_TYPE_ZERO_GRADIENT) pgrady1(i,k) = pgrady1(i,k) - correction_pgrad 
         IF (bcy2 == BC_TYPE_ZERO_GRADIENT) pgrady2(i,k) = pgrady2(i,k) + correction_pgrad 
       ENDDO ! i
       ENDDO ! k

       DO j=0,ny
       DO i=0,nx
         IF (bcz1 == BC_TYPE_ZERO_GRADIENT) pgradz1(i,j) = pgradz1(i,j) - correction_pgrad 
         IF (bcz2 == BC_TYPE_ZERO_GRADIENT) pgradz2(i,j) = pgradz2(i,j) + correction_pgrad 
       ENDDO ! i
       ENDDO ! j

    ENDIF ! bcx1

!   pgradflux = zero

!     DO k = 1,nz-1
!     DO j = 1,ny-1
!     DO i = 1,nx-1
!       pgxw     = pgradx1(j,k)*ium(i,j,k) + (oned - ium(i,j,k))*(pres(i,j,k)  -pres(i-1,j,k))*dxcinv(i)
!       pgxe     = pgradx2(j,k)*iup(i,j,k) + (oned - iup(i,j,k))*(pres(i+1,j,k)-pres(i,j,k)  )*dxcinv(i+1)
!       pgys     = pgrady1(i,k)*jum(i,j,k) + (oned - jum(i,j,k))*(pres(i,j,k)  -pres(i,j-1,k))*dycinv(j)
!       pgyn     = pgrady2(i,k)*jup(i,j,k) + (oned - jup(i,j,k))*(pres(i,j+1,k)-pres(i,j,k)  )*dycinv(j+1)
!       pgzb     = pgradz1(i,j)*kum(i,j,k) + (oned - kum(i,j,k))*(pres(i,j,k)  -pres(i,j,k-1))*dzcinv(k)
!       pgzf     = pgradz2(i,j)*kup(i,j,k) + (oned - kup(i,j,k))*(pres(i,j,k+1)-pres(i,j,k)  )*dzcinv(k+1)

!       pgradflux = pgradflux +                            &
!                ( -pgxw*dy(j)*dz(k)   &
!                  +pgxe*dy(j)*dz(k)   &
!                  -pgys*dx(i)*dz(k)   &
!                  +pgyn*dx(i)*dz(k)   &
!                  -pgzb*dx(i)*dy(j)   &
!                  +pgzf*dx(i)*dy(j) ) &
!                *(oned-REAL(iblank(i,j,k),KIND=CGREAL))

!     ENDDO
!     ENDDO
!     ENDDO

!     IF ( MOD(ntime,nmonitor) == 0 ) THEN
!       PRINT*,' check:pgrad flux     = ',pgradflux
!     END IF ! ntime

   END SUBROUTINE GCM_enforce_p_compatibility_outter
!-------------------------------------------------------------------------

