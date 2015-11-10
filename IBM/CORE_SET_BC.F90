!-------------------------------------------------
!  SUBROUTINE set_bc()
!-------------------------------------------------

!-------------------------------------------------
! bcflag = 1  => dirichlet bc
! bcflag = 2  => neumann bc
! bcflag = 3  => pulsatile bc
!
!     |-------|-------|-------|-------|-------
!     |       |       |       |       |
!     |   o   |   o   |   o   |   o   |   o
!     |       |  bcy  |  bcy  |       |
!     |-------|---+---|-------|-------|-------
!     |       |*******|*******|       |
!     |   obcx+*******|*******|bcxo   |   o
!     |       |*******|*******|       |
!     |-------|-------|-------|-------|-------
!
!
!-------------------------------------------------
SUBROUTINE set_bc()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE blasius_profile

    IMPLICIT NONE
    INTEGER             :: i,j,k

! Internal Boundaries Conditions
!  Need to set iup for GCM to satisfy mass conservation on Stair-step boundary
!        CALL write_dump_debug('bcxu',0,bcxu)

      IF (boundary_formulation == SSM_METHOD) CALL SSM_set_bc_internal()

      IF (boundary_formulation == GCM_METHOD) THEN
          CALL GCM_SetBodyInterceptValues()
          CALL GCM_SetBodyInterceptValuesFresh()
      ENDIF

! Outer boundary conditions
!
      IF (bcx1.eq.BC_TYPE_USER_SPECIFIED) THEN
           CALL blasius_velocity
      END IF

      DO k=0,nz
      DO j=0,ny

! left boundary
        i = 1
        SELECT CASE (bcx1)
          CASE (BC_TYPE_DIRICHLET)      ! dirichlet bc
            bcxu(i,j,k) = ux1
            bcxv(i,j,k) = vx1
            bcxw(i,j,k) = wx1
          CASE (BC_TYPE_ZERO_GRADIENT)  ! outflow bc ( zero gradient ; explicit)
            bcxu(i,j,k) = u(i,j,k)
            bcxv(i,j,k) = v(i,j,k)
            bcxw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PULSATILE_INFLOW)
            bcxu(i,j,k) = ux1*sin(twod*pi*freq_ux1*time)
            bcxv(i,j,k) = vx1*sin(twod*pi*freq_vx1*time)
            bcxw(i,j,k) = wx1*sin(twod*pi*freq_wx1*time)
          CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
            bcxu(i,j,k) = zero
            bcxv(i,j,k) = v(i,j,k)
            bcxw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PERIODIC)       ! periodic bc  (explicit & dirty implementation)
            bcxu(i,j,k) = half*( u(i,j,k) + u(nxc,j,k) )
            bcxv(i,j,k) = half*( v(i,j,k) + v(nxc,j,k) )
            bcxw(i,j,k) = half*( w(i,j,k) + w(nxc,j,k) )
          CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
            ux1 = u_blasius(j)
            bcxu(i,j,k) = ux1
            bcxv(i,j,k) = vx1
            bcxw(i,j,k) = wx1
          CASE (BC_TYPE_SHEAR)       ! shear bc
            bcxu(i,j,k) = ux1*y(j)/yout
            bcxv(i,j,k) = vx1
            bcxw(i,j,k) = wx1
        END SELECT

! right boundary
        i = nxc
        SELECT CASE (bcx2)
          CASE (BC_TYPE_DIRICHLET)      ! dirichlet bc
            bcxu(i,j,k) = ux2
            bcxv(i,j,k) = vx2
            bcxw(i,j,k) = wx2
          CASE (BC_TYPE_ZERO_GRADIENT)  ! outflow bc ( zero gradient ; explicit)
            bcxu(i,j,k) = u(i,j,k)
            bcxv(i,j,k) = v(i,j,k)
            bcxw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PULSATILE_INFLOW)

            bcxu(i,j,k) = ux2*sin(twod*pi*freq_ux2*time)

            if (time>tstart_gust) then
               bcxu(i,j,k) = ux2*abs( sin(2.0_CGREAL*pi*freq_ux2*time) )   !Changed by Wanh for gust purpose
            else
               bcxu(i,j,k) = 0
            endif

            bcxv(i,j,k) = vx2*sin(twod*pi*freq_vx2*time)
            bcxw(i,j,k) = wx2*sin(twod*pi*freq_wx2*time)
          CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
            bcxu(i,j,k) = zero
            bcxv(i,j,k) = v(i,j,k)
            bcxw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PERIODIC)       ! periodic bc
            bcxu(i,j,k) = half*( u(i,j,k) + u(1,j,k) )
            bcxv(i,j,k) = half*( v(i,j,k) + v(1,j,k) )
            bcxw(i,j,k) = half*( w(i,j,k) + w(1,j,k) )
          CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
            CALL blasius_velocity
            ux2 = u_blasius(j)
            bcxu(i,j,k) = ux2
            bcxv(i,j,k) = vx2
            bcxw(i,j,k) = wx2
          CASE (BC_TYPE_SHEAR)          ! shear bc
            bcxu(i,j,k) = ux2*y(j)/yout
            bcxv(i,j,k) = vx2
            bcxw(i,j,k) = wx2
        END SELECT

      ENDDO ! j
      ENDDO ! k

      DO k=0,nz
      DO i=0,nx

! bottom boundary
        j = 1
        SELECT CASE (bcy1)
          CASE (BC_TYPE_DIRICHLET)             ! dirichlet bc
            bcyu(i,j,k) = uy1
            bcyv(i,j,k) = vy1
            bcyw(i,j,k) = wy1
          CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
            bcyu(i,j,k) = u(i,j,k)
            bcyv(i,j,k) = v(i,j,k)
            bcyw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PULSATILE_INFLOW)
            bcyu(i,j,k) = uy1*sin(twod*pi*freq_uy1*time)
            bcyv(i,j,k) = vy1*sin(twod*pi*freq_vy1*time)
            bcyw(i,j,k) = wy1*sin(twod*pi*freq_wy1*time)
          CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
            bcyu(i,j,k) = u(i,j,k)
            bcyv(i,j,k) = zero
            bcyw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PERIODIC)       ! periodic bc
            bcyu(i,j,k) = half*( u(i,j,k) + u(i,nyc,k) )
            bcyv(i,j,k) = half*( v(i,j,k) + v(i,nyc,k) )
            bcyw(i,j,k) = half*( w(i,j,k) + w(i,nyc,k) )
          CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
            CALL blasius_velocity
            uy1 = u_blasius(j)
            bcyu(i,j,k) = uy1
            bcyv(i,j,k) = vy1
            bcyw(i,j,k) = wy1
          CASE (BC_TYPE_SHEAR)                ! shear bc
            bcyu(i,j,k) = uy1
            bcyv(i,j,k) = vy1*x(i)/xout
            bcyw(i,j,k) = wy1
        END SELECT

! top boundary
        j = nyc
        SELECT CASE (bcy2)
          CASE (BC_TYPE_DIRICHLET)             ! dirichlet bc
            bcyu(i,j,k) = uy2
            bcyv(i,j,k) = vy2
            bcyw(i,j,k) = wy2
          CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
            bcyu(i,j,k) = u(i,j,k)
            bcyv(i,j,k) = v(i,j,k)
            bcyw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PULSATILE_INFLOW)
            bcyu(i,j,k) = uy2*sin(twod*pi*freq_uy2*time)
            bcyv(i,j,k) = vy2*sin(twod*pi*freq_vy2*time)
            bcyw(i,j,k) = wy2*sin(twod*pi*freq_wy2*time)
          CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
            bcyu(i,j,k) = u(i,j,k)
            bcyv(i,j,k) = zero
            bcyw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PERIODIC)       ! periodic bc
            bcyu(i,j,k) = half*( u(i,j,k) + u(i,1,k) )
            bcyv(i,j,k) = half*( v(i,j,k) + v(i,1,k) )
            bcyw(i,j,k) = half*( w(i,j,k) + w(i,1,k) )
          CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
            CALL blasius_velocity
            uy2 = u_blasius(j)
            bcyu(i,j,k) = uy2
            bcyv(i,j,k) = vy2
            bcyw(i,j,k) = wy2
          CASE (BC_TYPE_SHEAR)             ! shear bc
            bcyu(i,j,k) = uy2
            bcyv(i,j,k) = vy2*x(i)/xout
            bcyw(i,j,k) = wy2
        END SELECT

      ENDDO ! i
      ENDDO ! k

      DO j=0,ny
      DO i=0,nx

! front boundary
        k = 1
        SELECT CASE (bcz1)
          CASE (BC_TYPE_DIRICHLET)             ! diriclet bc
            bczu(i,j,k) = uz1
            bczv(i,j,k) = vz1
            bczw(i,j,k) = wz1
          CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
            bczu(i,j,k) = u(i,j,k)
            bczv(i,j,k) = v(i,j,k)
            bczw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PULSATILE_INFLOW)
            bczu(i,j,k) = uz1*sin(twod*pi*freq_uz1*time)
            bczv(i,j,k) = vz1*sin(twod*pi*freq_vz1*time)
            bczw(i,j,k) = wz1*sin(twod*pi*freq_wz1*time)
          CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
            bczu(i,j,k) = u(i,j,k)
            bczv(i,j,k) = v(i,j,k)
            bczw(i,j,k) = 0.0
          CASE (BC_TYPE_PERIODIC)       ! periodic bc
            bczu(i,j,k) = half*( u(i,j,k) + u(i,j,nzc) )
            bczv(i,j,k) = half*( v(i,j,k) + v(i,j,nzc) )
            bczw(i,j,k) = half*( w(i,j,k) + w(i,j,nzc) )
          CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
            CALL blasius_velocity
            uz1 = u_blasius(j)
            bczu(i,j,k) = uz1
            bczv(i,j,k) = vz1
            bczw(i,j,k) = wz1
          CASE (BC_TYPE_SHEAR)          ! shear bc
            bczu(i,j,k) = uz1
            bczv(i,j,k) = vz1
            bczw(i,j,k) = wz1*y(j)/yout
        END SELECT

! back boundary
        k = nzc
        SELECT CASE (bcz2)
          CASE (BC_TYPE_DIRICHLET)             ! diriclet bc
            bczu(i,j,k) = uz2
            bczv(i,j,k) = vz2
            bczw(i,j,k) = wz2
          CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
            bczu(i,j,k) = u(i,j,k)
            bczv(i,j,k) = v(i,j,k)
            bczw(i,j,k) = w(i,j,k)
          CASE (BC_TYPE_PULSATILE_INFLOW)
            bczu(i,j,k) = uz2*sin(twod*pi*freq_uz2*time)
            bczv(i,j,k) = vz2*sin(twod*pi*freq_vz2*time)
            bczw(i,j,k) = wz2*sin(twod*pi*freq_wz2*time)
          CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
            bczu(i,j,k) = u(i,j,k)
            bczv(i,j,k) = v(i,j,k)
            bczw(i,j,k) = zero
          CASE (BC_TYPE_PERIODIC)       ! periodic bc
            bczu(i,j,k) = half*( u(i,j,k) + u(i,j,1) )
            bczv(i,j,k) = half*( v(i,j,k) + v(i,j,1) )
            bczw(i,j,k) = half*( w(i,j,k) + w(i,j,1) )
          CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
            CALL blasius_velocity
            uz2 = u_blasius(j)
            bczu(i,j,k) = uz2
            bczv(i,j,k) = vz2
            bczw(i,j,k) = wz2
          CASE (BC_TYPE_SHEAR)          ! shear bc
            bczu(i,j,k) = uz2
            bczv(i,j,k) = vz2
            bczw(i,j,k) = wz2*y(j)/yout
        END SELECT

      ENDDO ! i
      ENDDO ! j

!        CALL write_dump_debug('bcxu',1,bcxu)
!        CALL write_dump_debug('bcyv',1,bcyv)
!        CALL write_dump_debug('bczw',1,bczw)

! adjust BC to satisfy global mas conservation
!    If(MOD(ntime,5)==0) Then
      CALL enforce_global_mass_consv()
!    END IF

! set values for outer ghost points
      CALL set_outer_ghost_vel()

   END SUBROUTINE set_bc
!-------------------------------------------------------------------------

! compute mass flux at all boundaries and adjust outflow BC so as to satisfy
! global mass conservation.

   SUBROUTINE enforce_global_mass_consv()

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

       DO k=1,nzc
       DO j=1,nyc
       DO i=1,nxc
         massflux = massflux +                            &
                  ( -bcxu(i,j,k)*ium(i,j,k)*(1-iMarkm(i,j,k))*dy(j)*dz(k)   &
                    +bcxu(i,j,k)*iup(i,j,k)*(1-iMarkp(i,j,k))*dy(j)*dz(k)   &
                    -bcyv(i,j,k)*jum(i,j,k)*(1-jMarkm(i,j,k))*dx(i)*dz(k)   &
                    +bcyv(i,j,k)*jup(i,j,k)*(1-jMarkp(i,j,k))*dx(i)*dz(k)   &
                    -bczw(i,j,k)*kum(i,j,k)*(1-kMarkm(i,j,k))*dx(i)*dy(j)   &
                    +bczw(i,j,k)*kup(i,j,k)*(1-kMarkp(i,j,k))*dx(i)*dy(j) ) &
                 *REAL(1-iblank(i,j,k),KIND=CGREAL)
       ENDDO ! i
       ENDDO ! j
       ENDDO ! k

        correction_vel =-massflux/outflow_area

       IF ( MOD(ntime,nmonitor) == 0 ) THEN
         print *, '*****************'
         PRINT*,' SET_BC:massflux       = ',massflux
         PRINT*,' SET_BC:outflow_area   = ',outflow_area
         PRINT*,' SET_BC:correction_vel = ',correction_vel
         print *, '*****************'
         write(5669,*) nTime,massflux,correction_vel
       END IF ! ntime

       DO k=0,nz
       DO j=0,ny
         IF (bcx1 == BC_TYPE_ZERO_GRADIENT) bcxu(1,j,k)   = bcxu(1,j,k)   - correction_vel
         IF (bcx2 == BC_TYPE_ZERO_GRADIENT) bcxu(nxc,j,k) = bcxu(nxc,j,k) + correction_vel
       ENDDO ! j
       ENDDO ! k

       DO k=0,nz
       DO i=0,nx
         IF (bcy1 == BC_TYPE_ZERO_GRADIENT) bcyv(i,1,k)   = bcyv(i,1,k)   - correction_vel
         IF (bcy2 == BC_TYPE_ZERO_GRADIENT) bcyv(i,nyc,k) = bcyv(i,nyc,k) + correction_vel
       ENDDO ! i
       ENDDO ! k

       DO j=0,ny
       DO i=0,nx
         IF (bcz1 == BC_TYPE_ZERO_GRADIENT) bczw(i,j,1)   = bczw(i,j,1)   - correction_vel
         IF (bcz2 == BC_TYPE_ZERO_GRADIENT) bczw(i,j,nzc) = bczw(i,j,nzc) + correction_vel
       ENDDO ! i
       ENDDO ! j

    ENDIF ! bcx1

   END SUBROUTINE enforce_global_mass_consv
!-------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Purpose: Enforcing pressure boundary conditions for LSOR and MG method.
!
! Input: pressure values
!
! Output: enforced pressure values.
! Written by H. Dong
!------------------------------------------------------------------------------

   SUBROUTINE enforce_p_periodic(pres)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE pressure_arrays
    USE grid_arrays
    USE GCM_arrays

    IMPLICIT NONE

!    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT(INOUT) :: pres
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1) :: pres

    INTEGER             :: i,j,k,n

    IF (bcx1  .EQ.  BC_TYPE_PERIODIC .AND. &
        bcx2  .EQ.  BC_TYPE_PERIODIC) THEN

          DO k=0,nz
             DO j=0,ny
                 pres(0,j,k)  = pres(nxc,j,k)
                 pres(nx,j,k) = pres(1,j,k)
             ENDDO
          ENDDO

     END IF

     IF (bcy1 .EQ. BC_TYPE_PERIODIC .AND. &
         bcy2 .EQ. BC_TYPE_PERIODIC) THEN

         DO K = 0, nz
           DO I = 0, nx
                 pres(i,0,k)  = pres(i,nyc,k)
                 pres(i,ny,k) = pres(i,1,k)
            ENDDO
         ENDDO

     END IF

     IF (ndim .EQ. DIM_3D) THEN
         IF (bcz1 .EQ. BC_TYPE_PERIODIC .AND. &
             bcz2 .EQ. BC_TYPE_PERIODIC) THEN
          DO J = 0, ny
            DO I = 0, nx
                 pres(i,j,0)  = pres(i,j,nzc)
                 pres(i,j,nz) = pres(i,j,1)
            END DO
          END DO
         END IF
     END IF

   END SUBROUTINE enforce_p_periodic
!-------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Purpose: Enforcing velocity boundary conditions.
!
! Input: No
!
! Output: No
! Written by H. Dong
!------------------------------------------------------------------------------

   SUBROUTINE enforce_u_periodic
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER             :: i,j,k,n

    IF ( bcx1 == BC_TYPE_PERIODIC  .OR.  &
         bcx2 == BC_TYPE_PERIODIC  .OR.  &
         bcy1 == BC_TYPE_PERIODIC  .OR.  &
         bcy2 == BC_TYPE_PERIODIC  .OR.  &
         bcz1 == BC_TYPE_PERIODIC  .OR.  &
         bcz2 == BC_TYPE_PERIODIC       ) THEN

! remove outer boundary condition correction before modifying BC
       CALL remove_rhs_adjust_bc()

      IF (bcx1 == BC_TYPE_PERIODIC .AND. &
          bcx2 == BC_TYPE_PERIODIC) THEN
         DO k=0,nz
         DO j=0,ny
               bcxu(1,j,k) = half*( u(1,j,k) + u(nxc,j,k) )
               bcxv(1,j,k) = half*( v(1,j,k) + v(nxc,j,k) )
               bcxw(1,j,k) = half*( w(1,j,k) + w(nxc,j,k) )

               bcxu(nxc,j,k) =  bcxu(1,j,k)
               bcxv(nxc,j,k) =  bcxv(1,j,k)
               bcxw(nxc,j,k) =  bcxw(1,j,k)
         ENDDO ! j
         ENDDO ! k
      END IF

      IF (bcy1 == BC_TYPE_PERIODIC .AND. &
          bcy2 == BC_TYPE_PERIODIC) THEN
         DO k=0,nz
         DO i=0,nx
               bcyu(i,1,k) = half*( u(i,1,k) + u(i,nyc,k) )
               bcyv(i,1,k) = half*( v(i,1,k) + v(i,nyc,k) )
               bcyw(i,1,k) = half*( w(i,1,k) + w(i,nyc,k) )
               bcyu(i,nyc,k) = bcyu(i,1,k)
               bcyv(i,nyc,k) = bcyv(i,1,k)
               bcyw(i,nyc,k) = bcyw(i,1,k)
         ENDDO ! i
         ENDDO ! k
      END IF

     IF (bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC) THEN
         DO j=0,ny
         DO i=0,nx
              bczu(i,j,1) = half*( u(i,j,1) + u(i,j,nzc) )
              bczv(i,j,1) = half*( v(i,j,1) + v(i,j,nzc) )
              bczw(i,j,1) = half*( w(i,j,1) + w(i,j,nzc) )

              bczu(i,j,nzc) = bczu(i,j,1)
              bczv(i,j,nzc) = bczv(i,j,1)
              bczw(i,j,nzc) = bczw(i,j,1)
         ENDDO ! i
         ENDDO ! j
      END IF

! readjust outer boundary condition correction after modifying BC
       CALL rhs_readjust_bc()

    ENDIF ! bcx1

   END SUBROUTINE enforce_u_periodic

!------------------------------------------------------------------------------
!
! Purpose: Removing IUP, IUM, JUP, JUM, KUP, KUM for pressure periodic conditions.
!
! Input: No.
!
! Output: No.
! Written by H. Dong
!------------------------------------------------------------------------------

   SUBROUTINE remove_up_um

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,m

      IF (bcx1 == BC_TYPE_PERIODIC .AND. &
          bcx2 == BC_TYPE_PERIODIC) THEN
        DO j=0,ny
        DO k=0,nz
           ium(1,j,k)   = 0
           iup(nxc,j,k) = 0
        ENDDO
        ENDDO
      END IF

      IF (bcy1 == BC_TYPE_PERIODIC .AND. &
          bcy2 == BC_TYPE_PERIODIC) THEN
          DO i=0,nx
          DO k=0,nz
             jum(i,1,k)   = 0
             jup(i,nyc,k) = 0
          ENDDO
          ENDDO
      END IF

     IF (bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC) THEN
         DO i=0,nx
         DO j=0,ny
           kum(i,j,1)   = 0
           kup(i,j,nzc) = 0
         ENDDO
         ENDDO
     END IF

   END SUBROUTINE remove_up_um


!------------------------------------------------------------------------------
!
! Purpose: Adding IUP, IUM, JUP, JUM, KUP, KUM back for pressure periodic conditions.
!
! Input: No.
!
! Output: No.
! Written by H. Dong
!------------------------------------------------------------------------------

   SUBROUTINE add_up_um

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,m

      IF (bcx1 == BC_TYPE_PERIODIC .AND. &
          bcx2 == BC_TYPE_PERIODIC) THEN
        DO j=0,ny
        DO k=0,nz
           ium(1,j,k)   = 1
           iup(nxc,j,k) = 1
        ENDDO
        ENDDO
      END IF

      IF (bcy1 == BC_TYPE_PERIODIC .AND. &
          bcy2 == BC_TYPE_PERIODIC) THEN
          DO i=0,nx
          DO k=0,nz
             jum(i,1,k)   = 1
             jup(i,nyc,k) = 1
          ENDDO
          ENDDO
      END IF

     IF (bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC) THEN
         DO i=0,nx
         DO j=0,ny
           kum(i,j,1)   = 1
           kup(i,j,nzc) = 1
         ENDDO
         ENDDO
     END IF

   END SUBROUTINE add_up_um

   SUBROUTINE set_pressure_dirichlet_bc
!-----------------------------------------------|
! Purpose: Setting IUP, IUM, JUP, JUM, KUP, KUM |
! to be -1 at the outer boundary for Dirichlet|
! pressure condition.                           |
!                                               |
! Written by Z. Liang                           |
!-----------------------------------------------|

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,m

      IF (pbcx1 == PBC_DIRICHLET) THEN
        DO k=0,nz
        DO j=0,ny
           ium(1,  j,k) = -1
        ENDDO
        ENDDO
      END IF

      IF (pbcx2 == PBC_DIRICHLET) THEN
        DO k=0,nz
        DO j=0,ny
           iup(nxc,j,k) = -1
        ENDDO
        ENDDO
      END IF

      IF (pbcy1 == PBC_DIRICHLET) THEN
          DO k=0,nz
          DO i=0,nx
             jum(i,  1,k) = -1
          ENDDO
          ENDDO
      END IF

      IF (pbcy2 == PBC_DIRICHLET) THEN
          DO k=0,nz
          DO i=0,nx
             jup(i,nyc,k) = -1
          ENDDO
          ENDDO
      END IF

     IF (pbcz1 == PBC_DIRICHLET) THEN
         DO j=0,ny
         DO i=0,nx
           kum(i,j,  1) = -1
         ENDDO
         ENDDO
     END IF

     IF (pbcz2 == PBC_DIRICHLET) THEN
         DO j=0,ny
         DO i=0,nx
           kup(i,j,nzc) = -1
         ENDDO
         ENDDO
     END IF

   END SUBROUTINE set_pressure_dirichlet_bc
!------------------------------------------------------------------------------
   SUBROUTINE unset_pressure_dirichlet_bc
!-----------------------------------------------|
! Purpose: Restore IUP, IUM, JUP, JUM, KUP, KUM |
! to be 1 at the out boundary                   |
!                                               |
! Written by Z. Liang                           |
!-----------------------------------------------|

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,m

      IF (pbcx1 == PBC_DIRICHLET) THEN
        DO k=0,nz
        DO j=0,ny
           ium(1,  j,k) = 1
        ENDDO
        ENDDO
      END IF
      IF (pbcx2 == PBC_DIRICHLET) THEN
        DO k=0,nz
        DO j=0,ny
           iup(nxc,j,k) = 1
        ENDDO
        ENDDO
      END IF

      IF (pbcy1 == PBC_DIRICHLET) THEN
          DO k=0,nz
          DO i=0,nx
             jum(i,  1,k) = 1
          ENDDO
          ENDDO
      END IF
      IF (pbcy2 == PBC_DIRICHLET) THEN
          DO k=0,nz
          DO i=0,nx
             jup(i,nyc,k) = 1
          ENDDO
          ENDDO
      END IF

     IF (pbcz1 == PBC_DIRICHLET) THEN
         DO j=0,ny
         DO i=0,nx
           kum(i,j,  1) = 1
         ENDDO
         ENDDO
     END IF
     IF (pbcz2 == PBC_DIRICHLET) THEN
         DO j=0,ny
         DO i=0,nx
           kup(i,j,nzc) = 1
         ENDDO
         ENDDO
     END IF

   END SUBROUTINE unset_pressure_dirichlet_bc

