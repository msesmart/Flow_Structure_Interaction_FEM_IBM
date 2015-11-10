!-------------------------------------------------
!  SUBROUTINE solve_ad()
!  SUBROUTINE itsolv_ad(var,r)
!  SUBROUTINE itsolv_ad_x(var,r)
!  SUBROUTINE itsolv_ad_y(var,r)
!  SUBROUTINE itsolv_ad_z(var,r)
!  SUBROUTINE calc_residual_ad(var,r,resm)
!-------------------------------------------------

!-------------------------------------------------
!-------------------------------------------------
   SUBROUTINE solve_ad()

    USE global_parameters
    USE flow_parameters
    USE flow_parameters
    USE flow_arrays
    USE multiuse_arrays

    IMPLICIT NONE

!... Local Variables

    INTEGER           :: iter,i,j,k,loc(3),loc1(3),loc2(3),loc3(3)
    REAL(KIND=CGREAL) :: maxresu,maxresv,maxresw,maxres
!******************************************************************************

    iter      = 0
    maxres    = 1.E10_CGREAL 

! DEBUG
!    print *, 'before calc max(u), max,min(nlu)=', maxval(u), maxval(nlu),minval(nlu)
!    print *, 'before calc max(v), max,min(nlv)=', maxval(v), maxval(nlv),minval(nlv)
!    print *, 'before calc max,min(bcxu)=', maxval(bcxu), minval(bcxu)
!    print *, 'before calc max,min(bcyv)=', maxval(bcyv), minval(bcyv)
! END DEBUG

    maxresw=zero
    CALL calc_residual_ad(u,nlu,maxresu,loc1)
    CALL calc_residual_ad(v,nlv,maxresv,loc2)
    IF (ndim == DIM_3D) CALL calc_residual_ad(w,nlw,maxresw,loc3)
	  maxres = maxresu
	  IF (ABS(maxresv) > ABS(maxres)) maxres = maxresv
	  IF (ABS(maxresw) > ABS(maxres)) maxres = maxresw

! DEBUG
!    print *, 'after calc max(u), min(u)=', maxval(u), minval(u)
!    print *, 'after calc max(v), min(v)=', maxval(v), minval(v)
!    print *, 'after calc max,min(nlu)=', maxval(nlu),minval(nlu)
!    print *, 'after calc max,min(nlv)=', maxval(nlv),minval(nlv)
!    print *, 'after calc maxresu, maxresv=', maxresu,maxresv
!    print *, ' '
! END DEBUG

    SELECT CASE (boundary_formulation)
    CASE (NO_INTERNAL_BOUNDARY:SSM_METHOD)

!------This section was commented out and rewritten by H. Luo---------
!
!!$! solve u*
!!$          iter      = 0
!!$          maxres    = 1.E10_CGREAL 
!!$          DO WHILE ((iter .LT. iterMax_ad) .AND. (maxres .GT. restol_ad))
!!$
!!$	     CALL itsolv_ad(u,nlu)
!!$
!!$! new changes
!!$             CALL enforce_u_periodic
!!$	     CALL calc_residual_ad(u,nlu,maxres)
!!$	     iter=iter+1
!!$          ENDDO
!!$          IF ( iter .EQ. iterMax_ad .AND. maxres .GT. restol_ad ) THEN
!!$            PRINT*,'Velocity did not converge in ',iterMax_ad,' iterations'
!!$            PRINT*,'Final residual = ',maxres
!!$          ELSE
!!$            IF (MOD(ntime,nmonitor) == 0) THEN
!!$              PRINT*,'velocity convergence : ',iter,maxres
!!$            END IF ! ntime
!!$          ENDIF
!!$
!!$! solve v*
!!$	  iter=0
!!$          maxres    = 1.E10_CGREAL 
!!$          DO WHILE ((iter .LT. iterMax_ad) .AND. (maxres .GT. restol_ad))
!!$
!!$	     CALL itsolv_ad(v,nlv)
!!$
!!$! new changes
!!$             CALL enforce_u_periodic
!!$
!!$	     CALL calc_residual_ad(v,nlv,maxres)
!!$	     iter=iter+1
!!$          ENDDO
!!$          IF ( iter .EQ. iterMax_ad .AND. maxres .GT. restol_ad ) THEN
!!$            PRINT*,'Velocity did not converge in ',iterMax_ad,' iterations'
!!$            PRINT*,'Final residual = ',maxres
!!$          ELSE
!!$            IF (MOD(ntime,nmonitor) == 0) THEN
!!$              PRINT*,'velocity convergence : ',iter,maxres
!!$            END IF ! ntime
!!$          ENDIF
!!$
!!$! solve w*
!!$	  IF (ndim == DIM_3D) THEN
!!$	    iter=0
!!$            maxres    = 1.E10_CGREAL 
!!$            DO WHILE ((iter .LT. iterMax_ad) .AND. (maxres .GT. restol_ad))
!!$
!!$	       CALL itsolv_ad(w,nlw)
!!$
!!$! new changes
!!$            CALL enforce_u_periodic
!!$
!!$	       CALL calc_residual_ad(w,nlw,maxres)
!!$	       iter=iter+1
!!$            ENDDO
!!$            IF ( iter .EQ. iterMax_ad .AND. maxres .GT. restol_ad ) THEN
!!$              PRINT*,'Velocity did not converge in ',iterMax_ad,' iterations'
!!$              PRINT*,'Final residual = ',maxres
!!$            ELSE
!!$              IF (MOD(ntime,nmonitor) == 0) THEN
!!$                PRINT*,'velocity convergence : ',iter,maxres
!!$              END IF ! ntime
!!$            ENDIF
!!$	  ENDIF
      SELECT CASE (ad_solver_type)
      CASE (AD_SOLVER_TYPE_LSOR)
      CASE (AD_SOLVER_TYPE_MSIP)
        IF (ndim == DIM_2D) THEN
          CALL COEFF_AD_MSIP_2D(.TRUE.)
          CALL COEFF_LU_MSIP_2D()
        ELSE
          CALL COEFF_AD_MSIP_3D(.TRUE.)
          CALL COEFF_LU_MSIP_3D()
        END IF
      END SELECT

      DO WHILE ((iter .LT. iterMax_ad) .AND. (maxres .GT. restol_ad))

        SELECT CASE (ad_solver_type)
        CASE (AD_SOLVER_TYPE_LSOR)
          CALL itsolv_ad(u,nlu)
          CALL itsolv_ad(v,nlv)

          IF (ndim == DIM_3D) CALL itsolv_ad(w,nlw)
        CASE (AD_SOLVER_TYPE_MSIP)
          CALL itsolv_ad_MSIP(u,nlu)
          CALL itsolv_ad_MSIP(v,nlv)

          IF (ndim == DIM_3D) CALL itsolv_ad_MSIP(w,nlw)
        END SELECT
              
        CALL enforce_u_periodic

        IF(advec_scheme == CRANK_NICOLSON2) THEN
              
				  CALL face_vel()

          SELECT CASE (ad_solver_type)
          CASE (AD_SOLVER_TYPE_LSOR)
          CASE (AD_SOLVER_TYPE_MSIP)
            IF (ndim == DIM_2D) THEN
              CALL COEFF_AD_MSIP_2D(.FALSE.)
              CALL COEFF_LU_MSIP_2D()
            ELSE
              CALL COEFF_AD_MSIP_3D(.FALSE.)
              CALL COEFF_LU_MSIP_3D()
            END IF
          END SELECT
                
        ENDIF
              
        maxresw=zero
        SELECT CASE (ad_solver_type)
        CASE (AD_SOLVER_TYPE_LSOR)
          CALL calc_residual_ad(u,nlu,maxresu,loc1)
          CALL calc_residual_ad(v,nlv,maxresv,loc2)
          IF (ndim == DIM_3D) &
          CALL calc_residual_ad(w,nlw,maxresw,loc3)
        CASE (AD_SOLVER_TYPE_MSIP)
          CALL calc_residual_MSIP(u,nlu,maxresu,loc1)
          CALL calc_residual_MSIP(v,nlv,maxresv,loc2)
				  IF (ndim == DIM_3D) &
				  CALL calc_residual_MSIP(w,nlw,maxresw,loc3)
        END SELECT
                  
        maxres = maxresu
			  loc=loc1
			  IF (ABS(maxresv) > ABS(maxres)) THEN
				  maxres = maxresv
				  loc=loc2
			  ENDIF
        IF (ABS(maxresw) > ABS(maxres)) THEN
          maxres = maxresw
				  loc=loc3
        ENDIF
        iter = iter + 1

        IF ( iter .EQ. iterMax_ad .AND. maxres .GT. restol_ad ) THEN
            PRINT*,'Velocity did not converge in ',iterMax_ad,' iterations'
            PRINT*,'Final residual = ',maxres
        ELSE
            IF (MOD(ntime,nmonitor) == 0) THEN
!                    PRINT*,'velocity convergence SSM: ',iter,maxres,maxresu,maxresv
!                    WRITE(STDOUT,'(1X,A,3X,I4,4(1X,1PE19.11))') 'velocity convergence : ',iter,maxres,maxresu,maxresv,maxresw
              WRITE(STDOUT,'(1X,A,3X,I4,3(1X,1PE17.10),3(2X,I4))') 'velocity convergence : ',iter,maxresu,maxresv,maxresw,loc
            END IF ! ntime
        ENDIF

      ENDDO

      SELECT CASE (ad_solver_type)
      CASE (AD_SOLVER_TYPE_LSOR)
      CASE (AD_SOLVER_TYPE_MSIP)
        DEALLOCATE(LU,CA)
      END SELECT

	  CASE (GCM_METHOD)
      maxres    = 1.E10_CGREAL 
      DO WHILE ((iter .LT. iterMax_ad) .AND. (maxres .GT. restol_ad))

        CALL itsolv_ad(u,nlu)
        CALL itsolv_ad(v,nlv)

        IF (ndim == DIM_3D) CALL itsolv_ad(w,nlw)

        CALL GCM_GhostCell_Vel()

        ! Note that computation of face_v is redundant
        ! when the zero gradient condition is specified for the velocity
        ! because face_v is updated in GCM_enforce_global_mass_consv(),
        ! in GCM_GhostCell_Vel(). If that is the case, the time scheme
        ! is CRANK_NICOLSON2 even if advec_scheme is set to be 
        ! CRANK_NICOLSON1.
        ! -- H. Luo
        IF(advec_scheme == CRANK_NICOLSON2) CALL face_vel()

        CALL calc_residual_ad(u,nlu,maxresu,loc1)
        CALL calc_residual_ad(v,nlv,maxresv,loc2)

        maxres = maxresu
        IF (ABS(maxresv) > ABS(maxres)) maxres = maxresv
        IF (ndim == DIM_3D) THEN
          CALL calc_residual_ad(w,nlw,maxresw,loc3)
	        IF (ABS(maxresw) > ABS(maxres)) maxres = maxresw
        ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	    write(890,*)iter,maxresu,maxresv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        iter = iter + 1

        IF ( iter .EQ. iterMax_ad .AND. maxres .GT. restol_ad ) THEN
            PRINT*,'Velocity did not converge in ',iterMax_ad,' iterations'
            PRINT*,'Final residual = ',maxres
        ELSE
            IF (MOD(ntime,nmonitor) == 0) THEN
              PRINT*,'velocity convergence GCM: ',iter,maxres
            END IF ! ntime
        ENDIF
      ENDDO ! do while

    END SELECT



  END SUBROUTINE solve_ad

!----------------------------------------------------------
! currently coded as Line SOR with Gauss Siedel as smoother
!----------------------------------------------------------
   SUBROUTINE itsolv_ad(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

    
    IF (advec_scheme == ADAMS_BASHFORTH2) THEN
! Line solve in the x-direction
    CALL itsolv_ad_x_AB(var,r)

! Line solve in the y-direction    
    CALL itsolv_ad_y_AB(var,r)

! Line solve in the z-direction    
    IF (nDim == DIM_3D) CALL itsolv_ad_z_AB(var,r)

    ELSE
! Line solve in the x-direction
    CALL itsolv_ad_x(var,r)

! Line solve in the y-direction    
    CALL itsolv_ad_y(var,r)

! Line solve in the z-direction    
    IF (nDim == DIM_3D) CALL itsolv_ad_z(var,r)
    
    ENDIF
    
   END SUBROUTINE itsolv_ad
!----------------------------------------------

   SUBROUTINE itsolv_ad_x(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

!... Loop Variables

    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: omega_ad, half_dt
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2

    INTEGER :: iDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

!******************************************************************************

    iDirection = 1
    omega_ad = oned

    ! H. Luo start-------------------
    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = half * dt
    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
       half_dt  = zero
    ENDIF
    ! H. Luo end---------------------

    CALL enforce_u_periodic

! this loop is to transfer to RHS all the weights that do not fit on the 
! x-tdma

! Line solve in the x-direction
    
    DO k=1,nzc
    DO j=1,nyc

       JP   = j+1 
       JM   = j-1 
 
       KP   = k+1 
       KM   = k-1
 
      DO i=1,nxc

       riblank   = oned-(REAL(iblank(i,j,k),KIND=CGREAL))

       amx(i) = amx_ad(i,j,k)
       apx(i) = apx_ad(i,j,k)
       acx(i) =- ( amx(i) + apx(i) )      

       amy(j) = amy_ad(i,j,k)
       apy(j) = apy_ad(i,j,k)
       acy(j) =- ( amy(j) + apy(j) )     

       amz(k) = amz_ad(i,j,k)
       apz(k) = apz_ad(i,j,k)
       acz(k) =- ( amz(k) + apz(k) )    

!--------Added by H. Luo------------------------------------------------
! Note that this part only takes effect when the adection terms were
! treated with the CN scheme, since otherwise half_dt = 0.

       tmp1   = oned - fx(i  );   tmp2 = fx(i+1) 
       amx(i) = amx(i) - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
       apx(i) = apx(i) + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i)

       tmp1   =         fx(i  ) *(1 - ium(i,j,k))
       tmp2   = (oned - fx(i+1))*(1 - iup(i,j,k))
       acx(i) = acx(i) + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i  ,j,k)*tmp1)*dxinv(i)

       tmp1   = oned - fy(j  );   tmp2   =        fy(j+1)
       amy(j) = amy(j) - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
       apy(j) = apy(j) + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)

       tmp1   =         fy(j  ) *(1 - jum(i,j,k))
       tmp2   = (oned - fy(j+1))*(1 - jup(i,j,k))
       acy(j) = acy(j) + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j  ,k)*tmp1)*dyinv(j)

       tmp1   = oned - fz(k  );   tmp2   =        fz(k+1)
       amz(k) = amz(k) - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
       apz(k) = apz(k) + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)

       tmp1   =         fz(k  ) *(1 - kum(i,j,k))
       tmp2   = (oned - fz(k+1))*(1 - kup(i,j,k))
       acz(k) = acz(k) + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k  )*tmp1)*dzinv(k) 
!--------End adding-------------------------------------------------------

       rhs(i) = r(i,j,k) - ( var(i,JM,k)*amy(j)*(1-jum(i,j,k))  &
                            +var(i,JP,k)*apy(j)*(1-jup(i,j,k))  &
                            +var(i,j,KM)*amz(k)*(1-kum(i,j,k))  &
                            +var(i,j,KP)*apz(k)*(1-kup(i,j,k)) )

       amx(i) = amx(i)*riblank*(1-ium(i,j,k))
       apx(i) = apx(i)*riblank*(1-iup(i,j,k))
       acx(i) = oned + ( acx(i)+acy(j)+acz(k) ) *riblank   &
                            *(oned-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))
       rhs(i) = rhs(i)*riblank + var(i,j,k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL) 
              
       IF ( boundary_formulation == GCM_METHOD .AND. &
            fresh_cell(i,j,k) == 1                   ) THEN
         iFr = i; jFr = j; kFr = k;
         CALL GCM_correct_rhs_ad(iDirection,iFr,jFr,kFr,nx,rhs,var)
       END iF
       
      ENDDO ! i

!      IF (bcx1 == BC_TYPE_PERIODIC .AND. &
!          bcx2 == BC_TYPE_PERIODIC) THEN
!          beta_p = amx(1)
!          alpha_p = apx(nxc)

!          CALL cyclic_tdma(amx,acx,apx,rhs,dummy,alpha_p,beta_p,1,nxc)
!      ELSE

        IF (bcx1 == BC_TYPE_PERIODIC .AND. &
            bcx2 == BC_TYPE_PERIODIC) THEN
           rhs(1) = rhs(1) - var(nxc,j,k)*amx(1)
           rhs(nxc) = rhs(nxc) - var(1,j,k)*apx(nxc)
        END IF

         CALL tdma(amx,acx,apx,rhs,dummy,1,nxc)
!      END IF

      DO i=1,nxc
!      var(i,j,k) = var(i,j,k) + omega_ad*(dummy(i)-var(i,j,k))
       var(i,j,k) = dummy(i)
      ENDDO

    ENDDO ! j
    ENDDO ! k

   END SUBROUTINE itsolv_ad_x
!----------------------------------------------
!
!----------------------------------------------
   SUBROUTINE itsolv_ad_y(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

!... Loop Variables

    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: omega_ad, half_dt
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2

    INTEGER :: jDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

!******************************************************************************

    jDirection = 2
    omega_ad = oned

    ! H. Luo start-------------------
    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = half * dt
    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
       half_dt  = zero
    ENDIF
    ! H. Luo end---------------------

    CALL enforce_u_periodic 

! Line solve in the y-direction

    DO k=1,nzc
    DO i=1,nxc

       IP   = i + 1 
       IM   = i - 1 
 
       KP   = k + 1 
       KM   = k - 1 

      DO j=1,nyc

       riblank   = oned-(REAL(iblank(i,j,k),KIND=CGREAL))

       amx(i) = amx_ad(i,j,k)
       apx(i) = apx_ad(i,j,k)
       acx(i) =- ( amx(i) + apx(i) )                           

       amy(j) = amy_ad(i,j,k)
       apy(j) = apy_ad(i,j,k)
       acy(j) =- ( amy(j) + apy(j) )                          

       amz(k) = amz_ad(i,j,k)
       apz(k) = apz_ad(i,j,k)
       acz(k) =- ( amz(k) + apz(k) )   

!--------Added by H. Luo------------------------------------------------
! Note that this part only takes effect when the adection terms were
! treated with the CN scheme, since otherwise half_dt = 0.

       tmp1   = oned - fx(i  );   tmp2 = fx(i+1)
       amx(i) = amx(i) - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
       apx(i) = apx(i) + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i)

       tmp1   =         fx(i  ) *(1 - ium(i,j,k))
       tmp2   = (oned - fx(i+1))*(1 - iup(i,j,k))
       acx(i) = acx(i) + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i  ,j,k)*tmp1)*dxinv(i)

       tmp1   = oned - fy(j  );   tmp2   =        fy(j+1)
       amy(j) = amy(j) - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
       apy(j) = apy(j) + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)

       tmp1   =         fy(j  ) *(1 - jum(i,j,k))
       tmp2   = (oned - fy(j+1))*(1 - jup(i,j,k))
       acy(j) = acy(j) + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j  ,k)*tmp1)*dyinv(j)

       tmp1   = oned - fz(k  );   tmp2   =        fz(k+1)
       amz(k) = amz(k) - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
       apz(k) = apz(k) + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)

       tmp1   =         fz(k  ) *(1 - kum(i,j,k))
       tmp2   = (oned - fz(k+1))*(1 - kup(i,j,k))
       acz(k) = acz(k) + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k  )*tmp1)*dzinv(k)
!--------End adding-------------------------------------------------------
                      
       rhs(j) = r(i,j,k) -(  var(IM,j,k)*amx(i)*(1-ium(i,j,k))  &
                            +var(IP,j,k)*apx(i)*(1-iup(i,j,k))  &
                            +var(i,j,KM)*amz(k)*(1-kum(i,j,k))  &
                            +var(i,j,KP)*apz(k)*(1-kup(i,j,k)) )

       amy(j) = amy(j)*riblank*(1-jum(i,j,k))
       apy(j) = apy(j)*riblank*(1-jup(i,j,k))
       acy(j) = oned + ( acx(i)+acy(j)+acz(k) ) * riblank &
                            *(oned-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))
       rhs(j) = rhs(j)*riblank + var(i,j,k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL) 
       
       IF ( boundary_formulation == GCM_METHOD .AND. &
            fresh_cell(i,j,k) == 1                   ) THEN
         iFr = i; jFr = j; kFr = k;
         CALL GCM_correct_rhs_ad(jDirection,iFr,jFr,kFr,ny,rhs,var)
       ENDIF ! fresh_cell
       
      ENDDO ! j

!      IF (bcy1 == BC_TYPE_PERIODIC .AND. &
!          bcy2 == BC_TYPE_PERIODIC) THEN
!          beta_p = amy(1)
!          alpha_p = apy(nyc)

!          CALL cyclic_tdma(amy,acy,apy,rhs,dummy,alpha_p,beta_p,1,nyc)
!      ELSE

        IF (bcy1 == BC_TYPE_PERIODIC .AND. &
            bcy2 == BC_TYPE_PERIODIC) THEN
           rhs(1) = rhs(1) - var(i,nyc,k)*amy(1)
           rhs(nyc) = rhs(nyc) - var(i,1,k)*apy(nyc)
        END IF

        CALL tdma(amy,acy,apy,rhs,dummy,1,nyc)

!      END IF

      DO j=1,nyc
!      var(i,j,k) = var(i,j,k) + omega_ad*(dummy(j)-var(i,j,k))
       var(i,j,k) = dummy(j)
      ENDDO ! j

    ENDDO ! i
    ENDDO ! k

   END SUBROUTINE itsolv_ad_y
!----------------------------------------------------------

!----------------------------------------------------------
   SUBROUTINE itsolv_ad_z(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

!... Loop Variables
    
    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: omega_ad, half_dt
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2

    INTEGER :: kDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

!******************************************************************************

    kDirection = 3
    omega_ad = oned

    ! H. Luo start-------------------
    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = half * dt
    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
       half_dt  = zero
    ENDIF
    ! H. Luo end---------------------

    CALL enforce_u_periodic 
    
! Line solver in the z-direction
    DO i=1,nxc
    DO j=1,nyc
       IP   = i + 1 
       IM   = i - 1 
 
       JP   = j + 1 
       JM   = j - 1 

      DO k=1,nzc

       riblank   = oned-(REAL(iblank(i,j,k),KIND=CGREAL))

       amx(i) = amx_ad(i,j,k)
       apx(i) = apx_ad(i,j,k)
       acx(i) =- ( amx(i) + apx(i) )                           

       amy(j) = amy_ad(i,j,k)
       apy(j) = apy_ad(i,j,k)
       acy(j) =- ( amy(j) + apy(j) )                          

       amz(k) = amz_ad(i,j,k)
       apz(k) = apz_ad(i,j,k)
       acz(k) =- ( amz(k) + apz(k) )                         

!--------Added by H. Luo------------------------------------------------
! Note that this part only takes effect when the adection terms were
! treated with the CN scheme, since otherwise half_dt = 0.

       tmp1   = oned - fx(i  );   tmp2 = fx(i+1)
       amx(i) = amx(i) - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
       apx(i) = apx(i) + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i)

       tmp1   =         fx(i  ) *(1 - ium(i,j,k))
       tmp2   = (oned - fx(i+1))*(1 - iup(i,j,k))
       acx(i) = acx(i) + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i  ,j,k)*tmp1)*dxinv(i)

       tmp1   = oned - fy(j  );   tmp2   =        fy(j+1)
       amy(j) = amy(j) - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
       apy(j) = apy(j) + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)

       tmp1   =         fy(j  ) *(1 - jum(i,j,k))
       tmp2   = (oned - fy(j+1))*(1 - jup(i,j,k))
       acy(j) = acy(j) + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j  ,k)*tmp1)*dyinv(j)

       tmp1   = oned - fz(k  );   tmp2   =        fz(k+1)
       amz(k) = amz(k) - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
       apz(k) = apz(k) + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)

       tmp1   =         fz(k  ) *(1 - kum(i,j,k))
       tmp2   = (oned - fz(k+1))*(1 - kup(i,j,k))
       acz(k) = acz(k) + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k  )*tmp1)*dzinv(k)
!--------End adding-------------------------------------------------------

       rhs(k) = r(i,j,k) - var(i,JM,k)*amy(j)*(1-jum(i,j,k))  &
                         - var(i,JP,k)*apy(j)*(1-jup(i,j,k))  &
                         - var(IM,j,k)*amx(i)*(1-ium(i,j,k))  &
                         - var(IP,j,k)*apx(i)*(1-iup(i,j,k))

       amz(k) = amz(k)*riblank*(oned-kum(i,j,k))
       apz(k) = apz(k)*riblank*(oned-kup(i,j,k))
       acz(k) = oned + ( acx(i)+acy(j)+acz(k) ) * riblank &
                            *(oned-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL)) !  <--- transpose this properly RM 
       rhs(k) = rhs(k)*riblank + var(i,j,k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL)   !<--- transpose this properly RM
       
       
        IF ( boundary_formulation == GCM_METHOD .AND. &
             fresh_cell(i,j,k) == 1                   ) THEN
          iFr = i; jFr = j; kFr = k;
          CALL GCM_correct_rhs_ad(kDirection,iFr,jFr,kFr,ny,rhs,var)
        ENDIF ! fresh_cell
       
      ENDDO ! k

!     IF (bcz1 == BC_TYPE_PERIODIC .AND. &
!         bcz2 == BC_TYPE_PERIODIC) THEN
!        beta_p = amz(1)
!        alpha_p = apz(nzc)

!        CALL cyclic_tdma(amz,acz,apz,rhs,dummy,alpha_p,beta_p,1,nzc)
!     ELSE

        IF (bcz1 == BC_TYPE_PERIODIC .AND. &
            bcz2 == BC_TYPE_PERIODIC) THEN
           rhs(1) = rhs(1) - var(i,j,nzc)*amz(1)
           rhs(nzc) = rhs(nzc) - var(i,j,1)*apz(nzc)
        END IF

        CALL tdma(amz,acz,apz,rhs,dummy,1,nzc)
!     END IF

      DO k=1,nzc
!      var(i,j,k) = var(i,j,k) + omega_ad*(dummy(k)-var(i,j,k))
       var(i,j,k) = dummy(k)
      ENDDO !k

    ENDDO !j 
    ENDDO !i
    
   END SUBROUTINE itsolv_ad_z
!----------------------------------------------------------
! 
!----------------------------------------------------------
   SUBROUTINE calc_residual_ad(var,r,resm,loc)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays  ! H. Luo
 
    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)  ::var,r
    REAL(KIND=CGREAL),                                   INTENT (OUT) ::resm
    INTEGER,           DIMENSION(3),                    INTENT (OUT) ::loc

!... Loop Variables

    INTEGER            :: i,j,k

!... Local Variables
    
    INTEGER            :: iFr,jFr,kFr
    INTEGER            :: IP,IM,JP,JM,KP,KM

    REAL(KIND=CGREAL)  :: res,killForGCMFresh
    REAL(KIND=CGREAL)  :: half_dt, tmp1, tmp2

!******************************************************************************
    loc = 0

    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = half * dt
    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
       half_dt  = zero
    ENDIF

    CALL enforce_u_periodic

    resm = zero

    DO k=1,nzc
    DO j=1,nyc
    DO i=1,nxc

       IP   = i + 1 
       IM   = i - 1 
 
       JP   = j + 1 
       JM   = j - 1 
 
       KP   = k + 1 
       KM   = k - 1 

!--------Added by H. Luo------------------------------------------------

       amx(i) = amx_ad(i,j,k)
       apx(i) = apx_ad(i,j,k)
       acx(i) =- ( amx(i) + apx(i) )

       amy(j) = amy_ad(i,j,k)
       apy(j) = apy_ad(i,j,k)
       acy(j) =- ( amy(j) + apy(j) )

       amz(k) = amz_ad(i,j,k)
       apz(k) = apz_ad(i,j,k)
       acz(k) =- ( amz(k) + apz(k) )

! Note that this part only takes effect when the adection terms were
! treated with the CN scheme, since otherwise half_dt = 0.

       tmp1   = oned - fx(i  );   tmp2 = fx(i+1)
       amx(i) = amx(i) - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
       apx(i) = apx(i) + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i)

       tmp1   =         fx(i  ) *(1 - ium(i,j,k))
       tmp2   = (oned - fx(i+1))*(1 - iup(i,j,k))
       acx(i) = acx(i) + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i  ,j,k)*tmp1)*dxinv(i)

       tmp1   = oned - fy(j  );   tmp2   =        fy(j+1)
       amy(j) = amy(j) - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
       apy(j) = apy(j) + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)

       tmp1   =         fy(j  ) *(1 - jum(i,j,k))
       tmp2   = (oned - fy(j+1))*(1 - jup(i,j,k))
       acy(j) = acy(j) + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j  ,k)*tmp1)*dyinv(j)

       tmp1   = oned - fz(k  );   tmp2   =        fz(k+1)
       amz(k) = amz(k) - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
       apz(k) = apz(k) + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)

       tmp1   =         fz(k  ) *(1 - kum(i,j,k))
       tmp2   = (oned - fz(k+1))*(1 - kup(i,j,k))
       acz(k) = acz(k) + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k  )*tmp1)*dzinv(k)

!--------End adding-------------------------------------------------------

       killForGCMFresh = (oned-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))

       !res    = r(i,j,k) - var(i,j,k)*( oned                                          &
       !                                - ( (amx_ad(i,j,k)+apx_ad(i,j,k))                    &
       !                                   +(amy_ad(i,j,k)+apy_ad(i,j,k))                    &
       !                                   +(amz_ad(i,j,k)+apz_ad(i,j,k)) )*killForGCMFresh )&
       !                  - var(IM,j,k)*amx_ad(i,j,k)*(oned-ium(i,j,k))               &
       !                  - var(IP,j,k)*apx_ad(i,j,k)*(oned-iup(i,j,k))               &
       !                  - var(i,JM,k)*amy_ad(i,j,k)*(oned-jum(i,j,k))               &
       !                  - var(i,JP,k)*apy_ad(i,j,k)*(oned-jup(i,j,k))               &
       !                  - var(i,j,KM)*amz_ad(i,j,k)*(oned-kum(i,j,k))               &
       !                  - var(i,j,KP)*apz_ad(i,j,k)*(oned-kup(i,j,k))

       res    = r(i,j,k) - var(i,j,k)*( oned                                         &
                                       + ( acx(i) + acy(j) + acz(k) )*killForGCMFresh )    &
                         - var(IM,j,k)*amx(i)*(1-ium(i,j,k))               &
                         - var(IP,j,k)*apx(i)*(1-iup(i,j,k))               &
                         - var(i,JM,k)*amy(j)*(1-jum(i,j,k))               &
                         - var(i,JP,k)*apy(j)*(1-jup(i,j,k))               &
                         - var(i,j,KM)*amz(k)*(1-kum(i,j,k))               &
                         - var(i,j,KP)*apz(k)*(1-kup(i,j,k))

       IF ( boundary_formulation == GCM_METHOD .AND. &
            fresh_cell(i,j,k) == 1                   ) THEN
         iFr = i; jFr = j; kFr = k;
         CALL GCM_correct_res_ad(iFr,jFr,kFr,var,res)
       ENDIF ! fresh_cell
       
       res    = res*REAL(1-iblank(i,j,k),KIND=CGREAL)

       IF (ABS(res) > resm ) THEN
          resm = ABS(res)
         loc(1) = i
         loc(2) = j
         loc(3) = k
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !print*, 'Location of the max. residual:', iLoc, jLoc, kLoc

  END SUBROUTINE calc_residual_ad  
!----------------------------------------------------------

!----------------------------------------------------------
! MSIP as smoother
!----------------------------------------------------------
   SUBROUTINE itsolv_ad_MSIP(var,r)

    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

    IF (nDim == DIM_2D) THEN
       CALL itsolv_MSIP_2D(var,r,restol_ad)
    ELSE
       CALL itsolv_MSIP_3D(var,r,restol_ad)
    ENDIF
    
   END SUBROUTINE itsolv_ad_MSIP
!----------------------------------------------

   SUBROUTINE itsolv_ad_x_AB(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

!... Loop Variables

    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: omega_ad, half_dt
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2

    INTEGER :: iDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

!******************************************************************************

    iDirection = 1
    omega_ad = oned

    CALL enforce_u_periodic

! this loop is to transfer to RHS all the weights that do not fit on the 
! x-tdma

! Line solve in the x-direction
    
    DO k=1,nzc
    DO j=1,nyc

       JP   = j+1 
       JM   = j-1 
 
       KP   = k+1 
       KM   = k-1
 
      DO i=1,nxc

       riblank   = oned-(REAL(iblank(i,j,k),KIND=CGREAL))

       amx(i) = amx_ad(i,j,k)
       apx(i) = apx_ad(i,j,k)
       acx(i) =- ( amx(i) + apx(i) )      

       amy(j) = amy_ad(i,j,k)
       apy(j) = apy_ad(i,j,k)
       acy(j) =- ( amy(j) + apy(j) )     

       amz(k) = amz_ad(i,j,k)
       apz(k) = apz_ad(i,j,k)
       acz(k) =- ( amz(k) + apz(k) )    

       rhs(i) = r(i,j,k) - ( var(i,JM,k)*amy(j)*(1-jum(i,j,k))  &
                            +var(i,JP,k)*apy(j)*(1-jup(i,j,k))  &
                            +var(i,j,KM)*amz(k)*(1-kum(i,j,k))  &
                            +var(i,j,KP)*apz(k)*(1-kup(i,j,k)) )

       amx(i) = amx(i)*riblank*(1-ium(i,j,k))
       apx(i) = apx(i)*riblank*(1-iup(i,j,k))
       acx(i) = oned + ( acx(i)+acy(j)+acz(k) ) *riblank   &
                            *(oned-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))
       rhs(i) = rhs(i)*riblank + var(i,j,k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL) 
              
       IF ( boundary_formulation == GCM_METHOD .AND. &
            fresh_cell(i,j,k) == 1                   ) THEN
         iFr = i; jFr = j; kFr = k;
         CALL GCM_correct_rhs_ad(iDirection,iFr,jFr,kFr,nx,rhs,var)
       END iF
       
      ENDDO ! i

!      IF (bcx1 == BC_TYPE_PERIODIC .AND. &
!          bcx2 == BC_TYPE_PERIODIC) THEN
!          beta_p = amx(1)
!          alpha_p = apx(nxc)

!          CALL cyclic_tdma(amx,acx,apx,rhs,dummy,alpha_p,beta_p,1,nxc)
!      ELSE

        IF (bcx1 == BC_TYPE_PERIODIC .AND. &
            bcx2 == BC_TYPE_PERIODIC) THEN
           rhs(1) = rhs(1) - var(nxc,j,k)*amx(1)
           rhs(nxc) = rhs(nxc) - var(1,j,k)*apx(nxc)
        END IF

         CALL tdma(amx,acx,apx,rhs,dummy,1,nxc)
!      END IF

      DO i=1,nxc
!      var(i,j,k) = var(i,j,k) + omega_ad*(dummy(i)-var(i,j,k))
       var(i,j,k) = dummy(i)
      ENDDO

    ENDDO ! j
    ENDDO ! k

   END SUBROUTINE itsolv_ad_x_AB
!----------------------------------------------
!
!----------------------------------------------
   SUBROUTINE itsolv_ad_y_AB(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

!... Loop Variables

    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: omega_ad, half_dt
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2

    INTEGER :: jDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

!******************************************************************************

    jDirection = 2
    omega_ad = oned

    CALL enforce_u_periodic 

! Line solve in the y-direction

    DO k=1,nzc
    DO i=1,nxc

       IP   = i + 1 
       IM   = i - 1 
 
       KP   = k + 1 
       KM   = k - 1 

      DO j=1,nyc

       riblank   = oned-(REAL(iblank(i,j,k),KIND=CGREAL))

       amx(i) = amx_ad(i,j,k)
       apx(i) = apx_ad(i,j,k)
       acx(i) =- ( amx(i) + apx(i) )                           

       amy(j) = amy_ad(i,j,k)
       apy(j) = apy_ad(i,j,k)
       acy(j) =- ( amy(j) + apy(j) )                          

       amz(k) = amz_ad(i,j,k)
       apz(k) = apz_ad(i,j,k)
       acz(k) =- ( amz(k) + apz(k) )   

       rhs(j) = r(i,j,k) -(  var(IM,j,k)*amx(i)*(1-ium(i,j,k))  &
                            +var(IP,j,k)*apx(i)*(1-iup(i,j,k))  &
                            +var(i,j,KM)*amz(k)*(1-kum(i,j,k))  &
                            +var(i,j,KP)*apz(k)*(1-kup(i,j,k)) )

       amy(j) = amy(j)*riblank*(1-jum(i,j,k))
       apy(j) = apy(j)*riblank*(1-jup(i,j,k))
       acy(j) = oned + ( acx(i)+acy(j)+acz(k) ) * riblank &
                            *(oned-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))
       rhs(j) = rhs(j)*riblank + var(i,j,k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL) 
       
       IF ( boundary_formulation == GCM_METHOD .AND. &
            fresh_cell(i,j,k) == 1                   ) THEN
         iFr = i; jFr = j; kFr = k;
         CALL GCM_correct_rhs_ad(jDirection,iFr,jFr,kFr,ny,rhs,var)
       ENDIF ! fresh_cell
       
      ENDDO ! j

!      IF (bcy1 == BC_TYPE_PERIODIC .AND. &
!          bcy2 == BC_TYPE_PERIODIC) THEN
!          beta_p = amy(1)
!          alpha_p = apy(nyc)

!          CALL cyclic_tdma(amy,acy,apy,rhs,dummy,alpha_p,beta_p,1,nyc)
!      ELSE

        IF (bcy1 == BC_TYPE_PERIODIC .AND. &
            bcy2 == BC_TYPE_PERIODIC) THEN
           rhs(1) = rhs(1) - var(i,nyc,k)*amy(1)
           rhs(nyc) = rhs(nyc) - var(i,1,k)*apy(nyc)
        END IF

        CALL tdma(amy,acy,apy,rhs,dummy,1,nyc)

!      END IF

      DO j=1,nyc
!      var(i,j,k) = var(i,j,k) + omega_ad*(dummy(j)-var(i,j,k))
       var(i,j,k) = dummy(j)
      ENDDO ! j

    ENDDO ! i
    ENDDO ! k

   END SUBROUTINE itsolv_ad_y_AB
!----------------------------------------------------------

!----------------------------------------------------------
   SUBROUTINE itsolv_ad_z_AB(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

!... Loop Variables
    
    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: omega_ad, half_dt
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2

    INTEGER :: kDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

!******************************************************************************

    kDirection = 3
    omega_ad = oned

    CALL enforce_u_periodic 
    
! Line solver in the z-direction
    DO i=1,nxc
    DO j=1,nyc
       IP   = i + 1 
       IM   = i - 1 
 
       JP   = j + 1 
       JM   = j - 1 

      DO k=1,nzc

       riblank   = oned-(REAL(iblank(i,j,k),KIND=CGREAL))

       amx(i) = amx_ad(i,j,k)
       apx(i) = apx_ad(i,j,k)
       acx(i) =- ( amx(i) + apx(i) )                           

       amy(j) = amy_ad(i,j,k)
       apy(j) = apy_ad(i,j,k)
       acy(j) =- ( amy(j) + apy(j) )                          

       amz(k) = amz_ad(i,j,k)
       apz(k) = apz_ad(i,j,k)
       acz(k) =- ( amz(k) + apz(k) )                         

       rhs(k) = r(i,j,k) - var(i,JM,k)*amy(j)*(1-jum(i,j,k))  &
                         - var(i,JP,k)*apy(j)*(1-jup(i,j,k))  &
                         - var(IM,j,k)*amx(i)*(1-ium(i,j,k))  &
                         - var(IP,j,k)*apx(i)*(1-iup(i,j,k))

       amz(k) = amz(k)*riblank*(oned-kum(i,j,k))
       apz(k) = apz(k)*riblank*(oned-kup(i,j,k))
       acz(k) = oned + ( acx(i)+acy(j)+acz(k) ) * riblank &
                            *(oned-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL)) !  <--- transpose this properly RM 
       rhs(k) = rhs(k)*riblank + var(i,j,k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL)   !<--- transpose this properly RM
       
       
        IF ( boundary_formulation == GCM_METHOD .AND. &
             fresh_cell(i,j,k) == 1                   ) THEN
          iFr = i; jFr = j; kFr = k;
          CALL GCM_correct_rhs_ad(kDirection,iFr,jFr,kFr,ny,rhs,var)
        ENDIF ! fresh_cell
       
      ENDDO ! k

!     IF (bcz1 == BC_TYPE_PERIODIC .AND. &
!         bcz2 == BC_TYPE_PERIODIC) THEN
!        beta_p = amz(1)
!        alpha_p = apz(nzc)

!        CALL cyclic_tdma(amz,acz,apz,rhs,dummy,alpha_p,beta_p,1,nzc)
!     ELSE

        IF (bcz1 == BC_TYPE_PERIODIC .AND. &
            bcz2 == BC_TYPE_PERIODIC) THEN
           rhs(1) = rhs(1) - var(i,j,nzc)*amz(1)
           rhs(nzc) = rhs(nzc) - var(i,j,1)*apz(nzc)
        END IF

        CALL tdma(amz,acz,apz,rhs,dummy,1,nzc)
!     END IF

      DO k=1,nzc
!      var(i,j,k) = var(i,j,k) + omega_ad*(dummy(k)-var(i,j,k))
       var(i,j,k) = dummy(k)
      ENDDO !k

    ENDDO !j 
    ENDDO !i
    
   END SUBROUTINE itsolv_ad_z_AB
!----------------------------------------------------------
! 
!----------------------------------------------------------
subroutine identify_hybrid_cell         !add by yan
USE global_parameters
use boundary_arrays
USE flow_parameters
USE flow_arrays
USE multiuse_arrays
USE grid_arrays
use hybrid_cell_arrays

IMPLICIT NONE
integer :: i,j,k
integer :: im, ip, jm, jp, km, kp
!integer :: hybrid_mark(0:nx+1,0:ny+1,0:nz+1)

hybrid_mark=0

if(nbody_solid/=0)then
    do k=1,nzc
    do j=1,nyc
    do i=1,nxc
        iM = MAX(i-1,1)
        iP = MIN(i+1,nxc)
        jM = MAX(j-1,1)
        jP = MIN(j+1,nyc)
        kM = MAX(k-1,1)
        kP = MIN(k+1,nzc)
        if(iblank(i,j,k)/=iblank(ip,j,k).or.&
           iblank(i,j,k)/=iblank(im,j,k).or.&
           iblank(i,j,k)/=iblank(i,jp,k).or.&
           iblank(i,j,k)/=iblank(i,jm,k).or.&
           iblank(i,j,k)/=iblank(i,j,kp).or.&
           iblank(i,j,k)/=iblank(i,j,km))then
            if(iblank(i,j,k)==0)then
                hybrid_mark(i,j,k)=1
            end if   
        end if
    end do
    end do
    end do
end if



if(nbody_membrane/=0)then
    do k=1,nzc
    do j=1,nyc
    do i=1,nxc
        if(hybridMarkMemb(i,j,k)==1.or.hybridMarkMemb(i,j,k)==-1)then
            hybrid_mark(i,j,k)=1
        end if
    end do
    end do
    end do
end if



!    if(ntime==2)then
!    write(568,*) 'VARIABLES="X","Y","hybridMarkMemb","hybrid_mark"'
!    write(568,*) 'ZONE F=POINT, I=',nxc,' , J=',nyc
!    do j=1,nyc
!    do i=1,nxc
!        write(568,*) xc(i),yc(j),hybridMarkMemb(i,j,1),hybrid_mark(i,j,1)
!    end do
!    end do
!    end if

nhybrid=sum(hybrid_mark)
write(*,*) 'number of hybrid cells = ',nhybrid

end subroutine identify_hybrid_cell

!------------------------------------------------------------
!
!------------------------------------------------------------
subroutine update_hybrid_velocity           !add by yan
USE global_parameters
use boundary_arrays
USE flow_parameters
USE flow_arrays
USE multiuse_arrays
USE grid_arrays
use hybrid_cell_arrays

IMPLICIT NONE

integer :: i,j,k,n,iRow,ii,jj,kk
integer :: body_dist_min, closest_marker
integer :: h_mark(3)
!type(vector) :: vec
real(cgreal) :: temp_u1(0:nx+1,0:ny+1,0:nz+1), temp_v1(0:nx+1,0:ny+1,0:nz+1), temp_w1(0:nx+1,0:ny+1,0:nz+1)
real(cgreal) :: temp_u2(0:nx+1,0:ny+1,0:nz+1), temp_v2(0:nx+1,0:ny+1,0:nz+1), temp_w2(0:nx+1,0:ny+1,0:nz+1)
real(cgreal) :: uhb,vhb,whb,alpha1,alpha2,beta
!real(cgreal) :: temp_u3(0:nx+1,0:ny+1,0:nz+1), temp_v3(0:nx+1,0:ny+1,0:nz+1), temp_w3(0:nx+1,0:ny+1,0:nz+1)
real(cgreal) :: dist,distx,disty,distz,distb,ratiox,ratioy,ratioz
real(cgreal) :: xBIbeta,yBIbeta,zBIbeta
real(cgreal) :: xBIalpha1,yBIalpha1,zBIalpha1
real(cgreal) :: xBIalpha2,yBIalpha2,zBIalpha2
real(cgreal) :: xBIalpha3,yBIalpha3,zBIalpha3

temp_u1=u
temp_v1=v
temp_w1=w
u_bak=u
v_bak=v
w_bak=w
temp_u2=zero
temp_v2=zero
temp_w2=zero




do n=1,nhybrid

    xBIbeta=zero
    yBIbeta=zero
    zBIbeta=zero

    xBIalpha1=zero
    yBIalpha1=zero
    zBIalpha1=zero
    xBIalpha2=zero
    yBIalpha2=zero
    zBIalpha2=zero
    xBIalpha3=zero
    yBIalpha3=zero
    zBIalpha3=zero

    dist=zero
    distx=zero
    disty=zero
    distz=zero
    distb=zero

    i=ihybrid(n)
    j=jhybrid(n)
    k=khybrid(n)
    h_mark=(/i,j,k/)
    if(boundary_formulation == GCM_METHOD)then
        uhb=zero
        vhb=zero
        whb=zero
        do iRow=1,iRMax
            ii = iCIndex(n) + inccI(iRow)
            jj = jCIndex(n) + inccJ(iRow)
            kk = kCIndex(n) + inccK(iRow)

            IF ( ii /= i .OR. jj /= j .OR. kk /= k) THEN
                uhb = uhb + coeffD(iRow,n)* u(ii,jj,kk)
                vhb = vhb + coeffD(iRow,n)* v(ii,jj,kk)
                whb = whb + coeffD(iRow,n)* w(ii,jj,kk)
            ELSE
                uhb = uhb + coeffD(iRow,n)* uBIntercept(n)
                vhb = vhb + coeffD(iRow,n)* vBIntercept(n)
                whb = whb + coeffD(iRow,n)* wBIntercept(n)
            ENDIF ! ii

        end do
    
        temp_u2(i,j,k)=uhb
        temp_v2(i,j,k)=vhb
        temp_w2(i,j,k)=whb

        call calculate_dist(h_mark,dist)

        u(i,j,k)=temp_u1(i,j,k)*(oned-dist)+temp_u2(i,j,k)*dist
        v(i,j,k)=temp_v1(i,j,k)*(oned-dist)+temp_v2(i,j,k)*dist
        w(i,j,k)=temp_w1(i,j,k)*(oned-dist)+temp_w2(i,j,k)*dist
    else

!        if(iup(i,j,k)==1.or.ium(i,j,k)==1)then
!            temp_u2(i,j,k)=bcxu(i,j,k)
!            temp_v2(i,j,k)=bcxv(i,j,k)
!            temp_w2(i,j,k)=bcxw(i,j,k)
!        end if
!
!        if(jup(i,j,k)==1.or.jum(i,j,k)==1)then
!            temp_u2(i,j,k)=bcyu(i,j,k)
!            temp_v2(i,j,k)=bcyv(i,j,k)
!            temp_w2(i,j,k)=bcyw(i,j,k)
!        end if
!
!        if(ndim==dim_3d)then
!            if(kup(i,j,k)==1.or.kum(i,j,k)==1)then
!                temp_u2(i,j,k)=bczu(i,j,k)
!                temp_v2(i,j,k)=bczv(i,j,k)
!                temp_w2(i,j,k)=bczw(i,j,k)
!            end if
!        end if


!        call calculate_dist(h_mark,dist)


        CALL hybrid_Calc_BodyIntercept_Unstruc( i, j, k, xc(i), yc(j), zc(k),    &
                                        xBIbeta, yBIbeta, zBIbeta, closestElementHC(n) )
        distb=sqrt((xBIbeta-xc(i))**2+(yBIbeta-yc(j))**2+(zBIbeta-zc(k))**2)
        if(iup(i,j,k)==1) then
            CALL hybrid_Calc_BodyIntercept_Unstruc( i+1, j, k, xc(i+1), yc(j), zc(k),    &
                                        xBIalpha1, yBIalpha1, zBIalpha1, closestElementR1(n) )
            distx=sqrt((xBIalpha1-xc(i+1))**2+(yBIalpha1-yc(j))**2+(zBIalpha1-zc(k))**2)
        end if
        if(ium(i,j,k)==1) then
            CALL hybrid_Calc_BodyIntercept_Unstruc( i-1, j, k, xc(i-1), yc(j), zc(k),    &
                                        xBIalpha1, yBIalpha1, zBIalpha1, closestElementR1(n) )
            distx=sqrt((xBIalpha1-xc(i-1))**2+(yBIalpha1-yc(j))**2+(zBIalpha1-zc(k))**2)
        end if
        if(jup(i,j,k)==1) then
            CALL hybrid_Calc_BodyIntercept_Unstruc( i, j+1, k, xc(i), yc(j+1), zc(k),    &
                                        xBIalpha2, yBIalpha2, zBIalpha2, closestElementR1(n) )
            disty=sqrt((xBIalpha2-xc(i))**2+(yBIalpha2-yc(j+1))**2+(zBIalpha2-zc(k))**2)
        end if
        if(jum(i,j,k)==1) then
            CALL hybrid_Calc_BodyIntercept_Unstruc( i, j-1, k, xc(i), yc(j-1), zc(k),    &
                                        xBIalpha2, yBIalpha2, zBIalpha2, closestElementR1(n) )
            disty=sqrt((xBIalpha2-xc(i))**2+(yBIalpha2-yc(j-1))**2+(zBIalpha2-zc(k))**2)
        end if
        if(ndim==dim_3d)then
            if(kup(i,j,k)==1) then
                CALL hybrid_Calc_BodyIntercept_Unstruc( i, j, k+1, xc(i), yc(j), zc(k+1),    &
                                            xBIalpha3, yBIalpha3, zBIalpha3, closestElementR1(n) )
                distz=sqrt((xBIalpha3-xc(i))**2+(yBIalpha3-yc(j))**2+(zBIalpha3-zc(k+1))**2)
            end if
            if(kum(i,j,k)==1) then
                CALL hybrid_Calc_BodyIntercept_Unstruc( i, j, k-1, xc(i), yc(j), zc(k-1),    &
                                            xBIalpha3, yBIalpha3, zBIalpha3, closestElementR1(n) )
                distz=sqrt((xBIalpha3-xc(i))**2+(yBIalpha3-yc(j))**2+(zBIalpha3-zc(k-1))**2)
            end if
        end if


        ratiox=distx/(distx+distb)
        ratioy=disty/(disty+distb)

        if(ndim==dim_3d)then
            ratioz=distz/(distz+distb)
        end if
!        ratiox=1.0d0
!        ratioy=1.0d0
!        ratioz=1.0d0
!        if(n==33)then
!            write(*,*) 'n=29'
!            write(*,*) 'n=29'
!        end if

!        if(iup(i,j,k)==1) temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-u(i-1,j,k)
!        if(ium(i,j,k)==1) temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-u(i+1,j,k)
!        if(jup(i,j,k)==1) temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-v(i,j-1,k)
!        if(jum(i,j,k)==1) temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-v(i,j+1,k)
!        if(ndim==dim_3d)then
!            if(kup(i,j,k)==1) temp_w2(i,j,k)=2.0d0*bczw(i,j,k)-w(i,j,k-1)
!            if(kum(i,j,k)==1) temp_w2(i,j,k)=2.0d0*bczw(i,j,k)-w(i,j,k+1)
!        end if

        if((iup(i,j,k)==1.or.ium(i,j,k)==1).and.(jup(i,j,k)==0.and.jum(i,j,k)==0))then
            if(iup(i,j,k)==1)then
                !temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-u(i-1,j,k)
                temp_v2(i,j,k)=2.0d0*bcxv(i,j,k)-v(i,j+1,k)
                !print*, bcxu(i,j,k),bcxv(i,j,k)
            else
                !temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-u(i+1,j,k)
                temp_v2(i,j,k)=2.0d0*bcxv(i,j,k)-v(i,j-1,k)
                !print*, bcxu(i,j,k),bcxv(i,j,k)
            end if
            !u(i,j,k)=temp_u1(i,j,k)*(oned-ratiox)+temp_u2(i,j,k)*ratiox
            v(i,j,k)=temp_v1(i,j,k)*(oned-ratiox)+temp_v2(i,j,k)*ratiox
!            u(i,j,k)=temp_u1(i,j,k)*(oned-dist)+temp_u2(i,j,k)*dist
!            v(i,j,k)=temp_v1(i,j,k)*(oned-dist)+temp_v2(i,j,k)*dist
        end if
        if((jup(i,j,k)==1.or.jum(i,j,k)==1).and.(iup(i,j,k)==0.and.ium(i,j,k)==0))then
            if(jup(i,j,k)==1)then
                !temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-v(i+1,j-1,k)
                temp_u2(i,j,k)=2.0d0*bcyu(i,j,k)-u(i+1,j,k)
                !print*, bcyv(i,j,k),bcyu(i,j,k)
            else
                !temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-v(i-1,j+1,k)
                temp_u2(i,j,k)=2.0d0*bcyu(i,j,k)-u(i-1,j,k)
                !print*, bcyv(i,j,k),bcyu(i,j,k)
            end if
            u(i,j,k)=temp_u1(i,j,k)*(oned-ratioy)+temp_u2(i,j,k)*ratioy
            !v(i,j,k)=temp_v1(i,j,k)*(oned-ratioy)+temp_v2(i,j,k)*ratioy
!            u(i,j,k)=temp_u1(i,j,k)*(oned-dist)+temp_u2(i,j,k)*dist
!            v(i,j,k)=temp_v1(i,j,k)*(oned-dist)+temp_v2(i,j,k)*dist
        end if
!        if((iup(i,j,k)==1.or.ium(i,j,k)==1).and.(jup(i,j,k)==1.or.jum(i,j,k)==1))then
!            if(iup(i,j,k)==1.and.jup(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-u(i-1,j,k)
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-v(i,j-1,k)
!            else if(iup(i,j,k)==1.and.jum(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-u(i-1,j,k)
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-v(i,j+1,k)
!            else if(ium(i,j,k)==1.and.jup(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-u(i+1,j,k)
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-v(i,j-1,k)
!            else if(ium(i,j,k)==1.and.jum(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-u(i+1,j,k)
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-v(i,j+1,k)
!            end if
!            
!            u(i,j,k)=temp_u1(i,j,k)*(oned-ratiox)+temp_u2(i,j,k)*ratiox
!            v(i,j,k)=temp_v1(i,j,k)*(oned-ratioy)+temp_v2(i,j,k)*ratioy
!!            u(i,j,k)=temp_u1(i,j,k)*(oned-dist)+temp_u2(i,j,k)*dist
!!            v(i,j,k)=temp_v1(i,j,k)*(oned-dist)+temp_v2(i,j,k)*dist
!        end if
!
!        if(ndim==dim_3d)then
!            w(i,j,k)=temp_w1(i,j,k)*(oned-ratioz)+temp_w2(i,j,k)*ratioz
!!            w(i,j,k)=temp_w1(i,j,k)*(oned-dist)+temp_w2(i,j,k)*dist
!        end if

!        if((iup(i,j,k)==1.or.ium(i,j,k)==1).and.(jup(i,j,k)==0.or.jum(i,j,k)==0))then
!            if(iup(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-half*(u(i-1,j,k)+u(i,j+1,k))
!                temp_v2(i,j,k)=2.0d0*bcxv(i,j,k)-half*(v(i-1,j,k)+v(i,j+1,k))
!            else
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-half*(u(i+1,j,k)+u(i,j+1,k))
!                temp_v2(i,j,k)=2.0d0*bcxv(i,j,k)-half*(v(i+1,j,k)+v(i,j+1,k))
!            end if
!            u(i,j,k)=temp_u1(i,j,k)*(oned-ratiox)+temp_u2(i,j,k)*ratiox
!            v(i,j,k)=temp_v1(i,j,k)*(oned-ratiox)+temp_v2(i,j,k)*ratiox
!        end if
!        if((jup(i,j,k)==1.or.jum(i,j,k)==1).and.(iup(i,j,k)==0.or.ium(i,j,k)==0))then
!            if(jup(i,j,k)==1)then
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-half*(v(i,j-1,k)+v(i-1,j,k))
!                temp_u2(i,j,k)=2.0d0*bcyu(i,j,k)-half*(u(i,j-1,k)+u(i-1,j,k))
!            else
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-half*(v(i,j+1,k)+v(i-1,j,k))
!                temp_u2(i,j,k)=2.0d0*bcyu(i,j,k)-half*(u(i,j+1,k)+u(i-1,j,k))
!            end if
!            u(i,j,k)=temp_u1(i,j,k)*(oned-ratioy)+temp_u2(i,j,k)*ratioy
!            !u(i,j,k)=zero
!            v(i,j,k)=temp_v1(i,j,k)*(oned-ratioy)+temp_v2(i,j,k)*ratioy
!        end if
!        if((iup(i,j,k)==1.or.ium(i,j,k)==1).and.(jup(i,j,k)==1.or.jum(i,j,k)==1))then
!            if(iup(i,j,k)==1.and.jup(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-half*(u(i-1,j,k)+u(i,j-1,k))
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-half*(v(i-1,j,k)+v(i,j-1,k))
!            else if(iup(i,j,k)==1.and.jum(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-half*(u(i-1,j,k)+u(i,j+1,k))
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-half*(v(i-1,j,k)+v(i,j+1,k))
!            else if(ium(i,j,k)==1.and.jup(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-half*(u(i+1,j,k)+u(i,j-1,k))
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-half*(v(i+1,j,k)+v(i,j-1,k))
!            else if(ium(i,j,k)==1.and.jum(i,j,k)==1)then
!                temp_u2(i,j,k)=2.0d0*bcxu(i,j,k)-half*(u(i+1,j,k)+u(i,j+1,k))
!                temp_v2(i,j,k)=2.0d0*bcyv(i,j,k)-half*(v(i+1,j,k)+v(i,j+1,k))
!            end if
!            
!            u(i,j,k)=temp_u1(i,j,k)*(oned-ratiox)+temp_u2(i,j,k)*ratiox
!            v(i,j,k)=temp_v1(i,j,k)*(oned-ratioy)+temp_v2(i,j,k)*ratioy
!        end if
!
!        if(ndim==dim_3d)then
!            w(i,j,k)=temp_w1(i,j,k)*(oned-ratioz)+temp_w2(i,j,k)*ratioz
!        end if
        




    end if

end do



end subroutine update_hybrid_velocity
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine calculate_dist(h_mark,dist)          !add by yan
USE global_parameters
use boundary_arrays
USE flow_parameters
USE flow_arrays
USE multiuse_arrays
USE grid_arrays

IMPLICIT NONE

integer,intent(in) :: h_mark(3)
real(cgreal),intent(out) :: dist
real(cgreal) :: dip,dim,djp,djm,dkp,dkm
!type(vector) :: vip,vim,vjp,vjm,vkp,vkm
real(cgreal) :: vip(3),vim(3),vjp(3),vjm(3),vkp(3),vkm(3)
integer :: lip,lim,ljp,ljm,lkp,lkm
integer :: sum_lable
integer :: body_dist_min, closest_marker
!integer :: i,j,k
lip=0
lim=0
ljp=0
ljm=0
lkp=0
lkm=0
sum_lable=0

if(iblank(h_mark(1)+1,h_mark(2),h_mark(3))==1.or.hybridMarkMemb(h_mark(1)+1,h_mark(2),h_mark(3))*hybridMarkMemb(h_mark(1),h_mark(2),h_mark(3))==-1)then
    lip=1
    vip(1)=xc(h_mark(1)+1)
    vip(2)=yc(h_mark(2))
    vip(3)=zc(h_mark(3))
end if

if(iblank(h_mark(1)-1,h_mark(2),h_mark(3))==1.or.hybridMarkMemb(h_mark(1)-1,h_mark(2),h_mark(3))*hybridMarkMemb(h_mark(1),h_mark(2),h_mark(3))==-1)then
    lim=1
    vim(1)=xc(h_mark(1)-1)
    vim(2)=yc(h_mark(2))
    vim(3)=zc(h_mark(3))
end if

if(iblank(h_mark(1),h_mark(2)+1,h_mark(3))==1.or.hybridMarkMemb(h_mark(1),h_mark(2)+1,h_mark(3))*hybridMarkMemb(h_mark(1),h_mark(2),h_mark(3))==-1)then
    ljp=1
    vjp(1)=xc(h_mark(1))
    vjp(2)=yc(h_mark(2)+1)
    vjp(3)=zc(h_mark(3))
end if

if(iblank(h_mark(1),h_mark(2)-1,h_mark(3))==1.or.hybridMarkMemb(h_mark(1),h_mark(2)-1,h_mark(3))*hybridMarkMemb(h_mark(1),h_mark(2),h_mark(3))==-1)then
    ljm=1
    vjm(1)=xc(h_mark(1))
    vjm(2)=yc(h_mark(2)-1)
    vjm(3)=zc(h_mark(3))
end if

if(ndim==dim_3d)then
    if(iblank(h_mark(1),h_mark(2),h_mark(3)+1)==1.or.hybridMarkMemb(h_mark(1),h_mark(2),h_mark(3)+1)*hybridMarkMemb(h_mark(1),h_mark(2),h_mark(3))==-1)then
        lkp=1
        vkp(1)=xc(h_mark(1))
        vkp(2)=yc(h_mark(2))
        vkp(3)=zc(h_mark(3)+1)
    end if

    if(iblank(h_mark(1),h_mark(2),h_mark(3)-1)==1.or.hybridMarkMemb(h_mark(1),h_mark(2),h_mark(3)-1)*hybridMarkMemb(h_mark(1),h_mark(2),h_mark(3))==-1)then
        lkm=1
        vkm(1)=xc(h_mark(1))
        vkm(2)=yc(h_mark(2))
        vkm(3)=zc(h_mark(3)-1)
    end if
end if

sum_lable=lip+lim+ljp+ljm+lkp+lkm

if(ndim==dim_2d)then
    if(sum_lable>=3)then
        write(*,*) 'Error! More than two ghost nodes are found in 2D case!'
        stop
    end if
else if(ndim==dim_3d)then
    if(sum_lable>=4)then
        write(*,*) 'Error! More than three ghost nodes are found in 3D case!'
        stop
    end if
end if

if(lip==0)then
    dip=zero
else
    call find_closest_element(vip,body_dist_min,closest_marker)
    dip=sqrt((xbodymarker(body_dist_min,closest_marker)-vip(1))**twod+(ybodymarker(body_dist_min,closest_marker)-vip(2))**twod+&
        (zbodymarker(body_dist_min,closest_marker)-vip(3))**twod)
end if

if(lim==0)then
    dim=zero
else
    call find_closest_element(vim,body_dist_min,closest_marker)
    dim=sqrt((xbodymarker(body_dist_min,closest_marker)-vim(1))**twod+(ybodymarker(body_dist_min,closest_marker)-vim(2))**twod+&
        (zbodymarker(body_dist_min,closest_marker)-vim(3))**twod)
end if

if(ljp==0)then
    djp=zero
else
    call find_closest_element(vjp,body_dist_min,closest_marker)
    djp=sqrt((xbodymarker(body_dist_min,closest_marker)-vjp(1))**twod+(ybodymarker(body_dist_min,closest_marker)-vjp(2))**twod+&
        (zbodymarker(body_dist_min,closest_marker)-vjp(3))**twod)
end if

if(ljm==0)then
    djm=zero
else
    call find_closest_element(vjm,body_dist_min,closest_marker)
    djm=sqrt((xbodymarker(body_dist_min,closest_marker)-vjm(1))**twod+(ybodymarker(body_dist_min,closest_marker)-vjm(2))**twod+&
        (zbodymarker(body_dist_min,closest_marker)-vjm(3))**twod)
end if

if(ndim==dim_3d)then
    if(lkp==0)then
        dkp=zero
    else
        call find_closest_element(vkp,body_dist_min,closest_marker)
        dkp=sqrt((xbodymarker(body_dist_min,closest_marker)-vkp(1))**twod+(ybodymarker(body_dist_min,closest_marker)-vkp(2))**twod+&
            (zbodymarker(body_dist_min,closest_marker)-vkp(3))**twod)
    end if

    if(lkm==0)then
        dkm=zero
    else
        call find_closest_element(vkm,body_dist_min,closest_marker)
        dkm=sqrt((xbodymarker(body_dist_min,closest_marker)-vkm(1))**twod+(ybodymarker(body_dist_min,closest_marker)-vkm(2))**twod+&
            (zbodymarker(body_dist_min,closest_marker)-vkm(3))**twod)
    end if
end if

if(ndim==dim_2d)then
    dist=sqrt((dip/dxc(h_mark(1)+1))**twod+(dim/dxc(h_mark(1)))**twod+(djp/dyc(h_mark(2)+1))**twod+(djm/dyc(h_mark(2)))**twod)
else if(ndim==dim_3d)then
    dist=sqrt((dip/dxc(h_mark(1)+1))**twod+(dim/dxc(h_mark(1)))**twod+(djp/dyc(h_mark(2)+1))**twod+(djm/dyc(h_mark(2)))**twod+(dkp/dzc(h_mark(3)+1))**twod+(dkm/dzc(h_mark(3)))**twod)
end if


end subroutine calculate_dist
!-------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------
subroutine hybrid_cell_interpolation                !add by yan
USE global_parameters
USE flow_parameters
USE flow_arrays
USE grid_arrays
USE boundary_arrays
USE GCM_arrays
USE hybrid_cell_arrays
USE unstructured_surface_arrays
    
IMPLICIT NONE

INTEGER :: i,iBody,iRow,j,k,m,n
INTEGER :: iG, jG, kG, nbdr, iCIndx, jCIndx, kCIndx , iCIndxS, jCIndxS
INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
INTEGER :: iM, iP, jM, jP, kM, kP
INTEGER :: iRange, jRange, kRange
INTEGER :: node1,node2,node3

REAL(KIND=CGREAL) :: cosTheta,dsIntercept,sinTheta,               &
                    xBI,xBIN,xGC,xIP,xIPS,yBI,yBIN,yGC,yIP,yIPS,zBI,zBIN,zGC,zIP, &
                    minProbeLengthShear,slopeX, slopeY, slopeZ, maxDelta

REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: coeffDirc, coeffNeum

ALLOCATE( coeffDirc(iRMax) )
ALLOCATE( coeffNeum(iRMax) )

nhybrid=sum(hybrid_mark)

if(ntime>ntime_start+1)then
    call hybrid_mem_deallocate
end if

call hybrid_mem_allocate

nbdr=0
do k=1,nzc
do j=1,nyc
do i=1,nxc
    if(hybrid_mark(i,j,k)==1)then
        nbdr=nbdr+1
        ihybrid(nbdr)=i
        jhybrid(nbdr)=j
        khybrid(nbdr)=k
    end if
end do
end do
end do

do n=1,nhybrid
    iG = ihybrid(n) 
    jG = jhybrid(n)
    kG = khybrid(n)
      
    iCIndex(n) = -1
    jCIndex(n) = -1
    kCIndex(n) = -1
    
    iBody = 1
    
    xGC = xc(iG)
    yGC = yc(jG)
    zGC = zc(kG)
    
    CALL hybrid_Calc_BodyIntercept_Unstruc( iG, jG, kG, xGC, yGC, zGC,    &
                                        xBI, yBI, zBI, closestElementHC(n) )

    xBIntercept(n) = xBI
    yBIntercept(n) = yBI
    zBIntercept(n) = zBI

    dsIntercept = SQRT( (xGC-xBI)**2 + (yGC-yBI)**2 + (zGC-zBI)**2 )
    
    IF ( dsIntercept > SQRT(dxc(iG)**2 + dyc(jG)**2 + dzc(kG)**2) ) THEN
        PRINT*,' Normal intercept for hybrid cells might not be correct!'
        IF (dsIntercept/SQRT(dxc(iG)**2 + dyc(jG)**2 + dzc(kG)**2) > 2.0_CGREAL) THEN
            PRINT*,'Intercept is too long'
            STOP
        ENDIF
    ENDIF

    if(xGC<=xBI)then
        iCIndex(n)=iG-1
    else
        iCIndex(n)=iG
    end if
    if(yGC<=yBI)then
        jCIndex(n)=jG-1
    else
        jCIndex(n)=jG
    end if
    if(zGC<=zBI)then
        kCIndex(n)=kG-1
    else
        kCIndex(n)=kG
    end if

    iCIndx = iCIndex(n)
    jCIndx = jCIndex(n)
    kCIndx = kCIndex(n)
 
    xBIN = triElemNormx(iBody,closestElementHC(n))
    yBIN = triElemNormy(iBody,closestElementHC(n))
    zBIN = triElemNormz(iBody,closestElementHC(n))

    CALL hybrid_Calc_vanMatrixDN( iG, jG, kG, iCIndx, jCIndx, kCIndx,             &
                            xGC, yGC, zGC, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                            coeffDirc, coeffNeum      )
 
    coeffD(1:iRMax,n) = coeffDirc(1:iRMax)
    coeffN(1:iRMax,n) = coeffNeum(1:iRMax)

    CALL hybrid_calc_BIVelocity_Unstruc( iG, jG, kG, xGC, yGC, zGC, closestElementHC(n), &
                                           uBIntercept(n),        &
                                           vBIntercept(n),        &
                                           wBIntercept(n)         )

end do

deallocate(coeffDirc)
deallocate(coeffNeum)
!call update_hybrid_velocity


end subroutine hybrid_cell_interpolation
!--------------------------------------------------------------
!
!--------------------------------------------------------------
SUBROUTINE hybrid_Calc_VanMatrixDN( iG, jG, kG, iCIndx,jCIndx, kCIndx,      &            !add by yan
                                  xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                                  coeffDirc, coeffNeum      )
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    !USE gcm_arrays
    use hybrid_cell_arrays

    IMPLICIT NONE

!... parameters variables

    INTEGER, INTENT(IN)           :: iG, jG, kG, iCIndx,jCIndx, kCIndx
    REAL(KIND=CGREAL), INTENT(IN) :: xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN
    REAL(KIND=CGREAL), DIMENSION(iRMax), INTENT(OUT) :: coeffDirc, &
                                                          coeffNeum

!... loop variables

    INTEGER :: i,j,k,iRow
    INTEGER :: info

!... local variables
    
    REAL(KIND=CGREAL) :: rCond, xC1,xC2,xC3, xN1,xN2,xN3

!*****************************************************************************************
  
!   |-------|-------|---/---|-------|--         N : Nth ghost point
!   |   ii  |  iii  |  *    |       |           * : markers
!   |   0...|...O   | / .   |   .   |           O : other nodes used in bilinear interpolation
!   |   .   |   .   |*      |       |           + : probe tip (Image Point) 
!   |---.---|--+.---/-------|-------|--
!   |   .   |   .  *|       |       |
!   |   0...|. .O / |   N   |   .   |
!   |   i   |  iv*  |       |       |
!   |-------| --/ --|-------|-------|--

! interpolant      U = a X X X  + b X X  + c X X  + d X X
!                         1 2 3      1 2      1 3      2 3
!
!                    + e X  + f X + g X  + h
!                         1      2     3
!
!
!         [  X X     X     X   1  ]  [   ]     [     ] 
!      i  [   1 2     1     2     ]  [ a ]     [ U   ]
!         [                       ]  [   ]     [  i  ]
!         [  X X     X     X   1  ]  [   ]     [     ]
!      ii [   1 2     1     2     ]  [ b ]     [ U   ]
!         [                       ]  [   ]  =  [  ii ]
!     iii [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ c ]     [ U   ]
!         [                       ]  [   ]     [  iii]
!     iv  [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ d ]     [ U   ]
!         [                       ]  [   ]     [  iv ]
!
!
!   Van Matrix For Dirichlet conditions at Intersection Point (N)
!
!         [  X X     X     X   1  ]  [   ]     [     ] 
!      i  [   1 2     1     2     ]  [ a ]     [ U   ]
!         [                       ]  [   ]     [  i  ]
!         [  X X     X     X   1  ]  [   ]     [     ]
!      ii [   1 2     1     2     ]  [ b ]     [ U   ]
!         [                       ]  [   ]  =  [  ii ]
!     iii [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ c ]     [ U   ]
!         [                       ]  [   ]     [  iii]
!      N  [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ d ]     [ U   ]
!         [                       ]  [   ]     [  N  ]
!
!   Van Matrix For Neumann conditions at Intersection point (N)
!    B1 = n_x, B2 = n_y (components of normal vectors)
!    F_m = value of normal derivative 
!
!         [  X X           X     X   1  ]  [   ]     [     ] 
!      i  [   1 2           1     2     ]  [ a ]     [ U   ]
!         [                             ]  [   ]     [  i  ]
!         [  X X           X     X   1  ]  [   ]     [     ]
!      ii [   1 2           1     2     ]  [ b ]     [ U   ]
!         [                             ]  [   ]  =  [  ii ]
!     iii [  X X           X     X   1  ]  [   ]     [     ]
!         [   1 2           1     2     ]  [ c ]     [ U   ]
!         [                             ]  [   ]     [  iii]
!      N  [  B X  + B X    B     B   0  ]  [   ]     [     ]
!         [   1 2    2  1   1     2     ]  [ d ]     [ F   ]
!         [                             ]  [   ]     [  N  ]
!

    DO iRow= 1, iRMax
      i  = iCIndx + inccI(iRow)
      j  = jCIndx + inccJ(iRow)
      k  = kCIndx + inccK(iRow)

      xC1 = xc(i)
      xC2 = yc(j)
      xC3 = zc(k)

!-- Construct Vandermonde Matrices

!--- Dirichlet conditions for velocity field

      vanMD(iRow,1) = xC1*xC2*xC3
      vanMD(iRow,2) = xC1*xC2
      vanMD(iRow,3) = xC1*xC3
      vanMD(iRow,4) = xC2*xC3
      vanMD(iRow,5) = xC1
      vanMD(iRow,6) = xC2
      vanMD(iRow,7) = xC3
      vanMD(iRow,8) = oned

!--- Neumann conditions for pressure field


      vanMN(iRow,1) = xC1*xC2*xC3
      vanMN(iRow,2) = xC1*xC2
      vanMN(iRow,3) = xC1*xC3
      vanMN(iRow,4) = xC2*xC3
      vanMN(iRow,5) = xC1
      vanMN(iRow,6) = xC2
      vanMN(iRow,7) = xC3
      vanMN(iRow,8) = oned

!-- Correct For Ghost node part of cell formation, switch to Body Intercept point

      IF ( i==iG .AND. j == jG  .AND. k== kG) THEN
        xC1 = xBI
        xC2 = yBI
        xC3 = zBI
        xN1 = xBIN
        xN2 = yBIN
        xN3 = zBIN

        vanMD(iRow,1) = xC1*xC2*xC3
        vanMD(iRow,2) = xC1*xC2
        vanMD(iRow,3) = xC1*xC3
        vanMD(iRow,4) = xC2*xC3
        vanMD(iRow,5) = xC1
        vanMD(iRow,6) = xC2
        vanMD(iRow,7) = xC3
        vanMD(iRow,8) = oned

        vanMN(iRow,1) = xN1*xC2*XC3 + xN2*xC1*XC3 + xN3*XC1*XC2
        vanMN(iRow,2) = xN1*xC2 + xN2*xC1
        vanMN(iRow,3) = xN1*xC3 + xN3*xC1
        vanMN(iRow,4) = xN2*xC3 + xN3*xC2
        vanMN(iRow,5) = xN1
        vanMN(iRow,6) = xN2
        vanMN(iRow,7) = xN3
        vanMN(iRow,8) = zero

      ENDIF ! i
    ENDDO ! iRow		

! Compute inverse of Vandermonde Matrices

    CALL DGETRF(8, 8, vanMD,8,iPvtt, info) 
    CALL DGETRI(8, vanMD,8,iPvtt,workk, 8, info) 

    CALL DGETRF(8, 8, vanMN,8,iPvtt, info)
    CALL DGETRI(8, vanMN,8,iPvtt,workk, 8, info)

! Load Coeff-Matrices

    DO iRow = 1, iRMax
      coeffDirc(iRow) = vanMD(1,iRow)*xIP*yIP*zIP  &
                         + vanMD(2,iRow)*xIP*yIP      &
                         + vanMD(3,iRow)*xIP*zIP      &
                         + vanMD(4,iRow)*yIP*zIP      &
                         + vanMD(5,iRow)*xIP          &
                         + vanMD(6,iRow)*yIP          &
                         + vanMD(7,iRow)*zIP          &
                         + vanMD(8,iRow)

      coeffNeum(iRow) = vanMN(1,iRow)*xIP*yIP*zIP  &
                         + vanMN(2,iRow)*xIP*yIP      &
                         + vanMN(3,iRow)*xIP*zIP      &
                         + vanMN(4,iRow)*yIP*zIP      &
                         + vanMN(5,iRow)*xIP          &
                         + vanMN(6,iRow)*yIP          &
                         + vanMN(7,iRow)*zIP          &
                         + vanMN(8,iRow)
    ENDDO ! iRow 

  END SUBROUTINE hybrid_Calc_VanMatrixDN

!--------------------------------------------------------------
!
!--------------------------------------------------------------
subroutine hybrid_mem_allocate              !add by yan
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE hybrid_cell_arrays
USE unstructured_surface_arrays
    
IMPLICIT NONE

iRMax=8

ALLOCATE(ihybrid(1:nhybrid))
ALLOCATE(jhybrid(1:nhybrid)) 
ALLOCATE(khybrid(1:nhybrid)) 
ALLOCATE(closestElementHC(nhybrid))
ALLOCATE(closestElementR1(nhybrid))
ALLOCATE(closestElementR2(nhybrid))
ALLOCATE(closestElementR3(nhybrid))
ALLOCATE(xBIntercept(nhybrid))
ALLOCATE(yBIntercept(nhybrid))
ALLOCATE(zBIntercept(nhybrid))
ALLOCATE(uBIntercept(nhybrid))
ALLOCATE(vBIntercept(nhybrid))
ALLOCATE(wBIntercept(nhybrid))
ALLOCATE(coeffD(iRMax,nhybrid))
ALLOCATE(coeffN(iRMax,nhybrid))
ALLOCATE(iCIndex(nhybrid))
ALLOCATE(jCIndex(nhybrid))
ALLOCATE(kCIndex(nhybrid))

ihybrid            = 0
jhybrid            = 0 
khybrid            = 0

coeffD          = zero
coeffN          = zero

xBIntercept     = zero
yBIntercept     = zero 
zBIntercept     = zero 
uBIntercept     = zero
vBIntercept     = zero 
wBIntercept     = zero 

closestElementHC   = zero

end subroutine hybrid_mem_allocate
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
subroutine hybrid_mem_allocate_static              !add by yan
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE hybrid_cell_arrays
USE unstructured_surface_arrays
    
IMPLICIT NONE

iRMax=8

    ALLOCATE(inccI(iRMax))
    ALLOCATE(inccJ(iRMax))
    ALLOCATE(inccK(iRMax))
    ALLOCATE(iPvtt(iRMax))
    ALLOCATE(workk(iRMax))
    ALLOCATE(vanMD(iRMax,iRMax))
    ALLOCATE(vanMN(iRMax,iRMax))

! These allow us to define stencil image point
! Assumed clockwise from lower left corner.

      inccI(1) = 0
      inccJ(1) = 0
      inccK(1) = 0
      inccI(2) = 0
      inccJ(2) = 1
      inccK(2) = 0
      inccI(3) = 1
      inccJ(3) = 1
      inccK(3) = 0
      inccI(4) = 1
      inccJ(4) = 0
      inccK(4) = 0

      inccI(5) = 0  !----?????????
      inccJ(5) = 0
      inccK(5) = 1
      inccI(6) = 0
      inccJ(6) = 1
      inccK(6) = 1
      inccI(7) = 1
      inccJ(7) = 1
      inccK(7) = 1
      inccI(8) = 1
      inccJ(8) = 0
      inccK(8) = 1  ! ----????????

end subroutine hybrid_mem_allocate_static
!-------------------------------------------------------------
!
!-------------------------------------------------------------
SUBROUTINE hybrid_calc_BIVelocity_Unstruc( iGBI, jGBI, kGBI, xGBI, yGBI, zGBI, closestElementGBI,          &
                                          uGBI, vGBI, wGBI )

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE gcm_arrays
    use hybrid_cell_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... parameters

    INTEGER,           INTENT(IN)  :: iGBI, jGBI, kGBI, closestElementGBI 
    REAL(KIND=CGREAL), INTENT(IN)  :: xGBI, yGBI,zGBI
    REAL(KIND=CGREAL), INTENT(OUT) :: uGBI, vGBI, wGBI
!   
!... loop variables

    INTEGER :: i

!... local variables

    INTEGER                           :: iBody,node1,node2,node3,nMarker
    INTEGER                           :: info
    REAL(KIND=CGREAL)                 :: cX, cY, cZ, cC, rCond
    REAL(KIND=CGREAL), DIMENSION(4,4) :: vanTri
    REAL(KIND=CGREAL), DIMENSION(4)   :: rhsTri

! use the following approach
!
!  u = a x + b y + c z + d 
!  (a,b,c,d) determined by using four conditions
!   u = u(i) at ith node for i=1,3
!   GRAD(u) . n = 0  where n is normal to plane of triangle.
!  
!******************************************************************************

    iBody   = 1
    nMarker = nPtsBodyMarker(iBody)
!  
!   assume linear variation of velocity across element and then compute value at intercept.

    node1   = triElemNeig(iBody,1,closestElementGBI)
    node2   = triElemNeig(iBody,2,closestElementGBI)
    node3   = triElemNeig(iBody,3,closestElementGBI)
!
    vanTri(1,1) = xBodyMarker(iBody,node1) 
    vanTri(1,2) = yBodyMarker(iBody,node1) 
    vanTri(1,3) = zBodyMarker(iBody,node1) 
    vanTri(1,4) = oned

    vanTri(2,1) = xBodyMarker(iBody,node2) 
    vanTri(2,2) = yBodyMarker(iBody,node2) 
    vanTri(2,3) = zBodyMarker(iBody,node2) 
    vanTri(2,4) = oned

    vanTri(3,1) = xBodyMarker(iBody,node3) 
    vanTri(3,2) = yBodyMarker(iBody,node3) 
    vanTri(3,3) = zBodyMarker(iBody,node3) 
    vanTri(3,4) = oned

    vanTri(4,1) = triElemNormx(iBody,closestElementGBI) 
    vanTri(4,2) = triElemNormy(iBody,closestElementGBI) 
    vanTri(4,3) = triElemNormz(iBody,closestElementGBI) 
    vanTri(4,4) = zero

    CALL DGETRF(4, 4, vanTri,4,iPvtt, info)
    CALL DGETRI(4, vanTri,4,iPvtt,workk, 4, info)

! compute uGBI
    rhsTri(1) = uBodyMarker(iBody,node1)
    rhsTri(2) = uBodyMarker(iBody,node2)
    rhsTri(3) = uBodyMarker(iBody,node3)
    rhsTri(4) = zero

    cX = zero  
    cY = zero  
    cZ = zero  
    cC = zero  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO

    uGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

! compute vGBI
    rhsTri(1) = vBodyMarker(iBody,node1)
    rhsTri(2) = vBodyMarker(iBody,node2)
    rhsTri(3) = vBodyMarker(iBody,node3)
    rhsTri(4) = zero

    cX = zero  
    cY = zero  
    cZ = zero  
    cC = zero  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO

    vGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

! compute wGBI
    rhsTri(1) = wBodyMarker(iBody,node1)
    rhsTri(2) = wBodyMarker(iBody,node2)
    rhsTri(3) = wBodyMarker(iBody,node3)
    rhsTri(4) = zero

    cX = zero  
    cY = zero  
    cZ = zero  
    cC = zero  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO
  
    wGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

  END SUBROUTINE hybrid_calc_BIVelocity_Unstruc
!-------------------------------------------------------------
!
!-------------------------------------------------------------
subroutine hybrid_mem_deallocate            !add by yan
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE hybrid_cell_arrays
USE unstructured_surface_arrays
    
IMPLICIT NONE

DEALLOCATE(ihybrid)
DEALLOCATE(jhybrid) 
DEALLOCATE(khybrid) 
DEALLOCATE(closestElementHC)
DEALLOCATE(closestElementR1)
DEALLOCATE(closestElementR2)
DEALLOCATE(closestElementR3)
DEALLOCATE(xBIntercept)
DEALLOCATE(yBIntercept)
DEALLOCATE(zBIntercept)
DEALLOCATE(coeffD)
DEALLOCATE(coeffN)


DEALLOCATE(uBIntercept)
DEALLOCATE(vBIntercept)
DEALLOCATE(wBIntercept)

DEALLOCATE(iCIndex)
DEALLOCATE(jCIndex)
DEALLOCATE(kCIndex)


end subroutine hybrid_mem_deallocate



!------------------------------------------
  SUBROUTINE hybrid_calc_bodyIntercept_Unstruc(iGP, jGP, kGP, xGP, yGP, zGP, xBI, yBI, zBI , closestElementGP)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE gcm_arrays
    USE grid_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... parameters

    INTEGER,           INTENT(IN)  :: iGP, jGP, kGP
    REAL(KIND=CGREAL), INTENT(IN)  :: xGP, yGP, zGP
    INTEGER,           INTENT(OUT) :: closestElementGP
    REAL(KIND=CGREAL), INTENT(OUT) :: xBI, yBI, zBI
    
!... loop variables

    INTEGER :: iEdge,m,n,nc

!... local variables
 
    INTEGER, PARAMETER       :: NSIZE = 1000
    INTEGER, PARAMETER       :: MSIZE = 20
    INTEGER, DIMENSION(:),ALLOCATABLE   :: NeighElemInd

    REAL(KIND=CGREAL), DIMENSION(:),ALLOCATABLE   :: distMarker

    INTEGER                  :: nCheck
    INTEGER,DIMENSION(1)     :: iDummy(1)
    INTEGER                  :: iBody,numNeighElement
    INTEGER                  :: elemInd,node1,node2,node3,nMarker,iErr
    INTEGER                  :: nEdges,nodeCV,nodeA,nodeB
    INTEGER                  :: cElementGP(1:MSIZE),closestNodeGP(1:MSIZE)
    INTEGER                  :: shortestProbe,cMarker

    REAL(KIND=CGREAL)        :: cNormal,dMin,dsIntercept,xM,yM,zM
    REAL(KIND=CGREAL)        :: area123,areaDiff,distBIElem,distBIElemMin
    REAL(KIND=CGREAL)        :: epsiArea,distGPEI,distGPEIMin
    REAL(KIND=CGREAL)        :: xBITemp, yBITemp, zBITemp
    REAL(KIND=CGREAL)        :: xCV, yCV, zCV
    REAL(KIND=CGREAL)        :: xEI, yEI, zEI
    REAL(KIND=CGREAL)        :: xEITemp, yEITemp, zEITemp
    REAL(KIND=CGREAL)        :: vec01x, vec01y, vec01z
    REAL(KIND=CGREAL)        :: vec12x, vec12y, vec12z
    REAL(KIND=CGREAL)        :: magnitude12,magnitude12Inv,projectedLength
    REAL(KIND=CGREAL)        :: xBIT(1:MSIZE),yBIT(1:MSIZE),zBIT(1:MSIZE),dist(1:MSIZE)
    REAL(KIND=CGREAL)        :: distInside

!******************************************************************************


    iBody   = 1
    nMarker = nPtsBodyMarker(iBody)


! NCheck:  Number of closesest nodes to check
! high values of this variable increases robustness of procedure 
! and also CPU time for finding body intercept.

    nCheck = 3 

    IF (nCheck > MSIZE) THEN
       PRINT*,'nCheck in GCM_calc_bodyIntercept_Unstruc is limited to', MSIZE
       PRINT*,'Increase array size'
       STOP
    ENDIF


! ============================================================================
!   Allocate local array
! ============================================================================

    ALLOCATE(distMarker(nMarker),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &  
       'search_vertex_dotNorm: Memory Allocation Error for distMarker'
      STOP
    ENDIF ! ierr 
    ALLOCATE(NeighElemInd(NSIZE),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &  
       'search_vertex_dotNorm: Memory Allocation Error for NeighElemInd'
      STOP
    ENDIF ! ierr 

! ============================================================================
! Get closestMarker for generic point 
! ============================================================================

    dMin = 1.0E+5_CGREAL

    DO m = 1, nMarker
 
        xM = xBodyMarker(iBody,m)
        yM = yBodyMarker(iBody,m)
        zM = zBodyMarker(iBody,m)
          
        distMarker(m) = (xM-xGP)**2 + (yM-yGP)**2 + (zM-zGP)**2
	  
    ENDDO ! m 

    DO nc = 1,NCheck
	iDummy                        = MINLOC(distmarker(1:nMarker))
	closestNodeGP(nc)             = iDummy(1)
        distmarker(closestNodeGP(nc)) = 1.0E20_CGREAL
    ENDDO 

 !print*,closestNodeGP

! ============================================================================
! Find elements that share closest node/marker
! ============================================================================

    DO nc = 1,nCheck


    numNeighElement = 0
    DO m=1,totNumTriElem(iBody)
       IF ( triElemNeig(iBody,1,m) == closestNodeGP(nc) .OR. &
            triElemNeig(iBody,2,m) == closestNodeGP(nc) .OR. &
            triElemNeig(iBody,3,m) == closestNodeGP(nc) ) THEN
          numNeighElement               = numNeighElement + 1
          NeighElemInd(numNeighElement) = m
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if (closestnodeGP == 4312) then
!print*,m
! endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ENDIF
    ENDDO

! ============================================================================
!   Trap error if array NeighElemenInd overflows
! ============================================================================

    IF ( numNeighElement > NSIZE ) THEN
      WRITE(STDOUT,*) &
       'GCM_calc_bodyIntercept_Unstruc: Memory Overflow Error for NeighElemInd'
      WRITE(STDOUT,*) ' Allocated size = ',NSIZE
      WRITE(STDOUT,*) ' Current size   = ',numNeighElement
      WRITE(STDOUT,*) ' Aborting Run'

      STOP
    ENDIF ! NeighElemInd

! ============================================================================
!   Determine which element contains normal intercept
! ============================================================================

    closestElementGP = 0

    epsiArea = 1.0E-04_CGREAL
    distBIElemMin = 1.0E16_CGREAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if (closestnodeGP == 4312) then
!print*,numNeighElement,(NeighElemInd(n),n=1,numNeighElement)
! endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO n = 1,numNeighElement

     elemInd = NeighElemInd(n)

     node1   = triElemNeig(iBody,1,elemInd)
     node2   = triElemNeig(iBody,2,elemInd)
     node3   = triElemNeig(iBody,3,elemInd)
      
! ******************************************************************************
!    Check if BITemp is located inside triangle of surface element
!     through area difference
! ******************************************************************************
    
     CALL check_BIInsideTriangle(iBody,elemInd,node1,node2,node3,xGP,yGP,zGP,&
                                 xBITemp,yBITemp,zBITemp,area123,areaDiff,distInside)
    
!--------------------------------------------------------------------------
!    Select closest Elem and BI coordinates:
!     If BI falls inside the element use that
!     Else Base the selection on the minimum distance
!       between BI and either the norm to closest side or vertices of side
!--------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if (closestnodeGP == 4312) then
!print*,ABS(areaDiff),epsiArea*area123,elemInd
! endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     IF ( ABS(areaDiff) < epsiArea*area123) THEN
        xBI = xBITemp
        yBI = yBITemp
        zBI = zBITemp
        closestElementGP = elemInd
        GOTO 999
     ELSE 

        CALL calc_BIOutsideTriangle(iBody,elemInd,node1,node2,node3,xGP,yGP,zGP, closestNodeGP(nc), &
                                    xBITemp,yBITemp,zBITemp,distBIElem)
       
        IF (distBIElem <= distBIElemMin) THEN
          distBIElemMin = distBIElem
          closestElementGP = elemInd
          xBI = xBITemp
          yBI = yBITemp
          zBI = zBITemp
        ENDIF ! distBIElem
     ENDIF ! areaDiff

!DEBUG
!    IF (iGP == 126 .and. jGP == 23 .and. kGP == 12) then
!       WRITE(355,*)'ZONE'
!       WRITE(355,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!       WRITE(355,*)xBodyMarker(iBody,node2),yBodyMarker(iBody,node2),zBodyMarker(iBody,node2)
!       WRITE(355,*)xBodyMarker(iBody,node3),yBodyMarker(iBody,node3),zBodyMarker(iBody,node3)
!       WRITE(355,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!       WRITE(355,*)'ZONE'
!       WRITE(355,*)xGP,yGP,zGP
!       WRITE(355,*)xBItemp,yBItemp,zBItemp
!       WRITE(356,*)areadiff,closestElementGP
!    ENDIF
!DEBUG

    ENDDO ! n

! ============================================================================
!   Compute coordinates of Body Intercept in a robust manner
!    for the case where the temporary BI is located outside 
!    all the surface elements
!    1. Load coordinates of closest vertex (CV)
! ============================================================================

    xCV = xBodyMarker(iBody,closestNodeGP(nc))
    yCV = yBodyMarker(iBody,closestNodeGP(nc))
    zCV = zBodyMarker(iBody,closestNodeGP(nc))

    distGPEIMin = 1.0E+16_CGREAL

! ============================================================================
!    2. Detemine the indices of the 2 vertices connected to CV
! ============================================================================

    node1   = triElemNeig(iBody,1,closestElementGP)
    node2   = triElemNeig(iBody,2,closestElementGP)
    node3   = triElemNeig(iBody,3,closestElementGP)

    IF ( node1 == closestNodeGP(nc) ) THEN
      nodeCV = node1
      nodeA  = node2
      nodeB  = node3
    ELSEIF ( node2 == closestNodeGP(nc) ) THEN
      nodeCV = node2
      nodeA  = node3
      nodeB  = node1
    ELSEIF ( node3 == closestNodeGP(nc) ) THEN
      nodeCV = node3
      nodeA  = node1
      nodeB  = node2
    END IF ! node1 

! ============================================================================
!   3. Compute edge01 (CV-->GP), edge12 (CV-->A), edge13 (CV-->B) vectors
!      Project vector GP-CV onto CV-A or CV-B to find temporary edge intercept
!      If the projectedLength is < 0 or > edgeLength, EI is outside
!      Else EI is inside then compute its Location and distance to GP
! ============================================================================

     vec01x = xGP - xCV
     vec01y = yGP - yCV
     vec01z = zGP - zCV

     nEdges = 2

     DO iEdge = 1, nEdges
      SELECT CASE(iEdge)
        CASE(1)
         node1 = nodeA
        CASE(2)
         node1 = nodeB
      END SELECT ! iEdge

      vec12x = xBodyMarker(iBody,node1) - xCV
      vec12y = yBodyMarker(iBody,node1) - yCV
      vec12z = zBodyMarker(iBody,node1) - zCV

      magnitude12 = SQRT(vec12x**2 + vec12y**2 + vec12z**2)

      magnitude12Inv = oned/magnitude12
     
      vec12x = vec12x*magnitude12Inv
      vec12y = vec12y*magnitude12Inv
      vec12z = vec12z*magnitude12Inv

      projectedLength = vec01x*vec12x +vec01y*vec12y +vec01z*vec12z 
 
!----------------------------------------------------------------------------
!     Edge-Intercept (EI) point is outside Edge if pL < 0 or pL>magnitude12
!      else EI is inside and compute coordinates and distance
!      No need to take SQRT for distGPEI to save computations
!      Load EI into temporary value
!----------------------------------------------------------------------------

      IF ( projectedLength < zero .OR. &
           projectedLength > magnitude12     ) THEN
        xEITemp = xCV
        yEITemp = yCV
        zEITemp = zCV
      ELSE
       xEITemp = xCV + projectedLength*vec12x
       yEITemp = yCV + projectedLength*vec12y
       zEITemp = zCV + projectedLength*vec12z
      END IF ! projectedLength

!----------------------------------------------------------------------------
!     Find mininum value of |GP-EI| and corresponding EI
!----------------------------------------------------------------------------

      distGPEI = (xEITemp-xGP)**2.0_CGREAL &
               + (yEITemp-yGP)**2.0_CGREAL &
               + (zEITemp-zGP)**2.0_CGREAL 
  
      IF ( distGPEI < distGPEIMin ) THEN
        distGPEIMin = distGPEI
        xEI = xEITemp
        yEI = yEITemp
        zEI = zEITemp
      END IF ! distGPEI
     END DO ! iEdge

     xBI = xEI 
     yBI = yEI
     zBI = zEI

!DEBUG
!    IF (iGP == 126 .and. jGP == 23 .and. kGP == 12) then
!       WRITE(355,*)'ZONE'
!       WRITE(355,*)xGP,yGP,zGP
!       WRITE(355,*)xBI,yBI,zBI
!    ENDIF
!DEBUG

! TEMPORARY
!    WRITE(365,*)xBI,yBI,zBI
! END TEMPORARY

999 CONTINUE

  xBIT(nc) = xBI
  yBIT(nc) = yBI
  zBIT(nc) = zBI
  cElementGP(nc) = closestElementGP

  ENDDO  ! nc

  DO nc = 1,nCheck
     dist(nc) =  (xBIT(nc)-xGP)**2.0_CGREAL &
               + (yBIT(nc)-yGP)**2.0_CGREAL &
               + (zBIT(nc)-zGP)**2.0_CGREAL
  ENDDO

  iDummy           = MINLOC(dist(1:nCheck))
  shortestProbe    = iDummy(1)
  xBI              = xBIT(shortestProbe)
  yBI              = yBIT(shortestProbe)
  zBI              = zBIT(shortestProbe)
  closestElementGP = cElementGP(shortestProbe)

  DEALLOCATE(NeighElemInd)
  DEALLOCATE(distMarker)

  END SUBROUTINE hybrid_calc_bodyIntercept_Unstruc