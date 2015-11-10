 !-----------------------------------------------------------------------------
SUBROUTINE time_step_viscous

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE stat_arrays
    USE stat_vort_arrays
    USE fea_unstructure_surface
    USE body_dynamics
    USE nlold_arrays
    USE iso_c_binding
#ifdef PETSC
    use petsc_mod
#endif
    USE turb_parameters

    IMPLICIT NONE

    !... Local variables
    INTEGER              :: iBody,i,j,k,n
    REAL(KIND=CGREAL)	  :: sum
    REAL(KIND=CGREAL)     :: starttime, endtime
	REAL(KIND=CGREAL)     :: timer(7) ! 1: AD Eqn; 2: Poission Eqn; 3: Moving boundary; 4: FSI; 5: output; 6: Total
	REAL(KIND=CGREAL)     :: cumulate_time(7) ! 1: AD Eqn; 2: Poission Eqn;
    INTEGER              :: clock1, clock2, clock_rate
    INTEGER              :: clock01, clock02, clock_rate0
    INTEGER              :: clock11, clock12, clock_rate1
    INTEGER              :: Count_Converged

!******************************************************************************
!  initiate bodyMarkerVel(3*nPtsBodyMarker)
    DO i=1,nPtsBodyMarker(1)
        bodyMarkerVel(3*i)=0.0
        bodyMarkerVel(3*i-1)=0.0
        bodyMarkerVel(3*i-2)=0.0
    END DO

! Start time stepping loop
    cumulate_time=zero
	timer=zero
! Given ntime_start = 0
    DO ntime = ntime_start+1,ntime_start+no_tsteps        ! move solution from n --> n+1
      call system_clock(clock01)

! This is the time after this step is done, knowing the initial time is from 0. ntime = 1 corresponds to time (0+dt)
      time = time + dt
      Converged_FSI(:) = .FALSE.
      Count_Converged = 0

      IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) THEN
        WRITE(*,*)'======================================================='
        WRITE(*,'(A,I6,A,F15.7,A,I6)') 'NTIME = ',ntime,',  Time = ',time,',  NTIME TO GO = ',no_tsteps+ntime_start-ntime
      ENDIF

      IF (FSI_on) THEN
         FSI_CONVERGE = .False.
         FSI_ITERATION = 0
         yBodyMarkerOld(:,:) = yBodyMarker(:,:)
         struc_olddisp(:,:) = struc_disp(:,:)
         print *, 'coremax/min struc_disp(i,6) =', maxval(struc_disp(:,6)),minval(struc_disp(:,6))
      ENDIF

      IF (boundary_motion_type(1) >= FEA_FLOW_STRUC_INTERACTION) THEN
         niterFS = 0
      ENDIF

!------------------------------------------------------------------------------
!     Compute turbulent viscosity
!     Compute AD coefficients for non-moving boundaries
!     Set BC for viscosity
!     Add molecular viscosity
!     Note: nuTot is computed at n-time step level
!------------------------------------------------------------------------------

      IF ( turbActive == ACTIVE ) THEN
        IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) WRITE(*,*) 'Entering TURB_CalcVisc '
        CALL TURB_Calcvisc()

        IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) WRITE(*,*) 'Entering TURB_Visc_set_bc '
        CALL TURB_Visc_set_bc()

        DO k =0, nz+1
        DO j =0, ny+1
        DO i =0, nx+1
          viscTot(i,j,k) = viscTot(i,j,k) +reInv
          bcxvisc(i,j,k) = bcxvisc(i,j,k) +reInv
          bcyvisc(i,j,k) = bcyvisc(i,j,k) +reInv
          bczvisc(i,j,k)= bczvisc(i,j,k) +reInv
        END DO ! i
        END DO ! j
        END DO ! k

        IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) WRITE(*,*) 'Entering TURB_Visc_SetBoundCells '
        CALL TURB_Visc_SetBoundCells()

		IF ( boundary_motion /= MOVING_BOUNDARY ) THEN
          IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) WRITE(*,*) 'Entering set_solve_ad()'
          CALL set_solve_ad()
		END IF ! boundary_motion
      END IF   ! turbActive

!------------------------------------------------------------------------------
!     Compute advection-diffusion terms
!     Note: NL should be at n-time step level, hence they are computed
!           BEFORE moving the body boundary
!------------------------------------------------------------------------------
      CALL rhs_advec_diff()                              ! compute NLU   and VIS

!------------------------------------------------------------------------------
!     Move Body Boundary and compute coefficients
!------------------------------------------------------------------------------

      IF ( boundary_motion == MOVING_BOUNDARY .AND. ntime >= 1) THEN
		call system_clock(clock1)

        PRINT*,'CALL move_boundary()'
        if(boundary_motion_type(1) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
            boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING )then
            call move_boundary_defor() ! Added by G. Liu
        else
           CALL move_boundary() ! move boundary to (n+1)
        endif
        PRINT*,'CALL set_solve_ad()'
        CALL set_solve_ad()


        IF ( turbActive == ACTIVE ) THEN
          WRITE(*,*) 'Entering TURB_Visc_SetBoundCells '
          CALL TURB_Visc_SetBoundCells()
        END IF ! turbActive

        SELECT CASE (pp_solver_type)
		CASE (PP_SOLVER_TYPE_PETSC)
#ifdef PETSC
		  CALL petsc_setup_coeff_solver()
#else
		  PRINT*,'USER ERROR: PETSC Solver Active in Input Deck '
		  PRINT*,'  Code not compiled with PETSC Flag'
		  PRINT*,'  Code will stop and exit '
		  PRINT*,'  Either set solver to MG or compile with PETSC=1'
		  STOP
#endif

		END SELECT ! it_solver

		call system_clock(clock2, clock_rate)
		timer(3)=DBLE(clock2-clock1)/DBLE(clock_rate)
		IF ( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) &
		  WRITE(*,*) '***** Time(sec) for moving boundary:', timer(3)

      END IF ! boundary_motion

      DO iBody = 1,nBody
        IF (WALL_TYPE(iBody) == POROUS_OR_SLIP)THEN
          CALL wall_velocity(iBody)
        ENDIF
      END DO

      IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) THEN
        WRITE(*,*) 'AREAX1; AREAX2 = ',areax1,areax2
        WRITE(*,*) 'AREAY1; AREAY2 = ',areay1,areay2
        WRITE(*,*) 'AREAZ1; AREAZ2 = ',areaz1,areaz2
      END IF

!     WRITE(*,*) 'Entering FreshCell_CalcExpWeight '
      CALL FreshCell_CalcExpWeight()
      WRITE(*,*) '*********vega_interpolateindexratio_c '
      CALL vega_interpolateindexratio_c(markerInterpolateIndex,markerInterpolateRatio) ! Added by CJ Yuan July.14.2015
!------------------------------------------------------------------------------
!     Set viscosity for fresh cells
!------------------------------------------------------------------------------
!
!      IF ( turbActive == ACTIVE ) THEN
!        IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) WRITE(*,*) 'Entering TURB_Visc_SetFreshCell '
!        CALL TURB_Visc_SetFreshCell()
!      END IF ! turbActive

!------------------------------------------------------------------------------
!     Compute weights and update RHS for Fresh cells
!------------------------------------------------------------------------------
      CALL FreshCell_UpdateRhs()

      IF(boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT .or. boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
         boundary_motion_type(1) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
         boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING )THEN ! Added by Geng
        DO n=1,nbody
          xcent_prev(n)=xcent(n)
          ycent_prev(n)=ycent(n)
          zcent_prev(n)=zcent(n)
        ENDDO
        quat_prev=quat_iter
      ENDIF ! DYNAMICS_COUPLED_QUAT, DYNAMICS_COUPLED_MofI_QUAT, DYNAMICS_COUPLED_SWIMMING

      !call write_dump_debug2d_3vars('qPre',niterFS,U,V,p) ! liu debug
101   CALL set_bc()    ! fill boundary arrays with (u)^n+1
      IF (FSI_on .or. boundary_motion_type(1) >= FEA_FLOW_STRUC_INTERACTION) CALL face_vel()  !Added by Wanh Modified by CJ Yuan
      CALL rhs_adjust_bc()
!------------------------------------------------------------------------------
!     Adjust RHS for 2D computations
!------------------------------------------------------------------------------
      IF( nDim == DIM_2D) CALL rhs_adjust2D
!------------------------------------------------------------------------------
!     Solve intermediate velocity field using Thomas algorithm
!------------------------------------------------------------------------------
      call system_clock(clock1)
      WRITE(*,*) 'Entering SOLVE_AD '
      if(ntime==2 .and. niterFS == 0)then
          write(5678,*)  'VARIABLES="X","Y","iblank","iup","ium","jup","jum","kup","kum"'
          write(5678,*)  'ZONE F=POINT, I=',nxc,', J=',nyc
          k=1
          do j=1,nyc
          do i=1,nxc
            write(5678,5678) xc(i),yc(j),iblank(i,j,k),iup(i,j,k),ium(i,j,k),jup(i,j,k),jum(i,j,k),kup(i,j,k),kum(i,j,k)
          end do
          end do
      end if

      if(pressure_osc_velocity.or.pressure_osc_pressure)then
        call identify_hybrid_cell()
      end if

      IF(.not.dryRun)THEN
        if(pressure_osc_velocity)THEN
          call solve_ad()
          call hybrid_cell_interpolation
        ELSE
          call solve_ad()
        END IF
      END IF

      call system_clock(clock2, clock_rate)
      timer(1)=DBLE(clock2-clock1)/DBLE(clock_rate)

      IF ( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) &
        WRITE(*,*) '***** Time(sec) for solving advection-diffusion eq:', timer(1)
!------------------------------------------------------------------------------
!     Compute face velocities
!------------------------------------------------------------------------------
      WRITE(*,*) 'Entering FACE_VEL '
      CALL face_vel()
      if(pressure_osc_velocity)then
        call faceVelCorret()
        call update_hybrid_velocity
        CALL face_vel()
        call faceVelCorret()
      end if
      IF (cure_pressure_oscillations) CALL correct_face_vel()
!------------------------------------------------------------------------------
!     Compute RHS for the Poisson Pressure Equation (PPE)
!------------------------------------------------------------------------------
     WRITE(*,*) 'Entering RHS_POISSON '
 105     CALL rhs_poisson(sum)
      IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) THEN
         PRINT*, 'Sum of Poisson RHS = ',sum
      END IF
!------------------------------------------------------------------------------
!    Solve the Poisson Pressure Equation (PPE)
!------------------------------------------------------------------------------

      call system_clock(clock1)
      WRITE(*,*) 'Entering SOLVE_POISSON '
      if(.not.dryRun) CALL solve_poisson()
      call system_clock(clock2, clock_rate)
      timer(2)=DBLE(clock2-clock1)/DBLE(clock_rate)
      if(pressure_osc_velocity)then
          u=u_bak
          v=v_bak
          w=w_bak
      end if
      IF ( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) &
        WRITE(*,*) '***** Time(sec) for solving poisson eq:', timer(2)

!------------------------------------------------------------------------------
!    Correct velocity field and update pressure
!------------------------------------------------------------------------------
      WRITE(*,*) 'Entering CORRECT_VEL '
      if(.not.dryRun) CALL correct_vel()
      WRITE(*,*) 'Entering UPDATE_PRESSURE '
      if(.not.dryRun) CALL update_pressure()
      IF( nDim == DIM_2D) CALL vel_adjust2D
!------------------------------------------------------------------------------
!     Solve for structure, added by Wanh 05/05/10
!------------------------------------------------------------------------------
!      call system_clock(clock1)

      IF (FSI_on) THEN
         p = pPrime
         DO iBody = 1,Nbody
            WRITE(*,*) 'Body',iBody,'entering fea_getpressure'
            CALL fea_getpressure
            CALL fea_structure(iBody)
            IF (TRNS) THEN
               CALL fea_converge(iBody)
            ELSE
               CALL fea_converge_static(iBody)
            ENDIF
         ENDDO
!        CALL write_minmaxvals()
         IF (NTIME == 1) FSI_CONVERGE = .TRUE.
         IF (.NOT. FSI_CONVERGE) GOTO 101
         print *, 'core, after fsi, max/min struc_disp(i,6) =', maxval(struc_disp(:,6)),minval(struc_disp(:,6))
      ENDIF   !End if FSI_on

      IF (NTIME>0 .and. boundary_motion_type(1) >= FEA_FLOW_STRUC_INTERACTION) THEN
         p = pPrime
         DO iBody = 1,nBody
            IF (boundary_motion_type(1) == FEA_FLOW_STRUC_INTERACTION) THEN ! add CJ Yuan 05/26/2015
               IF (nBody_solid>0) THEN
                  CALL drag_lift_solid()
               ELSE IF (nBody_membrane>0) THEN
                  CALL drag_lift_membrane()
               ENDIF
            ELSE IF (boundary_motion_type(1) == DYNAMICS_COUPLED) THEN
               IF (nBody_solid>0) THEN
                  CALL MofI_CofG(iBody)
				  call out_MoI_CoM_debug
                  CALL drag_lift_solid()
               ELSE IF (nBody_membrane>0) THEN
                  CALL MofI_CofG_Mem_full_FBI(iBody)
                  CALL drag_lift_membrane()
               ENDIF

            ELSE IF (boundary_motion_type(1) == PARTIAL_DYNAMICS_COUPLED) THEN
               IF (nBody_membrane>0) THEN
                  CALL MofI_CofG_Mem(iBody)
                  CALL drag_lift_membrane_partial_dynamics(iBody)
               ENDIF

            ELSE if(boundary_motion_type(1) == BIO_DYNAMICS_COUPLED)then
                if(ibody==1.and.nBody_solid>0)then
                    call MofI_CofG(1)
                    call drag_lift_solid()
                end if
                if(ibody==1.and.nBody_membrane>0)then
                    call drag_lift_membrane()
                end if

            ELSE IF(boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT)THEN
                IF( nbody_solid>0 )then
                    call MofI_CofG_quat(iBody)     ! Get the position of the center of mass
                    call out_MoI_CoM_debug
                    call drag_lift_solid()
                    call moment_rot(iBody)       ! transform the moments from the inertia frame to the non-inertia frame
                else if( nbody_membrane>0)then
                    call MofI_CofG_quat(iBody)
                    call drag_lift_membrane()
                    call moment_rot(iBody)
                endif

            ELSE IF(boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT)THEN  ! Added by G. Liu
                IF( nbody_solid>0 )then
                    call MofI_CofG_transform(iBody)  ! Get the position of the center of mass, Moment of inertia transformation
                    call out_MoI_CoM_debug           ! debug
                    call drag_lift_solid()
!                    call moment_rot(iBody)       ! transform the moments from the inertia frame to the non-inertia frame
                else if( nbody_membrane>0)then
                    call MofI_CofG_transform(iBody)
                    call drag_lift_membrane()
                endif

            ELSE IF(boundary_motion_type(ibody) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
                    boundary_motion_type(iBody) == DYNAMICS_COUPLED_SWIMMING )THEN  ! Added by G. Liu
                IF( nbody_solid>0 )then
                    call MofI_CofG(iBody)            ! can only deal with one body
                    call out_MoI_CoM_debug           ! debug
                    call drag_lift_solid()
                else if( nbody_membrane>0)then
                    call MofI_CofG_Mem_full_FBI(iBody)
                    call drag_lift_membrane()
                endif
            ENDIF

            IF(boundary_motion_type(1)==BIO_DYNAMICS_COUPLED)THEN
                IF(ibody==1)THEN
                    call dynamics_motion(1)
                END IF
            ELSE IF(boundary_motion_type(1) == FEA_FLOW_STRUC_INTERACTION)THEN !Added by CJ Yuan
                print *, '****** Vega_deforamtion_c'
                call vega_deformation_c(markerPressure,markerInterpolateVelocity,bodyMarkerVel) !Added by CJ Yuan
            ELSE
                CALL dynamics_motion(iBody)
            END IF

            IF (boundary_motion_type(1) == DYNAMICS_COUPLED .or. boundary_motion_type(1) == BIO_DYNAMICS_COUPLED .OR. &
                boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT .or. boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT) THEN
               CALL compute_marker_vel(iBody) ! find the velocity here by CJ Yuan
            ELSE IF (boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
                     boundary_motion_type(iBody) == DYNAMICS_COUPLED_SWIMMING              )THEN ! Added by G. Liu
               CALL compute_marker_vel(iBody)
               call Add_defor_vel(iBody)

            ELSE IF (boundary_motion_type(1) == PARTIAL_DYNAMICS_COUPLED) THEN
               CALL compute_marker_vel_section(iBody)
            ENDIF

            if(boundary_motion_type(iBody)==BIO_FOLLOWED_DYNAMICS_COUPLED)then
                write(*,*) 'Skip dynamics convergence check for Body #',iBody
                Count_Converged=Count_Converged+1
            else if(boundary_motion_type(iBody)==FEA_FLOW_STRUC_INTERACTION) then
                CALL vega_markerVel_convergenceCheck  ! Added By CJ Yuan
            else
                CALL FSConvergeCheck(iBody)
            end if

            IF (Converged_FSI(iBody)) THEN
               WRITE (*,*) 'Body',iBody,' is converged.'
               Count_Converged = Count_Converged+1
               print *, 'Count_Converged=', Count_Converged
            ENDIF

         ENDDO

         IF (abs(Count_Converged - nBody) < 0.1) THEN
            WRITE (*,*) 'All body converged at ntime =', ntime, 'vega_reNewBodyPosition_c'
            CALL vega_reNewBodyPosition_c()
!            if(boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT .or. boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT)then
            if( boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT )then
              do n=1,nbody   ! quat_modify
                xcent(n) = xcent_prev(n) + dt*vxcent(n)
                ycent(n) = ycent_prev(n) + dt*vycent(n)
                zcent(n) = zcent_prev(n) + dt*vzcent(n)
              enddo
            endif
         ELSE
            niterFS = niterFS + 1
            WRITE (*,*) 'Iteration between F-S: ', niterFs

            IF (niterFS<niterFS_max) THEN
               Count_Converged = 0
               GOTO 101
            ELSE
               WRITE (*,*) 'Fluid-Solid convergence failed.'
!               STOP
            ENDIF
         ENDIF
      ENDIF   !End if (boundary_motion_type)
!------------------------------------------------------------------------------
!     Monitor output
!------------------------------------------------------------------------------
      call system_clock(clock1)

! output the position, velocity and acceleration of the body. Added by G. Liu
    IF (boundary_motion_type(1) == DYNAMICS_COUPLED .or. boundary_motion_type(1) == BIO_DYNAMICS_COUPLED .or. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT .or. boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING  ) THEN
      if(ntime.eq.1)then
        open(2012,file='kinematics.dat')
        write(2012,*)'time   acc_xCG acc_yCG  acc_zCG   vxcent  &
                   vycent   vzcent   xcent  ycent  zcent  angvx  angvy   angvz'
      else
        open(2012,file='kinematics.dat',access='append')
      endif
      write(2012,2011)time,acc_xCG(1),acc_yCG(1),acc_zCG(1), &
                      vxcent(1),vycent(1),vzcent(1), &
                      xcent(1),ycent(1),zcent(1),    &
                      angvx(1),angvy(1),angvz(1)
2011  format(13(1x,e13.6))
      close(2012)
    ENDIF  ! boundary_motion_type(1)

      IF ( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) THEN
        CALL write_monitor()
        IF (turbActive == ACTIVE) CALL TURB_write_monitor()
      ENDIF

      !IF ( nStat > STATS_NONE ) THEN
      IF ( ntime > STATS_SUM) THEN

        !IF ( MOD(ntime,nstat) == 0 .OR. ntime==ntime_start+1 ) THEN
         CALL calc_statistics(0)
         write(*,*) '++++++++++++++++++++++++SUM++++++++++++++++++++++++'
         write(*,*) ntime, no_tsteps
         write(*,*) '++++++++++++++++++++++++SUM++++++++++++++++++++++++'
!         pause
!         CALL calc_statistics_vorticity(0)
         statCtr = statCtr + 1
         statCtrv = statCtrv + 1
        !ENDIF
      END IF ! nStat

      IF ( nmonitor_probe_liftdrag > STATS_NONE ) THEN
        IF ( MOD(ntime,nmonitor_probe_liftdrag) == 0 .OR. ntime==ntime_start+1 ) THEN
          CALL write_probe_files()

          SELECT CASE(boundary_formulation)
            CASE (SSM_METHOD)
              IF(nBody_solid > 0) THEN
                CALL drag_lift_solid()
              END IF

              IF(nBody_membrane > 0) THEN
                CALL drag_lift_membrane()
                IF (NTIME>0 .and. boundary_motion_type(1) == PARTIAL_DYNAMICS_COUPLED) THEN
                   DO iBody = 1,nBody
                      CALL drag_lift_membrane_partial_dynamics(iBody)
                   ENDDO
                ENDIF
              END IF

            CASE (GCM_METHOD)
              CALL GCM_drag_lift()

           END SELECT ! boundary_formulation
         ENDIF
      END IF

      IF ( MOD(ntime,ndump)==0)THEN
         CALL write_dump()
      ENDIF

      IF ( MOD(ntime,nrestart) == 0 .OR. ntime==ntime_start+no_tsteps) THEN
        CALL write_restart()
      ENDIF

      IF ( ntime == no_tsteps) THEN
          CALL calc_statistics(1)
          write(*,*) '!!!!!!!!!!!!!!!!!!!!SUM_END!!!!!!!!!!!!!!!!!!!!'
          write(*,*) ntime, no_tsteps
          write(*,*) '!!!!!!!!!!!!!!!!!!!!SUM_END!!!!!!!!!!!!!!!!!!!!'
!          pause
!        IF ( nStat > STATS_NONE ) CALL calc_statistics_vorticity(1)
      ENDIF

      call system_clock(clock2, clock_rate)
      timer(5)=DBLE(clock2-clock1)/DBLE(clock_rate)

      call system_clock(clock02, clock_rate0)
      timer(6)=DBLE(clock02-clock01)/DBLE(clock_rate0)

      cumulate_time=cumulate_time+timer

      IF ( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) THEN
        WRITE(*,*) '***** Time(sec) for output:', timer(5)
        WRITE(*,*) '***** Time(sec) for current time step:', timer(6)
        WRITE(*,*)
        write(*,*) 'Cumulate time(sec) for moving boundary:', cumulate_time(3)
        write(*,*) 'Cumulate time(sec) for ad-diff eq:', cumulate_time(1)
        write(*,*) 'Cumulate time(sec) for poisson eq:', cumulate_time(2)
        write(*,*) 'Cumulate time(sec) for output:', cumulate_time(5)
        write(*,*) 'Cumulate time(sec) till current time step:', cumulate_time(6)
        write(1003,*) 'No.',ntime,'Cumulate time(sec) for moving boundary:', cumulate_time(3)
        write(1003,*) 'No.',ntime,'Cumulate time(sec) for ad-diff eq:', cumulate_time(1)
        write(1003,*) 'No.',ntime,'Cumulate time(sec) for poisson eq:', cumulate_time(2)
        write(1003,*) 'No.',ntime,'Cumulate time(sec) for output:', cumulate_time(5)
        write(1003,*) 'No.',ntime,'Cumulate time(sec) till current time step:', cumulate_time(6)
        write(1004,'(6(2X,E19.11))') time, timer(1), timer(2), timer(3), timer(5), timer(6)
        write(1005,'(6(2X,E19.11))') time, cumulate_time(1),cumulate_time(2),cumulate_time(3),cumulate_time(5),cumulate_time(6)
      ENDIF

      if(optimization.and.optStop)then
        write(*,*) 'Opt: Number of evaluations are enough!'
        call write_restart()
        stop
      end if
    ENDDO ! ntime
       5678 format(2f16.8,7i2)
   END SUBROUTINE time_step_viscous

SUBROUTINE interpolate_pressure_velocity(markerInterpolateIndex,markerInterpolateRatio,&
                            markerInterpolatePressure,markerInterpolateVelocity)

END SUBROUTINE interpolate_pressure_velocity
