!---------------------------------------
!  SUBROUTINE open_files()
!  SUBROUTINE read_inputs()
!  SUBROUTINE check_inputs()
!  SUBROUTINE read_inputs_bound()
!  SUBROUTINE check_inputs_bound()
!  SUBROUTINE init_simulation()
!  SUBROUTINE flow_stop()
!  SUBROUTINE init_flow()
!  SUBROUTINE read_restart_flow()
!  SUBROUTINE read_restart_body()
!  SUBROUTINE write_restart()
!  SUBROUTINE set_arrays_zero()
!  SUBROUTINE initialize_flowfield()
!  SUBROUTINE abort_vicar3d(abort_code)
!---------------------------------------
!
!-------------------------------------------------
! Code to simulate 3-D incompressible, unsteady, viscous
! flow in domain with multiple, complex embedded obstacles
!
! Only Cartesian grids (uniform or non-uniform) are allowed
!
! Obstacles boundaries are of a stairstep nature and obstacles
! are defined by a "iblank" variable which is 1 inside the obstacle
! 0 in the fluid.
!
! spatial discretization -  2nd order central difference
! temporal discretization-  2-step time-split scheme
!                           AB2 for non-linear term
!                           CN for for viscous terms
!
!     |-------|-------|-------|-------|-------
!     |       |       |       |       |
!     |   o   |   o   |   o   |   o   |   o
!     |       |       |       |       |
!     |-------|-------|-------|-------|-------
!     |       |*******|*******|*******|
!     |   o   +***1***|***1***|***1***|   o
!     |       |*******|*******|*******|
!     |-------|-------|-------|-------|-------
!     |       |*******|*******|       |
!     |   o   |***1***|***1***|   o   |   o
!     |       |*******|*******|       |
!     |-------|-------|-------|-------|-------
!     |       |*******|*******|       |
!     |   o   |***1***|***1***|   o   |   o
!     |       |*******|*******|       |
!     |-------|-------|-------|-------|-------
!     |       |       |       |       |
!     |   o   |   o   |   o   |   o   |   o
!     |       |       |       |       |
!     |-------|---+---|-------|-------|-------
!

!---------------------------------------------
! Main Program
!---------------------------------------------

PROGRAM VICAR3D
  USE global_parameters
  USE flow_parameters
  USE fea_unstructure_surface

!  USE flow_arrays
!  USE pressure_arrays
!  USE boundary_arrays
!  USE grid_arrays
!  USE multiuse_arrays
!  USE stat_arrays
!  USE stat_vort_arrays
!  USE body_dynamics
!  USE nlold_arrays

   USE iso_c_binding
#ifdef PETSC
   USE PETSC_MOD
#endif
   IMPLICIT NONE

   REAL(c_double),DIMENSION(3*numVertices):: bodyMarkerForce,bodyMarkerVel  ! Added by CJ Yuan July.17.2015
   REAL(c_double),DIMENSION(6*numVertices):: markerInterpolateRatio,markerInterpolateVelocity ! Added by CJ Yuan July.17.2015
   REAL(c_double),DIMENSION(2*numVertices):: markerInterpolatePressure ! Added by CJ Yuan July.17.2015
   INTEGER(c_int),DIMENSION(6*numVertices):: markerInterpolateIndex ! Added by CJ Yuan July.17.2015

! variables for time_step
    INTEGER              :: iBody,i,j,k,n
    REAL(KIND=CGREAL)	  :: sum
    REAL(KIND=CGREAL)     :: starttime, endtime
	REAL(KIND=CGREAL)     :: timer(7) ! 1: AD Eqn; 2: Poission Eqn; 3: Moving boundary; 4: FSI; 5: output; 6: Total
	REAL(KIND=CGREAL)     :: cumulate_time(7) ! 1: AD Eqn; 2: Poission Eqn;
    INTEGER              :: clock1, clock2, clock_rate
    INTEGER              :: clock01, clock02, clock_rate0
    INTEGER              :: clock11, clock12, clock_rate1
    INTEGER              :: Count_Converged

#ifdef PETSC
    PRINT*,'ENTERING PETSC_INIT'
    CALL petsc_init()
#endif

    PRINT*,'ENTERING OPEN_FILES'
    CALL open_files()

    PRINT*,'ENTERING READ_INPUTS'
    CALL read_inputs()

    PRINT*,'ENTERING CHECK_INPUTS'
    CALL check_inputs()

    PRINT*,'ENTERING READ_INPUTS_BOUND'
    CALL read_inputs_bound()

    PRINT*,'ENTERING CHECK_INPUTS_BOUND'
    CALL check_inputs_bound()

    PRINT*,'ENTERING TURB_CHECK_INPUTS'
    CALL turb_check_inputs()

    PRINT*,'ENTERING OPEN_DRAG_FILES'
    CALL open_drag_files()

    PRINT*,'ENTERING READ_PROBE_INPUT'
    CALL read_probe_inputs()

    PRINT*,'ENTERING INIT_SIMULATION'
    CALL init_simulation()

    IF (boundary_motion_type(1) == FEA_FLOW_STRUC_INTERACTION) THEN !Add By CJ Yuan June.15.2015
        CALL vega_FEM_initiate_c()  ! Modified by CJ Yuan July.13.2015
        CALL vega_interpolateindexratio_c(markerInterpolateIndex,markerInterpolateRatio) ! Added by CJ Yuan July.14.2015
    END IF
    IF(FSI_on) THEN    ! Added by Wanh 05/05/10
      PRINT*, 'ENTER FEA_INITIAL'
      CALL fea_initial()
    END IF

    PRINT*,'ENTERING OPEN_PROBE_FILES'
    CALL open_probe_files()

    PRINT*,'ENTERING TIME_STEP'
    IF (flow_type == VISCOUS_FLOW ) THEN
      !CALL time_step_viscous()
      cumulate_time=zero
	  timer=zero
	  DO ntime = ntime_start+1,ntime_start+no_tsteps
        CALL system_clock(clock01)
        time=time+dt
        Converged_FSI(:)=.FALSE.
        Count_Converged = 0
        IF(boundary_motion_type(1)>=FEA_FLOW_STRUC_INTERACTION)THEN
          niterFS = 0
        ENDIF
        CALL rhs_advec_diff()
        IF(boundary_motion==MOVING_BOUNDARY.AND.ntime>=1)THEN
		  CALL system_clock(clock1)
          PRINT*,'CALL move_boundary()'
          IF(boundary_motion_type(1) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
            boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING )THEN
            CALL move_boundary_defor() ! Added by G. Liu
          ELSE
            CALL move_boundary() ! move boundary to (n+1)
          ENDIF
          PRINT*,'CALL set_solve_ad()'
          CALL set_solve_ad()
          SELECT CASE(pp_solver_type)
            CASE(PP_SOLVER_TYPE_PETSC)
            #ifdef PETSC
		      CALL petsc_setup_coeff_solver()
            #else
              PRINT*,'USER ERROR: PETSC Solver Active in Input Deck '
		      PRINT*,'  Code not compiled with PETSC Flag'
		      PRINT*,'  Code will stop and exit '
              PRINT*,'  Either set solver to MG or compile with PETSC=1'
		      STOP
            #endif
          END SELECT !it_solver
          CALL system_clock(clock2, clock_rate)
		  timer(3)=DBLE(clock2-clock1)/DBLE(clock_rate)
		  IF(MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) &
		  WRITE(*,*) '***** Time(sec) for moving boundary:', timer(3)
        ENDIF !boundary_motion
        DO iBody=1,nBody
          IF(WALL_TYPE(iBody) == POROUS_OR_SLIP)THEN
            CALL wall_velocity(iBody)
          ENDIF
        ENDDO
        IF( MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) THEN
          WRITE(*,*) 'AREAX1; AREAX2 = ',areax1,areax2
          WRITE(*,*) 'AREAY1; AREAY2 = ',areay1,areay2
          WRITE(*,*) 'AREAZ1; AREAZ2 = ',areaz1,areaz2
        ENDIF
        CALL FreshCell_CalcExpWeight()
        WRITE(*,*) '*********vega_interpolateindexratio_c '
        CALL vega_interpolateindexratio_c(markerInterpolateIndex,markerInterpolateRatio) ! Added by CJ Yuan July.14.2015
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
        ENDIF
101     CALL set_bc()    ! fill boundary arrays with (u)^n+1
        IF(FSI_on .or. boundary_motion_type(1)>=FEA_FLOW_STRUC_INTERACTION) CALL face_vel()  !Added by Wanh Modified by CJ Yuan
        CALL rhs_adjust_bc()
        IF(nDim==DIM_2D) CALL rhs_adjust2D
        CALL system_clock(clock1)
        WRITE(*,*) 'Entering SOLVE_AD '
        IF(pressure_osc_velocity.or.pressure_osc_pressure) CALL identify_hybrid_cell()
        IF(.not.dryRun)THEN
          IF(pressure_osc_velocity)THEN
            CALL solve_ad()
            CALL hybrid_cell_interpolation
          ELSE
            CALL solve_ad()
          ENDIF
        ENDIF
        CALL system_clock(clock2, clock_rate)
        timer(1)=DBLE(clock2-clock1)/DBLE(clock_rate)
        IF(MOD(ntime,nmonitor) == 0 .OR. ntime==ntime_start+1 ) &
        WRITE(*,*) '***** Time(sec) for solving advection-diffusion eq:', timer(1)
        WRITE(*,*) 'Entering FACE_VEL '
        CALL face_vel()
        IF(pressure_osc_velocity)THEN
          CALL faceVelCorret()
          CALL update_hybrid_velocity
          CALL face_vel()
          CALL faceVelCorret()
        ENDIF
        IF(cure_pressure_oscillations) CALL correct_face_vel()
        WRITE(*,*) 'Entering RHS_POISSON '
105     CALL rhs_poisson(sum)



    ELSE
      CALL time_step_potential()
    ENDIF

    CALL flow_stop()

#ifdef PETSC
    CALL petsc_exit()
#endif

END PROGRAM VICAR3D

!---------------------------------------------
! File numbers 50 to 100 reserved for input/output files
!---------------------------------------------

   SUBROUTINE open_files()

    USE global_parameters
    USE flow_parameters
    USE turb_parameters

    IMPLICIT NONE

    pi = 4.0_CGREAL*ATAN(oned)

! set File Units numbers

    ifuInput         = 50
    ifuIblankIn      = 52
    ifuRstrtFlowIn   = 53
    ifuRstrtFlowOut  = 54
    ifuRstrtBodyIn   = 55
    ifuRstrtBodyOut  = 56
    ifuBodyIn        = 57
    ifuProbeIn       = 58
    ifuMarkerIn      = 59
    ifuUnstrucSurfIn = 60
    ifuUnstrucSurfOut= 61
    ifuBodyOut       = 62

    ifuDragOut      = 100
    ifuDragOutHinge = 120       !Added by Wanh
    ifuPrsbMomentRef= 130
    ifuFreshCellOut = 150
    ifuMarkerOut    = 160
    ifuStatOut      = 170
    ifuStatPlot     = 171
    ifuProbeOut     = 180
    ifuPowerMembrane2D=200

    ifuRstrtTurbIn  = 63
    ifuRstrtTurbOut = 64

!   Added by Wanh 05/05/10
    ifuFEA_input = 65
    ifuFEA_Bound_Cond = 66
    ifuBodyCGPath = 70          !Added by Wanh for (partially) dynamic coupling
    ifuRstrtDynamicsIn = 75     !Added by Wanh for (partially) dynamic coupling
    ifuRstrtDynamicsOut = 85     !Added by Wanh for (partially) dynamic coupling
    ifuDynaNonIner      = 86    ! Added by Geng for dynamic coupling based on non-inertia frame
    ifuforce_debug      = 230   ! Added by Geng for debug

    ifuDragOutZone  = 220       !Added by Chengyu

! Open files

    OPEN(ifuInput,         FILE='input.dat',             STATUS='OLD')
    OPEN(ifuIblankIn,      FILE='iblank_in.dat')
    OPEN(ifuRstrtFlowIn,   FILE='restart_flow_in.dat'                    ,FORM='UNFORMATTED')
    OPEN(ifuRstrtBodyIn,   FILE='restart_body_in.dat'                    ,FORM='UNFORMATTED')
    OPEN(ifuMarkerIn,      FILE='marker_in.dat',         STATUS='UNKNOWN')
    OPEN(ifuUnstrucSurfIn, FILE='unstruc_surface_in.dat',STATUS='UNKNOWN')
    OPEN(ifuUnstrucSurfOut,FILE='unstruc_surface_out.dat',STATUS='UNKNOWN')
    OPEN(ifuFreshCellOut,  FILE='freshcell_out.dat',     STATUS='UNKNOWN')

!   Added by Wanh 05/05/10
!    OPEN(ifuFEA_Input,FILE='fea_input.dat',     STATUS='UNKNOWN')
!    OPEN(ifuFEA_Bound_Cond,FILE='fea_boundcond.dat',     STATUS='UNKNOWN')

    OPEN(ifuRstrtDynamicsIn,  FILE='restart_dynamics_in.dat')
    open(4030,file='opt_evaluation_history.dat',status='replace')

! Set index for restarting
    indexRstrt = 1

   END SUBROUTINE open_files
!---------------------------------------------

   SUBROUTINE open_drag_files()

    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE
    INTEGER :: ibody, izoneMax

    CHARACTER*9          :: dragfile
    CHARACTER*5          :: powerfile
    CHARACTER*35         :: indragfile
    CHARACTER*35         :: inpowerfile
    character*19         :: fname

    IF (body_type /= CANONICAL) RETURN

!   Open drag and lift file

    print *, 'ifuDragOut,ifuDragOutHinge=', ifuDragOut,ifuDragOutHinge
    DO ibody = 1, nbody
      dragfile = TRIM("drag_lift")
      powerfile = TRIM("power")
      WRITE(indragfile,101) dragfile,ibody
101   FORMAT(a,'_body_',i3.3,'.dat')
      OPEN(UNIT=ifuDragOut+ibody-1,FILE=indragfile,FORM='formatted',ACTION="write")

      WRITE(inpowerfile,103) powerfile,ibody
103   FORMAT(a,'_body_',i3.3,'.dat')
      OPEN(UNIT=ifuPowerMembrane2D+ibody-1,FILE=inpowerfile,FORM='formatted',ACTION="write")

      !Added by Wanh for hinged body
      WRITE(indragfile,102) dragfile,ibody
102   FORMAT(a,'_hingedbody_',i3.3,'.dat')
      OPEN(UNIT=ifuDragOutHinge+ibody-1,FILE=indragfile,FORM='formatted',ACTION="write")

      IF(zoneSeparate) THEN
      !Added by Chengyu for calculating the force/power of different zones on a membrane surface
      DO izoneMax=1, zoneMax(ibody)
      WRITE(indragfile,104) dragfile,ibody, izoneMax
104       FORMAT(a,'_body_',i3.3,'_zone_',i3.3,'.dat')
          OPEN(UNIT=ifuDragOutZone+izoneMax-1,FILE=indragfile,FORM='formatted',ACTION="write")
      ENDDO !izoneMax
      ENDIF
    END DO ! ibody

    if(boundary_motion_type(1) >= DYNAMICS_COUPLED)then  ! Added by Geng for debug
      do ibody=1,nbody
        write(fname,105)'force_debug_',ibody,'.dat'
105     format(A12,I3.3,A4)
        OPEN(UNIT=ifuforce_debug+ibody-1,FILE=fname,FORM='formatted',ACTION="write")
      enddo
    endif

   END SUBROUTINE open_drag_files
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE read_inputs()

    USE global_parameters
    USE flow_parameters
    USE mg_parameters
    USE turb_global_parameters
    USE turb_parameters

    IMPLICIT NONE

    INTEGER :: n

    READ(ifuInput,*)
    READ(ifuInput,*) nread,ntime_skip_check,dryRun,ibkOut
    READ(ifuInput,*)
    READ(ifuInput,*) ndim, flow_type
    READ(ifuInput,*)
    READ(ifuInput,*) nx,ny,nz
    READ(ifuInput,*)
    READ(ifuInput,*) xgrid_unif,xout
    READ(ifuInput,*)
    READ(ifuInput,*) ygrid_unif,yout
    READ(ifuInput,*)
    READ(ifuInput,*) zgrid_unif,zout
    READ(ifuInput,*)
    READ(ifuInput,*) uinit,vinit,winit, vper !new change for 2D-3D perturbations
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcx1,ux1,vx1,wx1,freq_ux1,freq_vx1,freq_wx1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcx2,ux2,vx2,wx2,freq_ux2,freq_vx2,freq_wx2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcy1,uy1,vy1,wy1,freq_uy1,freq_vy1,freq_wy1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcy2,uy2,vy2,wy2,freq_uy2,freq_vy2,freq_wy2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcz1,uz1,vz1,wz1,freq_uz1,freq_vz1,freq_wz1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcz2,uz2,vz2,wz2,freq_uz2,freq_vz2,freq_wz2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) pbcx1, pppx1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) pbcx2, pppx2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) pbcy1, pppy1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) pbcy2, pppy2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) pbcz1, pppz1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) pbcz2, pppz2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) no_tsteps,nmonitor,ndump,nrestart,nstat,nmonitor_probe_liftdrag, STATS_SUM
    READ(ifuInput,*)
    READ(ifuInput,*) format_dump
    READ(ifuInput,*)
    READ(ifuInput,*) re,dt
    READ(ifuInput,*)
    READ(ifuInput,*) frac_step_type, advec_scheme
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) alfa
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) internal_boundary_present
    READ(ifuInput,*)
    READ(ifuInput,*) body_type
    READ(ifuInput,*)
    READ(ifuInput,*) boundary_formulation
    READ(ifuInput,*)
    READ(ifuInput,*) probeLengthNormalized
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) iterMax_ad,restol_ad,omega_adv, ad_solver_type
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) pp_solver_type, iRedBlack
    READ(ifuInput,*)
    READ(ifuInput,*) omega
    READ(ifuInput,*)
    READ(ifuInput,*) iterMax_Poisson,restol_Poisson,iterResPoisson
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) mgcyclex, mgcycley, mgcyclez
    READ(ifuInput,*)
    READ(ifuInput,*) iterFinest, iterInter, iterCoarsest
    READ(ifuInput,*)
    READ(ifuInput,*) mgLevels_X, mgLevels_Y, mgLevels_Z
    READ(ifuInput,*)
    READ(ifuInput,*) infoconv
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*)turbActive
    READ(ifuInput,*)
    READ(ifuInput,*)turbModel
    READ(ifuInput,*)
    READ(ifuInput,*)cSmagFix
    READ(ifuInput,*)
    READ(ifuInput,*)filterWidthModel
    READ(ifuInput,*)
    READ(ifuInput,*)testFilterDir(DIRX),testFilterDir(DIRY),testFilterDir(DIRZ)
    READ(ifuInput,*)
    READ(ifuInput,*)homogDir(DIRX),homogDir(DIRY),homogDir(DIRZ)
    READ(ifuInput,*)
    READ(ifuInput,*)fWidthRatio(DIRX),fWidthRatio(DIRY),fWidthRatio(DIRZ)
    READ(ifuInput,*)
    READ(ifuInput,*)turbLagrTimeScale
    READ(ifuInput,*)
    READ(ifuInput,*)AdditionalOutput, BinaryOutput
    READ(ifuInput,*)
    READ(ifuInput,*)Full_Coarsening, cure_pressure_oscillations
    READ(ifuInput,*)
    READ(ifuInput,*)pressure_osc_velocity, pressure_osc_pressure
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*)optimization
    READ(ifuInput,*)
    READ(ifuInput,*)evaluationInterval, errorPermitted

    nxc=nx-1
    nyc=ny-1
    nzc=nz-1

! print input -----------------------------------------------------------------

    CALL print_inputs()

    IF ( turbActive == ACTIVE ) CALL TURB_print_inputs

! set appropriate flags -------------------------------------------------------

    IF ( internal_boundary_present == INTR_BOUND_NONE ) THEN
      boundary_formulation = NO_INTERNAL_BOUNDARY
      body_type            = NONE
    END IF ! internal_boundary_present

    PRINT*,' Internal_Boundary_Present (0: INTR_BOUND_NONE, 1: INTR_BOUND_PRESENT)', &
             internal_boundary_present
    PRINT*,' Body_Type            (1: General, 2: Canonical) = ',body_type
    PRINT*,' Boundary_Formulation (1: SSM, 2: GCM)           = ',boundary_formulation

    IF ( internal_boundary_present == INTR_BOUND_NONE ) THEN
      bodyInterceptWeight = zero
      imagePointWeight    = zero
      gcmFlag             = zero
    END IF ! internal_boundary_present

    boundary_motion = FIXED_BOUNDARY
    IF (iRedBlack == 0) THEN
        TNcolorX = 1
        TNcolorY = 1
        TNcolorZ = 1
        iStep = 1
        jStep = 1
        kStep = 1
    ELSE IF(iRedBlack == 1) THEN
        TNcolorX = 2
        TNcolorY = 1
        TNcolorZ = 1

        iStep = 1
        jStep = 2
        kStep = (ndim-DIM_2D) + 1
    ELSE
        TNcolorX = 2
        TNcolorY = 2
        TNcolorZ = 2

        iStep = 2
        jStep = 2
        kStep = (ndim-DIM_2D) + 1
    END IF

   END SUBROUTINE read_inputs
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE check_inputs()

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE turb_global_parameters
    USE turb_parameters

    IMPLICIT NONE

    INTEGER :: iBody,n

    IF ( nread > 1 .OR. nread < 0 ) THEN
      PRINT*,'Incorrect Value for NREAD in Input',nread
      PRINT*,'Either use 0 or 1'
      STOP
    END IF ! nread

    IF ( ndim > DIM_3D .OR. ndim < DIM_2D ) THEN
      PRINT*,'Incorrect Value for NDIM in Input',ndim
      PRINT*,'Either use 2 or 3'
      STOP
    END IF ! ndim

    IF ( nx < 3 ) THEN
      PRINT*,'Incorrect Value for NX in Input',nx
      PRINT*,'Use at Least 3'
      STOP
    END IF ! nx

    IF ( ny < 3 ) THEN
      PRINT*,'Incorrect Value for NY in Input',ny
      PRINT*,'Use at Least 3'
      STOP
    END IF ! ny

    IF ( nz < 3 ) THEN
      PRINT*,'Incorrect Value for NZ in Input',nz
      PRINT*,'Use at Least 3'
      STOP
    END IF ! nz

    IF ( xgrid_unif < UNIFORM_GRID .AND. xgrid_unif > NONUNIFORM_GRID ) THEN
      PRINT*,'Incorrect Value for XGRID_UNIF in Input',xgrid_unif
      PRINT*,'Use either 1 (Uniform) or 2 (Non-Uniform)'
      STOP
    END IF ! xgrid_unif

    IF ( ygrid_unif < UNIFORM_GRID .AND. ygrid_unif > NONUNIFORM_GRID ) THEN
      PRINT*,'Incorrect Value for YGRID_UNIF in Input',ygrid_unif
      PRINT*,'Use either 1 (Uniform) or 2 (Non-Uniform)'
      STOP
    END IF ! ygrid_unif

    IF ( zgrid_unif < UNIFORM_GRID .AND. zgrid_unif > NONUNIFORM_GRID ) THEN
      PRINT*,'Incorrect Value for ZGRID_UNIF in Input',zgrid_unif
      PRINT*,'Use either 1 (Uniform) or 2 (Non-Uniform)'
      STOP
    END IF ! zgrid_unif

    IF ( nDim == DIM_2D .AND. zgrid_unif /= UNIFORM_GRID ) THEN
      PRINT*,'Incorrect Value for ZGRID_UNIF in Input',zgrid_unif
      PRINT*,'Use a value of 1 (Uniform) for 2D simulations'
      STOP
    END IF ! nDim

    IF ( pp_solver_type <= 0 .OR. pp_solver_type > PP_SOLVER_TYPE_MG_Point_Jacobi ) THEN
      PRINT*,'Incorrect Value for pp_solver_type in Input',pp_solver_type
      PRINT*,'Either use 1 (LSOR) or 2 (PETSC) or 3 (MG), OR 4 (SIP), OR 5 (MG-SIP), OR 6 (MSIP), OR 7 (MG-MSIP)'
      STOP
    END IF ! pp_solver_type

    IF (ABS(alfa)>1.0E-3_CGREAL) THEN
      Hybrid=.TRUE.
    ELSE
      Hybrid=.FALSE.
    END IF

    IF ( iterMax_ad <= 0 ) THEN
      PRINT*,'Incorrect Value for iterMax_ad in Input',iterMax_ad
      PRINT*,'Use at least 1'
      STOP
    END IF ! iterMax_ad

    IF ( iterMax_Poisson <= 0 ) THEN
      PRINT*,'Incorrect Value for IterMaxPETSC in Input',iterMax_Poisson
      PRINT*,'Use at least 1'
      STOP
    END IF ! iterMax_Poisson

    IF ( restol_ad <= zero .AND. restol_ad > oned) THEN
      PRINT*,'Incorrect Value for restol_ad in Input',restol_ad
      PRINT*,'Use Positive value and less than 1'
      STOP
    END IF ! restol_ad

    IF ( restol_Poisson <= zero .AND. restol_Poisson > oned) THEN
      PRINT*,'Incorrect Value for restol_Poisson in Input',restol_Poisson
      PRINT*,'Use Positive value and less than 1'
      STOP
    END IF ! restol_Poisson

    IF ( iterResPoisson <= 0 ) THEN
      PRINT*,'Incorrect Value for IterResPoisson in Input',iterResPoisson
      PRINT*,'Use at least 1'
      STOP
    END IF ! iterResPoisson

    IF ( internal_boundary_present < INTR_BOUND_NONE  .OR.    &
         internal_boundary_present > INTR_BOUND_PRESENT       ) THEN
      PRINT*,'Incorrect Value for INTERNAL_BOUNDARY_PRESENT in Input',&
              internal_boundary_present
      PRINT*,'Either use 0 (None) or 1 (Present)'
      STOP
    END IF ! internal_boundary_present

    IF ( format_dump < TECPLOT .OR. format_dump > FIELDVIEW ) THEN
      PRINT*,'Incorrect Value for FORMAT_DUMP in Input',format_dump
      PRINT*,'Either use 1 (Tecplot) or 2 (FieldView)'
      STOP
    END IF ! format_dump

    IF ( nStat < STATS_NONE ) THEN
      PRINT*,'Incorrect Value for NSTAT in Input',nstat
      PRINT*,'Either use either 0 or Positive value'
      STOP
    END IF ! ndim

    IF ( bcx1 < BC_TYPE_DIRICHLET .OR. &
         bcx1 > BC_TYPE_SHEAR          ) THEN
      PRINT*,'Incorrect Value for bcx1 in Input',bcx1
      PRINT*,'Bounds are 1 (Dirichlet) to 7 (Shear)'
      STOP
    END IF ! bcx1

    IF ( bcx2 < BC_TYPE_DIRICHLET .OR. &
         bcx2 > BC_TYPE_SHEAR          ) THEN
      PRINT*,'Incorrect Value for bcx2 in Input',bcx2
      PRINT*,'Bounds are 1 (Dirichlet) to 7 (Shear)'
      STOP
    END IF ! bcx2

    IF ( bcy1 < BC_TYPE_DIRICHLET .OR. &
         bcy1 > BC_TYPE_SHEAR          ) THEN
      PRINT*,'Incorrect Value for bcy1 in Input',bcy1
      PRINT*,'Bounds are 1 (Dirichlet) to 7 (Shear)'
      STOP
    END IF ! bcy1

    IF ( bcy2 < BC_TYPE_DIRICHLET .OR. &
         bcy2 > BC_TYPE_SHEAR          ) THEN
      PRINT*,'Incorrect Value for bcx2 in Input',bcy2
      PRINT*,'Bounds are 1 (Dirichlet) to 7 (Shear)'
      STOP
    END IF ! bcy2

    IF ( bcz1 < BC_TYPE_DIRICHLET .OR. &
         bcz1 > BC_TYPE_SHEAR          ) THEN
      PRINT*,'Incorrect Value for bcy1 in Input',bcz1
      PRINT*,'Bounds are 1 (Dirichlet) to 7 (Shear)'
      STOP
    END IF ! bcz1

    IF ( bcz2 < BC_TYPE_DIRICHLET .OR. &
         bcz2 > BC_TYPE_SHEAR          ) THEN
      PRINT*,'Incorrect Value for bcz2 in Input',bcz2
      PRINT*,'Bounds are 1 (Dirichlet) to 7 (Shear)'
      STOP
    END IF ! bcz2

    IF ( bcx1 == BC_TYPE_PERIODIC .AND. &
         bcx2 /= BC_TYPE_PERIODIC       ) THEN
      PRINT*,'bcx1 is a PERIODIC boundary while bcx2 is not ',bcx2
      PRINT*,'Incorrect Value for bcx2 in Input ',bcx2
      PRINT*,'Set bcx2 in Input to ',BC_TYPE_PERIODIC
      STOP
    END IF ! bcx1

    IF ( bcx1 /= BC_TYPE_PERIODIC .AND. &
         bcx2 == BC_TYPE_PERIODIC       ) THEN
      PRINT*,'bcx2 is a PERIODIC boundary while bcx1 is not ',bcx1
      PRINT*,'Incorrect Value for bcx1 in Input ',bcx1
      PRINT*,'Set bcx1 in Input to ',BC_TYPE_PERIODIC
      STOP
    END IF ! bcx1

    IF ( bcy1 == BC_TYPE_PERIODIC .AND. &
         bcy2 /= BC_TYPE_PERIODIC       ) THEN
      PRINT*,'bcy1 is a PERIODIC boundary while bcy2 is not ',bcy2
      PRINT*,'Incorrect Value for bcy2 in Input ',bcy2
      PRINT*,'Set bcy2 in Input to ',BC_TYPE_PERIODIC
      STOP
    END IF ! bcy1

    IF ( bcy1 /= BC_TYPE_PERIODIC .AND. &
         bcy2 == BC_TYPE_PERIODIC       ) THEN
      PRINT*,'bcy2 is a PERIODIC boundary while bcy1 is not ',bcy1
      PRINT*,'Incorrect Value for bcy1 in Input ',bcy1
      PRINT*,'Set bcy1 in Input to ',BC_TYPE_PERIODIC
      STOP
    END IF ! bcy1

    IF ( bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 /= BC_TYPE_PERIODIC       ) THEN
      PRINT*,'bcz1 is a PERIODIC boundary while bcz2 is not ',bcz2
      PRINT*,'Incorrect Value for bcx2 in Input ',bcz2
      PRINT*,'Set bcz2 in Input to ',BC_TYPE_PERIODIC
      STOP
    END IF ! bcz1

    IF ( bcz1 /= BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC       ) THEN
      PRINT*,'bcz2 is a PERIODIC boundary while bcz1 is not ',bcz1
      PRINT*,'Incorrect Value for bcx1 in Input ',bcz1
      PRINT*,'Set bcz1 in Input to ',BC_TYPE_PERIODIC
      STOP
    END IF ! bcz1


    IF (pp_solver_type .EQ. PP_SOLVER_TYPE_MG) THEN

     IF ( nx-1-2**mgLevels_X .LT. 2) THEN
        PRINT*, 'Too many grid levels in x direction have been set!'
        STOP
     END IF

     IF ( ny-1-2**mgLevels_Y .LT. 2) THEN
        PRINT*, 'Too many grid levels in y direction have been set!'
        STOP
     END IF

     IF (NDIM  .EQ.  DIM_3D) THEN
       IF ( nz-1-2**mgLevels_Z .LT. 2) THEN
          PRINT*, 'Too many grid levels in z direction have been set!'
          STOP
       END IF
     END IF

    END IF

    IF ( frac_step_type < NO_VAN_KAN .OR. &
         frac_step_type > VAN_KAN          ) THEN
      PRINT*,'Incorrect Value for frac_step_type in Input',frac_step_type
      PRINT*,'Bounds are 0 (no-van-kan) to 1 (van-kan)'
      STOP
    END IF ! frac_step_type

    IF ( turbActive < INACTIVE .OR. turbActive > ACTIVE ) THEN
      WRITE(STDOUT,'(A,2X,I2)') 'Incorrect Value for turbActive', turbActive
      WRITE(STDOUT,'(A)') 'Use turbActive: 0 (Inactive Turbulence Model)'
      WRITE(STDOUT,'(A)') '             or 1 (  Active Turbulence Model)'
      STOP
    ENDIF ! turbActive

    IF ( turbActive == ACTIVE ) CALL TURB_check_inputs

   END SUBROUTINE check_inputs
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE print_inputs()

    USE global_parameters
    USE flow_parameters
    USE turb_global_parameters
    USE turb_parameters

    IMPLICIT NONE

    PRINT*,' NREAD in Input = ',nread
    PRINT*,' NDIM in Input  = ',ndim

    PRINT*,' NX in Input = ',nx
    PRINT*,' NY in Input = ',ny
    PRINT*,' NZ in Input = ',nz

    PRINT*,' XGRID_UNIF in Input = ',xgrid_unif
    PRINT*,' YGRID_UNIF in Input = ',ygrid_unif
    PRINT*,' ZGRID_UNIF in Input = ',zgrid_unif

    PRINT*,' pp_solver_type in Input = ',pp_solver_type
    PRINT*,' iterMax_ad        in Input = ',iterMax_ad
    PRINT*,' iterMax_Poisson   in Input = ',iterMax_Poisson

    PRINT*,' restol_ad       in Input = ',restol_ad
    PRINT*,' restol_Poisson  in Input = ',restol_Poisson
    PRINT*,' IterResPoisson in Input = ',iterResPoisson

    PRINT*,' INTERNAL_BOUNDARY_PRESENT in Input = ',&
              internal_boundary_present

    PRINT*,' FORMAT_DUMP in Input = ',format_dump
    PRINT*,' NSTAT       in Input = ',nstat

    PRINT*,' frac_step_type         in Input = ',frac_step_type

    PRINT*,' TurbActive (0: Inactive, 1: Active)', turbActive

   END SUBROUTINE print_inputs
!---------------------------------------------
!---------------------------------------------
   SUBROUTINE read_inputs_bound()

    USE global_parameters
    USE flow_parameters
    USE unstructured_surface_arrays
    USE mg_parameters
    USE fea_unstructure_surface
    USE usr_module
    USE body_dynamics

    IMPLICIT NONE


    REAL(KIND=CGREAL)  :: readDummyR, Aratio
    INTEGER           :: n,readDummyInt
    INTEGER           :: ntimePerCycle_check
    INTEGER           :: I, J

    REAL(KIND=CGREAL) :: xcent_tmp,ycent_tmp,zcent_tmp          !Added by Wanh
    REAL(KIND=CGREAL) :: vxcent_tmp,vycent_tmp,vzcent_tmp       !Added by Wanh
    REAL(KIND=CGREAL) :: angvx_tmp,angvy_tmp,angvz_tmp          !Added by Wanh
    REAL(KIND=CGREAL) :: angle_quat,axis(3)                     !Added by Geng

    DO I=1, 32
      J=ISHFT(1,I-1)
      OUTPUT_PARA(I)=J
    END DO

    IF ( internal_boundary_present == INTR_BOUND_NONE) GOTO 999

! set appropriate weight factors for GCM method
!
!   U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)
!    gp                       ip      b
!
    bodyInterceptWeight = probeLengthNormalized /(probeLengthNormalized-oned)
    imagePointWeight    = oned/(probeLengthNormalized-oned)

    IF (boundary_formulation == GCM_METHOD) THEN
      PRINT*,'   probeLengthNormalized = ',probeLengthNormalized
      PRINT*,'   bodyInterceptWeight = ',bodyInterceptWeight
      PRINT*,'   imagePointWeight    = ',imagePointWeight
    ENDIF

    gcmFlag = REAL((NO_INTERNAL_BOUNDARY - boundary_formulation),KIND=CGREAL)  &
             *REAL((SSM_METHOD           - boundary_formulation) ,KIND=CGREAL) &
             /REAL((NO_INTERNAL_BOUNDARY - GCM_METHOD),KIND=CGREAL) &
             /REAL((SSM_METHOD           - GCM_METHOD),KIND=CGREAL)

    PRINT*,'   gcmFlag    = ',gcmFlag

    IF ( boundary_formulation /= NO_INTERNAL_BOUNDARY .AND. &
         body_type            == CANONICAL                  ) THEN
      OPEN(ifuBodyIn,FILE='canonical_body_in.dat',STATUS='UNKNOWN')
      OPEN(ifuBodyOut,FILE='canonical_body_out.dat',STATUS='UNKNOWN')
      READ(ifuBodyIn,*) nbody, nbody_solid, nbody_membrane,nSection, channel_flow, zoneSeparate !nSection added by Wanh for partially dynamic coupling
      PRINT*,'   nBody = ',nBody
      PRINT*,'   nBody_solid =', nBody_solid
      PRINT*,'   nBody_membrane =', nBody_membrane

! allocate memory -------------------------------------------------------------

      ALLOCATE(canonical_body_type(nBody))
      ALLOCATE(body_dim(nBody))
      ALLOCATE(boundary_motion_type(nBody))

! new arrary for mixed bodies

      ALLOCATE(Fort_Formatted(nBody))

      ALLOCATE(wall_type(nBody))
      ALLOCATE(nPtsBodyMarkerOrig(nbody))
      ALLOCATE(nPtsBodyMarker(nbody))
      ALLOCATE(rigidRef1(nbody))
      ALLOCATE(rigidRef2(nbody))
      ALLOCATE(rigidRef3(nbody))
      ALLOCATE(zoneMax(nbody))   !added by Chengyu
      ALLOCATE(totNumTriElem(nBody))
      ALLOCATE(unstruc_surface_type(nBody))

      ALLOCATE(n_phi(nbody))
      ALLOCATE(n_theta(nbody))

      ALLOCATE(radiusx(nbody))
      ALLOCATE(radiusy(nbody))
      ALLOCATE(radiusz(nbody))
      ALLOCATE(alpha(nbody))
      ALLOCATE(cosalpha(nbody))
      ALLOCATE(sinalpha(nbody))
      ALLOCATE(xcent(nSection*nbody))
      ALLOCATE(ycent(nSection*nbody))
      ALLOCATE(zcent(nSection*nbody))
      ALLOCATE(vxcent(nSection*nbody))
      ALLOCATE(vycent(nSection*nbody))
      ALLOCATE(vzcent(nSection*nbody))
      ALLOCATE(angvx(nSection*nbody))
      ALLOCATE(angvy(nSection*nbody))
      ALLOCATE(angvz(nSection*nbody))
!new arrays for rotation
      ALLOCATE(xcentinit(nbody))
      ALLOCATE(ycentinit(nbody))
      ALLOCATE(zcentinit(nbody))
      ALLOCATE(vxcentTrans(nbody))
      ALLOCATE(vycentTrans(nbody))
      ALLOCATE(vzcentTrans(nbody))
      ALLOCATE(ampx(nbody))
      ALLOCATE(ampy(nbody))
      ALLOCATE(ampz(nbody))
      ALLOCATE(freqx(nbody))
      ALLOCATE(freqy(nbody))
      ALLOCATE(freqz(nbody))

      ALLOCATE(phase(nbody))
      ALLOCATE(cosphase(nbody))
      ALLOCATE(sinphase(nbody))
      ALLOCATE(angvxinit(nbody))
      ALLOCATE(angvyinit(nbody))
      ALLOCATE(angvzinit(nbody))
      ALLOCATE(ampangx(nbody))
      ALLOCATE(ampangy(nbody))
      ALLOCATE(ampangz(nbody))
      ALLOCATE(freqangx(nbody))
      ALLOCATE(freqangy(nbody))
      ALLOCATE(freqangz(nbody))
      ALLOCATE(xcentConstr(nBody),ycentConstr(nBody),zcentConstr(nBody)) ! VEERA - Flow induced Motion
      ALLOCATE(density_solid(nBody))! VEERA - Flow induced Motion
      ALLOCATE(ntimePerCycle(nbody))
      Allocate(DoF_on(nbody,6))


!      ALLOCATE(bodyMarkerForce(nPtsMax*3)) !added by CJ Yuan
!      ALLOCATE(bodyMarkerVel(nPtsMax*3)) !added by CJ Yuan
!      ALLOCATE(markerInterpolateIndex(nPtsMax*6)) ! Added by CJ Yuan July.13.2015
!      ALLOCATE(markerInterpolateRatio(nPtsMax*6)) ! Added by CJ Yuan July.13.2015
!      ALLOCATE(markerPressure(nPtsMax*2)) ! Added by CJ Yuan July.13.2015
!      ALLOCATE(markerInterpolateVelocity(nPtsMax*6)) ! Added by CJ Yuan July.13.2015

! initialize values and arrays ------------------------------------------------

      ntimePerCycle  = 0  ! added for no-stop feature
      nPtsBodyMarker = 0
      unstruc_surface_type  = 0

      n_phi       = 0
      n_theta     = 0

      radiusx     = zero
      radiusy     = zero
      radiusz     = zero

      xcent       = zero
      ycent       = zero
      zcent       = zero

      vxcent      = zero
      vycent      = zero
      vzcent      = zero

      angvx       = zero
      angvy       = zero
      angvz       = zero

!new arrays for rotation

      xcentinit   = zero
      ycentinit   = zero
      zcentinit   = zero

      vxcentTrans = zero
      vycentTrans = zero
      vzcentTrans = zero

      ampx        = zero
      ampy        = zero
      ampz        = zero

      freqx       = zero
      freqy       = zero
      freqz       = zero

      angvxinit   = zero
      angvyinit   = zero
      angvzinit   = zero

      ampangx     = zero
      ampangy     = zero
      ampangz     = zero

      freqangx    = zero
      freqangy    = zero
      freqangz    = zero


      DO n = 1,nBody
        unstruc_surface_type(n) = SOLID_BODY
      ENDDO

      IF (nBody_Membrane > 0) THEN
         DO n = nBody_Solid+1, nBody
            unstruc_surface_type(n) = MEMBRANE
         ENDDO
      ENDIF

! read input file pertinent to internal boundary ------------------------------

      PRINT*,'Reading Canonical_Body_In.dat File'
      DO n=1,nbody

        IF (zoneSeparate) THEN
            READ(ifuBodyIn,*)canonical_body_type(n),body_dim(n),boundary_motion_type(n), zoneMax(n) !added by Chengyu
        ELSE
            READ(ifuBodyIn,*)canonical_body_type(n),body_dim(n),boundary_motion_type(n)
        ENDIF

        READ(ifuBodyIn,*)wall_type(n)

        PRINT*,'  canonical_body_type   = ',canonical_body_type(n)
        PRINT*,'  body_type             = ',body_dim(n)
        PRINT*,'  boundary_motion_type  = ',boundary_motion_type(n)
        PRINT*,'  wall_type             = ',wall_type(n)

        IF ( boundary_motion_type(n) /= STATIONARY )  boundary_motion = MOVING_BOUNDARY

!--------- checking ntimePerCycle -------- added by Haibo-----------
        IF ( boundary_motion_type(n) == PRESCRIBED .OR. boundary_motion_type(n) == BIO_FOLLOWED_DYNAMICS_COUPLED .or. &
             boundary_motion_type(n) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
             boundary_motion_type(n) == DYNAMICS_COUPLED_SWIMMING  ) THEN
           CALL read_marker_vel_check(n, ntimePerCycle_check)
           ntimePerCycle(n) = ntimePerCycle_check
        END IF
!--------------------------------------------------------------------

!	Added by Wanh 05/05/10  by C Yuan 05/26/2015
	IF ( boundary_motion_type(n) == FEA_FLOW_STRUC_INTERACTION) THEN
!           FSI_on = .TRUE.
!           CALL read_marker_vel_check(n, ntimePerCycle_check)
!           ntimePerCycle(n) = ntimePerCycle_check
        END IF

	IF ( boundary_motion_type(n) == PARTIAL_DYNAMICS_COUPLED) THEN
           CALL read_marker_vel_check(n, ntimePerCycle_check)
           ntimePerCycle(n) = ntimePerCycle_check
        END IF

        SELECT CASE (canonical_body_type(n))

          CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)

	     READ(ifuBodyIn,*)nPtsBodyMarker(n),readDummyInt
             READ(ifuBodyIn,*)radiusx(n),radiusy(n)
             READ(ifuBodyIn,*)xcent(n),ycent(n)
             READ(ifuBodyIn,*)alpha(n)
             READ(ifuBodyIn,*)vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
             READ(ifuBodyIn,*)ampx(n),ampy(n),ampz(n)
             READ(ifuBodyIn,*)freqx(n),freqy(n),freqz(n)
             READ(ifuBodyIn,*)angvx(n),angvy(n),angvz(n)
             READ(ifuBodyIn,*)phase(n)          ! phase angle advance of angular velocity over translational velocity oscillation
             READ(ifuBodyIn,*)ampangx(n),ampangy(n),ampangz(n)
             READ(ifuBodyIn,*)freqangx(n),freqangy(n),freqangz(n)
	     READ(ifuBodyIn,*)density_fluid, density_solid(n)
	     READ(ifuBodyIn,*)xcentConstr(n),ycentConstr(n),zcentConstr(n)

          CASE(ELLIPSOID)

             READ(ifuBodyIn,*)n_theta(n),readDummyInt
             READ(ifuBodyIn,*)radiusx(n),radiusy(n),radiusz(n)
             READ(ifuBodyIn,*)xcent(n),ycent(n),zcent(n)
             READ(ifuBodyIn,*)alpha(n)                                     ! not active
             READ(ifuBodyIn,*)vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
             READ(ifuBodyIn,*)ampx(n),ampy(n),ampz(n)
             READ(ifuBodyIn,*)freqx(n),freqy(n),freqz(n)
             READ(ifuBodyIn,*)angvx(n),angvy(n),angvz(n)
             READ(ifuBodyIn,*)phase(n)
             READ(ifuBodyIn,*)ampangx(n),ampangy(n),ampangz(n)             ! not active
             READ(ifuBodyIn,*)freqangx(n),freqangy(n),freqangz(n)          ! not active
	     READ(ifuBodyIn,*)density_fluid, density_solid(n)
	     READ(ifuBodyIn,*)xcentConstr(n),ycentConstr(n),zcentConstr(n)

             IF (boundary_formulation == GCM_METHOD) THEN
                  Aratio = radiusz(n)/radiusx(n)
             ELSE
                  Aratio = oned
             END IF

             IF ( MOD(n_theta(n),2) == 0 ) THEN
               n_phi(n) = IDNINT(half*REAL(n_theta(n),KIND=CGREAL)*Aratio)
             ELSE
               n_phi(n) = IDNINT(half*REAL(n_theta(n)-1,KIND=CGREAL)*Aratio)+1
             ENDIF ! n_theta

             nPtsBodyMarker(n) = n_theta(n)*n_phi(n)

           CASE(UNSTRUCTURED_SURFACE)

! Note: for unstruc surface, nPtsBodyMarker = total number of nodes
	    !Yan_canonical_body_in_parameters
             READ(ifuBodyIn,*)DoF_on(n,1),DoF_on(n,2),DoF_on(n,3),DoF_on(n,4),DoF_on(n,5),DoF_on(n,6)  ! Added by G. Liu
	         READ(ifuBodyIn,*)nPtsBodyMarker(n),totNumTriElem(n)           ! nPtsBodyMarker(n)
             READ(ifuBodyIn,*)rigidRef1(n),rigidRef2(n),rigidRef3(n)
             READ(ifuBodyIn,*)readDummyR,readDummyR, readDummyR  ! radiusx(n),radiusy(n)
             READ(ifuBodyIn,*)xcent(n),ycent(n),zcent(n)
 !            READ(ifuBodyIn,*)readDummyR             ! alpha(n)
             READ(ifuBodyIn,*)quat_init(:)            ! added by Geng
             axis(1)=quat_init(2)
             axis(2)=quat_init(3)
             axis(3)=quat_init(4)
             angle_quat=quat_init(1)*pi/180.0         ! added by Geng
             call quaternion_form(axis,angle_quat,quat_init) ! added by Geng
             READ(ifuBodyIn,*)vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
             READ(ifuBodyIn,*)ampx(n),ampy(n),ampz(n)
             READ(ifuBodyIn,*)freqx(n),freqy(n),freqz(n)
             READ(ifuBodyIn,*)angvx(n),angvy(n),angvz(n)
             READ(ifuBodyIn,*)phase(n)
             READ(ifuBodyIn,*)ampangx(n),ampangy(n),ampangz(n)
             READ(ifuBodyIn,*)freqangx(n),freqangy(n),freqangz(n)
	         READ(ifuBodyIn,*)density_fluid, density_solid(n)
 !            write(*,*)'debug density_solid',density_solid(n)
	         READ(ifuBodyIn,*)xcentConstr(n),ycentConstr(n),zcentConstr(n)


             !IF ( boundary_motion_type(n) == PARTIAL_DYNAMICS_COUPLED) THEN

!               READ(ifuBodyIn,*)DynamicCharLength(n)      !Added by Wanh for partially dynamic coupling

                READ(ifuBodyIn,*)hingemarker
!                IF (nSection == 1) hingemarker = 1      !Subjected to change later.
                READ(ifuBodyIn,*)thickoverlength
                READ(ifuBodyIn,*)Limiter_angle_amp
                READ(ifuBodyIn,*)A_x, A_theta
                READ(ifuBodyIn,*) FreqN_torsion, Freq_prsb        !Natural freq. of spring and Presb. motion freq.
	            IF (body_dim(n) > 2) READ(ifuBodyIn,*) DepthOverLength
             !ENDIF

             !IF ( boundary_motion_type(n) == DYNAMICS_COUPLED) THEN
             !   READ(ifuBodyIn,*)thickoverlength
             !ENDIF

             alpha(n)     = zero

        END SELECT ! canonical_body_type

! initialize variables

        xcentinit(n) = xcent(n)
        ycentinit(n) = ycent(n)
        zcentinit(n) = zcent(n)

        vxcent(n)    = vxcentTrans(n)
        vycent(n)    = vycentTrans(n)
        vzcent(n)    = vzcentTrans(n)

        angvxinit(n) = angvx(n)
        angvyinit(n) = angvy(n)
        angvzinit(n) = angvz(n)

        cosalpha(n)  = COS(alpha(n)*PI/180.0_CGREAL)
        sinalpha(n)  = SIN(alpha(n)*PI/180.0_CGREAL)

        cosphase(n)  = COS(phase(n)*PI/180.0_CGREAL)
        sinphase(n)  = SIN(phase(n)*PI/180.0_CGREAL)

	IF ( boundary_motion_type(n) == FEA_FLOW_STRUC_INTERACTION) THEN
           ALLOCATE (Lambdau_Aitken(nBody))
           ALLOCATE (Lambdav_Aitken(nBody))
           ALLOCATE (Lambdaw_Aitken(nBody))
           Lambdau_Aitken(:) = 0
           Lambdav_Aitken(:) = 0
           Lambdaw_Aitken(:) = 0
        ENDIF


! write variables

        PRINT*,' Body Number = ',n
	PRINT*,'   nPtsBodyMarker = ',nPtsBodyMarker(n)
        PRINT*,'   Canonical Body Type = ', canonical_body_type(n)

        SELECT CASE(canonical_body_type(n))

          CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)
            PRINT*,'   RadiusX, RadiusY           = ',radiusx(n),radiusy(n)
            PRINT*,'   XCent, YCent               = ',xcent(n),ycent(n)

          CASE(ELLIPSOID)
            PRINT*,'   RadiusX, RadiusY, RadiusZ  = ',radiusx(n),radiusy(n),radiusz(n)
            PRINT*,'   XCent, YCent, ZCent        = ',xcent(n),ycent(n),zcent(n)

          CASE(UNSTRUCTURED_SURFACE)
           PRINT*,'   totNumTriElem =',totNumTriElem(n)
           PRINT*,'   XCent, YCent, ZCent        = ',xcent(n),ycent(n),zcent(n)
        END SELECT ! canonical_body_type

        PRINT*,'   Alpha                         = ',alpha(n)
        PRINT*,'   VXCentTrans, VYCentTrans, VZCentTrans     = ',&
                   vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
        PRINT*,'   AmpX, AmpY, AmpZ             = ',ampx(n),ampy(n),ampz(n)
        PRINT*,'   FreqX, FreqY, FreqZ          = ',freqx(n),freqy(n),freqz(n)
        PRINT*,'   AngVx, AngVy, AngVz          = ',angvxinit(n),angvyinit(n),angvzinit(n)
        PRINT*,'   AmpAngX, AmpAngY, AmpAngZ    = ',ampangx(n),ampangy(n),ampangz(n)
        PRINT*,'   FreqAngX, FreqAngY, FreqAngZ = ',freqangx(n),freqangy(n),freqangz(n)
	PRINT*,'   Density fluid, Density Solid = ',density_fluid, density_solid(n)

      ENDDO ! n

      if(channel_flow)then                        !add by Yan, for BC of inlet and outlet of channel flow
          READ(ifuBodyIn,*) nGate
          allocate(bcTypeGate(nGate))
          allocate(bcGateU(nGate))
          allocate(bcGateV(nGate))
          allocate(bcGateW(nGate))
          allocate(freqGate(nGate))
          do i=1,nGate
              READ(ifuBodyIn,*) bcTypeGate(i),bcGateU(i),bcGateV(i),bcGateW(i),freqGate(i)
          end do
      end if

      READ(ifuBodyIn,*)Prsb_MomentRef

      IF (Prsb_MomentRef) THEN
        OPEN(ifuPrsbMomentRef,FILE='input_momentref.dat',STATUS='UNKNOWN')
      ENDIF

!     Added by Wanh for partially dynamic coupling
      ALLOCATE(vxcent_prev(nSection))
      ALLOCATE(vycent_prev(nSection))
      ALLOCATE(vzcent_prev(nSection))

      ALLOCATE(vxcent_iter(nSection))
      ALLOCATE(vycent_iter(nSection))
      ALLOCATE(vzcent_iter(nSection))

      ALLOCATE(acc_xCG(nSection))
      ALLOCATE(acc_yCG(nSection))
      ALLOCATE(acc_zCG(nSection))

      ALLOCATE(angvx_old(nSection*nBody))
      ALLOCATE(angvy_old(nSection*nBody))
      ALLOCATE(angvz_old(nSection*nBody))

      ALLOCATE(angvx_iter(nSection*nBody))
      ALLOCATE(angvy_iter(nSection*nBody))
      ALLOCATE(angvz_iter(nSection*nBody))

      ALLOCATE(Converged_FSI(nBody))
      ALLOCATE(DynamicCharLength(nSection))
      ALLOCATE(section_xcent(nBody,nSection))
      ALLOCATE(section_ycent(nBody,nSection))
      ALLOCATE(section_zcent(nBody,nSection))
      ALLOCATE(LambdaVel_Aitken(nSection))
      ALLOCATE(LambdaAngVel_Aitken(nSection))
      ALLOCATE(deltaVxcent_Aitken(nSection))
      ALLOCATE(deltaVycent_Aitken(nSection))
      ALLOCATE(deltaVzcent_Aitken(nSection))
      ALLOCATE(deltaAngvx_Aitken(nSection))
      ALLOCATE(deltaAngvy_Aitken(nSection))
      ALLOCATE(deltaAngvz_Aitken(nSection))

      angvx_old   = zero
      angvy_old   = zero
      angvz_old   = zero

    ENDIF ! boundary_formulation

    IF(boundary_motion_type(1)==DYNAMICS_COUPLED_QUAT .or. boundary_motion_type(1)==DYNAMICS_COUPLED_MofI_QUAT .or. &
       boundary_motion_type(1)==DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
       boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING )THEN  ! Added by Geng
      allocate(mass_in(nBody))
      allocate(Icm_in(nBody,6))
      allocate(xcent_prev(nbody),ycent_prev(nbody),zcent_prev(nbody))
      allocate(angvx_noniner(nbody),angvy_noniner(nbody),angvz_noniner(nbody))

      open(ifuDynaNonIner,file='mass_I_noniner.dat')
      do n=1,nbody
        read(ifuDynaNonIner,*)mass_in(n)
        read(ifuDynaNonIner,*)Icm_in(n,1),Icm_in(n,2),Icm_in(n,3)
        read(ifuDynaNonIner,*)Icm_in(n,4),Icm_in(n,5),Icm_in(n,6)
!        write(*,*)'liu debug:', mass_in(n),n,'in main.f90'
!        write(*,*)'liu debug:', Icm_in(n,3)
!        stop
      enddo
      close(ifuDynaNonIner)

!      do n=1,nbody
!        mass_in(n)=mass_in(n)*density_solid(n)
!        Icm_in(n,:)=Icm_in(n,:)*density_solid(n)
!        write(*,*)'mass_in for DYNAMICS_COUPLED_QUAT: body ',n, mass_in(n)
!        write(*,*)'MofI for DYNAMICS_COUPLED_QUAT: body',n, Icm_in(n,:)
!      enddo

    ENDIF ! DYNAMICS_COUPLED_QUAT

999 CONTINUE

   END SUBROUTINE read_inputs_bound
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE check_inputs_bound()

    USE global_parameters
    USE flow_parameters
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER :: iBody,n
    REAL(KIND=CGREAL) :: rEps

    IF ( internal_boundary_present == INTR_BOUND_NONE) GOTO 999

    IF ( body_type < GENERAL .OR. body_type > CANONICAL         ) THEN
      PRINT*,'Incorrect Value of BODY_TYPE in Input ', body_type
      STOP
    ENDIF ! internal_boundary_present

    IF ( boundary_formulation    < SSM_METHOD       .OR.    &
         boundary_formulation    > GCM_METHOD               ) THEN
      PRINT*,' Internal Boundary Present in Calculations '
      PRINT*,'  Incorrect Value for BOUNDARY_FORMULATION in Input',boundary_formulation
      PRINT*,'  Either use 1 (SSM: StairStep Method) or 2 (GCM: GhostCell Method)'
      STOP
    END IF ! boundary_formulation

    IF ( body_type            /= CANONICAL .AND.  &
         boundary_formulation == GCM_METHOD  ) THEN
      PRINT*,' Unsupported Option in VICAR3D '
      PRINT*,' Cannot run non-canonical body in GCM mode yet'
      PRINT*,' Use  boundary_formulation  = 1 (SSM_METHOD)'
      STOP
    END IF ! boundary_formulation

    DO n=1,nbody
      IF ( body_type == CANONICAL ) THEN
      IF ( canonical_body_type(n) < ELLIPTIC_CYLINDER .OR. &
           canonical_body_type(n) > UNSTRUCTURED_SURFACE          ) THEN
        PRINT*,' Incorrect Value for CANONICAL_BODY_TYPE in Input',canonical_body_type(n)
        PRINT*,'  Either use 1 (Elliptic Cylinder), 2 (General Cylinder) '
        PRINT*,'             3 (Ellipsoid), 4 (Unstructured Surface)'
        STOP
      END IF ! canonical_body_type
      END IF ! body_type
    ENDDO

    DO n=1,nbody
!    print *,' boundary_motion_type(n)=', n, boundary_motion_type(n)

      IF ( boundary_motion_type(n) < STATIONARY .OR. &
!           boundary_motion_type(n) > PRESCRIBED    ) THEN     !changed by Wanh
           boundary_motion_type(n) > DYNAMICS_COUPLED_SWIMMING) THEN
        PRINT*,' Incorrect Value for BOUNDARY_MOTION_TYPE in Input', boundary_motion_type(n)
        PRINT*,' Either use 0 (Stationary), 1 (Forced Motion), 2 (Flow Induced), 3(Prescribed),'
	PRINT*,'   4(FEA_FLOW_STRUC_INTERACTION), 5(PARTIAL_DYNAMICS_COUPLED) 6(DYNAMICS_COUPLED)'
        STOP
      END IF ! boundary_motion
    ENDDO

    DO n=1,nbody
      IF ( wall_type(n) < NONPOROUS_AND_NONSLIP .AND. &
           wall_type(n) >   POROUS_OR_SLIP            ) THEN
        PRINT*,' Incorrect Value for WALL_TYPE in Input',wall_type(n)
        PRINT*,'  Either use 0 (NonPorous and NonSlip), or 1 (Porous or Slip)'
        STOP
      END IF ! wall_type
    ENDDO

    IF (  boundary_motion < FIXED_BOUNDARY .OR.       &
          boundary_motion > MOVING_BOUNDARY            ) THEN
      PRINT*,' Incorrect Value for BOUNDARY_MOTION in Input',boundary_motion
      PRINT*,'  Should be either 1 (Fixed Boundary) or 2 (Moving Boundary) '
      STOP
    END IF ! boundary_formulation

    IF ( ABS(probeLengthNormalized-2.0_CGREAL) > EPSILON(rEps) ) THEN
      PRINT*,' Unconventional Value for probeLengthNormalized in Input',probeLengthNormalized
      PRINT*,' Typical Value is  2.0'
      STOP
    END IF ! probeLengthNormalized

    DO iBody = 1, nBody
      IF ( nPtsBodyMarker(iBody) < 2) THEN
       PRINT*,' Incorrect Value for nPtsBodyMarker in Input',nPtsBodyMarker(iBody)
       PRINT*,' for body', iBody
       PRINT*,'  Use at least 2 Points to Mark the Body'
       STOP
      END IF ! nPtsBodyMarker

    END DO ! iBody

    DO n=1,nbody
      IF ( canonical_body_type(n) == UNSTRUCTURED_SURFACE .AND. &
           totNumTriElem(n)       < 2                           ) THEN
        PRINT*,' Incorrect Value for TotNumTriElem of Unstructured Surface',totNumTriElem(n)
        PRINT*,'  Use At Least 2 Elements '
        STOP
      END IF ! canonical_body_type
    ENDDO

    IF (boundary_formulation == SSM_METHOD .AND. flow_type == POTENTIAL_FLOW) THEN
      DO n=1,nbody
        IF ( boundary_motion_type(n) > STATIONARY  ) THEN                              ! simple but required that pgradx1 etc be made into
          PRINT*,' Cannot run potential flow with moving SSM boundary'                 ! 3D arrays. Decided not to bother .. RM
          PRINT*,' Should try running GCM'
          STOP
        END IF ! boundary_motion
        ENDDO
    ENDIF

999 CONTINUE

   END SUBROUTINE check_inputs_bound
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE init_simulation()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

#ifdef PETSC
    USE PETSC_MOD
#endif
    USE mg_parameters
    USE mg_arrays
    USE turb_parameters

    IMPLICIT NONE
    INTEGER :: ibody, i, j, k

    CALL allocate_memory()

    CALL BOUNDARY_allocate_memory()

    IF ( turbActive == ACTIVE ) THEN
      WRITE(*,*) 'CALLING TURB_allocate_memory from MAIN.F90'
      CALL TURB_allocate_memory()
    ENDIF ! turbActive

    CALL set_arrays_zero()

    CALL make_grid()

    CALL metrics()

    CALL init_flow()

! set read iblank flag

    readIblankFlag = .FALSE.
!    IF ( boundary_motion == MOVING_BOUNDARY .AND. nRead == 1) readIblankFlag = .TRUE.

! initialize markers and iblank

    IF ( body_type == CANONICAL) THEN
      PRINT*,'CALLING initialize_marker from MAIN.F90'
      CALL initialize_marker()    ! set initial location and velocity of markers
    END IF ! boundary_formulation

    IF ( boundary_formulation /= NO_INTERNAL_BOUNDARY ) THEN

        CALL set_initial_iblank()

    END IF ! boundary_formulation

!    IF ( boundary_motion == MOVING_BOUNDARY .AND. nRead == 1) readIblankFlag = .TRUE.

! set boundary

    PRINT*,'ENTERING SET_BOUNDARY'
    CALL set_boundary()

!    CALL write_dump()

! reset flag for non initialization stage

    readIblankFlag = .FALSE.

! set boundary arrays for non-Restart simulation

    IF ( nRead == 0 ) THEN
      PRINT*,'ENTERING SET_BC'
      CALL set_bc()       ! fill boundary arrays with u^0
    ENDIF

    PRINT*,'ENTERING FreshCell_CalcExpWeight'
    CALL FreshCell_CalcExpWeight()

    IF (nRead == 1) THEN
      CALL enforce_global_mass_consv()
      CALL write_monitor()
      CALL divergence()
    ENDIF

    SELECT CASE (pp_solver_type)
      CASE ( PP_SOLVER_TYPE_PETSC )
#ifdef PETSC
         PRINT*,'ENTERING PETSC_SETUP_SOLVER'
         CALL petsc_setup_solver()
#else
         PRINT*,'USER ERROR: PETSC Solver Active in Input Deck '
         PRINT*,'  Code not compiled with PETSC Flag'
         PRINT*,'  Code will stop and exit '
         PRINT*,'  Either set solver to MG or compile with PETSC=1'
         STOP
#endif
      CASE ( PP_SOLVER_TYPE_MG, PP_SOLVER_TYPE_MG_Point_Jacobi)
         CALL MG_initial()
      CASE ( PP_SOLVER_TYPE_MG_SIP )
         CALL MG_initial()
      CASE ( PP_SOLVER_TYPE_MG_MSIP )
         CALL MG_initial()

    END SELECT ! pp_solver_type

    PRINT*,'ENTERING SET_SOLVE_AD'
    CALL set_solve_ad()

    print*,'advec_scheme = (1 for AB2; 2 for CN1; 3 for CN2) ', &
            advec_scheme

   END SUBROUTINE init_simulation
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE flow_stop()

    IMPLICIT NONE

!   CALL statistics(1)

    CALL deallocate_memory()

   END SUBROUTINE flow_stop
!-------------------------------------------

!-------------------------------------------
   SUBROUTINE init_flow()

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE stat_arrays
    USE turb_parameters

    IMPLICIT NONE

    reinv   = oned/re
    dtinv   = oned/dt

    statCtr = 0

    IF (nread == 1) then
      PRINT*,'ENTERING READ_RESTART_FLOW'
      CALL read_restart_flow()

      PRINT*,'ENTERING READ_RESTART_BODY'
      CALL read_restart_body()

      IF (boundary_formulation == SSM_METHOD) THEN
        IF (nBody_solid > 0) THEN
          CALL SSM_set_internal_iup_solid()
          CALL SSM_set_internal_area()
        END IF
        IF (nBody_membrane > 0) THEN
          CALL SSM_set_internal_iup_membrane()
        END IF
      END IF
    ELSE

      CALL initialize_flowfield()

      CALL spanwise_pert() ! new for spanwise perturbation by Haibo

      ntime_start = 0
      IF ( turbActive == ACTIVE .AND. turbModel == TURB_MODEL_DYNLAGR ) THEN
        WRITE(STDOUT,*) '   Entering TURB_InitPhiField '
        CALL TURB_InitPhiField()
      END IF ! turbActive

    ENDIF

   END SUBROUTINE init_flow
!-------------------------------------------

!-------------------------------------------
   SUBROUTINE read_restart_flow()

    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE nlold_arrays
    USE boundary_arrays
    USE grid_arrays
    USE turb_parameters
    USE body_dynamics
    USE usr_module

    IMPLICIT NONE

    INTEGER           :: nxx,nyy,nzz,n
    REAL(KIND=CGREAL) :: dtr,rer,rEps
    INTEGER           :: i,j,k

    READ(ifuRstrtFlowIn)nxx,nyy,nzz
    IF ((nxx /= nx) .or. (nyy /= ny) .or. (nzz /= nz) ) THEN
      PRINT*,' FATAL ERROR in READ_RESTART_FLOW'
      PRINT*,'  CODE WILL CEASE ITS EXECUTION '
      PRINT*,'  Inconsistent Number of Grid points between Input and Restart'
      PRINT*,'  NX, NY, NZ in Input file   - ',nx,ny,nz
      PRINT*,'  NX, NY, NZ in Restart file - ',nxx,nyy,nzz
      PRINT*,'  POTENTIAL SOLUTION: Check values in input.dat'
      CALL abort_vicar3d(10)
    ENDIF

! Read restart file
! Load all values for bc[x,y,z][u,v,w] arrays
!  this kernel is needed since bcxu represents n-level value
!   for gradient, symmetry and periodic conditions
!   as well as boundary marker field

    READ(ifuRstrtFlowIn)ntime_start,time,dtr,rer
    READ(ifuRstrtFlowIn)u(0:nx+1,0:ny+1,0:nz+1), &
                        v(0:nx+1,0:ny+1,0:nz+1), &
                        w(0:nx+1,0:ny+1,0:nz+1), &
                   face_u(0:nx+1,0:ny+1,0:nz+1), &
                   face_v(0:nx+1,0:ny+1,0:nz+1), &
                   face_w(0:nx+1,0:ny+1,0:nz+1), &
                        p(0:nx+1,0:ny+1,0:nz+1), &
                   nluold(0:nx+1,0:ny+1,0:nz+1), &
                   nlvold(0:nx+1,0:ny+1,0:nz+1), &
                   nlwold(0:nx+1,0:ny+1,0:nz+1), &
                     bcxu(0:nx+1,0:ny+1,0:nz+1), &   ! Left Boundary
                     bcxv(0:nx+1,0:ny+1,0:nz+1), &
                     bcxw(0:nx+1,0:ny+1,0:nz+1), &
                     bcyu(0:nx+1,0:ny+1,0:nz+1), &   ! Right Boundary
                     bcyv(0:nx+1,0:ny+1,0:nz+1), &
                     bcyw(0:nx+1,0:ny+1,0:nz+1), &
                     bczu(0:nx+1,0:ny+1,0:nz+1), &   ! Back Boundary
                     bczv(0:nx+1,0:ny+1,0:nz+1), &
                     bczw(0:nx+1,0:ny+1,0:nz+1), &
                       pgradx1(0:ny+1,0:nz+1),   &
                       pgradx2(0:ny+1,0:nz+1),   &
                       pgrady1(0:nx+1,0:nz+1),   &
                       pgrady2(0:nx+1,0:nz+1),   &
                       pgradz1(0:nx+1,0:ny+1),   &
                       pgradz2(0:nx+1,0:ny+1),   &
                       viscTot(0:nx+1,0:ny+1,0:nz+1), &
                       bcxVisc(0:nx+1,0:ny+1,0:nz+1), &
                       bcyVisc(0:nx+1,0:ny+1,0:nz+1), &
                       bczVisc(0:nx+1,0:ny+1,0:nz+1)

    pPrime(0:nx+1,0:ny+1,0:nz+1) = p(0:nx+1,0:ny+1,0:nz+1)

    IF (ABS( dtr - dt ) .GT. EPSILON(rEps) ) THEN
      PRINT*,'WARNING in READ_RESTART_FLOW: Time Step size changed from Restart'
      PRINT*,' DT in Input file   - ',dt
      PRINT*,' DT in Restart file - ',dtr
    ENDIF ! dtr

    IF (ABS( rer - re ) .GT. EPSILON(rEps) ) THEN
      PRINT*,'WARNING in READ_RESTART_FLOW: Reynolds Number changed from Restart'
      PRINT*,' Re in Input file   - ',re
      PRINT*,' Re in Restart file - ',rer
    ENDIF ! rer

    IF ( turbActive == ACTIVE ) THEN
      PRINT*,'ENTERING TURB_READ_RESTART'
     CALL TURB_read_restart()
   END IF ! turbActive

    CLOSE(ifuRstrtFlowIn)

    write(*,*) '%%%%%%%%%%%%%%ntime_restart =', ntime_start

   END SUBROUTINE read_restart_flow
!-------------------------------------------

!-------------------------------------------
   SUBROUTINE read_restart_body()

    USE flow_parameters
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    USE usr_module
    USE body_dynamics

    IMPLICIT NONE

!... Loop variables
    INTEGER :: n,m,j,i

!... Local variables
    INTEGER :: boundary_motion_typeR, canonical_body_typeR, body_DimR, nptsBodyMarkerOrigR, &
               nBodyR,nPtsBodyMarkerR, ntimeR,totNumTriElemR,wall_typeR, nFreshR

    REAL(KIND=CGREAL) :: alphaR,ampxR,ampyR,ampzR,angvxR,angvyR,angvzR,   &
                         ampangxR,ampangyR,ampangzR,freqxR,freqyR,freqzR, &
                         freqangxR,freqangyR,freqangzR,phaseR,            &
                         radiusxR,radiusyR,radiuszR,rEps,rER, timeR, &
                         dtR, vxCentTransR, vyCentTransR, vzCentTransR

! Read Restart File

    READ(ifuRstrtBodyIn)ntimeR,timeR,dtR,reR
    READ(ifuRstrtBodyIn)nBodyR


! Trap Error for inconsistent timestep count

    IF ( ntimeR /= ntime_start ) THEN
      PRINT*,' FATAL ERROR in READ_RESTART_BODY'
      PRINT*,'  CODE WILL CEASE ITS EXECUTION '
      PRINT*,'  Inconsistent Time Step'
      PRINT*,'  ntime_start Read in restart_flow_in.dat    = ',ntime_start
      PRINT*,'  ntime_start Read in restart_body_in.dat    = ',ntimeR
      PRINT*,'  POTENTIAL SOLUTION: Check if restart files properly copied'
      STOP
    ENDIF ! ntimeR

! Trap Error for inconsistent body count

    IF ( nBodyR /= nBody ) THEN
      PRINT*,' FATAL ERROR in READ_RESTART_BODY'
      PRINT*,'  CODE WILL CEASE ITS EXECUTION '
      PRINT*,'  Inconsistent Number of Body Marker'
      PRINT*,'  nBody Read in canonical_body_in.dat    = ',nBody
      PRINT*,'  nBody Read from Restart File, restart_body_in.dat, = ',nBodyR
      PRINT*,'  POTENTIAL SOLUTION: Check nBody in canonical_body_in.dat'
      STOP
    ENDIF ! nBodyR

! Loop over all internal boundaries
!  Hard trap for nptsBodyMarker, totNumTriElem
!  canonical_body_type, [x,y,z,u,v,w]bodyMarker
!  Soft trap for radius,alpha,vxcentTrans, etc..

    DO n=1,nBody

      READ(ifuRstrtBodyIn)nptsBodyMarkerOrigR,nptsBodyMarkerR
      READ(ifuRstrtBodyIn)canonical_body_typeR, body_DimR, boundary_motion_typeR
      READ(ifuRstrtBodyIn)wall_typeR
      READ(ifuRstrtBodyIn)radiusxR,radiusyR,radiuszR
      READ(ifuRstrtBodyIn)xcent(n),ycent(n),zcent(n)   ! overwrite value from input file
      READ(ifuRstrtBodyIn)alpha(n)
      READ(ifuRstrtBodyIn)vxcentTransR,vycentTransR,vzcentTransR
      READ(ifuRstrtBodyIn)ampxR,ampyR,ampzR
      READ(ifuRstrtBodyIn)freqxR,freqyR,freqzR
      READ(ifuRstrtBodyIn)angvxR,angvyR,angvzR
      READ(ifuRstrtBodyIn)phaseR
      READ(ifuRstrtBodyIn)ampangxR,ampangyR,ampangzR
      READ(ifuRstrtBodyIn)freqangxR,freqangyR,freqangzR

! Trap hard errors

      IF ( nPtsBodyMarkerOrigR /= nPtsBodyMarkerOrig(n) ) THEN
        PRINT*,' FATAL ERROR in READ_RESTART_BODY '
        PRINT*,'  CODE WILL CEASE ITS EXECUTION '
        PRINT*,'  Inconsistent Number of Body Marker Points for Body Number  =', n
        PRINT*,'  nPtsBodyMarker Read from canonical_body_in.dat = ',nPtsBodyMarkerOrig(n)
        PRINT*,'  nPtsBodyMarker Read from restart_body_in.dat   = ',nPtsBodyMarkerOrigR
        PRINT*,'  POTENTIAL SOLUTION: Check nPtsBodyMarker in canonical_body_in.dat'
        STOP
      ENDIF ! nptsBodyMarkerR

      IF ( canonical_body_typeR /= canonical_body_type(n) ) THEN
        PRINT*,' FATAL ERROR in READ_RESTART_BODY '
        PRINT*,'  CODE WILL CEASE ITS EXECUTION '
        PRINT*,'  Inconsistent Canonical Body Type for Body Number = ', n
        PRINT*,'  canonical_body_type Read from canonical_body_in.dat = ',canonical_body_type(n)
        PRINT*,'  canonical_body_type Read from restart_body_in.dat   = ',canonical_body_typeR
        PRINT*,'  POTENTIAL SOLUTION: Check canonical_body_type in canonical_body_in.dat'
        STOP
      ENDIF ! canonical_body_typeR

      IF ( body_DimR /= body_Dim(n) ) THEN
        PRINT*,' FATAL ERROR in READ_RESTART_BODY '
        PRINT*,'  CODE WILL CEASE ITS EXECUTION '
        PRINT*,'  Inconsistent Canonical Body Type for Body Number = ', n
        PRINT*,'  body_Dim Read from canonical_body_in.dat = ',body_Dim(n)
        PRINT*,'  body_Dim Read from restart_body_in.dat   = ',body_DimR
        PRINT*,'  POTENTIAL SOLUTION: Check body_Dim in canonical_body_in.dat'
        STOP
      ENDIF ! body_DimR

! - Read Body Markers -------------------------------------------------------------------------

      READ(ifuRstrtBodyIn)xBodyMarker(n,1:nPtsBodyMarker(n)), &
                          yBodyMarker(n,1:nPtsBodyMarker(n)), &
                          zBodyMarker(n,1:nPtsBodyMarker(n))

      IF (canonical_body_type(n) == ELLIPTIC_CYLINDER   .OR.  &
          canonical_body_type(n) == GENERAL_CYLINDER    .OR.  &
          canonical_body_type(n) == UNSTRUCTURED_SURFACE    ) THEN
        READ(ifuRstrtBodyIn)totNumTriElemR

        IF ( totNumTriElemR /= totNumTriElem(n) .AND. &
             canonical_body_type(n) == UNSTRUCTURED_SURFACE ) THEN
          PRINT*,' FATAL ERROR in READ_RESTART_BODY'
          PRINT*,'  CODE WILL CEASE ITS EXECUTION '
          PRINT*,'  Inconsistent Total Number of Triangular Elements for Body Number = ', n
          PRINT*,'  totNumTriElem Read from unstruc_surface_in.dat = ',totNumTriElem(n)
          PRINT*,'  totNumTriElem Read from restart_body_in.dat    = ',totNumTriElemR
          PRINT*,'  POTENTIAL SOLUTION: Check totNumTriElem in unstruc_surface_in.dat'
          STOP
        ENDIF ! totNumTriElemR

        READ(ifuRstrtBodyIn)triElemNeig(n,1,1:totNumTriElem(n)), &
                            triElemNeig(n,2,1:totNumTriElem(n)), &
                            triElemNeig(n,3,1:totNumTriElem(n))
        READ(ifuRstrtBodyIn)normDirFlag
      ENDIF ! canonical_body_type

      READ(ifuRstrtBodyIn)uBodyMarker(n,1:nPtsBodyMarker(n)), &
                          vBodyMarker(n,1:nPtsBodyMarker(n)), &
                          wBodyMarker(n,1:nPtsBodyMarker(n))

! Trap hard error for radius

      SELECT CASE( canonical_body_type(n) )

       CASE (ELLIPTIC_CYLINDER:GENERAL_CYLINDER)
          IF ( ABS( radiusxR - radiusx(n) ) .GT. EPSILON(rEps) .OR. &
               ABS( radiusyR - radiusy(n) ) .GT. EPSILON(rEps) ) THEN
            PRINT*,' ERROR in READ_RESTART_BODY-CODE WILL CEASE ITS EXECUTION '
            PRINT*,'  Inconsistent Radius for Body Number = ', n
            PRINT*,'  radiusx,  radiusy Read from canonical_body_in.dat = ',&
                      radiusx(n), radiusy(n)
            PRINT*,'  radiusx,  radiusy Read from restart_body_in.dat   = ',&
                      radiusxR, radiusyR
            PRINT*,'  POTENTIAL SOLUTION: Check radiusx, radiusy in canonical_body_in.dat'
            STOP
          ENDIF ! radiusxR

       CASE (ELLIPSOID)
          IF ( ABS( radiusxR - radiusx(n) ) .GT. EPSILON(rEps) .OR. &
               ABS( radiusyR - radiusy(n) ) .GT. EPSILON(rEps) .OR. &
               ABS( radiuszR - radiusz(n) ) .GT. EPSILON(rEps)      ) THEN
            PRINT*,' ERROR in READ_RESTART_BODY-CODE WILL CEASE ITS EXECUTION '
            PRINT*,'  Inconsistent Radius for Body Number = ', n
            PRINT*,'  radiusx,  radiusy, radiusz Read from canonical_body_in.dat = ',&
                      radiusx(n), radiusy(n), radiusz(n)
            PRINT*,'  radiusx,  radiusy, radiusz Read from restart_body_in.dat   = ',&
                      radiusxR, radiusyR, radiuszR
            PRINT*,'  POTENTIAL SOLUTION: Check radiusx, radiusy, radiusz in canonical_body_in.dat'
            STOP
          ENDIF ! radiusxR

      END SELECT ! canonical_body_type

! - Trap soft errors

      IF ( boundary_motion_typeR /=  boundary_motion_type(n) ) THEN
        PRINT*,' WARNING in READ_RESTART_BODY '
        PRINT*,'  boundary_motion_type has been modified in Input File for Body Number = ', n
        PRINT*,'  boundary_motion_type Read from canonical_body_in.dat = ', boundary_motion_type(n)
        PRINT*,'  boundary_motion_type Read from restart_body_in.dat   = ', boundary_motion_typeR
      ENDIF ! boundary_motion_typeR

      IF ( wall_typeR /=  wall_type(n) ) THEN
        PRINT*,' WARNING in READ_RESTART_BODY '
        PRINT*,'  wall_type has been modified in Input File for Body Number = ', n
        PRINT*,'  wall_type Read from canonical_body_in.dat = ', wall_type(n)
        PRINT*,'  wall_type Read from restart_body_in.dat   = ', wall_typeR
      ENDIF ! wall_typeR

      IF ( ABS( vxcentTransR - vxcentTrans(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( vycentTransR - vycentTrans(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( vzcentTransR - vzcentTrans(n) ) .GT. EPSILON(rEps)      ) THEN
        PRINT*,' WARNING in READ_RESTART_BODY '
        PRINT*,'  v[x,y,z]centTrans has been modified in Input File for Body Number = ', n
        PRINT*,'  vxcentTrans,vycentTrans,vzcentTrans Read from canonical_body_in.dat = ', &
                  vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
        PRINT*,'  vxcentTrans,vycentTrans,vzcentTrans Read from restart_body_in.dat   = ', &
                  vxcentTransR,vycentTransR,vzcentTransR
      ENDIF ! vxcentTransR

      IF ( ABS( ampxR - ampx(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( ampyR - ampy(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( ampzR - ampz(n) ) .GT. EPSILON(rEps)      ) THEN
        PRINT*,' WARNING in READ_RESTART_BODY '
        PRINT*,'  amp[x,y,z] has been modified in Input File for Body Number = ', n
        PRINT*,'  ampx,ampy,ampz Read from canonical_body_in.dat = ', ampx(n),ampy(n),ampz(n)
        PRINT*,'  ampx,ampy,ampz Read from restart_body_in.dat   = ', ampxR,ampyR,ampzR
      ENDIF ! ampxR

      IF ( ABS( freqxR - freqx(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( freqyR - freqy(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( freqzR - freqz(n) ) .GT. EPSILON(rEps)      ) THEN
        PRINT*,' WARNING in READ_RESTART_BODY '
        PRINT*,'  freq[x,y,z] has been modified in Input File for Body Number = ', n
        PRINT*,'  freqx,freqy,freqz Read from canonical_body_in.dat = ', freqx(n),freqy(n),freqz(n)
        PRINT*,'  freqx,freqy,freqz Read from restart_body_in.dat   = ', freqxR,freqyR,freqzR
      ENDIF ! freqxR

      IF ( ABS( phaseR - phase(n) ) .GT. EPSILON(rEps) ) THEN
        PRINT*,' WARNING in READ_RESTART_BODY '
        PRINT*,'  Phase has been modified in Input File for Body Number = ', n
        PRINT*,'  Phase Read from canonical_body_in.dat  = ', phase(n)
        PRINT*,'  Phase Read from restart_body_in.dat    = ', phaseR
      ENDIF ! phaseR

      IF ( ABS( freqangxR - freqangx(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( freqangyR - freqangy(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( freqangzR - freqangz(n) ) .GT. EPSILON(rEps)      ) THEN
        PRINT*,' WARNING in READ_RESTART_BODY '
        PRINT*,'  freqang[x,y,z] has been modified in Input File for Body Number = ', n
        PRINT*,'  freqangx,freqangy,freqangz Read from canonical_body_in.dat = ', &
                  freqangx(n),freqangy(n),freqangz(n)
        PRINT*,'  freqangx,freqangy,freqangz Read from restart_body_in.dat   = ', &
                  freqangxR,freqangyR,freqangzR
      ENDIF ! freqxR

      IF ( ABS( ampangxR - ampangx(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( ampangyR - ampangy(n) ) .GT. EPSILON(rEps) .OR. &
           ABS( ampangzR - ampangz(n) ) .GT. EPSILON(rEps)      ) THEN
        PRINT*,' WARNING in READ_RESTART_BODY '
        PRINT*,'  ampang[x,y,z] has been modified in Input File for Body Number = ', n
        PRINT*,'  ampangx,ampangy,ampangz Read from canonical_body_in.dat = ', &
                  ampangx(n),ampangy(n),ampangz(n)
        PRINT*,'  ampangx,ampangy,ampangz Read from restart_body_in.dat   = ', &
                  ampangxR,ampangyR,ampangzR
      ENDIF ! ampxR

    ENDDO ! n

! Read iblank and fresh_cell

    READ(ifuRstrtBodyIn)iblank(0:nx+1,0:ny+1,0:nz+1), &
                        fresh_cell(0:nx+1,0:ny+1,0:nz+1)

    IF ( boundary_formulation == GCM_METHOD )   &
      READ(ifuRstrtBodyIn)bodyNum(0:nx+1,0:ny+1,0:nz+1)
    READ(ifuRstrtBodyIn)nFreshR

    PRINT*,'Number of Fresh Cells read from restrt file = ',nFreshR

    CLOSE(ifuRstrtBodyIn)


    READ(ifuRstrtDynamicsIn,*) Limiter_angle_amp
    READ(ifuRstrtDynamicsIn,*) Limiter_On
    READ(ifuRstrtDynamicsIn,*) aero_moment, aero_moment_threshold
    READ(ifuRstrtDynamicsIn,*) moment_z, moment_z_pretime
    DO i = 2,nSection
       READ(ifuRstrtDynamicsIn,*) angvz(i)
!       READ(ifuRstrtDynamicsIn,*) angvz_prev(i)
    ENDDO

    CLOSE (ifuRstrtDynamicsIn)

   END SUBROUTINE read_restart_body
!-------------------------------------------

!-------------------------------------------
   SUBROUTINE write_restart()

    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE nlold_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    USE turb_global_parameters
    USE turb_parameters

    USE body_dynamics           !Added by Wanh
    USE usr_module              !Added by Wanh

    IMPLICIT NONE

    INTEGER :: n,i

    PRINT*,'Writing out restart file: ntime' ,ntime

! Open Files
!  Use 2 restart files to avoid clobbering problems.

    SELECT CASE(indexRstrt)
      CASE(1)
        OPEN(ifuRstrtFlowOut,FILE='restart_flow_out1.dat',FORM='UNFORMATTED')
        OPEN(ifuRstrtBodyOut,FILE='restart_body_out1.dat',FORM='UNFORMATTED')
        IF (turbActive == ACTIVE .AND. turbModel == TURB_MODEL_DYNLAGR) &
          OPEN(ifuRstrtTurbOut,FILE='restart_turb_out1.dat',FORM='UNFORMATTED')

        OPEN(ifuRstrtDynamicsOut,FILE='restart_dynamics_out1.dat')

        indexRstrt = 2

      CASE(2)
        OPEN(ifuRstrtFlowOut,FILE='restart_flow_out2.dat',FORM='UNFORMATTED')
        OPEN(ifuRstrtBodyOut,FILE='restart_body_out2.dat',FORM='UNFORMATTED')
        IF (turbActive == ACTIVE .AND. turbModel == TURB_MODEL_DYNLAGR) &
          OPEN(ifuRstrtTurbOut,FILE='restart_turb_out2.dat',FORM='UNFORMATTED')

        OPEN(ifuRstrtDynamicsOut,FILE='restart_dynamics_out2.dat')

        indexRstrt = 1

    END SELECT ! indexRstrt


! Flow Field

! Write restart file
! Load all values for bc[x,y,z][u,v,w] arrays
!  this kernel is needed since bcxu represents n-level value
!   for gradient, symmetry and periodic conditions
!   as well as boundary marker field

    WRITE(ifuRstrtFlowOut)nx,ny,nz
    WRITE(ifuRstrtFlowOut)ntime,time,dt,re
    WRITE(ifuRstrtFlowOut)u(0:nx+1,0:ny+1,0:nz+1), &
                          v(0:nx+1,0:ny+1,0:nz+1), &
                          w(0:nx+1,0:ny+1,0:nz+1), &
                     face_u(0:nx+1,0:ny+1,0:nz+1), &
                     face_v(0:nx+1,0:ny+1,0:nz+1), &
                     face_w(0:nx+1,0:ny+1,0:nz+1), &
                          p(0:nx+1,0:ny+1,0:nz+1), &
                     nluold(0:nx+1,0:ny+1,0:nz+1), &
                     nlvold(0:nx+1,0:ny+1,0:nz+1), &
                     nlwold(0:nx+1,0:ny+1,0:nz+1), &
                       bcxu(0:nx+1,0:ny+1,0:nz+1), &   ! Left Boundary
                       bcxv(0:nx+1,0:ny+1,0:nz+1), &
                       bcxw(0:nx+1,0:ny+1,0:nz+1), &
                       bcyu(0:nx+1,0:ny+1,0:nz+1), &   ! Right Boundary
                       bcyv(0:nx+1,0:ny+1,0:nz+1), &
                       bcyw(0:nx+1,0:ny+1,0:nz+1), &
                       bczu(0:nx+1,0:ny+1,0:nz+1), &   ! Back Boundary
                       bczv(0:nx+1,0:ny+1,0:nz+1), &
                       bczw(0:nx+1,0:ny+1,0:nz+1), &
                       pgradx1(0:ny+1,0:nz+1),     &
                       pgradx2(0:ny+1,0:nz+1),     &
                       pgrady1(0:nx+1,0:nz+1),     &
                       pgrady2(0:nx+1,0:nz+1),     &
                       pgradz1(0:nx+1,0:ny+1),     &
                       pgradz2(0:nx+1,0:ny+1),     &
                       viscTot(0:nx+1,0:ny+1,0:nz+1), &
                       bcxVisc(0:nx+1,0:ny+1,0:nz+1), &
                       bcyVisc(0:nx+1,0:ny+1,0:nz+1), &
                       bczVisc(0:nx+1,0:ny+1,0:nz+1)

! Body Field
!  Write restart file for moving boundary case only


       WRITE(ifuRstrtBodyOut)ntime,time,dt,re
       WRITE(ifuRstrtBodyOut)nBody

       DO n=1,nBody

         WRITE(ifuRstrtBodyOut)nptsBodyMarkerOrig(n),nptsBodyMarker(n)
         WRITE(ifuRstrtBodyOut)canonical_body_type(n),body_dim(n),boundary_motion_type(n)
         WRITE(ifuRstrtBodyOut)wall_type(n)
         WRITE(ifuRstrtBodyOut)radiusx(n),radiusy(n),radiusz(n)
         WRITE(ifuRstrtBodyOut)xcent(n),ycent(n),zcent(n)
         WRITE(ifuRstrtBodyOut)alpha(n)
         WRITE(ifuRstrtBodyOut)vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
         WRITE(ifuRstrtBodyOut)ampx(n),ampy(n),ampz(n)
         WRITE(ifuRstrtBodyOut)freqx(n),freqy(n),freqz(n)
         WRITE(ifuRstrtBodyOut)angvx(n),angvy(n),angvz(n)
         WRITE(ifuRstrtBodyOut)phase(n)
         WRITE(ifuRstrtBodyOut)ampangx(n),ampangy(n),ampangz(n)
         WRITE(ifuRstrtBodyOut)freqangx(n),freqangy(n),freqangz(n)

         WRITE(ifuRstrtBodyOut)xBodyMarker(n,1:nPtsBodyMarker(n)), &
                               yBodyMarker(n,1:nPtsBodyMarker(n)), &
                               zBodyMarker(n,1:nPtsBodyMarker(n))

         IF (canonical_body_type(n) == ELLIPTIC_CYLINDER   .OR.  &
             canonical_body_type(n) == GENERAL_CYLINDER    .OR.  &
             canonical_body_type(n) == UNSTRUCTURED_SURFACE    ) THEN
            WRITE(ifuRstrtBodyOut)totNumTriElem(n)
            WRITE(ifuRstrtBodyOut)triElemNeig(n,1,1:totNumTriElem(n)), &
                                  triElemNeig(n,2,1:totNumTriElem(n)), &
                                  triElemNeig(n,3,1:totNumTriElem(n))
            WRITE(ifuRstrtBodyOut)normDirFlag
         ENDIF ! canonical_body_type

         WRITE(ifuRstrtBodyOut)uBodyMarker(n,1:nPtsBodyMarker(n)), &
                               vBodyMarker(n,1:nPtsBodyMarker(n)), &
                               wBodyMarker(n,1:nPtsBodyMarker(n))
       ENDDO ! n

       WRITE(ifuRstrtBodyOut)iblank(0:nx+1,0:ny+1,0:nz+1), &
                         fresh_cell(0:nx+1,0:ny+1,0:nz+1)
       IF(boundary_formulation == GCM_METHOD)   &
             WRITE(ifuRstrtBodyOut)bodyNum(0:nx+1,0:ny+1,0:nz+1)
       WRITE(ifuRstrtBodyOut)nFresh

    IF ( turbActive == ACTIVE ) THEN
      PRINT*,'ENTERING TURB_WRITE_RESTART'
      CALL TURB_write_restart()
    END IF ! turbActive

    WRITE(ifuRstrtDynamicsOut,*) Limiter_angle_amp                      !Added by Wanh
    WRITE(ifuRstrtDynamicsOut,*) Limiter_On                             !Added by Wanh
    WRITE(ifuRstrtDynamicsOut,*) aero_moment, aero_moment_threshold     !Added by Wanh
    WRITE(ifuRstrtDynamicsOut,*) moment_z, moment_z_pretime             !Added by Wanh
    DO i = 2,nSection
       WRITE(ifuRstrtDynamicsOut,*) angvz(i)                            !Added by Wanh
!       WRITE(ifuRstrtDynamicsOut,*) angvz_prev(i)                            !Added by Wanh
    ENDDO

    CLOSE(ifuRstrtFlowOut)
    CLOSE(ifuRstrtBodyOut)
    CLOSE(ifuRstrtTurbOut)

    CLOSE(ifuRstrtDynamicsOut)                  !Added by Wanh

   END SUBROUTINE write_restart
!-------------------------------------------

!-------------------------------------------
   SUBROUTINE set_arrays_zero()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE nlold_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE stat_arrays

    IMPLICIT NONE

    time    = zero

    x       = zero
    y       = zero
    z       = zero

    dx      = zero
    dy      = zero
    dz      = zero

    dxc     = zero
    dyc     = zero
    dzc     = zero

    dxinv   = zero
    dyinv   = zero
    dzinv   = zero

    dxcinv  = zero
    dycinv  = zero
    dzcinv  = zero

    fx      = zero
    fy      = zero
    fz      = zero

    iblank      = 0
    fresh_cell  = 0
    exp_weight  = zero

    iup     = 0
    ium     = 0
    jup     = 0
    jum     = 0
    kup     = 0
    kum     = 0

    u       = zero
    v       = zero
    w       = zero

    face_u  = zero
    face_v  = zero
    face_w  = zero

    p       = zero
    pPrime  = zero
    pgradx1 = zero
    pgradx2 = zero
    pgrady1 = zero
    pgrady2 = zero
    pgradz1 = zero
    pgradz2 = zero

    nlu     = zero
    nlv     = zero
    nlw     = zero

    nluold  = zero
    nlvold  = zero
    nlwold  = zero

    uTilde  = zero
    vTilde  = zero
    wTilde  = zero

    bcxu    = zero
    bcyu    = zero
    bczu    = zero

    bcxv    = zero
    bcyv    = zero
    bczv    = zero

    bcxw    = zero
    bcyw    = zero
    bczw    = zero

    amx     = zero
    apx     = zero
    acx     = zero

    amy     = zero
    apy     = zero
    acy     = zero

    amz     = zero
    apz     = zero
    acz     = zero

    rhs     = zero
    dummy   = zero

    face1   = zero
    face2   = zero

    viscTot = zero
    bcxvisc = zero
    bcyvisc = zero
    bczvisc = zero

    IF ( nStat > STATS_NONE ) THEN
      uAv       = zero
      vAv       = zero
      wAv       = zero
      pAv       = zero
      uvAv      = zero
      uwAv      = zero
      vwAv      = zero
      uuAv      = zero
      vvAv      = zero
      wwAv      = zero
    END IF ! nStat


   END SUBROUTINE set_arrays_zero
!-------------------------------------------

!-------------------------------------------
   SUBROUTINE initialize_flowfield()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays

    IMPLICIT NONE

    u(0:,0:,0:)       = uinit
    v(0:,0:,0:)       = vinit
    w(0:,0:,0:)       = winit
    face_u(0:,0:,0:)  = uinit
    face_v(0:,0:,0:)  = vinit
    face_w(0:,0:,0:)  = winit

    p(0:,0:,0:)       = zero
    pPrime(0:,0:,0:)  = zero

    viscTot(0:,0:,0:) = reinv
    bcxvisc(0:,0:,0:) = reinv
    bcyvisc(0:,0:,0:) = reinv
    bczvisc(0:,0:,0:) = reinv

   END SUBROUTINE initialize_flowfield
!-------------------------------------------
!-------------------------------------------
   SUBROUTINE spanwise_pert()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays

    IMPLICIT NONE

    INTEGER :: I, J, K
    REAL(KIND=CGREAL) :: har1, har2

    DO K = 0, nz
       DO J = 0, ny
           DO I = 0, nx
                   CALL RANDOM_NUMBER(har1)
                   CALL RANDOM_NUMBER(har2)
                   u(i,j,k) = u(i,j,k) + vper*(har1-har2)
                   v(i,j,k) = v(i,j,k) + vper*(har1-har2)
                   w(i,j,k) = w(i,j,k) + vper*(har1-har2)
                   face_u(i,j,k) = face_u(i,j,k) + vper*(har1-har2)
                   face_v(i,j,k) = face_v(i,j,k) + vper*(har1-har2)
                   face_w(i,j,k) = face_w(i,j,k) + vper*(har1-har2)
            END DO
        END DO
    END DO

   END SUBROUTINE spanwise_pert
!-------------------------------------------

!-------------------------------------------
   SUBROUTINE abort_vicar3d(abort_code)

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: abort_code

!
    SELECT CASE (abort_code)
      CASE(10)
        PRINT*,'Grid sizes in restart and input files do not match'
      CASE(20)
        PRINT*,'Iblank cell sticking out on its own'
    END SELECT

    STOP

   END SUBROUTINE abort_vicar3d
!-------------------------------------------
