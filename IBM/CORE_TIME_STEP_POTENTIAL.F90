!-----------------------------------------------------------------------------
   SUBROUTINE time_step_potential()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE stat_arrays
#ifdef PETSC
    use petsc_mod
#endif
    USE turb_parameters
    
    IMPLICIT NONE

!... Local variables

    INTEGER              :: iBody,i,j,k
    REAL(KIND=CGREAL)    :: sum, phiTMax, phiTMin,xhat,yhat,vdot
    REAL(KIND=CGREAL)    :: starttime, endtime
    INTEGER              :: clock1, clock2, clock_rate
  
!******************************************************************************

! Start time stepping loop

    DO ntime = ntime_start+1,ntime_start+no_tsteps        ! move solution from n --> n+1

      time = time + dt

! starting with a zero guess for potential
      pPrime = zero

      IF( MOD(ntime,nmonitor) == 0 ) THEN
        WRITE(*,*)'======================================================='
        WRITE(*,'(A,I6,A,F15.7,A,I6)') 'NTIME = ',ntime,',  Time = ',time,',  NTIME TO GO = ',no_tsteps+ntime_start-ntime
      ENDIF

!------------------------------------------------------------------------------
!     Move Body Boundary and compute coefficients
!------------------------------------------------------------------------------

      IF ( boundary_motion == MOVING_BOUNDARY .AND. ntime > 1) THEN
        CALL move_boundary()   ! move boundary to (n+1)

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

      END IF ! boundary_motion

      DO iBody = 1,nBody
        IF (WALL_TYPE(iBody) == POROUS_OR_SLIP) THEN  
          CALL wall_velocity(iBody)
        ENDIF                                        
      END DO

      IF( MOD(ntime,nmonitor) == 0 ) THEN
        WRITE(*,*) 'AREAX1; AREAX2 = ',areax1,areax2
        WRITE(*,*) 'AREAY1; AREAY2 = ',areay1,areay2
        WRITE(*,*) 'AREAZ1; AREAZ2 = ',areaz1,areaz2
      END IF

!     WRITE(*,*) 'Entering SET_BC '
      CALL set_bc()    ! fill boundary arrays with (u)^n+1

! set boundary conditions for velocity at outer boundaries.      
      DO k = 1,nz-1
      DO j = 1,ny-1
        pgradx1(j,k) = bcxu(1,j,k)
        pgradx2(j,k) = bcxu(nx-1,j,k)
      ENDDO
      ENDDO

      DO k = 1,nz-1
      DO i = 1,nx-1
        pgrady1(i,k) = bcyv(i,1,k)
        pgrady2(i,k) = bcyv(i,ny-1,k)
      ENDDO
      ENDDO

      DO j = 1,ny-1
      DO i = 1,nx-1
        pgradz1(i,j) = bczw(i,j,1)
        pgradz2(i,j) = bczw(i,j,nz-1)
      ENDDO
      ENDDO

! Set RHS of Poisson Equation to Zero 
      nlu = zero

!------------------------------------------------------------------------------
!    Solve the Poisson Pressure Equation (PPE)
!------------------------------------------------------------------------------

      call system_clock(clock1)
      CALL solve_poisson()
      call system_clock(clock2, clock_rate)
      IF ( MOD(ntime,nmonitor) == 0 ) & 
        WRITE(*,*) '*****Total time for solving poisson eq is:', &
         REAL(clock2-clock1)/REAL(clock_rate)
	!endtime-starttime

! Compute potential velocity
      CALL gradient(pPrime)
      u   = nlu
      v   = nlv	
      w   = nlw	

! Compute d(phi)/dt term
      DO k=1,nz-1
      DO j=1,ny-1
      DO i=1,nx-1
        p(i,j,k)  = ( pPrime(i,j,k) - p(i,j,k) )/dt
        p(i,j,k)  = p(i,j,k)*REAL(1-iblank(i,j,k),KIND=CGREAL)
      ENDDO
      ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      write(701,*)'VARIABLES="X","Y","PHIC","PHITC","PHITE"'
!      write(701,*)'ZONE F=POINT, I=',nx-1,', J=',ny-1
!      k = 1
!      do j=1,ny-1
!      do i=1,nx-1
!         xhat = xc(i) - xcent(1)
!         yhat = yc(j) - ycent(1)
!         vdot = ampy(1)*2.0_CGREAL*PI*freqy(1)*COS(2.0_CGREAL*PI*freqy(1)*time)
!         nlv(i,j,k) = -(-vycent(1)**2 + vdot*yhat)* 0.25_CGREAL/(xhat**2 + yhat**2)
!         write(701,129)xc(i),yc(j),pPrime(i,j,k),p(i,j,k),nlv(i,j,k)*(oned - REAL(iblank(i,j,k),KIND=CGREAL))
!      enddo
!      enddo
!129	format(5(2x,e12.5))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      phiTMax = MAXVAL(p(1:nx-1,1:ny-1,1:nz-1))
      phiTMin = MINVAL(p(1:nx-1,1:ny-1,1:nz-1))

      PRINT*,'Max Value of PhiT = ',phiTMax
      PRINT*,'Min Value of PhiT = ',phiTMin
 
! compute pressure using Bernoulli Equation
      DO k=1,nz-1
      DO j=1,ny-1
      DO i=1,nx-1
        p(i,j,k) = ( half*(oned - u(i,j,k)*u(i,j,k) - v(i,j,k)*v(i,j,k) )   &
                     - p(i,j,k) )*REAL(1-iblank(i,j,k),KIND=CGREAL)
      ENDDO
      ENDDO
      ENDDO

!------------------------------------------------------------------------------
!    Monitor output
!------------------------------------------------------------------------------

      IF ( MOD(ntime,nmonitor) == 0 ) THEN
        CALL write_monitor()
      ENDIF

      IF ( nStat > STATS_NONE ) THEN
        IF ( MOD(ntime,nstat) == 0 ) THEN
         CALL calc_statistics(0)
         statCtr = statCtr + 1
        ENDIF
      END IF ! nStat

      IF ( MOD(ntime,nmonitor_probe_liftdrag) == 0 ) THEN
        CALL write_probe_files()

        SELECT CASE(boundary_formulation)

          CASE (SSM_METHOD)
            CALL drag_lift_potential()

          CASE (GCM_METHOD)
            CALL GCM_drag_lift_potential()

         END SELECT ! boundary_formulation
      ENDIF

      IF ( MOD(ntime,ndump) == 0 )    THEN
         CALL write_dump()
      ENDIF

      IF ( MOD(ntime,nrestart) == 0 .OR. ntime==ntime_start+no_tsteps) THEN
          CALL write_restart()

        IF ( nStat > STATS_NONE ) CALL calc_statistics(1)
      ENDIF

! Save current potential in p
      p  = pPrime 

      phiTMax = MAXVAL(p(1:nx-1,1:ny-1,1:nz-1))
      phiTMin = MINVAL(p(1:nx-1,1:ny-1,1:nz-1))

      PRINT*,'Max Value of PhiOld = ',phiTMax
      PRINT*,'Min Value of PhiOld = ',phiTMin


    ENDDO ! ntime
       
   END SUBROUTINE time_step_potential

