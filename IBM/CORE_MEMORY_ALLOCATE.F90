!--------------------------------------------------
!   SUBROUTINE allocate_memory()
!   SUBROUTINE deallocate_memory()
!--------------------------------------------------

   SUBROUTINE allocate_memory()

    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE nlold_arrays
    USE multiuse_arrays
    USE grid_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE stat_arrays
    USE usr_module ,ONLY : scx,scy,scz,scmx,scmy,scmz  ! VEERA  -- 4 Flow-Induced Motion - Boundary_marker_vel.F90
    USE body_dynamics, ONLY : nSection
    USE Pressure_Aitken_Array                   !Wanh

    IMPLICIT NONE
    
    INTEGER :: nmax
    INTEGER :: iErr

    ALLOCATE(x(0:nx+1)     ,y(0:ny+1)     ,z(0:nz+1)     )
    ALLOCATE(xc(0:nx+1)    ,yc(0:ny+1)    ,zc(0:nz+1)    )
    ALLOCATE(dx(0:nx+1)    ,dy(0:ny+1)    ,dz(0:nz+1)    )
    ALLOCATE(dxc(0:nx+1)   ,dyc(0:ny+1)   ,dzc(0:nz+1)   )
    ALLOCATE(dxinv(0:nx+1) ,dyinv(0:ny+1) ,dzinv(0:nz+1) )
    ALLOCATE(dxcinv(0:nx+1),dycinv(0:ny+1),dzcinv(0:nz+1))
    ALLOCATE(fx(0:nx+1)    ,fy(0:ny+1)    ,fz(0:nz+1)    )

    ALLOCATE(u(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(v(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(w(0:nx+1,0:ny+1,0:nz+1))

    ALLOCATE(u_bak(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(v_bak(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(w_bak(0:nx+1,0:ny+1,0:nz+1))

    ALLOCATE(face_u(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(face_v(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(face_w(0:nx+1,0:ny+1,0:nz+1))

!    ALLOCATE(Usign(0:nx+1,0:ny+1,0:nz+1))    !
!    ALLOCATE(Vsign(0:nx+1,0:ny+1,0:nz+1))    !  Added by Rupesh
!    ALLOCATE(Wsign(0:nx+1,0:ny+1,0:nz+1))    !

    ALLOCATE(p(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(pPrime(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(pgradx1(0:ny+1,0:nz+1))
    ALLOCATE(pgradx2(0:ny+1,0:nz+1))
    ALLOCATE(pgrady1(0:nx+1,0:nz+1))
    ALLOCATE(pgrady2(0:nx+1,0:nz+1))
    ALLOCATE(pgradz1(0:nx+1,0:ny+1))
    ALLOCATE(pgradz2(0:nx+1,0:ny+1))

    ALLOCATE(nlu(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(nlv(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(nlw(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(nlu0(0:nx+1,0:ny+1,0:nz+1)) !
    ALLOCATE(nlw0(0:nx+1,0:ny+1,0:nz+1)) !

    ALLOCATE(nluold(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(nlvold(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(nlwold(0:nx+1,0:ny+1,0:nz+1))

    IF (boundary_motion_type(1) == FEA_FLOW_STRUC_INTERACTION .OR.  &
        boundary_motion_type(1) == PARTIAL_DYNAMICS_COUPLED .OR.    &
        boundary_motion_type(1) == DYNAMICS_COUPLED .OR. &
        boundary_motion_type(1) == BIO_DYNAMICS_COUPLED.or. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT .OR. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. & 
        boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING) THEN
    ALLOCATE(nlu_FSI(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(nlv_FSI(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(nlw_FSI(0:nx+1,0:ny+1,0:nz+1))
    ENDIF

    ALLOCATE(uTilde(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(vTilde(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(wTilde(0:nx+1,0:ny+1,0:nz+1))

    ALLOCATE(bcxu(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(bcxv(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(bcxw(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(bcyu(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(bcyv(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(bcyw(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(bczu(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(bczv(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(bczw(0:nx+1,0:ny+1,0:nz+1))

    ALLOCATE(amx(0:nx+1),acx(0:nx+1),apx(0:nx+1))
    ALLOCATE(amy(0:ny+1),acy(0:ny+1),apy(0:ny+1))
    ALLOCATE(amz(0:nz+1),acz(0:nz+1),apz(0:nz+1))

    nmax = nx+ny+nz
    ALLOCATE(rhs(0:nmax+1),dummy(0:nmax+1))
    ALLOCATE(face1(0:nmax+1),face2(0:nmax+1))

    ALLOCATE(amx_ad(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(apx_ad(0:nx+1,0:ny+1,0:nz+1))
    
    ALLOCATE(amy_ad(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(apy_ad(0:nx+1,0:ny+1,0:nz+1))

    ALLOCATE(amz_ad(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(apz_ad(0:nx+1,0:ny+1,0:nz+1))

    IF (nStat > STATS_NONE ) THEN
      ALLOCATE(uav(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(vav(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(wav(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(pav(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(uvAv(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(uwAv(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(vwAv(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(uuAv(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(vvAv(0:nx+1,0:ny+1,0:nz+1))
      ALLOCATE(wwAv(0:nx+1,0:ny+1,0:nz+1))
    ENDIF ! nStat
    
    if(pressure_osc_velocity.or.pressure_osc_pressure)then                                 !add by Yan
        allocate(hybrid_mark(0:nx+1,0:ny+1,0:nz+1))
    end if

!  ALLOCATE(scx(nBody)) 	! VEERA  -- 4 Flow-Induced Motion - Boundary_marker_vel.F90
!  ALLOCATE(scy(nBody))  ! VEERA
!  ALLOCATE(scz(nBody))	! VEERA     
!  ALLOCATE(scmx(nBody)) ! VEERA
!  ALLOCATE(scmy(nBody)) ! VEERA
!  ALLOCATE(scmz(nBody)) ! VEERA
  
  ALLOCATE(scx(nSection*nBody))       !Changed by Wanh  
  ALLOCATE(scy(nSection*nBody))  
  ALLOCATE(scz(nSection*nBody))
  ALLOCATE(scmx(nSection*nBody)) 
  ALLOCATE(scmy(nSection*nBody))
  ALLOCATE(scmz(nSection*nBody))

  IF (pressureAitkenOn) THEN
     ALLOCATE(deltapPrime_Lplus1(0:nx+1,0:ny+1,0:nz+1))       !Added by Wanh for pressure Aitken acceleration.
     ALLOCATE( deltapPrime_prevL(0:nx+1,0:ny+1,0:nz+1))       !Added by Wanh for pressure Aitken acceleration.
     ALLOCATE( diff_deltapPrimeL(0:nx+1,0:ny+1,0:nz+1))       !Added by Wanh for pressure Aitken acceleration.
     ALLOCATE(pPrime_L(0:nx+1,0:ny+1,0:nz+1))                 !Added by Wanh for pressure Aitken acceleration.
  ENDIF


!------------------------------------------------------------------------------
!   Arrays pertinent to viscosity
!------------------------------------------------------------------------------

    ALLOCATE(viscTot(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for viscTot'
      STOP
    ENDIF ! iErr 

    ALLOCATE(bcxvisc(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcxvisc'
      STOP
    ENDIF ! iErr
      
    ALLOCATE(bcyvisc(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcyvisc'
      STOP
    ENDIF ! iErr

    ALLOCATE(bczvisc(0:nx+1,0:ny+1,0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bczvisc'
      STOP
    ENDIF ! iErr
        
  END SUBROUTINE allocate_memory
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE deallocate_memory()

    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE nlold_arrays
    USE multiuse_arrays
    USE grid_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE stat_arrays
    USE mg_arrays
    USE usr_module ,ONLY : scx,scy,scz,scmx,scmy,scmz  ! VEERA  -- 4 Flow-Induced Motion - Boundary_marker_vel.F90 
    USE Pressure_Aitken_Array                   !Wanh
    
    IMPLICIT NONE

    DEALLOCATE(x  ,y  ,z)
    DEALLOCATE(xc ,yc ,zc)
    DEALLOCATE(dx ,dy ,dz)
    DEALLOCATE(dxc,dyc,dzc)
    DEALLOCATE(dxinv ,dyinv ,dzinv)
    DEALLOCATE(dxcinv,dycinv,dzcinv)
    DEALLOCATE(fx    ,fy    ,fz)

    DEALLOCATE(iblank)
    DEALLOCATE(fresh_cell)
    DEALLOCATE(exp_weight)

    DEALLOCATE(iup)
    DEALLOCATE(ium)
    DEALLOCATE(jup)
    DEALLOCATE(jum)
    DEALLOCATE(kup)
    DEALLOCATE(kum)

    IF (Hybrid) THEN
      DEALLOCATE(iupp)    !
      DEALLOCATE(iumm)    !
      DEALLOCATE(jupp)    !  Added by Rupesh (used for 2nd Upwinding)
      DEALLOCATE(jumm)    !
      DEALLOCATE(kupp)    !
      DEALLOCATE(kumm)    !
    ENDIF

    DEALLOCATE(u)
    DEALLOCATE(v)
    DEALLOCATE(w)

    DEALLOCATE(u_bak)
    DEALLOCATE(v_bak)
    DEALLOCATE(w_bak)

    DEALLOCATE(face_u)
    DEALLOCATE(face_v)
    DEALLOCATE(face_w)

!    DEALLOCATE(Usign)  !
!    DEALLOCATE(Vsign)  !  Added by Rupesh (used for 2nd Upwinding)
!    DEALLOCATE(Wsign)  !

    DEALLOCATE(p)
    DEALLOCATE(pPrime)
    DEALLOCATE(pgradx1)
    DEALLOCATE(pgradx2)
    DEALLOCATE(pgrady1)
    DEALLOCATE(pgrady2)
    DEALLOCATE(pgradz1)
    DEALLOCATE(pgradz2)

    DEALLOCATE(nlu)
    DEALLOCATE(nlv)
    DEALLOCATE(nlw)

    DEALLOCATE(nluold)
    DEALLOCATE(nlvold)
    DEALLOCATE(nlwold)

    IF (boundary_motion_type(1) == FEA_FLOW_STRUC_INTERACTION .OR.  &
        boundary_motion_type(1) == PARTIAL_DYNAMICS_COUPLED .OR.    &
        boundary_motion_type(1) == DYNAMICS_COUPLED .OR. &
        boundary_motion_type(1) == BIO_DYNAMICS_COUPLED .or. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT .OR. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
        boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING) THEN
    DEALLOCATE(nlu_FSI)
    DEALLOCATE(nlv_FSI)
    DEALLOCATE(nlw_FSI)
    ENDIF

    DEALLOCATE(uTilde)
    DEALLOCATE(vTilde)
    DEALLOCATE(wTilde)

    DEALLOCATE(bcxu)
    DEALLOCATE(bcxv)
    DEALLOCATE(bcxw)
    DEALLOCATE(bcyu)
    DEALLOCATE(bcyv)
    DEALLOCATE(bcyw)
    DEALLOCATE(bczu)
    DEALLOCATE(bczv)
    DEALLOCATE(bczw)

    DEALLOCATE(amx,acx,apx)
    DEALLOCATE(amy,acy,apy)
    DEALLOCATE(amz,acz,apz)
    DEALLOCATE(rhs,dummy)
    DEALLOCATE(face1,face2)
    
    DEALLOCATE(amx_ad)
    DEALLOCATE(apx_ad)    
    DEALLOCATE(amy_ad)
    DEALLOCATE(apy_ad)
    DEALLOCATE(amz_ad)
    DEALLOCATE(apz_ad)

    DEALLOCATE(angvx, angvy, angvz)
    DEALLOCATE(angvx_old, angvy_old, angvz_old)

    IF (pp_solver_type == PP_SOLVER_TYPE_MG) THEN
      DEALLOCATE (mgrid_I, mgrid_J, mgrid_K)
      DEALLOCATE (dxcinv_MG, dxinv_MG)
      DEALLOCATE (dycinv_MG, dyinv_MG)
      DEALLOCATE (dzcinv_MG, dzinv_MG)
    END IF ! pp_solver_type

    IF ( nStat > STATS_NONE ) THEN
      DEALLOCATE(uav)
      DEALLOCATE(vav)
      DEALLOCATE(wav)
      DEALLOCATE(pav)
      DEALLOCATE(uuav)
      DEALLOCATE(uvav)
      DEALLOCATE(uwav)
      DEALLOCATE(vvav)
      DEALLOCATE(vwav)
      DEALLOCATE(wwav)
    END IF ! nStat
    
  DEALLOCATE(scx)  ! VEERA  -- 4 Flow-Induced Motion - Boundary_marker_vel.F90
  DEALLOCATE(scy)  ! VEERA
  DEALLOCATE(scz)  ! VEERA     
  DEALLOCATE(scmx) ! VEERA
  DEALLOCATE(scmy) ! VEERA
  DEALLOCATE(scmz) ! VEERA
  DEALLOCATE(xcentConstr)   ! VEERA
  DEALLOCATE(ycentConstr)   ! VEERA
  DEALLOCATE(zcentConstr)   ! VEERA
  DEALLOCATE(density_solid) ! VEERA

  IF (pressureAitkenOn) THEN
     DEALLOCATE(deltapPrime_Lplus1)       ! Added by Wanh for pressure Aitken
     DEALLOCATE(deltapPrime_prevL)        ! Added by Wanh for pressure Aitken
     DEALLOCATE(diff_deltapPrimeL)        ! Added by Wanh for pressure Aitken
     DEALLOCATE(pPrime_L)       ! Added by Wanh for pressure Aitken
  ENDIF

  if(pressure_osc_velocity.or.pressure_osc_pressure)then                                 !add by Yan
        deallocate(hybrid_mark)
  end if
    
  END SUBROUTINE deallocate_memory
!------------------------------------------------------------------------------
