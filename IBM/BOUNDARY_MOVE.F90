! Move body markers and centroids and set new velocity of marker points and centroid

SUBROUTINE move_boundary()

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface
USE usr_module
USE body_dynamics

IMPLICIT NONE

INTEGER :: i,n,ifortCent,iBody,iFort,m,j,k,i_flag,presb_flag,iSection
INTEGER :: time_momentref
INTEGER :: LE_marker, TE_marker, iMarker

! The following moment_z_pretime is for restart.
IF (NTIME == 1) THEN
  moment_z_pretime = 1e-20
ELSE
  moment_z_pretime = moment_z
ENDIF

! Update body velocity
DO iBody=1,nBody
  SELECT CASE (boundary_motion_type(iBody))
  CASE (FORCED)
    angvx_old(ibody) = angvx(ibody)
    angvy_old(ibody) = angvy(ibody)
    angvz_old(ibody) = angvz(ibody)
    CALL forced_motion(iBody)
    CALL compute_marker_vel(iBody)
! print *, 'max/min ubodymarker:', maxval(ubodymarker), minval(ubodymarker)
  CASE (FLOW_INDUCED)
    angvx_old(ibody) = angvx(ibody)
    angvy_old(ibody) = angvy(ibody)
    angvz_old(ibody) = angvz(ibody)
    CALL flow_induced_motion(iBody)
    CALL compute_marker_vel(iBody)
  CASE (PRESCRIBED)
    CALL read_marker_vel(iBody)
  CASE (FEA_FLOW_STRUC_INTERACTION) !Added by Wanh 05/05/10 Modified by CJ Yuan Jun.3.2015
    PRINT*,'   FEA_FLOW_STRUC_INTERACTION vega_vel_update'
    CALL vega_vel_update ! added by CJ Yuan Jun.3.2015

  CASE (PARTIAL_DYNAMICS_COUPLED)
!For partially prescribed and partially dynamic coupling motion, the
!vxent/vycent/angvz are cacluated based on the body with dynamic coupling oly.
    IF (niterFS < 1e-1) CALL read_prescribed_marker_vel(iBody)
!vxcent_prev(iBody) = vxcent(iBody)
!vycent_prev(iBody) = vycent(iBody)
    DO iSection = 2, Nsection
      vxcent_prev(iSection) = vxcent(iSection)
      vycent_prev(iSection) = vycent(iSection)
      vzcent_prev(iSection) = vzcent(iSection)
      angvx_old(iSection) = angvx(iSection)
      angvy_old(iSection) = angvy(iSection)
      angvz_old(iSection) = angvz(iSection)
    ENDDO
    CALL compute_marker_vel_section(iBody)

  CASE (DYNAMICS_COUPLED)
    vxcent_prev(iBody) = vxcent(iBody)
    vycent_prev(iBody) = vycent(iBody)
    vzcent_prev(iBody) = vzcent(iBody)
    angvx_old(ibody) = angvx(ibody)
    angvy_old(ibody) = angvy(ibody)
    angvz_old(ibody) = angvz(ibody)
    CALL compute_marker_vel(iBody)

  CASE ( DYNAMICS_COUPLED_QUAT )
    vxcent_prev(iBody) = vxcent(iBody)
    vycent_prev(iBody) = vycent(iBody)
    vzcent_prev(iBody) = vzcent(iBody)
    angvx_old(ibody) = angvx(ibody)
    angvy_old(ibody) = angvy(ibody)
    angvz_old(ibody) = angvz(ibody)
    CALL compute_marker_vel(iBody)

  CASE ( DYNAMICS_COUPLED_MofI_QUAT )
    vxcent_prev(iBody) = vxcent(iBody)
    vycent_prev(iBody) = vycent(iBody)
    vzcent_prev(iBody) = vzcent(iBody)
    angvx_old(ibody) = angvx(ibody)
    angvy_old(ibody) = angvy(ibody)
    angvz_old(ibody) = angvz(ibody)
    CALL compute_marker_vel(iBody)

  CASE (BIO_DYNAMICS_COUPLED)
    vxcent_prev(1) = vxcent(1)
    vycent_prev(1) = vycent(1)
    vzcent_prev(1) = vzcent(1)
    angvx_old(1) = angvx(1)
    angvy_old(1) = angvy(1)
    angvz_old(1) = angvz(1)
    CALL compute_marker_vel(iBody)
    DO i=1,nPtsBodyMarker(iBody)
    xBodyMarker(iBody,i) = xBodyMarker(iBody,i) + dt*uBodyMarker(iBody,i)
    yBodyMarker(iBody,i) = yBodyMarker(iBody,i) + dt*vBodyMarker(iBody,i)
    zBodyMarker(iBody,i) = zBodyMarker(iBody,i) + dt*wBodyMarker(iBody,i)
    ENDDO ! i

  CASE (BIO_FOLLOWED_DYNAMICS_COUPLED)
    call read_marker_vel(iBody)
    call correct_marker_vel(iBody)
    call compute_marker_vel(iBody)
!
! CASE (DYNAMICS_COUPLED_FALLING_DEFOR)
! call read_marker_vel(iBody) ! Read pure deforming velocity here.
! call update_xbodymarker_noniner(iBody) ! Update body shape in non-inertia frame using the last deforming velocity. Added by G. Liu
! call get_marker_vel_defor(iBody) ! Deforming velocity of the present and the last time step. Added by G. Liu
! call correct_marker_vel_defor(iBody) ! Change the coordinate of the deforming velocity
! call compute_marker_vel(iBody) !

  END SELECT
ENDDO

if(boundary_motion_type(1)==BIO_DYNAMICS_COUPLED)then
  DO i=1,nPtsBodyMarker(1)
    xBodyMarker(1,i) = xBodyMarker(1,i) - dt*uBodyMarker(1,i)
    yBodyMarker(1,i) = yBodyMarker(1,i) - dt*vBodyMarker(1,i)
    zBodyMarker(1,i) = zBodyMarker(1,i) - dt*wBodyMarker(1,i)
  ENDDO ! i
end if

!
! n+1 n
! x - x
! - i - i n+1
! ------------ = v
! dt - i

!...Update location of boundary points
!...Update centroid location

DO iBody=1,nBody
  SELECT CASE (ndim)
  CASE (DIM_2D)
    !IF (.not. FSI_on) THEN
    DO i=1,nPtsBodyMarker(iBody)
      xBodyMarker(iBody,i) = xBodyMarker(iBody,i) + dt*uBodyMarker(iBody,i)
      yBodyMarker(iBody,i) = yBodyMarker(iBody,i) + dt*vBodyMarker(iBody,i)
      zBodyMarker(iBody,i) = zBodyMarker(iBody,i) + dt*wBodyMarker(iBody,i)
    ENDDO ! i
    CALL UnstrucSurfInDomain(iBody)

    xcent(iBody) = xcent(iBody) + dt*vxcent(iBody)
    ycent(iBody) = ycent(iBody) + dt*vycent(iBody)
    zcent(iBody) = zcent(iBody)
    !updates the angle of attack for rotating bodies.
    alpha(iBody) = alpha(iBody) + dt*angvz(iBody)*180.0_CGREAL/PI
    cosalpha(iBody) = COS(alpha(iBody)*PI/180.0_CGREAL)
    sinalpha(iBody) = SIN(alpha(iBody)*PI/180.0_CGREAL)

  CASE (DIM_3D)
    DO i=1,nPtsBodyMarker(iBody)
      xBodyMarker(iBody,i) = xBodyMarker(iBody,i) + dt*uBodyMarker(iBody,i)
      yBodyMarker(iBody,i) = yBodyMarker(iBody,i) + dt*vBodyMarker(iBody,i)
      zBodyMarker(iBody,i) = zBodyMarker(iBody,i) + dt*wBodyMarker(iBody,i)
    ENDDO ! i
    CALL UnstrucSurfInDomain(iBody)

    xcent(iBody) = xcent(iBody) + dt*vxcent(iBody)
    ycent(iBody) = ycent(iBody) + dt*vycent(iBody)
    zcent(iBody) = zcent(iBody) + dt*vzcent(iBody)
    !updates the angle of attack for rotating bodies.
    alpha(iBody) = alpha(iBody) + dt*angvz(iBody)*180.0_CGREAL/PI
    cosalpha(iBody) = COS(alpha(iBody)*PI/180.0_CGREAL)
    sinalpha(iBody) = SIN(alpha(iBody)*PI/180.0_CGREAL)
  END SELECT ! canonical_body_type
ENDDO ! iBody

IF (Prsb_MomentRef) THEN
  READ(ifuPrsbMomentRef,*) time_momentref, Moment_refx, Moment_refy, Moment_refz
ENDIF

IF (MOD(ntime,nmonitor) == 0) THEN
  ifort = 264
  DO iBody = 1,nBody
    ifortCent = ifort+iBody
    IF (boundary_motion_type(iBody) == DYNAMICS_COUPLED.OR. boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT .or. &
      boundary_motion_type(iBody) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
      boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
      boundary_motion_type(iBody) == DYNAMICS_COUPLED_SWIMMING ) THEN
      write(ifortCent,123)ntime,time,xcent(iBody),ycent(iBody),zcent(iBody),&
      vxcent(iBody),vycent(iBody),vzcent(iBody),&
      angvx(iBody),angvy(iBody),angvz(iBody)
    ENDIF
    IF (boundary_motion_type(iBody) == PARTIAL_DYNAMICS_COUPLED) THEN
      IF (ndim == 2) THEN
        TE_marker = nPtsBodyMarker(iBody)/3
        LE_marker = 1
      ENDIF
      iMarker = TE_marker
      WRITE (ifortCent,124) time, xBodyMarker(iBody,iMarker), &
            yBodyMarker(iBody,iMarker), zBodyMarker(iBody,iMarker), &
             -atan(xBodyMarker(iBody,TE_marker)-xBodyMarker(iBody,LE_marker) / &
             yBodyMarker(iBody,TE_marker)-yBodyMarker(iBody,LE_marker) )
    ENDIF
  END DO ! iBody
END IF ! ntime

123 FORMAT(i6,10(2x,e14.7))
124 FORMAT(5(2x,e14.7))

  CALL calculate_arclength_norm_ds()
  CALL set_boundary()

END SUBROUTINE move_boundary
!
!---------------------------------------------------------------------------
SUBROUTINE compute_marker_vel(iBody)
!
! this is a second-order accurate algorithm for motion developed by
! R. Mittal that preserves length during rotation.

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface
USE body_dynamics

IMPLICIT NONE

INTEGER,INTENT(IN) :: iBody

INTEGER :: i,m,j,id,iSection
REAL(KIND=CGREAL) :: uB,vB,wB
REAL(KIND=CGREAL) :: temp_uB, temp_vB, temp_wB, temp_angvx, temp_angvy, temp_angvz
REAL(KIND=CGREAL) :: total_omega

!
!
! n+1 n+1 n+1/2
! v = v + v
! - i - c -P
!
!where
!
! n+1/2 -1[ n+1/2 ( n n ) ]
! v = A [ omega X ( x - r ) ]
! -P [ ---- ( - -c ) ]
!
! where
! 2 2 2 2
! [ 1+ox dt /4 -ozdt/2 - oxoydt /4 -oydt/2 + oxozdt /4 ]
! -1 1 | 2 2 2 2 |
! A = ---------------------|ozdt/2 - oxoydt /4 1+oy dt /4 -oxdt/2 - oyozdt /4 |
! 2 2 2 2 | 2 2 2 2 |
! 1+dt (ox + oy + oz ) |oydt/2 + oxozdt /4 oxdt/2 - oyozdt /4 1+oz dt /4 |
!
! and
! n n+1 n n+1
! [ 0 -(oz + oz )/2 (oy + oy )/2 ]
! n+1/2 | n n+1 n n+1 |
! omega = |(oz + oz )/2 0 -(ox + ox )/2 |
! | n n+1 n n+1 |
! [ -(oy + oy )/2 (ox + ox )/2 0 ]

if(boundary_motion_type(iBody)==BIO_FOLLOWED_DYNAMICS_COUPLED)then
temp_angvx = half*(angvx_old(1)+angvx(1))
temp_angvy = half*(angvy_old(1)+angvy(1))
temp_angvz = half*(angvz_old(1)+angvz(1))
total_omega = oned/(oned + &
0.25_CGREAL*dt**2*(angvx(1)**2 + angvy(1)**2 + angvz(1)**2))
DO i=1,nPtsBodyMarker(iBody)

temp_uB = ( temp_angvy*(zBodyMarker(iBody,i)-zcent(1)) &
- temp_angvz*(yBodyMarker(iBody,i)-ycent(1)) )

temp_vB = -( temp_angvx*(zBodyMarker(iBody,i)-zcent(1)) &
- temp_angvz*(xBodyMarker(iBody,i)-xcent(1)) )

temp_wB = ( temp_angvx*(yBodyMarker(iBody,i)-ycent(1)) &
- temp_angvy*(xBodyMarker(iBody,i)-xcent(1)) )

uB = vxcent(1) &
+ total_omega*( temp_uB*(oned + 0.25_CGREAL*dt**2*angvx(1)**2) &
+ temp_vB*( -half*dt*angvz(1) - 0.25_CGREAL*dt**2*angvx(1)*angvy(1)) &
+ temp_wB*( -half*dt*angvy(1) + 0.25_CGREAL*dt**2*angvx(1)*angvz(1)))

vB = vycent(1) &
+ total_omega*( temp_uB*(half*dt*angvz(1) - 0.25_CGREAL*dt**2*angvx(1)*angvy(1)) &
+ temp_vB*(oned + 0.25_CGREAL*dt**2*angvy(1)**2) &
+ temp_wB*( -half*dt*angvx(1) - 0.25_CGREAL*dt**2*angvy(1)*angvz(1)))

wB = vzcent(1) &
+ total_omega*( temp_uB*(half*dt*angvy(1) + 0.25_CGREAL*dt**2*angvx(1)*angvz(1)) &
+ temp_vB*( half*dt*angvx(1) - 0.25_CGREAL*dt**2*angvy(1)*angvz(1)) &
+ temp_wB*( oned + 0.25_CGREAL*dt**2*angvz(1)**2))

uBodyMarker(iBody,i) = uBodyMarkerDyn(iBody,i) + uB
vBodyMarker(iBody,i) = vBodyMarkerDyn(iBody,i) + vB
wBodyMarker(iBody,i) = wBodyMarkerDyn(iBody,i) + wB

ENDDO ! end loop of i

!else if(boundary_motion_type(iBody)==DYNAMICS_COUPLED_FALLING_DEFOR)THEN ! Added by G. Liu
!
! temp_angvx = half*(angvx_old(iBody)+angvx(iBody))
! temp_angvy = half*(angvy_old(iBody)+angvy(iBody))
! temp_angvz = half*(angvz_old(iBody)+angvz(iBody))
! total_omega = oned/(oned + &
! 0.25_CGREAL*dt**2*(angvx(iBody)**2 + angvy(iBody)**2 + angvz(iBody)**2))
!
! DO i=1,nPtsBodyMarker(iBody)
!
! temp_uB = ( temp_angvy*(zBodyMarker(iBody,i)-zcent(iBody)) &
! - temp_angvz*(yBodyMarker(iBody,i)-ycent(iBody)) )
!
! temp_vB = -( temp_angvx*(zBodyMarker(iBody,i)-zcent(iBody)) &
! - temp_angvz*(xBodyMarker(iBody,i)-xcent(iBody)) )
!
! temp_wB = ( temp_angvx*(yBodyMarker(iBody,i)-ycent(iBody)) &
! - temp_angvy*(xBodyMarker(iBody,i)-xcent(iBody)) )
!
! uB = vxcent(iBody) &
! + total_omega*( temp_uB*(oned + 0.25_CGREAL*dt**2*angvx(iBody)**2) &
! + temp_vB*( -half*dt*angvz(iBody) - 0.25_CGREAL*dt**2*angvx(iBody)*angvy(iBody)) &
! + temp_wB*( -half*dt*angvy(iBody) + 0.25_CGREAL*dt**2*angvx(iBody)*angvz(iBody)))
!
! vB = vycent(iBody) &
! + total_omega*( temp_uB*(half*dt*angvz(iBody) - 0.25_CGREAL*dt**2*angvx(iBody)*angvy(iBody)) &
! + temp_vB*(oned + 0.25_CGREAL*dt**2*angvy(iBody)**2) &
! + temp_wB*( -half*dt*angvx(iBody) - 0.25_CGREAL*dt**2*angvy(iBody)*angvz(iBody)))
!
! wB = vzcent(iBody) &
! + total_omega*( temp_uB*(half*dt*angvy(iBody) + 0.25_CGREAL*dt**2*angvx(iBody)*angvz(iBody)) &
! + temp_vB*( half*dt*angvx(iBody) - 0.25_CGREAL*dt**2*angvy(iBody)*angvz(iBody)) &
! + temp_wB*( oned + 0.25_CGREAL*dt**2*angvz(iBody)**2))
!
! uBodyMarker(iBody,i) = uB + uBodyMarkerDefor(iBody,i)
! vBodyMarker(iBody,i) = vB + vBodyMarkerDefor(iBody,i)
! wBodyMarker(iBody,i) = wB + wBodyMarkerDefor(iBody,i)

else

temp_angvx = half*(angvx_old(iBody)+angvx(iBody))
temp_angvy = half*(angvy_old(iBody)+angvy(iBody))
temp_angvz = half*(angvz_old(iBody)+angvz(iBody))
total_omega = oned/(oned + &
0.25_CGREAL*dt**2*(angvx(iBody)**2 + angvy(iBody)**2 + angvz(iBody)**2))

DO i=1,nPtsBodyMarker(iBody)
! DO id=1,nPtsBodyMarker(iBody)-nPrescribedMarker(iBody) !Changed by Wanh for partially dynamic coupling.
! Do iSection = 2, nSection ! section 1 is prescribed in default.
! i = DynamicMarker(iBody,iSection,id)
! Note x/y/zcent are based on the dynamic parts of the body.

temp_uB = ( temp_angvy*(zBodyMarker(iBody,i)-zcent(iBody)) &
- temp_angvz*(yBodyMarker(iBody,i)-ycent(iBody)) )

temp_vB = -( temp_angvx*(zBodyMarker(iBody,i)-zcent(iBody)) &
- temp_angvz*(xBodyMarker(iBody,i)-xcent(iBody)) )

temp_wB = ( temp_angvx*(yBodyMarker(iBody,i)-ycent(iBody)) &
- temp_angvy*(xBodyMarker(iBody,i)-xcent(iBody)) )

uB = vxcent(iBody) &
+ total_omega*( temp_uB*(oned + 0.25_CGREAL*dt**2*angvx(iBody)**2) &
+ temp_vB*( -half*dt*angvz(iBody) - 0.25_CGREAL*dt**2*angvx(iBody)*angvy(iBody)) &
+ temp_wB*( -half*dt*angvy(iBody) + 0.25_CGREAL*dt**2*angvx(iBody)*angvz(iBody)))

vB = vycent(iBody) &
+ total_omega*( temp_uB*(half*dt*angvz(iBody) - 0.25_CGREAL*dt**2*angvx(iBody)*angvy(iBody)) &
+ temp_vB*(oned + 0.25_CGREAL*dt**2*angvy(iBody)**2) &
+ temp_wB*( -half*dt*angvx(iBody) - 0.25_CGREAL*dt**2*angvy(iBody)*angvz(iBody)))

wB = vzcent(iBody) &
+ total_omega*( temp_uB*(half*dt*angvy(iBody) + 0.25_CGREAL*dt**2*angvx(iBody)*angvz(iBody)) &
+ temp_vB*( half*dt*angvx(iBody) - 0.25_CGREAL*dt**2*angvy(iBody)*angvz(iBody)) &
+ temp_wB*( oned + 0.25_CGREAL*dt**2*angvz(iBody)**2))

uBodyMarker(iBody,i) = uB
vBodyMarker(iBody,i) = vB
wBodyMarker(iBody,i) = wB
! ENDDO ! end loop of iSection

ENDDO ! end loop of i
end if

END SUBROUTINE compute_marker_vel
!---------------------------------------------------------------------------
SUBROUTINE read_marker_vel(iBody)

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays

IMPLICIT NONE

INTEGER,INTENT(IN) :: iBody

INTEGER :: m,ifort,ifortCent,nptsbm_temp
INTEGER :: i, ntime_check, nptsbm_skip ! added by Haibo
INTEGER :: istat1, istat2
REAL(KIND=CGREAL) :: temp_time,temp_dt

REAL(KIND=CGREAL) :: time_skip, dt_skip
REAL(KIND=CGREAL) :: uBodyMarker_skip,vBodyMarker_skip,wBodyMarker_skip

ifort = 40
ifortCent = ifort+iBody

ntime_check = mod(ntime-1, ntimePerCycle(iBody))
IF (ntime > 1 .AND. ntime_check == 0) THEN
close (ifortCent)
END IF

IF (ntime_skip_check == 1) THEN
ntime_skip = mod(ntime_start, ntimePerCycle(iBody))
END IF

IF (nread == 1 .AND. ntime_skip /= 0) THEN
close (ifortCent)
IF (Fort_Formatted(iBody)) THEN
DO i = 1, ntime_skip
READ(ifortCent,*) dt_skip,time_skip,nptsbm_skip
DO m=1,nPtsBodyMarker(iBody)
READ(ifortCent,*)uBodyMarker_skip,vBodyMarker_skip,wBodyMarker_skip
END DO ! m
END DO
ELSE
DO i = 1, ntime_skip
READ(ifortCent) dt_skip,time_skip,nptsbm_skip
DO m=1,nPtsBodyMarker(iBody)
READ(ifortCent)uBodyMarker_skip,vBodyMarker_skip,wBodyMarker_skip
END DO ! m
END DO
ENDIF
write (*,'(A,I2,A,I3,A,I4,A)') 'Velocity file for body',iBody,' of', nBody,' skipped',ntime_skip,' steps'

IF (iBody == nBody) THEN
ntime_skip = 0
ntime_skip_check = 0
END IF
END IF
!----------------------------------------------------
! print *, 'test in read_marker_vel'
IF (Fort_Formatted(iBody)) THEN
READ(ifortCent,*)temp_dt,temp_time,nptsbm_temp
ELSE
READ(ifortCent)temp_dt,temp_time,nptsbm_temp
ENDIF
PRINT*,temp_dt,temp_time,nptsbm_temp
PRINT*,time

IF (ABS(dt-temp_dt) > 1.0E-8) THEN
PRINT*,'DT in code and marker_vel files do not match'
PRINT*,dt,temp_dt
PRINT*,'Aborting Run'
STOP
ENDIF

IF (nptsbm_temp /= nPtsBodyMarker(iBody)) THEN
PRINT*,'nPtsBodyMarker in marker_vel file does not match code'
PRINT*,nPtsBodyMarker(iBody),nptsbm_temp
PRINT*,'Aborting Run'
STOP
ENDIF

IF (Fort_Formatted(iBody)) THEN
DO m=1,nPtsBodyMarker(iBody)
READ(ifortCent,*)uBodyMarker(iBody,m),vBodyMarker(iBody,m),wBodyMarker(iBody,m)
END DO ! m
ELSE
DO m=1,nPtsBodyMarker(iBody)
READ(ifortCent)uBodyMarker(iBody,m),vBodyMarker(iBody,m),wBodyMarker(iBody,m)
END DO ! m
ENDIF
END SUBROUTINE read_marker_vel
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE correct_marker_vel(iBody)

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE operation

IMPLICIT NONE

INTEGER,INTENT(IN) :: iBody
integer :: i
real(kind=cgreal) :: ang1,ang2,ang3
real(kind=cgreal) :: quat(4),qt1(4),qt2(4)
real(kind=cgreal) :: vecInit1(3),vecInit2(3),vec1(3),vec2(3),veltoRot(3),norm(3),normInit(3),axis1(3),axis2(3)

vecInit1(1)=xBodyMarkerInit(1,rigidRef3(1))-xBodyMarkerInit(1,rigidRef1(1))
vecInit1(2)=yBodyMarkerInit(1,rigidRef3(1))-yBodyMarkerInit(1,rigidRef1(1))
vecInit1(3)=zBodyMarkerInit(1,rigidRef3(1))-zBodyMarkerInit(1,rigidRef1(1))

vecInit2(1)=xBodyMarkerInit(1,rigidRef2(1))-xBodyMarkerInit(1,rigidRef1(1))
vecInit2(2)=yBodyMarkerInit(1,rigidRef2(1))-yBodyMarkerInit(1,rigidRef1(1))
vecInit2(3)=zBodyMarkerInit(1,rigidRef2(1))-zBodyMarkerInit(1,rigidRef1(1))

normInit=cross(vecInit2,vecInit1)

vec1(1)=xBodyMarker(1,rigidRef3(1))-xBodyMarker(1,rigidRef1(1))
vec1(2)=yBodyMarker(1,rigidRef3(1))-yBodyMarker(1,rigidRef1(1))
vec1(3)=zBodyMarker(1,rigidRef3(1))-zBodyMarker(1,rigidRef1(1))

vec2(1)=xBodyMarker(1,rigidRef2(1))-xBodyMarker(1,rigidRef1(1))
vec2(2)=yBodyMarker(1,rigidRef2(1))-yBodyMarker(1,rigidRef1(1))
vec2(3)=zBodyMarker(1,rigidRef2(1))-zBodyMarker(1,rigidRef1(1))

norm=cross(vec2,vec1)

if(ndim==dim_2d)then
if(vecInit1(2)>=zero)then
ang1=acos(vecInit1(1)/mo(vecInit1,3))
else
ang1=-acos(vecInit1(1)/mo(vecInit1,3))
end if

if(vec1(2)>=zero)then
ang2=acos(vec1(1)/mo(vec1,3))
else
ang2=-acos(vec1(1)/mo(vec1,3))
end if

call quaternion_form((/0.0d0,0.0d0,1.0d0/),-ang1,qt1)
call quaternion_form((/0.0d0,0.0d0,1.0d0/),ang2,qt2)
quat=qm(qt1,qt2)
do i=1,nPtsBodyMarker(iBody)
velToRot(1)=uBodyMarker(iBody,i)
velToRot(2)=vBodyMarker(iBody,i)
velToRot(3)=wBodyMarker(iBody,i)
call quaternion_rotation(velToRot,quat)
uBodyMarkerDyn(iBody,i)=velToRot(1)
vBodyMarkerDyn(iBody,i)=velToRot(2)
wBodyMarkerDyn(iBody,i)=velToRot(3)
end do

else if(ndim==dim_3d)then
if(abs(vecInit2(1)-vec2(1))>zero.or.abs(vecInit2(2)-vec2(2))>zero.or.abs(vecInit2(3)-vec2(3))>zero.or.&
abs(vecInit1(1)-vec1(1))>zero.or.abs(vecInit1(2)-vec1(2))>zero.or.abs(vecInit1(3)-vec1(3))>zero)then
ang1=vector_angle(vecInit2,vec2)
axis1=cross(vecInit2,vec2)/mo(cross(vecInit2,vec2),3)
call quaternion_form(axis1,ang1,qt1)
call quaternion_rotation(normInit,qt1)
ang2=vector_angle(normInit,norm)
axis2=cross(normInit,norm)/mo(cross(normInit,norm),3)
call quaternion_form(axis2,ang2,qt2)
quat=qm(qt1,qt2)
do i=1,nPtsBodyMarker(iBody)
velToRot(1)=uBodyMarker(iBody,i)
velToRot(2)=vBodyMarker(iBody,i)
velToRot(3)=wBodyMarker(iBody,i)
call quaternion_rotation(velToRot,quat)
uBodyMarkerDyn(iBody,i)=velToRot(1)
vBodyMarkerDyn(iBody,i)=velToRot(2)
wBodyMarkerDyn(iBody,i)=velToRot(3)
end do
else
uBodyMarkerDyn=uBodyMarker
vBodyMarkerDyn=vBodyMarker
wBodyMarkerDyn=wBodyMarker
end if

end if

end subroutine correct_marker_vel
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE read_marker_vel_check(iBody, ntime_check)

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays

IMPLICIT NONE

INTEGER,INTENT(IN) :: iBody
INTEGER,INTENT(OUT) :: ntime_check

INTEGER :: m,ifort,ifortCent,nptsbm_temp
REAL(KIND=CGREAL) :: temp_time,temp_dt
REAL(KIND=CGREAL) :: uBodyMarker_check,vBodyMarker_check,wBodyMarker_check

INTEGER :: istat,itmp

ifort = 40
ifortCent = ifort+iBody
IF (ifortCent > 80) THEN
PRINT*,'Too many bodies. Having problems with file numbers'
STOP
ENDIF

ntime_check = 0
CLOSE(ifortCent)
Fort_Formatted(iBody)=.FALSE.
READ(ifortCent, IOSTAT=istat)temp_dt,temp_time,nptsbm_temp
IF (istat >0 .or. abs(temp_dt-dt)>1E-6 ) THEN ! NEED TO REOPEN IN FORMATTED
Fort_Formatted(iBody)=.TRUE.
WRITE(*,'(A,I2,A)') 'OPEN FORT.',ifortCent,' IN FORMATTED'
ELSE
! write(*,*) istat, temp_dt,temp_time,nptsbm_temp
WRITE(*,'(A,I2,A)') 'OPEN FORT.',ifortCent,' IN UNFORMATTED'
ENDIF
CLOSE(ifortCent)
DO

IF (Fort_Formatted(iBody)) THEN
READ(ifortCent, *, IOSTAT=istat)temp_dt,temp_time,nptsbm_temp
IF (istat >0) THEN
PRINT *, 'ERROR when read_marker_vel_check in formatted'
STOP
ENDIF
IF (istat <0) EXIT
IF ( boundary_motion_type(iBody) == PARTIAL_DYNAMICS_COUPLED .OR. &
boundary_motion_type(iBody) == FEA_FLOW_STRUC_INTERACTION) THEN
DO m=1,nptsbm_temp
READ(ifortCent,*)itmp,uBodyMarker_check,vBodyMarker_check,wBodyMarker_check
END DO ! m
ELSE IF ( boundary_motion_type(iBody) == PRESCRIBED) THEN
DO m=1,nptsbm_temp
READ(ifortCent,*)uBodyMarker_check,vBodyMarker_check,wBodyMarker_check
END DO ! m
ENDIF
ELSE ! GOOD IN UNFORMATTED
READ(ifortCent, IOSTAT=istat)temp_dt,temp_time,nptsbm_temp
IF (istat >0) THEN
PRINT *, 'ERROR when read_marker_vel_check in Unformatted'
STOP
ENDIF
IF (istat <0) EXIT
IF ( boundary_motion_type(iBody) == PARTIAL_DYNAMICS_COUPLED .OR. &
boundary_motion_type(iBody) == FEA_FLOW_STRUC_INTERACTION) THEN
DO m=1,nptsbm_temp
READ(ifortCent)itmp,uBodyMarker_check,vBodyMarker_check,wBodyMarker_check
END DO ! m
ELSE IF ( boundary_motion_type(iBody) == PRESCRIBED) THEN
DO m=1,nptsbm_temp
READ(ifortCent)uBodyMarker_check,vBodyMarker_check,wBodyMarker_check
END DO ! m
ENDIF
ENDIF
ntime_check = ntime_check + 1

END DO

WRITE(*,*) 'ntimePerPeriod for body ', iBody , ' is ', ntime_check
close(ifortCent)

END SUBROUTINE read_marker_vel_check
!------------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE read_prescribed_marker_vel(n)
! This subroutine is copied from read_marker_vel() and modified accordingly,
! for the purpose of FSI, in which only a part of markers are
! given prescribed motion.
! Unlike the previous subroutine read_marker_vel(), the node number in the
! read-in "fort" files need to be specified on structure. To avoid uncessary changes
! in subroutine read_marker_vel() and format of data entries in "fort" files,
! read_prescribed_marker_vel() is used to read in the prescribed
! motion on specified nodes. //Wanh 05/05/10
!
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface
USE body_dynamics

IMPLICIT NONE

INTEGER,INTENT(IN) :: n
INTEGER :: prescribed_point, iSection, id

INTEGER :: m,im,ifort,ifortCent
INTEGER :: i, ntime_check, nptsbm_skip ! added by Haibo

REAL(KIND=CGREAL) :: temp_time,temp_dt

REAL(KIND=CGREAL) :: time_skip, dt_skip
REAL(KIND=CGREAL) :: uBodyMarker_skip,vBodyMarker_skip,wBodyMarker_skip

ifort = 40
ifortCent = ifort+n
IF (ifortCent > 80) PRINT*,'Might have problems with file numbers'

!--------- new added by Haibo ----------------------
ntime_check = mod(ntime-1, ntimePerCycle(n))

IF (ntime_skip_check == 1) THEN
ntime_skip = mod(ntime_start, ntimePerCycle(n))
END IF

write(*,*) 'ntime_check, ntime_skip =', ntime_check, ntime_skip, ifortCent, ntime

IF (ntime > 1.AND.ntime_check == 0) THEN
close (ifortCent)
END IF

IF (nread == 1.AND. ntime_skip /= 0) THEN
close (ifortCent)
print *, 'here ...'
DO i = 1, ntime_skip
READ(ifortCent,*)dt_skip,time_skip,nPrescribedMarker(n)
DO m=1,nPrescribedMarker(n)
READ(ifortCent,*)im,uBodyMarker_skip,vBodyMarker_skip,wBodyMarker_skip
END DO ! m
END DO

IF (n == nBody) THEN
ntime_skip = 0
ntime_skip_check = 0
END IF
END IF
!----------------------------------------------------
READ(ifortCent,*)temp_dt,temp_time,nPrescribedMarker(n)

IF (ABS(dt-temp_dt) > 1.0E-8) THEN
PRINT*,'DT in code and marker_vel files do not match.'
PRINT*,dt,temp_dt
PRINT*,'Aborting Run'
STOP
ENDIF

IF (FSI_ON) THEN
DO i=1,nPtsBodyMarker(n) !for simplicity, velocity on all markers are read in.
uBodyMarker(n,i) = struc_vel(i,1)
vBodyMarker(n,i) = struc_vel(i,2)
wBodyMarker(n,i) = struc_vel(i,3)
ENDDO
ENDIF

DO m=1,nPrescribedMarker(n)
READ(ifortCent,*) PrescribedMarker(n,m),uBodyMarker(n,PrescribedMarker(n,m)), &
vBodyMarker(n,PrescribedMarker(n,m)),wBodyMarker(n,PrescribedMarker(n,m))
END DO ! m

id = 0
DO i=1,nPtsBodyMarker(n)
prescribed_point = 0
DO m=1,nPrescribedMarker(n)
IF ( i == PrescribedMarker(n,m) ) THEN
prescribed_point = 1
EXIT
ENDIF
ENDDO

DO iSection = 2,nSection
!****** Temporarily 2 sections are assigned, need to changed for nSection>2.
IF (prescribed_point == 0) then
id = id+1
DynamicMarker(n,iSection,id) = i
ENDIF
SectionMarker(n,iSection) = id

ENDDO ! loop of iSection

ENDDO


END SUBROUTINE read_prescribed_marker_vel
!------------------------------------------------------------------------------

SUBROUTINE compute_marker_vel_section(n)
!
! this is a second-order accurate algorithm for motion developed by
! R. Mittal that preserves length during rotation.

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface
USE body_dynamics

IMPLICIT NONE

INTEGER,INTENT(IN) :: n

INTEGER :: i,m,j,id,iSection
REAL(KIND=CGREAL) :: uB,vB,wB
REAL(KIND=CGREAL) :: temp_uB, temp_vB, temp_wB, temp_angvx, temp_angvy, temp_angvz
REAL(KIND=CGREAL) :: total_omega
REAL(KIND=CGREAL) :: uref, vref, wref
REAL(KIND=CGREAL) :: xref, yref, zref
REAL(KIND=CGREAL) :: usection_cent, vsection_cent,wsection_cent

WRITE (*,*) 'Calling subroutine compute_marker_vel_section ...'

DO iSection = 2, nSection

xref = xBodyMarker(n,hingemarker)
yref = yBodyMarker(n,hingemarker)
zref = zBodyMarker(n,hingemarker)

uref = uBodyMarker(n,hingemarker)
vref = vBodyMarker(n,hingemarker)
wref = wBodyMarker(n,hingemarker)

print *, ' xref,yref,zref:', xref, yref,zref
print *, ' uref,vref,wref:', uref, vref,wref
print *, ' section_xcent,section_ycent:',section_xcent(n,iSection),section_ycent(n,iSection)

temp_angvx = half*(angvx_old(iSection)+angvx(iSection))
temp_angvy = half*(angvy_old(iSection)+angvy(iSection))
temp_angvz = half*(angvz_old(iSection)+angvz(iSection))

print *, ' angvz(iSection) = ', angvz(iSection)
print *, ' temp_angvz = ', temp_angvz
print *, ' section_ycent(n,iSection) = ', section_ycent(n,iSection)

total_omega = oned/(oned + &
0.25_CGREAL*dt**2*(angvx(iSection)**2 + angvy(iSection)**2 + angvz(iSection)**2))

usection_cent = uref - temp_angvz*(section_ycent(n,iSection)-yref)
vsection_cent = vref + temp_angvz*(section_xcent(n,iSection)-xref)
wsection_cent = 0.0

DO id=1,SectionMarker(n,iSection)
i = DynamicMarker(n,iSection,id)

temp_uB = ( temp_angvy*(zBodyMarker(n,i)-section_zcent(n,iSection)) &
- temp_angvz*(yBodyMarker(n,i)-section_ycent(n,iSection)) )

temp_vB = -( temp_angvx*(zBodyMarker(n,i)-section_zcent(n,iSection)) &
- temp_angvz*(xBodyMarker(n,i)-section_xcent(n,iSection)) )

temp_wB = ( temp_angvx*(yBodyMarker(n,i)-section_ycent(n,iSection)) &
- temp_angvy*(xBodyMarker(n,i)-section_xcent(n,iSection)) )

uB = usection_cent &
+ total_omega*( temp_uB*(oned + 0.25_CGREAL*dt**2*angvx(iSection)**2) &
+ temp_vB*( -half*dt*angvz(iSection) - 0.25_CGREAL*dt**2*angvx(iSection)*angvy(iSection)) &
+ temp_wB*( -half*dt*angvy(iSection) + 0.25_CGREAL*dt**2*angvx(iSection)*angvz(iSection)))

vB = vsection_cent &
+ total_omega*( temp_uB*(half*dt*angvz(iSection) - 0.25_CGREAL*dt**2*angvx(iSection)*angvy(iSection)) &
+ temp_vB*(oned + 0.25_CGREAL*dt**2*angvy(iSection)**2) &
+ temp_wB*( -half*dt*angvx(iSection) - 0.25_CGREAL*dt**2*angvy(iSection)*angvz(iSection)))

wB = wsection_cent &
+ total_omega*( temp_uB*(half*dt*angvy(iSection) + 0.25_CGREAL*dt**2*angvx(iSection)*angvz(iSection)) &
+ temp_vB*( half*dt*angvx(iSection) - 0.25_CGREAL*dt**2*angvy(iSection)*angvz(iSection)) &
+ temp_wB*( oned + 0.25_CGREAL*dt**2*angvz(iSection)**2))

uBodyMarker(n,i) = uB
vBodyMarker(n,i) = vB
wBodyMarker(n,i) = wB
ENDDO ! end loop of iSection

print *, 'max/min ubodymarker:', maxval(uBodyMarker), minval(uBodyMarker)
print *, 'max/min vbodymarker:', maxval(vBodyMarker), minval(vBodyMarker)

ENDDO ! end loop of id

END SUBROUTINE compute_marker_vel_section
!====================================================================================================
SUBROUTINE move_boundary_defor()
! Added by G. Liu. This subroutine is used in "motion tyep = DYNAMICS_COUPLED_FALLING_DEFOR" only.
! update internal boundary (with deformation) position
! First rotate the body without deforming motion, and then add the deforming motion
!
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface
USE usr_module
USE body_dynamics

IMPLICIT NONE

character fname*25
INTEGER :: i,n,ifortCent,iBody,iFort,m,j,k,i_flag,presb_flag,iSection
real(cgreal) :: uBodyMarker_temp(nBody,nPtsMax),vBodyMarker_temp(nBody,nPtsMax),wBodyMarker_temp(nBody,nPtsMax)


do iBody=1,nBody
call compute_marker_vel(iBody) !
xBodymarker(iBody,:) = xBodyMarker(iBody,:) + dt*uBodyMarker(iBody,:)
yBodymarker(iBody,:) = yBodyMarker(iBody,:) + dt*vBodyMarker(iBody,:)
uBodyMarker_temp(iBody,:) = uBodyMarker(iBody,:)
vBodyMarker_temp(iBody,:) = vBodyMarker(iBody,:)
if(ndim == DIM_3D)then
zBodymarker(iBody,:) = zBodyMarker(iBody,:) + dt*wBodyMarker(iBody,:)
wBodyMarker_temp(iBody,:) = wBodyMarker(iBody,:)
endif
call get_quat_defor(iBody)
write(*,*)'get_quat_defor done!'
call read_marker_vel(iBody)
uBodyMarkerDefor(iBody,:) = uBodyMarker(iBody,:)
vBodyMarkerDefor(iBody,:) = vBodyMarker(iBody,:)
wBodyMarkerDefor(iBody,:) = wBodyMarker(iBody,:)
call rotate_defor_vel(quat_trans,iBody) ! uBodyMarkerDefor -> uBodyMarkerDyn (non-inertia frame -> inertia frame)
call update_bodyshape_iner(iBody)
call update_bodyshape_noniner(iBody)

uBodyMarker(iBody,:) = uBodyMarker_temp(iBody,:) + uBodyMarkerDyn(iBody,:)
vBodyMarker(iBody,:) = vBodyMarker_temp(iBody,:) + vBodyMarkerDyn(iBody,:)
wBodyMarker(iBody,:) = wBodyMarker_temp(iBody,:) + wBodyMarkerDyn(iBody,:)

enddo ! iBody

Do iBody=1, nBody
vxcent_prev(iBody) = vxcent(iBody)
vycent_prev(iBody) = vycent(iBody)
vzcent_prev(iBody) = vzcent(iBody)
angvx_old(ibody) = angvx(ibody)
angvy_old(ibody) = angvy(ibody)
angvz_old(ibody) = angvz(ibody)
ENDDO

! write(*,*)'liu debug in move_boundary_defor: stop'
! stop

CALL calculate_arclength_norm_ds()

CALL set_boundary()

endsubroutine move_boundary_defor

!=============================================================================
subroutine rotate_defor_vel(quat,iBody)
!
!
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface
USE usr_module
USE body_dynamics

real(cgreal) :: quat(4) ! input
real(cgreal) :: vec(3)
integer i,m

do i=1,nPtsBodyMarker(iBody)
vec(1) = uBodyMarkerDefor(iBody,i)
vec(2) = vBodyMarkerDefor(iBody,i)
vec(3) = wBodyMarkerDefor(iBody,i)
call quaternion_rotation(vec,quat)
uBodyMarkerDyn(iBody,i) = vec(1) ! deforming velocity in initia frame
vBodyMarkerDyn(iBody,i) = vec(2)
wBodyMarkerDyn(iBody,i) = vec(3)
enddo


endsubroutine

!==============================================================================
subroutine update_bodyshape_iner(iBody)
! Added by G. Liu
!
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
use body_dynamics
use usr_module

integer :: i

xBodyMarker(iBody,:) = xBodyMarker(iBody,:) + dt*uBodyMarkerDyn(iBody,:)
yBodyMarker(iBody,:) = yBodyMarker(iBody,:) + dt*vBodyMarkerDyn(iBody,:)
zBodyMarker(iBody,:) = zBodyMarker(iBody,:) + dt*wBodyMarkerDyn(iBody,:)

end subroutine
!=================================================================
subroutine update_bodyshape_noniner(iBody)
! if ntime == 0, then uBodyMarkerDefor_prvs = 0.0
!
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
use body_dynamics
use usr_module

integer :: i

xBodyMarkerNonIner(iBody,:) = xBodyMarkerNonIner(iBody,:) + dt*uBodyMarkerDefor(iBody,:)
yBodyMarkerNonIner(iBody,:) = yBodyMarkerNonIner(iBody,:) + dt*vBodyMarkerDefor(iBody,:)
zBodyMarkerNonIner(iBody,:) = zBodyMarkerNonIner(iBody,:) + dt*wBodyMarkerDefor(iBody,:)


end subroutine
!======================================================================================
subroutine Add_defor_vel(iBody)
!Added by G. Liu
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface
USE usr_module
USE body_dynamics

integer :: iBody

uBodyMarker(iBody,:) = uBodyMarker(iBody,:) + uBodyMarkerDyn(iBody,:)
vBodyMarker(iBody,:) = vBodyMarker(iBody,:) + vBodyMarkerDyn(iBody,:)
wBodyMarker(iBody,:) = wBodyMarker(iBody,:) + wBodyMarkerDyn(iBody,:)

endsubroutine
!=====================================================================================
