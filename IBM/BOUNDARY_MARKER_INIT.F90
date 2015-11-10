!------------------------------------------------------------------------------
SUBROUTINE initialize_marker()

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
use body_dynamics

IMPLICIT NONE


!...Loop variables
INTEGER :: i,j,iBody
INTEGER :: k, nnt

!...Local variables
INTEGER :: m, nBodyMarkerIn,nBodyRstrt,n, &
nPtsBodyMarkerIn, totNumTriElemIn
INTEGER, DIMENSION(nBody) :: body_type_orig

REAL(KIND=CGREAL) :: dMin,dxMin,dyMin,dzMin, &
bodyResolution, bodyResolutionNormalized, &
rad,theta,phi, xTemp, yTemp
LOGICAL :: readMarkerFlag

! Set read Marker flag

PRINT*,'SETTING UP CANONICAL BODIES in INITIALIZE_MARKER '

IF ( nread == 1) GOTO 1000

! Save copy of canonical body type

body_type_orig(1:nBody) = canonical_body_type(1:nBody)

! Initialize for no restart

uBodyMarker = zero
vBodyMarker = zero
wBodyMarker = zero

xBodyMarker = zero
yBodyMarker = zero
zBodyMarker = zero

! Select appropriate body type
! Setting up marker locations

OPEN(UNIT=145,file='marker_unstruc_tri.dat',STATUS='UNKNOWN')
IF(zoneSeparate) then
OPEN(UNIT=245,file='marker_unstruc_tri_zoneColor.dat',STATUS='UNKNOWN') !added by Chengyu for checking the zone numbers of the unstruc_surface file
ENDIF

DO iBody = 1, nBody

SELECT CASE (canonical_body_type(iBody))

CASE(ELLIPTIC_CYLINDER)
PRINT*,' SETTING UP ELLIPTIC CYLINDER'

!-- Test resolution for GCM
IF (boundary_formulation == GCM_METHOD) THEN
dxMin = MINVAL(dx(1:nxc))
dyMin = MINVAL(dy(1:nyc))
dMin = MIN(dxMin,dyMin)
bodyResolution = PI*( radiusx(iBody)+radiusy(iBody) ) / &
REAL(nPtsBodyMarker(iBody),KIND=CGREAL)
bodyResolutionNormalized = dMin/bodyResolution

PRINT*,' dx Min ',dxMin
PRINT*,' dy Min ',dyMin
PRINT*,' d Min ',dMin
PRINT*,' Body Resolution',bodyResolution
PRINT*,' Current Normalized Resolution for Body ',bodyResolutionNormalized
IF (bodyResolutionNormalized < 2.0_CGREAL ) THEN
PRINT*,' Ideal Normalized Resolution Should be at LEAST 2 .. aborting'
STOP
ENDIF
ENDIF

DO m = 1, nPtsBodyMarkerOrig(iBody)
theta = (REAL(m,KIND=CGREAL)-oned)*2.0_CGREAL*PI/ &
REAL(nPtsBodyMarkerOrig(iBody),KIND=CGREAL) 
xTemp = radiusx(iBody)*COS(theta)
yTemp = radiusy(iBody)*SIN(theta)
xBodyMarker(iBody,m) = xcent(iBody) + xTemp*cosalpha(iBody) - yTemp*sinalpha(iBody) 
yBodyMarker(iBody,m) = ycent(iBody) + xTemp*sinalpha(iBody) + yTemp*cosalpha(iBody) 
zBodyMarker(iBody,m) = z(1)
ENDDO ! m

CALL extend_cylinder_3D(iBody)

CASE(GENERAL_CYLINDER)
PRINT*,' SETTING UP GENERAL CYLINDER'
READ(ifuMarkerIn,*) nBodyMarkerIn
IF ( nBodyMarkerIn /= nPtsBodyMarkerOrig(iBody) ) THEN
PRINT*,'Init_Marker: Inconsistent body_in.dat and marker_in.dat files for body = ', iBody
PRINT*,' Reading in body_in.dat nPtsBodyMarker = ',nPtsBodyMarkerOrig(iBody)
PRINT*,' Reading from marker_in.dat nBodyMarkerIn = ', nBodyMarkerIn
CALL abort_vicar3d(10)
ENDIF ! nBodyMarkerIn

PRINT*,' iBody ', iBody
PRINT*,' nBodyMarkerIn ', nBodyMarkerIn
PRINT*,' nPtsBodyMarker ', nPtsBodyMarker(iBody)
DO m = 1, nPtsBodyMarkerOrig(iBody)
READ(ifuMarkerIn,*)xBodyMarker(iBody,m),yBodyMarker(iBody,m)
zBodyMarker(iBody,m) = z(1)
ENDDO

CALL extend_cylinder_3D(iBody)

CASE(ELLIPSOID)
PRINT*,' SETTING UP ELLIPSOID'

m = 0

DO i = 1, n_phi(iBody)
phi = ( REAL(i,KIND=CGREAL)-oned )*PI/ &
(REAL(n_phi(iBody),KIND=CGREAL)-oned)
IF (i .EQ. 1 .OR. i .EQ. n_phi(iBody)) THEN
theta = zero
m = m + 1

xTemp = radiusx(iBody)*COS(theta)*SIN(phi)
yTemp = radiusy(iBody)*SIN(theta)*SIN(phi)
xBodyMarker(iBody,m) = xcent(iBody) + xTemp*cosalpha(iBody) - yTemp*sinalpha(iBody) 
yBodyMarker(iBody,m) = ycent(iBody) + xTemp*sinalpha(iBody) + yTemp*cosalpha(iBody)
zBodyMarker(iBody,m) = zcent(iBody) + radiusz(iBody)*COS(phi)
ELSE
DO j = 1,n_theta(iBody)
theta = (REAL(j,KIND=CGREAL)-oned)*2.0_CGREAL*PI/ &
REAL(n_theta(iBody),KIND=CGREAL)
m = m + 1

xTemp = radiusx(iBody)*COS(theta)*SIN(phi)
yTemp = radiusy(iBody)*SIN(theta)*SIN(phi)
xBodyMarker(iBody,m) = xcent(iBody) + xTemp*cosalpha(iBody) - yTemp*sinalpha(iBody)
yBodyMarker(iBody,m) = ycent(iBody) + xTemp*sinalpha(iBody) + yTemp*cosalpha(iBody)
zBodyMarker(iBody,m) = zcent(iBody) + radiusz(iBody)*COS(phi)
ENDDO ! j
END IF
ENDDO ! i

nPtsBodyMarker(iBody) = m
write(*,*) 'Total number of marker points are: ', m
j = 0
i = 1

DO k = 1, n_theta(iBody)
j = j + 1
triElemNeig(iBody,1,j) = i
triElemNeig(iBody,2,j) = i + k
IF (i+k .EQ. n_theta(ibody) + 1) THEN
triElemNeig(iBody,3,j) = i + 1
ELSE 
triElemNeig(iBody,3,j) = i + k + 1
END IF
END DO ! end k 

nnt = 1

DO i = 2, m-n_theta(iBody)-1

j = j + 1
triElemNeig(iBody,1,j) = i
triElemNeig(iBody,2,j) = i + n_theta(iBody)
! write (*,*) 'i, n_theta =', i, nnt*n_theta(ibody)+1

IF ( i .EQ. nnt*n_theta(ibody)+1 ) THEN
triElemNeig(iBody,3,j) = i + 1 
ELSE
triElemNeig(iBody,3,j) = i + n_theta(iBody) + 1
END IF

j = j + 1 
triElemNeig(iBody,1,j) = i
IF ( i .EQ. nnt*n_theta(ibody) + 1 ) THEN
triElemNeig(iBody,2,j) = i + 1
triElemNeig(iBody,3,j) = i + 1 - n_theta(iBody)
nnt = nnt + 1 
ELSE
triElemNeig(iBody,2,j) = i + n_theta(iBody) + 1
triElemNeig(iBody,3,j) = i + 1
END IF

ENDDO ! i

DO i = m-n_theta(iBody), m-1

j = j + 1
triElemNeig(iBody,1,j) = i
triElemNeig(iBody,2,j) = m
IF (i+1 .EQ. m) THEN
triElemNeig(iBody,3,j) = m-n_theta(iBody) 
ELSE
triElemNeig(iBody,3,j) = i + 1
END IF
END DO ! end i

totNumTriElem(iBody) = j

write(*,*) 'Total # of elements: ', j

! END IF

CASE(UNSTRUCTURED_SURFACE)
PRINT*,' SETTING UP UNSTRUCTURED SURFACE'

IF ( boundary_motion /= MOVING_BOUNDARY .OR. nread /= 1) THEN
READ(ifuUnstrucSurfIn,*)
READ(ifuUnstrucSurfIn,*) nPtsBodyMarkerIn, totNumTriElemIn
print *, 'nPtsBodyMarkerIn, totNumTriElemIn in BOUNDARY_MARKER_INIT is:', nPtsBodyMarkerIn, totNumTriElemIn
READ(ifuUnstrucSurfIn,*)

IF ( nPtsBodyMarkerIn /= nPtsBodyMarker(iBody) ) THEN
PRINT*,'Init_Marker: Inconsistent canonical_body_in.dat and unstruc_surface_in.dat files for body = ', iBody
PRINT*,' Reading in canonical_body_in.dat nPtsBodyMarker = ',nPtsBodyMarker(iBody)
PRINT*,' Reading from unstruc_surface_in.dat nPtsBodyMarkerIn = ',nPtsBodyMarkerIn
CALL abort_vicar3d(10)
ENDIF ! nPtsBodyMarkerIn

IF ( totNumTriElemIn /= totNumTriElem(iBody) ) THEN
PRINT*,'Init_Marker: Inconsistent canonical_body_in.dat and unstruc_surface_in.dat files for body = ', iBody
PRINT*,' Reading in canonical_body_in.dat totNumTriElem = ', totNumTriElem(iBody)
PRINT*,' Reading from unstruc_surface_in.dat totNumTriElemIn = ', totNumTriElemIn
CALL abort_vicar3d(10)
ENDIF ! totNumTriElemIn

DO m=1,nPtsBodyMarker(iBody)
if(.not.channel_flow .AND. .not.zoneSeparate)then
READ(ifuUnstrucSurfIn,*) i,xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m)
elseif(zoneSeparate) then
READ(ifuUnstrucSurfIn,*) i,xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m),zoneMarker(iBody,m) !add zoneMarker array, added by Chengyu 2/12/14
elseif(channel_flow) then
READ(ifuUnstrucSurfIn,*) i,xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m),gateLabel(iBody,m) !add by Yan for channel flow BC
end if
ENDDO
zoneMax(iBody)=maxval(zoneMarker) !add zoneMax array, added by Chengyu 2/12/14

READ(ifuUnstrucSurfIn,*)
DO j=1,totNumTriElem(iBody)
READ(ifuUnstrucSurfIn,*) i,triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
ENDDO
READ(ifuUnstrucSurfIn,*)
READ(ifuUnstrucSurfIn,*)pointOutsideBodyX(iBody),pointOutsideBodyY(iBody),pointOutsideBodyZ(iBody)
ENDIF

if(boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_SWIMMING )then ! Added by G. Liu
xBodyMarkerNonIner(iBody,:)=xBodyMarker(iBody,:)
yBodyMarkerNonIner(iBody,:)=yBodyMarker(iBody,:)
zBodyMarkerNonIner(iBody,:)=zBodyMarker(iBody,:)
call Marker_move_rot(iBody,xcent(ibody),ycent(ibody),zcent(ibody),quat_init) ! translation and rotation of the marker points. Non-inertial -> inertial
quat_iter=quat_init
quat_phys=quat_init
endif

! IF( boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR )THEN ! Added by G. Liu (deforming velocity translation)
! call Defor_vel_rot(iBody,quat_init)
! ENDIF

CALL WriteSurfInTecplot(iBody)

CALL UnstrucSurfInDomain(iBody)

xBodyMarkerInit(iBody,:)=xBodyMarker(iBody,:)
yBodyMarkerInit(iBody,:)=yBodyMarker(iBody,:)
zBodyMarkerInit(iBody,:)=zBodyMarker(iBody,:)

END SELECT ! canonical_body_type

ENDDO ! iBody

! set up initial marker velocities

DO iBody = 1, nBody

SELECT CASE (canonical_body_type(iBody))

CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)

SELECT CASE (boundary_motion_type(iBody))

CASE (STATIONARY)
SELECT CASE (wall_type(iBody))
CASE(POROUS_OR_SLIP)
CALL wall_velocity(iBody)
CASE(NONPOROUS_AND_NONSLIP)
DO m = 1, nPtsBodyMarkerOrig(iBody)
uBodyMarker(iBody,m) = zero
vBodyMarker(iBody,m) = zero
wBodyMarker(iBody,m) = zero
ENDDO
END SELECT ! wall_type

CASE (FORCED)
DO m = 1, nPtsBodyMarkerOrig(iBody)
uBodyMarker(iBody,m) = vxcent(iBody) &
+ ( -angvz(iBody)*(yBodyMarker(iBody,m)-ycent(iBody)) )
vBodyMarker(iBody,m) = vycent(iBody) &
- ( -angvz(iBody)*(xBodyMarker(iBody,m)-xcent(iBody)) )
wBodyMarker(iBody,m) = vzcent(iBody) 
ENDDO

CASE (FLOW_INDUCED:PRESCRIBED)
DO m = 1, nPtsBodyMarkerOrig(iBody)
uBodyMarker(iBody,m) = zero
vBodyMarker(iBody,m) = zero
wBodyMarker(iBody,m) = zero
ENDDO


END SELECT ! boundary_motion_type

CALL extend_cylinder_vel_3d(iBody)

CASE(ELLIPSOID:UNSTRUCTURED_SURFACE)

PRINT*,' SETTING UP MOVING ELLIPSOID OR UNSTRUCTURED SURFACE : SSM'

! .. set marker velocity
SELECT CASE (boundary_motion_type(iBody))

CASE (STATIONARY)
SELECT CASE (wall_type(iBody))
CASE(POROUS_OR_SLIP)
CALL wall_velocity(iBody)
CASE(NONPOROUS_AND_NONSLIP)
DO m = 1, nPtsBodyMarker(iBody)
uBodyMarker(iBody,m) = zero
vBodyMarker(iBody,m) = zero
wBodyMarker(iBody,m) = zero
ENDDO
END SELECT ! wall_type

CASE (FORCED)
DO m = 1, nPtsBodyMarker(iBody)
uBodyMarker(iBody,m) = vxcent(iBody) &
+ ( angvy(iBody)*(zBodyMarker(iBody,m)-zcent(iBody)) &
- angvz(iBody)*(yBodyMarker(iBody,m)-ycent(iBody)) )
vBodyMarker(iBody,m) = vycent(iBody) &
- ( angvx(iBody)*(zBodyMarker(iBody,m)-zcent(iBody)) &
- angvz(iBody)*(xBodyMarker(iBody,m)-xcent(iBody)) )
wBodyMarker(iBody,m) = vzcent(iBody) &
+ ( angvx(iBody)*(yBodyMarker(iBody,m)-ycent(iBody)) &
- angvy(iBody)*(xBodyMarker(iBody,m)-xcent(iBody)) )
ENDDO

CASE (FLOW_INDUCED:PRESCRIBED)
DO m = 1, nPtsBodyMarker(iBody)
uBodyMarker(iBody,m) = zero
vBodyMarker(iBody,m) = zero
wBodyMarker(iBody,m) = zero
ENDDO

CASE (BIO_FOLLOWED_DYNAMICS_COUPLED)
DO m = 1, nPtsBodyMarkerOrig(iBody)
uBodyMarker(iBody,m) = zero
vBodyMarker(iBody,m) = zero
wBodyMarker(iBody,m) = zero
ENDDO

END SELECT ! boundary_motion_type

END SELECT ! canonical_body_type

ENDDO ! iBody

WRITE(ifuBodyOut,*)nbody

! write body markers into Output file
DO iBody = 1, nBody

! writing out unstructured surface
WRITE(ifuUnstrucSurfOut,*)
WRITE(ifuUnstrucSurfOut,*) nPtsBodyMarker(iBody), totNumTriElem(iBody)
WRITE(ifuUnstrucSurfOut,*)
DO m=1,nPtsBodyMarker(iBody)
WRITE(ifuUnstrucSurfOut,*) m,xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m)
ENDDO
WRITE(ifuUnstrucSurfOut,*)
DO j=1,totNumTriElem(iBody)
WRITE(ifuUnstrucSurfOut,*) j,triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
ENDDO
WRITE(ifuUnstrucSurfOut,*)
WRITE(ifuUnstrucSurfOut,*)pointOutsideBodyX(iBody),pointOutsideBodyY(iBody),pointOutsideBodyZ(iBody)

! writing out body parameter file
IF(canonical_body_type(iBody) <= GENERAL_CYLINDER) THEN
WRITE(*,*) '2D BODY IS NOT SUPPORTED ANYMORE!'
! The following 5 lines are commented out by Wanh temporarily.
! STOP
! canonical_body_type(iBody) = 4
! body_dim(iBody) = 2
! radiusz(iBody) = zero
! zcent(iBody) = zero

ELSE IF (canonical_body_type(iBody) == ELLIPSOID .AND. &
boundary_formulation == GCM_METHOD) THEN
canonical_body_type(iBody) = 4
ENDIF

WRITE(ifuBodyOut,*)canonical_body_type(iBody),body_dim(iBody),boundary_motion_type(iBody)
WRITE(ifuBodyOut,*)wall_type(iBody)
WRITE(ifuBodyOut,*)nPtsBodyMarker(iBody),totNumTriElem(iBody)
WRITE(ifuBodyOut,*)radiusx(iBody),radiusy(iBody),radiusz(iBody)
WRITE(ifuBodyOut,*)xcent(iBody),ycent(iBody),zcent(iBody)
WRITE(ifuBodyOut,*)alpha(iBody)
WRITE(ifuBodyOut,*)vxcentTrans(iBody),vycentTrans(iBody),vzcentTrans(iBody)
WRITE(ifuBodyOut,*)ampx(iBody),ampy(iBody),ampz(iBody)
WRITE(ifuBodyOut,*)freqx(iBody),freqy(iBody),freqz(iBody)
WRITE(ifuBodyOut,*)angvx(iBody),angvy(iBody),angvz(iBody)
WRITE(ifuBodyOut,*)phase(iBody)
WRITE(ifuBodyOut,*)ampangx(iBody),ampangy(iBody),ampangz(iBody)
WRITE(ifuBodyOut,*)freqangx(iBody),freqangy(iBody),freqangz(iBody)
WRITE(ifuBodyOut,*)density_fluid, density_solid(iBody)
WRITE(ifuBodyOut,*)xcentConstr(iBody),ycentConstr(iBody),zcentConstr(iBody) 

ENDDO

DO iBody = 1, nBody
! IF(canonical_body_type(iBody) <= GENERAL_CYLINDER) THEN
! The following 11 lines are commented out by Wanh temporarily.
! IF(body_type_orig(iBody) <= GENERAL_CYLINDER .OR. &
! (body_type_orig(iBody) == ELLIPSOID .AND. &
! boundary_formulation == GCM_METHOD) ) THEN

! PRINT*,'######################CODE NEEDS TO BE RERUN ##########################################'
! PRINT*,'Body parametrs have been written out in : canonical_body_out.dat'
! PRINT*,'Cylinder Surface(s) have been written out in : unstructured_surface_out.dat'
! PRINT*,'Following needs to be done:'
! PRINT*,'(1) Rename unstruc_surface_out.dat TO unstruc_surface_in.dat'
! PRINT*,'(2) Rename canonical_body_out.dat TO canonical_body_in.dat'
! PRINT*,'(3) Run code again'
! STOP
! ENDIF
ENDDO

1000 CONTINUE

OPEN(ifuMarkerOut,FILE='marker_out.dat',STATUS='UNKNOWN')
DO iBody = 1, nBody
WRITE(ifuMarkerOut,*) 'Body:',iBody
DO m = 1, nPtsBodyMarkerOrig(iBody)
WRITE(ifuMarkerOut,*) xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m) 
ENDDO
ENDDO

PRINT*, 'CALL calculate_arclength_norm_ds() in initialize_marker'
CALL calculate_arclength_norm_ds()

CLOSE(145)

END SUBROUTINE initialize_marker 
!-----------------------------------------------------------------

! extend cylinder to make pseudo-3D body

SUBROUTINE extend_cylinder_3D(iBody)

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays

IMPLICIT NONE

INTEGER, INTENT(IN) :: iBody

!...Loop variables
INTEGER :: i,j,k,m,mInd,mp, &
nPtsBodyMarkerIn, totNumTriElemIn

INTEGER :: nodesMax

REAL(KIND=CGREAL) :: dist_beg_end_marker,ds_marker

! Determine if body if open or closed
dist_beg_end_marker = SQRT( (xBodyMarker(iBody,1)-xBodyMarker(iBody,nPtsBodymarkerOrig(iBody)))**2 &
+(yBodyMarker(iBody,1)-yBodyMarker(iBody,nPtsBodymarkerOrig(iBody)))**2 )
ds_marker = SQRT( (xBodyMarker(iBody,2)-xBodyMarker(iBody,1))**2 &
+(yBodyMarker(iBody,2)-yBodyMarker(iBody,1))**2 ) 
IF (dist_beg_end_marker > 2.0_CGREAL*ds_marker) THEN
PRINT*,'Body seems to be an "open" body'
nodesMax = nPtsBodymarkerOrig(iBody)-1
totNumTriElem(iBody) = 2*(nPtsBodyMarkerOrig(iBody)-1)*(nz-1)
ELSE
nodesMax = nPtsBodymarkerOrig(iBody)
ENDIF

! assign nodal cooridinates
DO k=2,nz
DO m=1,nPtsBodymarkerOrig(iBody)
mInd = (k-1)*nPtsBodymarkerOrig(iBody)+m
xBodyMarker(iBody,mInd) = xBodyMarker(iBody,m)
yBodyMarker(iBody,mInd) = yBodyMarker(iBody,m)
zBodyMarker(iBody,mInd) = z(k)
ENDDO
ENDDO

! create neighbor list for elements
j = 0
DO k = 1, nz-1
DO m = 1, nodesMax

j = j+1

mp = m+1
IF (m == nPtsBodyMarkerOrig(iBody)) mp = 1

triElemNeig(iBody,1,j) = (k-1)*nPtsBodymarkerOrig(iBody) + m
triElemNeig(iBody,2,j) = k*nPtsBodymarkerOrig(iBody) + m
triElemNeig(iBody,3,j) = (k-1)*nPtsBodymarkerOrig(iBody) + mp

j = j+1
triElemNeig(iBody,1,j) = (k-1)*nPtsBodymarkerOrig(iBody) + mp
triElemNeig(iBody,2,j) = k*nPtsBodymarkerOrig(iBody) + m
triElemNeig(iBody,3,j) = k*nPtsBodymarkerOrig(iBody) + mp

ENDDO
ENDDO

CALL WriteSurfInTecplot(iBody)

pointOutsideBodyX(iBody) =-10.0_CGREAL
pointOutsideBodyy(iBody) =-10.0_CGREAL
pointOutsideBodyz(iBody) =-10.0_CGREAL

END SUBROUTINE extend_cylinder_3D

!------------------------------------------------------------------------------
!-----------------------------------------------------------------

! extend cylinder to make pseudo-3D body

SUBROUTINE extend_cylinder_vel_3D(iBody)

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays

IMPLICIT NONE

INTEGER, INTENT(IN) :: iBody

!...Loop variables
INTEGER :: i,j,k,m,mInd

DO k=2,nz
DO m=1,nPtsBodymarkerOrig(iBody)
mInd = (k-1)*nPtsBodymarkerOrig(iBody)+m
uBodyMarker(iBody,mInd) = uBodyMarker(iBody,m)
vBodyMarker(iBody,mInd) = vBodyMarker(iBody,m)
wBodyMarker(iBody,mInd) = wBodyMarker(iBody,m)
ENDDO
ENDDO

END SUBROUTINE extend_cylinder_vel_3D

!------------------------------------------------------------------------------

SUBROUTINE UnstrucSurfInDomain(iBody)

USE global_parameters
USE flow_parameters
USE boundary_arrays

IMPLICIT NONE

INTEGER, INTENT(IN) :: iBody

!...Loop variables
INTEGER :: m

!--Check whether the whole body is not in flow domain

DO m=1,nPtsBodyMarker(iBody)
IF (xBodyMarker(iBody,m)>0 .AND. xBodyMarker(iBody,m)<XOUT .AND. &
yBodyMarker(iBody,m)>0 .AND. yBodyMarker(iBody,m)<YOUT .AND. &
zBodyMarker(iBody,m)>0 .AND. ZBodyMarker(iBody,m)<ZOUT) RETURN

ENDDO

WRITE(*,'(A,I2,A)') 'The whole unstructure body(No.',iBody,') is outside flow domain.'
STOP

END SUBROUTINE UnstrucSurfInDomain

!------------------------------------------------------------------------------

SUBROUTINE WriteSurfInTecplot(iBody)

USE global_parameters
USE flow_parameters
USE boundary_arrays
USE unstructured_surface_arrays

IMPLICIT NONE

INTEGER, INTENT(IN) :: iBody

!...Loop variables
INTEGER :: m, j

!--Write surface mesh data in Tecplot Format for checking

WRITE(145,*) 'TITLE="3D TRIANGULAR SURFACE DATA"'
WRITE(145,*) 'VARIABLES="X","Y","Z"'
WRITE(145,*) 'ZONE N=',nPtsBodyMarker(iBody),',E=',totNumTriElem(iBody),'F=FEPOINT, ET=TRIANGLE'
DO m=1,nPtsBodyMarker(iBody)
WRITE(145,*) xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m)
ENDDO
DO j=1,totNumTriElem(iBody)
WRITE(145,*) triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
ENDDO

!--Write surface mesh data in Tecplot Format for checking the zone numbers !added by Chengyu
WRITE(245,*) 'TITLE="3D TRIANGULAR SURFACE DATA"'
WRITE(245,*) 'VARIABLES="X","Y","Z","Color" '
WRITE(245,*) 'ZONE N=',nPtsBodyMarker(iBody),',E=',totNumTriElem(iBody),'DATAPACKING=POINT, ZONETYPE=FETriangle'

DO m=1,nPtsBodyMarker(iBody)
WRITE(245,1000) xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m), zoneMarker(iBody,m)
ENDDO
DO j=1,totNumTriElem(iBody)
WRITE(245,*) triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
ENDDO

1000 format (3E18.6,2X,I6)
END SUBROUTINE WriteSurfInTecplot
!=================================================================================

SUBROUTINE Marker_move_rot(iBody,x1,y1,z1,quat) ! Added by G. Liu

USE global_parameters
USE flow_parameters
USE boundary_arrays
USE unstructured_surface_arrays
use body_dynamics

IMPLICIT NONE

INTEGER, INTENT(IN) :: iBody
REAL(KIND=CGREAL), INTENT(IN) :: x1,y1,z1,quat(4)

!...Loop variables
INTEGER :: m, j
REAL(KIND=CGREAL) :: vec1(3),vec2(3)


DO m=1,nPtsBodyMarker(iBody)
vec1(1)=xBodyMarkerNonIner(iBody,m)
vec1(2)=yBodyMarkerNonIner(iBody,m)
vec1(3)=zBodyMarkerNonIner(iBody,m)
call quaternion_rotation(vec1,quat)
xBodyMarker(iBody,m) = vec1(1) + x1
yBodyMarker(iBody,m) = vec1(2) + y1
zBodyMarker(iBody,m) = vec1(3) + z1
! xBodyMarkerNonIner(iBody,m) = vec1(1)
! yBodyMarkerNonIner(iBody,m) = vec1(2)
! zBodyMarkerNonIner(iBody,m) = vec1(3)
ENDDO

! DO m=1,nPtsBodyMarker(iBody)
! xBodyMarker(iBody,m) = xBodyMarkerNonIner(iBody,m) + x1
! yBodyMarker(iBody,m) = yBodyMarkerNonIner(iBody,m) + y1
! zBodyMarker(iBody,m) = zBodyMarkerNonIner(iBody,m) + z1
! ENDDO

END SUBROUTINE Marker_move_rot
