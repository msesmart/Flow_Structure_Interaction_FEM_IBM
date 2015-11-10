!---------------------------------------------------------------
! subroutine enforces general mean + sinusoidal component on translational
! and angular velocity

SUBROUTINE forced_motion(bodyNumber)

USE global_parameters
USE flow_parameters

IMPLICIT NONE

INTEGER,INTENT (IN) ::bodyNumber


vxcent(bodyNumber) = vxcentTrans(bodyNumber) &
+ ampx(bodyNumber)*SIN(2.0_CGREAL*PI*freqx(bodyNumber)*time)
vycent(bodyNumber) = vycentTrans(bodyNumber) &
+ ampy(bodyNumber)*SIN(2.0_CGREAL*PI*freqy(bodyNumber)*time)
vzcent(bodyNumber) = vzcentTrans(bodyNumber) &
+ ampz(bodyNumber)*SIN(2.0_CGREAL*PI*freqz(bodyNumber)*time)

angvx(bodyNumber) = angvxinit(bodyNumber) &
+ ampangx(bodyNumber) &
*( SIN(2.0_CGREAL*PI*freqangx(bodyNumber)*time)*cosphase(bodyNumber) &
+COS(2.0_CGREAL*PI*freqangx(bodyNumber)*time)*sinphase(bodyNumber) )
angvy(bodyNumber) = angvyinit(bodyNumber) &
+ ampangy(bodyNumber) &
*( SIN(2.0_CGREAL*PI*freqangy(bodyNumber)*time)*cosphase(bodyNumber) &
+COS(2.0_CGREAL*PI*freqangy(bodyNumber)*time)*sinphase(bodyNumber) )
angvz(bodyNumber) = angvzinit(bodyNumber) &
+ ampangz(bodyNumber) &
*( SIN(2.0_CGREAL*PI*freqangz(bodyNumber)*time)*cosphase(bodyNumber) &
+COS(2.0_CGREAL*PI*freqangz(bodyNumber)*time)*sinphase(bodyNumber) )

print *, 'angvz(bodyNumber) = ', angvz(bodyNumber)
END SUBROUTINE forced_motion

!---------------------------------------------------------------
! subroutine to specify wall velocity for slip/porous walls

SUBROUTINE wall_velocity(bodyNumber)

USE global_parameters
USE flow_parameters
USE boundary_arrays
IMPLICIT NONE

INTEGER,INTENT (IN) ::bodyNumber

INTEGER ::m,mMin,mMax,mInd,k
REAL(KIND=CGREAL) ::uBodyMarkerRel,vBodyMarkerRel,wBodyMarkerRel

mMin = 95
mMax = 105

IF (boundary_motion /= MOVING_BOUNDARY) THEN ! Rupeshs additions start here
uBodyMarker(bodyNumber,1:nPtsBodyMarker(bodyNumber)) = zero
vBodyMarker(bodyNumber,1:nPtsBodyMarker(bodyNumber)) = zero
wBodyMarker(bodyNumber,1:nPtsBodyMarker(bodyNumber)) = zero
ENDIF

DO m = mMin, mMax

uBodyMarkerRel = zero
! uBodyMarkerRel = -1.0*sin(2.0*(4.0*atan(1.0))*0.32*time) ! Uo*sin(wt)
vBodyMarkerRel = SIN(2.0_CGREAL*PI*time/2.0_CGREAL)
wBodyMarkerRel = zero
uBodyMarker(bodyNumber,m) = uBodyMarker(bodyNumber,m) + uBodyMarkerRel
vBodyMarker(bodyNumber,m) = vBodyMarker(bodyNumber,m) + vBodyMarkerRel
wBodyMarker(bodyNumber,m) = wBodyMarker(bodyNumber,m) + wBodyMarkerRel
WRITE(*,123) bodyNumber,m,uBodyMarker(bodyNumber,m)
123 FORMAT ('uBodyMarker(',I1,',',I4,') = ',E18.7)

ENDDO ! Rupeshs additions end here


! Extending velocity across span for 2D body
IF ( ndim == DIM_2D .AND. ntime > 0 ) THEN

DO k=2,nz
DO m=mMin,mMax
mInd = (k-1)*nPtsBodymarker(bodyNumber)/nz + m
uBodyMarker(bodyNumber,mInd) = uBodyMarker(bodyNumber,m)
vBodyMarker(bodyNumber,mInd) = vBodyMarker(bodyNumber,m)
wBodyMarker(bodyNumber,mInd) = wBodyMarker(bodyNumber,m)
ENDDO
ENDDO

ENDIF

END SUBROUTINE wall_velocity

!---------------------------------------------------------------

SUBROUTINE flow_induced_motion(bodyNumber)

USE global_parameters
USE flow_parameters
USE grid_arrays, ONLY : xc,yc,zc,dx,dy,dz
USE boundary_arrays, ONLY : Iblank
Use usr_module
IMPLICIT NONE

INTEGER :: i,j,k
INTEGER,INTENT (IN) ::bodyNumber

REAL :: det
REAL,DIMENSION(1:3) :: rhside

! density_solid = 10.0_CGREAL
! density_fluid = oned

density_ratio = density_solid(bodyNumber) / density_fluid

lScale = oned !=> Length scale
vScale = oned !=> Velocity scale

CALL MofI_CofG(bodyNumber) ! Moment of Inertia and Center of Gravity

! bodyforce_moment

! IF ( body_dim(bodyNumber) == 2 ) THEN
IF ( ndim == DIM_2D ) THEN

force_x = 2.0_CGREAL * sCx(bodyNumber) / ( density_fluid * lScale* (vScale)**2 *zout)
force_y = 2.0_CGREAL * sCy(bodyNumber) / ( density_fluid * lScale* (vScale)**2 *zout)
force_z = 2.0_CGREAL * sCz(bodyNumber) / ( density_fluid * lScale* (vScale)**2 *zout)
moment_x =2.0_CGREAL * scmx(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3)*zout)
moment_y =2.0_CGREAL * scmy(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3)*zout)
moment_z =2.0_CGREAL * scmz(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3)*zout)

ELSE

force_x = 2.0_CGREAL * sCx(bodyNumber) / ( density_fluid * lScale* (vScale)**2 )
force_y = 2.0_CGREAL * sCy(bodyNumber) / ( density_fluid * lScale* (vScale)**2 )
force_z = 2.0_CGREAL * sCz(bodyNumber) / ( density_fluid * lScale* (vScale)**2 )
moment_x =2.0_CGREAL * scmx(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3))
moment_y =2.0_CGREAL * scmy(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3))
moment_z =2.0_CGREAL * scmz(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3))

ENDIF

non_dim_volume = volume / ( lScale**3 )

vxcent_prev(bodyNumber) = vxcent(bodyNumber)
vycent_prev(bodyNumber) = vycent(bodyNumber)
vzcent_prev(bodyNumber) = vzcent(bodyNumber)

angvx_prev=angvx(bodyNumber)
angvy_prev=angvy(bodyNumber)
angvz_prev=angvz(bodyNumber)

! Linear Momentum Equ..(Fx i + Fy j + Fy k) = m[(du/dx)i + (dv/dy - g)j + (dw/dz)k]

! vxcent(bodyNumber)= xcentConstr(bodyNumber)*( vxcent_prev + dt*force_x*0.5/(non_dim_volume*density_ratio) )
!
! vycent(bodyNumber)= ycentConstr(bodyNumber)*( vycent_prev + dt*force_y*0.5/(non_dim_volume*density_ratio) &
! - dt*(1-1/density_ratio)/vScale**2 )
!
! vzcent(bodyNumber)= zcentConstr(bodyNumber)*( vzcent_prev + dt*force_z*0.5/(non_dim_volume*density_ratio) )

! Wanh changed vxcent_prev to vxcent_prev(bodyNumber)
vxcent(bodyNumber)= xcentConstr(bodyNumber)*( vxcent_prev(bodyNumber) + dt*force_x*0.5/(non_dim_volume*density_ratio) )

vycent(bodyNumber)= ycentConstr(bodyNumber)*( vycent_prev(bodyNumber) + dt*force_y*0.5/(non_dim_volume*density_ratio) &
- dt*(1-1/density_ratio)/vScale**2 )

vzcent(bodyNumber)= zcentConstr(bodyNumber)*( vzcent_prev(bodyNumber) + dt*force_z*0.5/(non_dim_volume*density_ratio) )

! Angular Momentum Equ..
! n+1 n n+1
!T = d([I]{omega}) /dt => ([I]{w}) - ([I]{w}) = dt * T

!Approximation used: n n+1 n n+1
!T = d([I]{omega}) /dt => [I] ({w} - {w} ) = dt * T
!
! n n+1
! rhside = ([I]{w}) + dt * T

rhside(1) = nonDimM_I(1,1)*angvx_prev + nonDimM_I(1,2)*angvy_prev + nonDimM_I(1,3)*angvz_prev + dt*moment_x
rhside(2) = nonDimM_I(2,1)*angvx_prev + nonDimM_I(2,2)*angvy_prev + nonDimM_I(2,3)*angvz_prev + dt*moment_y
rhside(3) = nonDimM_I(3,1)*angvx_prev + nonDimM_I(3,2)*angvy_prev + nonDimM_I(3,3)*angvz_prev + dt*moment_z

det = nonDimM_I(1,1)*(nonDimM_I(2,2)*nonDimM_I(3,3) - nonDimM_I(2,3)*nonDimM_I(3,2)) - &
nonDimM_I(1,2)*(nonDimM_I(2,1)*nonDimM_I(3,3) - nonDimM_I(2,3)*nonDimM_I(3,1)) + &
nonDimM_I(1,3)*(nonDimM_I(2,1)*nonDimM_I(3,2) - nonDimM_I(3,1)*nonDimM_I(2,2))

invMI(1,1) = (nonDimM_I(2,2)*nonDimM_I(3,3)-nonDimM_I(2,3)*nonDimM_I(3,2))/det
invMI(1,2) = (nonDimM_I(1,3)*nonDimM_I(3,2)-nonDimM_I(1,2)*nonDimM_I(3,3))/det
invMI(1,3) = (nonDimM_I(1,2)*nonDimM_I(2,3)-nonDimM_I(1,3)*nonDimM_I(2,2))/det

invMI(2,1) = (nonDimM_I(2,3)*nonDimM_I(3,1)-nonDimM_I(2,1)*nonDimM_I(3,3))/det
invMI(2,2) = (nonDimM_I(1,1)*nonDimM_I(3,3)-nonDimM_I(1,3)*nonDimM_I(3,1))/det
invMI(2,3) = (nonDimM_I(1,3)*nonDimM_I(2,1)-nonDimM_I(1,1)*nonDimM_I(2,3))/det

invMI(3,1) = (nonDimM_I(2,1)*nonDimM_I(3,2)-nonDimM_I(3,1)*nonDimM_I(2,2))/det
invMI(3,2) = (nonDimM_I(1,2)*nonDimM_I(3,1)-nonDimM_I(1,1)*nonDimM_I(3,2))/det
invMI(3,3) = (nonDimM_I(1,1)*nonDimM_I(2,2)-nonDimM_I(1,2)*nonDimM_I(2,1))/det


angvx(bodyNumber)= invMI(1,1)*rhside(1) + invMI(1,2)*rhside(2) + invMI(1,3)*rhside(3)
angvy(bodyNumber)= invMI(2,1)*rhside(1) + invMI(2,2)*rhside(2) + invMI(2,3)*rhside(3)
angvz(bodyNumber)= invMI(3,1)*rhside(1) + invMI(3,2)*rhside(2) + invMI(3,3)*rhside(3)

write(1001,*)'--------------------'
write(1001,*)'Volume=', volume
write(1001,*)'x,y,z', xcent(bodyNumber), ycent(bodyNumber), zcent(bodyNumber)
DO I = 1,3
write(1001,*) (nonDimM_I(I,J),J=1,3)
ENDDO
write(1002,200)time,vxcent(bodyNumber),vycent(bodyNumber),vzcent(bodyNumber),angvx(bodyNumber), &
angvy(bodyNumber),angvz(bodyNumber),force_x,force_y,force_z,moment_x,moment_y,moment_z

200 FORMAT(13(1x,E12.5))
END SUBROUTINE flow_induced_motion

!---------------------------------------------------------------
SUBROUTINE MofI_CofG(bodyNumber)

USE global_parameters
USE flow_parameters
USE grid_arrays, ONLY : xc,yc,zc,dx,dy,dz
USE boundary_arrays, ONLY : Iblank
USE usr_module
use body_dynamics, ONLY : niterFS

Implicit None

INTEGER :: i,j,k, bodyNumber

volume = zero
I_XX =zero
I_YY =zero
I_ZZ =zero
I_XY =zero
I_YZ =zero
I_XZ =zero

xcent(bodyNumber) = zero
ycent(bodyNumber) = zero
zcent(bodyNumber) = zero

! Moment of inertia at the previous time step

if((ntime .gt. 1) .and. (niterFS .eq. 0))then
do i=1,3
do j=1,3
nonDimM_I_prvs(i,j)=nonDimM_I(i,j)
enddo
enddo
endif


! Centroid and Volume

IF (ndim == DIM_2D) THEN
k = 2
DO j =2, ny-1
DO i =2, nx-1
volume = volume + dx(i)*dy(j)*Iblank(i,j,k)
xcent(bodyNumber) = xcent(bodyNumber) + xc(i)*dx(i)*dy(j)*Iblank(i,j,k)
ycent(bodyNumber) = ycent(bodyNumber) + yc(j)*dx(i)*dy(j)*Iblank(i,j,k)
ENDDO
ENDDO

ELSE
DO k =2, nz-1
DO j =2, ny-1
DO i =2, nx-1
volume = volume + dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
xcent(bodyNumber) = xcent(bodyNumber) + xc(i)*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
ycent(bodyNumber) = ycent(bodyNumber) + yc(j)*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
zcent(bodyNumber) = zcent(bodyNumber) + zc(k)*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
ENDDO
ENDDO
ENDDO

ENDIF

xcent(bodyNumber) = xcent(bodyNumber)/volume
ycent(bodyNumber) = ycent(bodyNumber)/volume
zcent(bodyNumber) = zcent(bodyNumber)/volume

! Moment of Inertia Tensor

DO k =2, nz-1
DO j =2, ny-1
DO i =2, nx-1
I_XX = I_XX + ((yc(j)-ycent(bodyNumber))**2 + (zc(k)-zcent(bodyNumber))**2)*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
I_YY = I_YY + ((xc(i)-xcent(bodyNumber))**2 + (zc(k)-zcent(bodyNumber))**2)*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
I_ZZ = I_ZZ + ((yc(j)-ycent(bodyNumber))**2 + (xc(i)-xcent(bodyNumber))**2)*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
I_XY = I_XY + ( (xc(i)-xcent(bodyNumber)) * (yc(j)-ycent(bodyNumber)) )*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
I_YZ = I_YZ + ( (zc(k)-zcent(bodyNumber)) * (yc(j)-ycent(bodyNumber)) )*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
I_XZ = I_XZ + ( (xc(i)-xcent(bodyNumber)) * (zc(k)-zcent(bodyNumber)) )*dx(i)*dy(j)*dz(k)*Iblank(i,j,k)
ENDDO
ENDDO
ENDDO

I_XX = I_XX * density_solid(bodyNumber)
I_YY = I_YY * density_solid(bodyNumber)
I_ZZ = I_ZZ * density_solid(bodyNumber)
I_XY = I_XY * density_solid(bodyNumber)
I_YZ = I_YZ * density_solid(bodyNumber)
I_XZ = I_XZ * density_solid(bodyNumber)

IF (ndim == DIM_2D) THEN ! Area Moment of Inertia
I_XX =zero
I_YY =zero
I_ZZ =zero
I_XY =zero
I_YZ =zero
I_XZ =zero
k = 2
DO j =2, ny-1
DO i =2, nx-1
I_ZZ = I_ZZ+((yc(j)-ycent(bodyNumber))**2 + (xc(i)-xcent(bodyNumber))**2)*dx(i)*dy(j)*Iblank(i,j,k)
ENDDO
ENDDO

I_ZZ = I_ZZ * density_solid(bodyNumber) !This line is added by Wanh on 04/10/11

ENDIF

! nonDimM_I(1,1) =2.0_CGREAL *I_XX/(density_fluid*lScale**5)
! nonDimM_I(1,2) =2.0_CGREAL *I_XY/(density_fluid*lScale**5)
! nonDimM_I(1,3) =2.0_CGREAL *I_XZ/(density_fluid*lScale**5)

! nonDimM_I(2,1) =2.0_CGREAL *I_XY/(density_fluid*lScale**5)
! nonDimM_I(2,2) =2.0_CGREAL *I_YY/(density_fluid*lScale**5)
! nonDimM_I(2,3) =2.0_CGREAL *I_YZ/(density_fluid*lScale**5)

! nonDimM_I(3,1) =2.0_CGREAL *I_XZ/(density_fluid*lScale**5)
! nonDimM_I(3,2) =2.0_CGREAL *I_YZ/(density_fluid*lScale**5)
! nonDimM_I(3,3) =2.0_CGREAL *I_ZZ/(density_fluid*lScale**5)

! Do not know why to multiply 2.0 in the above nonDimM_I, because
! non-dimensionalized by 0.5*density_fluid*L^5?

nonDimM_I(1,1) = I_XX/(density_fluid*lScale**5) !removed 2.0 by Wanh.
nonDimM_I(1,2) = -I_XY/(density_fluid*lScale**5)
nonDimM_I(1,3) = -I_XZ/(density_fluid*lScale**5)

nonDimM_I(2,1) = -I_XY/(density_fluid*lScale**5)
nonDimM_I(2,2) = I_YY/(density_fluid*lScale**5)
nonDimM_I(2,3) = -I_YZ/(density_fluid*lScale**5)

nonDimM_I(3,1) = -I_XZ/(density_fluid*lScale**5)
nonDimM_I(3,2) = -I_YZ/(density_fluid*lScale**5)
nonDimM_I(3,3) = I_ZZ/(density_fluid*lScale**5) !removed 2.0 by Wanh.

if(ntime.eq.1 .and. niterFS .eq. 0)then
do i=1,3
do j=1,3
nonDimM_I_prvs(i,j)=nonDimM_I(i,j)
enddo
enddo
endif

END SUBROUTINE MofI_CofG

!==========================================================
subroutine MofI_CofG_quat(iBody)
! Mass and moment of inertia are read from mass_I_noniner.dat
! Added by Geng

USE global_parameters
USE flow_parameters
USE usr_module
USE body_dynamics
USE boundary_arrays

integer :: iBody
integer :: n,i,j,k

volume = 0.0
nonDimM_I = 0.0
n=iBody
write(*,*)'liu debug in MofI_CofG_quat',iBody,n

if(ntime.gt.1 .and. niterFS .eq. 0)then
do i=1,3
do j=1,3
nonDimM_I_prvs(i,j)=nonDimM_I(i,j)
enddo
enddo
endif


IF (ndim == DIM_2D) THEN
! xcent(n)=xcent_prev(n)+dt*vxcent(n) ! quat_modify
! ycent(n)=ycent_prev(n)+dt*vycent(n)
xcent(n)=xcent_prev(n)
ycent(n)=ycent_prev(n)
zcent(n)=0.0
volume = mass_in(n)
nonDimM_I(3,3) = Icm_in(n,3)*density_solid(bodyNumber)/density_fluid
ELSE

! xcent(n)=xcent_prev(n)+dt*vxcent(n) ! quat_modify
! ycent(n)=ycent_prev(n)+dt*vycent(n)
! zcent(n)=zcent_prev(n)+dt*vzcent(n)

xcent(n)=xcent_prev(n) !+dt*vxcent(n)
ycent(n)=ycent_prev(n) !+dt*vycent(n)
zcent(n)=zcent_prev(n) !+dt*vzcent(n)

volume = mass_in(iBody)
nonDimM_I(1,1) = Icm_in(iBody,1)*density_solid(iBody)/density_fluid
nonDimM_I(1,2) =-Icm_in(iBody,4)*density_solid(iBody)/density_fluid
nonDimM_I(1,3) =-Icm_in(iBody,6)*density_solid(iBody)/density_fluid
nonDimM_I(2,1) =-Icm_in(iBody,4)*density_solid(iBody)/density_fluid
nonDimM_I(2,2) = Icm_in(iBody,2)*density_solid(iBody)/density_fluid
nonDimM_I(2,3) =-Icm_in(iBody,5)*density_solid(iBody)/density_fluid
nonDimM_I(3,1) =-Icm_in(iBody,6)*density_solid(iBody)/density_fluid
nonDimM_I(3,2) =-Icm_in(iBody,5)*density_solid(iBody)/density_fluid
nonDimM_I(3,3) = Icm_in(iBody,3)*density_solid(iBody)/density_fluid
ENDIF

if(ntime.eq.1 .and. niterFS .eq. 0)then
do i=1,3
do j=1,3
nonDimM_I_prvs(i,j)=nonDimM_I(i,j) !???? why ntime ==1
enddo
enddo
endif

! write(*,*)'liu debug', volume,mass_in(n),n, 'in subroutine MofI_CofG_quat'
! stop
endsubroutine MofI_CofG_quat
!======================================================================================
subroutine MofI_CofG_transform(iBody)
! Mass and moment of inertia are read from mass_I_noniner.dat
! The moment of inertia will be transformed to the inertia frame in this subroutine
! Added by Geng

USE global_parameters
USE flow_parameters
USE usr_module
USE body_dynamics
USE boundary_arrays

integer :: iBody
integer :: n,i,j,k
REAL(KIND=CGREAL) :: A_trans(3,3),AT_trans(3,3)

if(ntime.gt.1 .and. niterFS .eq. 0)then
do i=1,3
do j=1,3
nonDimM_I_prvs(i,j)=nonDimM_I(i,j)
enddo
enddo
endif


volume = 0.0
nonDimM_I = 0.0
n=iBody



! liu debug
if(ntime.gt.1) then
write(*,*)'********* nonDiM_I_prvs in MofI_CofG_transform ************'
write(*,*)'ntime,niterFS:',ntime,niterFS
write(*,*)'nonDimM_I_prvs:'
write(*,*)nonDimM_I_prvs(1,:)
write(*,*)nonDimM_I_prvs(2,:)
write(*,*)nonDimM_I_prvs(3,:)
write(*,*)'nonDimM_I:'
write(*,*)nonDimM_I(1,:)
write(*,*)nonDimM_I(2,:)
write(*,*)nonDimM_I(3,:)
write(*,*)'***********************************************************'
endif

IF (ndim == DIM_2D) THEN
! xcent(n)=xcent_prev(n)+dt*vxcent(n) ! quat_modify
! ycent(n)=ycent_prev(n)+dt*vycent(n)
xcent(n) = 0.5*( xBodyMarker(n,rigidRef1(n)) + xBodyMarker(n,rigidRef2(n)) ) ! matched com
ycent(n) = 0.5*( yBodyMarker(n,rigidRef1(n)) + yBodyMarker(n,rigidRef2(n)) )
zcent(n)=0.0
volume = mass_in(n)
nonDimM_I(3,3) = Icm_in(n,3)*density_solid(iBody)/density_fluid
! write(*,*)'Icm_in(n,3)',Icm_in(n,3),density_solid(iBody),density_fluid
! write(*,*)'Iz=',nonDimM_I(3,3)
! stop
ELSE

! xcent(n)=xcent_prev(n) ! liu debug
! ycent(n)=ycent_prev(n)
! zcent(n)=zcent_prev(n)

xcent(n) = 0.5*( xBodyMarker(n,rigidRef1(n)) + xBodyMarker(n,rigidRef2(n)) ) ! matched com
ycent(n) = 0.5*( yBodyMarker(n,rigidRef1(n)) + yBodyMarker(n,rigidRef2(n)) )
zcent(n) = 0.5*( zBodyMarker(n,rigidRef1(n)) + zBodyMarker(n,rigidRef2(n)) )


volume = mass_in(iBody)
nonDimM_I(1,1) = Icm_in(iBody,1)*density_solid(iBody)/density_fluid
nonDimM_I(1,2) =-Icm_in(iBody,4)*density_solid(iBody)/density_fluid
nonDimM_I(1,3) =-Icm_in(iBody,6)*density_solid(iBody)/density_fluid

nonDimM_I(2,1) =-Icm_in(iBody,4)*density_solid(iBody)/density_fluid
nonDimM_I(2,2) = Icm_in(iBody,2)*density_solid(iBody)/density_fluid
nonDimM_I(2,3) =-Icm_in(iBody,5)*density_solid(iBody)/density_fluid

nonDimM_I(3,1) =-Icm_in(iBody,6)*density_solid(iBody)/density_fluid
nonDimM_I(3,2) =-Icm_in(iBody,5)*density_solid(iBody)/density_fluid
nonDimM_I(3,3) = Icm_in(iBody,3)*density_solid(iBody)/density_fluid
ENDIF

! write(*,*) 'liu debug for I_ori in BOUNDARY_MARKER_VEL'
! write(*,*) nonDimM_I(1,:)
! write(*,*) nonDimM_I(2,:)
! write(*,*) nonDimM_I(3,:)

call get_trans_quat ! matched quaternion

if(ntime ==1 .and. niterFS .eq.0)then
open(1205,file='quaternion_debug.dat')
else
open(1205,file='quaternion_debug.dat',access='append')
endif
write(1205,1206)ntime,niterFS,quat_trans(:)
close(1205)
1206 format(I5, I4, 4(1x,e16.8))

! call get_trans_matrix(quat_prev,A_trans)
call get_trans_matrix(quat_trans,A_trans) ! matched quaternion
call get_transpose(A_trans,AT_trans)
call get_I_inertia(A_trans,nonDimM_I,AT_trans)

! write(*,*) 'liu debug for I_transform in BOUNDARY_MARKER_VEL'
! write(*,*) nonDimM_I(1,:)
! write(*,*) nonDimM_I(2,:)
! write(*,*) nonDimM_I(3,:)

! stop


if(ntime.eq.1 .and. niterFS .eq. 0)then
do i=1,3
do j=1,3
nonDimM_I_prvs(i,j)=nonDimM_I(i,j) !???? why ntime ==1
enddo
enddo
endif

! write(*,*)'liu debug', volume,mass_in(n),n, 'in subroutine MofI_CofG_quat'
! stop
endsubroutine MofI_CofG_transform

!============================================================
SUBROUTINE fea_structure_update(iBody)
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface
INTEGER :: iBody,i
DO i=1,nPtsBodyMarker(iBody)
uBodyMarker(iBody,i) = struc_vel(i,1)
vBodyMarker(iBody,i) = struc_vel(i,2)
wBodyMarker(iBody,i) = struc_vel(i,3)
xBodyMarker(iBody,i) = xBodyMarker(iBody,i)+struc_disp(i,1)
yBodyMarker(iBody,i) = yBodyMarker(iBody,i)+struc_disp(i,2)
zBodyMarker(iBody,i) = zBodyMarker(iBody,i)+struc_disp(i,3)
! Rotation have not taken into account yet
ENDDO
END SUBROUTINE fea_structure_update

!============================================================

SUBROUTINE vega_vel_update
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE fea_unstructure_surface

INTEGER :: iBody,i

iBody=1
DO i=1,nPtsBodyMarker(iBody)
uBodyMarker(iBody,i) = bodyMarkerVel(i*3-2)
vBodyMarker(iBody,i) = bodyMarkerVel(i*3-1)
wBodyMarker(iBody,i) = bodyMarkerVel(i*3)
!xBodyMarker(iBody,i) = xBodyMarker(iBody,i)+struc_disp(i,1)
!yBodyMarker(iBody,i) = yBodyMarker(iBody,i)+struc_disp(i,2)
!zBodyMarker(iBody,i) = zBodyMarker(iBody,i)+struc_disp(i,3)
! Rotation have not taken into account yet
ENDDO

END SUBROUTINE vega_vel_update
!-----------------------------------------------------------------------------------------------------
! Added by Wanh for partially or fully dynamic coupling of membrane
SUBROUTINE MofI_CofG_Mem(bodyNumber)

! For rigid body only.
! The membrane is considered as a rectangle

USE global_parameters
USE flow_parameters
USE grid_arrays, ONLY : xc,yc,zc,dx,dy,dz
USE boundary_arrays, ONLY : Iblank,xBodyMarker,yBodyMarker,zBodyMarker
USE usr_module
USE body_dynamics

Implicit None

INTEGER :: i,j,k,iSection,bodyNumber,id
REAL(KIND=CGREAL) :: xsec_min, xsec_max,ysec_min,ysec_max,zsec_min,zsec_max,section_length

LScale = 1.0

volume = zero
I_XX =zero
I_YY =zero
I_ZZ =zero
I_XY =zero
I_YZ =zero
I_XZ =zero

! Centroid and Volume

IF (nSection>=2 .and. boundary_motion_type(bodyNumber) == PARTIAL_DYNAMICS_COUPLED) THEN
xcent = zero
ycent = zero
zcent = zero

DO iSection = 2, nSection
xsec_min = 1e10
xsec_max = -1e10
ysec_min = 1e10
ysec_max = -1e10
zsec_min = 1e10
zsec_max = -1e10

DO id = 1, SectionMarker(bodyNumber,iSection)
xsec_min = min(xsec_min,xBodyMarker(bodyNumber,DynamicMarker(bodyNumber,iSection,id)))
xsec_max = max(xsec_max,xBodyMarker(bodyNumber,DynamicMarker(bodyNumber,iSection,id)))
ysec_min = min(ysec_min,yBodyMarker(bodyNumber,DynamicMarker(bodyNumber,iSection,id)))
ysec_max = max(ysec_max,yBodyMarker(bodyNumber,DynamicMarker(bodyNumber,iSection,id)))
zsec_min = min(zsec_min,zBodyMarker(bodyNumber,DynamicMarker(bodyNumber,iSection,id)))
zsec_max = max(zsec_max,zBodyMarker(bodyNumber,DynamicMarker(bodyNumber,iSection,id)))
ENDDO

section_xcent(bodyNumber,iSection) = 0.5*(xsec_min+xsec_max)
section_ycent(bodyNumber,iSection) = 0.5*(ysec_min+ysec_max)
section_zcent(bodyNumber,iSection) = 0.5*(zsec_min+zsec_max)

!print *, 'section_xcent(bodyNumber,iSection)=', section_xcent(bodyNumber,iSection)
!print *, 'section_ycent(bodyNumber,iSection)=', section_ycent(bodyNumber,iSection)
!print *, 'section_zcent(bodyNumber,iSection)=', section_zcent(bodyNumber,iSection)

section_length = sqrt((xsec_max-xsec_min)**2+(ysec_max-ysec_min)**2)
print *, 'section_length = ', section_length

! Moment of Inertia Tensor
!Here section with equal length is assumed.

I_ZZ = (density_solid(bodyNumber)*thickoverlength*section_length**3)/3.0

volume = thickoverlength*section_length*1.0

IF (ndim == DIM_3D) THEN
I_ZZ = (density_solid(bodyNumber)*DepthOverLength*thickoverlength*section_length**3)/3.0
volume = thickoverlength*section_length*DepthOverLength
ENDIF

nonDimM_I(3,3) =oned *I_ZZ/(density_fluid*LScale**5)
ENDDO !end loop of iSection
ENDIF !nSection>=2 and partially prescribed motion

IF (nSection==1) THEN
xsec_min = 1e10
xsec_max = -1e10
ysec_min = 1e10
ysec_max = -1e10
zsec_min = 1e10
zsec_max = -1e10

! xcent(bodyNumber) = (MINVAL(xBodyMarker(bodyNumber,:)+MINVAL(xBodyMarker(bodyNumber,:))/2.0_CGREAL

IF (ndim == DIM_2D) THEN
volume = thickoverlength*LScale*1.0
ELSE IF (ndim == DIM_3D) THEN
volume = thickoverlength*LScale*(DepthOverLength*LScale)
ENDIF

IF (boundary_motion_type(bodyNumber) == DYNAMICS_COUPLED.or.boundary_motion_type(bodyNumber) == BIO_DYNAMICS_COUPLED) THEN
I_ZZ = (density_solid(bodyNumber)*volume*LScale**2)/12.0
ELSEIF (boundary_motion_type(bodyNumber) == PARTIAL_DYNAMICS_COUPLED) THEN
I_ZZ = (density_solid(bodyNumber)*volume*LScale**2)/3.0
ENDIF

print *, ' I_zz, thickoverlenght, density_solid =', I_zz, thickoverlength, density_solid(bodyNumber)
nonDimM_I(3,3) =oned *I_ZZ/(density_fluid*LScale**5)

ENDIF !nSection=1

END SUBROUTINE MofI_CofG_Mem
!===========================================================================
SUBROUTINE MofI_CofG_Mem_Full_FBI(bodyNumber)
! Caculate the mass, center of gravity and moment of inertia of a membrane
! For DYNAMICS_COUPLED option
! Added by Geng

USE global_parameters
USE flow_parameters
USE grid_arrays, ONLY : xc,yc,zc,dx,dy,dz
USE boundary_arrays
USE usr_module
USE body_dynamics
USE unstructured_surface_arrays

Implicit None

INTEGER :: i,j,k,bodyNumber,id,iBody,i1,i2,i3
REAL(KIND=CGREAL) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,aa,bb,cc,pp,ss,xtc,ytc,ztc,totarea

iBody=bodyNumber
volume = zero
totarea = zero
I_XX =zero
I_YY =zero
I_ZZ =zero
I_XY =zero
I_YZ =zero
I_XZ =zero

xcent = zero
ycent = zero
zcent = zero

! Moment of inertia at the previous time step

if((ntime .gt. 1) .and. (niterFS .eq. 0))then
  do i=1,3
  do j=1,3
    nonDimM_I_prvs(i,j)=nonDimM_I(i,j)
  enddo
  enddo
endif

! Centroid and Volume

IF (ndim == DIM_2D) THEN
  write(*,*)'The present code cannot deal with a full-FBI problem for a 2D membrane. Stop!'
  stop
ELSE
  DO j=1,totNumTriElem(iBody)
    i1 = TriElemNeig(iBody,1,j)
    i2 = TriElemNeig(iBody,2,j)
    i3 = TriElemNeig(iBody,3,j)
    x1 = xBodyMarker(iBody,i1)
    y1 = yBodyMarker(iBody,i1)
    z1 = zBodyMarker(iBody,i1)
    x2 = xBodyMarker(iBody,i2)
    y2 = yBodyMarker(iBody,i2)
    z2 = zBodyMarker(iBody,i2)
    x3 = xBodyMarker(iBody,i3)
    y3 = yBodyMarker(iBody,i3)
    z3 = zBodyMarker(iBody,i3)
    xtc = (x1+x2+x3)/3.0
    ytc = (y1+y2+y3)/3.0
    ztc = (z1+z2+z3)/3.0

    aa = dsqrt((x1-x2)**2.0+(y1-y2)**2.0+(z1-z2)**2.0)
    bb = dsqrt((x1-x3)**2.0+(y1-y3)**2.0+(z1-z3)**2.0)
    cc = dsqrt((x2-x3)**2.0+(y2-y3)**2.0+(z2-z3)**2.0)
    pp = 0.5*(aa+bb+cc)
    ss = dsqrt(pp*(pp-aa)*(pp-bb)*(pp-cc))
    totarea = totarea + ss

    xcent(iBody) = xcent(iBody) + xtc*ss
    ycent(iBody) = ycent(iBody) + ytc*ss
    zcent(iBody) = zcent(iBody) + ztc*ss
  ENDDO ! j loop

  volume = totarea*thickoverlength
  xcent(bodyNumber) = xcent(bodyNumber)/totarea
  ycent(bodyNumber) = ycent(bodyNumber)/totarea
  zcent(bodyNumber) = zcent(bodyNumber)/totarea

  !--- Moment of inertia -------------
  DO j=1,totNumTriElem(iBody)
    i1 = TriElemNeig(iBody,1,j)
    i2 = TriElemNeig(iBody,2,j)
    i3 = TriElemNeig(iBody,3,j)
    x1 = xBodyMarker(iBody,i1)
    y1 = yBodyMarker(iBody,i1)
    z1 = zBodyMarker(iBody,i1)
    x2 = xBodyMarker(iBody,i2)
    y2 = yBodyMarker(iBody,i2)
    z2 = zBodyMarker(iBody,i2)
    x3 = xBodyMarker(iBody,i3)
    y3 = yBodyMarker(iBody,i3)
    z3 = zBodyMarker(iBody,i3)
    xtc = (x1+x2+x3)/3.0
    ytc = (y1+y2+y3)/3.0
    ztc = (z1+z2+z3)/3.0
    aa = dsqrt((x1-x2)**2.0+(y1-y2)**2.0+(z1-z2)**2.0)
    bb = dsqrt((x1-x3)**2.0+(y1-y3)**2.0+(z1-z3)**2.0)
    cc = dsqrt((x2-x3)**2.0+(y2-y3)**2.0+(z2-z3)**2.0)
    pp = 0.5*(aa+bb+cc)
    ss = dsqrt(pp*(pp-aa)*(pp-bb)*(pp-cc))*thickoverlength

    I_xx = I_xx + ( (ytc- ycent(iBody))**2.0 + (ztc-zcent(iBody))**2.0 )*ss
    I_yy = I_yy + ( (xtc- xcent(iBody))**2.0 + (ztc-zcent(iBody))**2.0 )*ss
    I_zz = I_zz + ( (xtc- xcent(iBody))**2.0 + (ytc-ycent(iBody))**2.0 )*ss
    I_xy = I_xy + ( (ytc- ycent(iBody))*(xtc-xcent(iBody)) )*ss
    I_xz = I_xz + ( (ztc- zcent(iBody))*(xtc-xcent(iBody)) )*ss
    I_yz = I_yz + ( (ztc- zcent(iBody))*(ytc-ycent(iBody)) )*ss
  ENDDO
  I_XX = I_XX * density_solid(bodyNumber)
  I_YY = I_YY * density_solid(bodyNumber)
  I_ZZ = I_ZZ * density_solid(bodyNumber)
  I_XY = I_XY * density_solid(bodyNumber)
  I_YZ = I_YZ * density_solid(bodyNumber)
  I_XZ = I_XZ * density_solid(bodyNumber)

  ! non-dimensionalized by 0.5*density_fluid*L^5?

!  nonDimM_I(1,1) = I_XX/(density_fluid*lScale**5) !removed 2.0 by Wanh.
!  nonDimM_I(1,2) = -I_XY/(density_fluid*lScale**5)
!  nonDimM_I(1,3) = -I_XZ/(density_fluid*lScale**5)
!
!  nonDimM_I(2,1) = -I_XY/(density_fluid*lScale**5)
!  nonDimM_I(2,2) = I_YY/(density_fluid*lScale**5)
!  nonDimM_I(2,3) = -I_YZ/(density_fluid*lScale**5)
!
!  nonDimM_I(3,1) = -I_XZ/(density_fluid*lScale**5)
!  nonDimM_I(3,2) = -I_YZ/(density_fluid*lScale**5)
!  nonDimM_I(3,3) = I_ZZ/(density_fluid*lScale**5) !removed 2.0 by Wanh.


  nonDimM_I(1,1) =  I_XX
  nonDimM_I(1,2) = -I_XY
  nonDimM_I(1,3) = -I_XZ

  nonDimM_I(2,1) = -I_XY
  nonDimM_I(2,2) =  I_YY
  nonDimM_I(2,3) = -I_YZ

  nonDimM_I(3,1) = -I_XZ
  nonDimM_I(3,2) = -I_YZ
  nonDimM_I(3,3) =  I_ZZ


ENDIF

if(ntime.eq.1 .and. niterFS .eq. 0)then
  do i=1,3
  do j=1,3
    nonDimM_I_prvs(i,j)=nonDimM_I(i,j)
  enddo
  enddo
endif

! write the file (liu debug)
if(ntime.eq.1 .and. niterFS.eq.0)then
  open(4001,file='m_I_debug.dat')
  write(4001,*)'ntime,niterFS,volume,Ixx,Iyy,Izz,-Ixy,-Ixz,-Iyz'
else
  open(4001,file='m_I_debug.dat',access='append')
endif
write(4001,4002)ntime,niterFS,volume,nonDimM_I(1,1),nonDimM_I(2,2),nonDimM_I(3,3),nonDimM_I(1,2),nonDimM_I(1,3),nonDimM_I(2,3)
4002 format(I5,I3,7(1x,e12.5))
close(4001)

END SUBROUTINE MofI_CofG_Mem_Full_FBI
!=============================================================================
SUBROUTINE dynamics_motion(iBody)
! Calculate the acceleration and velocity of the center of gravity.
USE global_parameters
USE flow_parameters
USE usr_module
USE body_dynamics
USE boundary_arrays

! USE GCM_arrays, ONLY : iPvt, work

IMPLICIT NONE

INTEGER :: iBody, iSection,i
INTEGER :: info
INTEGER, PARAMETER :: iRowMax=8
integer :: niter_fb

INTEGER :: iPvt(iRowMax)

REAL(KIND=CGREAL) :: non_dim_mass,forceXtmp,forceYtmp,forceZtmp,momentXtmp,momentYtmp,momentZtmp
REAL(KIND=CGREAL) :: LScaleCube
REAL(KIND=CGREAL) :: h,M
REAL(KIND=CGREAL) :: Istar, Fr_square, Fr_Belmonte, U0_Belmonte
REAL(KIND=CGREAL) :: tscale_down, tscale_pendulum
REAL(KIND=CGREAL) :: ycentcur,vycentcur,acc_yCGcur
! REAL(KIND=CGREAL) :: vxcent_wo_relax,vycent_wo_relax,vzcent_wo_relax
REAL(KIND=CGREAL) :: deltaVxcent_Aitken_prev,deltaVycent_Aitken_prev,deltaVzcent_Aitken_prev
REAL(KIND=CGREAL) :: dAitkenU,dAitkenV,dAitkenW
! REAL(KIND=CGREAL) :: angvx_wo_relax,angvy_wo_relax,angvz_wo_relax
REAL(KIND=CGREAL) :: deltaAngvx_Aitken_prev,deltaAngvy_Aitken_prev,deltaAngvz_Aitken_prev
REAL(KIND=CGREAL) :: dAitkenAngx,dAitkenAngy,dAitkenAngz
REAL(KIND=CGREAL) :: Aitken_nominator,Aitken_denominator
REAL(KIND=CGREAL) :: angRHS(BODY_DIM3),angLHS(BODY_DIM3,BODY_DIM3),Lmatrx(BODY_DIM3,BODY_DIM3)
REAL(KIND=CGREAL) :: Iomg1,Iomg2,Iomg3
REAL(KIND=CGREAL) :: work(iRowMax)
REAL(KIND=CGREAL), PARAMETER :: NUM_ZERO = 1E-30
REAL(KIND=CGREAL) :: freq
REAL(KIND=CGREAL) :: K_theta, theta, theta_limit, d_theta
REAL(KIND=CGREAL) :: tip_x, tip_y, hinge_x, hinge_y
REAL(KIND=CGREAL) :: uref, vref
REAL(KIND=CGREAL) :: accx_o, accy_o
REAL(KIND=CGREAL) :: theta_1,omg_1,angzacc_1,length_1, length_section
REAL(KIND=CGREAL) :: noninertial_moment
REAL(KIND=CGREAL) :: dx_center2hinge, dy_center2hinge
REAL(KIND=CGREAL) :: omega_in(3) ! Added by Geng

character fname*17

density_ratio = density_solid(iBody)/density_fluid

LScale = 1.0
LScaleCube = LScale*LScale*LScale

h = thickoverlength*LScale

d_theta = 2*PI/180.0 ! 2-degree


!if(ntime.gt.1) then
!write(*,*)'********* nonDiM_I_prvs in dynamics_motion ************'
!write(*,*)'ntime,niterFS:',ntime,niterFS
!write(*,*)nonDimM_I_prvs(1,:)
!write(*,*)nonDimM_I_prvs(2,:)
!write(*,*)nonDimM_I_prvs(3,:)
!endif

IF (boundary_motion_type(iBody) == DYNAMICS_COUPLED .or. boundary_motion_type(iBody) == BIO_DYNAMICS_COUPLED .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT .or. boundary_motion_type(iBody) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR ) THEN

VScale = sqrt( 2.0_CGREAL*(density_ratio-1)*h*GRAVITY_ACC )
ELSE IF (boundary_motion_type(iBody) == PARTIAL_DYNAMICS_COUPLED) THEN
! for character freq 1
freq = Freq_prsb
! VScale = 2*(atan(1.0)*4.0)*freq*A_x
VScale = 1
ENDIF


WRITE(*,*) '------------------------------------------------------'
WRITE(*,*) 'Calling dynamics_motion ...'
WRITE(*,*) ' iBody = ', iBody

! if(ntime==1 .and. niterFS==0)then
! open(1028,file='volume.dat')
! else
! open(1028,file='volume.dat',access='append')
! endif
! write(1028,*)'volume=',volume,'ntime=',ntime,'niterFS=',niterFS
! close(1028)

IF (boundary_motion_type(iBody) == DYNAMICS_COUPLED .or. boundary_motion_type(iBody) == BIO_DYNAMICS_COUPLED .or. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT .OR. boundary_motion_type(iBody) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_SWIMMING ) THEN

 IF(boundary_motion_type(iBody) /= DYNAMICS_COUPLED_QUAT .and. boundary_motion_type(iBody) /= DYNAMICS_COUPLED_MofI_QUAT) THEN
  IF (nbody_solid >0) THEN
   CALL MofI_CofG(iBody) ! Moment of Inertia, volume, and Center of Gravity
  ELSE IF (nbody_membrane >0) THEN
   CALL MofI_CofG_Mem_full_FBI(iBody) ! Moment of Inertia, volume, and Center of Gravity
  ENDIF

  print *, ' Dimensional I (area) =', nonDimM_I(3,3)*(LScale**5) &
                                       *density_fluid/density_solid(iBody)
  print *, ' nonDimensional I =', nonDimM_I(3,3)
  endif

non_dim_volume = volume/LScaleCube
non_dim_mass = density_ratio*non_dim_volume !Non-dimensionlzed Body Mass
write(*,*)'liu debug', volume,non_dim_mass !

IF (ntime < 2) THEN

  U0_Belmonte = sqrt((density_ratio-1)*GRAVITY_ACC*non_dim_volume) !Belmonte 1998
  tscale_down = LScale/U0_Belmonte
  tscale_pendulum = sqrt(density_ratio/(density_ratio-1)*LScale/GRAVITY_ACC)
  Fr_Belmonte = sqrt(non_dim_mass)

  WRITE(*,*) ' Density ratio of Body:', iBody, ' is:', density_ratio
  WRITE(*,*) ' VScale = ', VScale, ' U0_Belmonte = ', U0_Belmonte
  WRITE(*,*) ' Fr_Belmonte =', Fr_Belmonte

  !Transfer I from MofI_CofG to the Istar, defined in Andersen et al (2005 I)
  WRITE(*,*) ' Dimensional I (area) =', nonDimM_I(3,3)*(LScale**5) &
                                        *density_fluid/density_solid(iBody)
ENDIF

print *, ' scy(iBody) =', scy(iBody)
print *, ' scmz(iBody) =', scmz(iBody)

IF (ndim == DIM_2D) THEN

  if(boundary_motion_type(iBody) == BIO_DYNAMICS_COUPLED)then
    forceXtmp=zero
    forceYtmp=zero
    momentZtmp=zero
    do i=1,nBody
      forceXtmp=forceXtmp+scx(i)
      forceYtmp=forceYtmp+scy(i)
      momentZtmp=momentZtmp+scmz(i)
    end do
    force_x = forceXtmp/zout/non_dim_volume/density_ratio
    force_y = forceYtmp/zout/non_dim_volume/density_ratio - 0.5/thickoverlength/density_ratio
    moment_z = momentZtmp/zout/LScaleCube


! else if(boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT)THEN ! Added by Geng
! force_x = scx(iBody)/zout/mass_in(iBody)
! force_y = scy(iBody)/zout/mass_in(iBody) - 0.5/thickoverlength/density_ratio
! moment_z = scmz(iBody)/zout/LScaleCube

  else
    force_x = scx(iBody)/zout/non_dim_volume/density_ratio
    if(boundary_motion_type(iBody) .ge. DYNAMICS_COUPLED_SWIMMING)then
       force_y = scy(iBody)/zout/non_dim_volume/density_ratio
    else
       force_y = scy(iBody)/zout/non_dim_volume/density_ratio - 0.5/thickoverlength/density_ratio
    endif
    moment_z = scmz(iBody)/zout/LScaleCube
    !aero_moment = scmz(iBody)/zout/LScaleCube/VScale/VScale
  end if


!    ! liu debug
!    if(ntime==1 .and. niterFS == 0)then
!      open(999,file='niter0_debug.dat')
!    else
!      open(999,file='niter0_debug.dat',status='old',access='append')
!    endif
!      write(999,*)'ntime=',ntime,niterFS
!      write(999,*)scx(iBody),scy(iBody),scmz(iBody)
!    close(999)

ELSE
 if(boundary_motion_type(iBody) == BIO_DYNAMICS_COUPLED)then

 forceXtmp=zero
 forceYtmp=zero
 forceZtmp=zero
 momentXtmp=zero
 momentYtmp=zero
 momentZtmp=zero
 do i=1,nBody
   forceXtmp=forceXtmp+scx(i)
   forceYtmp=forceYtmp+scy(i)
   forceZtmp=forceZtmp+scz(i)
   momentXtmp=momentXtmp+scmx(i)
   momentYtmp=momentYtmp+scmy(i)
   momentZtmp=momentZtmp+scmz(i)
 end do

 force_x = forceXtmp/non_dim_mass
 force_y = forceYtmp/non_dim_mass - 0.5/thickoverlength/density_ratio !?????
 force_z = forceZtmp/non_dim_mass
 moment_x = momentXtmp/LScaleCube
 moment_y = momentYtmp/LScaleCube

 !aero_moment = scmz(iBody)/LScaleCube/VScale/VScale
 aero_moment = momentZtmp/LScaleCube
 moment_z = aero_moment

! else if(boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT)then ! Added by Geng
! force_x = scx(iBody)/mass_in(iBody)
! force_y = scy(iBody)/mass_in(iBody) - 0.5/thickoverlength/density_ratio
! force_z = scz(iBody)/mass_in(iBody)
! moment_x = scmx(iBody)/LScaleCube
! moment_y = scmy(iBody)/LScaleCube
!
! aero_moment = scmz(iBody)/LScaleCube
! moment_z = aero_moment

 else

  force_x = scx(iBody)/non_dim_mass
  if(boundary_motion_type(iBody) .ge. DYNAMICS_COUPLED_SWIMMING)then
    force_y = scy(iBody)/non_dim_mass
  else
    force_y = scy(iBody)/non_dim_mass - 0.5/thickoverlength/density_ratio
  endif
  force_z = scz(iBody)/non_dim_mass
  moment_x = scmx(iBody)/LScaleCube
  moment_y = scmy(iBody)/LScaleCube

  !aero_moment = scmz(iBody)/LScaleCube/VScale/VScale
  aero_moment = scmz(iBody)/LScaleCube
  moment_z = aero_moment
 end if
ENDIF

acc_xCG(iBody) = force_x
acc_yCG(iBody) = force_y
acc_zCG(iBody) = force_z

vxcent_iter(iBody) = vxcent(iBody)
vycent_iter(iBody) = vycent(iBody)
vzcent_iter(iBody) = vzcent(iBody)

angvx_iter(iBody) = angvx(iBody)
angvy_iter(iBody) = angvy(iBody)
angvz_iter(iBody) = angvz(iBody)

SELECT CASE (ndim)

CASE (DIM_2D)
  if(DoF_on(ibody,1) == 1)then  ! Added by G. Liu, for DoF control, DoF: Degree of Freedom
     vxcent_wo_relax = vxcent_prev(iBody) + acc_xCG(iBody)*dt
  else
     vxcent_wo_relax = zero
  endif

  if(DoF_on(ibody,2) == 1)then
     vycent_wo_relax = vycent_prev(iBody) + acc_yCG(iBody)*dt
  else
     vycent_wo_relax = zero
  endif
  vzcent_wo_relax = zero

  angvx_wo_relax = zero
  angvy_wo_relax = zero
  if(DoF_on(ibody,6) == 1)then
    angvz_wo_relax = angvz_old(iBody) + moment_z/nonDimM_I(3,3)*dt
  else
    angvz_wo_relax = zero
  endif
  ! print *, 'moment_z, nonDimM_I(3,3)=', moment_z, nonDimM_I(3,3)
  WRITE (STDOUT,'(A,1X,1PE20.10,1X,1PE20.10)') ' W/o relax angvz,vycent:', angvz_wo_relax,vycent_wo_relax

CASE (DIM_3D)
  if(DoF_on(ibody,1) == 1)then
    vxcent_wo_relax = vxcent_prev(iBody) + acc_xCG(iBody)*dt
  else
    vxcent_wo_relax = zero
  endif

  if(DoF_on(ibody,2) == 1)then
    vycent_wo_relax = vycent_prev(iBody) + acc_yCG(iBody)*dt
  else
    vycent_wo_relax = zero
  endif

  if(DoF_on(ibody,3) == 1)then
    vzcent_wo_relax = vzcent_prev(iBody) + acc_zCG(iBody)*dt
  else
    vzcent_wo_relax = zero
  endif

  Iomg1 = nonDimM_I(1,1)*angvx_old(iBody) + nonDimM_I(1,2)*angvy_old(iBody) &
        + nonDimM_I(1,3)*angvz_old(iBody)
  Iomg2 = nonDimM_I(2,1)*angvx_old(iBody) + nonDimM_I(2,2)*angvy_old(iBody) &
        + nonDimM_I(2,3)*angvz_old(iBody)
  Iomg3 = nonDimM_I(3,1)*angvx_old(iBody) + nonDimM_I(3,2)*angvy_old(iBody) &
        + nonDimM_I(3,3)*angvz_old(iBody)

IF(boundary_motion_type(iBody) == DYNAMICS_COUPLED .or. boundary_motion_type(iBody) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR .or. boundary_motion_type(iBody) == DYNAMICS_COUPLED_SWIMMING )THEN ! Added by Geng
  AngvIntg = 3
ELSEIF(boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT .OR. &
  boundary_motion_type(iBody) == BIO_DYNAMICS_COUPLED )THEN
  AngvIntg = 2
ENDIF ! Added by Geng

IF (AngvIntg == 0) THEN
  Lmatrx(1,2) = Iomg3*angvy_old(iBody)*dt
  Lmatrx(1,3) = Iomg2*angvz_old(iBody)*dt
  Lmatrx(2,3) = Iomg1*angvz_old(iBody)*dt
  Lmatrx(2,1) = Iomg3*angvx_old(iBody)*dt
  Lmatrx(3,1) = Iomg2*angvx_old(iBody)*dt
  Lmatrx(3,2) = Iomg1*angvy_old(iBody)*dt

  angRHS(1) = Iomg1 + moment_x*dt - Lmatrx(1,2) + Lmatrx(1,3)
  angRHS(2) = Iomg2 + moment_y*dt - Lmatrx(2,3) + Lmatrx(2,1)
  angRHS(3) = Iomg3 + moment_z*dt - Lmatrx(3,1) + Lmatrx(3,2)

  angLHS(1,1) = nonDimM_I(1,1)
  angLHS(1,2) = nonDimM_I(1,2)
  angLHS(1,3) = nonDimM_I(1,3)
  angLHS(2,1) = nonDimM_I(2,1)
  angLHS(2,2) = nonDimM_I(2,2)
  angLHS(2,3) = nonDimM_I(2,3)
  angLHS(3,1) = nonDimM_I(3,1)
  angLHS(3,2) = nonDimM_I(3,2)
  angLHS(3,3) = nonDimM_I(3,3)

  angvx_wo_relax = zero
  angvy_wo_relax = zero
  if(DoF_on(iBody,6)==1)then ! added by G. Liu
    angvz_wo_relax = angvz_old(iBody) + moment_z/nonDimM_I(3,3)*dt
  else
    angvz_wo_relax = zero
  endif

ELSE IF (AngvIntg == 1) THEN
  Iomg1 = nonDimM_I(1,1)*angvx_old(iBody) + nonDimM_I(1,2)*angvy_old(iBody) &
  + nonDimM_I(1,3)*angvz_old(iBody)
  Iomg2 = nonDimM_I(2,1)*angvx_old(iBody) + nonDimM_I(2,2)*angvy_old(iBody) &
  + nonDimM_I(2,3)*angvz_old(iBody)
  Iomg3 = nonDimM_I(3,1)*angvx_old(iBody) + nonDimM_I(3,2)*angvy_old(iBody) &
  + nonDimM_I(3,3)*angvz_old(iBody)

  Lmatrx(1,2) = 0.5*Iomg3*dt
  Lmatrx(1,3) = 0.5*Iomg2*dt
  Lmatrx(2,1) = Lmatrx(1,2)
  Lmatrx(2,3) = 0.5*Iomg1*dt
  Lmatrx(3,1) = Lmatrx(1,3)
  Lmatrx(3,2) = Lmatrx(2,3)

  angRHS(1) = Iomg1 + moment_x*dt - angvy_old(iBody)*Lmatrx(1,2) + angvz_old(iBody)*Lmatrx(1,3)
  angRHS(2) = Iomg2 + moment_y*dt - angvz_old(iBody)*Lmatrx(2,3) + angvx_old(iBody)*Lmatrx(2,1)
  angRHS(3) = Iomg3 + moment_z*dt - angvx_old(iBody)*Lmatrx(3,1) + angvy_old(iBody)*Lmatrx(3,2)

  angLHS(1,1) = nonDimM_I(1,1)
  angLHS(1,2) = nonDimM_I(1,2) + Lmatrx(1,2)
  angLHS(1,3) = nonDimM_I(1,3) - Lmatrx(1,3)
  angLHS(2,1) = nonDimM_I(2,1) - Lmatrx(2,1)
  angLHS(2,2) = nonDimM_I(2,2)
  angLHS(2,3) = nonDimM_I(2,3) + Lmatrx(2,3)
  angLHS(3,1) = nonDimM_I(3,1) + Lmatrx(3,1)
  angLHS(3,2) = nonDimM_I(3,2) - Lmatrx(3,2)
  angLHS(3,3) = nonDimM_I(3,3)

ELSE IF (AngvIntg == 2) THEN
  Lmatrx(1,1) = ( nonDimM_I(3,1)*angvy_old(iBody) - nonDimM_I(2,1)*angvz_old(iBody) )*dt
  Lmatrx(2,2) = ( nonDimM_I(1,2)*angvz_old(iBody) - nonDimM_I(3,2)*angvx_old(iBody) )*dt
  Lmatrx(3,3) = ( nonDimM_I(2,3)*angvx_old(iBody) - nonDimM_I(1,3)*angvy_old(iBody) )*dt

  Lmatrx(1,2) = ( nonDimM_I(3,2)*angvy_old(iBody) + nonDimM_I(3,3)*angvz_old(iBody) )*dt
  Lmatrx(1,3) = ( nonDimM_I(2,2)*angvy_old(iBody) + nonDimM_I(2,3)*angvz_old(iBody) )*dt
  Lmatrx(2,3) = ( nonDimM_I(1,1)*angvx_old(iBody) + nonDimM_I(1,3)*angvz_old(iBody) )*dt
  Lmatrx(2,1) = ( nonDimM_I(3,1)*angvx_old(iBody) + nonDimM_I(3,3)*angvz_old(iBody) )*dt
  Lmatrx(3,1) = ( nonDimM_I(2,1)*angvx_old(iBody) + nonDimM_I(2,2)*angvy_old(iBody) )*dt
  Lmatrx(3,2) = ( nonDimM_I(1,1)*angvx_old(iBody) + nonDimM_I(1,2)*angvy_old(iBody) )*dt

  angRHS(1) = Iomg1 + moment_x*dt - angvy_old(iBody)*Lmatrx(1,2) + angvz_old(iBody)*Lmatrx(1,3)
  angRHS(2) = Iomg2 + moment_y*dt - angvz_old(iBody)*Lmatrx(2,3) + angvx_old(iBody)*Lmatrx(2,1)
  angRHS(3) = Iomg3 + moment_z*dt - angvx_old(iBody)*Lmatrx(3,1) + angvy_old(iBody)*Lmatrx(3,2)

  angLHS(1,1) = nonDimM_I(1,1) + Lmatrx(1,1)
  angLHS(1,2) = nonDimM_I(1,2)
  angLHS(1,3) = nonDimM_I(1,3)
  angLHS(2,1) = nonDimM_I(2,1)
  angLHS(2,2) = nonDimM_I(2,2) + Lmatrx(2,2)
  angLHS(2,3) = nonDimM_I(2,3)
  angLHS(3,1) = nonDimM_I(3,1)
  angLHS(3,2) = nonDimM_I(3,2)
  angLHS(3,3) = nonDimM_I(3,3) + Lmatrx(3,3)

ELSE IF (AngvIntg == 3) THEN
  Iomg1 = nonDimM_I_prvs(1,1)*angvx_old(ibody)+nonDimM_I_prvs(1,2)*angvy_old(ibody)+nonDimM_I_prvs(1,3)*angvz_old(ibody)
  Iomg2 = nonDimM_I_prvs(2,1)*angvx_old(ibody)+nonDimM_I_prvs(2,2)*angvy_old(ibody)+nonDimM_I_prvs(2,3)*angvz_old(ibody)
  Iomg3 = nonDimM_I_prvs(3,1)*angvx_old(ibody)+nonDimM_I_prvs(3,2)*angvy_old(ibody)+nonDimM_I_prvs(3,3)*angvz_old(ibody)

  angRHS(1) = Iomg1 + moment_x*dt
  angRHS(2) = Iomg2 + moment_y*dt
  angRHS(3) = Iomg3 + moment_z*dt

  angLHS(1,1) = nonDimM_I(1,1)
  angLHS(1,2) = nonDimM_I(1,2)
  angLHS(1,3) = nonDimM_I(1,3)
  angLHS(2,1) = nonDimM_I(2,1)
  angLHS(2,2) = nonDimM_I(2,2)
  angLHS(2,3) = nonDimM_I(2,3)
  angLHS(3,1) = nonDimM_I(3,1)
  angLHS(3,2) = nonDimM_I(3,2)
  angLHS(3,3) = nonDimM_I(3,3)

ENDIF ! end of if (AngvIntg == 0/1/2/3)

CALL DGETRF(3, 3, angLHS,3,iPvt, info)
CALL DGETRI(3, angLHS,3,iPvt,work, 3, info)

if(DoF_on(iBody,4) == 1)then
  angvx_wo_relax = angLHS(1,1)*angRHS(1) + angLHS(1,2)*angRHS(2) + angLHS(1,3)*angRHS(3)
else
  angvx_wo_relax = zero
endif

if(DoF_on(iBody,5) == 1)then
  angvy_wo_relax = angLHS(2,1)*angRHS(1) + angLHS(2,2)*angRHS(2) + angLHS(2,3)*angRHS(3)
else
  angvy_wo_relax = zero
endif

if(DoF_on(iBody,6) == 1)then
  angvz_wo_relax = angLHS(3,1)*angRHS(1) + angLHS(3,2)*angRHS(2) + angLHS(3,3)*angRHS(3)
else
  angvz_wo_relax = zero
endif

END SELECT !body_dim

CALL Aitken_dyn(iBody)

WRITE (STDOUT,'(4X,A,1X,1PE20.10,1X,1PE20.10)') 'vxcent, vycent in dynamics:', vxcent(iBody),vycent(iBody)
WRITE (STDOUT,'(4X,A,1X,1PE20.10)') 'angvz in dynamics:', angvz(iBody)

! liu debug
if(ntime .eq. 1 .and. niterFS .eq. 0)then
open(1212,file='dynamics_debug.dat')
else
open(1212,file='dynamics_debug.dat',access='append')
endif
write(1212,*) '**************** ntime',ntime,' niterFS',niterFS,' *******************'
write(1212,*) 'I'
write(1212,*) nonDimM_I(1,:)
write(1212,*) nonDimM_I(2,:)
write(1212,*) nonDimM_I(3,:)
write(1212,*) 'I_prev:'
write(1212,*) nonDimM_I_prvs(1,:)
write(1212,*) nonDimM_I_prvs(2,:)
write(1212,*) nonDimM_I_prvs(3,:)
write(1212,*) 'moment in x, y, z'
write(1212,*) moment_x,moment_y,moment_z
write(1212,*) 'angv_old:', angvx_old(ibody),angvy_old(ibody),angvz_old(ibody)
write(1212,*) 'angv_rlx:', angvx_wo_relax,angvy_old(ibody),angvz_old(ibody)
write(1212,*) 'angv :', angvx(iBody),angvy(iBody),angvz(iBody)
write(1212,*) ' '
close(1212)

! IF DYNAMICS_COUPLED_QUAT, angv|non-inertial -> angv|inertial

if(boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT)then
  call angv_noniner2iner(iBody)
  omega_in(1)=angvx_noniner(iBody)
  omega_in(2)=angvy_noniner(iBody)
  omega_in(3)=angvz_noniner(iBody)
  call quat_update(quat_prev,omega_in,1,quat_iter) ! Quaternion is updated based on the angular velocity in inertia frame, so frameflag=1
endif

if(boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT)then
  omega_in(1) = angvx(iBody)
  omega_in(2) = angvy(iBody)
  omega_in(3) = angvz(iBody)
  call quat_update(quat_prev,omega_in,2,quat_iter) ! Quaternion is updated based on the angular velocity in inertia frame, so frameflag=2
endif

write(fname,101)'kine_body_',iBody,'.dat'
101 format(A10,I3.3,A4)

IF (NTIME == 1 .and. niterFS==0)then
open(ifuBodyCGPath+iBody,file=fname)
WRITE (ifuBodyCGPath+iBody,'(A)') 'ntime niterFS acc_xCG acc_yCG acc_zCG vxcent &
vycent vzcent xcent ycent zcent angvx angvy angvz'
ELSE
open(ifuBodyCGPath+iBody,file=fname,access='append')
ENDIF
WRITE (ifuBodyCGPath+iBody, 121) ntime,niterFS,acc_xCG(iBody),acc_yCG(iBody),acc_zCG(iBody), &
vxcent(iBody),vycent(iBody),vzcent(iBody), &
xcent(iBody),ycent(iBody),zcent(iBody), &
angvx(iBody),angvy(iBody),angvz(iBody)
close(ifuBodyCGPath+iBody)
ENDIF ! end of DYNAMICS (nBody_solid)

!================================
! Partially dynamics coupling   !
!================================

IF (boundary_motion_type(iBody) == PARTIAL_DYNAMICS_COUPLED) THEN
IF (nbody_solid >0) THEN
CALL MofI_CofG(iBody) ! Moment of Inertia, volume, and Center of Gravity
ELSE IF (nbody_membrane >0) THEN
CALL MofI_CofG_Mem(iBody) ! Moment of Inertia, volume, and Center of Gravity
ENDIF

print *, ' Dimensional I (area) =', nonDimM_I(3,3)*(LScale**5) &
*density_fluid/density_solid(iBody)
print *, ' nonDimensional I =', nonDimM_I(3,3)

!Note volume is only of the portion which experiences FSI. Portion under
!prescribed portion is not included.
non_dim_volume = volume/LScaleCube
non_dim_mass = density_ratio*non_dim_volume !Non-dimensionlzed Body Mass

! K_theta = 0 !Spring constant
K_theta = 4*PI*PI*FreqN_torsion*FreqN_torsion*nonDimM_I(3,3)
hinge_x = xBodyMarker(iBody,hingemarker)
hinge_y = yBodyMarker(iBody,hingemarker)
print *, 'A_x =', A_x
accx_o = -2*PI*PI*freq*freq*A_x*cos(2*PI*freq*ntime*dt) !LE where prescribed motion added.
accy_o = 0
theta_1 = A_theta*sin(2*PI*freq*ntime*dt) !Link on which motion is prescribed
omg_1 = A_theta*2*PI*freq*cos(2*PI*freq*ntime*dt)
angzacc_1 = -4*PI*PI*A_theta*freq*freq*sin(2*PI*freq*ntime*dt)
length_1 = sqrt( (xBodyMarker(iBody,hingemarker)-xBodyMarker(iBody,1))**2 + &
(yBodyMarker(iBody,hingemarker)-yBodyMarker(iBody,1))**2)

accx_hinge = accx_o + angzacc_1*length_1*cos(theta_1) - omg_1*omg_1*length_1*sin(theta_1)
accy_hinge = accy_o + angzacc_1*length_1*sin(theta_1) + omg_1*omg_1*length_1*cos(theta_1)
!uref = uBodyMarker(iBody,hingemarker)
uref = uBodyMarker(iBody,1) !Note this is different from HingedPlate4_Constraint2.

tip_x = xBodyMarker(iBody,INT(nPtsBodyMarker(iBody)/nz))
tip_y = yBodyMarker(iBody,INT(nPtsBodyMarker(iBody)/nz))
print *, ' hingemarker = ', hingemarker
print *, ' hinge_x,hinge_y =', hinge_x, hinge_y
print *, ' tip_x,tip_y =', tip_x, tip_y

IF (NTIME <= 2) THEN
Limiter_change_direction = .False.
ENDIF

IF ( abs(hinge_x-tip_x)<1e-10 ) THEN
theta = 0
ELSE
theta = -atan( (tip_x-hinge_x)/(tip_y-hinge_y) )
ENDIF

print *, ' theta =', theta
print *, ' yBodyMarker(iBody,hingemarker)=', yBodyMarker(iBody,hingemarker)
print *, ' section_ycent(iBody,2)=', section_ycent(iBody,2)
print *, ' ntime, Limiter_On=', NTIME, Limiter_On

! The forces are non-dimensioned by rho_f * U**2 * L**2
! The moments are non-dimensioned by rho_f * U**2 * L**3

DO iSection = 2,nSection
dx_center2hinge = section_xcent(iBody,iSection)-xBodyMarker(iBody,hingemarker)
dy_center2hinge = section_ycent(iBody,iSection)-yBodyMarker(iBody,hingemarker)
dy_center2hinge = -dy_center2hinge ! add negative since positive dy_center2hinge wanted.
print *, 'dx_center2hinge, dy_center2hinge = ', dx_center2hinge, dy_center2hinge
print *, 'accx_hinge, accy_hinge =', accx_hinge, accy_hinge

IF (ndim == DIM_2D) THEN
!moment w.r.t. the hinge point
!The follwing "-" before scmz is just for temporary
!moment_z = scmz(iSection)/zout/LScaleCube &
!- 1.0/Fr_square*(section_xcent(iBody,iSection)-xBodyMarker(iBody,hingemarker))*non_dim_volume

!Use character time T=1/f=1, character speed U=L/T=1
!Note scmz is divided by VScale^2
aero_moment = scmz(iSection)/zout/LScaleCube/VScale/VScale
grav_moment = - (density_ratio-1)*GRAVITY_ACC*dx_center2hinge*non_dim_volume
noninertial_moment = -density_ratio*non_dim_volume*dy_center2hinge*accx_hinge &
-density_ratio*non_dim_volume*dx_center2hinge*accy_hinge

moment_z = aero_moment + grav_moment + noninertial_moment

ELSE
aero_moment = scmz(iSection)/LScaleCube/VScale/VScale
grav_moment = - (density_ratio-1)*GRAVITY_ACC &
* (section_xcent(iBody,iSection)-xBodyMarker(iBody,hingemarker))*non_dim_volume
noninertial_moment = -density_ratio*non_dim_volume*dy_center2hinge*accx_hinge &
-density_ratio*non_dim_volume*dx_center2hinge*accy_hinge
moment_z = aero_moment + grav_moment + noninertial_moment

ENDIF !end of ndim

print *, ' dabs(angvz(iSection)) from prev. step:', dabs(angvz(iSection))

IF ( MOD(ntime,ntimePerCycle(iBody)) == 1 &
.OR. MOD(ntime,ntimePerCycle(iBody)) == ntimePerCycle(iBody)/2+1 &
! Limiter_On changes from T to F
.OR. (moment_z_pretime*moment_z<=0 .AND. dabs(angvz(iSection))<1e-10 .AND. &
Limiter_change_direction .eqv. .False. ) ) THEN
Limiter_On = .False.
! print *, ' in here 1?'
ENDIF


IF (Limiter_change_direction) THEN
IF ( ( MOD(ntime,ntimePerCycle(iBody))>ntimePerCycle(iBody)/4 .and. &
MOD(ntime,ntimePerCycle(iBody))<ntimePerCycle(iBody)/2 ) &
.or. MOD(ntime,ntimePerCycle(iBody))>ntimePerCycle(iBody)/4*3 ) THEN
Limiter_change_direction = .False.
ENDIF
ENDIF

print *, ' aero moment, threshold =', aero_moment, aero_moment_threshold
print *, ' grav moment =', grav_moment
print *, ' noninertial moment =', noninertial_moment
print *, ' Limiter_On is ', Limiter_On
print *, ' niterFS = ', niterFS
print *, ' moment_z, moment_z_pretime = ', moment_z, moment_z_pretime
print *, ' theta and limit:', theta, Limiter_angle_amp

angvz_iter(iSection) = angvz(iSection)

SELECT CASE (ndim)

CASE (DIM_2D)
angvx_wo_relax = zero
angvy_wo_relax = zero
angvz_wo_relax = angvz_old(iSection) + (moment_z - K_theta*theta)/nonDimM_I(3,3)*dt/2/PI/freq/freq

CASE (DIM_3D)
angvx_wo_relax = zero
angvy_wo_relax = zero
angvz_wo_relax = angvz_old(iSection) + (moment_z - K_theta*theta)/nonDimM_I(3,3)*dt/2/PI/freq/freq

END SELECT !body_dim

!CALL Aitken_hingemembrane(iBody,iSection,theta,d_theta) !add by yan
CALL Aitken_hingemembrane(iBody,iSection)

print *, ' angvz(iSection) in dynamics : ', angvz(iSection)

write(fname,101)'kine_body_',iBody,'.dat'

IF (NTIME == 1 .and. niterFS==0)then
open(ifuBodyCGPath+iBody,file=fname)
WRITE (ifuBodyCGPath+iBody,'(A)') 'ntime niterFS acc_xCG acc_yCG acc_zCG vxcent &
vycent vzcent xcent ycent zcent angvx angvy angvz'
ELSE
open(ifuBodyCGPath+iBody,file=fname,access='append')
ENDIF
WRITE (ifuBodyCGPath+iBody, 121) ntime,niterFS,acc_xCG(iBody),acc_yCG(iBody),acc_zCG(iBody), &
vxcent(iBody),vycent(iBody),vzcent(iBody), &
xcent(iBody),ycent(iBody),zcent(iBody), &
angvx(iBody),angvy(iBody),angvz(iBody)
close(ifuBodyCGPath+iBody)


! IF (NTIME == 1)then
!
! open(ifuBodyCGPath+iBody,FILE=
! WRITE (ifuBodyCGPath+iBody,'(A)') 'time acc_xCG acc_yCG acc_zCG vxcent &
! vycent vzcent xcent ycent zcent angvx angvy angvz'
! ENDIF
! WRITE (ifuBodyCGPath+iBody, 121) dt*ntime,acc_xCG(iSection),acc_yCG(iSection), acc_zCG(iSection), &
! vxcent(iSection),vycent(iSection),vzcent(iSection), &
! xcent(iSection),ycent(iSection),zcent(iSection), &
! angvx(iSection),angvy(iSection),angvz(iSection)
ENDDO ! end loop of iSection.

ENDIF !end of nBody_membrane

121 FORMAT (I5,1x,I4,1x,12(E13.6,1x))

END SUBROUTINE dynamics_motion
!======================================================================================================
subroutine dynamics_motion_quat(iBody)
USE global_parameters
USE flow_parameters
USE usr_module
USE body_dynamics
USE boundary_arrays

integer :: iBody




ENDSUBROUTINE dynamics_motion_quat
! -------------------------------------------------------------------------------------------------------
SUBROUTINE Aitken_dyn(iBody)
USE body_dynamics
USE flow_parameters, ONLY : vxcent,vycent,vzcent, &
angvx,angvy,angvz, &
angvx_old,angvy_old,angvz_old, &
ntime,ntimePerCycle
USE usr_module

IMPLICIT NONE

INTEGER :: iBody

REAL(KIND=CGREAL) :: deltaVxcent_Aitken_prev,deltaVycent_Aitken_prev,deltaVzcent_Aitken_prev
REAL(KIND=CGREAL) :: dAitkenU,dAitkenV,dAitkenW
REAL(KIND=CGREAL) :: deltaAngvx_Aitken_prev,deltaAngvy_Aitken_prev,deltaAngvz_Aitken_prev
REAL(KIND=CGREAL) :: dAitkenAngx,dAitkenAngy,dAitkenAngz
REAL(KIND=CGREAL) :: Aitken_nominator,Aitken_denominator
REAL(KIND=CGREAL), PARAMETER :: NUM_ZERO = 1E-30

IF (AITKEN) THEN

IF (niterFS <1e-1) THEN !First iteration
deltaVxcent_Aitken(iBody) = vxcent_prev(iBody) - vxcent_wo_relax
deltaVycent_Aitken(iBody) = vycent_prev(iBody) - vycent_wo_relax
deltaVzcent_Aitken(iBody) = vzcent_prev(iBody) - vzcent_wo_relax

deltaAngvx_Aitken(iBody) = angvx_old(iBody) - angvx_wo_relax
deltaAngvy_Aitken(iBody) = angvy_old(iBody) - angvy_wo_relax
deltaAngvz_Aitken(iBody) = angvz_old(iBody) - angvz_wo_relax

LambdaVel_Aitken(iBody) = LambdaVel_Aitken_Init
LambdaAngVel_Aitken(iBody) = LambdaAngVel_Aitken_Init

vxcent(iBody) = LambdaVel_Aitken(iBody)*vxcent_prev(iBody) + (1-LambdaVel_Aitken(iBody))*vxcent_wo_relax
vycent(iBody) = LambdaVel_Aitken(iBody)*vycent_prev(iBody) + (1-LambdaVel_Aitken(iBody))*vycent_wo_relax
vzcent(iBody) = LambdaVel_Aitken(iBody)*vzcent_prev(iBody) + (1-LambdaVel_Aitken(iBody))*vzcent_wo_relax
angvx(iBody) = LambdaAngVel_Aitken(iBody)*angvx_old(iBody) + (1-LambdaAngVel_Aitken(iBody))*angvx_wo_relax
angvy(iBody) = LambdaAngVel_Aitken(iBody)*angvy_old(iBody) + (1-LambdaAngVel_Aitken(iBody))*angvy_wo_relax
angvz(iBody) = LambdaAngVel_Aitken(iBody)*angvz_old(iBody) + (1-LambdaAngVel_Aitken(iBody))*angvz_wo_relax

ELSE
deltaVxcent_Aitken_prev = deltaVxcent_Aitken(iBody)
deltaVycent_Aitken_prev = deltaVycent_Aitken(iBody)
deltaVzcent_Aitken_prev = deltaVzcent_Aitken(iBody)
deltaAngvx_Aitken_prev = deltaAngvx_Aitken(iBody)
deltaAngvy_Aitken_prev = deltaAngvy_Aitken(iBody)
deltaAngvz_Aitken_prev = deltaAngvz_Aitken(iBody)

! write(*,*) 'Aitken_dyn_TEST1',deltaVxcent_Aitken_prev
! write(*,*) 'Aitken_dyn_TEST1',deltaVycent_Aitken_prev
! write(*,*) 'Aitken_dyn_TEST1',deltaVzcent_Aitken_prev
! write(*,*) 'Aitken_dyn_TEST1',deltaAngvx_Aitken_prev
! write(*,*) 'Aitken_dyn_TEST1',deltaAngvy_Aitken_prev
! write(*,*) 'Aitken_dyn_TEST1',deltaAngvz_Aitken_prev

deltaVxcent_Aitken(iBody) = vxcent_iter(iBody) - vxcent_wo_relax
deltaVycent_Aitken(iBody) = vycent_iter(iBody) - vycent_wo_relax
deltaVzcent_Aitken(iBody) = vzcent_iter(iBody) - vzcent_wo_relax

! write(*,*) 'Aitken_dyn_TEST2',deltaVxcent_Aitken(iBody)
! write(*,*) 'Aitken_dyn_TEST2',deltaVycent_Aitken(iBody)
! write(*,*) 'Aitken_dyn_TEST2',deltaVzcent_Aitken(iBody)


dAitkenU = deltaVxcent_Aitken_prev - deltaVxcent_Aitken(iBody)
dAitkenV = deltaVycent_Aitken_prev - deltaVycent_Aitken(iBody)
dAitkenW = deltaVzcent_Aitken_prev - deltaVzcent_Aitken(iBody)

! write(*,*) 'Aitken_dyn_TEST3',dAitkenU
! write(*,*) 'Aitken_dyn_TEST3',dAitkenV
! write(*,*) 'Aitken_dyn_TEST3',dAitkenW

Aitken_nominator = dAitkenU*deltaVxcent_Aitken(iBody) + dAitkenV*deltaVycent_Aitken(iBody) &
+ dAitkenW*deltaVzcent_Aitken(iBody)

Aitken_denominator = dAitkenU*dAitkenU + dAitkenV*dAitkenV + dAitkenW*dAitkenW

! write(*,*) 'Aitken_dyn_TEST4',Aitken_nominator
! write(*,*) 'Aitken_dyn_TEST4',Aitken_denominator

! print *, 'ibody:',ibody, 'Aitken nom/denom:', Aitken_nominator,Aitken_nominator
! print *, 'vx/ycent_wo_relax:', vxcent_wo_relax,vycent_wo_relax
! print *, 'deltaVx/ycent_Aitken_prev:', deltaVxcent_Aitken_prev,deltaVycent_Aitken_prev

IF (abs(Aitken_denominator)>NUM_ZERO) THEN
LambdaVel_Aitken(iBody) = LambdaVel_Aitken(iBody) + (LambdaVel_Aitken(iBody)-1)* &
Aitken_nominator/Aitken_denominator
ENDIF
IF (abs(LambdaVel_Aitken(iBody)-1.0)<1e-5) LambdaVel_Aitken(iBody) = LambdaVel_Aitken_Init

! WRITE(*,*) 'UnderRelaxation factor:', 1-LambdaVel_Aitken(iBody)

vxcent(iBody) = LambdaVel_Aitken(iBody)*vxcent_iter(iBody) + (1-LambdaVel_Aitken(iBody))*vxcent_wo_relax
vycent(iBody) = LambdaVel_Aitken(iBody)*vycent_iter(iBody) + (1-LambdaVel_Aitken(iBody))*vycent_wo_relax
vzcent(iBody) = LambdaVel_Aitken(iBody)*vzcent_iter(iBody) + (1-LambdaVel_Aitken(iBody))*vzcent_wo_relax
! print *, 'vx/y/zcent:', vxcent(iBody), vycent(iBody),vzcent(iBody)

! write(*,*) 'Aitken_dyn_TEST5',vxcent(iBody)
! write(*,*) 'Aitken_dyn_TEST5',vycent(iBody)
! write(*,*) 'Aitken_dyn_TEST5',vzcent(iBody)

deltaAngvx_Aitken(iBody) = angvx_iter(iBody) - angvx_wo_relax
deltaAngvy_Aitken(iBody) = angvy_iter(iBody) - angvy_wo_relax
deltaAngvz_Aitken(iBody) = angvz_iter(iBody) - angvz_wo_relax

! write(*,*) 'Aitken_dyn_TEST6',deltaAngvx_Aitken(iBody)
! write(*,*) 'Aitken_dyn_TEST6',deltaAngvy_Aitken(iBody)
! write(*,*) 'Aitken_dyn_TEST6',deltaAngvz_Aitken(iBody)

dAitkenAngx = deltaAngvx_Aitken_prev - deltaAngvx_Aitken(iBody)
dAitkenAngy = deltaAngvy_Aitken_prev - deltaAngvy_Aitken(iBody)
dAitkenAngz = deltaAngvz_Aitken_prev - deltaAngvz_Aitken(iBody)

! write(*,*) 'Aitken_dyn_TEST7',dAitkenAngx
! write(*,*) 'Aitken_dyn_TEST7',dAitkenAngy
! write(*,*) 'Aitken_dyn_TEST7',dAitkenAngz

Aitken_nominator = dAitkenAngx*deltaAngvx_Aitken(iBody) + dAitkenAngy*deltaAngvy_Aitken(iBody) &
+ dAitkenAngz*deltaAngvz_Aitken(iBody)

Aitken_denominator = dAitkenAngx*dAitkenAngx + dAitkenAngy*dAitkenAngy + dAitkenAngz*dAitkenAngz

! write(*,*) 'Aitken_dyn_TEST8',Aitken_nominator
! write(*,*) 'Aitken_dyn_TEST8',Aitken_denominator

IF (abs(Aitken_denominator)>NUM_ZERO) THEN
LambdaAngVel_Aitken(iBody) = LambdaAngVel_Aitken(iBody) + (LambdaAngVel_Aitken(iBody)-1)* &
Aitken_nominator/Aitken_denominator
ENDIF

IF (abs(LambdaAngVel_Aitken(iBody)-1.0)<1e-5) LambdaAngVel_Aitken(iBody) = LambdaAngVel_Aitken_Init

angvx(iBody) = LambdaAngVel_Aitken(iBody)*angvx_iter(iBody) + (1-LambdaAngVel_Aitken(iBody))*angvx_wo_relax
angvy(iBody) = LambdaAngVel_Aitken(iBody)*angvy_iter(iBody) + (1-LambdaAngVel_Aitken(iBody))*angvy_wo_relax
angvz(iBody) = LambdaAngVel_Aitken(iBody)*angvz_iter(iBody) + (1-LambdaAngVel_Aitken(iBody))*angvz_wo_relax

! write(*,*) 'Aitken_dyn_TEST9',angvx(iBody)
! write(*,*) 'Aitken_dyn_TEST9',angvy(iBody)
! write(*,*) 'Aitken_dyn_TEST9',angvz(iBody)


ENDIF !endif of niterFS

ELSE IF (.NOT. AITKEN) THEN

vxcent(iBody) = vxcent_iter(iBody) + LambdaVel_Aitken_Init* &
(vxcent_wo_relax - vxcent_iter(iBody))
vycent(iBody) = vycent_iter(iBody) + LambdaVel_Aitken_Init* &
(vycent_wo_relax - vycent_iter(iBody))
vzcent(iBody) = vzcent_iter(iBody) + LambdaVel_Aitken_Init* &
(vzcent_wo_relax - vzcent_iter(iBody))
angvx(iBody) = angvx_iter(iBody) + LambdaAngVel_Aitken_Init* &
(angvx_wo_relax - angvx_iter(iBody) )
angvy(iBody) = angvy_iter(iBody) + LambdaAngVel_Aitken_Init* &
(angvy_wo_relax - angvy_iter(iBody) )
angvz(iBody) = angvz_iter(iBody) + LambdaAngVel_Aitken_Init* &
(angvz_wo_relax - angvz_iter(iBody) )

! The above is actually for one body, since moment_z is not an array.
ENDIF

END SUBROUTINE Aitken_dyn
! -------------------------------------------------------------------------------------------------------
SUBROUTINE Aitken_hingemembrane(iBody,iSection)

USE body_dynamics
USE flow_parameters, ONLY : vxcent,vycent,vzcent, &
angvx,angvy,angvz, &
angvx_old,angvy_old,angvz_old, &
ntime,ntimePerCycle
USE usr_module

IMPLICIT NONE

INTEGER :: iBody, iSection

! REAL(KIND=CGREAL) :: vxcent_wo_relax,vycent_wo_relax,vzcent_wo_relax
REAL(KIND=CGREAL) :: deltaVxcent_Aitken_prev,deltaVycent_Aitken_prev,deltaVzcent_Aitken_prev
REAL(KIND=CGREAL) :: dAitkenU,dAitkenV,dAitkenW
! REAL(KIND=CGREAL) :: angvx_wo_relax,angvy_wo_relax,angvz_wo_relax
REAL(KIND=CGREAL) :: deltaAngvx_Aitken_prev,deltaAngvy_Aitken_prev,deltaAngvz_Aitken_prev
REAL(KIND=CGREAL) :: dAitkenAngx,dAitkenAngy,dAitkenAngz
REAL(KIND=CGREAL) :: Aitken_nominator,Aitken_denominator
REAL(KIND=CGREAL) :: theta, theta_limit, d_theta
REAL(KIND=CGREAL), PARAMETER :: NUM_ZERO = 1E-30

IF (AITKEN) THEN
write(*,*) "HAHAHAHAHAAHAHAH,niterFS=",niterFS

IF (niterFS <1e-1) THEN !First iteration
deltaVxcent_Aitken(iSection) = zero
deltaVycent_Aitken(iSection) = zero
deltaAngvx_Aitken(iSection) = zero
deltaAngvy_Aitken(iSection) = zero
deltaAngvz_Aitken(iSection) = angvz_old(iSection) - angvz_wo_relax

LambdaVel_Aitken(iSection) = LambdaVel_Aitken_Init
LambdaAngVel_Aitken(iSection) = LambdaAngVel_Aitken_Init

angvz(iSection) = LambdaAngVel_Aitken(iSection)*angvz_old(iSection) &
+ (1-LambdaAngVel_Aitken(iSection))*angvz_wo_relax

print *, ' iSection=',iSection, ' AngVel_aitken=', LambdaAngVel_Aitken(iSection)
print *, ' angvz_old, wo_relax, angvz:',angvz_old(iSection),angvz_wo_relax, angvz(iSection)

IF (.NOT. Limiter_On) THEN

IF (dabs(theta-Limiter_angle_amp) < d_theta &
.AND. angvz(iSection)*angvz_old(iSection)>NUM_ZERO ) THEN
! The plate is accelerating, and approaching to the limiter.
Limiter_On = .TRUE.
print *, ' Limiter_On is ', Limiter_On
print *, ' theta and limit:', theta, Limiter_angle_amp
print *, ' here 2'
theta_limit = -1.0 * theta_limit
Limiter_change_direction = .True.
Limiter_angle_amp = -1.0*Limiter_angle_amp
ENDIF

ENDIF

IF ( MOD(ntime,ntimePerCycle(iBody)) <= 20 &
.OR.( MOD(ntime,ntimePerCycle(iBody)) > ntimePerCycle(iBody)/2 .AND. &
MOD(ntime,ntimePerCycle(iBody)) <= ntimePerCycle(iBody)/2+20 ) ) THEN
Limiter_On = .False.
print *, 'here 4'
ENDIF

IF (Limiter_On) angvz(iSection) = 0

ELSE
deltaVxcent_Aitken_prev = deltaVxcent_Aitken(iSection)
deltaVycent_Aitken_prev = deltaVycent_Aitken(iSection)
deltaAngvz_Aitken_prev = deltaAngvz_Aitken(iSection)
deltaAngvx_Aitken_prev = deltaAngvx_Aitken(iSection)
deltaAngvy_Aitken_prev = deltaAngvy_Aitken(iSection)

write(*,*) 'TEST,1:::::',deltaVxcent_Aitken_prev
write(*,*) 'TEST,1:::::',deltaVycent_Aitken_prev
write(*,*) 'TEST,1:::::',deltaAngvz_Aitken_prev

dAitkenU = deltaVxcent_Aitken_prev - deltaVxcent_Aitken(iSection)
dAitkenV = deltaVycent_Aitken_prev - deltaVycent_Aitken(iSection)
! dAitkenW = deltaVzcent_Aitken_prev - deltaVzcent_Aitken(iSection)
write(*,*) 'TEST,2:::::',dAitkenU
write(*,*) 'TEST,2:::::',dAitkenV

deltaVxcent_Aitken(iSection) = zero
deltaVycent_Aitken(iSection) = zero
deltaAngvz_Aitken(iSection) = angvz_iter(iSection) - angvz_wo_relax
deltaAngvx_Aitken(iSection) = zero
deltaAngvy_Aitken(iSection) = zero
write(*,*) 'TEST,3:::::',deltaAngvz_Aitken(iSection)

dAitkenAngx = deltaAngvx_Aitken_prev - deltaAngvx_Aitken(iSection)
dAitkenAngy = deltaAngvy_Aitken_prev - deltaAngvy_Aitken(iSection)
dAitkenAngz = deltaAngvz_Aitken_prev - deltaAngvz_Aitken(iSection)

write(*,*) 'TEST,4:::::',dAitkenAngx
write(*,*) 'TEST,4:::::',dAitkenAngy
write(*,*) 'TEST,4:::::',dAitkenAngz

Aitken_nominator = dAitkenAngx*deltaAngvx_Aitken(iSection) + dAitkenAngy*deltaAngvy_Aitken(iSection) &
+ dAitkenAngz*deltaAngvz_Aitken(iSection)

write(*,*) 'TEST,5:::::',Aitken_nominator

Aitken_denominator = dAitkenAngx*dAitkenAngx + dAitkenAngy*dAitkenAngy + dAitkenAngz*dAitkenAngz

write(*,*) 'TEST,6:::::',Aitken_denominator

IF (abs(Aitken_denominator)>NUM_ZERO) THEN
LambdaAngVel_Aitken(iSection) = LambdaAngVel_Aitken(iSection) + (LambdaAngVel_Aitken(iSection)-1)* &
Aitken_nominator/Aitken_denominator

write(*,*) 'TEST,7:::::',LambdaAngVel_Aitken(iSection)
ENDIF

IF (abs(LambdaAngVel_Aitken(iSection)-1.0)<1e-5) LambdaAngVel_Aitken(iSection) = LambdaAngVel_Aitken_Init

angvz(iSection) = LambdaAngVel_Aitken(iSection)*angvz_iter(iSection) + (1-LambdaAngVel_Aitken(iSection))*angvz_wo_relax
write(*,*) 'TEST,8:::::',angvz(iSection)
IF (Limiter_On) angvz(iSection) = 0

ENDIF !endif of niterFS

ELSE IF (.NOT. AITKEN) THEN

IF (abs(aero_moment)<aero_moment_threshold) THEN
angvz(iSection) = angvz_iter(iSection) + LambdaAngVel_Aitken_Init* &
(angvz_wo_relax - angvz_iter(iSection) )
ELSE
angvz(iSection) = 0
ENDIF

! The above is actually for one body, since moment_z is not an array.
ENDIF

END SUBROUTINE Aitken_hingemembrane
! -------------------------------------------------------------------------------------------------------
SUBROUTINE vega_markerVel_convergenceCheck  ! Added by CJ Yuan
  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE usr_module
  USE body_dynamics
  IMPLICIT NONE
  INTEGER :: i,iBody
  REAL(KIND=CGREAL) :: max_diff

  max_diff=0.0
  iBody=1
  DO i=1,nPtsBodyMarker(1)
    IF(abs(bodyMarkerVel(i*3-2)-uBodyMarker(iBody,i))>max_diff)THEN
      max_diff=abs(bodyMarkerVel(i*3-2)-uBodyMarker(iBody,i))
    END IF
    IF(abs(bodyMarkerVel(i*3-1)-vBodyMarker(iBody,i))>max_diff)THEN
      max_diff=abs(bodyMarkerVel(i*3-1)-vBodyMarker(iBody,i))
    END IF
    IF(abs(bodyMarkerVel(i*3)-wBodyMarker(iBody,i))>max_diff)THEN
      max_diff=abs(bodyMarkerVel(i*3)-wBodyMarker(iBody,i))
    END IF
  ENDDO

  IF(max_diff<=1.0e-3)THEN
    Converged_FSI(iBody) = .True.
  ELSE
    WRITE (*,*) 'Vega FSI Not Converged',max_diff
    DO i=1,nPtsBodyMarker(1)
      uBodyMarker(iBody,i)=0.5*uBodyMarker(iBody,i)+0.5*bodyMarkerVel(i*3-2)
      vBodyMarker(iBody,i)=0.5*vBodyMarker(iBody,i)+0.5*bodyMarkerVel(i*3-1)
      wBodyMarker(iBody,i)=0.5*wBodyMarker(iBody,i)+0.5*bodyMarkerVel(i*3)
      bodyMarkerVel(i*3-2)=uBodyMarker(iBody,i)
      bodyMarkerVel(i*3-1)=vBodyMarker(iBody,i)
      bodyMarkerVel(i*3)=wBodyMarker(iBody,i)
    ENDDO
  END IF
  !Converged_FSI(iBody) = .True. ! Add for weak FSI set(No convergence iteration)

END SUBROUTINE vega_markerVel_convergenceCheck

! -------------------------------------------------------------------------------------------------------
SUBROUTINE FSConvergeCheck(iBody)

USE flow_parameters, ONLY : ntime,vxcent,vycent,vzcent,angvx,angvy,angvz,ndim,boundary_motion_type
USE global_parameters, ONLY : DIM_2D, DIM_3D
USE body_dynamics
USE usr_module

IMPLICIT NONE

INTEGER :: iBody,iSection
LOGICAL :: converge_log

REAL(KIND=CGREAL) :: vxcentdiff,vycentdiff,vzcentdiff
REAL(KIND=CGREAL) :: angvxdiff, angvydiff, angvzdiff
REAL(KIND=CGREAL) :: converge_cretia, converge_cretia_ang
REAL(KIND=CGREAL), PARAMETER :: ZEROVEL = 1e-6

WRITE (*,*) 'Calling FSConvergeCheck ...'
IF (ntime>5) THEN
  converge_cretia = Convg_Cretia/5.0
  converge_cretia_ang = Convg_Cretia_Ang/5.0
ELSE
  converge_cretia = Convg_Cretia
  converge_cretia_ang = Convg_Cretia_Ang
ENDIF

converge_log = .False.

SELECT CASE (ndim)

  CASE (DIM_2D)
    IF (boundary_motion_type(iBody) == DYNAMICS_COUPLED .or. boundary_motion_type(iBody) == BIO_DYNAMICS_COUPLED.or. &
    boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT .OR. boundary_motion_type(iBody) == DYNAMICS_COUPLED_MofI_QUAT .OR.&
    boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
    boundary_motion_type(iBody) == DYNAMICS_COUPLED_SWIMMING) THEN
      IF ( abs(vxcent(iBody)) < ZEROVEL .and. abs(vxcent_iter(iBody)) < ZEROVEL ) THEN
         vxcentdiff = 0.0
      ELSE
         vxcentdiff = abs( (vxcent(iBody)-vxcent_iter(iBody))/vxcent(iBody) )
      ENDIF

      IF ( abs(vycent(iBody)) < ZEROVEL .and. abs(vycent_iter(iBody)) < ZEROVEL ) THEN
        vycentdiff = 0.0
      ELSE
        vycentdiff = abs( (vycent(iBody)-vycent_iter(iBody))/vycent(iBody) )
      ENDIF

      IF ( abs(angvz(iBody)) < ZEROVEL .and. abs(angvz_iter(iBody)) < ZEROVEL) THEN
        angvzdiff = 0.0
      ELSE
        angvzdiff = abs( (angvz(iBody)-angvz_iter(iBody))/angvz(iBody) )
      ENDIF

    ENDIF !end of boundary_motion_type == DYNAMICS_COUPLED

    IF (boundary_motion_type(iBody) == PARTIAL_DYNAMICS_COUPLED) THEN
      DO iSection = 2, nSection
        IF (NTIME>1 .and. abs(angvz(iSection))<1e-20) THEN
          Converged_FSI(iBody) = .True.
          WRITE (*,*) ' Body:', iBody, 'at Ntime =', ntime, ' niterFS =',niterFS
          RETURN
        ENDIF

        IF ( abs(angvz(iSection)) < ZEROVEL .and. abs(angvz_iter(iSection)) < ZEROVEL ) THEN
          print *, ' angvz(iSection)=', angvz(iSection)
          angvzdiff = 0.0
        ELSE
          angvzdiff = abs( (angvz(iSection)-angvz_iter(iSection))/angvz(iSection) )
        ENDIF
! converge_log might be an array for section number >2.
      ENDDO ! end loop of iSection
    ENDIF !end of boundary_motion_type == PARTIAL_DYNAMICS_COUPLED

  converge_log = ( vxcentdiff<converge_cretia .and. &
  vycentdiff<converge_cretia .and. &
  angvzdiff<converge_cretia_ang )
  print *, ' converge_log=', converge_log

CASE (DIM_3D)
IF (boundary_motion_type(iBody) == DYNAMICS_COUPLED.or. boundary_motion_type(iBody) == BIO_DYNAMICS_COUPLED .or. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_QUAT .OR. boundary_motion_type(iBody) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
boundary_motion_type(iBody) == DYNAMICS_COUPLED_SWIMMING) THEN
IF ( abs(vxcent(iBody)) < ZEROVEL .and. abs(vxcent_iter(iBody)) < ZEROVEL ) THEN
vxcentdiff = 0.0
ELSE
vxcentdiff = abs( (vxcent(iBody)-vxcent_iter(iBody))/vxcent(iBody) )
ENDIF

IF ( abs(vycent(iBody)) < ZEROVEL .and. abs(vycent_iter(iBody)) < ZEROVEL ) THEN
vycentdiff = 0.0
ELSE
vycentdiff = abs( (vycent(iBody)-vycent_iter(iBody))/vycent(iBody) )
ENDIF

IF ( abs(vzcent(iBody)) < ZEROVEL .and. abs(vzcent_iter(iBody)) < ZEROVEL ) THEN
vzcentdiff = 0.0
ELSE
vzcentdiff = abs( (vzcent(iBody)-vzcent_iter(iBody))/vzcent(iBody) )
ENDIF

IF ( abs(angvx(iBody)) < ZEROVEL .and. abs(angvx_iter(iBody)) < ZEROVEL ) THEN
angvxdiff = 0.0
ELSE
angvxdiff = abs( (angvx(iBody)-angvx_iter(iBody))/angvx(iBody) )
ENDIF

IF ( abs(angvy(iBody)) < ZEROVEL .and. abs(angvy_iter(iBody)) < ZEROVEL ) THEN
angvydiff = 0.0
ELSE
angvydiff = abs( (angvy(iBody)-angvy_iter(iBody))/angvy(iBody) )
ENDIF

ENDIF !end of boundary_motion_type == DYNAMICS_COUPLED

IF (boundary_motion_type(iBody) == PARTIAL_DYNAMICS_COUPLED) THEN
DO iSection = 2, nSection

IF (NTIME>1 .and. abs(angvz(iSection))<1e-20) THEN
Converged_FSI(iBody) = .True.
WRITE (*,*) ' Body:', iBody, 'at Ntime =', ntime, ' niterFS =',niterFS
RETURN
ENDIF

IF ( abs(angvz(iSection)) < ZEROVEL .and. abs(angvz_iter(iSection)) < ZEROVEL ) THEN
angvzdiff = 0.0
ELSE
angvzdiff = abs( (angvz(iSection)-angvz_iter(iSection))/angvz(iSection) )
ENDIF

ENDDO ! end loop of iSection
ENDIF !end of boundary_motion_type == PARTIAL_DYNAMICS_COUPLED

converge_log = ( vxcentdiff<converge_cretia .and. &
vycentdiff<converge_cretia .and. &
vzcentdiff<converge_cretia .and. &
angvxdiff<converge_cretia_ang .and. &
angvydiff<converge_cretia_ang .and. &
angvzdiff<converge_cretia_ang )

print *, ' converge_log=', converge_log
write(*,*) 'vel_cretia=',converge_cretia
write(*,*) 'ang_cretia=',converge_cretia_ang
write(*,*) 'vx_diff=',vxcentdiff
write(*,*) 'vy_diff=',vycentdiff
write(*,*) 'vz_diff=',vzcentdiff
write(*,*) 'angvx_diff=',angvxdiff
write(*,*) 'angvy_diff=',angvydiff
write(*,*) 'angvz_diff=',angvzdiff

WRITE(*,*)'liu debug in FSConvergeCheck',vycent(iBody),vycent_iter(iBody)
! stop

END SELECT !body_dim

IF (converge_log) THEN
Converged_FSI(iBody) = .True.
ELSE
Converged_FSI(iBody) = .False.
ENDIF

WRITE (*,*) ' Body:', iBody, 'at Ntime =', ntime, ' niterFS =',niterFS
IF (ndim == DIM_2D) THEN
! WRITE (*,*) 'FSI vx/y/centdiff:', vxcentdiff,vycentdiff
WRITE (*,*) ' angvzdiff:', angvzdiff
ENDIF

IF (ndim == DIM_3D) THEN
! WRITE (*,*) ' FSI angvx/y/zdiff:',angvxdiff, angvydiff, angvzdiff
! WRITE (*,*) ' angvx,angvx_iter:', angvx(iBody), angvx_iter(iBody)
! WRITE (*,*) ' angvy,angvy_iter:', angvy(iBody), angvy_iter(iBody)
WRITE (*,*) ' angvzdiff:', angvzdiff
ENDIF

END SUBROUTINE FSConvergeCheck


! Added by Wanh end
! -------------------------------------------------------------------------------------------------------
subroutine angv_noniner2iner(iBody)
! Added by Geng

USE global_parameters
USE flow_parameters
USE usr_module
USE body_dynamics
USE boundary_arrays

integer :: iBody
REAL(KIND=CGREAL) :: vec1(3)



vec1(1) = angvx(ibody)
vec1(2) = angvy(ibody)
vec1(3) = angvz(ibody)


angvx_noniner(iBody)=vec1(1)
angvy_noniner(iBody)=vec1(2)
angvz_noniner(iBody)=vec1(3)

call quaternion_rotation(vec1,quat_prev) ! quat_modify quat_iter -> quat_prev
angvx(ibody) = vec1(1)
angvy(ibody) = vec1(2)
angvz(ibody) = vec1(3)

endsubroutine angv_noniner2iner
!================================================
subroutine quat_update(quat_in,omega_in,frameflag,quat_out)
! update quaternion according to the angular velocity in body frame
! Added by Geng
!---------------------------
! dq/dt=0.5*omega*q
! q(n+1) = q(n)+dq/dt*dt
!---------------------------
use flow_parameters
use operation

REAL(KIND=CGREAL) :: quat_in(4),omega_in(3),quat_out(4)
REAL(KIND=CGREAL) :: q1(4),norm_q
integer :: frameflag ! 1: omega_in is in the non-inertia frame. 2: inerita frame
integer :: m

q1(1)=0.0
do m=2,4
q1(m)=omega_in(m-1)*0.5
enddo

if(frameflag ==1)then
quat_out=qm(quat_in,q1)
elseif(frameflag ==2)then
quat_out=qm(q1,quat_in)
else
write(*,*)'Frameflag should be either 1 or 2, stop!'
stop
endif

do m=1,4
quat_out(m) = quat_out(m)*dt + quat_in(m)
enddo
norm_q=dsqrt(quat_out(1)**2+quat_out(2)**2+quat_out(3)**2+quat_out(4)**2)
quat_out(:)=quat_out(:)/norm_q

end subroutine quat_update
!===========================================================================
! call get_trans_matrix(quat_prev,A_trans)
! call get_transpose(A_trans,AT_trans)
! call get_I_inertia(A_trans,nonDimM_I,AT_trans)
!=================================================================
subroutine get_trans_matrix(quat_trans,A_trans)
use flow_parameters
use operation

REAL(KIND=CGREAL) :: quat_trans(4) ! IN
REAL(KIND=CGREAL) :: A_trans(3,3) ! OUT
REAL(KIND=CGREAL) :: q1(0:3),mtrx(3,3)

integer i,j
do i=1,4
q1(i-1)=quat_trans(i)
enddo

mtrx(1,1) = 1.d0-2.d0*(q1(2)**2.d0+q1(3)**2.d0)
mtrx(1,2) = 2.d0*(q1(1)*q1(2)+q1(0)*q1(3))
mtrx(1,3) = 2.d0*(q1(1)*q1(3)-q1(0)*q1(2))
mtrx(2,1) = 2.d0*(q1(1)*q1(2)-q1(0)*q1(3))
mtrx(2,2) = 1.d0-2.d0*(q1(1)**2.d0+q1(3)**2.d0)
mtrx(2,3) = 2.d0*(q1(3)*q1(2)+q1(0)*q1(1))
mtrx(3,1) = 2.d0*(q1(1)*q1(3)+q1(0)*q1(2))
mtrx(3,2) = 2.d0*(q1(3)*q1(2)-q1(0)*q1(1))
mtrx(3,3) = 1.d0-2.d0*(q1(1)**2.d0+q1(2)**2.d0) ! this is code is from samane's Quaternions.f90

! A_trans = mtrx
call get_transpose(mtrx,A_trans)

endsubroutine
!================================================================
subroutine get_transpose(A_trans,AT_trans)
use flow_parameters
use operation

integer :: i,j
REAL(KIND=CGREAL) :: A_trans(3,3) ! IN
REAL(KIND=CGREAL) :: AT_trans(3,3) ! OUT


do i=1,3
do j=1,3
AT_trans(j,i) = A_trans(i,j)
enddo
enddo

endsubroutine
!================================================================
subroutine get_I_inertia(A_trans,MoI,AT_trans) ! Added by G. Liu I_non_inertia -> I_inertia (I: moment of inertia)
use flow_parameters
use operation

integer :: i,j,m
REAL(KIND=CGREAL) :: A_trans(3,3), AT_trans(3,3) ! IN
REAL(KIND=CGREAL) :: MoI(3,3) ! IN/OUT
REAL(KIND=CGREAL) :: mtemp1(3,3),mtemp2(3,3)

mtemp1=MoI
MoI=zero

call matrix_multi(A_trans,mtemp1,mtemp2)
call matrix_multi(mtemp2,AT_trans,MoI)

endsubroutine
!================================================================
subroutine matrix_multi(a,b,c) ! Added by G. Liu
use flow_parameters
use operation

integer :: n,i,j,m ! n: IN
real(KIND=CGREAL) :: a(3,3),b(3,3) ! IN
REAL(KIND=CGREAL) :: c(3,3) ! OUT

c=0.0
n=3
do i=1,n
do j=1,n
do m=1,n
c(i,j)=c(i,j)+a(i,m)*b(m,j)
enddo
enddo
enddo

endsubroutine
!==============================================================
subroutine get_trans_quat ! Added by G. Liu
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE operation
use body_dynamics

REAL(KIND=CGREAL) :: vec1(3),vec2(3),vecInit1(3),vecInit2(3)
REAL(KIND=CGREAL) :: centnew(3),centInit(3)
REAL(KIND=CGREAL) :: norm(3),normInit(3),axis1(3),axis2(3)
real(kind=cgreal) :: ang1,ang2,ang3
real(kind=cgreal) :: quat(4),qt1(4),qt2(4)


centnew(1) = 0.5*( xBodyMarker(1,rigidRef1(1)) + xBodyMarker(1,rigidRef2(1)) ) ! matched com
centnew(2) = 0.5*( yBodyMarker(1,rigidRef1(1)) + yBodyMarker(1,rigidRef2(1)) )
centnew(3) = 0.5*( zBodyMarker(1,rigidRef1(1)) + zBodyMarker(1,rigidRef2(1)) )

vec1(1) = xBodyMarker(1,rigidRef1(1))-centnew(1)
vec1(2) = yBodyMarker(1,rigidRef1(1))-centnew(2)
vec1(3) = zBodyMarker(1,rigidRef1(1))-centnew(3)

vec2(1) = xBodyMarker(1,rigidRef3(1))-centnew(1)
vec2(2) = yBodyMarker(1,rigidRef3(1))-centnew(2)
vec2(3) = zBodyMarker(1,rigidRef3(1))-centnew(3)

centInit(1) = 0.5*( xBodyMarkerNonIner(1,rigidRef1(1)) + xBodyMarkerNonIner(1,rigidRef2(1)) ) ! matched com
centInit(2) = 0.5*( yBodyMarkerNonIner(1,rigidRef1(1)) + yBodyMarkerNonIner(1,rigidRef2(1)) )
centInit(3) = 0.5*( zBodyMarkerNonIner(1,rigidRef1(1)) + zBodyMarkerNonIner(1,rigidRef2(1)) )

vecInit1(1) = xBodyMarkerNonIner(1,rigidRef1(1))-centInit(1)
vecInit1(2) = yBodyMarkerNonIner(1,rigidRef1(1))-centInit(2)
vecInit1(3) = zBodyMarkerNonIner(1,rigidRef1(1))-centInit(3)

vecInit2(1) = xBodyMarkerNonIner(1,rigidRef3(1))-centInit(1)
vecInit2(2) = yBodyMarkerNonIner(1,rigidRef3(1))-centInit(2)
vecInit2(3) = zBodyMarkerNonIner(1,rigidRef3(1))-centInit(3)

write(*,*) ' debug in get_trans_quat'
write(*,*)vecInit1(:)
write(*,*)vecInit2(:)
write(*,*)vec1(:)
write(*,*)vec2(:)


normInit=cross(vecInit2,vecInit1)
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
else if(ndim==dim_3d)then
if(abs(vecInit2(1)-vec2(1))>zero .or. abs(vecInit2(2)-vec2(2))>zero .or. abs(vecInit2(3)-vec2(3))>zero.or.&
abs(vecInit1(1)-vec1(1))>zero.or.abs(vecInit1(2)-vec1(2))>zero.or.abs(vecInit1(3)-vec1(3))>zero)then
ang1=vector_angle(vecInit2,vec2)
axis1=cross(vecInit2,vec2)/mo(cross(vecInit2,vec2),3)
write(*,*)'ang1,axis1:',ang1,axis1
call quaternion_form(axis1,ang1,qt1)
call quaternion_rotation(normInit,qt1)
if(abs(normInit(1)-norm(1))>zero .or. abs(normInit(2)-norm(2))>zero .or. abs(normInit(3)-norm(3))>zero)then
ang2=vector_angle(normInit,norm)
if(mo(cross(normInit,norm),3) .le. 1.e-15)then
axis2=(/0.0,0.0,1.0/)
else
axis2=cross(normInit,norm)/mo(cross(normInit,norm),3)
endif
write(*,*)'ang2,axis2:',ang2,axis2
call quaternion_form(axis2,ang2,qt2)
else
qt2=(/0.0,0.0,0.0,1.0/)
endif
quat=qm(qt1,qt2)

else
quat=quat_prev
endif
write(*,*)'quat',quat
endif
quat_trans=quat

! stop
endsubroutine
!==============================================================
!==============================================================
subroutine get_quat_defor(iBody) ! Added by G. Liu
! output: quat_trans
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
USE operation
use body_dynamics

REAL(KIND=CGREAL) :: vec1(3),vec2(3),vecInit1(3),vecInit2(3)
REAL(KIND=CGREAL) :: centnew(3),centInit(3)
REAL(KIND=CGREAL) :: norm(3),normInit(3),axis1(3),axis2(3)
real(kind=cgreal) :: ang1,ang2,ang3
real(kind=cgreal) :: quat(4),qt1(4),qt2(4)


! centnew(1) = 0.5*( xBodyMarker(1,rigidRef1(1)) + xBodyMarker(1,rigidRef2(1)) ) ! matched com
! centnew(2) = 0.5*( yBodyMarker(1,rigidRef1(1)) + yBodyMarker(1,rigidRef2(1)) )
! centnew(3) = 0.5*( zBodyMarker(1,rigidRef1(1)) + zBodyMarker(1,rigidRef2(1)) )

vec1(1) = xBodyMarker(iBody,rigidRef2(iBody)) - xBodyMarker(iBody,rigidRef1(iBody))
vec1(2) = yBodyMarker(iBody,rigidRef2(iBody)) - yBodyMarker(iBody,rigidRef1(iBody))
vec1(3) = zBodyMarker(iBody,rigidRef2(iBody)) - zBodyMarker(iBody,rigidRef1(iBody))

vec2(1) = xBodyMarker(iBody,rigidRef3(iBody)) - xBodyMarker(iBody,rigidRef1(iBody))
vec2(2) = yBodyMarker(iBody,rigidRef3(iBody)) - yBodyMarker(iBody,rigidRef1(iBody))
vec2(3) = zBodyMarker(iBody,rigidRef3(iBody)) - zBodyMarker(iBody,rigidRef1(iBody))


vecInit1(1) = xBodyMarkerNonIner(iBody,rigidRef2(iBody)) - xBodyMarkerNonIner(iBody,rigidRef1(iBody))
vecInit1(2) = yBodyMarkerNonIner(iBody,rigidRef2(iBody)) - yBodyMarkerNonIner(iBody,rigidRef1(iBody))
vecInit1(3) = zBodyMarkerNonIner(iBody,rigidRef2(iBody)) - zBodyMarkerNonIner(iBody,rigidRef1(iBody))

vecInit2(1) = xBodyMarkerNonIner(iBody,rigidRef3(iBody)) - xBodyMarkerNonIner(iBody,rigidRef1(iBody))
vecInit2(2) = yBodyMarkerNonIner(iBody,rigidRef3(iBody)) - yBodyMarkerNonIner(iBody,rigidRef1(iBody))
vecInit2(3) = zBodyMarkerNonIner(iBody,rigidRef3(iBody)) - zBodyMarkerNonIner(iBody,rigidRef1(iBody))


write(*,*) ' debug in get_quat_defor'
write(*,*)vecInit1(:)
write(*,*)vecInit2(:)
write(*,*)vec1(:)
write(*,*)vec2(:)


normInit=cross(vecInit2,vecInit1)
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
else if(ndim==dim_3d)then
if(abs(vecInit2(1)-vec2(1))>zero .or. abs(vecInit2(2)-vec2(2))>zero .or. abs(vecInit2(3)-vec2(3))>zero .or. &
abs(vecInit1(1)-vec1(1))>zero .or. abs(vecInit1(2)-vec1(2))>zero .or. abs(vecInit1(3)-vec1(3))>zero)then
ang1=vector_angle(vecInit2,vec2)
axis1=cross(vecInit2,vec2)/mo(cross(vecInit2,vec2),3)
write(*,*)'ang1,axis1:',ang1,axis1
call quaternion_form(axis1,ang1,qt1)
call quaternion_rotation(normInit,qt1)
if(abs(normInit(1)-norm(1))>zero .or. abs(normInit(2)-norm(2))>zero .or. abs(normInit(3)-norm(3))>zero)then
ang2=vector_angle(normInit,norm)
if(mo(cross(normInit,norm),3) .le. 1.e-15)then
axis2=(/0.0,0.0,1.0/)
else
axis2=cross(normInit,norm)/mo(cross(normInit,norm),3)
endif
write(*,*)'ang2,axis2:',ang2,axis2
call quaternion_form(axis2,ang2,qt2)
else
qt2=(/0.0,0.0,0.0,1.0/)
endif
quat=qm(qt1,qt2)

else
quat=(/0.0,0.0,0.0,1.0/)
endif

endif
write(*,*)'quat for coordinate transformation:',quat
quat_trans=quat

! stop
endsubroutine

!======================================================================

subroutine out_MoI_CoM_debug ! Added by G. Liu
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays
use body_dynamics
use usr_module

if(ntime .eq.1 .and. niterFS .eq. 0)then
open(1207,file='mass_I_debug.dat')
else
open(1207,file='mass_I_debug.dat',access='append')
endif

write(1207,1208)ntime,niterFS,volume,nonDimM_I(1,1),nonDimM_I(2,2),nonDimM_I(3,3),nonDimM_I(1,2),nonDimM_I(1,3),nonDimM_I(2,3), &
xcent(1),ycent(1),zcent(1)
1208 format(I5,I4,10(1x,e15.6))
close(1207)

endsubroutine
