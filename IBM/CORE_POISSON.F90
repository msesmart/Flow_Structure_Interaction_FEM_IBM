!--------------------------------------------
! SUBROUTINE rhs_poisson(sum) 
! SUBROUTINE solve_poisson()
! SUBROUTINE itsolv(var,r)
! SUBROUTINE calc_residual(var,r,resm)
!--------------------------------------------
SUBROUTINE find_closest_element(vec,body_dist_min,closest_marker)

USE global_parameters
USE flow_parameters
USE flow_arrays
USE boundary_arrays
USE grid_arrays
USE gcm_arrays
USE unstructured_surface_arrays
USE body_dynamics

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(3) :: Vec
INTEGER :: closest_marker,body_dist_min
INTEGER :: iBody, m
INTEGER, DIMENSION(1) :: iclose
REAL(KIND=CGREAL) :: dist_min,dist_min_body
REAL(KIND=CGREAL) :: distX, distY, distZ
REAL(KIND=CGREAL) :: dist_sqr(nPtsMax)

dist_min = 1.0E10_CGREAL
body_dist_min = 0

DO iBody = 1, nBody

SELECT CASE(ndim)
CASE(DIM_2D)
DO m=1,nPtsBodyMarker(iBody)
distX = (vec(1)-xBodyMarker(iBody,m))
distY = (vec(2)-yBodyMarker(iBody,m))
dist_sqr(m) = distX**2 + distY**2
ENDDO ! m

CASE(DIM_3D)
DO m=1,nPtsBodyMarker(iBody)
distX = (vec(1)-xBodyMarker(iBody,m))
distY = (vec(2)-yBodyMarker(iBody,m))
distZ = (vec(3)-zBodyMarker(iBody,m))
dist_sqr(m) = distX**2 + distY**2 + distZ**2

ENDDO ! m

END SELECT ! canonical_body_type

dist_min_body = MINVAL(dist_sqr(1:nPtsBodyMarker(iBody)))
IF ( dist_min_body <= dist_min ) THEN
dist_min = dist_min_body
iclose = MINLOC(dist_sqr(1:nPtsBodyMarker(iBody)))
closest_marker = iclose(1)
body_dist_min = iBody
ENDIF ! dist_min_body

ENDDO ! iBody

END SUBROUTINE find_closest_element
!-------------------------------------------------------------------------------
SUBROUTINE find_closest_element_modified(vec,body_dist_min,closest_marker,body_start,body_end)

USE global_parameters
USE flow_parameters
USE flow_arrays
USE boundary_arrays
USE grid_arrays
USE gcm_arrays
USE unstructured_surface_arrays
USE body_dynamics

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(3) :: Vec
INTEGER :: closest_marker,body_dist_min
INTEGER :: iBody, m, body_start,body_end
INTEGER, DIMENSION(1) :: iclose
REAL(KIND=CGREAL) :: dist_min,dist_min_body
REAL(KIND=CGREAL) :: distX, distY, distZ
REAL(KIND=CGREAL) :: dist_sqr(nPtsMax)

dist_min = 1.0E10_CGREAL
body_dist_min = 0

DO iBody = body_start, body_end

SELECT CASE(ndim)
CASE(DIM_2D)
DO m=1,nPtsBodyMarker(iBody)
distX = (vec(1)-xBodyMarker(iBody,m))
distY = (vec(2)-yBodyMarker(iBody,m))
dist_sqr(m) = distX**2 + distY**2
ENDDO ! m

CASE(DIM_3D)
DO m=1,nPtsBodyMarker(iBody)
distX = (vec(1)-xBodyMarker(iBody,m))
distY = (vec(2)-yBodyMarker(iBody,m))
distZ = (vec(3)-zBodyMarker(iBody,m))
dist_sqr(m) = distX**2 + distY**2 + distZ**2

ENDDO ! m

END SELECT ! canonical_body_type

dist_min_body = MINVAL(dist_sqr(1:nPtsBodyMarker(iBody)))
IF ( dist_min_body <= dist_min ) THEN
dist_min = dist_min_body
iclose = MINLOC(dist_sqr(1:nPtsBodyMarker(iBody)))
closest_marker = iclose(1)
body_dist_min = iBody
ENDIF ! dist_min_body

ENDDO ! iBody

END SUBROUTINE find_closest_element_modified
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

SUBROUTINE rhs_poisson(sum) 

USE global_parameters
USE flow_parameters
USE flow_arrays
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays

IMPLICIT NONE

REAL(KIND=CGREAL), INTENT(OUT) :: sum

REAL(KIND=CGREAL) :: sum_rhs_ghost
INTEGER :: i,j,k
REAL(KIND=CGREAL) :: AAe, AAw, AAn, AAs, AAf, AAb, AAg, Ag
REAL(KIND=CGREAL) :: Be, Bw, Bn, Bs, Bf, Bb, B
REAL(KIND=CGREAL) :: face_ue, face_uw, face_vn, face_vs, face_wf, face_wb
INTEGER :: closest_marker,body_dist_min
REAL(KIND=CGREAL), DIMENSION(3) :: Ub
! REAL(KIND=CGREAL) :: taylor_series_matching

!******************************************************************************

! rhs = [ d(U)/dx + d(V)/dy + d(W)/dz ] / dt

sum = zero
sum_rhs_ghost = zero
! nlu0(0:,0:,0:) = zero ! debug nlu
! nlw(0:,0:,0:) = zero ! debug nlu
nlu0(0:,0:,0:) = zero ! debug nlu
nlw0(0:,0:,0:) = zero ! debug nlu

DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
! nlu0(i,j,k) = (dtinv)* &
IF (cure_pressure_oscillations .AND. ivc(i,j,k)>0) THEN
AAe=face(cell(ivc(i,j,k))%F_ip)%a
AAw=face(cell(ivc(i,j,k))%F_im)%a
AAn=face(cell(ivc(i,j,k))%F_jp)%a
AAs=face(cell(ivc(i,j,k))%F_jm)%a
AAf=face(cell(ivc(i,j,k))%F_kp)%a
AAb=face(cell(ivc(i,j,k))%F_km)%a

!! VOF differential form
! AAe=oned
! AAw=oned
! AAn=oned
! AAs=oned
! AAf=oned
! AAb=oned
! doesn't work :(

AAg=face(cell(ivc(i,j,k))%F_slice)%a

call find_closest_element(face(cell(ivc(i,j,k))%F_slice)%centroid, body_dist_min, closest_marker)

Ub(1)=uBodyMarker(body_dist_min,closest_marker)
Ub(2)=vBodyMarker(body_dist_min,closest_marker)
Ub(3)=wBodyMarker(body_dist_min,closest_marker)

! IF (iblank(i,j,k)==0) THEN
! face_ue=face_u(i+1,j,k)
! face_uw=face_u(i,j,k)
! face_vn=face_v(i,j+1,k)
! face_vs=face_v(i,j,k)
! face_wf=face_w(i,j,k+1)
! face_wb=face_w(i,j,k)
! ELSE
! IF (iblank(i+1,j,k)==0) then
! face_ue=face_u(i+1,j,k)
! ELSE
! face_ue=Ub%x(1)
! ENDIF
! 
! IF (iblank(i-1,j,k)==0) then
! face_uw=face_u(i,j,k)
! ELSE
! face_uw=Ub%x(1)
! ENDIF
! 
! IF (iblank(i,j+1,k)==0) then
! face_vn=face_v(i,j+1,k)
! ELSE
! face_vn=Ub%x(2)
! ENDIF
! 
! IF (iblank(i,j-1,k)==0) then
! face_vs=face_v(i,j,k)
! ELSE 
! face_vs=Ub%x(2)
! ENDIF
! 
! IF (ndim==DIM_3D) THEN
! 
! IF (iblank(i,j,k+1)==0) then
! face_wf=face_w(i,j,k+1)
! ELSE
! face_wf=Ub%x(3)
! ENDIF
! IF (iblank(i,j,k-1)==0) then
! face_wb=face_w(i,j,k)
! ELSE
! face_wb=Ub%x(3)
! ENDIF
! ELSE
! face_wf=Ub%x(3)
! face_wb=Ub%x(3)
! ENDIF 
! ENDIF

IF (cell(ivc(i,j,k))%F_ip>0) THEN
face_ue=face(cell(ivc(i,j,k))%F_ip)%vel
ELSE
face_ue=face_u(i+1,j,k)
ENDIF
IF (cell(ivc(i,j,k))%F_im>0) THEN
face_uw=face(cell(ivc(i,j,k))%F_im)%vel
ELSE
face_uw=face_u(i,j,k)
ENDIF
IF (cell(ivc(i,j,k))%F_jp>0) THEN
face_vn=face(cell(ivc(i,j,k))%F_jp)%vel
ELSE
face_vn=face_v(i,j+1,k)
ENDIF
IF (cell(ivc(i,j,k))%F_jm>0) THEN
face_vs=face(cell(ivc(i,j,k))%F_jm)%vel
ELSE
face_vs=face_v(i,j,k)
ENDIF
IF (ndim==DIM_3D) THEN
IF (cell(ivc(i,j,k))%F_kp>0) THEN
face_wf=face(cell(ivc(i,j,k))%F_kp)%vel
ELSE
face_wf=face_w(i,j,k+1)
ENDIF
IF (cell(ivc(i,j,k))%F_km>0) THEN
face_wb=face(cell(ivc(i,j,k))%F_km)%vel
ELSE
face_wb=face_w(i,j,k)
ENDIF
ELSE
face_wf=face_w(i,j,k+1)
face_wb=face_w(i,j,k)
ENDIF

nlu0(i,j,k) = (dtinv)* &
( ( face_ue*AAe - face_uw*AAw )*dxinv(i) &
+( face_vn*AAn - face_vs*AAs )*dyinv(j) &
+( face_wf*AAf - face_wb*AAb )*dzinv(k) &
+DOT_PRODUCT(Ub, cell(ivc(i,j,k))%slice_normal)*AAg)
write(1001,*) i, j, IBLANK(I,J,K), nlu0(i,j,k), (dtinv)* &
( ( face_u(i+1,j,k) - face_u(i,j,k) )*dxinv(i) &
+( face_v(i,j+1,k) - face_v(i,j,k) )*dyinv(j) &
+( face_w(i,j,k+1) - face_w(i,j,k) )*dzinv(k) )

ELSE
nlu0(i,j,k) = (dtinv)* &
( ( face_u(i+1,j,k) - face_u(i,j,k) )*dxinv(i) &
+( face_v(i,j+1,k) - face_v(i,j,k) )*dyinv(j) &
+( face_w(i,j,k+1) - face_w(i,j,k) )*dzinv(k) )
nlu0(i,j,k) = nlu0(i,j,k)*REAL(1-iblank(i,j,k),KIND=CGREAL)
ENDIF
sum = sum + nlu0(i,j,k)*dx(i)*dy(j)*dz(k)
ENDDO
ENDDO
ENDDO

! CALL write_dump_debug('U ',0,U)
! CALL write_dump_debug('V ',0,V)
! CALL write_dump_debug('bcxu',0,bcxu)
! CALL write_dump_debug('bcyu',0,bcyu)
! CALL write_dump_debug('bcxv',0,bcxv)
! CALL write_dump_debug('bcyv',0,bcyv)
! CALL write_dump_debug('fU ',0,face_U)
! CALL write_dump_debug('fV ',0,face_V)
! CALL write_dump_debug('fW ',0,face_W)
! CALL write_dump_debug('nlu ',1,nlu)

! CALL write_dump()

IF (cure_pressure_oscillations) THEN
DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
IF (iblank(i,j,k)==1 .AND. ivc(i,j,k)>0) THEN

IF (iblank(i+1,j,k)==0 .AND. ivc(i+1,j,k)>0) THEN
IF (cell(ivc(i+1,j,k))%volumn>half) THEN
Be=cell(ivc(i,j,k))%slice_normal(1)**2*face(cell(ivc(i,j,k))%F_ip)%a
ENDIF
ELSE
Be=zero
ENDIF
IF (iblank(i-1,j,k)==0 .AND. ivc(i-1,j,k)>0) THEN
IF (cell(ivc(i-1,j,k))%volumn>half) THEN
Bw=cell(ivc(i,j,k))%slice_normal(1)**2*face(cell(ivc(i,j,k))%F_im)%a
ENDIF
ELSE
Bw=zero
ENDIF

IF (iblank(i,j+1,k)==0 .AND. ivc(i,j+1,k)>0) THEN
IF (cell(ivc(i,j+1,k))%volumn>half) THEN
Bn=cell(ivc(i,j,k))%slice_normal(2)**2*face(cell(ivc(i,j,k))%F_jp)%a
ENDIF
ELSE
Bn=zero
ENDIF
IF (iblank(i,j-1,k)==0 .AND. ivc(i,j-1,k)>0) THEN
IF (cell(ivc(i,j-1,k))%volumn>half) THEN
Bs=cell(ivc(i,j,k))%slice_normal(2)**2*face(cell(ivc(i,j,k))%F_jm)%a
ENDIF
ELSE
Bs=zero
ENDIF

IF (iblank(i,j,k+1)==0 .AND. ivc(i,j,k+1)>0) THEN
IF (cell(ivc(i,j,k+1))%volumn>half) THEN
Bf=cell(ivc(i,j,k))%slice_normal(3)**2*face(cell(ivc(i,j,k))%F_kp)%a
ENDIF
ELSE
Bf=zero
ENDIF
IF (iblank(i,j,k-1)==0 .AND. ivc(i,j,k-1)>0) THEN
IF (cell(ivc(i,j,k-1))%volumn>half) THEN
Bb=cell(ivc(i,j,k))%slice_normal(3)**2*face(cell(ivc(i,j,k))%F_km)%a
ENDIF
ELSE
Bb=zero
ENDIF

B = Be + Bw + Bn + Bs + Bf + Bb

IF (iblank(i+1,j,k)==0 .AND. ivc(i+1,j,k)>0) THEN
IF (cell(ivc(i+1,j,k))%volumn>half) THEN
nlu0(i+1,j,k) = nlu0(i+1,j,k) + Be / B * nlu0(i,j,k)
ENDIF
ENDIF
IF (iblank(i-1,j,k)==0 .AND. ivc(i-1,j,k)>0) THEN
IF (cell(ivc(i-1,j,k))%volumn>half) THEN
nlu0(i-1,j,k) = nlu0(i-1,j,k) + Bw / B * nlu0(i,j,k)
ENDIF
ENDIF

IF (iblank(i,j+1,k)==0 .AND. ivc(i,j+1,k)>0) THEN
IF (cell(ivc(i,j+1,k))%volumn>half) THEN
nlu0(i,j+1,k) = nlu0(i,j+1,k) + Bn / B * nlu0(i,j,k)
ENDIF
ENDIF
IF (iblank(i,j-1,k)==0 .AND. ivc(i,j-1,k)>0) THEN
IF (cell(ivc(i,j-1,k))%volumn>half) THEN
nlu0(i,j-1,k) = nlu0(i,j-1,k) + Bs / B * nlu0(i,j,k)
ENDIF
ENDIF

IF (iblank(i,j,k+1)==0 .AND. ivc(i,j,k+1)>0) THEN
IF (cell(ivc(i,j,k+1))%volumn>half) THEN
nlu0(i,j,k+1) = nlu0(i,j,k+1) + Bf / B * nlu0(i,j,k)
ENDIF
ENDIF
IF (iblank(i,j,k-1)==0 .AND. ivc(i,j,k-1)>0) THEN
IF (cell(ivc(i,j,k-1))%volumn>half) THEN
nlu0(i,j,k-1) = nlu0(i,j,k-1) + Bb / B * nlu0(i,j,k)
ENDIF
ENDIF

nlu0(i,j,k) = zero

ENDIF
ENDDO
ENDDO
ENDDO
ENDIF

print *, 'in rhs poisson, max/min nlu:', maxval(nlu0),minval(nlu0)

!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
! DO k = 1,nzc
! DO j = 1,nyc
! DO i = 1,nxc
! if (ghostcellmark(i,j,k) == 1) then
! sum_rhs_ghost = sum_rhs_ghost + &
! ( ( face_u(i+1,j,k) - face_u(i,j,k) )*dxinv(i) &
! +( face_v(i,j+1,k) - face_v(i,j,k) )*dyinv(j) &
! +( face_w(i,j,k+1) - face_w(i,j,k) )*dzinv(k) )*dx(i)*dy(j)*dz(k)
! endif
! ENDDO
! ENDDO
! ENDDO
! write(863,*)'sum_rhs_ghost =',sum_rhs_ghost
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
END SUBROUTINE rhs_poisson
!-------------------------------------------------------------------------------
SUBROUTINE subtracting_average_pressure()
USE global_parameters
USE flow_parameters
USE grid_arrays
USE pressure_arrays
USE boundary_arrays

IMPLICIT NONE

INTEGER :: i,j,k, ii, jj, kk, i1, k1, k2
REAL(KIND=CGREAL) :: p1
REAL(KIND=CGREAL) :: pTotal, volTotal, pAverage

!------- new part for subtracting average pressure ----------
pTotal = zero
volTotal = zero
DO K = 1, nzc
DO J = 1, nyc
DO I = 1, nxc
pTotal = pTotal + pPrime(i,j,k)*dx(i)*dy(j)*dz(k)*(1-iblank(i,j,k))
volTotal = volTotal + dx(i)*dy(j)*dz(k)*(1-iblank(i,j,k))
END DO
END DO
END DO
pAverage = pTotal/volTotal

WRITE(*,'(A,1PE15.7,A,1PE15.7)') 'Subtract average pressure ',pAverage,' in volumn', volTotal

DO K = 1, nzc
DO J = 1, nyc
DO I = 1, nxc
IF (iblank(i,j,k)==1) THEN
pPrime(i,j,k) = zero ! For dead cells
ELSE
pPrime(i,j,k) = pPrime(i,j,k) - pAverage
ENDIF
END DO
END DO
END DO

DO K = 1, nz-1
DO J = 1, ny-1
DO I = 1, nx-1
IF (fresh_cell(i,j,k)==1) THEN
i1=0
p1=zero
if (nDim .EQ. DIM_3D) then
k1=k-1
k2=k+1
else
k1=k
k2=k
endif

do kk=k1, k2
do jj=j-1, j+1
do ii=i-1, i+1
if (iblank(ii,jj,kk)/=1 .and. fresh_cell(ii,jj,kk)/=1) then
i1=i1+1
p1=p1+pPrime(ii,jj,kk)
endif
END DO
END DO
END DO

pPrime(i,j,k) = p1/i1 ! interpolate for fresh cells
ENDIF
END DO
END DO
END DO

END SUBROUTINE subtracting_average_pressure

!-------------------------------------------------------------------------------
SUBROUTINE solve_poisson()

USE global_parameters
USE flow_parameters
USE flow_arrays
USE pressure_arrays
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays
USE mg_parameters
USE mg_ARRAYS
USE fea_unstructure_surface, ONLY : FSI_ITERATION
USE Pressure_Aitken_Array

IMPLICIT NONE

!... Local variables
EXTERNAL MG_itsolv
EXTERNAL MG_itsolv_Point_Jacobi
EXTERNAL MG_itsolv_MSIP_2D, MG_itsolv_MSIP_3D

EXTERNAL MG_Residual
EXTERNAL MG_Residual_MSIP

EXTERNAL MG_Precondition_MSIP_2D, MG_Precondition_MSIP_3D
EXTERNAL MG_Precondition_MSIP_Conformal_Mapping

INTEGER :: iter,i,j,k,loc(3)
INTEGER :: iErr
REAL(KIND=CGREAL) :: maxres
REAL(KIND=CGREAL) :: compTime

!******************************************************************************

IF (pbcx1 == PBC_DIRICHLET .OR. &
pbcy1 == PBC_DIRICHLET .OR. &
pbcz1 == PBC_DIRICHLET) THEN

CALL set_pressure_dirichlet_bc()

ELSE

CALL subtracting_average_pressure()

END IF

!---------------------------------------------------
!--------------------------------------------------------------------
! Remove the ups and ums for periodic condition in certain direction
! This only works for LSOR and MG method, not for PETSC method
!-------------------------------------------------------------------
IF (pp_solver_type .NE. PP_SOLVER_TYPE_PETSC) THEN
IF (bcx1 == BC_TYPE_PERIODIC .OR. &
bcy1 == BC_TYPE_PERIODIC .OR. &
bcz1 == BC_TYPE_PERIODIC) THEN
CALL remove_up_um
END IF
END IF

iter = 0
maxres = 1.0E10_CGREAL 

IF (pressureAitkenOn) THEN !Added by Wanh for Aitken acceleration
pPrime_L(:,:,:) = pPrime(:,:,:)
ENDIF

IF (MOD(ntime,nmonitor) == 0) THEN
CALL calc_residual(pPrime,nlu0,maxres,loc)
WRITE(STDOUT,'(A,2(2X,1PE15.7),3(2X,I3))') &
'Initial Residual Entering Poisson Solver', &
maxres,restol_Poisson,loc(1),loc(2),loc(3)
END IF ! ntime 

SELECT CASE(pp_solver_type)

CASE (PP_SOLVER_TYPE_LSOR)
DO WHILE ((iter .LT. iterMax_Poisson) .AND. (maxres .GT. restol_Poisson))
CALL itsolv(pPrime,nlu0)
IF (boundary_formulation == GCM_METHOD) CALL GCM_correct_pressure()
IF (bcx1 == BC_TYPE_PERIODIC .OR. & 
bcy1 == BC_TYPE_PERIODIC .OR. &
bcz1 == BC_TYPE_PERIODIC) THEN 
CALL enforce_p_periodic(pPrime)
END IF
CALL calc_residual(pPrime,nlu0,maxres,loc)
iter = iter + 1
WRITE(STDOUT,'(A,I4,2X,1PE15.7,3(2X,I4))')&
'LSOR: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3)
ENDDO
IF ( iter .EQ. iterMax_Poisson .AND. maxres .GT. restol_Poisson ) THEN
PRINT*,'LSOR: Pressure did not converge in ',iterMax_Poisson,' iterations'
PRINT*,'Final residual = ',maxres,' at',loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
ELSE
IF (MOD(ntime,nmonitor) == 0) THEN
WRITE(STDOUT,'(A,I4,2X,1PE15.7,4(2X,I4))')&
'LSOR: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
END IF ! ntime
ENDIF

CASE (PP_SOLVER_TYPE_PETSC)

#ifdef PETSC

! --- Enter PETSc Solver
! For inactive pressure compatibility, apply Usual Iterative process
! For active pressure compatibility, exit PETSC solver at each iteration

IF ( boundary_formulation == NONE .OR. &
boundary_formulation == SSM_METHOD ) THEN
CALL petsc_solver(pPrime,nlu0)

ELSE IF ( boundary_formulation == GCM_METHOD ) THEN
DO WHILE ((iter .LT. iterMax_Poisson) .AND. (maxres .GT. restol_Poisson)) 
CALL petsc_solver(pPrime,nlu0)
CALL GCM_correct_pressure()

iter = iter + 1
CALL calc_residual(pPrime,nlu0,maxres,loc)
! PRINT*,'PETSC: Iter, Maxres = ',iter,maxres,maxres*dt
ENDDO ! iter
ENDIF ! boundary_formulation

IF ( iter .EQ. iterMax_Poisson .AND. maxres .GT. restol_Poisson ) THEN
PRINT*,'PETSC: Pressure Convergence at NTIME = ',ntime
PRINT*,'PETSC: Pressure did not converge in ',iterMax_Poisson,' iterations'
PRINT*,'Final residual = ',maxres,' at',loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
ELSE
IF (MOD(ntime,nmonitor) == 0) THEN
WRITE(STDOUT,'(A,I4,2X,1PE15.7,4(2X,I4))')&
'PETSC: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
END IF ! ntime
ENDIF

#else
PRINT*,'USER ERROR: PETSC Solver Active in Input Deck '
PRINT*,' Code not compiled with PETSC Flag'
PRINT*,' Code will stop and exit '
PRINT*,' Either set solver to MG or compile with PETSC=1'
STOP
#endif

CASE (PP_SOLVER_TYPE_MG, PP_SOLVER_TYPE_MG_Point_Jacobi) 

DO WHILE ((iter .LT. iterMax_Poisson) .AND. (maxres .GT. restol_Poisson))

IF (pp_solver_type==PP_SOLVER_TYPE_MG) THEN
CALL mg_solver(pPrime,nlu0,MG_itsolv,MG_Residual,compTime)
ELSE
CALL mg_solver(pPrime,nlu0,MG_itsolv_Point_Jacobi,MG_Residual,compTime)
ENDIF
IF (bcx1 == BC_TYPE_PERIODIC .OR. & 
bcy1 == BC_TYPE_PERIODIC .OR. &
bcz1 == BC_TYPE_PERIODIC) THEN 
CALL enforce_p_periodic(pPrime)
END IF

IF (boundary_formulation == GCM_METHOD) CALL GCM_correct_pressure()

CALL calc_residual(pPrime,nlu0,maxres,loc)
iter = iter + 1
IF (MOD(ntime,nmonitor) == 0) &
WRITE(STDOUT,'(A,I4,2X,1PE15.7,3(2X,I4),2X,F7.3)')&
'MG-LSOR: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3),compTime

IF (pressureAitkenOn) CALL Pressure_Aitken(loc,iter)

ENDDO

IF ( iter .EQ. iterMax_Poisson .AND. maxres .GT. restol_Poisson ) THEN
PRINT*,'MG-LSOR: Pressure Convergence at NTIME = ',ntime
PRINT*,'MG-LSOR: Pressure did not converge in ',iterMax_Poisson,' iterations'
PRINT*,'Final residual = ',maxres,' at',loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
ELSE
IF (MOD(ntime,nmonitor) == 0) THEN
WRITE(STDOUT,'(A,I4,2X,1PE15.7,4(2X,I4))')&
'MG-LSOR: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
END IF ! ntime
ENDIF

!******************************************************************************
CASE (PP_SOLVER_TYPE_SIP) 
IF (ndim == DIM_2D) THEN
CALL COEFF_PRESSURE_MSIP_2D()
CALL PRECONDITION_SIP_2D()
ELSE
WRITE(STDOUT,*) 'SIP 3D is not available.'
WRITE(STDOUT,*) 'Please select other algorithms.'
STOP
END IF

DO WHILE ((iter .LT. iterMax_Poisson) .AND. (maxres .GT. restol_Poisson))
IF (ndim == DIM_3D) THEN
CALL itsolv_SIP_3D(pPrime,nlu0)
ELSE
CALL itsolv_SIP_2D(pPrime,nlu0)
END IF

! IF (pbcx1 == PBC_NEUMANN .AND. pbcx2 == PBC_NEUMANN .AND. &
! pbcy1 == PBC_NEUMANN .AND. pbcy2 == PBC_NEUMANN .AND. &
! pbcz1 == PBC_NEUMANN .AND. pbcz2 == PBC_NEUMANN) &
! CALL subtracting_average_pressure()

IF (boundary_formulation == GCM_METHOD) CALL GCM_correct_pressure()

CALL calc_residual(pPrime,nlu0,maxres,loc)
iter = iter + 1 
WRITE(STDOUT,'(A,I4,2X,1PE15.7,3(2X,I4))')&
'SIP: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3)
ENDDO

DEALLOCATE(LU,CA)

IF ( iter .EQ. iterMax_Poisson .AND. maxres .GT. restol_Poisson ) THEN
PRINT*,'SIP: Pressure did not converge in ',iterMax_Poisson,' iterations'
PRINT*,'Final residual = ',maxres,' at',loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
ELSE
IF (MOD(ntime,nmonitor) == 0) THEN
WRITE(STDOUT,'(A,I4,2X,1PE15.7,4(2X,I4))')&
'SIP: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
END IF ! ntime
ENDIF

!******************************************************************************
CASE (PP_SOLVER_TYPE_MSIP) 
IF (ndim == DIM_2D) THEN
CALL COEFF_PRESSURE_MSIP_2D()
CALL COEFF_LU_MSIP_2D()
ELSE
CALL COEFF_PRESSURE_MSIP_3D()
CALL COEFF_LU_MSIP_3D()
END IF

DO WHILE ((iter .LT. iterMax_Poisson) .AND. (maxres .GT. restol_Poisson))
IF (ndim == DIM_3D) THEN
CALL itsolv_MSIP_3D(pPrime,nlu0,restol_Poisson)
ELSE
CALL itsolv_MSIP_2D(pPrime,nlu0,restol_Poisson)
END IF

! IF (pbcx1 == PBC_NEUMANN .AND. pbcx2 == PBC_NEUMANN .AND. &
! pbcy1 == PBC_NEUMANN .AND. pbcy2 == PBC_NEUMANN .AND. &
! pbcz1 == PBC_NEUMANN .AND. pbcz2 == PBC_NEUMANN) &
! CALL subtracting_average_pressure()

IF (boundary_formulation == GCM_METHOD) CALL GCM_correct_pressure()

CALL calc_residual(pPrime,nlu0,maxres,loc)
iter = iter + 1 
WRITE(STDOUT,'(A,I4,2X,1PE15.7,3(2X,I4))')&
'MSIP: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3)
ENDDO

DEALLOCATE(LU,CA)

IF ( iter .EQ. iterMax_Poisson .AND. maxres .GT. restol_Poisson ) THEN
PRINT*,'MSIP: Pressure did not converge in ',iterMax_Poisson,' iterations'
PRINT*,'Final residual = ',maxres,' at',loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
ELSE
IF (MOD(ntime,nmonitor) == 0) THEN
WRITE(STDOUT,'(A,I4,2X,1PE15.7,4(2X,I4))')&
'MSIP: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
END IF ! ntime
ENDIF

!******************************************************************************
CASE (PP_SOLVER_TYPE_MG_MSIP) 
CALL MG_Memory_Allocation_iblank(mgLevels_X,ICOORD)

IF (Full_Coarsening) THEN

CALL MG_Prepare_Iblank_FMG

ELSE

CALL MG_Memory_Allocation_iblank(mgLevels_Y,JCOORD)
! IF (Conformal_Mapping) THEN
! CALL MG_Prepare_Iblank_Conformal_Mapping(ICOORD)
! CALL MG_Prepare_Iblank_Conformal_Mapping(JCOORD)
! ELSE
CALL MG_Prepare_Iblank(ICOORD)
CALL MG_Prepare_Iblank(JCOORD)
! ENDIF
IF (ndim == DIM_3D) THEN
CALL MG_Memory_Allocation_iblank(mgLevels_Z,KCOORD)
! IF (Conformal_Mapping) THEN
! CALL MG_Prepare_Iblank_Conformal_Mapping(KCOORD)
! ELSE
CALL MG_Prepare_Iblank(KCOORD)
! ENDIF
END IF

ENDIF

! IF (Conformal_Mapping) THEN
! CALL MG_PREPARE_MSIP(MG_Precondition_MSIP_Conformal_Mapping)
! ELSE
IF (ndim == DIM_2D) THEN
CALL MG_PREPARE_MSIP(MG_Precondition_MSIP_2D)
ELSE
CALL MG_PREPARE_MSIP(MG_Precondition_MSIP_3D)
ENDIF
! ENDIF

DO WHILE ((iter .LT. iterMax_Poisson) .AND. (maxres .GT. restol_Poisson))

IF (Full_Coarsening) THEN

IF (ndim == DIM_2D) THEN
CALL mg_solver_FMG(pPrime,nlu0,MG_itsolv_MSIP_2D,MG_Residual_MSIP,compTime)
ELSE
CALL mg_solver_FMG(pPrime,nlu0,MG_itsolv_MSIP_3D,MG_Residual_MSIP,compTime)
ENDIF

ELSE
IF (ndim == DIM_2D) THEN
CALL mg_solver(pPrime,nlu0,MG_itsolv_MSIP_2D,MG_Residual_MSIP,compTime)
ELSE
CALL mg_solver(pPrime,nlu0,MG_itsolv_MSIP_3D,MG_Residual_MSIP,compTime)
ENDIF

ENDIF 
! IF (pbcx1 == PBC_NEUMANN .AND. pbcx2 == PBC_NEUMANN .AND. &
! pbcy1 == PBC_NEUMANN .AND. pbcy2 == PBC_NEUMANN .AND. &
! pbcz1 == PBC_NEUMANN .AND. pbcz2 == PBC_NEUMANN) &
! CALL subtracting_average_pressure()

IF (boundary_formulation == GCM_METHOD) CALL GCM_correct_pressure()

CALL calc_residual_MG_MSIP(pPrime,nlu0,maxres,loc)
iter = iter + 1 
WRITE(STDOUT,'(A,I4,2X,1PE15.7,3(2X,I4),2X,1F8.3,1X)')&
'MG-MSIP: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3),compTime
! WRITE(STDOUT,'(2X,1F8.3)') compTime

IF (pressureAitkenOn) CALL Pressure_Aitken(loc,iter)
! STOP

ENDDO

CALL MG_Memory_Deallocation_iblank(mgLevels_X,ICOORD)
DO i=1, mgLevels_X
DEALLOCATE(MGX(i)%CA, MGX(i)%LU)
ENDDO

IF (Full_Coarsening) THEN
ELSE

CALL MG_Memory_Deallocation_iblank(mgLevels_Y,JCOORD)
DO i=2, mgLevels_Y
DEALLOCATE(MGY(i)%CA, MGY(i)%LU)
ENDDO
IF (ndim == DIM_3D) THEN
CALL MG_Memory_Deallocation_iblank(mgLevels_Z,KCOORD)
DO i=2, mgLevels_Z
DEALLOCATE(MGZ(i)%CA, MGZ(i)%LU)
ENDDO
END IF

ENDIF

IF ( iter .EQ. iterMax_Poisson .AND. maxres .GT. restol_Poisson ) THEN
PRINT*,'MG-MSIP: Pressure did not converge in ',iterMax_Poisson,' iterations'
PRINT*,'Final residual = ',maxres,' at',loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
ELSE
IF (MOD(ntime,nmonitor) == 0) THEN
WRITE(STDOUT,'(A,I4,2X,1PE15.7,4(2X,I4))')&
'MG-MSIP: Pressure Convergence : ',iter,maxres,loc(1),loc(2),loc(3),&
fresh_cell(loc(1),loc(2),loc(3))
END IF ! ntime
ENDIF

END SELECT

!--------------------------------------------------------------------
! Add the ups and ums back for periodic condition in certain direction
! This only works for LSOR and MG method, not for PETSC method
!-------------------------------------------------------------------
IF (pp_solver_type .NE. PP_SOLVER_TYPE_PETSC) THEN
IF (bcx1 == BC_TYPE_PERIODIC .OR. & 
bcy1 == BC_TYPE_PERIODIC .OR. &
bcz1 == BC_TYPE_PERIODIC) THEN 
CALL add_up_um 
END IF
END IF

IF (pbcx1 == PBC_DIRICHLET .OR. & !
pbcy1 == PBC_DIRICHLET .OR. & ! Added by H. Luo
pbcz1 == PBC_DIRICHLET) THEN !

CALL unset_pressure_dirichlet_bc()
END IF

END SUBROUTINE solve_poisson

!----------------------------------------------------------
! Line SOR with Gauss Siedel as smoother
!----------------------------------------------------------
SUBROUTINE itsolv(var,r)

USE global_parameters
USE flow_parameters
USE flow_arrays
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays
USE solver_arrays
USE GCM_arrays

IMPLICIT NONE

!... parameters

REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT (IN OUT) ::var
REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT (IN) ::r

!... Local variables

REAL(KIND=CGREAL) :: VolCell
INTEGER :: i,j,k
INTEGER :: iBody, iRow, iG,jG, iNode, jNode, n
INTEGER :: loc(3)
REAL(KIND=CGREAL) :: maxres 


!******************************************************************************
! Note.. homogeneous Neumann BC have been hardcoded here

CALL enforce_p_periodic(var)

! Line solve in the x-direction
DO k=1,nzc
DO j=1,nyc
DO i=1,nxc
amx(i) = dxcinv(i) *dxinv(i)*(1 - ium(i,j,k) )
apx(i) = dxcinv(i+1)*dxinv(i)*(1 - iup(i,j,k) )
acx(i) = - ( amx(i) + apx(i) )

amy(j) = dycinv(j) *dyinv(j)*(1 - jum(i,j,k) )
apy(j) = dycinv(j+1)*dyinv(j)*(1 - jup(i,j,k) )
acy(j) = - ( amy(j) + apy(j) )

amz(k) = dzcinv(k) *dzinv(k)*(1 - kum(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL)
apz(k) = dzcinv(k+1)*dzinv(k)*(1 - kup(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL)
acz(k) = - ( amz(k) + apz(k) )


rhs(i) = r(i,j,k) - var(i,j-1,k)*amy(j) &
- var(i,j+1,k)*apy(j) &
- var(i,j,k-1)*amz(k) &
- var(i,j,k+1)*apz(k) &
+ nlw0(i,j,k)

amx(i) = amx(i)*REAL(1-iblank(i,j,k),KIND=CGREAL)
apx(i) = apx(i)*REAL(1-iblank(i,j,k),KIND=CGREAL)
acx(i) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank(i,j,k),KIND=CGREAL) &
+ REAL(iblank(i,j,k),KIND=CGREAL) 
rhs(i) = rhs(i)*REAL(1-iblank(i,j,k),KIND=CGREAL) &
+ REAL(ghostcellMark(i,j,k),KIND=CGREAL)*var(i,j,k)

ENDDO ! i 

IF (bcx1 == BC_TYPE_PERIODIC .AND. bcx2 == BC_TYPE_PERIODIC) THEN
rhs(1) = rhs(1) - var(nxc,j,k)*amx(1)
rhs(nxc) = rhs(nxc) - var(1,j,k)*apx(nxc)
END IF

IF (pbcx1 == PBC_DIRICHLET .and. pbcx2 == PBC_DIRICHLET) THEN !
rhs(1 ) = rhs(1 ) - pppx1 * amx(1 ) !
rhs(nxc) = rhs(nxc) - pppx2 * apx(nxc) ! Added by H. Luo
END IF !

CALL tdma(amx,acx,apx,rhs,dummy,1,nxc)

DO i=1,nxc
var(i,j,k) = var(i,j,k) + omega*(dummy(i)-var(i,j,k))
ENDDO

IF (pbcx1 == PBC_DIRICHLET .and. pbcx2 == PBC_DIRICHLET) THEN !
var(0, j,k) = pppx1 !
var(nx,j,k) = pppx2 ! Added by H. Luo
END IF

ENDDO ! j 
ENDDO ! k

CALL enforce_p_periodic(var)

! Line solve in the y-direction
DO k=1,nzc
DO i=1,nxc
DO j=1,nyc
amx(i) = dxcinv(i) *dxinv(i)*(1 - ium(i,j,k) )
apx(i) = dxcinv(i+1)*dxinv(i)*(1 - iup(i,j,k) )
acx(i) = - ( amx(i) + apx(i) )

amy(j) = dycinv(j) *dyinv(j)*(1 - jum(i,j,k) )
apy(j) = dycinv(j+1)*dyinv(j)*(1 - jup(i,j,k) )
acy(j) = - ( amy(j) + apy(j) )

amz(k) = dzcinv(k) *dzinv(k)*(1 - kum(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL)
apz(k) = dzcinv(k+1)*dzinv(k)*(1 - kup(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL)
acz(k) = - ( amz(k) + apz(k) )

rhs(j) = r(i,j,k) - var(i-1,j,k)*amx(i) & 
- var(i+1,j,k)*apx(i) & 
- var(i,j,k-1)*amz(k) &
- var(i,j,k+1)*apz(k) &
+ nlw0(i,j,k)

!-- Modify rhs and coefficients for Ghost cells in GCM

amy(j) = amy(j)*REAL(1-iblank(i,j,k),KIND=CGREAL)
apy(j) = apy(j)*REAL(1-iblank(i,j,k),KIND=CGREAL)
acy(j) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank(i,j,k),KIND=CGREAL) &
+ REAL(iblank(i,j,k),KIND=CGREAL) 
rhs(j) = rhs(j)*REAL(1-iblank(i,j,k),KIND=CGREAL) &
+ REAL(ghostcellMark(i,j,k),KIND=CGREAL)*var(i,j,k)
ENDDO

IF (bcy1 == BC_TYPE_PERIODIC .AND. bcy2 == BC_TYPE_PERIODIC) THEN
rhs(1) = rhs(1) - var(i,nyc,k)*amy(1)
rhs(nyc) = rhs(nyc) - var(i,1,k)*apy(nyc)
END IF

IF (pbcy1 == PBC_DIRICHLET .and. pbcy2 == PBC_DIRICHLET) THEN !
rhs(1 ) = rhs(1 ) - pppy1 * amy(1 ) !
rhs(nyc) = rhs(nyc) - pppy2 * apy(nyc) ! Added by H. Luo
END IF !

CALL tdma(amy,acy,apy,rhs,dummy,1,nyc)

DO j=1,nyc
var(i,j,k) = var(i,j,k) + omega*(dummy(j)-var(i,j,k))
ENDDO

ENDDO
ENDDO

! CALL calc_residual(var, r, maxres, loc)
! write(*,*) 'y maximum residual is: ', maxres


! Line solver in the z-direction
IF (ndim == DIM_3D) THEN
CALL enforce_p_periodic(var)
DO j=1,nyc
DO i=1,nxc
DO k=1,nzc
amx(i) = dxcinv(i) *dxinv(i)*(1 - ium(i,j,k) )
apx(i) = dxcinv(i+1)*dxinv(i)*(1 - iup(i,j,k) )
acx(i) = - ( amx(i) + apx(i) )

amy(j) = dycinv(j) *dyinv(j)*(1 - jum(i,j,k) )
apy(j) = dycinv(j+1)*dyinv(j)*(1 - jup(i,j,k) )
acy(j) = - ( amy(j) + apy(j) )

amz(k) = dzcinv(k) *dzinv(k)*(1 - kum(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL)
apz(k) = dzcinv(k+1)*dzinv(k)*(1 - kup(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL)
acz(k) = - ( amz(k) + apz(k) )

rhs(k) = r(i,j,k) - var(i,j-1,k)*amy(j) &
- var(i,j+1,k)*apy(j) &
- var(i-1,j,k)*amx(i) &
- var(i+1,j,k)*apx(i) &
+ nlw0(i,j,k)

amz(k) = amz(k)*REAL(1-iblank(i,j,k),KIND=CGREAL)
apz(k) = apz(k)*REAL(1-iblank(i,j,k),KIND=CGREAL)
acz(k) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank(i,j,k),KIND=CGREAL) &
+ REAL(iblank(i,j,k),KIND=CGREAL) 
rhs(k) = rhs(k)*REAL(1-iblank(i,j,k),KIND=CGREAL) &
+ REAL(ghostcellMark(i,j,k),KIND=CGREAL)*var(i,j,k)
ENDDO ! k

IF (bcz1 == BC_TYPE_PERIODIC .AND. bcz2 == BC_TYPE_PERIODIC) THEN
rhs(1) = rhs(1) - var(i,j,nzc)*amz(1)
rhs(nzc) = rhs(nzc) - var(i,j,1)*apz(nzc)
END IF

IF (pbcz1 == PBC_DIRICHLET .and. pbcz2 == PBC_DIRICHLET) THEN !
rhs(1 ) = rhs(1 ) - pppz1 * amz(1 ) !
rhs(nzc) = rhs(nzc) - pppz2 * apz(nzc) ! Added by H. Luo
END IF 

CALL tdma(amz,acz,apz,rhs,dummy,1,nzc)

DO k=1,nzc
var(i,j,k) = var(i,j,k) + omega*(dummy(k)-var(i,j,k))
ENDDO ! k

ENDDO ! i 
ENDDO ! j
ENDIF ! ndim

! CALL calc_residual(var, r, maxres, loc)
! write(*,*) 'z maximum residual is: ', maxres

END SUBROUTINE itsolv
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE calc_residual(var,r,resm,loc)

USE global_parameters
USE flow_parameters
USE flow_arrays
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays
USE solver_arrays
USE GCM_arrays

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT (IN) ::var,r
REAL(KIND=CGREAL), INTENT (OUT) ::resm
INTEGER, DIMENSION(3), INTENT (OUT) ::loc

INTEGER :: i,j,k
INTEGER :: iG,jG,iBody,iRow,n
REAL(KIND=CGREAL) :: res
REAL(KIND=CGREAL) :: bmx,bpx,bcx,bc
REAL(KIND=CGREAL) :: bmy,bpy,bcy
REAL(KIND=CGREAL) :: bmz,bpz,bcz
REAL(KIND=CGREAL) :: AAe, AAw, AAn, AAs, AAf, AAb

!******************************************************************************

loc=0
resm = zero

DO k=1,nzc
DO j=1,nyc
DO i=1,nxc

IF (cure_pressure_oscillations .AND. iblank(i,j,k)==0 .AND. ivc(i,j,k)>0) THEN
AAe=face(cell(ivc(i,j,k))%F_ip)%a
AAw=face(cell(ivc(i,j,k))%F_im)%a
AAn=face(cell(ivc(i,j,k))%F_jp)%a
AAs=face(cell(ivc(i,j,k))%F_jm)%a
AAf=face(cell(ivc(i,j,k))%F_kp)%a
AAb=face(cell(ivc(i,j,k))%F_km)%a

bmx = dxcinv(i) *dxinv(i)*(1 - ium(i,j,k) ) * AAw
bpx = dxcinv(i+1)*dxinv(i)*(1 - iup(i,j,k) ) * AAe
bcx = - ( bmx + bpx )

bmy = dycinv(j) *dyinv(j)*(1 - jum(i,j,k) ) * AAs
bpy = dycinv(j+1)*dyinv(j)*(1 - jup(i,j,k) ) * AAn
bcy = - ( bmy + bpy )

bmz = dzcinv(k) *dzinv(k)*(1 - kum(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL) * AAb
bpz = dzcinv(k+1)*dzinv(k)*(1 - kup(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL) * AAf
bcz = - ( bmz + bpz )

bc = (bcx+bcy+bcz)

ELSE
bmx = dxcinv(i) *dxinv(i)*(1 - ium(i,j,k) )
bpx = dxcinv(i+1)*dxinv(i)*(1 - iup(i,j,k) )
bcx = - ( bmx + bpx )

bmy = dycinv(j) *dyinv(j)*(1 - jum(i,j,k) )
bpy = dycinv(j+1)*dyinv(j)*(1 - jup(i,j,k) )
bcy = - ( bmy + bpy )

bmz = dzcinv(k) *dzinv(k)*(1 - kum(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL)
bpz = dzcinv(k+1)*dzinv(k)*(1 - kup(i,j,k) ) &
*REAL((ndim - DIM_2D),KIND=CGREAL)
bcz = - ( bmz + bpz )

bmx = bmx*REAL(1-iblank(i,j,k),KIND=CGREAL)
bpx = bpx*REAL(1-iblank(i,j,k),KIND=CGREAL)

bmy = bmy*REAL(1-iblank(i,j,k),KIND=CGREAL)
bpy = bpy*REAL(1-iblank(i,j,k),KIND=CGREAL)

bmz = bmz*REAL(1-iblank(i,j,k),KIND=CGREAL)
bpz = bpz*REAL(1-iblank(i,j,k),KIND=CGREAL)

bc = (bcx+bcy+bcz)*REAL(1-iblank(i,j,k),KIND=CGREAL) &
+ REAL(iblank(i,j,k),KIND=CGREAL) 
ENDIF
res = r(i,j,k) - var(i,j,k)*bc &
- var(i-1,j,k)*bmx &
- var(i+1,j,k)*bpx &
- var(i,j-1,k)*bmy &
- var(i,j+1,k)*bpy &
- var(i,j,k-1)*bmz &
- var(i,j,k+1)*bpz & 
! - dxinv(i)*( iup(i,j,k)*pgradx2(j,k) - ium(i,j,k)*pgradx1(j,k) ) &
! *(1.0_CGREAL - pot_flag(i,j,k))&
! - dyinv(j)*( jup(i,j,k)*pgrady2(i,k) - jum(i,j,k)*pgrady1(i,k) ) &
! *(1.0_CGREAL - pot_flag(i,j,k))&
! - dzinv(k)*( kup(i,j,k)*pgradz2(i,j) - kum(i,j,k)*pgradz1(i,j) ) &
! *(1.0_CGREAL - pot_flag(i,j,k))
+ nlw0(i,j,k)
res = res*REAL(1-iblank(i,j,k),KIND=CGREAL)
IF (ABS(res) > resm ) THEN
resm = ABS(res)
loc(1) = i
loc(2) = j
loc(3) = k
ENDIF 

ENDDO
ENDDO
ENDDO

END SUBROUTINE calc_residual
!!------------------------------------------------------------------------------- 
#ifdef PETSC
SUBROUTINE petsc_solver(var,r)

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays
! use petsc_mod
USE GCM_arrays

IMPLICIT NONE

INTEGER :: i,j,k
INTEGER :: ind, dimPetsc, iterations
LOGICAL :: converged

REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,0:nz+1),INTENT (IN OUT) ::var
REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,0:nz+1),INTENT (IN) ::r
REAL(KIND=CGREAL),DIMENSION((nxc)*(nyc)*(nzc)) :: f, uu

!******************************************************************************

ind = 0

DO k=1,nzc
DO j=1,nyc
DO i=1,nxc
ind = ind + 1
f(ind) = r(i,j,k) &
- dxinv(i)*( iup(i,j,k)*pgradx2(j,k) - ium(i,j,k)*pgradx1(j,k) ) &
- dyinv(j)*( jup(i,j,k)*pgrady2(i,k) - jum(i,j,k)*pgrady1(i,k) ) &
- dzinv(k)*( kup(i,j,k)*pgradz2(i,j) - kum(i,j,k)*pgradz1(i,j) ) 

f(ind) = f(ind)*REAL(1-iblank(i,j,k),KIND=CGREAL) &
+ REAL(ghostcellMark(i,j,k),KIND=CGREAL)*var(i,j,k)

uu(ind) = var(i,j,k) 
ENDDO ! i
ENDDO ! j
ENDDO ! k

dimPetsc = (nxc)*(nyc)*(nzc)

! call petsc_run_solver(dimPetsc, uu, f, iterations, converged) !commented out by Wanh temporarily

! IF (MOD(ntime,nmonitor) == 0) THEN
! IF (converged) then
! PRINT*, 'Solution converged in ', iterations, ' iterations'
! ELSE
! PRINT*, 'Solution did not converge in ', iterations, ' iterations'
! ENDIF ! converged
! ENDIF ! ntime

ind=0
DO k=1,nzc
DO j=1,nyc
DO i=1,nxc
ind = ind + 1
var(i,j,k) = uu(ind) 
ENDDO ! i 
ENDDO ! j
ENDDO ! k

END SUBROUTINE petsc_solver
#endif

SUBROUTINE Pressure_Aitken(loc,iter)

USE global_parameters
USE flow_parameters
USE global_parameters
USE pressure_arrays
USE Pressure_Aitken_Array

IMPLICIT NONE

INTEGER :: iter
INTEGER, DIMENSION(3) ::loc
INTEGER, PARAMETER :: findAlphaSize = 40

INTEGER :: i,j,k,imin,imax,jmin,jmax,kmin,kmax

REAL(KIND=CGREAL) :: deltapPrime(0:nx+1,0:ny+1,0:nz+1)
REAL(KIND=CGREAL) :: InnerProduct_Ratio, ALPHA_PressureAitken
REAL(KIND=CGREAL) :: Ratio_nominator, Ratio_denominator
REAL(KIND=CGREAL) :: LambdaP_prev

Ratio_nominator = 0_CGREAL
Ratio_denominator = 0_CGREAL

IF (ndim == DIM_3D) THEN
! kmin = loc(3)-findAlphaSize
! kmax = loc(3)+findAlphaSize
kmin = 0
kmax = nz+1
! kmin = loc(3)-findAlphaSize
! kmax = loc(3)+findAlphaSize
ELSE
kmin = loc(3)
kmax = loc(3)
ENDIF

! imin = loc(1)-findAlphaSize
! imax = loc(1)+findAlphaSize
! jmin = loc(2)-findAlphaSize
! jmax = loc(2)+findAlphaSize
imin = 0
imax = nx+1
jmin = 0
jmax = ny+1

! Calculation of under-relaxtion coefficient 
DO k = kmin, kmax
DO j = jmin, jmax
DO i = imin, imax
deltapPrime_Lplus1(i,j,k) = pPrime_L(i,j,k)-pPrime(i,j,k) !pPrime is value of w/o iteration
diff_deltapPrimeL(i,j,k) = deltapPrime_prevL(i,j,k)-deltapPrime_Lplus1(i,j,k)
deltapPrime_prevL(i,j,k) = deltapPrime_Lplus1(i,j,k)
Ratio_nominator = Ratio_nominator + diff_deltapPrimeL(i,j,k)*deltapPrime_Lplus1(i,j,k)
Ratio_denominator = Ratio_denominator + diff_deltapPrimeL(i,j,k)*diff_deltapPrimeL(i,j,k)
ENDDO
ENDDO
ENDDO

IF (iter == 1) THEN
LambdaP = LambdaP_Init
ELSE
InnerProduct_Ratio = Ratio_nominator/Ratio_denominator
print *, 'Before update LambdaP, InnerProduct_Ratio =', LambdaP, InnerProduct_Ratio
print *, 'Ratio_nom, denom =', Ratio_nominator, Ratio_denominator
LambdaP = LambdaP+(LambdaP-1)*InnerProduct_Ratio
IF (abs(LambdaP-1.0)<1e-5) THEN
LambdaP = LambdaP_Init
ELSE IF (LambdaP<0) THEN
LambdaP = 0
ENDIF
LambdaP_prev = LambdaP
! print *, 'After update LambdaP = ', LambdaP
ENDIF

ALPHA_PressureAitken = 1 - LambdaP
! write (*,*) 'Alpha for pressure Aitken acceleration is: ', ALPHA_PressureAitken
! write (*,*) ' '

DO k = kmin, kmax
DO j = jmin, jmax
DO i = imin, imax
pPrime(i,j,k) = pPrime_L(i,j,k)-ALPHA_PressureAitken*deltapPrime_Lplus1(i,j,k)
pPrime_L(i,j,k) = pPrime(i,j,k)
ENDDO
ENDDO
ENDDO

! pPrime = pPrime_L - ALPHA_PressureAitken*deltapPrime_Lplus1
! pPrime_L = pPrime

END SUBROUTINE Pressure_Aitken


!-------------------------------------------------------------------------------

SUBROUTINE Setup_Poisson_Tester_RHS(Problem)

USE global_parameters
USE flow_parameters
USE flow_arrays
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays

IMPLICIT NONE

REAL(KIND=CGREAL) :: sum, x1

INTEGER :: i,j,k,problem

!******************************************************************************

! rhs = [ d(U)/dx + d(V)/dy + d(W)/dz ] / dt

sum = zero
nlu0 = zero

DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
! nlu0(i,j,k) = (dtinv)* & 
! ( ( sin(2*pi*x(i+1)/xout)*sin(2*pi*yc(j)/yout) - sin(2*pi*x(i)/xout)*sin(2*pi*yc(j)/yout) )*dxinv(i) &
! +( sin(2*pi*xc(i)/xout)*sin(2*pi*y(j+1)/yout) - sin(2*pi*xc(i)/xout)*sin(2*pi*y(j)/yout) )*dyinv(j) )
! nlu0(i,j,k) = -3.d0 * cos(xc(i))*cos(yc(j))*cos(zc(k))
nlu0(i,j,k) = -2.d0 * cos(xc(i))*cos(yc(j))
! nlu0(i,j,k) = nlu0(i,j,k)*REAL(1-iblank(i,j,k),KIND=CGREAL)
sum = sum + nlu0(i,j,k)*dx(i)*dy(j)*dz(k)
ENDDO
ENDDO
ENDDO

! CALL write_dump_debug('nlu0',0,nlu0)

if (pp_solver_type >3) then 
DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
IF (J==1) THEN
END IF
IF (J==NYc) THEN
END IF
IF (I==1) THEN
END IF
IF (I==NXc) THEN
END IF
ENDDO
ENDDO
ENDDO
end if

WRITE(STDOUT,'(5X,A,1X,2(1PE21.12))') 'In rhs poisson, max/min nlu0:', maxval(nlu0), minval(nlu0)
WRITE(STDOUT,'(5X,A,1X,1PE21.12)') 'Sum Of Poisson RHS =',sum

END SUBROUTINE Setup_Poisson_Tester_RHS

!-------------------------------------------------------------------------------

SUBROUTINE Poisson_Tester()

USE global_parameters
USE flow_parameters
USE flow_arrays
USE pressure_arrays
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays
USE mg_parameters

IMPLICIT NONE
! 
EXTERNAL MG_Precondition_MSIP_2D, MG_Precondition_MSIP_3D
EXTERNAL MG_itsolv_MSIP_2D, MG_itsolv_MSIP_3D
EXTERNAL MG_Residual_MSIP

EXTERNAL MG_itsolv
EXTERNAL MG_itsolv_Point_Jacobi

EXTERNAL MG_Residual

INTEGER :: iter, M
INTEGER :: loc(3), loc2(3)

INTEGER :: clock1, clock2, clock_rate
REAL(KIND=CGREAL) :: maxres, maxres1, maxres2, maxres3
REAL(KIND=CGREAL) :: compTime
! 

CALL Setup_Poisson_Tester_RHS(3)

IF (pbcx1 == PBC_DIRICHLET .OR. &
pbcy1 == PBC_DIRICHLET .OR. &
pbcz1 == PBC_DIRICHLET .OR. &
pbcx2 == PBC_DIRICHLET .OR. &
pbcy2 == PBC_DIRICHLET .OR. &
pbcz2 == PBC_DIRICHLET) THEN

CALL set_pressure_dirichlet_bc()

ELSE

CALL subtracting_average_pressure()

END IF

call system_clock(clock1)

SELECT CASE(pp_solver_type)

CASE (PP_SOLVER_TYPE_SIP)
IF (ndim == DIM_2D) THEN
CALL COEFF_PRESSURE_MSIP_2D()
CALL PRECONDITION_SIP_2D()
ELSE
WRITE(STDOUT,*) 'SIP 3D is not available.'
WRITE(STDOUT,*) 'Please select other algorithms.'
STOP
END IF
CASE (PP_SOLVER_TYPE_MSIP) 
IF (ndim == DIM_2D) THEN
CALL COEFF_PRESSURE_MSIP_2D()
CALL COEFF_LU_MSIP_2D()
ELSE
CALL COEFF_PRESSURE_MSIP_3D()
CALL COEFF_LU_MSIP_3D()
END IF
CASE (PP_SOLVER_TYPE_MG_MSIP) 
CALL MG_Memory_Allocation_iblank(mgLevels_X,ICOORD)

IF (Full_Coarsening) THEN
CALL MG_Prepare_Iblank_FMG
ELSE
CALL MG_Memory_Allocation_iblank(mgLevels_Y,JCOORD)
CALL MG_Prepare_Iblank(ICOORD)
CALL MG_Prepare_Iblank(JCOORD)
IF (ndim == DIM_3D) THEN
CALL MG_Memory_Allocation_iblank(mgLevels_Z,KCOORD)
CALL MG_Prepare_Iblank(KCOORD)
END IF

ENDIF
IF (ndim == DIM_2D) THEN
CALL MG_PREPARE_MSIP(MG_Precondition_MSIP_2D)
ELSE
CALL MG_PREPARE_MSIP(MG_Precondition_MSIP_3D)
ENDIF
END SELECT

do M=1, 1

iter = 0
maxres = 1.0E10_CGREAL 
pPrime=zero
DO WHILE ((iter .LT. iterMax_Poisson) .AND. (maxres .GT. restol_Poisson))

SELECT CASE(pp_solver_type)

CASE (PP_SOLVER_TYPE_LSOR)
CALL itsolv(pPrime,nlu0)
CALL calc_residual(pPrime,nlu0,maxres,loc)
CASE (PP_SOLVER_TYPE_MG)
CALL mg_solver(pPrime,nlu0,MG_itsolv,MG_Residual,compTime)
CALL calc_residual(pPrime,nlu0,maxres,loc)
CASE (PP_SOLVER_TYPE_MG_Point_Jacobi) 
CALL mg_solver(pPrime,nlu0,MG_itsolv_Point_Jacobi,MG_Residual,compTime)
CALL calc_residual(pPrime,nlu0,maxres,loc)
CASE (PP_SOLVER_TYPE_SIP) 
IF (ndim == DIM_3D) THEN
CALL itsolv_SIP_3D(pPrime,nlu0)
ELSE
CALL itsolv_SIP_2D(pPrime,nlu0)
END IF
CALL calc_residual(pPrime,nlu0,maxres,loc)
CASE (PP_SOLVER_TYPE_MSIP) 
IF (ndim == DIM_3D) THEN
CALL itsolv_MSIP_3D(pPrime,nlu0,restol_Poisson)
ELSE
CALL itsolv_MSIP_2D(pPrime,nlu0,restol_Poisson)
END IF
CALL calc_residual(pPrime,nlu0,maxres,loc)
CASE (PP_SOLVER_TYPE_MG_MSIP) 
IF (Full_Coarsening) THEN

IF (ndim == DIM_2D) THEN
CALL mg_solver_FMG(pPrime,nlu0,MG_itsolv_MSIP_2D,MG_Residual_MSIP,compTime)
ELSE
CALL mg_solver_FMG(pPrime,nlu0,MG_itsolv_MSIP_3D,MG_Residual_MSIP,compTime)
ENDIF

ELSE
IF (ndim == DIM_2D) THEN
CALL mg_solver(pPrime,nlu0,MG_itsolv_MSIP_2D,MG_Residual_MSIP,compTime)
ELSE
CALL mg_solver(pPrime,nlu0,MG_itsolv_MSIP_3D,MG_Residual_MSIP,compTime)
ENDIF

ENDIF 
! CALL calc_residual_MG_MSIP(pPrime,nlu,maxres,loc)
CALL calc_residual(pPrime,nlu0,maxres,loc)

END SELECT
iter = iter + 1
write(1001,*) iter, maxres
write(1002,*) iter, compTime
call system_clock(clock2, clock_rate)
WRITE(1003,*) maxres, REAL(clock2-clock1)/REAL(clock_rate)
IF (MOD(ntime,nmonitor) == 0) &
WRITE(STDOUT,'(5X,I5.5,1X,E19.11,1X,A,I4.4,A,I4.4,A,I4.4,A)') &
iter,maxres,'(',loc(1),',',loc(2),',',loc(3),')'
CALL calc_grid_convergence_test(pPrime,nlu0,maxres3,loc2)
CALL calc_grid_convergence_test1(pPrime,nlu0,maxres1)
CALL calc_grid_convergence_test2(pPrime,nlu0,maxres2)
WRITE(1004,'(4(2X,E19.11))') dx(1), maxres1, maxres2, maxres3
ENDDO
enddo
CALL write_dump_debug('p   ',0,pPrime)

call system_clock(clock2, clock_rate)
WRITE(*,*) '*****Total time for solving poisson eq is:', &
REAL(clock2-clock1)/REAL(clock_rate)

stop

END SUBROUTINE Poisson_Tester

!-------------------------------------------------------------------------------
SUBROUTINE calc_grid_convergence_test(var,r,resm,loc)

USE global_parameters
USE flow_parameters
USE boundary_arrays
USE grid_arrays
USE MG_parameters
USE MG_arrays

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT (IN) ::var,r
REAL(KIND=CGREAL), INTENT (OUT) ::resm
INTEGER, DIMENSION(3), INTENT (OUT) ::loc

REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

INTEGER :: i,j,k, MZ
INTEGER :: iG,jG
REAL(KIND=CGREAL) :: res
REAL(KIND=CGREAL) :: bmx,bpx,bcx,bc
REAL(KIND=CGREAL) :: bmy,bpy,bcy
REAL(KIND=CGREAL) :: bmz,bpz,bcz

!******************************************************************************
loc = 0

AS=>MGX(1)%CA(1,:,:,:)
AW=>MGX(1)%CA(2,:,:,:)
AP=>MGX(1)%CA(3,:,:,:)
AE=>MGX(1)%CA(4,:,:,:)
AN=>MGX(1)%CA(5,:,:,:)
IF (ndim == DIM_3D) THEN
AF=>MGX(1)%CA(6,:,:,:)
AB=>MGX(1)%CA(7,:,:,:)
END IF

resm = zero

IF (ndim == DIM_3D) THEN
MZ=nzc
ELSE
MZ=1
END IF

DO k=1,mz 
DO j=1,nyc
DO i=1,nxc

res = VAR(i,j,k) - cos(xc(i))*cos(yc(j))

IF (ABS(res) > resm ) THEN
resm = ABS(res)
loc(1) = i
loc(2) = j
loc(3) = k
ENDIF 

ENDDO
ENDDO
ENDDO
END SUBROUTINE calc_grid_convergence_test
!!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE calc_grid_convergence_test1(var,r,resm)

USE global_parameters
USE flow_parameters
USE boundary_arrays
USE grid_arrays
USE MG_parameters
USE MG_arrays

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT (IN) ::var,r
REAL(KIND=CGREAL), INTENT (OUT) ::resm

REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

INTEGER :: i,j,k, MZ
INTEGER :: iG,jG
REAL(KIND=CGREAL) :: res
REAL(KIND=CGREAL) :: bmx,bpx,bcx,bc
REAL(KIND=CGREAL) :: bmy,bpy,bcy
REAL(KIND=CGREAL) :: bmz,bpz,bcz

!******************************************************************************
AS=>MGX(1)%CA(1,:,:,:)
AW=>MGX(1)%CA(2,:,:,:)
AP=>MGX(1)%CA(3,:,:,:)
AE=>MGX(1)%CA(4,:,:,:)
AN=>MGX(1)%CA(5,:,:,:)
IF (ndim == DIM_3D) THEN
AF=>MGX(1)%CA(6,:,:,:)
AB=>MGX(1)%CA(7,:,:,:)
END IF

res = zero

IF (ndim == DIM_3D) THEN
MZ=nzc
ELSE
MZ=1
END IF

DO k=1,mz 
DO j=1,nyc
DO i=1,nxc

res = res + abs((VAR(i,j,k) - cos(xc(i))*cos(yc(j))))

ENDDO
ENDDO
ENDDO

resm=res/nxc/nyc/mz

END SUBROUTINE calc_grid_convergence_test1
!!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE calc_grid_convergence_test2(var,r,resm)

USE global_parameters
USE flow_parameters
USE boundary_arrays
USE grid_arrays
USE MG_parameters
USE MG_arrays

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1), INTENT (IN) ::var,r
REAL(KIND=CGREAL), INTENT (OUT) ::resm

REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

INTEGER :: i,j,k, MZ
INTEGER :: iG,jG
REAL(KIND=CGREAL) :: res
REAL(KIND=CGREAL) :: bmx,bpx,bcx,bc
REAL(KIND=CGREAL) :: bmy,bpy,bcy
REAL(KIND=CGREAL) :: bmz,bpz,bcz

!******************************************************************************
AS=>MGX(1)%CA(1,:,:,:)
AW=>MGX(1)%CA(2,:,:,:)
AP=>MGX(1)%CA(3,:,:,:)
AE=>MGX(1)%CA(4,:,:,:)
AN=>MGX(1)%CA(5,:,:,:)
IF (ndim == DIM_3D) THEN
AF=>MGX(1)%CA(6,:,:,:)
AB=>MGX(1)%CA(7,:,:,:)
END IF

res = zero

IF (ndim == DIM_3D) THEN
MZ=nzc
ELSE
MZ=1
END IF

DO k=1,mz 
DO j=1,nyc
DO i=1,nxc

res = res + (VAR(i,j,k) - cos(xc(i))*cos(yc(j)))**2

ENDDO
ENDDO
ENDDO

resm=sqrt(res)/nxc/nyc/mz
END SUBROUTINE calc_grid_convergence_test2
!!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE GCM_correct_pressure() !r,mx,my,mz)
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays
USE pressure_arrays

IMPLICIT NONE

! INTEGER :: mx,my,mz
! REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (INOUT) :: r

INTEGER :: I,J,K
!
CALL set_outer_ghost_pres(pPrime, nx+1, ny+1, nz+1)
CALL GCM_p_set_bc_internal(pPrime, nx+1, ny+1, nz+1)
CALL GCM_enforce_p_compatibility(pPrime, nx+1, ny+1, nz+1)

DO k = 1, nzc
DO j = 1, nyc
DO i = 1, nxc
nlw0(i,j,k) = - dxinv(i)*( iup(i,j,k)*pgradx2(j,k) - ium(i,j,k)*pgradx1(j,k) )*(oned - pot_flag(i,j,k)) &
- dyinv(j)*( jup(i,j,k)*pgrady2(i,k) - jum(i,j,k)*pgrady1(i,k) )*(oned - pot_flag(i,j,k)) &
- dzinv(k)*( kup(i,j,k)*pgradz2(i,j) - kum(i,j,k)*pgradz1(i,j) )*(oned - pot_flag(i,j,k))
END DO
END DO
END DO

END SUBROUTINE GCM_correct_pressure
!-------------------------------------------------------------------------------

