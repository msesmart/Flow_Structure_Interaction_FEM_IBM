!---------------------------------
! SUBROUTINE rhs_advec_diff() 
! SUBROUTINE rhs_adjust_bc()
!---------------------------------


!-------------------------------------------------------------------------------
! n { n n-1} { n n n }
! RHS = u - dt { 3/2 NL -1/2 NL } + dt(1/2Re) {am u + ac u + ap u }
! { i i } { i i-1 i i i i+1} 
!-------------------------------------------------------------------------------

SUBROUTINE rhs_advec_diff() 

USE global_parameters
USE flow_parameters
USE flow_arrays
USE grid_arrays
USE boundary_arrays
USE multiuse_arrays
USE nlold_arrays
USE pressure_arrays
USE solver_ad_arrays

IMPLICIT NONE

!------------------------------------------------------------------------------
! Non-linear convection term
!------------------------------------------------------------------------------

CALL rhs_advec(u,nlu,bcxu,bcyu,bczu) 
CALL rhs_advec(v,nlv,bcxv,bcyv,bczv) 
IF ( nDim == DIM_3D ) THEN
CALL rhs_advec(w,nlw,bcxw,bcyw,bczw)
END IF ! nDim

!------------------------------------------------------------------------------
! Load Adams-Bashforth terms
!------------------------------------------------------------------------------

CALL rhs_loadAB(nlu,nluOld)
CALL rhs_loadAB(nlv,nlvOld)
IF ( nDim == DIM_3D ) THEN
CALL rhs_loadAB(nlw,nlwOld)
END IF ! nDim

!------------------------------------------------------------------------------ 
! Diffusion terms 
!------------------------------------------------------------------------------

CALL rhs_diff(u,nlu,bcxu,bcyu,bczu) 
CALL rhs_diff(v,nlv,bcxv,bcyv,bczv) 
IF ( nDim == DIM_3D ) THEN
CALL rhs_diff(w,nlw,bcxw,bcyw,bczw)
END IF ! nDim

!------------------------------------------------------------------------------
! Invoke Van-kan formalism
!------------------------------------------------------------------------------

IF (frac_step_type == VAN_KAN) CALL rhs_vankan

IF (boundary_motion_type(1) == FEA_FLOW_STRUC_INTERACTION .OR. &
boundary_motion_type(1) == PARTIAL_DYNAMICS_COUPLED .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED .OR. &
boundary_motion_type(1) == BIO_DYNAMICS_COUPLED .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING ) THEN
! Save nlu, nlv, nlw for interation of FSI, the following 3 lines are added by Wanh 06/24/11
nlu_FSI = nlu
nlv_FSI = nlv
nlw_FSI = nlw
ENDIF

!------------------------------------------------------------------------------
CONTAINS

SUBROUTINE rhs_advec(vel,nlvel,bcxvel,bcyvel,bczvel) 

USE global_parameters
USE flow_parameters
USE flow_arrays
USE grid_arrays
USE boundary_arrays

IMPLICIT NONE

!... parameters

REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN) :: vel
REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN) :: bcxvel
REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN) :: bcyvel
REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN) :: bczvel
REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(OUT) :: nlvel

!... loop variables

INTEGER :: i,j,k
INTEGER :: imm,jmm,kmm ! Added by Rupesh
INTEGER :: ipp,jpp,kpp ! used in 2nd Upwinding


!... local variables
REAL(KIND=CGREAL) :: edxWeightm1, edxWeightm2, edxWeightp1, edxWeightp2
REAL(KIND=CGREAL) :: wdxWeightm1, wdxWeightm2, wdxWeightp1, wdxWeightp2

REAL(KIND=CGREAL) :: ndxWeightm1, ndxWeightm2, ndxWeightp1, ndxWeightp2
REAL(KIND=CGREAL) :: sdxWeightm1, sdxWeightm2, sdxWeightp1, sdxWeightp2

REAL(KIND=CGREAL) :: bdxWeightm1, bdxWeightm2, bdxWeightp1, bdxWeightp2
REAL(KIND=CGREAL) :: fdxWeightm1, fdxWeightm2, fdxWeightp1, fdxWeightp2

REAL(KIND=CGREAL) :: Usign, UsignP, Vsign, VsignP, Wsign, WsignP

REAL(KIND=CGREAL) :: vele,velw,veln,vels,velf,velb
REAL(KIND=CGREAL) :: vele_Up,velw_Up,veln_Up,vels_Up,velf_Up,velb_Up ! Added by Rupesh 
REAL(KIND=CGREAL) :: tempvel

!******************************************************************************

!---------------------------------------------------------------------------------------------------
! Convective terms:
! nlvel = d(vel*U)/dx + d(vel*V)/dy + d(vel*W)/dz
!---------------------------------------------------------------------------------------------------
! The following decription was added by Rupesh and pertains to 2nd Upwind scheme for 
! convected face velocities. 
! ___________________________________________________
! | | | | | |
! | | | | | |
! | o | o w+ o +e o | o |
! | WW | W | P | E | EE |
! |__________|_________|_________|_________|__________|
!
! |<--S_uu-->|<--S_u-->|<--S_c-->|<--S_d-->|<--S_dd-->| 1.Flow: Left --> Right
! 
! |<--S_dd-->|<--S_d-->|<--S_c-->|<--S_u-->|<--S_uu-->| 2.Flow: Left <-- Right 
!
! Subscripts 'u' and 'd' in S refer to upstream and downstream resply.
!
! LOGIC: 1. If the flow is from left to right across a face (say, face e), then
!
! | { S_u + 2S_c} { S_c }
! u_e| = {-----------}u_P - {-----------}u_W
! |up { S_u + S_c } { S_u + S_c } 
!
! 2. If the flow is from right to left across a face (say, face e), then
!
! | { S_uu + 2S_u} { S_u }
! u_e| = {------------}u_E - {------------}u_EE
! |up { S_uu + S_u } { S_uu + S_u } 
! 
! - It should be noted that for u_w, the above formulae are still valid, provided the stencil 
! is offset by one cell.
! - These formulae are derived from: u_face = u + (Grad u)(dot)(dS), where 
!
! 'u' and 'Grad u' are cellcentered value and its gradient in the upstream cell resply, 
! and dS is the displacement vector from the upstream cell-centroid to the face centroid.
! 'Grad u' is approximated by upwind differencing based on the direction of the wind.
!
!----------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! x-direction 
!------------------------------------------------------------------------------

DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
vele = ( fx(i+1) *vel(i+1,j,k) &
+ ( oned-fx(i+1) ) *vel(i,j,k) )*(1-iup(i,j,k)) &
+ bcxvel(i,j,k)*iup(i,j,k)

velw = ( fx(i) *vel(i,j,k) &
+ ( oned-fx(i) ) *vel(i-1,j,k) )*(1-ium(i,j,k)) &
+ bcxvel(i,j,k)*ium(i,j,k)

IF (Hybrid) THEN

imm = MAX(i-2,0) 
ipp = MIN(i+2,nx)

USign = SIGN(oned,face_u(i,j,k))
USignP = SIGN(oned,face_u(i+1,j,k))

!
! 2nd Upwind differencing added by Rupesh
!
edxWeightm1 = ( dxinv(i-1)+twod*dxinv(i) ) / &
( dxinv(i-1)+dxinv(i) )
edxWeightm2 = dxinv(i) / ( dxinv(i-1)+dxinv(i) )
edxWeightp1 = ( dxinv(ipp)+twod*dxinv(i+1) ) / &
( dxinv(ipp)+dxinv(i+1) )
edxWeightp2 = dxinv(i+1) / ( dxinv(ipp)+dxinv(i+1) )
!
wdxWeightm1 = ( dxinv(imm)+twod*dxinv(i-1) ) / &
( dxinv(imm)+dxinv(i-1) )
wdxWeightm2 = dxinv(i-1) / ( dxinv(imm)+dxinv(i-1) )
wdxWeightp1 = ( dxinv(i+1)+twod*dxinv(i) ) / &
( dxinv(i+1)+dxinv(i) )
wdxWeightp2 = dxinv(i) / ( dxinv(i+1)+dxinv(i) )
!
vele_Up = ( ( oned+UsignP ) &
*( iumm(i,j,k)*vele + (1-iumm(i,j,k) )* &
( edxWeightm1*vel(i,j,k) - edxWeightm2*vel(i-1,j,k) ) ) &
+ ( oned-UsignP ) &
*( iupp(i,j,k)*vele + (1-iupp(i,j,k) )* &
( edxWeightp1*vel(i+1,j,k) - edxWeightp1*vel(ipp,j,k) ) ) &
)*half*(1-iup(i,j,k)) + bcxvel(i,j,k)*iup(i,j,k)
!
velw_Up = ( ( oned+Usign ) &
*( iumm(i,j,k)*velw + (1-iumm(i,j,k) )* &
( wdxWeightm1*vel(i-1,j,k) - wdxWeightm2*vel(imm,j,k) ) ) &
+ ( oned-Usign ) &
*( iupp(i,j,k)*velw + (1-iupp(i,j,k) )* &
( wdxWeightp1*vel(i,j,k) - wdxWeightp2*vel(i+1,j,k) ) ) &
)*half*(1-ium(i,j,k)) + bcxvel(i,j,k)*ium(i,j,k)
!
! Original CDS definition commented by Rupesh
! 
nlvel(i,j,k) = ( ( (oned - alfa)*vele + alfa*vele_Up )* &
face_u(i+1,j,k) & 
- ( (oned - alfa)*velw + alfa*velw_Up )* &
face_u(i,j,k) &
)*dxinv(i) ! Hybrid definition added by Rupesh
ELSE
nlvel(i,j,k) = ( vele*face_u(i+1,j,k) -velw*face_u(i,j,k) )*dxinv(i) 
ENDIF ! HYBRID

ENDDO ! i
ENDDO ! j
ENDDO ! k

!------------------------------------------------------------------------------
! y-direction
!------------------------------------------------------------------------------

DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
veln = ( fy(j+1) *vel(i,j+1,k) &
+ ( oned-fy(j+1) ) *vel(i,j,k) )*(1-jup(i,j,k)) &
+ bcyvel(i,j,k)*jup(i,j,k)

vels = ( fy(j) *vel(i,j,k) &
+ ( oned-fy(j) ) *vel(i,j-1,k) )*(1-jum(i,j,k)) &
+ bcyvel(i,j,k)*jum(i,j,k)

IF (Hybrid) THEN

VSign = SIGN(oned,face_v(i,j,k))
VSignP = SIGN(oned,face_v(i,j+1,k))

jmm = MAX(j-2,0)
jpp = MIN(j+2,ny)

!
! 2nd Upwind differencing added by Rupesh
!
ndxWeightm1 = ( dyinv(j-1)+twod*dyinv(j) ) / &
( dyinv(j-1)+dyinv(j) )
ndxWeightm2 = dyinv(j) / ( dyinv(j-1)+dyinv(j) )
ndxWeightp1 = ( dyinv(jpp)+twod*dyinv(j+1) ) / &
( dyinv(jpp)+dyinv(j+1) )
ndxWeightp2 = dyinv(j+1) / ( dyinv(jpp)+dyinv(j+1) )
!
sdxWeightm1 = ( dyinv(jmm)+twod*dyinv(j-1) ) / &
( dyinv(jmm)+dyinv(j-1) )
sdxWeightm2 = dyinv(j-1) / ( dyinv(jmm)+dyinv(j-1) ) 
sdxWeightp1 = ( dyinv(j+1)+twod*dyinv(j) ) / &
( dyinv(j+1)+dyinv(j) )
sdxWeightp2 = dyinv(j) / ( dyinv(j+1)+dyinv(j) )
!
veln_Up = ( ( oned+VsignP ) &
*( jumm(i,j,k)*veln + (1-jumm(i,j,k) )* &
( ndxWeightm1*vel(i,j,k) - ndxWeightm2*vel(i,j-1,k) ) ) &
+ ( oned-VsignP ) &
*( jupp(i,j,k)*veln + (1-jupp(i,j,k) )* &
( ndxWeightp1*vel(i,j+1,k) - ndxWeightp2*vel(i,jpp,k) ) ) &
)*half*(1-jup(i,j,k)) + bcyvel(i,j,k)*jup(i,j,k)
!
vels_Up = ( ( oned+Vsign ) &
*( jumm(i,j,k)*vels + (1-jumm(i,j,k) )* &
( sdxWeightm1*vel(i,j-1,k) - sdxWeightm2*vel(i,jmm,k) ) ) &
+ ( oned-Vsign ) &
*( jupp(i,j,k)*vels + (1-jupp(i,j,k) )* &
( sdxWeightp1*vel(i,j,k) - sdxWeightp2*vel(i,j+1,k) ) ) &
)*half*(1-jum(i,j,k)) + bcyvel(i,j,k)*jum(i,j,k)
!
! Original CDS definition commented by Rupesh
!
nlvel(i,j,k) = nlvel(i,j,k) + &
( ( (oned - alfa)*veln + alfa*veln_Up )* &
face_v(i,j+1,k) &
- ( (oned - alfa)*vels + alfa*vels_Up )* &
face_v(i,j,k) &
)*dyinv(j) ! Hybrid definition added by Rupesh
ELSE
nlvel(i,j,k) = nlvel(i,j,k) & 
+ ( veln*face_v(i,j+1,k) -vels*face_v(i,j,k) )*dyinv(j)
ENDIF ! HYBRID

ENDDO ! i
ENDDO ! j
ENDDO ! k

!------------------------------------------------------------------------------
! z-direction
!------------------------------------------------------------------------------
IF ( nDim == DIM_3D ) THEN 
DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
velf = ( fz(k+1) *vel(i,j,k+1) &
+ ( oned-fz(k+1) ) *vel(i,j,k) )*(1-kup(i,j,k)) &
+ bczvel(i,j,k)*kup(i,j,k)

velb = ( fz(k) *vel(i,j,k) &
+ ( oned-fz(k) ) *vel(i,j,k-1) )*(1-kum(i,j,k)) &
+ bczvel(i,j,k)*kum(i,j,k)

IF (Hybrid) THEN

WSign = SIGN(oned,face_w(i,j,k))
WSignP = SIGN(oned,face_w(i,j,k+1))

kmm = MAX(k-2,0)
kpp = MIN(k+2,nz)

!
! 2nd Upwind differencing added by Rupesh
!
fdxWeightm1 = ( dzinv(k-1)+twod*dzinv(k) ) / &
( dzinv(k-1)+dzinv(k) )
fdxWeightm2 = dzinv(k) / ( dzinv(k-1)+dzinv(k) )
fdxWeightp1 = ( dzinv(kpp)+twod*dzinv(k+1) ) / &
( dzinv(kpp)+dzinv(k+1) )
fdxWeightp2 = dzinv(k+1) / ( dzinv(kpp)+dzinv(k+1) )
!
bdxWeightm1 = ( dzinv(kmm)+twod*dzinv(k-1) ) / &
( dzinv(kmm)+dzinv(k-1) )
bdxWeightm2 = dzinv(k-1) / ( dzinv(kmm)+dzinv(k-1) )
bdxWeightp1 = ( dzinv(k+1)+twod*dzinv(k) ) / &
( dzinv(k+1)+dzinv(k) )
bdxWeightp2 = dzinv(k) / ( dzinv(k+1)+dzinv(k) )
!
velf_Up = ( ( oned+WsignP ) &
*( kumm(i,j,k)*velf + (1-kumm(i,j,k) )* &
( fdxWeightm1*vel(i,j,k) - fdxWeightm2*vel(i,j,k-1) ) ) &
+ ( oned-WsignP ) &
*( kupp(i,j,k)*velf + (1-kupp(i,j,k) )* &
( fdxWeightp1*vel(i,j,k+1) - fdxWeightp2*vel(i,j,kpp) ) ) &
)*half*(1-kup(i,j,k)) + bczvel(i,j,k)*kup(i,j,k)
!
velb_Up = ( ( oned+Wsign ) &
*( kumm(i,j,k)*velb + (1-kumm(i,j,k) )* &
( bdxWeightm1*vel(i,j,k-1) - bdxWeightm2*vel(i,j,kmm) ) ) & 
+ ( oned-Wsign ) &
*( kupp(i,j,k)*velb + (1-kupp(i,j,k) )* &
( bdxWeightp1*vel(i,j,k) - bdxWeightp2*vel(i,j,k+1) ) ) &
)*half*(1-kum(i,j,k)) + bczvel(i,j,k)*kum(i,j,k)
!
! Original CDS definition commented by Rupesh
!
nlvel(i,j,k) = nlvel(i,j,k) + &
( ( (oned - alfa)*velf + alfa*velf_Up )* &
face_w(i,j,k+1) &
- ( (oned - alfa)*velb + alfa*velb_Up )* &
face_w(i,j,k) &
)*dzinv(k) ! Hybrid definition added by Rupesh 
ELSE
nlvel(i,j,k) = nlvel(i,j,k) & 
+ ( velf*face_w(i,j,k+1) -velb*face_w(i,j,k) )*dzinv(k)
ENDIF ! HYBRID

ENDDO ! i
ENDDO ! j
ENDDO ! k
END IF ! nDim

END SUBROUTINE rhs_advec
!------------------------------------------------------------------------------

SUBROUTINE rhs_loadAB(nlvel,nlvelOld) 

USE global_parameters
USE flow_parameters
USE boundary_arrays

IMPLICIT NONE

!... parameters

REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(INOUT) :: nlvel,nlvelOld

!... loop variables

INTEGER :: i,j,k

!... local variables

REAL(KIND=CGREAL) :: tempvel 

!******************************************************************************

! The CN scheme for the advection terms were added by H. Luo
IF (advec_scheme == ADAMS_BASHFORTH2) THEN 
DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
tempvel = exp_weight(i,j,k) *nlvel(i,j,k) &
- ( exp_weight(i,j,k)-oned ) *nlvelOld(i,j,k)
nlvelOld(i,j,k) = nlvel(i,j,k)
nlvel(i,j,k) = tempvel
ENDDO ! i
ENDDO ! j
ENDDO ! k 
ELSEIF (advec_scheme == CRANK_NICOLSON1 .or. &
advec_scheme == CRANK_NICOLSON2) THEN
DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
nlvel(i,j,k) = nlvel(i,j,k) * half
ENDDO ! i
ENDDO ! j
ENDDO ! k 
ENDIF

END SUBROUTINE rhs_loadAB
!------------------------------------------------------------------------------

SUBROUTINE rhs_diff(vel,nlvel,bcxvel,bcyvel,bczvel) 

USE global_parameters
USE flow_parameters
USE flow_arrays
USE grid_arrays
USE boundary_arrays
USE solver_ad_arrays

IMPLICIT NONE

!... parameters

REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN) :: vel
REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN) :: bcxvel
REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN) :: bcyvel
REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN) :: bczvel
REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(INOUT) :: nlvel

!... loop variables

INTEGER :: i,j,k

!... local variables

REAL(KIND=CGREAL) :: acx,amx,apx,acy,amy,apy,acz,amz,apz
REAL(KIND=CGREAL) :: diffvel, rnDim
REAL(KIND=CGREAL) :: nuE,nuW,nuS,nuN,nuF,nuB

!******************************************************************************

!------------------------------------------------------------------------------
! Diffusive terms
! nlvel = d/dx_j [(1/Re+nut) d(vel)/dx_j ]
!
! Note: Since amx_ad coefficients have already a negative sign in them (-1/2 dt)
! diffvel needs only to be subtracted in computing nlvel term
!------------------------------------------------------------------------------

rnDim = REAL((ndim - DIM_2D),KIND=CGREAL) 

DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
nuE = ( fx(i+1) *viscTot(i+1,j,k) &
+ ( oned-fx(i+1) ) *viscTot(i,j,k) )*(1-iup(i,j,k)) &
+ bcxvisc(i,j,k)*iup(i,j,k)

nuW = ( fx(i) *viscTot(i,j,k) &
+ ( oned-fx(i) ) *viscTot(i-1,j,k) )*(1-ium(i,j,k)) &
+ bcxvisc(i,j,k)*ium(i,j,k)

nuN = ( fy(j+1) *viscTot(i,j+1,k) &
+ ( oned-fy(j+1) ) *viscTot(i,j,k) )*(1-jup(i,j,k)) &
+ bcyvisc(i,j,k)*jup(i,j,k)

nuS = ( fy(j) *viscTot(i,j,k) &
+ ( oned-fy(j) ) *viscTot(i,j-1,k) )*(1-jum(i,j,k)) &
+ bcyvisc(i,j,k)*jum(i,j,k)

nuF = ( fz(k+1) *viscTot(i,j,k+1) &
+ ( oned-fz(k+1) ) *viscTot(i,j,k) )*(1-kup(i,j,k)) &
+ bczvisc(i,j,k)*kup(i,j,k)

nuB = ( fz(k) *viscTot(i,j,k) &
+ ( oned-fz(k) ) *viscTot(i,j,k-1) )*(1-kum(i,j,k)) &
+ bczvisc(i,j,k)*kum(i,j,k)

amx = ( dxcinv(i)*(1-ium(i,j,k)) &
+dxinv(i)*ium(i,j,k)*twod )*dxinv(i)
apx = ( dxcinv(i+1)*(1-iup(i,j,k)) &
+dxinv(i)*iup(i,j,k)*twod )*dxinv(i)
amx = amx*nuW
apx = apx*nuE
acx = - ( amx + apx )

amy = ( dycinv(j)*(1-jum(i,j,k)) &
+dyinv(j)*jum(i,j,k)*twod )*dyinv(j)
apy = ( dycinv(j+1)*(1-jup(i,j,k)) &
+dyinv(j)*jup(i,j,k)*twod )*dyinv(j)
amy = amy * nuS
apy = apy * nuN
acy = - ( amy + apy )

amz = ( dzcinv(k)*(1-kum(i,j,k)) &
+dzinv(k)*kum(i,j,k)*twod )*dzinv(k)*rnDim
apz = ( dzcinv(k+1)*(1-kup(i,j,k)) &
+dzinv(k)*kup(i,j,k)*twod )*dzinv(k)*rnDim
amz = amz * nuB
apz = apz * nuF
acz = - ( amz + apz )

diffvel = amx*( vel(i-1,j,k)*(1-ium(i,j,k)) &
+ bcxvel(i,j,k)*ium(i,j,k) ) &
+ acx* vel(i,j,k) &
+ apx*( vel(i+1,j,k)*(1-iup(i,j,k)) &
+ bcxvel(i,j,k)*iup(i,j,k) ) &
+ amy*( vel(i,j-1,k)*(1-jum(i,j,k)) &
+ bcyvel(i,j,k)*jum(i,j,k) ) &
+ acy* vel(i,j,k) &
+ apy*( vel(i,j+1,k)*(1-jup(i,j,k)) &
+ bcyvel(i,j,k)*jup(i,j,k) ) &
+ amz*( vel(i,j,k-1)*(1-kum(i,j,k)) &
+ bczvel(i,j,k)*kum(i,j,k) ) &
+ acz* vel(i,j,k) &
+ apz*( vel(i,j,k+1)*(1-kup(i,j,k)) &
+ bczvel(i,j,k)*kup(i,j,k) ) 

nlvel(i,j,k) = vel(i,j,k) - dt*nlvel(i,j,k) + half*dt*diffvel
ENDDO ! i
ENDDO ! j
ENDDO ! k

END SUBROUTINE rhs_diff
!------------------------------------------------------------------------------

SUBROUTINE rhs_vankan

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE pressure_arrays
USE multiuse_arrays

IMPLICIT NONE

!... loop variables

INTEGER :: i,j,k

!... local variables

REAL(KIND=CGREAL) :: pe,pw,pn,ps,pf,pb,pgx,pgy,pgz

!******************************************************************************

!------------------------------------------------------------------------------ 
! add pressure gradient
!------------------------------------------------------------------------------

DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc

pe = ( fx(i+1)*p(i+1,j,k) + (oned-fx(i+1))*p(i,j,k) )*(1-iup(i,j,k)) &
+ p(i,j,k)*iup(i,j,k)

pw = ( fx(i) *p(i,j,k) + (oned-fx(i)) *p(i-1,j,k) )*(1-ium(i,j,k)) &
+ p(i,j,k)*ium(i,j,k)

pn = ( fy(j+1)*p(i,j+1,k) + (oned-fy(j+1))*p(i,j,k) )*(1-jup(i,j,k)) &
+ p(i,j,k)*jup(i,j,k)

ps = ( fy(j) *p(i,j,k) + (oned-fy(j)) *p(i,j-1,k) )*(1-jum(i,j,k)) &
+ p(i,j,k)*jum(i,j,k)

pf = ( fz(k+1)*p(i,j,k+1) + (oned-fz(k+1))*p(i,j,k) )*(1-kup(i,j,k)) &
+ p(i,j,k)*kup(i,j,k)

pb = ( fz(k) *p(i,j,k) + (oned-fz(k)) *p(i,j,k-1) )*(1-kum(i,j,k)) &
+ p(i,j,k)*kum(i,j,k)

pgx= (pe-pw)*dxinv(i)
pgy= (pn-ps)*dyinv(j)
pgz= (pf-pb)*dzinv(k)

nlu(i,j,k) = nlu(i,j,k) - dt*pgx*REAL(1-iblank(i,j,k),KIND=CGREAL)
nlv(i,j,k) = nlv(i,j,k) - dt*pgy*REAL(1-iblank(i,j,k),KIND=CGREAL)
nlw(i,j,k) = nlw(i,j,k) - dt*pgz*REAL(1-iblank(i,j,k),KIND=CGREAL)

ENDDO
ENDDO
ENDDO

END SUBROUTINE rhs_vankan
!------------------------------------------------------------------------------

END SUBROUTINE rhs_advec_diff 
!-------------------------------------------------------------------------------
SUBROUTINE rhs_adjust_bc()

USE global_parameters
USE flow_parameters
USE flow_arrays
USE boundary_arrays
USE grid_arrays
USE multiuse_arrays
USE solver_ad_arrays
USE GCM_arrays
USE nlold_arrays

IMPLICIT NONE

INTEGER :: i,j,k,n,iFr,jFr,kFr,iRow
REAL(KIND=CGREAL) :: amxd,apxd,acxd
REAL(KIND=CGREAL) :: amyd,apyd,acyd
REAL(KIND=CGREAL) :: amzd,apzd,aczd
REAL(KIND=CGREAL) :: rnDim
REAL(KIND=CGREAL) :: nuE,nuW,nuS,nuN,nuF,nuB

REAL(KIND=CGREAL) :: half_dt, tmp1, tmp2, tmp3 ! H. Luo

!--- H. Luo start---------------

IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
half_dt = half * dt
ELSEIF (advec_scheme == ADAMS_BASHFORTH2) THEN
half_dt = zero
ENDIF
!print*,'rhs_adjust_bc: half_dt = ', half_dt 
!--- H. Luo end-----------------

rnDim = REAL((ndim - DIM_2D),KIND=CGREAL) 

DO k=1,nzc
DO j=1,nyc
DO i=1,nxc
nuE = ( fx(i+1) *viscTot(i+1,j,k) &
+ ( oned-fx(i+1) ) *viscTot(i,j,k) )*(1-iup(i,j,k)) &
+ bcxvisc(i,j,k)*iup(i,j,k)

nuW = ( fx(i) *viscTot(i,j,k) &
+ ( oned-fx(i) ) *viscTot(i-1,j,k) )*(1-ium(i,j,k)) &
+ bcxvisc(i,j,k)*ium(i,j,k)

nuN = ( fy(j+1) *viscTot(i,j+1,k) &
+ ( oned-fy(j+1) ) *viscTot(i,j,k) )*(1-jup(i,j,k)) &
+ bcyvisc(i,j,k)*jup(i,j,k)

nuS = ( fy(j) *viscTot(i,j,k) &
+ ( oned-fy(j) ) *viscTot(i,j-1,k) )*(1-jum(i,j,k)) &
+ bcyvisc(i,j,k)*jum(i,j,k)

nuF = ( fz(k+1) *viscTot(i,j,k+1) &
+ ( oned-fz(k+1) ) *viscTot(i,j,k) )*(1-kup(i,j,k)) &
+ bczvisc(i,j,k)*kup(i,j,k)

nuB = ( fz(k) *viscTot(i,j,k) &
+ ( oned-fz(k) ) *viscTot(i,j,k-1) )*(1-kum(i,j,k)) &
+ bczvisc(i,j,k)*kum(i,j,k)

amxd = ( dxcinv(i)*(1-ium(i,j,k)) &
+dxinv(i)*ium(i,j,k)*twod )*dxinv(i)
apxd = ( dxcinv(i+1)*(1-iup(i,j,k)) &
+dxinv(i)*iup(i,j,k)*twod )*dxinv(i)

amyd = ( dycinv(j)*(1-jum(i,j,k)) &
+dyinv(j)*jum(i,j,k)*twod )*dyinv(j)
apyd = ( dycinv(j+1)*(1-jup(i,j,k)) &
+dyinv(j)*jup(i,j,k)*twod )*dyinv(j)

amzd = ( dzcinv(k)*(1-kum(i,j,k)) &
+dzinv(k)*kum(i,j,k)*twod )*dzinv(k)
apzd = ( dzcinv(k+1)*(1-kup(i,j,k)) &
+dzinv(k)*kup(i,j,k)*twod )*dzinv(k)

amxd = - (0.50_CGREAL*dt*nuW)*amxd
apxd = - (0.50_CGREAL*dt*nuE)*apxd

amyd = - (0.50_CGREAL*dt*nuS)*amyd
apyd = - (0.50_CGREAL*dt*nuN)*apyd

amzd = - (0.50_CGREAL*dt*nuB)*amzd*rnDim 
apzd = - (0.50_CGREAL*dt*nuF)*apzd*rnDim 

!--------Added by H. Luo-------------------------------------------
tmp1 = (oned - fx(i ))*(1 - ium(i,j,k)) + ium(i,j,k)
tmp2 = fx(i+1) *(1 - iup(i,j,k)) + iup(i,j,k)
amxd = amxd - half_dt * face_u(i ,j,k)* tmp1 *dxinv(i)
apxd = apxd + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i) 

tmp1 = (oned - fy(j ))*(1 - jum(i,j,k)) + jum(i,j,k)
tmp2 = fy(j+1) *(1 - jup(i,j,k)) + jup(i,j,k)
amyd = amyd - half_dt * face_v(i,j ,k)* tmp1 *dyinv(j)
apyd = apyd + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)

tmp1 = (oned - fz(k ))*(1 - kum(i,j,k)) + kum(i,j,k)
tmp2 = fz(k+1) *(1 - kup(i,j,k)) + kup(i,j,k)
amzd = amzd - half_dt * face_w(i,j,k )* tmp1 *dzinv(k)
apzd = apzd + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)
!--------End adding------------------------------------------------


!print *, amxd, apxd, amyd, apyd, amzd, apzd
!------------------------------------------------------------------------------
! Take care of fresh cells in SSM
!------------------------------------------------------------------------------

IF ( boundary_motion == MOVING_BOUNDARY .AND. &
boundary_formulation==SSM_METHOD ) THEN

amxd = amxd* REAL(1-fresh_cell(i,j,k),KIND=CGREAL) &
- REAL(fresh_cell(i,j,k),KIND=CGREAL)*ium(i,j,k)*(twod*dxcinv(i) )**sidw
apxd = apxd* REAL(1-fresh_cell(i,j,k),KIND=CGREAL) &
- REAL(fresh_cell(i,j,k),KIND=CGREAL)*iup(i,j,k)*(twod*dxcinv(i+1))**sidw

amyd = amyd* REAL(1-fresh_cell(i,j,k),KIND=CGREAL) &
- REAL(fresh_cell(i,j,k),KIND=CGREAL)*jum(i,j,k)*(twod*dycinv(j) )**sidw
apyd = apyd* REAL(1-fresh_cell(i,j,k),KIND=CGREAL) &
- REAL(fresh_cell(i,j,k),KIND=CGREAL)*jup(i,j,k)*(twod*dycinv(j+1))**sidw

amzd = amzd* REAL(1-fresh_cell(i,j,k),KIND=CGREAL) &
- REAL(fresh_cell(i,j,k),KIND=CGREAL)*kum(i,j,k)*(twod*dzcinv(k) )**sidw
apzd = apzd* REAL(1-fresh_cell(i,j,k),KIND=CGREAL) &
- REAL(fresh_cell(i,j,k),KIND=CGREAL)*kup(i,j,k)*(twod*dzcinv(k+1))**sidw

amzd = amzd*rnDim 
apzd = apzd*rnDim 

ENDIF

!print *, amxd, apxd, amyd, apyd, amzd, apzd
! if (fresh_cell(i,j,k) == 1) then
! print *, 'i,j,nlu/v=', i,j,nlu(i,j,k),nlv(i,j,k)
! endif

IF (boundary_motion_type(1) == FEA_FLOW_STRUC_INTERACTION .OR. &
boundary_motion_type(1) == PARTIAL_DYNAMICS_COUPLED .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED .OR. &
boundary_motion_type(1) == BIO_DYNAMICS_COUPLED .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED_QUAT .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED_MofI_QUAT .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED_FALLING_DEFOR .OR. &
boundary_motion_type(1) == DYNAMICS_COUPLED_SWIMMING) THEN
tmp1=nlu_FSI(i,j,k)
tmp2=nlv_FSI(i,j,k)
tmp3=nlw_FSI(i,j,k)
ELSE
tmp1=nlu(i,j,k)
tmp2=nlv(i,j,k)
tmp3=nlw(i,j,k)
ENDIF

nlu(i,j,k) = TMP1 &
- amxd*bcxu(i,j,k)*ium(i,j,k) &
- apxd*bcxu(i,j,k)*iup(i,j,k) &
- amyd*bcyu(i,j,k)*jum(i,j,k) &
- apyd*bcyu(i,j,k)*jup(i,j,k) &
- amzd*bczu(i,j,k)*kum(i,j,k) &
- apzd*bczu(i,j,k)*kup(i,j,k)

nlv(i,j,k) = TMP2 &
- amxd*bcxv(i,j,k)*ium(i,j,k) &
- apxd*bcxv(i,j,k)*iup(i,j,k) &
- amyd*bcyv(i,j,k)*jum(i,j,k) &
- apyd*bcyv(i,j,k)*jup(i,j,k) &
- amzd*bczv(i,j,k)*kum(i,j,k) &
- apzd*bczv(i,j,k)*kup(i,j,k) 

nlw(i,j,k) = TMP3 &
- amxd*bcxw(i,j,k)*ium(i,j,k) &
- apxd*bcxw(i,j,k)*iup(i,j,k) &
- amyd*bcyw(i,j,k)*jum(i,j,k) &
- apyd*bcyw(i,j,k)*jup(i,j,k) &
- amzd*bczw(i,j,k)*kum(i,j,k) &
- apzd*bczw(i,j,k)*kup(i,j,k)

ENDDO
ENDDO
ENDDO

!------------------------------------------------------------------------------
! Take care of fresh cells in GCM
!------------------------------------------------------------------------------

IF ( boundary_motion == MOVING_BOUNDARY .AND. &
boundary_formulation == GCM_METHOD .AND. nFresh > 0 ) THEN

DO n=1,nFresh
iFr = iFresh(n)
jFr = jFresh(n)
kFr = kFresh(n)
DO iRow = 1,iRowMax
i = iFreshCellIndex(n) + incI(iRow)
j = jFreshCellIndex(n) + incJ(iRow)
k = kFreshCellIndex(n) + incK(iRow)
IF ( i==iFr .AND. j==jFr .AND. k==kFr) THEN
nlu(iFr,jFr,kFr) = coeffGCMFreshD(iRow,n)*uBodyInterceptFresh(n)
nlv(iFr,jFr,kFr) = coeffGCMFreshD(iRow,n)*vBodyInterceptFresh(n)
nlw(iFr,jFr,kFr) = coeffGCMFreshD(iRow,n)*wBodyInterceptFresh(n)
ENDIF
ENDDO
ENDDO
ENDIF

END SUBROUTINE rhs_adjust_bc 
!------------------------------------------------------------------------------- 

!-------------------------------------------------------------------------------
! SUBROUTINE rhs_adjust2D() 
!-------------------------------------------------------------------------------

SUBROUTINE rhs_adjust2D() 

USE global_parameters
USE flow_parameters
USE multiuse_arrays

IMPLICIT NONE

INTEGER :: i,j,k

! Set RHS for momentum equations for 2D calculations

! copy k=1 plane to other planes

DO k = 2,nzc
DO j = 1,nyc
DO i = 1,nxc
nlu(i,j,k) = nlu(i,j,1)
nlv(i,j,k) = nlv(i,j,1)
ENDDO ! i
ENDDO ! j
ENDDO ! k

! zero w-component

DO k = 1,nzc
DO j = 1,nyc
DO i = 1,nxc
nlw(i,j,k) = zero
ENDDO ! i
ENDDO ! j
ENDDO ! k

END SUBROUTINE rhs_adjust2D
!-------------------------------------------------------------------------------
