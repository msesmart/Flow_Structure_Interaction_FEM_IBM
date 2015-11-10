MODULE global_parameters

IMPLICIT NONE

!INTEGER, PARAMETER :: CGREAL = SELECTED_REAL_KIND(P=15,R=307) ! double precision float 64bit
INTEGER, PARAMETER :: CGREAL =8
INTEGER(c_int) numVertices=6000 ! Added by CJ Yuan July.17.2015

REAL(KIND=CGREAL), PARAMETER :: zero = 0.0_CGREAL
REAL(KIND=CGREAL), PARAMETER :: oned = 1.0_CGREAL
REAL(KIND=CGREAL), PARAMETER :: twod = 2.0_CGREAL
REAL(KIND=CGREAL), PARAMETER :: half = 0.5_CGREAL

INTEGER, PARAMETER :: VISCOUS_FLOW = 1, &
POTENTIAL_FLOW = 2

INTEGER, PARAMETER :: UNIFORM_GRID = 1, &
NONUNIFORM_GRID = 2

INTEGER, PARAMETER :: BC_TYPE_DIRICHLET = 1, &
BC_TYPE_ZERO_GRADIENT = 2, &
BC_TYPE_PULSATILE_INFLOW = 3, &
BC_TYPE_SYMMETRY = 4, &
BC_TYPE_PERIODIC = 5, &
BC_TYPE_USER_SPECIFIED = 6, &
BC_TYPE_SHEAR = 7

INTEGER, PARAMETER :: AD_SOLVER_TYPE_LSOR = 0, &
AD_SOLVER_TYPE_MSIP = 1

INTEGER, PARAMETER :: PP_SOLVER_TYPE_LSOR = 1, &
PP_SOLVER_TYPE_PETSC = 2, &
PP_SOLVER_TYPE_MG = 3, &
PP_SOLVER_TYPE_SIP = 4, &
PP_SOLVER_TYPE_MG_SIP = 5, &
PP_SOLVER_TYPE_MSIP = 6, &
PP_SOLVER_TYPE_MG_MSIP = 7, &
PP_SOLVER_TYPE_MG_Point_Jacobi = 8

INTEGER, PARAMETER :: FIXED_BOUNDARY = 1, &
MOVING_BOUNDARY = 2

INTEGER, PARAMETER :: STATIONARY = 0, &
FORCED = 1, &
FLOW_INDUCED = 2, &
PRESCRIBED = 3, &
FEA_FLOW_STRUC_INTERACTION = 4, & !Added by Wanh by CJ Yuan
PARTIAL_DYNAMICS_COUPLED = 5, & !Added by Wanh
DYNAMICS_COUPLED = 6, & !Added by Wanh
BIO_DYNAMICS_COUPLED = 7, & !Added by Yan
BIO_FOLLOWED_DYNAMICS_COUPLED = 8, & !Added by Yan
DYNAMICS_COUPLED_QUAT = 9, & !Added by Geng (for rigid falling body)
DYNAMICS_COUPLED_MofI_QUAT = 10, & !Added by Geng (for rigid falling body)
DYNAMICS_COUPLED_FALLING_DEFOR = 11, & !Added by Geng (for falling body with active deforming motion)
DYNAMICS_COUPLED_SWIMMING = 12 !Added by Geng (for free swimming, 1 body, either solid body or membrane)

INTEGER, PARAMETER :: PBC_DIRICHLET = 1, &!added by H. Luo
PBC_NEUMANN = 2 !

INTEGER, PARAMETER :: ADAMS_BASHFORTH2 = 1, &!added by H. Luo
CRANK_NICOLSON1 = 2, &!
CRANK_NICOLSON2 = 3 !

INTEGER, PARAMETER :: NONPOROUS_AND_NONSLIP = 0, &
POROUS_OR_SLIP = 1

INTEGER, PARAMETER :: NONE = 0, &
GENERAL = 1, &
CANONICAL = 2

INTEGER, PARAMETER :: ELLIPTIC_CYLINDER = 1, &
GENERAL_CYLINDER = 2, &
ELLIPSOID = 3, &
UNSTRUCTURED_SURFACE = 4

INTEGER, PARAMETER :: OPEN_MEMBRANE = 1, &
CLOSED_MEMBRANE = 2

INTEGER, PARAMETER :: SOLID_BODY = 1, &
MEMBRANE = 2

INTEGER, PARAMETER :: BODY_DIM2 = 2, &
BODY_DIM3 = 3

INTEGER, PARAMETER :: INTR_BOUND_NONE = 0, &
INTR_BOUND_PRESENT = 1

INTEGER, PARAMETER :: NO_VAN_KAN = 0, &
VAN_KAN = 1

INTEGER, PARAMETER :: TECPLOT = 1, &
FIELDVIEW = 2

INTEGER, PARAMETER :: IBLANK_READ = 1

INTEGER, PARAMETER :: DIM_2D = 2, &
DIM_3D = 3

INTEGER, PARAMETER :: IBLANK_USED = 1

INTEGER, PARAMETER :: NO_INTERNAL_BOUNDARY = 0, &
SSM_METHOD = 1, &
GCM_METHOD = 2

INTEGER, PARAMETER :: INVISCID = 1

INTEGER, PARAMETER :: ICOORD = 1, &
JCOORD = 2, &
KCOORD = 3

INTEGER, PARAMETER :: STATS_NONE = 0

INTEGER, PARAMETER :: STDOUT = 6

INTEGER, PARAMETER :: INACTIVE = 0, &
ACTIVE = 1, &
ERR_NONE = 0

INTEGER, PARAMETER :: OUTPUT_IBLANK = 1, &
OUTPUT_FRESH_CELL = 2, &
OUTPUT_GHOSTCELLMARK = 4, &
OUTPUT_BODYNUM = 8, &
OUTPUT_IUM = 16, &
OUTPUT_IUP = 32, &
OUTPUT_JUM = 64, &
OUTPUT_JUP = 128, &
OUTPUT_KUM = 256, &
OUTPUT_KUP = 512, &
OUTPUT_GHOSTCELLSOLID =1024, &
OUTPUT_GHOSTCELLMEMB =2048, &
OUTPUT_IBLANK_SOLID =4096, &
OUTPUT_IBLANK_MEMB =8192

LOGICAL :: cure_pressure_oscillations

LOGICAL :: pressure_osc_velocity !add by Yan

LOGICAL :: pressure_osc_pressure !add by Yan

real(kind=cgreal) :: errorPermitted !add by Yan

integer :: evaluationInterval,optCount=0 !add by Yan

logical :: optStop=.false.,optimization

LOGICAL :: Full_Coarsening

REAL(KIND=CGREAL) :: sumV4=zero,sumV7=zero,sumV10=zero !add by Yan
REAL(KIND=CGREAL) :: V4Pre=zero,V4=zero,V7Pre=zero,V7=zero,V10Pre=zero,V10=zero !add by Yan

LOGICAL, PARAMETER :: Conformal_Mapping = .FALSE.

LOGICAL, PARAMETER :: pressureAitkenOn = .FALSE.

LOGICAL, PARAMETER :: Aitken = .TRUE.

TYPE :: CELL_picar3d
INTEGER :: I, J, K
INTEGER :: F_ip, F_im
INTEGER :: F_jp, F_jm
INTEGER :: F_kp, F_km
INTEGER :: F_slice
REAL(KIND=CGREAL), DIMENSION(3) :: slice_normal
REAL(KIND=CGREAL) :: volumn !, volumn_fraction
END TYPE CELL_picar3d

TYPE :: FACE_picar3d
! REAL(KIND=CGREAL) :: area
REAL(KIND=CGREAL) :: a
REAL(KIND=CGREAL) :: vel
REAL(KIND=CGREAL), DIMENSION(3) :: centroid
END TYPE FACE_picar3d

CONTAINS

FUNCTION CROSS(v1,v2)
REAL(KIND=CGREAL), DIMENSION(3) :: CROSS
REAL(KIND=CGREAL), DIMENSION(3), INTENT(IN) :: V1, V2
CROSS(1) = V1(2)*V2(3) - V1(3)*V2(2)
CROSS(2) = V1(3)*V2(1) - V1(1)*V2(3)
CROSS(3) = V1(1)*V2(2) - V1(2)*V2(1)
END FUNCTION CROSS

END MODULE global_parameters
!------------------------------------------------------

MODULE flow_parameters

USE global_parameters
use iso_c_binding
IMPLICIT NONE

REAL(KIND=CGREAL), PARAMETER :: sidw = 2.0_CGREAL ! parameter used in IDW interpolation
REAL(KIND=CGREAL) :: pi

INTEGER :: nread
INTEGER :: ndim
INTEGER :: flow_type
INTEGER :: nx,ny,nz
INTEGER :: nxc,nyc,nzc
INTEGER :: xgrid_unif,ygrid_unif,zgrid_unif
INTEGER :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
INTEGER :: no_tsteps,nmonitor,ndump,nrestart,nstat,nmonitor_probe_liftdrag
INTEGER :: STATS_SUM !added by chengyu
INTEGER :: format_dump
INTEGER :: ntime,ntime_start, ntime_skip, ntime_skip_check
INTEGER :: pp_solver_type
INTEGER :: ad_solver_type
INTEGER :: mlev,iwmg,mlw,iterMax_ad
INTEGER :: boundary_motion, nbody, nbody_solid, nbody_membrane, nGate !add by Yan
INTEGER :: nBody_Active !Added by Wanh for multi-boody in FSI
INTEGER :: boundary_formulation
INTEGER :: iterMax_Poisson
INTEGER :: body_type
INTEGER :: iterResPoisson
INTEGER :: ifuBodyIn, ifuInput, ifuIblankIn, ifuMarkerIn, ifuProbeIn, &
ifuRstrtFlowIn, ifuRstrtBodyIn, ifuUnstrucSurfIn, ifuUnstrucSurfOut, &
ifuBodyOut
INTEGER :: ifuDragOut, ifuFreshCellOut, ifuMarkerOut, ifuProbeOut, &
ifuRstrtFlowOut , ifuRstrtBodyOut, ifuStatOut, ifuStatPlot, ifuPowerMembrane2D
INTEGER :: ifuDragOutHinge !Added by Wanh
INTEGER :: ifuRstrtDynamicsIn !Added by Wanh
INTEGER :: ifuRstrtDynamicsOut !Added by Wanh
INTEGER :: ifuDragOutZone !Added by Chengyu
INTEGER :: ifuPrsbMomentRef
INTEGER :: ifuDynaNonIner !Added by Geng
INTEGER :: ifuforce_debug !Added by Geng for debug

INTEGER :: ifuFEA_input, ifuFEA_Bound_Cond
INTEGER :: ifuBodyCGPath !Added by Wanh for partially dynamic coupling
INTEGER :: internal_boundary_present,nPtsMax
INTEGER :: indexRstrt, indexStat, indexStatVort
INTEGER :: frac_step_type
INTEGER :: BinaryOutput, AdditionalOutput
LOGICAL, DIMENSION(:), ALLOCATABLE :: Fort_Formatted
INTEGER :: OUTPUT_PARA(32)

REAL(KIND=CGREAL) :: ktotalPre,kRigidPre
REAL(KIND=CGREAL) :: cpwInertiaPre,cpwRigidPre,cpwAeroPre,cpwAero

INTEGER, DIMENSION(:), ALLOCATABLE :: nPtsBodyMarkerOrig,canonical_body_type, &
body_dim,boundary_motion_type,wall_type, &
unstruc_surface_type,&
bcTypeGate !add by Yan

INTEGER, DIMENSION(:), ALLOCATABLE :: nPtsBodyMarker,rigidRef1,rigidRef2,rigidRef3

INTEGER, DIMENSION(:), ALLOCATABLE :: n_theta,n_phi

INTEGER, DIMENSION(:), ALLOCATABLE :: ntimePerCycle ! added by Haibo

INTEGER, DIMENSION(:), ALLOCATABLE :: zoneMax !added by Chengyu

INTEGER, ALLOCATABLE :: DoF_on(:,:)  ! added by G. Liu (read it from Canonical_body_in.dat)

LOGICAL :: readIblankFlag
LOGICAL :: channel_flow,dryRun,ibkOut !add by Yan
LOGICAL :: zoneSeparate !added by Chengyu

REAL(KIND=CGREAL) :: xout,yout,zout
REAL(KIND=CGREAL) :: uinit,vinit,winit
REAL(KIND=CGREAL) :: ux1,ux2,vx1,vx2,wx1,wx2
REAL(KIND=CGREAL) :: uy1,uy2,vy1,vy2,wy1,wy2
REAL(KIND=CGREAL) :: uz1,uz2,vz1,vz2,wz1,wz2
REAL(KIND=CGREAL) :: freq_ux1,freq_vx1,freq_wx1
REAL(KIND=CGREAL) :: freq_ux2,freq_vx2,freq_wx2
REAL(KIND=CGREAL) :: freq_uy1,freq_vy1,freq_wy1
REAL(KIND=CGREAL) :: freq_uy2,freq_vy2,freq_wy2
REAL(KIND=CGREAL) :: freq_uz1,freq_vz1,freq_wz1
REAL(KIND=CGREAL) :: freq_uz2,freq_vz2,freq_wz2
REAL(KIND=CGREAL) :: tstart_gust !Added by Wanh

REAL(KIND=CGREAL) :: re,dt,reinv,dtinv
REAL(KIND=CGREAL) :: restol_ad,restol_Poisson,omega, omega_adv
REAL(KIND=CGREAL) :: time
REAL(KIND=CGREAL) :: bodyInterceptWeight, imagePointWeight, probeLengthNormalized, &
gcmFlag

REAL(KIND=CGREAL) :: areax1,areax2, &
areay1,areay2, &
areaz1,areaz2

REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: radiusx,radiusy,radiusz
REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: alpha,cosalpha,sinalpha
REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: phase,cosphase,sinphase

REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: vxcent,vycent,vzcent, &
angvx,angvy,angvz, &
xcent,ycent,zcent
! for new motion -- Added by Haibo
REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angvx_old,angvy_old,angvz_old

REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: xcentinit,ycentinit,zcentinit, &
vxcentTrans,vycentTrans,vzcentTrans, &
ampx,ampy,ampz,freqx,freqy,freqz

REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angvxinit,angvyinit,angvzinit, &
ampangx,ampangy,ampangz, &
freqangx,freqangy,freqangz
REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: bcGateU,bcGateV,bcGateW,freqGate

REAL(KIND=CGREAL) :: area_left,area_right, &
area_bot,area_top, &
area_back,area_front, &
outflow_area

REAL(KIND=CGREAL) :: alfa ! Weighting factor for hybrid scheme - added by Rupesh
LOGICAL :: Hybrid

REAL(KIND=CGREAL) :: vper ! for 2D-3D initial random perturbations

! For Flow-Induced Motion. --- Added by veera

REAL(KIND=CGREAL) ,DIMENSION(:),ALLOCATABLE :: xcentConstr,ycentConstr,zcentConstr ! Centroid Constraint Flag

REAL(KIND=CGREAL) :: density_fluid

REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: density_solid

INTEGER :: pbcx1,pbcx2, pbcy1,pbcy2, pbcz1,pbcz2 ! Added by H. Luo
REAL(KIND=CGREAL) :: pppx1,pppx2, pppy1,pppy2, pppz1,pppz2 !
INTEGER :: advec_scheme ! added for implicit scheme

! For reading in body center info. for prescribed motion
LOGICAL :: Prsb_MomentRef !T for reading in moment ref point (e.g. body center location) from input_momentref.dat
REAL(KIND=CGREAL) :: Moment_refx, Moment_refy, Moment_refz

!REAL(c_double), DIMENSION(:) , ALLOCATABLE ::bodyMarkerForce,bodyMarkerVel ! Added by CJ Yuan Jun.2.2015
!REAL(c_double), DIMENSION(:) , ALLOCATABLE ::markerInterpolateRatio,markerPressure ! Added by CJ Yuan July.13.2015
!REAL(c_double), DIMENSION(:) , ALLOCATABLE ::markerInterpolateVelocity ! Added by CJ Yuan July.13.2015
!INTEGER (c_int), DIMENSION(:), ALLOCATABLE ::markerInterpolateIndex ! Added by CJ Yuan July.13.2015

END MODULE flow_parameters
!------------------------------------------------------

MODULE grid_arrays

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: x,y,z,xc,yc,zc
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dx,dy,dz
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxinv,dyinv,dzinv
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxc,dyc,dzc
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxcinv,dycinv,dzcinv
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: fx,fy,fz

END MODULE grid_arrays
!------------------------------------------------------

MODULE flow_arrays

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: u,v,w
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: u_bak,v_bak,w_bak !add by yan
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: face_u,face_v,face_w
! REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Usign,Vsign,Wsign ! For 2nd Upwinding - Added by Rupesh
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxu,bcxv,bcxw
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcyu,bcyv,bcyw
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bczu,bczv,bczw
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: viscTot
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxvisc,bcyvisc,bczvisc

INTEGER :: CUTTING_FACE_N
INTEGER :: CUTTING_EDGE_N

INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: iex, iey, iez
INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: ifx, ify, ifz

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: edge_CUT !, area_face
END MODULE flow_arrays
!------------------------------------------------------

MODULE boundary_arrays

USE global_parameters

IMPLICIT NONE

INTEGER :: nFresh
INTEGER, DIMENSION(:), ALLOCATABLE :: num_iblank
INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: zoneMarker !added by Chengyu

INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: iblank,bcBlank
INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: iblank_solid, iblank_memb, gateTest
INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: ghostCellMark
INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: hybridMarkMemb
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: conflictCell,conflictBCi,conflictBCj,conflictBCk
INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: ghostCellMemb, ghostCellSolid

INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: hybrid_mark !add by yan

INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: fresh_cell

INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: iup,ium,jup,jum,kup,kum
INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: iMarkM,iMarkP,jMarkM,jMarkP,kMarkM,kMarkP
INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: iupp,iumm,jupp,jumm,kupp,kumm ! For 2nd Upwinding - Added by Rupesh
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: exp_weight
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: pot_flag

REAL(KIND=CGREAL), DIMENSION(:,:) , ALLOCATABLE :: xBodyMarker,yBodyMarker,zBodyMarker, &
uBodyMarker,vBodyMarker,wBodyMarker, &
uBodyMarkerDyn,vBodyMarkerDyn,wBodyMarkerDyn, &! Added by Yan, (in DYNAMICS_COUPLED_FALLING_DEFOR, they are deforming velocity in inertia frame)
uBodyMarkerDefor,vBodyMarkerDefor,wBodyMarkerDefor, & ! Added by G. Liu, deforming velocity in non-inertia frame
uBodyMarkerDefor_prvs,vBodyMarkerDefor_prvs,wBodyMarkerDefor_prvs, & ! Added by G. Liu
uBodyMarkerDefor_prvs_iner,vBodyMarkerDefor_prvs_iner,wBodyMarkerDefor_prvs_iner, & ! Added by G. Liu, the deforming velocity of previous time step in inertia frame
axBodyMarker,ayBodyMarker,azBodyMarker, &
sBodyMarker,dsBodyMarker,xNormBodyMarker, &
yNormBodyMarker, &
xBodyMarkerInit,yBodyMarkerInit,zBodyMarkerInit, &
xBodyMarkerNoniner,yBodyMarkerNoniner,zBodyMarkerNoniner

Integer, DIMENSION(:,:) , ALLOCATABLE :: gateLabel

REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: pgradx1,pgradx2
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: pgrady1,pgrady2
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: pgradz1,pgradz2

INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: iblank_IN
INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: ivc

INTEGER :: CUTTING_CELL_N

TYPE (FACE_picar3d), DIMENSION(:), POINTER :: face
TYPE (CELL_picar3d), DIMENSION(:), POINTER :: cell
! TYPE :: CELL_VOLUMN_
! INTEGER :: I, J, K
! REAL(KIND=CGREAL) :: VOLUMN
! END TYPE CELL_VOLUMN_
!
! TYPE (CELL_VOLUMN_), DIMENSION(:), POINTER :: CELL_V
END MODULE boundary_arrays
!------------------------------------------------------

MODULE multiuse_arrays

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nlu,nlv,nlw,nlu0,nlw0
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: uTilde,vTilde,wTilde

REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: LU, CA
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: LU2D, CA2D

END MODULE multiuse_arrays
!------------------------------------------------------

MODULE pressure_arrays

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: p, pPrime

END MODULE pressure_arrays
!------------------------------------------------------

MODULE nlold_arrays

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nluOld,nlvOld,nlwOld
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nlu_FSI, &
nlv_FSI, &
nlw_FSI

END MODULE nlold_arrays
!------------------------------------------------------

MODULE solver_arrays

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amx,apx,acx
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amy,apy,acy
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amz,apz,acz
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: rhs,dummy
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: face1, face2

END MODULE solver_arrays
!------------------------------------------------------

MODULE solver_ad_arrays

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amx_ad,apx_ad
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amy_ad,apy_ad
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amz_ad,apz_ad

END MODULE solver_ad_arrays
!------------------------------------------------------

MODULE GCM_arrays

USE global_parameters

IMPLICIT NONE

INTEGER :: iRowMax, nGhost
INTEGER, DIMENSION(:) , ALLOCATABLE :: incI, incJ, incK, iPvt
INTEGER, DIMENSION(:) , ALLOCATABLE :: closestMarker, &
iGhost,jGhost,kGhost, &
iCellIndex,jCellIndex,kCellIndex

INTEGER, DIMENSION(:) , ALLOCATABLE :: iGhostP,jGhostP,kGhostP !Added 05/21/10

INTEGER, DIMENSION(:) , ALLOCATABLE :: iFresh,jFresh,kFresh &
,iFreshCellIndex,jFreshCellIndex,kFreshCellIndex
INTEGER, DIMENSION(:) , ALLOCATABLE :: closestMarkerFresh
INTEGER, DIMENSION(:) , ALLOCATABLE :: closestElementFresh
INTEGER, DIMENSION(:) , ALLOCATABLE :: iBodyRank
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: bodyNum
INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: zoneCheck !added by Chengyu
INTEGER, DIMENSION(:) , ALLOCATABLE :: iCellIndexS, jCellIndexS, kCellIndexS

REAL(KIND=CGREAL), DIMENSION(2) :: det
REAL(KIND=CGREAL), DIMENSION(:) , ALLOCATABLE :: work
REAL(KIND=CGREAL), DIMENSION(:) , ALLOCATABLE :: closestMarkerRatio, &
probeLength, &
xBodyInterceptTang, &
yBodyInterceptTang, &
zBodyInterceptTang, &
xBodyInterceptNorm, &
yBodyInterceptNorm, &
zBodyInterceptNorm, &
xBodyIntercept, &
yBodyIntercept, &
zBodyIntercept, &
uBodyIntercept, &
vBodyIntercept, &
wBodyIntercept, &
pBodyIntercept, &
dpdnBodyIntercept, &
dpdtBodyIntercept, &
xImagePoint, &
yImagePoint, &
zImagePoint, &
fluxGhost

REAL(KIND=CGREAL), DIMENSION(:) , ALLOCATABLE :: xBIG,yBIG,ZBIG

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: closestMarkerRatioFresh, &
xBodyInterceptFresh, &
yBodyInterceptFresh, &
zBodyInterceptFresh, &
uBodyInterceptFresh, &
vBodyInterceptFresh, &
wBodyInterceptFresh

REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: coeffGCMD, coeffGCMN, &
vanMatrixD, vanMatrixN, &
coeffGCMFreshD,dSFaceProject, &
xBodyCentroid,yBodyCentroid, &
sBodyCentroid, &
xCentroidTang,yCentroidTang, &
xCentroidNorm,yCentroidNorm

REAL(KIND=CGREAL), DIMENSION(:) , ALLOCATABLE :: probeLengthS,probeLengthNormalizedS, &
imagePointWeightS, &
xImagePointS,yImagePointS,zImagePointS

REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: coeffGCMDS, coeffGCMNS, &
vanMatrixDS, vanMatrixNS

END MODULE GCM_arrays
!------------------------------------------------------
module hybrid_cell_arrays
USE global_parameters
implicit none

INTEGER :: iRMax, nhybrid
INTEGER, DIMENSION(:) , ALLOCATABLE :: closestMarker, &
ihybrid,jhybrid,khybrid, &
iCIndex,jCIndex,kCIndex
INTEGER, DIMENSION(:), ALLOCATABLE :: closestElementHC,closestElementR1,closestElementR2,closestElementR3
REAL(KIND=CGREAL), DIMENSION(:) , ALLOCATABLE :: xBIntercept, &
yBIntercept, &
zBIntercept, &
uBIntercept, &
vBIntercept, &
wBIntercept
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: coeffD, coeffN

REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: vanMD, vanMN, &
coeffFreshD,dSFProject, &
xBCentroid,yBCentroid, &
sBCentroid, &
xCenTang,yCenTang, &
xCenNorm,yCenNorm
INTEGER, DIMENSION(:) , ALLOCATABLE :: inccI, inccJ, inccK, iPvtt
REAL(KIND=CGREAL), DIMENSION(:) , ALLOCATABLE :: workk

end module hybrid_cell_arrays
!------------------------------------------------------

MODULE unstructured_surface_arrays

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemNormX,triElemNormY, &
triElemNormZ,triElemArea
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemCentX,triElemCentY, &
triElemCentZ
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemTang1X,triElemTang1Y, &
triElemTang1Z
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemTang2X,triElemTang2Y, &
triElemTang2Z
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: pointOutsideBodyX,pointOutsideBodyY, &
pointOutsideBodyZ,surfArea

INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: triElemNeig
INTEGER, DIMENSION(:), ALLOCATABLE :: totNumTriElem
INTEGER, DIMENSION(:), ALLOCATABLE :: closestElementGC,cElementG

REAL (KIND=CGREAL) :: normDirFlag

END MODULE unstructured_surface_arrays
!------------------------------------------------------

MODULE probe_parameters

IMPLICIT NONE

INTEGER :: nProbe
INTEGER, DIMENSION(:),ALLOCATABLE :: iProbe, jProbe, kProbe

END MODULE probe_parameters
!------------------------------------------------------
!------------------------------------------------------
MODULE stat_arrays

USE global_parameters

IMPLICIT NONE

INTEGER :: statCtr
REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uAv,vAv,wAv,pAv
REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uvAv,vwAv,uwAv
REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uuAv,vvAv,wwAv

END MODULE stat_arrays
!------------------------------------------------------

MODULE blasius_profile

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL) :: cavity_H,slot_H,ddratio,d,delta,uinf
REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:) :: eta,u_blasius
INTEGER :: l,junk,i_start

END MODULE blasius_profile
!------------------------------------------------------

MODULE mg_parameters

USE global_parameters
USE flow_parameters

IMPLICIT NONE

INTEGER :: mgLevels_X, mgLevels_Y, mgLevels_Z
INTEGER :: mgcyclex, mgcycley, mgcyclez, infoconv, incrlev
INTEGER :: iterFinest, iterInter, iterCoarsest
INTEGER :: ittt1, nCount

INTEGER :: iRedBlack, TNcolorX, TNcolorY, TNcolorZ, iStep, jStep, kStep
!new for Redblack LSOR

END MODULE mg_parameters


MODULE mg_arrays

USE global_parameters

IMPLICIT NONE

INTEGER, DIMENSION(:), ALLOCATABLE :: mgrid_I, mgrid_J, mgrid_K
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ghostcellMark_MG
INTEGER(1), DIMENSION(:,:,:), POINTER :: iblank_MG

TYPE :: MGtype
INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: iblank
REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: CA
REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: LU
END TYPE MGtype

TYPE (MGtype), DIMENSION(:), ALLOCATABLE, TARGET :: MGX, MGY, MGZ

REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: rhs_mg, phi_mg

! REAL(KIND=CGREAL),DIMENSION(:,:), ALLOCATABLE :: pgradx1_MG, pgradx2_MG
! REAL(KIND=CGREAL),DIMENSION(:,:), ALLOCATABLE :: pgrady1_MG, pgrady2_MG
! REAL(KIND=CGREAL),DIMENSION(:,:), ALLOCATABLE :: pgradz1_MG, pgradz2_MG

INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: ium_mg, iup_mg, &
jum_mg, jup_mg, kum_mg, kup_mg

REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: &
dxcinv_mg,dycinv_mg,dzcinv_mg
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: &
dxinv_mg,dyinv_mg,dzinv_mg
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: &
dxc_mg,dyc_mg,dzc_mg
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: &
dx_mg,dy_mg,dz_mg
REAL(KIND=CGREAL),DIMENSION(:,:), ALLOCATABLE :: x_mg, xc_mg, y_mg, &
yc_mg,z_mg, zc_mg

END MODULE mg_arrays

!------------------------------------------------------
MODULE usr_module

USE global_parameters

REAL(KIND=CGREAL) :: density_ratio
REAL(KIND=CGREAL) :: lScale,vScale
REAL(KIND=CGREAL) :: non_dim_volume,volume
REAL(KIND=CGREAL) :: I_XX,I_YY,I_ZZ,I_XY,I_YZ,I_XZ
! REAL(KIND=CGREAL) :: vxcent_prev,vycent_prev,vzcent_prev
REAL(KIND=CGREAL) :: angvx_prev, angvy_prev, angvz_prev
REAL(KIND=CGREAL) :: moment_x,moment_y,moment_z
REAL(KIND=CGREAL) :: force_x,force_y,force_z
REAL(KIND=CGREAL) :: moment_z_pretime
REAL(KIND=CGREAL),DIMENSION(1:3,1:3) ::nonDimM_I, invMI,nonDimM_I_prvs
REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: scx,scy,scz
REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: scmx,scmy,scmz
! Added by Wanh for FSI
REAL(KIND=CGREAL),DIMENSION(:,:),ALLOCATABLE :: uBodyMarker_iter,vBodyMarker_iter,wBodyMarker_iter
REAL(KIND=CGREAL),DIMENSION(:,:),ALLOCATABLE :: uBodyMarker_worelax,vBodyMarker_worelax,wBodyMarker_worelax
! Added by Wanh for partially dynamic coupling
REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: vxcent_prev,vycent_prev,vzcent_prev !Changed by Wanh
! REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: angvx_prev, angvy_prev, angvz_prev !Changed by Wanh
REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: vxcent_iter,vycent_iter,vzcent_iter !Added by Wanh
REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: angvx_iter, angvy_iter, angvz_iter !Added by Wanh
REAL(KIND=CGREAL) :: vxcent_wo_relax,vycent_wo_relax,vzcent_wo_relax
REAL(KIND=CGREAL) :: angvx_wo_relax,angvy_wo_relax,angvz_wo_relax

END MODULE usr_module

!------------------------------------------------------
MODULE stat_vort_arrays

USE global_parameters

IMPLICIT NONE

INTEGER :: statCtrv
REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: oxAv,oyAv,ozAv
REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: oxoxAv,oyoyAv,ozozAv

END MODULE stat_vort_arrays

!------------------------------------------------------
! Added for structure by Wanh 05/05/10
MODULE fea_unstructure_surface
USE global_parameters

LOGICAL :: FSI_on
LOGICAL :: FSI_CONVERGE
LOGICAL :: TRNS
INTEGER :: FSI_ITERATION
REAL(KIND=CGREAL) :: FSIConvg_Criteria = 1e-3

INTEGER, DIMENSION(:),ALLOCATABLE :: JBC
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: PROPERTY !At most 2 type of materials
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: LOAD,dLOAD,WK
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: ACC,VEL,DISP
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: struc_disp,struc_vel
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: struc_disp_iter,struc_vel_iter
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: struc_disp_worelax, struc_vel_worelax
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: struc_olddisp

INTEGER, DIMENSION(:), ALLOCATABLE :: ELMatType
INTEGER, DIMENSION(:), ALLOCATABLE :: nPrescribedMarker
INTEGER, DIMENSION(:,:), ALLOCATABLE :: PrescribedMarker

INTEGER :: ELETYPE, NMATP
INTEGER :: NBC
INTEGER :: NEQ,IBAND,MAXSTIFF
INTEGER :: nodeDoF=6

INTEGER, DIMENSION(:), ALLOCATABLE :: IPROF,IPROF2,NLOC
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: STF,MSS
REAL(KIND=CGREAL), DIMENSION(:,:) , ALLOCATABLE :: pBodyMarker,pxBodyMarker, &
pyBodyMarker,pzBodyMarker

REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemP
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: yBodymarkerOld(:,:)

REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: xBodyMarker0,yBodyMarker0,zBodyMarker0

REAL(KIND=CGREAL) :: massratio

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: Lambdau_Aitken,Lambdav_Aitken,Lambdaw_Aitken

REAL(KIND=CGREAL), PARAMETER :: Lambdau_Aitken_Init = 0.3,Lambdav_Aitken_Init = 0.3,Lambdaw_Aitken_Init = 0.3
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: deltau_Aitken, deltav_Aitken, deltaw_Aitken
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: deltau_Aitken_prev,deltav_Aitken_prev,deltaw_Aitken_prev

END MODULE fea_unstructure_surface


! Added by Wanh
MODULE body_dynamics

USE global_parameters

IMPLICIT NONE

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: Ipaxis_x, Ipaxis_y, Ipaxis_z
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: Rho_fs !density ratio fluid/solid
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: xCG,yCG,zCG
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: uCG,vCG,wCG
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: acc_xCG,acc_yCG,acc_zCG
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: omg_x,omg_y,omg_z
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: DynamicCharLength

REAL(KIND=CGREAL), PARAMETER :: GRAVITY_ACC = 9.81
REAL(KIND=CGREAL), PARAMETER :: Convg_Cretia = 5e-4, Convg_Cretia_Ang = 5e-3

REAL(KIND=CGREAL) :: thickoverlength !for 2D body only
REAL(KIND=CGREAL) :: DepthOverLength !for 3D body only

REAL(KIND=CGREAL),DIMENSION(:,:),ALLOCATABLE :: section_vxcent,section_vycent,section_vzcent, &
section_angvx,section_angvy,section_angvz, &
section_xcent,section_ycent,section_zcent
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: LambdaVel_Aitken,LambdaAngVel_Aitken
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: deltaVxcent_Aitken,deltaVycent_Aitken,deltaVzcent_Aitken
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: deltaAngvx_Aitken,deltaAngvy_Aitken,deltaAngvz_Aitken

REAL(KIND=CGREAL) :: Char_I !Moment of Inertia
REAL(KIND=CGREAL) :: Fn
REAL(KIND=CGREAL) :: Rho_f, Rho_s
REAL(KIND=CGREAL) :: FreqN_torsion, Freq_prsb
REAL(KIND=CGREAL), PARAMETER :: LambdaVel_Aitken_Init = 0.2, LambdaAngVel_Aitken_Init = 0.2
REAL(KIND=CGREAL) :: Quat_init(4), Quat_phys(4),Quat_iter(4),Quat_prev(4),quat_trans(4)

REAL(KIND=CGREAL) :: aero_moment,grav_moment,aero_moment_threshold,grav_moment_threshold
REAL(KIND=CGREAL) :: uref_threshold
REAL(KIND=CGREAL) :: accx_hinge, accy_hinge
REAL(KIND=CGREAL) :: A_x, A_theta

REAL(KIND=CGREAL) :: Limiter_angle_amp

INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: DynamicMarker !nBody,nSection,markers
INTEGER, DIMENSION(:,:), ALLOCATABLE :: SectionMarker !nBody,nSection
INTEGER :: nSection, hingemarker
INTEGER :: niterFS
INTEGER, PARAMETER :: niterFS_max = 30
! INTEGER, PARAMETER :: AngvIntg = 3 !Used for scheme selection of angular velocity integration.
INTEGER :: AngvIntg = 3 !Used for scheme selection of angular velocity integration.

LOGICAL, DIMENSION(:), ALLOCATABLE :: Converged_FSI
LOGICAL :: Limiter_On, Limiter_change_direction

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: mass_in ! mass_in(nbody) Added by Geng
REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: Icm_in ! Icm_in(nbody,6) Added by Geng
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: xcent_prev,ycent_prev,zcent_prev ! xcent_prev(nbody) Added by Geng
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: angvx_noniner,angvy_noniner,angvz_noniner ! angvx_noniner(nbody) Added by Geng

END MODULE body_dynamics


! For Aitken acceleration for pressure
MODULE Pressure_Aitken_Array

USE global_parameters

REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: deltapPrime_Lplus1
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: deltapPrime_prevL
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: diff_deltapPrimeL
REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: pPrime_L
REAL(KIND=CGREAL) :: LambdaP
REAL(KIND=CGREAL), PARAMETER :: LambdaP_Init = 0.2

END MODULE Pressure_Aitken_Array

! Added by Wanh end


module operation
use global_parameters
contains
function mo(a,dim)
implicit none
integer :: dim
real*8 :: a(dim)
real*8 :: mo

integer :: i
real*8 :: sum

sum=0.0d0
do i=1,dim
sum=sum+a(i)*a(i)
end do
mo=sqrt(sum)
return

end function mo

!-----------------------------------
!
!-----------------------------------

real(8) function dot(a,b)
implicit none

real(8) :: a(3),b(3)

dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
return
end function dot

!-----------------------------------
!
!-----------------------------------

function qm(a,b)
implicit none

real*8 :: a(4),b(4)
real*8 :: qm(4)

qm(1)=a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4)
qm(2)=a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3)
qm(3)=a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2)
qm(4)=a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)

end function qm
!-----------------------------------
!
!-----------------------------------

function projection(a,n)
implicit none
real(8) :: a(3),n(3)
real(8) :: projection(3)

projection=cross(n,cross(a,n))
return
end function projection
!-----------------------------------
!
!-----------------------------------

real(8) function vector_angle(a,b)
implicit none
real(8) :: a(3),b(3),c,d,pi=4.0d0*atan(1.0d0)

c=dot(a,b)/(mo(a,3)*mo(b,3))
d=mo(cross(a,b),3)/(mo(a,3)*mo(b,3))

if(c<1.0d0.and.c>-1.0d0)then
vector_angle=acos(c)
else if(c>0.0d0)then
vector_angle=0.0d0
else if(c<0.0d0)then
vector_angle=pi
end if

if(d>=-pi/2.0d0.and.d<0.0d0)then
vector_angle=2.0d0*pi-vector_angle
end if

return
end function vector_angle
!-----------------------------------
!
!-----------------------------------
function min_abs(a,b)
implicit none

real*8 :: a,b
real*8 :: min_abs

if(abs(a)<=abs(b))then
min_abs=a
else
min_abs=b
end if

return

end function min_abs

!-----------------------------------
!
!-----------------------------------
function max_abs(a,b)
implicit none

real*8 :: a,b
real*8 :: max_abs

if(abs(a)>=abs(b))then
max_abs=a
else
max_abs=b
end if

return

end function max_abs




end module operation
