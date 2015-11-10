 ! TURB parameters ------------------------------------------------------------

   MODULE turb_global_parameters
   
    IMPLICIT NONE

    INTEGER, PARAMETER :: TURB_NOMODEL       = 0, &  ! No Turbulence Model
                          TURB_MODEL_DYNSMAG = 1, &  ! LES dynamic Smagorinsky
                          TURB_MODEL_DYNLAGR = 2     ! LES dynamic Lagrangian
    
    INTEGER, PARAMETER :: TURB_FWIDTH_CUBEROOT = 1, &
                          TURB_FWIDTH_GEOM     = 2
    
    INTEGER, PARAMETER :: TURB_LAGR_TSCALE_TSL = 1, &
                          TURB_LAGR_TSCALE_JFM = 2
  
    INTEGER, PARAMETER :: S11 = 1, &                 ! symm tensor components
                          S12 = 2, & 
                          S13 = 3, &
                          S22 = 4, &
                          S23 = 5, &
                          S33 = 6, &
                          TENSOR_SYM_NELM = 6        ! elm number of symm tensor  
 
    INTEGER, PARAMETER :: DIRX = 1, &                
                          DIRY = 2, &
                          DIRZ = 3
 
    INTEGER, PARAMETER :: GRADX = 1, &                
                          GRADY = 2, &
                          GRADZ = 3
                          
    INTEGER, PARAMETER :: ZERO_DIR  = 0, &
                          ONE_DIR   = 1, &
                          TWO_DIR   = 2, &
                          THREE_DIR = 3
               
   END MODULE turb_global_parameters
!------------------------------------------------------------------------------

   MODULE turb_parameters
   
    USE global_parameters
    
    IMPLICIT NONE
    
    INTEGER :: turbActive
    INTEGER :: filterWidthModel, turbModel, turbLagrTimeScale
    INTEGER :: ifuRstrtTurbIn, ifuRstrtTurbOut
    INTEGER :: indexRstrtTurb
    INTEGER, DIMENSION(3) :: homogDir, testFilterDir
    
    REAL(KIND=CGREAL)               :: oneThird, twoSqrt, twoThird
    REAL(KIND=CGREAL)               :: cSmagFix
    REAL(KIND=CGREAL), DIMENSION(3) :: fWidthRatio
  
   END MODULE turb_parameters
!------------------------------------------------------------------------------

   MODULE turb_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
    REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: mij, lij, sij
    REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: qField,qTestField
    REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: phiLM, phiMM  

   END MODULE turb_arrays
!------------------------------------------------------------------------------

