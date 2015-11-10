!-------------------------------------------------------------------------------

SUBROUTINE PRECONDITION_SIP_2D()
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: LW, LS, LP, UN, UE
    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: AW, AS, AP, AN, AE

    INTEGER :: I,J,K, II,JJ
    REAL(KIND=CGREAL) :: P1, P2, BETA
!    
    LS=>LU(:,:,1,1)
    LW=>LU(:,:,1,2)
    LP=>LU(:,:,1,3)
    UE=>LU(:,:,1,4)
    UN=>LU(:,:,1,5)
    
    AS=>CA(:,:,1,1)
    AW=>CA(:,:,1,2)
    AP=>CA(:,:,1,3)
    AE=>CA(:,:,1,4)
    AN=>CA(:,:,1,5)

    K=1
    
    BETA=OMEGA
    UE=zero
    UN=zero
    
    DO J=1, nyc
    JJ=J+1
    DO I=1, nxc
    II=I+1
    IF (IBLANK(I,J,K)==1) CYCLE
        LS(II,JJ)=AS(I,J)/(1.+BETA*UE(II,JJ-1))     ! b 
        LW(II,JJ)=AW(I,J)/(1.+BETA*UN(II-1,JJ))     ! c
        P1=BETA*LS(II,JJ)*UE(II,JJ-1)
        P2=BETA*LW(II,JJ)*UN(II-1,JJ)
        LP(II,JJ)=1.d0/(AP(I,J)+P1+P2-LS(II,JJ)*UN(II,JJ-1)-LW(II,JJ)*UE(II-1,JJ))
        UE(II,JJ)=(AE(I,J)-P1)*LP(II,JJ)            ! e
        UN(II,JJ)=(AN(I,J)-P2)*LP(II,JJ)            ! f
    END DO
    END DO

END SUBROUTINE PRECONDITION_SIP_2D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE PRECONDITION_SIP_3D()
!    USE global_parameters
!    USE flow_parameters
!    USE boundary_arrays
!    USE multiuse_arrays
!
END SUBROUTINE PRECONDITION_SIP_3D
!-------------------------------------------------------------------------------
   
!-------------------------------------------------------------------------------
SUBROUTINE COEFF_AD_MSIP_2D(NEW)
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_ad_arrays

    IMPLICIT NONE

    LOGICAL :: NEW
    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: AW, AS, AP, AN, AE

    INTEGER :: I,J,K, II,JJ
    
    REAL(KIND=CGREAL) :: half_dt
    REAL(KIND=CGREAL) :: tmp1, tmp2
    REAL(KIND=CGREAL) :: acx,amx,apx,acy,amy,apy

!    
    
    IF (NEW) THEN
        ALLOCATE(LU(0:NX+1,0:NY+1,1,7))
        ALLOCATE(CA(NX+1,NY+1,1,5))
    ENDIF
    
    AS=>CA(:,:,1,1)
    AW=>CA(:,:,1,2)
    AP=>CA(:,:,1,3)
    AE=>CA(:,:,1,4)
    AN=>CA(:,:,1,5)

    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = half * dt
    K=1
    DO J=1, nyc
    DO I=1, nxc
    IF (iblank(I,J,K)==1) THEN
        ae(i,j)=zero
        aw(i,j)=zero
        an(i,j)=zero
        as(i,j)=zero
        ap(i,j)=oned
    ELSE
		amx = amx_ad(i,j,k)
		apx = apx_ad(i,j,k)
		acx =- ( amx + apx )      

		amy = amy_ad(i,j,k)
		apy = apy_ad(i,j,k)
		acy =- ( amy + apy )     

		tmp1   = oned - fx(i  );   tmp2 = fx(i+1) 
		amx = amx - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
		apx = apx + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i)   

		tmp1   =         fx(i  ) *(1 - ium(i,j,k))
		tmp2   = (oned - fx(i+1))*(1 - iup(i,j,k))
		acx = acx + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i  ,j,k)*tmp1)*dxinv(i)

		tmp1   = oned - fy(j  );   tmp2   =        fy(j+1)
		amy = amy - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
		apy = apy + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)

		tmp1   =         fy(j  ) *(1 - jum(i,j,k))
		tmp2   = (oned - fy(j+1))*(1 - jup(i,j,k))
		acy = acy + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j  ,k)*tmp1)*dyinv(j)

        ae(i,j)=apx*(1 - iup(i,j,k) )
        aw(i,j)=amx*(1 - ium(i,j,k) )
        an(i,j)=apy*(1 - jup(i,j,k) )
        as(i,j)=amy*(1 - jum(i,j,k) )
        ap(i,j)=oned + acx + acy

    ENDIF
    
    END DO
    END DO

    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
    K=1
    DO J=1, nyc
    DO I=1, nxc
    IF (iblank(I,J,K)==1) THEN
      ae(i,j)=zero
      aw(i,j)=zero
      an(i,j)=zero
      as(i,j)=zero
      ap(i,j)=oned
    ELSE
		  amx = amx_ad(i,j,k)
		  apx = apx_ad(i,j,k)
		  acx =- ( amx + apx )      

		  amy = amy_ad(i,j,k)
		  apy = apy_ad(i,j,k)
		  acy =- ( amy + apy )     

      ae(i,j)=apx*(1 - iup(i,j,k) )
      aw(i,j)=amx*(1 - ium(i,j,k) )
      an(i,j)=apy*(1 - jup(i,j,k) )
      as(i,j)=amy*(1 - jum(i,j,k) )
      ap(i,j)=oned + acx + acy

    ENDIF
    
    END DO
    END DO

    ENDIF


END SUBROUTINE COEFF_AD_MSIP_2D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE COEFF_AD_MSIP_3D(NEW)
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_ad_arrays

    IMPLICIT NONE

    LOGICAL :: NEW
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER :: I,J,K, II,JJ,KK

    REAL(KIND=CGREAL) :: half_dt
    REAL(KIND=CGREAL) :: tmp1, tmp2
    REAL(KIND=CGREAL) :: acx,amx,apx,acy,amy,apy,acz,amz,apz

!
    
    IF (NEW) THEN
        ALLOCATE(LU(0:NX+1,0:NY+1,0:NZ+1,13))
        ALLOCATE(CA(NX+1,NY+1,NZ+1,7))
    ENDIF

    AS=>CA(:,:,:,1)
    AW=>CA(:,:,:,2)
    AP=>CA(:,:,:,3)
    AE=>CA(:,:,:,4)
    AN=>CA(:,:,:,5)
    AF=>CA(:,:,:,6)
    AB=>CA(:,:,:,7)

    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = half * dt
    DO K=1, nzc
    DO J=1, nyc
    DO I=1, nxc
    IF (iblank(I,J,K)==1) THEN
        ae(i,j,K)=zero
        aw(i,j,K)=zero
        an(i,j,K)=zero
        as(i,j,K)=zero
        af(i,j,K)=zero
        ab(i,j,K)=zero
        ap(i,j,K)=oned
    ELSE
		amx = amx_ad(i,j,k)
		apx = apx_ad(i,j,k)
		acx =- ( amx + apx )      

		amy = amy_ad(i,j,k)
		apy = apy_ad(i,j,k)
		acy =- ( amy + apy )     

		amz = amz_ad(i,j,k)
		apz = apz_ad(i,j,k)
		acz =- ( amz + apz )    

		tmp1   = oned - fx(i  );   tmp2 = fx(i+1) 
		amx = amx - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
		apx = apx + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i)   

		tmp1   =         fx(i  ) *(1 - ium(i,j,k))
		tmp2   = (oned - fx(i+1))*(1 - iup(i,j,k))
		acx = acx + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i  ,j,k)*tmp1)*dxinv(i)

		tmp1   = oned - fy(j  );   tmp2   =        fy(j+1)
		amy = amy - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
		apy = apy + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)

		tmp1   =         fy(j  ) *(1 - jum(i,j,k))
		tmp2   = (oned - fy(j+1))*(1 - jup(i,j,k))
		acy = acy + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j  ,k)*tmp1)*dyinv(j)

		tmp1   = oned - fz(k  );   tmp2   =        fz(k+1)
		amz = amz - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
		apz = apz + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)

		tmp1   =         fz(k  ) *(1 - kum(i,j,k))
		tmp2   = (oned - fz(k+1))*(1 - kup(i,j,k))
		acz = acz + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k  )*tmp1)*dzinv(k) 

        ae(i,j,K)=apx*(1 - iup(i,j,k) )
        aw(i,j,K)=amx*(1 - ium(i,j,k) )
        an(i,j,K)=apy*(1 - jup(i,j,k) )
        as(i,j,K)=amy*(1 - jum(i,j,k) )
        af(i,j,K)=apz*(1 - kup(i,j,k) )
        ab(i,j,K)=amz*(1 - kum(i,j,k) )
        ap(i,j,K)=oned + acx + acy + acz

    ENDIF
    
    END DO
    END DO
    END DO

    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN

    DO K=1, nzc
    DO J=1, nyc
    DO I=1, nxc
    IF (iblank(I,J,K)==1) THEN
        ae(i,j,K)=zero
        aw(i,j,K)=zero
        an(i,j,K)=zero
        as(i,j,K)=zero
        af(i,j,K)=zero
        ab(i,j,K)=zero
        ap(i,j,K)=oned
    ELSE
		amx = amx_ad(i,j,k)
		apx = apx_ad(i,j,k)
		acx =- ( amx + apx )      

		amy = amy_ad(i,j,k)
		apy = apy_ad(i,j,k)
		acy =- ( amy + apy )     

		amz = amz_ad(i,j,k)
		apz = apz_ad(i,j,k)
		acz =- ( amz + apz )    

        ae(i,j,K)=apx*(1 - iup(i,j,k) )
        aw(i,j,K)=amx*(1 - ium(i,j,k) )
        an(i,j,K)=apy*(1 - jup(i,j,k) )
        as(i,j,K)=amy*(1 - jum(i,j,k) )
        af(i,j,K)=apz*(1 - kup(i,j,k) )
        ab(i,j,K)=amz*(1 - kum(i,j,k) )
        ap(i,j,K)=oned + acx + acy + acz

    ENDIF
    
    END DO
    END DO
    END DO

    ENDIF


END SUBROUTINE COEFF_AD_MSIP_3D
!-------------------------------------------------------------------------------
   
!-------------------------------------------------------------------------------
SUBROUTINE COEFF_PRESSURE_MSIP_2D()
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: AW, AS, AP, AN, AE

    INTEGER :: I,J,K, II,JJ
    REAL(KIND=CGREAL) :: PHI1, PHI4, BETA
!    
    
    ALLOCATE(LU(0:NX+1,0:NY+1,1,7))
    ALLOCATE(CA(NX+1,NY+1,1,5))
    
    AS=>CA(:,:,1,1)
    AW=>CA(:,:,1,2)
    AP=>CA(:,:,1,3)
    AE=>CA(:,:,1,4)
    AN=>CA(:,:,1,5)

    K=1
    
    DO J=1, nyc
    DO I=1, nxc
    IF (iblank(I,J,K)==1) THEN
        ae(i,j)=zero
        aw(i,j)=zero
        an(i,j)=zero
        as(i,j)=zero
        ap(i,j)=oned
    ELSE
        ae(i,j)=dxcinv(i+1)*dxinv(i)
        aw(i,j)=dxcinv(i)  *dxinv(i)
        an(i,j)=dycinv(j+1)*dyinv(j)
        as(i,j)=dycinv(j)  *dyinv(j)
        ap(i,j)=-( ae(i,j) + aw(i,j) + an(i,j) +  as(i,j))
        
        ap(i,j)=ap(i,j) + ae(i,j)*iup(i,j,k) + aw(i,j)*ium(i,j,k) &
                        + an(i,j)*jup(i,j,k) + as(i,j)*jum(i,j,k)
        ae(i,j)=ae(i,j)*(1 - babs(iup(i,j,k)) )
        aw(i,j)=aw(i,j)*(1 - babs(ium(i,j,k)) )
        an(i,j)=an(i,j)*(1 - babs(jup(i,j,k)) )
        as(i,j)=as(i,j)*(1 - babs(jum(i,j,k)) )
    ENDIF
    END DO
    END DO

END SUBROUTINE COEFF_PRESSURE_MSIP_2D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE COEFF_PRESSURE_MSIP_3D()
!    USE global_parameters
!    USE flow_parameters
!    USE grid_arrays
!    USE boundary_arrays
!    USE multiuse_arrays
!
!    IMPLICIT NONE
!
!    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: AW, AS, AP, AN, AE
!
!    INTEGER :: I,J,K, II,JJ
!    REAL(KIND=CGREAL) :: PHI1, PHI4, BETA
!!
END SUBROUTINE COEFF_PRESSURE_MSIP_3D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE COEFF_LU_MSIP_2D()
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: B, C, D, E, F, G, H
    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: AW, AS, AP, AN, AE

    INTEGER :: I,J,K, II,JJ
    REAL(KIND=CGREAL) :: PHI1, PHI4, BETA
!    
    B=>LU(:,:,1,1)
    C=>LU(:,:,1,2)
    D=>LU(:,:,1,3)
    E=>LU(:,:,1,4)
    F=>LU(:,:,1,5)
    G=>LU(:,:,1,6)
    H=>LU(:,:,1,7)
    
    AS=>CA(:,:,1,1)
    AW=>CA(:,:,1,2)
    AP=>CA(:,:,1,3)
    AE=>CA(:,:,1,4)
    AN=>CA(:,:,1,5)

    BETA=OMEGA
    F=zero
    G=zero
    H=zero
    
    K=1
    DO J=1, nyc
    JJ=J+1
    DO I=1, nxc
    II=I+1
    IF (IBLANK(I,J,K)==1) CYCLE
        B(II,JJ)=AS(I,J)/(1.D0-BETA*F(II,JJ-1)*F(II+1,JJ-1))      ! b 
        C(II,JJ)=-B(II,JJ)*F(II,JJ-1)                             ! c
        D(II,JJ)=(AW(I,J)-B(II,JJ)*G(II,JJ-1))/(1+2.D0*BETA*G(II-1,JJ))
        PHI1=C(II,JJ)*F(II+1,JJ-1)
        PHI4=D(II,JJ)*G(II-1,JJ)
        E(II,JJ)=1.D0/(AP(I,J)-B(II,JJ)*H(II,JJ-1)-C(II,JJ)*G(II+1,JJ-1)-D(II,JJ)*F(II-1,JJ)+2.D0*BETA*(PHI1+PHI4))            ! e
        F(II,JJ)=(AE(I,J)-C(II,JJ)*H(II+1,JJ-1)-2.D0*BETA*PHI1)*E(II,JJ)            ! f
        G(II,JJ)=-D(II,JJ)*H(II-1,JJ)*E(II,JJ)  
        H(II,JJ)=(AN(I,J)-BETA*PHI4)*E(II,JJ)  
    END DO
    END DO

END SUBROUTINE COEFF_LU_MSIP_2D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE COEFF_LU_MSIP_3D()
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: A, B, C, D, E, F, G, H, P, R, S, U, V
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER :: I,J,K, II,JJ,KK
    REAL(KIND=CGREAL) :: PHI1, PHI2, PHI3, PHI4, BETA
    REAL(KIND=CGREAL) :: PHI5, PHI6, PHI7, PHI8
    REAL(KIND=CGREAL) :: PHI9, PHI10, PHI11, PHI12

    B=>LU(:,:,:,1)
    C=>LU(:,:,:,2)
    D=>LU(:,:,:,3)
    E=>LU(:,:,:,4)
    F=>LU(:,:,:,5)
    G=>LU(:,:,:,6)
    H=>LU(:,:,:,7)
    A=>LU(:,:,:,8)
    P=>LU(:,:,:,9)
    R=>LU(:,:,:,10)
    S=>LU(:,:,:,11)
    U=>LU(:,:,:,12)
    V=>LU(:,:,:,13)
    
    AS=>CA(:,:,:,1)
    AW=>CA(:,:,:,2)
    AP=>CA(:,:,:,3)
    AE=>CA(:,:,:,4)
    AN=>CA(:,:,:,5)
    AF=>CA(:,:,:,6)
    AB=>CA(:,:,:,7)

    BETA=OMEGA

    H=zero
    P=zero
    R=zero
    S=zero
    U=zero
    V=zero

    DO KK=1, nzc
    K=KK+1
    DO JJ=1, nyc
    J=JJ+1
    DO II=1, nxc
    I=II+1
    IF (iblank(II,JJ,KK)==1) CYCLE

        A(I,J,K)=AB(II,JJ,KK)/(1+BETA*(P(I,J,K-1)-H(I,J,K-1)*(H(I+1,J,K-1)+R(I+1,J,K-1))    &
              -(R(I,J,K-1)-P(I+1,J,K-1)*H(I,J,K-1))*(H(I,J+1,K-1)+P(I,J+1,K-1)+R(I,J+1,K-1))))
        B(I,J,K)=-A(I,J,K)*H(I,J,K-1)
        C(I,J,K)=-A(I,J,K)*R(I,J,K-1)-B(I,J,K)*P(I+1,J,K-1)
        D(I,J,K)=(AS(II,JJ,KK)-A(I,J,K)*S(I,J,K-1)    &
                    +BETA*( (H(I+1,J-1,K)+2*S(I+1,J-1,K)+V(I+1,J-1,K))*B(I,J,K)*S(I+1,J,K-1) &
                           -S(I-1,J,K)*(AW(II,JJ,KK)-A(I,J,K)*U(I,J,K-1)))) / &
                  (1+BETA*(2*S(I,J-1,K)+U(I,J-1,K)-S(I-1,J,K)*P(I,J-1,K)-H(I,J-1,K)*(H(I+1,J-1,K)+2*S(I+1,J-1,K)+V(I+1,J-1,K)))) 
                   
        E(I,J,K)=-B(I,J,K)*S(I+1,J,K-1)-D(I,J,K)*H(I,J-1,K)
        F(I,J,K)=(AW(II,JJ,KK)-A(I,J,K)*U(I,J,K-1)-D(I,J,K)*P(I,J-1,K)  &
                  -BETA*(A(I,J,K)*P(I,J,K-1)+C(I,J,K)*P(I,J+1,K-1)+D(I,J,K)*U(I,J-1,K))) /  &
                  (1+BETA*(2*P(I-1,J,K)+S(I-1,J,K)+2*U(I-1,J,K)))

        PHI1=B(I,J,K)*H(I+1,J,K-1)
        PHI2=A(I,J,K)*P(I,J,K-1)
        PHI3=B(I,J,K)*R(I+1,J,K-1)+C(I,J,K)*H(I,J+1,K-1)
        PHI4=C(I,J,K)*P(I,J,K-1)
        PHI5=C(I,J,K)*R(I,J+1,K-1)
        PHI6=E(I,J,K)*H(I+1,J-1,K)
        PHI7=F(I,J,K)*P(I-1,J,K)
        PHI8=D(I,J,K)*S(I,J-1,K)
        PHI9=E(I,J,K)*S(I+1,J-1,K)
        PHI10=D(I,J,K)*U(I,J-1,K)+F(I,J,K)*S(I-1,J,K)
        PHI11=E(I,J,K)*V(I+1,J-1,K)
        PHI12=F(I,J,K)*U(I-1,J,K)

        G(I,J,K)=1.D0/(AP(II,JJ,KK)-A(I,J,K)*V(I,J,K-1)-B(I,J,K)*U(I+1,J,K-1)-C(I,J,K)*S(I,J+1,K-1) &
                 -D(I,J,K)*R(I,J-1,K)-E(I,J,K)*P(I+1,J-1,K)-F(I,J,K)*H(I-1,J,K) &
                 +BETA*(2*(PHI1+PHI2+PHI3)+3*PHI4+2*(PHI5+PHI6+PHI7+PHI8)+3*PHI9+2*(PHI10+PHI11+PHI12)))
        H(I,J,K)=(AE(II,JJ,KK)-B(I,J,K)*V(I+1,J,K-1)-E(I,J,K)*R(I+1,J-1,K)  &
                -BETA*(2*PHI1+PHI3+2*PHI6+PHI9+PHI11))*G(I,J,K)
        P(I,J,K)=(-C(I,J,K)*U(I,J+1,K-1)-F(I,J,K)*R(I-1,J,K))*G(I,J,K)
        R(I,J,K)=(AN(II,JJ,KK)-C(I,J,K)*V(I,J+1,K-1)-BETA*(PHI2+PHI3+2*PHI4+2*PHI5+PHI7))*G(I,J,K)
        S(I,J,K)=(-D(I,J,K)*V(I,J-1,K)-E(I,J,K)*U(I+1,J-1,K))*G(I,J,K)
        U(I,J,K)=-F(I,J,K)*V(I-1,J,K)*G(I,J,K)
        V(I,J,K)=(AF(II,JJ,KK)-BETA*(PHI8+PHI9+PHI10+PHI11+PHI12))*G(I,J,K)
        
    END DO ! II
    END DO ! JJ
    END DO ! KK
    
END SUBROUTINE COEFF_LU_MSIP_3D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE itsolv_SIP_2D(var,r)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE pressure_arrays

    IMPLICIT NONE

!... parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

!... Local variables

    INTEGER :: i,j,K,II,JJ
    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: LW, LS, LP, UN, UE
    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: AW, AS, AP, AN, AE
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1) :: RES
    REAL(KIND=CGREAL) :: RESN 
    
    LS=>LU(:,:,1,1)
    LW=>LU(:,:,1,2)
    LP=>LU(:,:,1,3)
    UE=>LU(:,:,1,4)
    UN=>LU(:,:,1,5)
    
    AS=>CA(:,:,1,1)
    AW=>CA(:,:,1,2)
    AP=>CA(:,:,1,3)
    AE=>CA(:,:,1,4)
    AN=>CA(:,:,1,5)


    RES=zero
    K=1
    
    RESN=zero
    DO J=1,nyc
    JJ=J+1
    DO I=1,nxc
    II=I+1
    IF (IBLANK(I,J,K)==1) THEN
      RES(I,J)=zero
    ELSE
      RES(I,J)=R(I,J,K)+nlw(I,J,K)-(AP(I,J)*VAR(I,J,K)+AN(I,J)*VAR(I,J+1,K)+AS(I,J)*VAR(I,J-1,K)+AE(I,J)*VAR(I+1,J,K)+AW(I,J)*VAR(I-1,J,K))
!      RESN=RESN+ABS(RES(I,J))
      RESN=MAX(RESN,ABS(RES(I,J)))
      RES(I,J)=(RES(I,J)-LS(II,JJ)*RES(I,J-1)-LW(II,JJ)*RES(I-1,J))*LP(II,JJ)
    ENDIF
    END DO
    END DO

    IF(RESN .LT. restol_Poisson) RETURN

    DO J=nyc,1,-1
    JJ=J+1
    DO I=nxc,1,-1
    II=I+1
    IF (IBLANK(I,J,K)==1) THEN
      RES(I,J)=zero
      IF (boundary_formulation == GCM_METHOD .AND. ghostcellMark(i,j,k)==1) THEN
      ELSE
        VAR(I,J,K)=zero
      ENDIF
    ELSE
      RES(I,J)=RES(I,J)-UN(II,JJ)*RES(I,J+1)-UE(II,JJ)*RES(I+1,J)
      VAR(I,J,K)=VAR(I,J,K)+RES(I,J)
    ENDIF
    END DO
    END DO
    
    RESN=zero
    DO J=nyc,1,-1
    JJ=J+1
    DO I=1,nxc
    II=I+1
    IF (IBLANK(I,J,K)==1) THEN
      RES(I,J)=zero
    ELSE
      RES(I,J)=R(I,J,K)+nlw(I,J,K)-(AP(I,J)*VAR(I,J,K)+AN(I,J)*VAR(I,J+1,K)+AS(I,J)*VAR(I,J-1,K)+AE(I,J)*VAR(I+1,J,K)+AW(I,J)*VAR(I-1,J,K))
!      RESN=RESN+ABS(RES(I,J))
      RESN=MAX(RESN,ABS(RES(I,J)))
      RES(I,J)=(RES(I,J)-LS(II,JJ)*RES(I,J+1)-LW(II,JJ)*RES(I-1,J))*LP(II,JJ)
    ENDIF
    END DO
    END DO
    
    IF(RESN .LT. restol_Poisson) RETURN

    DO J=1,nyc
    JJ=J+1
    DO I=nxc,1,-1
    II=I+1
    IF (IBLANK(I,J,K)==1) THEN
      RES(I,J)=zero
      IF (boundary_formulation == GCM_METHOD .AND. ghostcellMark(i,j,k)==1) THEN
      ELSE
        VAR(I,J,K)=zero
      ENDIF
    ELSE
      RES(I,J)=RES(I,J)-UN(II,JJ)*RES(I,J-1)-UE(II,JJ)*RES(I+1,J)
      VAR(I,J,K)=VAR(I,J,K)+RES(I,J)
    ENDIF
    END DO
    END DO

    DO I=K+1, nzc
    VAR(:,:,I)=VAR(:,:,K)
    END DO

END SUBROUTINE itsolv_SIP_2D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE itsolv_SIP_3D(var,r)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE pressure_arrays

    IMPLICIT NONE

!... parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::r

!... Local variables

END SUBROUTINE itsolv_SIP_3D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE itsolv_MSIP_2D(var,RHS,restol)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

!... parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::RHS
    REAL(KIND=CGREAL) :: restol
    
!... Local variables

    INTEGER :: i,j,K,II,JJ
    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: B, C, D, E, F, G, H
    REAL(KIND=CGREAL), DIMENSION(:,:), POINTER :: AW, AS, AP, AN, AE
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1) :: RES
    REAL(KIND=CGREAL) :: RESN 
    
    B=>LU(:,:,1,1)
    C=>LU(:,:,1,2)
    D=>LU(:,:,1,3)
    E=>LU(:,:,1,4)
    F=>LU(:,:,1,5)
    G=>LU(:,:,1,6)
    H=>LU(:,:,1,7)
    
    AS=>CA(:,:,1,1)
    AW=>CA(:,:,1,2)
    AP=>CA(:,:,1,3)
    AE=>CA(:,:,1,4)
    AN=>CA(:,:,1,5)

    RES=zero
    RESN=zero
    
    K=1
    DO J=1,nyc
    JJ=J+1
    DO I=1,nxc
    II=I+1
    IF (IBLANK(I,J,K)==1) THEN
      RES(I,J)=zero
    ELSE
      RES(I,J)=RHS(I,J,K)-(AP(I,J)*VAR(I,J,K)+AN(I,J)*VAR(I,J+1,K)+AS(I,J)*VAR(I,J-1,K)+AE(I,J)*VAR(I+1,J,K)+AW(I,J)*VAR(I-1,J,K))
!      RESN=RESN+ABS(RES(I,J))
      RESN=MAX(RESN,ABS(RES(I,J)))
      RES(I,J)=(RES(I,J)-B(II,JJ)*RES(I,J-1)-C(II,JJ)*RES(I+1,J-1)-D(II,JJ)*RES(I-1,J))*E(II,JJ)
    ENDIF
    END DO
    END DO

    IF(RESN .LT. restol) RETURN

    DO J=nyc,1,-1
    JJ=J+1
    DO I=nxc,1,-1
    II=I+1
    IF (IBLANK(I,J,K)==1) THEN
      RES(I,J)=zero
      VAR(I,J,K)=zero
    ELSE
      RES(I,J)=RES(I,J)-F(II,JJ)*RES(I+1,J)-G(II,JJ)*RES(I-1,J+1)-H(II,JJ)*RES(I,J+1)
      VAR(I,J,K)=VAR(I,J,K)+RES(I,J)
    ENDIF
    END DO
    END DO
    
    DO I=K+1, nzc
    VAR(:,:,I)=VAR(:,:,K)
    END DO

END SUBROUTINE itsolv_MSIP_2D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE itsolv_MSIP_3D(var,RHS,restol)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

!... parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)     ::RHS
    REAL(KIND=CGREAL) :: restol

!... Local variables

    INTEGER :: i,j,k, II, JJ, KK
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: A, B, C, D, E, F, G, H, P, R, S, U, V
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1) :: RES
    REAL(KIND=CGREAL) :: RESN 

    B=>LU(:,:,:,1)
    C=>LU(:,:,:,2)
    D=>LU(:,:,:,3)
    E=>LU(:,:,:,4)
    F=>LU(:,:,:,5)
    G=>LU(:,:,:,6)
    H=>LU(:,:,:,7)
    A=>LU(:,:,:,8)
    P=>LU(:,:,:,9)
    R=>LU(:,:,:,10)
    S=>LU(:,:,:,11)
    U=>LU(:,:,:,12)
    V=>LU(:,:,:,13)
    
    AS=>CA(:,:,:,1)
    AW=>CA(:,:,:,2)
    AP=>CA(:,:,:,3)
    AE=>CA(:,:,:,4)
    AN=>CA(:,:,:,5)
    AF=>CA(:,:,:,6)
    AB=>CA(:,:,:,7)

    RES=zero
    RESN=zero
    
    DO K=1,nzc
    KK=K+1
    DO J=1,nyc
    JJ=J+1
    DO I=1,nxc
    II=I+1
      IF (iblank(I,J,K)==1) THEN
        RES(I,J,K)=zero
      ELSE
        RES(I,J,K)=RHS(I,J,K)-(AP(I,J,K)*VAR(I,J,K)+  &
          AE(I,J,K)*VAR(I+1,J,K)+AW(I,J,K)*VAR(I-1,J,K)+  &
          AN(I,J,K)*VAR(I,J+1,K)+AS(I,J,K)*VAR(I,J-1,K)+  &
          AF(I,J,K)*VAR(I,J,K+1)+AB(I,J,K)*VAR(I,J,K-1))
!        RESN=RESN+ABS(RES(I,J,K))
        RESN=MAX(RESN,ABS(RES(I,J,K)))
        RES(I,J,K)=(RES(I,J,K)-A(II,JJ,KK)*RES(I,J,K-1)-B(II,JJ,KK)*RES(I+1,J,K-1)-C(II,JJ,KK)*RES(I,J+1,K-1)-  &
                               D(II,JJ,KK)*RES(I,J-1,K)-E(II,JJ,KK)*RES(I+1,J-1,K)-F(II,JJ,KK)*RES(I-1,J,K))*G(II,JJ,KK)
      ENDIF
    END DO
    END DO
    END DO

    IF(RESN .LT. restol) RETURN
    
    DO K=nzc,1,-1
    KK=K+1
    DO J=nyc,1,-1
    JJ=J+1
    DO I=nxc,1,-1
    II=I+1
      IF (iblank(I,J,K)==1) THEN
        VAR(I,J,K)=zero
        RES(I,J,K)=zero
      ELSE
        RES(I,J,K)=RES(I,J,K)-H(II,JJ,KK)*RES(I+1,J,K  )-P(II,JJ,KK)*RES(I-1,J+1,K)-R(II,JJ,KK)*RES(I,J+1,K)-   &
                              S(II,JJ,KK)*RES(I,J-1,K+1)-U(II,JJ,KK)*RES(I-1,J,K+1)-V(II,JJ,KK)*RES(I,J,K+1)
        VAR(I,J,K)=VAR(I,J,K)+RES(I,J,K)
      ENDIF
    END DO
    END DO
    END DO

END SUBROUTINE itsolv_MSIP_3D
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE calc_residual_MSIP(var,r,resm,loc)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)  ::var,r
    REAL(KIND=CGREAL),                                   INTENT (OUT) ::resm
    INTEGER,           DIMENSION(3),                     INTENT (OUT) ::loc

    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER              :: i,j,k, MZ
    INTEGER              :: iG,jG
    REAL(KIND=CGREAL)    :: res
    REAL(KIND=CGREAL) :: bmx,bpx,bcx,bc
    REAL(KIND=CGREAL) :: bmy,bpy,bcy
    REAL(KIND=CGREAL) :: bmz,bpz,bcz

!******************************************************************************
    loc = 0
    
    AS=>CA(:,:,:,1)
    AW=>CA(:,:,:,2)
    AP=>CA(:,:,:,3)
    AE=>CA(:,:,:,4)
    AN=>CA(:,:,:,5)
    IF (ndim == DIM_3D) THEN
    AF=>CA(:,:,:,6)
    AB=>CA(:,:,:,7)
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
        IF (iblank(i,j,k)==1) cycle
        bmx =   aw(i,j,K)
        bpx =   ae(i,j,K)

        bmy =   as(i,j,K)
        bpy =   an(i,j,K)

        IF (ndim == DIM_3D) THEN
        bmz =   AB(I,J,K)
        bpz =   AF(I,J,K)
        ELSE
        bmz=zero
        bpz=zero
        END IF
        bc =   AP(I,J,K)

       res    = r(i,j,k) - var(i,j,k)*bc  &
                         - var(i-1,j,k)*bmx                &
                         - var(i+1,j,k)*bpx                &
                         - var(i,j-1,k)*bmy                &
                         - var(i,j+1,k)*bpy                &
                         - var(i,j,k-1)*bmz                &
                         - var(i,j,k+1)*bpz                 
                               
       IF (ABS(res) > resm ) THEN
         resm = ABS(res)
         loc(1) = i
         loc(2) = j
         loc(3) = k
       ENDIF 

    ENDDO
    ENDDO
    ENDDO
END SUBROUTINE calc_residual_MSIP
!!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  SUBROUTINE MG_PREPARE_MSIP(MG_Precondition_MSIP_)
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mg_parameters
    USE mg_arrays    

    IMPLICIT NONE
    EXTERNAL MG_Precondition_MSIP_
    
    INTEGER :: I, ILEVEL, iErr
!    
    IF (Full_Coarsening) THEN

    DO i=1, mgLevels_X
        IF (ndim == DIM_3D) THEN
            ALLOCATE( MGX(i)%CA(7,mgrid_I(i)+1,mgrid_I(i)+1,nz),STAT=ierr )
            ALLOCATE( MGX(i)%LU(0:mgrid_I(i)+1,0:mgrid_I(i)+1,0:nz,13),STAT=ierr )
        ELSE
            ALLOCATE( MGX(i)%CA(5,mgrid_I(i)+1,mgrid_I(i)+1,1),STAT=ierr )
            ALLOCATE( MGX(i)%LU(0:mgrid_I(i)+1,0:mgrid_I(i)+1,1,7),STAT=ierr )
        ENDIF
    END DO

    ELSE

    DO i=1, mgLevels_X
        IF (ndim == DIM_3D) THEN
            ALLOCATE( MGX(i)%CA(7,mgrid_I(i)+1,ny,nz),STAT=ierr )
            ALLOCATE( MGX(i)%LU(0:mgrid_I(i)+1,0:ny,0:nz,13),STAT=ierr )
        ELSE
            ALLOCATE( MGX(i)%CA(5,mgrid_I(i)+1,ny,1),STAT=ierr )
            ALLOCATE( MGX(i)%LU(0:mgrid_I(i)+1,0:ny,1,7),STAT=ierr )
        ENDIF
    END DO
        
    DO i=2, mgLevels_Y
        IF (ndim == DIM_3D) THEN
            ALLOCATE( MGY(i)%CA(7,nx,mgrid_J(i)+1,nz),STAT=ierr )
            ALLOCATE( MGY(i)%LU(0:nx,0:mgrid_J(i)+1,0:nz,13),STAT=ierr )
        ELSE
            ALLOCATE( MGY(i)%CA(5,nx,mgrid_J(i)+1,1),STAT=ierr )
            ALLOCATE( MGY(i)%LU(0:nx,0:mgrid_J(i)+1,1,7),STAT=ierr )
        ENDIF
    END DO

    IF (ndim == DIM_3D) THEN
        DO i=2, mgLevels_Z
            ALLOCATE( MGZ(i)%CA(7,nx,NY,mgrid_K(i)+1),STAT=ierr )
            ALLOCATE( MGZ(i)%LU(0:nx,0:NY,0:mgrid_K(i)+1,13),STAT=ierr )
        END DO
    ENDIF
    
    ENDIF
    
    CALL MG_Allocate_Memory_IUP

    IF (Full_Coarsening) THEN
      DO ILEVEL=1, mgLevels_X
     
        IF (ILEVEL>1) THEN
          iblank_MG => MGX(ILEVEL)%iblank
          CALL MG_Prepare( ILEVEL, mgrid_I(iLEVEL), mgrid_I(iLEVEL), NZC)
            
        ELSE
          iblank_MG => iblank
          CALL MG_Prepare_BC(ILEVEL)
            
        ENDIF
           
        CALL MG_Precondition_MSIP_(ILEVEL, 1, 1, mgrid_I(iLEVEL), mgrid_I(iLEVEL), NZC)

      END DO ! ILEVEL
        
    ELSE
!
    DO ILEVEL=1, mgLevels_X
 
      IF (ILEVEL>1) THEN
        iblank_MG => MGX(ILEVEL)%iblank
        CALL MG_Prepare( ILEVEL, mgrid_I(iLEVEL), NYC, NZC)
        
      ELSE
        iblank_MG => iblank
        CALL MG_Prepare_BC(ILEVEL)
        
      ENDIF
       
      CALL MG_Precondition_MSIP_(ILEVEL, 1, 1, mgrid_I(iLEVEL), nyc, NZC)

    END DO ! ILEVEL
    
!    
    DO ILEVEL=2, mgLevels_Y
    
      iblank_MG => MGY(ILEVEL)%iblank
        
      CALL MG_Prepare( ILEVEL, NX-1, mgrid_J(iLEVEL), NZC)
       
      CALL MG_Precondition_MSIP_(1, ILEVEL, 1, NXC, mgrid_J(iLEVEL), NZC)

    END DO ! ILEVEL
    
!    
    IF (ndim == DIM_3D) THEN
    
      DO ILEVEL=2, mgLevels_Z
    
        iblank_MG => MGZ(ILEVEL)%iblank

        CALL MG_Prepare( ILEVEL, NXC, NYC, mgrid_K(iLEVEL))
       
        CALL MG_Precondition_MSIP_(1, 1, ILEVEL, NXC, NYC, mgrid_K(iLEVEL))

      END DO ! ILEVEL
    
    ENDIF
    ENDIF
    
    DEALLOCATE(ium_MG, iup_MG, jum_MG, jup_MG, kum_MG, kup_MG)
    
END SUBROUTINE MG_PREPARE_MSIP
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE MG_Precondition_MSIP_2D(nLevX, nLevY, nLevZ, MX, MY, MZ)
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mg_parameters
    USE mg_arrays    

    IMPLICIT NONE

    INTEGER :: nLevX, nLevY, nLevZ, MX, MY, MZ
    
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: B, C, D, E, F, G, H
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE

    INTEGER :: I,J,K, II,JJ,KK
    REAL(KIND=CGREAL) :: PHI1, PHI4, BETA

    REAL(KIND=CGREAL) :: AAe, AAw, AAn, AAs

    IF (nLevY==1) THEN
        B=>MGX(nLevX)%LU(:,:,:,1)
        C=>MGX(nLevX)%LU(:,:,:,2)
        D=>MGX(nLevX)%LU(:,:,:,3)
        E=>MGX(nLevX)%LU(:,:,:,4)
        F=>MGX(nLevX)%LU(:,:,:,5)
        G=>MGX(nLevX)%LU(:,:,:,6)
        H=>MGX(nLevX)%LU(:,:,:,7)
        
        AS=>MGX(nLevX)%CA(1,:,:,:)
        AW=>MGX(nLevX)%CA(2,:,:,:)
        AP=>MGX(nLevX)%CA(3,:,:,:)
        AE=>MGX(nLevX)%CA(4,:,:,:)
        AN=>MGX(nLevX)%CA(5,:,:,:)

    ELSE
        B=>MGY(nLevY)%LU(:,:,:,1)
        C=>MGY(nLevY)%LU(:,:,:,2)
        D=>MGY(nLevY)%LU(:,:,:,3)
        E=>MGY(nLevY)%LU(:,:,:,4)
        F=>MGY(nLevY)%LU(:,:,:,5)
        G=>MGY(nLevY)%LU(:,:,:,6)
        H=>MGY(nLevY)%LU(:,:,:,7)
        
        AS=>MGY(nLevY)%CA(1,:,:,:)
        AW=>MGY(nLevY)%CA(2,:,:,:)
        AP=>MGY(nLevY)%CA(3,:,:,:)
        AE=>MGY(nLevY)%CA(4,:,:,:)
        AN=>MGY(nLevY)%CA(5,:,:,:)

    ENDIF

    K=1
    ! Naumann B.C.on ANY boundary
    DO J=1, MY
    DO I=1, MX
	  IF (iblank_MG(I,J,K)==1) THEN
		  ae(i,j,K)=zero
		  aw(i,j,K)=zero
		  an(i,j,K)=zero
		  as(i,j,K)=zero
		  ap(i,j,K)=oned
	  ELSE
        IF (cure_pressure_oscillations .AND. iblank_MG(i,j,k)==0 .AND. ivc(i,j,k)>0 .AND. nLevX==1 .AND. nLevY==1 .AND. nLevZ==1) THEN
          AAe=face(cell(ivc(i,j,k))%F_ip)%a
          AAw=face(cell(ivc(i,j,k))%F_im)%a
          AAn=face(cell(ivc(i,j,k))%F_jp)%a
          AAs=face(cell(ivc(i,j,k))%F_jm)%a

          ae(i,j,K)=dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*AAe
          aw(i,j,K)=dxcinv_mg(i,  nLevX)*dxinv_mg(i,nLevX)*AAw
          an(i,j,K)=dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*AAn
          as(i,j,K)=dycinv_mg(j,  nLevY)*dyinv_mg(j,nLevY)*AAs
          
        ELSE
	  IF (Full_Coarsening) THEN
		ae(i,j,K)=dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)
		aw(i,j,K)=dxcinv_mg(i,  nLevX)*dxinv_mg(i,nLevX)
		an(i,j,K)=dycinv_mg(j+1,nLevX)*dyinv_mg(j,nLevX)
		as(i,j,K)=dycinv_mg(j,  nLevX)*dyinv_mg(j,nLevX)
                
	  ELSE
		ae(i,j,K)=dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)
		aw(i,j,K)=dxcinv_mg(i,  nLevX)*dxinv_mg(i,nLevX)
		an(i,j,K)=dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)
		as(i,j,K)=dycinv_mg(j,  nLevY)*dyinv_mg(j,nLevY)
                
	  ENDIF

        ENDIF
        ap(i,j,K)=-( ae(i,j,K) + aw(i,j,K) + an(i,j,K) +  as(i,j,K))

        ap(i,j,K)=ap(i,j,K) + ae(i,j,K)*iup_MG(i,j,k) + aw(i,j,K)*ium_MG(i,j,k) &
                            + an(i,j,K)*jup_MG(i,j,k) + as(i,j,K)*jum_MG(i,j,k)

        ae(i,j,K)=ae(i,j,K)*(1 - babs(iup_MG(i,j,k) ))
        aw(i,j,K)=aw(i,j,K)*(1 - babs(ium_MG(i,j,k) ))
        an(i,j,K)=an(i,j,K)*(1 - babs(jup_MG(i,j,k) ))
        as(i,j,K)=as(i,j,K)*(1 - babs(jum_MG(i,j,k) ))
          
      ENDIF
    END DO
    END DO
    
    BETA=OMEGA
    F=zero
    G=zero
    H=zero
    
    KK=1
    DO J=1, MY
    JJ=J+1
    DO I=1, MX
    II=I+1
    IF (iblank_MG(I,J,K)==1) CYCLE

        B(II,JJ,KK)=AS(I,J,K)/(1.D0-BETA*F(II,JJ-1,KK)*F(II+1,JJ-1,KK))      ! b 
        C(II,JJ,KK)=-B(II,JJ,KK)*F(II,JJ-1,KK)                             ! c
        D(II,JJ,KK)=(AW(I,J,K)-B(II,JJ,KK)*G(II,JJ-1,KK))/(1+2.D0*BETA*G(II-1,JJ,KK))

        PHI1=C(II,JJ,KK)*F(II+1,JJ-1,KK)
        PHI4=D(II,JJ,KK)*G(II-1,JJ,KK)

        E(II,JJ,KK)=1.D0/(AP(I,J,K)-B(II,JJ,KK)*H(II,JJ-1,KK)-C(II,JJ,KK)*G(II+1,JJ-1,KK)- &
                                    D(II,JJ,KK)*F(II-1,JJ,KK)+2.D0*BETA*(PHI1+PHI4))            ! e
        F(II,JJ,KK)=(AE(I,J,K)-C(II,JJ,KK)*H(II+1,JJ-1,KK)-2.D0*BETA*PHI1)*E(II,JJ,KK)            ! f
        G(II,JJ,KK)=-D(II,JJ,KK)*H(II-1,JJ,KK)*E(II,JJ,KK)  
        H(II,JJ,KK)=(AN(I,J,K)-BETA*PHI4)*E(II,JJ,KK)  

    END DO ! I
    END DO ! J
    
END SUBROUTINE MG_Precondition_MSIP_2D
!---------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE MG_Precondition_MSIP_3D(nLevX, nLevY, nLevZ, MX, MY, MZ)
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mg_parameters
    USE mg_arrays    

    IMPLICIT NONE

    INTEGER :: nLevX, nLevY, nLevZ, MX, MY, MZ
    
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: A, B, C, D, E, F, G, H, P, R, S, U, V
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER :: I,J,K, II,JJ,KK
    REAL(KIND=CGREAL) :: PHI1, PHI2, PHI3, PHI4, BETA
    REAL(KIND=CGREAL) :: PHI5, PHI6, PHI7, PHI8
    REAL(KIND=CGREAL) :: PHI9, PHI10, PHI11, PHI12

    REAL(KIND=CGREAL) :: AAe, AAw, AAn, AAs, AAf, AAb

    IF (nLevY==1 .AND. nLevZ==1) THEN
        B=>MGX(nLevX)%LU(:,:,:,1)
        C=>MGX(nLevX)%LU(:,:,:,2)
        D=>MGX(nLevX)%LU(:,:,:,3)
        E=>MGX(nLevX)%LU(:,:,:,4)
        F=>MGX(nLevX)%LU(:,:,:,5)
        G=>MGX(nLevX)%LU(:,:,:,6)
        H=>MGX(nLevX)%LU(:,:,:,7)
        A=>MGX(nLevX)%LU(:,:,:,8)
        P=>MGX(nLevX)%LU(:,:,:,9)
        R=>MGX(nLevX)%LU(:,:,:,10)
        S=>MGX(nLevX)%LU(:,:,:,11)
        U=>MGX(nLevX)%LU(:,:,:,12)
        V=>MGX(nLevX)%LU(:,:,:,13)
        
        AS=>MGX(nLevX)%CA(1,:,:,:)
        AW=>MGX(nLevX)%CA(2,:,:,:)
        AP=>MGX(nLevX)%CA(3,:,:,:)
        AE=>MGX(nLevX)%CA(4,:,:,:)
        AN=>MGX(nLevX)%CA(5,:,:,:)
        AF=>MGX(nLevX)%CA(6,:,:,:)
        AB=>MGX(nLevX)%CA(7,:,:,:)

    ELSE IF (NLEVX==1 .AND. nLevZ==1) THEN
        B=>MGY(nLevY)%LU(:,:,:,1)
        C=>MGY(nLevY)%LU(:,:,:,2)
        D=>MGY(nLevY)%LU(:,:,:,3)
        E=>MGY(nLevY)%LU(:,:,:,4)
        F=>MGY(nLevY)%LU(:,:,:,5)
        G=>MGY(nLevY)%LU(:,:,:,6)
        H=>MGY(nLevY)%LU(:,:,:,7)
        A=>MGY(nLevY)%LU(:,:,:,8)
        P=>MGY(nLevY)%LU(:,:,:,9)
        R=>MGY(nLevY)%LU(:,:,:,10)
        S=>MGY(nLevY)%LU(:,:,:,11)
        U=>MGY(nLevY)%LU(:,:,:,12)
        V=>MGY(nLevY)%LU(:,:,:,13)
        
        AS=>MGY(nLevY)%CA(1,:,:,:)
        AW=>MGY(nLevY)%CA(2,:,:,:)
        AP=>MGY(nLevY)%CA(3,:,:,:)
        AE=>MGY(nLevY)%CA(4,:,:,:)
        AN=>MGY(nLevY)%CA(5,:,:,:)
        AF=>MGY(nLevY)%CA(6,:,:,:)
        AB=>MGY(nLevY)%CA(7,:,:,:)
    ELSE
        B=>MGZ(nLevZ)%LU(:,:,:,1)
        C=>MGZ(nLevZ)%LU(:,:,:,2)
        D=>MGZ(nLevZ)%LU(:,:,:,3)
        E=>MGZ(nLevZ)%LU(:,:,:,4)
        F=>MGZ(nLevZ)%LU(:,:,:,5)
        G=>MGZ(nLevZ)%LU(:,:,:,6)
        H=>MGZ(nLevZ)%LU(:,:,:,7)
        A=>MGZ(nLevZ)%LU(:,:,:,8)
        P=>MGZ(nLevZ)%LU(:,:,:,9)
        R=>MGZ(nLevZ)%LU(:,:,:,10)
        S=>MGZ(nLevZ)%LU(:,:,:,11)
        U=>MGZ(nLevZ)%LU(:,:,:,12)
        V=>MGZ(nLevZ)%LU(:,:,:,13)
        
        AS=>MGZ(nLevZ)%CA(1,:,:,:)
        AW=>MGZ(nLevZ)%CA(2,:,:,:)
        AP=>MGZ(nLevZ)%CA(3,:,:,:)
        AE=>MGZ(nLevZ)%CA(4,:,:,:)
        AN=>MGZ(nLevZ)%CA(5,:,:,:)
        AF=>MGZ(nLevZ)%CA(6,:,:,:)
        AB=>MGZ(nLevZ)%CA(7,:,:,:)
    ENDIF
    
    ! Naumann B.C.on ANY boundary
    DO K=1, MZ
    DO J=1, MY
    DO I=1, MX
    IF (iblank_MG(I,J,K)==1) THEN
        ae(i,j,K)=zero
        aw(i,j,K)=zero
        an(i,j,K)=zero
        as(i,j,K)=zero
        af(i,j,K)=zero
        ab(i,j,K)=zero
        ap(i,j,K)=oned
    ELSE
      
        IF (cure_pressure_oscillations .AND. iblank_MG(i,j,k)==0 .AND. ivc(i,j,k)>0 .AND. nLevX==1 .AND. nLevY==1 .AND. nLevZ==1) THEN
          AAe=face(cell(ivc(i,j,k))%F_ip)%a
          AAw=face(cell(ivc(i,j,k))%F_im)%a
          AAn=face(cell(ivc(i,j,k))%F_jp)%a
          AAs=face(cell(ivc(i,j,k))%F_jm)%a
          AAf=face(cell(ivc(i,j,k))%F_kp)%a
          AAb=face(cell(ivc(i,j,k))%F_km)%a
!          AAg=face(cell(ivc(i,j,k))%F_slice)%a

          ae(i,j,K)=dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)*AAe
          aw(i,j,K)=dxcinv_mg(i,  nLevX)*dxinv_mg(i,nLevX)*AAw
          an(i,j,K)=dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)*AAn
          as(i,j,K)=dycinv_mg(j,  nLevY)*dyinv_mg(j,nLevY)*AAs
          af(i,j,K)=dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)*AAf
          ab(i,j,K)=dzcinv_mg(k,  nLevZ)*dzinv_mg(k,nLevZ)*AAb
          
        ELSE
          ae(i,j,K)=dxcinv_mg(i+1,nLevX)*dxinv_mg(i,nLevX)
          aw(i,j,K)=dxcinv_mg(i,  nLevX)*dxinv_mg(i,nLevX)
          an(i,j,K)=dycinv_mg(j+1,nLevY)*dyinv_mg(j,nLevY)
          as(i,j,K)=dycinv_mg(j,  nLevY)*dyinv_mg(j,nLevY)
          af(i,j,K)=dzcinv_mg(k+1,nLevZ)*dzinv_mg(k,nLevZ)
          ab(i,j,K)=dzcinv_mg(k,  nLevZ)*dzinv_mg(k,nLevZ)
          
        ENDIF
        
        ap(i,j,K)=-( ae(i,j,K) + aw(i,j,K) + an(i,j,K) +  as(i,j,K) + af(i,j,K) +  ab(i,j,K) )

        ap(i,j,K)=ap(i,j,K) + ae(i,j,K)*iup_MG(i,j,k) + aw(i,j,K)*ium_MG(i,j,k) &
                            + an(i,j,K)*jup_MG(i,j,k) + as(i,j,K)*jum_MG(i,j,k) &
                            + af(i,j,K)*kup_MG(i,j,k) + ab(i,j,K)*kum_MG(i,j,k)
                            
        ae(i,j,K)=ae(i,j,K)*(1 - babs(iup_MG(i,j,k) ))
        aw(i,j,K)=aw(i,j,K)*(1 - babs(ium_MG(i,j,k) ))
        an(i,j,K)=an(i,j,K)*(1 - babs(jup_MG(i,j,k) ))
        as(i,j,K)=as(i,j,K)*(1 - babs(jum_MG(i,j,k) ))
        af(i,j,K)=af(i,j,K)*(1 - babs(kup_MG(i,j,k) ))
        ab(i,j,K)=ab(i,j,K)*(1 - babs(kum_MG(i,j,k) ))
          
    ENDIF
    END DO
    END DO
    END DO
	
    BETA=OMEGA

    H=zero
    P=zero
    R=zero
    S=zero
    U=zero
    V=zero

    DO KK=1, MZ
    K=KK+1
    DO JJ=1, MY
    J=JJ+1
    DO II=1, MX
    I=II+1
    IF (iblank_MG(II,JJ,KK)==1) CYCLE

        A(I,J,K)=AB(II,JJ,KK)/(1+BETA*(P(I,J,K-1)-H(I,J,K-1)*(H(I+1,J,K-1)+R(I+1,J,K-1))    &
              -(R(I,J,K-1)-P(I+1,J,K-1)*H(I,J,K-1))*(H(I,J+1,K-1)+P(I,J+1,K-1)+R(I,J+1,K-1))))
        B(I,J,K)=-A(I,J,K)*H(I,J,K-1)
        C(I,J,K)=-A(I,J,K)*R(I,J,K-1)-B(I,J,K)*P(I+1,J,K-1)
        D(I,J,K)=(AS(II,JJ,KK)-A(I,J,K)*S(I,J,K-1)    &
                    +BETA*( (H(I+1,J-1,K)+2*S(I+1,J-1,K)+V(I+1,J-1,K))*B(I,J,K)*S(I+1,J,K-1) &
                           -S(I-1,J,K)*(AW(II,JJ,KK)-A(I,J,K)*U(I,J,K-1)))) / &
                  (1+BETA*(2*S(I,J-1,K)+U(I,J-1,K)-S(I-1,J,K)*P(I,J-1,K)-H(I,J-1,K)*(H(I+1,J-1,K)+2*S(I+1,J-1,K)+V(I+1,J-1,K)))) 
                   
        E(I,J,K)=-B(I,J,K)*S(I+1,J,K-1)-D(I,J,K)*H(I,J-1,K)
        F(I,J,K)=(AW(II,JJ,KK)-A(I,J,K)*U(I,J,K-1)-D(I,J,K)*P(I,J-1,K)  &
                  -BETA*(A(I,J,K)*P(I,J,K-1)+C(I,J,K)*P(I,J+1,K-1)+D(I,J,K)*U(I,J-1,K))) /  &
                  (1+BETA*(2*P(I-1,J,K)+S(I-1,J,K)+2*U(I-1,J,K)))

        PHI1=B(I,J,K)*H(I+1,J,K-1)
        PHI2=A(I,J,K)*P(I,J,K-1)
        PHI3=B(I,J,K)*R(I+1,J,K-1)+C(I,J,K)*H(I,J+1,K-1)
        PHI4=C(I,J,K)*P(I,J,K-1)
        PHI5=C(I,J,K)*R(I,J+1,K-1)
        PHI6=E(I,J,K)*H(I+1,J-1,K)
        PHI7=F(I,J,K)*P(I-1,J,K)
        PHI8=D(I,J,K)*S(I,J-1,K)
        PHI9=E(I,J,K)*S(I+1,J-1,K)
        PHI10=D(I,J,K)*U(I,J-1,K)+F(I,J,K)*S(I-1,J,K)
        PHI11=E(I,J,K)*V(I+1,J-1,K)
        PHI12=F(I,J,K)*U(I-1,J,K)

        G(I,J,K)=1.D0/(AP(II,JJ,KK)-A(I,J,K)*V(I,J,K-1)-B(I,J,K)*U(I+1,J,K-1)-C(I,J,K)*S(I,J+1,K-1) &
                 -D(I,J,K)*R(I,J-1,K)-E(I,J,K)*P(I+1,J-1,K)-F(I,J,K)*H(I-1,J,K) &
                 +BETA*(2*(PHI1+PHI2+PHI3)+3*PHI4+2*(PHI5+PHI6+PHI7+PHI8)+3*PHI9+2*(PHI10+PHI11+PHI12)))
        H(I,J,K)=(AE(II,JJ,KK)-B(I,J,K)*V(I+1,J,K-1)-E(I,J,K)*R(I+1,J-1,K)  &
                -BETA*(2*PHI1+PHI3+2*PHI6+PHI9+PHI11))*G(I,J,K)
        P(I,J,K)=(-C(I,J,K)*U(I,J+1,K-1)-F(I,J,K)*R(I-1,J,K))*G(I,J,K)
        R(I,J,K)=(AN(II,JJ,KK)-C(I,J,K)*V(I,J+1,K-1)-BETA*(PHI2+PHI3+2*PHI4+2*PHI5+PHI7))*G(I,J,K)
        S(I,J,K)=(-D(I,J,K)*V(I,J-1,K)-E(I,J,K)*U(I+1,J-1,K))*G(I,J,K)
        U(I,J,K)=-F(I,J,K)*V(I-1,J,K)*G(I,J,K)
        V(I,J,K)=(AF(II,JJ,KK)-BETA*(PHI8+PHI9+PHI10+PHI11+PHI12))*G(I,J,K)
        
    END DO ! II
    END DO ! JJ
    END DO ! KK
    
END SUBROUTINE MG_Precondition_MSIP_3D
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE MG_itsolv_MSIP_2D(var,RHS, nLevX,nLevY,nLevZ,mx,my,mz) 

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nLevX, nLevY, nLevZ, mx, my, mz
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)     :: RHS
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (INOUT)  :: var

    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: B, C, D, E, F, G, H
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE

    INTEGER :: i,j,k, II, JJ, KK
    REAL(KIND=CGREAL), DIMENSION(0:MX,0:MY) :: RES
    REAL(KIND=CGREAL) :: RESN 

    IF (nLevY==1) THEN
        B=>MGX(nLevX)%LU(:,:,:,1)
        C=>MGX(nLevX)%LU(:,:,:,2)
        D=>MGX(nLevX)%LU(:,:,:,3)
        E=>MGX(nLevX)%LU(:,:,:,4)
        F=>MGX(nLevX)%LU(:,:,:,5)
        G=>MGX(nLevX)%LU(:,:,:,6)
        H=>MGX(nLevX)%LU(:,:,:,7)
        
        AS=>MGX(nLevX)%CA(1,:,:,:)
        AW=>MGX(nLevX)%CA(2,:,:,:)
        AP=>MGX(nLevX)%CA(3,:,:,:)
        AE=>MGX(nLevX)%CA(4,:,:,:)
        AN=>MGX(nLevX)%CA(5,:,:,:)

    ELSE
        B=>MGY(nLevY)%LU(:,:,:,1)
        C=>MGY(nLevY)%LU(:,:,:,2)
        D=>MGY(nLevY)%LU(:,:,:,3)
        E=>MGY(nLevY)%LU(:,:,:,4)
        F=>MGY(nLevY)%LU(:,:,:,5)
        G=>MGY(nLevY)%LU(:,:,:,6)
        H=>MGY(nLevY)%LU(:,:,:,7)
        
        AS=>MGY(nLevY)%CA(1,:,:,:)
        AW=>MGY(nLevY)%CA(2,:,:,:)
        AP=>MGY(nLevY)%CA(3,:,:,:)
        AE=>MGY(nLevY)%CA(4,:,:,:)
        AN=>MGY(nLevY)%CA(5,:,:,:)
    
    ENDIF

    RES=zero
    RESN=zero
    
    K=1
    KK=1
    DO J=1,MY-1
    JJ=J+1
    DO I=1,MX-1
    II=I+1
      IF (iblank_MG(I,J,K)==1) THEN
        RES(I,J)=zero
      ELSE
        RES(I,J)=RHS(I,J,K)-(AP(I,J,K)*VAR(I,J,K)+	&
                             AN(I,J,K)*VAR(I,J+1,K)+&
                             AS(I,J,K)*VAR(I,J-1,K)+&
                             AE(I,J,K)*VAR(I+1,J,K)+&
                             AW(I,J,K)*VAR(I-1,J,K))
!        RESN=RESN+ABS(RES(I,J,K))
        RES(I,J)=(RES(I,J)-B(II,JJ,KK)*RES(I,J-1)-C(II,JJ,KK)*RES(I+1,J-1)-D(II,JJ,KK)*RES(I-1,J))*E(II,JJ,KK)
      ENDIF
    END DO
    END DO

!    IF(RESN .LT. restol_Poisson) RETURN
    
    DO J=MY-1,1,-1
    JJ=J+1
    DO I=MX-1,1,-1
    II=I+1
      IF (iblank_MG(I,J,K)==1) THEN
        RES(I,J)=zero
        VAR(I,J,K)=zero
      ELSE
        RES(I,J)=RES(I,J)-F(II,JJ,KK)*RES(I+1,J)-G(II,JJ,KK)*RES(I-1,J+1)-H(II,JJ,KK)*RES(I,J+1)
        VAR(I,J,K)=VAR(I,J,K)+RES(I,J)
      ENDIF
    END DO
    END DO
    
    DO I=K+1, MZ-1
    VAR(:,:,I)=VAR(:,:,K)
    END DO
    
!    CALL write_dump_debug_mg('var ',10,var,nLevX,nLevY,nLevZ,0,mx,0,my,0,mz)
!    CALL write_dump_debug_mg('res ',nLevX+10,res,nLevX,nLevY,nLevZ,0,mx,0,my,1,1)

END SUBROUTINE MG_itsolv_MSIP_2D
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE MG_itsolv_MSIP_3D(var,RHS, nLevX,nLevY,nLevZ,mx,my,mz) 

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nLevX, nLevY, nLevZ, mx, my, mz
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)     :: RHS
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (INOUT)  :: var

    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: A, B, C, D, E, F, G, H, P, R, S, U, V
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER :: i,j,k, II, JJ, KK
    REAL(KIND=CGREAL), DIMENSION(0:MX,0:MY,0:MZ) :: RES
    REAL(KIND=CGREAL) :: RESN 

    IF (nLevY==1 .AND. nLevZ==1) THEN
        B=>MGX(nLevX)%LU(:,:,:,1)
        C=>MGX(nLevX)%LU(:,:,:,2)
        D=>MGX(nLevX)%LU(:,:,:,3)
        E=>MGX(nLevX)%LU(:,:,:,4)
        F=>MGX(nLevX)%LU(:,:,:,5)
        G=>MGX(nLevX)%LU(:,:,:,6)
        H=>MGX(nLevX)%LU(:,:,:,7)
        P=>MGX(nLevX)%LU(:,:,:,9)
        A=>MGX(nLevX)%LU(:,:,:,8)
        R=>MGX(nLevX)%LU(:,:,:,10)
        S=>MGX(nLevX)%LU(:,:,:,11)
        U=>MGX(nLevX)%LU(:,:,:,12)
        V=>MGX(nLevX)%LU(:,:,:,13)
        
        AS=>MGX(nLevX)%CA(1,:,:,:)
        AW=>MGX(nLevX)%CA(2,:,:,:)
        AP=>MGX(nLevX)%CA(3,:,:,:)
        AE=>MGX(nLevX)%CA(4,:,:,:)
        AN=>MGX(nLevX)%CA(5,:,:,:)
        AF=>MGX(nLevX)%CA(6,:,:,:)
        AB=>MGX(nLevX)%CA(7,:,:,:)

    ELSE IF (nLevX==1 .AND. nLevZ==1) THEN
        B=>MGY(nLevY)%LU(:,:,:,1)
        C=>MGY(nLevY)%LU(:,:,:,2)
        D=>MGY(nLevY)%LU(:,:,:,3)
        E=>MGY(nLevY)%LU(:,:,:,4)
        F=>MGY(nLevY)%LU(:,:,:,5)
        G=>MGY(nLevY)%LU(:,:,:,6)
        H=>MGY(nLevY)%LU(:,:,:,7)
        A=>MGY(nLevY)%LU(:,:,:,8)
        P=>MGY(nLevY)%LU(:,:,:,9)
        R=>MGY(nLevY)%LU(:,:,:,10)
        S=>MGY(nLevY)%LU(:,:,:,11)
        U=>MGY(nLevY)%LU(:,:,:,12)
        V=>MGY(nLevY)%LU(:,:,:,13)
       
        AS=>MGY(nLevY)%CA(1,:,:,:)
        AW=>MGY(nLevY)%CA(2,:,:,:)
        AP=>MGY(nLevY)%CA(3,:,:,:)
        AE=>MGY(nLevY)%CA(4,:,:,:)
        AN=>MGY(nLevY)%CA(5,:,:,:)
        AF=>MGY(nLevY)%CA(6,:,:,:)
        AB=>MGY(nLevY)%CA(7,:,:,:)

    ELSE
        B=>MGZ(nLevZ)%LU(:,:,:,1)
        C=>MGZ(nLevZ)%LU(:,:,:,2)
        D=>MGZ(nLevZ)%LU(:,:,:,3)
        E=>MGZ(nLevZ)%LU(:,:,:,4)
        F=>MGZ(nLevZ)%LU(:,:,:,5)
        G=>MGZ(nLevZ)%LU(:,:,:,6)
        H=>MGZ(nLevZ)%LU(:,:,:,7)
        A=>MGZ(nLevZ)%LU(:,:,:,8)
        P=>MGZ(nLevZ)%LU(:,:,:,9)
        R=>MGZ(nLevZ)%LU(:,:,:,10)
        S=>MGZ(nLevZ)%LU(:,:,:,11)
        U=>MGZ(nLevZ)%LU(:,:,:,12)
        V=>MGZ(nLevZ)%LU(:,:,:,13)
        
        AS=>MGZ(nLevZ)%CA(1,:,:,:)
        AW=>MGZ(nLevZ)%CA(2,:,:,:)
        AP=>MGZ(nLevZ)%CA(3,:,:,:)
        AE=>MGZ(nLevZ)%CA(4,:,:,:)
        AN=>MGZ(nLevZ)%CA(5,:,:,:)
        AF=>MGZ(nLevZ)%CA(6,:,:,:)
        AB=>MGZ(nLevZ)%CA(7,:,:,:)
   
    ENDIF

    RES=zero
    RESN=zero
    
    DO K=1,MZ-1
    KK=K+1
    DO J=1,MY-1
    JJ=J+1
    DO I=1,MX-1
    II=I+1
      IF (iblank_MG(I,J,K)==1) THEN
        RES(I,J,K)=zero
      ELSE
        RES(I,J,K)=RHS(I,J,K)-(AP(I,J,K)*VAR(I,J,K)+  &
          AE(I,J,K)*VAR(I+1,J,K)+AW(I,J,K)*VAR(I-1,J,K)+  &
          AN(I,J,K)*VAR(I,J+1,K)+AS(I,J,K)*VAR(I,J-1,K)+  &
          AF(I,J,K)*VAR(I,J,K+1)+AB(I,J,K)*VAR(I,J,K-1))
!        RESN=RESN+ABS(RES(I,J,K))
        RES(I,J,K)=(RES(I,J,K)-A(II,JJ,KK)*RES(I,J,K-1)-B(II,JJ,KK)*RES(I+1,J,K-1)-C(II,JJ,KK)*RES(I,J+1,K-1)-  &
                               D(II,JJ,KK)*RES(I,J-1,K)-E(II,JJ,KK)*RES(I+1,J-1,K)-F(II,JJ,KK)*RES(I-1,J,K))*G(II,JJ,KK)
      ENDIF
    END DO
    END DO
    END DO

!    IF(RESN .LT. restol_Poisson) RETURN
    
    DO K=MZ-1,1,-1
    KK=K+1
    DO J=MY-1,1,-1
    JJ=J+1
    DO I=MX-1,1,-1
    II=I+1
      IF (iblank_MG(I,J,K)==1) THEN
        VAR(I,J,K)=zero
        RES(I,J,K)=zero
      ELSE
        RES(I,J,K)=RES(I,J,K)-H(II,JJ,KK)*RES(I+1,J,K  )-P(II,JJ,KK)*RES(I-1,J+1,K)-R(II,JJ,KK)*RES(I,J+1,K)-   &
                              S(II,JJ,KK)*RES(I,J-1,K+1)-U(II,JJ,KK)*RES(I-1,J,K+1)-V(II,JJ,KK)*RES(I,J,K+1)
        VAR(I,J,K)=VAR(I,J,K)+RES(I,J,K)
      ENDIF
    END DO
    END DO
    END DO

END SUBROUTINE MG_itsolv_MSIP_3D
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE MG_Residual_MSIP(var,R,resL,nLevX,nLevY,nLevZ,mx, my, mz) 

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nLevX, nLevY, nLevZ, mx, my, mz

    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (INOUT)  :: resL
    REAL(KIND=CGREAL), DIMENSION(0:mx,0:my,0:mz), INTENT (IN)     :: var, R
!
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER            :: i,j,k,iErr,MK
    INTEGER            :: iG,jG,iBody,iRow,n
    REAL(KIND=CGREAL)  :: res, resMax
    REAL(KIND=CGREAL)  :: bmx,bpx,bcx,bc
    REAL(KIND=CGREAL)  :: bmy,bpy,bcy
    REAL(KIND=CGREAL)  :: bmz,bpz,bcz

    IF (nLevY==1 .AND. nLevZ==1) THEN
        AS=>MGX(nLevX)%CA(1,:,:,:)
        AW=>MGX(nLevX)%CA(2,:,:,:)
        AP=>MGX(nLevX)%CA(3,:,:,:)
        AE=>MGX(nLevX)%CA(4,:,:,:)
        AN=>MGX(nLevX)%CA(5,:,:,:)
        IF (ndim == DIM_3D) THEN
            AF=>MGX(nLevX)%CA(6,:,:,:)
            AB=>MGX(nLevX)%CA(7,:,:,:)
        ENDIF
    ELSE IF (NLEVX==1 .AND. nLevZ==1) THEN
        AS=>MGY(nLevY)%CA(1,:,:,:)
        AW=>MGY(nLevY)%CA(2,:,:,:)
        AP=>MGY(nLevY)%CA(3,:,:,:)
        AE=>MGY(nLevY)%CA(4,:,:,:)
        AN=>MGY(nLevY)%CA(5,:,:,:)
        IF (ndim == DIM_3D) THEN
            AF=>MGY(nLevY)%CA(6,:,:,:)
            AB=>MGY(nLevY)%CA(7,:,:,:)
        ENDIF
    ELSE
        AS=>MGZ(nLevZ)%CA(1,:,:,:)
        AW=>MGZ(nLevZ)%CA(2,:,:,:)
        AP=>MGZ(nLevZ)%CA(3,:,:,:)
        AE=>MGZ(nLevZ)%CA(4,:,:,:)
        AN=>MGZ(nLevZ)%CA(5,:,:,:)
        AF=>MGZ(nLevZ)%CA(6,:,:,:)
        AB=>MGZ(nLevZ)%CA(7,:,:,:)
    ENDIF

    resMax = zero

    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
       CALL enforce_p_periodic(var) 
    ENDIF ! bcx1

    IF (ndim == DIM_3D) THEN
        MK=MZ-1
    ELSE
        MK=1
    ENDIF
!   ------------------------------------
    DO k = 1, mK
    DO j = 1, my-1
    DO i = 1, mx-1

        bmx =   AW(I,J,K)
        bpx =   AE(I,J,K)
!            bcx = - ( bmx + bpx ) 
     
        bmy =   AS(I,J,K)
        bpy =   AN(I,J,K)
!            bcy = - ( bmy + bpy ) 
        IF (ndim == DIM_3D) THEN
            bmz =   AB(I,J,K)
            bpz =   AF(I,J,K)
!            bcz = - ( bmz + bpz ) 
        ELSE
            BMZ=zero
            BPZ=zero
!            BCZ=zero
        ENDIF
 
        bc =   AP(I,J,K)*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )  &
               + REAL(iblank_MG(i,j,k),KIND=CGREAL)

        res    = R(i,j,k) - var(i,j,k)*bc  &
                     - var(i-1,j,k)*bmx                &
                     - var(i+1,j,k)*bpx                &
                     - var(i,j-1,k)*bmy                &
                     - var(i,j+1,k)*bpy                &
                     - var(i,j,k-1)*bmz                &
                     - var(i,j,k+1)*bpz

       resL(i,j,k)  = res*(oned-REAL(iblank_MG(i,j,k),KIND=CGREAL) )
!       if ((i==1 .and. j==1) .or. (i==2 .and. j==1) .or. (i==3 .and. j==1)) then
!       print *, resL(i,j,k),iblank_MG(i,j,k)
!      print *, bc, var(i,j,k), - var(i,j,k)*bc
!      print *, bmx, var(i-1,j,k), - var(i-1,j,k)*bmx
!      print *, bpx, var(i+1,j,k), - var(i+1,j,k)*bpx
!      print *, bmy, var(i,j-1,k), - var(i,j-1,k)*bmy
!      print *, bpy, var(i,j+1,k), - var(i,j+1,k)*bpy
!      print *, bmz, var(i,j,k-1), - var(i,j,k-1)*bmz
!      print *, bpz, var(i,j,k+1), - var(i,j,k+1)*bpz
!      print *, res, r(i,j,k), - var(i,j,k)*bc  &
!                     - var(i-1,j,k)*bmx                &
!                     - var(i+1,j,k)*bpx                &
!                     - var(i,j-1,k)*bmy                &
!                     - var(i,j+1,k)*bpy                &
!                     - var(i,j,k-1)*bmz                &
!                     - var(i,j,k+1)*bpz
!      print *
!       endif
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!    CALL write_dump_debug2d('AS  ',0,AS)
!    CALL write_dump_debug2d('AW  ',0,AW)
!    CALL write_dump_debug2d('AP  ',0,AP)
!    CALL write_dump_debug2d('AE  ',0,AE)
!    CALL write_dump_debug2d('AN  ',0,AN)
!	  CALL write_dump_debug_MG('resL',1,resL, nLevX,nLevY,nLevZ, mx-1,my-1, 0,mx,0,my,0,mz)
!    stop

    IF (infoconv == 1) THEN
      resMax = MAXVAL(ABS(resL))
        WRITE(*,100) nLevX, nLevY, nLevZ, resmax
    ENDIF ! infoconv

100   FORMAT('MG: residual check : ',1x,3I6,2x,E19.11)
END SUBROUTINE MG_Residual_MSIP
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE calc_residual_MG_MSIP(var,r,resm,loc)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,0:nz+1),  INTENT (IN)  ::var,r
    REAL(KIND=CGREAL),                                   INTENT (OUT) ::resm
    INTEGER,           DIMENSION(3),                     INTENT (OUT) ::loc

    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER              :: i,j,k, MZ
    INTEGER              :: iG,jG
    REAL(KIND=CGREAL)    :: res
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
        IF (iblank(i,j,k)==1) cycle
        bmx =   aw(i,j,K)
        bpx =   ae(i,j,K)

        bmy =   as(i,j,K)
        bpy =   an(i,j,K)

        IF (ndim == DIM_3D) THEN
        bmz =   AB(I,J,K)
        bpz =   AF(I,J,K)
        ELSE
        bmz=zero
        bpz=zero
        END IF
        bc =   AP(I,J,K)

       res    = r(i,j,k) - var(i,j,k)*bc  &
                         - var(i-1,j,k)*bmx                &
                         - var(i+1,j,k)*bpx                &
                         - var(i,j-1,k)*bmy                &
                         - var(i,j+1,k)*bpy                &
                         - var(i,j,k-1)*bmz                &
                         - var(i,j,k+1)*bpz                 
                               
       IF (ABS(res) > resm ) THEN
         resm = ABS(res)
         loc(1) = i
         loc(2) = j
         loc(3) = k
       ENDIF 

    ENDDO
    ENDDO
    ENDDO
END SUBROUTINE calc_residual_MG_MSIP
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE GCM_MG_Precondition_MSIP_2D()
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE grid_arrays
    USE mg_arrays    

    IMPLICIT NONE

    INTEGER :: MX, MY, MZ
    
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: B, C, D, E, F, G, H
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE

    INTEGER :: I,J,K, II,JJ,KK
    REAL(KIND=CGREAL) :: PHI1, PHI4, BETA

    REAL(KIND=CGREAL) :: AAe, AAw, AAn, AAs

    MX=nxc
    MY=nyc
    MZ=nzc

    B=>MGX(1)%LU(:,:,:,1)
    C=>MGX(1)%LU(:,:,:,2)
    D=>MGX(1)%LU(:,:,:,3)
    E=>MGX(1)%LU(:,:,:,4)
    F=>MGX(1)%LU(:,:,:,5)
    G=>MGX(1)%LU(:,:,:,6)
    H=>MGX(1)%LU(:,:,:,7)
        
    AS=>MGX(1)%CA(1,:,:,:)
    AW=>MGX(1)%CA(2,:,:,:)
    AP=>MGX(1)%CA(3,:,:,:)
    AE=>MGX(1)%CA(4,:,:,:)
    AN=>MGX(1)%CA(5,:,:,:)

    K=1
    ! Naumann B.C.on ANY boundary
    DO J=1, MY
    DO I=1, MX
	  IF (iblank(I,J,K)==1) THEN
		  ae(i,j,K)=zero
		  aw(i,j,K)=zero
		  an(i,j,K)=zero
		  as(i,j,K)=zero
		  ap(i,j,K)=oned
	  ELSE
      IF (cure_pressure_oscillations .AND. iblank_MG(i,j,k)==0 .AND. ivc(i,j,k)>0 ) THEN
        AAe=face(cell(ivc(i,j,k))%F_ip)%a
        AAw=face(cell(ivc(i,j,k))%F_im)%a
        AAn=face(cell(ivc(i,j,k))%F_jp)%a
        AAs=face(cell(ivc(i,j,k))%F_jm)%a

        ae(i,j,K)=dxcinv(i+1)*dxinv(i)*AAe
        aw(i,j,K)=dxcinv(i  )*dxinv(i)*AAw
        an(i,j,K)=dycinv(j+1)*dyinv(j)*AAn
        as(i,j,K)=dycinv(j  )*dyinv(j)*AAs
          
      ELSE
	      IF (Full_Coarsening) THEN
		      ae(i,j,K)=dxcinv(i+1)*dxinv(i)
		      aw(i,j,K)=dxcinv(i  )*dxinv(i)
		      an(i,j,K)=dycinv(j+1)*dyinv(j)
		      as(i,j,K)=dycinv(j  )*dyinv(j)
                
	      ELSE
		      ae(i,j,K)=dxcinv(i+1)*dxinv(i)
		      aw(i,j,K)=dxcinv(i  )*dxinv(i)
		      an(i,j,K)=dycinv(j+1)*dyinv(j)
		      as(i,j,K)=dycinv(j  )*dyinv(j)
                
	      ENDIF

      ENDIF
      ap(i,j,K)=-( ae(i,j,K) + aw(i,j,K) + an(i,j,K) +  as(i,j,K))

      ap(i,j,K)=ap(i,j,K) + ae(i,j,K)*iup(i,j,k) + aw(i,j,K)*ium(i,j,k) &
                          + an(i,j,K)*jup(i,j,k) + as(i,j,K)*jum(i,j,k)

      ae(i,j,K)=ae(i,j,K)*(1 - babs(iup(i,j,k) ))
      aw(i,j,K)=aw(i,j,K)*(1 - babs(ium(i,j,k) ))
      an(i,j,K)=an(i,j,K)*(1 - babs(jup(i,j,k) ))
      as(i,j,K)=as(i,j,K)*(1 - babs(jum(i,j,k) ))
          
    ENDIF
    END DO
    END DO
    
    BETA=OMEGA
    F=zero
    G=zero
    H=zero
    
    KK=1
    DO J=1, MY
      JJ=J+1
      DO I=1, MX
        II=I+1
        IF (iblank(I,J,K)==1) CYCLE

        B(II,JJ,KK)=AS(I,J,K)/(1.D0-BETA*F(II,JJ-1,KK)*F(II+1,JJ-1,KK))      ! b 
        C(II,JJ,KK)=-B(II,JJ,KK)*F(II,JJ-1,KK)                             ! c
        D(II,JJ,KK)=(AW(I,J,K)-B(II,JJ,KK)*G(II,JJ-1,KK))/(1+2.D0*BETA*G(II-1,JJ,KK))

        PHI1=C(II,JJ,KK)*F(II+1,JJ-1,KK)
        PHI4=D(II,JJ,KK)*G(II-1,JJ,KK)

        E(II,JJ,KK)=1.D0/(AP(I,J,K)-B(II,JJ,KK)*H(II,JJ-1,KK)-C(II,JJ,KK)*G(II+1,JJ-1,KK)- &
                                    D(II,JJ,KK)*F(II-1,JJ,KK)+2.D0*BETA*(PHI1+PHI4))            ! e
        F(II,JJ,KK)=(AE(I,J,K)-C(II,JJ,KK)*H(II+1,JJ-1,KK)-2.D0*BETA*PHI1)*E(II,JJ,KK)            ! f
        G(II,JJ,KK)=-D(II,JJ,KK)*H(II-1,JJ,KK)*E(II,JJ,KK)  
        H(II,JJ,KK)=(AN(I,J,K)-BETA*PHI4)*E(II,JJ,KK)  

      END DO ! I
    END DO ! J
    
END SUBROUTINE GCM_MG_Precondition_MSIP_2D
!---------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE GCM_MG_Precondition_MSIP_3D()
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE grid_arrays
    USE mg_arrays    

    IMPLICIT NONE

    INTEGER :: MX, MY, MZ
    
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: A, B, C, D, E, F, G, H, P, R, S, U, V
    REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: AW, AS, AP, AN, AE, AF, AB

    INTEGER :: I,J,K, II,JJ,KK
    REAL(KIND=CGREAL) :: PHI1, PHI2, PHI3, PHI4, BETA
    REAL(KIND=CGREAL) :: PHI5, PHI6, PHI7, PHI8
    REAL(KIND=CGREAL) :: PHI9, PHI10, PHI11, PHI12

    REAL(KIND=CGREAL) :: AAe, AAw, AAn, AAs, AAf, AAb

    MX=nxc
    MY=nyc
    MZ=nzc

    B=>MGX(1)%LU(:,:,:,1)
    C=>MGX(1)%LU(:,:,:,2)
    D=>MGX(1)%LU(:,:,:,3)
    E=>MGX(1)%LU(:,:,:,4)
    F=>MGX(1)%LU(:,:,:,5)
    G=>MGX(1)%LU(:,:,:,6)
    H=>MGX(1)%LU(:,:,:,7)
    A=>MGX(1)%LU(:,:,:,8)
    P=>MGX(1)%LU(:,:,:,9)
    R=>MGX(1)%LU(:,:,:,10)
    S=>MGX(1)%LU(:,:,:,11)
    U=>MGX(1)%LU(:,:,:,12)
    V=>MGX(1)%LU(:,:,:,13)
        
    AS=>MGX(1)%CA(1,:,:,:)
    AW=>MGX(1)%CA(2,:,:,:)
    AP=>MGX(1)%CA(3,:,:,:)
    AE=>MGX(1)%CA(4,:,:,:)
    AN=>MGX(1)%CA(5,:,:,:)
    AF=>MGX(1)%CA(6,:,:,:)
    AB=>MGX(1)%CA(7,:,:,:)
    
    ! Naumann B.C.on ANY boundary
    DO K=1, MZ
    DO J=1, MY
    DO I=1, MX
    IF (iblank(I,J,K)==1) THEN
        ae(i,j,K)=zero
        aw(i,j,K)=zero
        an(i,j,K)=zero
        as(i,j,K)=zero
        af(i,j,K)=zero
        ab(i,j,K)=zero
        ap(i,j,K)=oned
    ELSE
      
        IF (cure_pressure_oscillations .AND. iblank_MG(i,j,k)==0 .AND. ivc(i,j,k)>0) THEN
          AAe=face(cell(ivc(i,j,k))%F_ip)%a
          AAw=face(cell(ivc(i,j,k))%F_im)%a
          AAn=face(cell(ivc(i,j,k))%F_jp)%a
          AAs=face(cell(ivc(i,j,k))%F_jm)%a
          AAf=face(cell(ivc(i,j,k))%F_kp)%a
          AAb=face(cell(ivc(i,j,k))%F_km)%a
!          AAg=face(cell(ivc(i,j,k))%F_slice)%a

          ae(i,j,K)=dxcinv(i+1)*dxinv(i)*AAe
          aw(i,j,K)=dxcinv(i  )*dxinv(i)*AAw
          an(i,j,K)=dycinv(j+1)*dyinv(j)*AAn
          as(i,j,K)=dycinv(j  )*dyinv(j)*AAs
          af(i,j,K)=dzcinv(k+1)*dzinv(k)*AAf
          ab(i,j,K)=dzcinv(k  )*dzinv(k)*AAb
          
        ELSE
          ae(i,j,K)=dxcinv(i+1)*dxinv(i)
          aw(i,j,K)=dxcinv(i  )*dxinv(i)
          an(i,j,K)=dycinv(j+1)*dyinv(j)
          as(i,j,K)=dycinv(j  )*dyinv(j)
          af(i,j,K)=dzcinv(k+1)*dzinv(k)
          ab(i,j,K)=dzcinv(k  )*dzinv(k)
          
        ENDIF
        
        ap(i,j,K)=-( ae(i,j,K) + aw(i,j,K) + an(i,j,K) +  as(i,j,K) + af(i,j,K) +  ab(i,j,K) )

        ap(i,j,K)=ap(i,j,K) + ae(i,j,K)*iup(i,j,k) + aw(i,j,K)*ium(i,j,k) &
                            + an(i,j,K)*jup(i,j,k) + as(i,j,K)*jum(i,j,k) &
                            + af(i,j,K)*kup(i,j,k) + ab(i,j,K)*kum(i,j,k)
                            
        ae(i,j,K)=ae(i,j,K)*(1 - babs(iup(i,j,k) ))
        aw(i,j,K)=aw(i,j,K)*(1 - babs(ium(i,j,k) ))
        an(i,j,K)=an(i,j,K)*(1 - babs(jup(i,j,k) ))
        as(i,j,K)=as(i,j,K)*(1 - babs(jum(i,j,k) ))
        af(i,j,K)=af(i,j,K)*(1 - babs(kup(i,j,k) ))
        ab(i,j,K)=ab(i,j,K)*(1 - babs(kum(i,j,k) ))
          
    ENDIF
    END DO
    END DO
    END DO
	
    BETA=OMEGA

    H=zero
    P=zero
    R=zero
    S=zero
    U=zero
    V=zero

    DO KK=1, MZ
    K=KK+1
    DO JJ=1, MY
    J=JJ+1
    DO II=1, MX
    I=II+1
    IF (iblank(II,JJ,KK)==1) CYCLE

        A(I,J,K)=AB(II,JJ,KK)/(1+BETA*(P(I,J,K-1)-H(I,J,K-1)*(H(I+1,J,K-1)+R(I+1,J,K-1))    &
              -(R(I,J,K-1)-P(I+1,J,K-1)*H(I,J,K-1))*(H(I,J+1,K-1)+P(I,J+1,K-1)+R(I,J+1,K-1))))
        B(I,J,K)=-A(I,J,K)*H(I,J,K-1)
        C(I,J,K)=-A(I,J,K)*R(I,J,K-1)-B(I,J,K)*P(I+1,J,K-1)
        D(I,J,K)=(AS(II,JJ,KK)-A(I,J,K)*S(I,J,K-1)    &
                    +BETA*( (H(I+1,J-1,K)+2*S(I+1,J-1,K)+V(I+1,J-1,K))*B(I,J,K)*S(I+1,J,K-1) &
                           -S(I-1,J,K)*(AW(II,JJ,KK)-A(I,J,K)*U(I,J,K-1)))) / &
                  (1+BETA*(2*S(I,J-1,K)+U(I,J-1,K)-S(I-1,J,K)*P(I,J-1,K)-H(I,J-1,K)*(H(I+1,J-1,K)+2*S(I+1,J-1,K)+V(I+1,J-1,K)))) 
                   
        E(I,J,K)=-B(I,J,K)*S(I+1,J,K-1)-D(I,J,K)*H(I,J-1,K)
        F(I,J,K)=(AW(II,JJ,KK)-A(I,J,K)*U(I,J,K-1)-D(I,J,K)*P(I,J-1,K)  &
                  -BETA*(A(I,J,K)*P(I,J,K-1)+C(I,J,K)*P(I,J+1,K-1)+D(I,J,K)*U(I,J-1,K))) /  &
                  (1+BETA*(2*P(I-1,J,K)+S(I-1,J,K)+2*U(I-1,J,K)))

        PHI1=B(I,J,K)*H(I+1,J,K-1)
        PHI2=A(I,J,K)*P(I,J,K-1)
        PHI3=B(I,J,K)*R(I+1,J,K-1)+C(I,J,K)*H(I,J+1,K-1)
        PHI4=C(I,J,K)*P(I,J,K-1)
        PHI5=C(I,J,K)*R(I,J+1,K-1)
        PHI6=E(I,J,K)*H(I+1,J-1,K)
        PHI7=F(I,J,K)*P(I-1,J,K)
        PHI8=D(I,J,K)*S(I,J-1,K)
        PHI9=E(I,J,K)*S(I+1,J-1,K)
        PHI10=D(I,J,K)*U(I,J-1,K)+F(I,J,K)*S(I-1,J,K)
        PHI11=E(I,J,K)*V(I+1,J-1,K)
        PHI12=F(I,J,K)*U(I-1,J,K)

        G(I,J,K)=1.D0/(AP(II,JJ,KK)-A(I,J,K)*V(I,J,K-1)-B(I,J,K)*U(I+1,J,K-1)-C(I,J,K)*S(I,J+1,K-1) &
                 -D(I,J,K)*R(I,J-1,K)-E(I,J,K)*P(I+1,J-1,K)-F(I,J,K)*H(I-1,J,K) &
                 +BETA*(2*(PHI1+PHI2+PHI3)+3*PHI4+2*(PHI5+PHI6+PHI7+PHI8)+3*PHI9+2*(PHI10+PHI11+PHI12)))
        H(I,J,K)=(AE(II,JJ,KK)-B(I,J,K)*V(I+1,J,K-1)-E(I,J,K)*R(I+1,J-1,K)  &
                -BETA*(2*PHI1+PHI3+2*PHI6+PHI9+PHI11))*G(I,J,K)
        P(I,J,K)=(-C(I,J,K)*U(I,J+1,K-1)-F(I,J,K)*R(I-1,J,K))*G(I,J,K)
        R(I,J,K)=(AN(II,JJ,KK)-C(I,J,K)*V(I,J+1,K-1)-BETA*(PHI2+PHI3+2*PHI4+2*PHI5+PHI7))*G(I,J,K)
        S(I,J,K)=(-D(I,J,K)*V(I,J-1,K)-E(I,J,K)*U(I+1,J-1,K))*G(I,J,K)
        U(I,J,K)=-F(I,J,K)*V(I-1,J,K)*G(I,J,K)
        V(I,J,K)=(AF(II,JJ,KK)-BETA*(PHI8+PHI9+PHI10+PHI11+PHI12))*G(I,J,K)
        
    END DO ! II
    END DO ! JJ
    END DO ! KK
    
END SUBROUTINE GCM_MG_Precondition_MSIP_3D
!---------------------------------------------------------------------

