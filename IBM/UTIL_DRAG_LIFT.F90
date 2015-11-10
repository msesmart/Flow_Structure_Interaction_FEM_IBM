
SUBROUTINE drag_lift_solid()
!   Compute the Lift and Drag coefficients for n-Bodies in the flow
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE usr_module ,ONLY : scx,scy,scz,scmx,scmy,scmz	!VEERA
    USE body_dynamics

    IMPLICIT NONE

    INTEGER              :: i,j,k
    INTEGER              :: ibody
    INTEGER              :: m,mG

    REAL(KIND=CGREAL)    :: amx,apx,acx
    REAL(KIND=CGREAL)    :: amy,apy,acy
    REAL(KIND=CGREAL)    :: amz,apz,acz

    REAL(KIND=CGREAL)    :: cxp,cyp,czp,cxs,cys,czs
    REAL(KIND=CGREAL)    :: cmxp,cmyp,cmzp,cmxs,cmys,cmzs
    REAL(KIND=CGREAL)    :: cx,cy,cz
    REAL(KIND=CGREAL)    :: cmx,cmy,cmz
    REAL(KIND=CGREAL)    :: cpw
    REAL(KIND=CGREAL)    :: distMin,dist
    REAL(KIND=CGREAL)    :: xp,yp,zp,xs,ys,zs
    REAL(KIND=CGREAL)    :: xc_xcent_min,xc_xcent_max,yc_ycent_min,yc_ycent_max
    REAL(KIND=CGREAL)    :: xp_min,xp_max,yp_min,yp_max

    CHARACTER*9          :: dragfile
    CHARACTER*25         :: indragfile
    CHARACTER*8         :: bmfFileName ! for debug

    DO iBody = 1, nBody_solid

!---- Pressure Term

      cxp  = zero
      cyp  = zero
      czp  = zero
      cmxp = zero
      cmyp = zero
      cmzp = zero
      cpw =  zero          !Power coefficient by Meliha
      cxs  = zero
      cys  = zero
      czs  = zero
      cmxs = zero
      cmys = zero
      cmzs = zero

      xc_xcent_min = 1e+10
      xc_xcent_max = -1e-10
      yc_ycent_min = 1e+10
      yc_ycent_max = -1e-10
      xp_min = 1e+10
      xp_max = -1e-10
      yp_min = 1e+10
      yp_max = -1e-10

      DO i=1,nPtsBodyMarker(iBody)
        bodyMarkerForce(3*i-2)=0.0
        bodyMarkerForce(3*i-1)=0.0
        bodyMarkerForce(3*i)=0.0
      ENDDO

      DO k = 1, nz-1
      DO j = 1, ny-1
      DO i = 1, nx-1

        distMin = 1.0E8_CGREAL

        IF ( bodyNum(i-1,j  ,k  ) == iBody .OR. &
             bodyNum(i+1,j  ,k  ) == iBody .OR. &
             bodyNum(i  ,j+1,k  ) == iBody .OR. &
             bodyNum(i  ,j-1,k  ) == iBody .OR. &
             bodyNum(i  ,j  ,k+1) == iBody .OR. &
             bodyNum(i  ,j  ,k-1) == iBody      ) THEN
            ! The following do-loop of m is added for getting mG in order to calculate power.
            DO m = 1,nPtsBodyMarker(iBody)
               dist = (xc(i) - xBodyMarker(iBody,m))**2 &
                     +(yc(j) - yBodyMarker(iBody,m))**2 &
                     +(zc(k) - zBodyMarker(iBody,m))**2
               IF ( dist <= distMin ) THEN
                  distMin = dist
                  mG      = m
               ENDIF
            ENDDO ! m

            xp = -p(i,j,k)*ium(i,j,k)*iblank(i-1,j,k)*dy(j)*dz(k) &
	             +p(i,j,k)*iup(i,j,k)*iblank(i+1,j,k)*dy(j)*dz(k)
	        yp = -p(i,j,k)*jum(i,j,k)*iblank(i,j-1,k)*dx(i)*dz(k) &
	             +p(i,j,k)*jup(i,j,k)*iblank(i,j+1,k)*dx(i)*dz(k)
	        zp = -p(i,j,k)*kum(i,j,k)*iblank(i,j,k-1)*dx(i)*dy(j) &
	             +p(i,j,k)*kup(i,j,k)*iblank(i,j,k+1)*dx(i)*dy(j)
            cxp = cxp + xp
	        cyp = cyp + yp
	        czp = czp + zp

            bodyMarkerForce(3*mG-2)=bodyMarkerForce(3*mG-2)+xp
            bodyMarkerForce(3*mG-1)=bodyMarkerForce(3*mG-1)+yp
            bodyMarkerForce(3*mG)=bodyMarkerForce(3*mG)+zp

          IF (.NOT. Prsb_MomentRef) THEN
             cmxp = cmxp + ( yc(j) - ycent(iBody) )*zp
             cmyp = cmyp - ( xc(i) - xcent(iBody) )*zp
             cmzp = cmzp + ( xc(i) - xcent(iBody) )*yp &
                         - ( yc(j) - ycent(iBody) )*xp

             IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                cmxp = cmxp - ( zc(k) - zcent(iBody) )*yp
                cmyp = cmyp + ( zc(k) - zcent(iBody) )*xp
             ENDIF

          ELSEIF (Prsb_MomentRef) THEN
             cmxp = cmxp + ( yc(j) - Moment_refy )*zp
             cmyp = cmyp - ( xc(i) - Moment_refx )*zp
             cmzp = cmzp + ( xc(i) - Moment_refx )*yp &
                         - ( yc(j) - Moment_refy )*xp

             IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                cmxp = cmxp - ( zc(k) - Moment_refz )*yp
                cmyp = cmyp + ( zc(k) - Moment_refz )*xp
             ENDIF
          ENDIF
!----     Shear Stress Term

          amx = dxinv(i)*ium(i,j,k)*2.0_CGREAL
          apx = dxinv(i)*iup(i,j,k)*2.0_CGREAL

          amy = dyinv(j)*jum(i,j,k)*2.0_CGREAL
          apy = dyinv(j)*jup(i,j,k)*2.0_CGREAL

          amz = dzinv(k)*kum(i,j,k)*2.0_CGREAL
          apz = dzinv(k)*kup(i,j,k)*2.0_CGREAL

          xs = reinv * dx(i)* dz(k) *                              &
                    ( amy*(    u(i,j,k) -bcyu(i,j,k) )*iblank(i,j-1,k)   &
                    - apy*( bcyu(i,j,k) -   u(i,j,k) )*iblank(i,j+1,k) ) &
                    +reinv * dx(i)* dy(j) *                              &
                    ( amz*(    u(i,j,k) -bczu(i,j,k) )*iblank(i,j,k-1)   &
                    - apz*( bczu(i,j,k) -   u(i,j,k) )*iblank(i,j,k+1) )

          ys = reinv * dy(j)* dz(k) *                              &
                    ( amx*(    v(i,j,k) -bcxv(i,j,k) )*iblank(i-1,j,k)   &
                    - apx*( bcxv(i,j,k) -   v(i,j,k) )*iblank(i+1,j,k) ) &
                    +reinv * dx(i)* dy(j) *                              &
                    ( amz*(    v(i,j,k) -bczv(i,j,k) )*iblank(i,j,k-1)   &
                    - apz*( bczv(i,j,k) -   v(i,j,k) )*iblank(i,j,k+1) )

          zs = reinv * dy(j)* dz(k) *                              &
                    ( amx*(    w(i,j,k) -bcxw(i,j,k) )*iblank(i-1,j,k)   &
                    - apx*( bcxw(i,j,k) -   w(i,j,k) )*iblank(i+1,j,k) ) &
                    +reinv * dx(i)* dz(k) *                              &
                    ( amy*(    w(i,j,k) -bcyw(i,j,k) )*iblank(i,j-1,k)   &
                    - apy*( bcyw(i,j,k) -   w(i,j,k) )*iblank(i,j+1,k) )

          cxs = cxs + xs
          cys = cys + ys
          czs = czs + zs

          bodyMarkerForce(3*mG-2)=xs+bodyMarkerForce(3*mG-2)
          bodyMarkerForce(3*mG-1)=ys+bodyMarkerForce(3*mG-1)
          bodyMarkerForce(3*mG)=zs+bodyMarkerForce(3*mG)

          IF (.NOT. Prsb_MomentRef) THEN
             cmxs = cmxs + ( yc(j) - ycent(iBody) )*zs
             cmys = cmys - ( xc(i) - xcent(iBody) )*zs
             cmzs = cmzs + ( xc(i) - xcent(iBody) )*ys &
                         - ( yc(j) - ycent(iBody) )*xs

             IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                cmxs = cmxs - ( zc(k) - zcent(iBody) )*ys
                cmys = cmys + ( zc(k) - zcent(iBody) )*xs
             ENDIF
          ELSEIF (Prsb_MomentRef) THEN
             cmxs = cmxs + ( yc(j) - Moment_refy )*zs
             cmys = cmys - ( xc(i) - Moment_refx )*zs
             cmzs = cmzs + ( xc(i) - Moment_refx )*ys &
                         - ( yc(j) - Moment_refy )*xs

             IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                cmxs = cmxs - ( zc(k) - Moment_refz )*ys
                cmys = cmys + ( zc(k) - Moment_refz )*xs
             ENDIF

          ENDIF

          if (xc(i) - xcent(iBody) < xc_xcent_min) xc_xcent_min = xc(i) - xcent(iBody)
          if (xc(i) - xcent(iBody) > xc_xcent_max) xc_xcent_max = xc(i) - xcent(iBody)
          if (yc(j) - ycent(iBody) < yc_ycent_min) yc_ycent_min = yc(j) - ycent(iBody)
          if (yc(j) - ycent(iBody) > yc_ycent_max) yc_ycent_max = yc(j) - ycent(iBody)

          if (xp > xp_max) xp_max = xp
          if (xp < xp_min) xp_min = xp
          if (yp > yp_max) yp_max = yp
          if (yp < yp_min) yp_min = yp

! Power
          cpw = cpw + (xp+xs)*(uBodyMarker(ibody,mG)-uinit) + &
                      (yp+ys)*(vBodyMarker(ibody,mG)-vinit) + &
                      (zp+zs)*(wBodyMarker(ibody,mG)-winit)
        ENDIF ! bodyNum

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

      scx(ibody)  = cxp  + cxs         !VEERA ... Fed to Flow Induced Motion subroutine
      scy(ibody)  = cyp  + cys         !
      scz(ibody)  = czp  + czs         !
      scmx(ibody) = cmxp + cmxs        !
      scmy(ibody) = cmyp + cmys        !
      scmz(ibody) = cmzp + cmzs        !VEERA ... Fed to Flow Induced Motion subroutine

!--- Construct Components
      cxp = 2.0_CGREAL*cxp
      cyp = 2.0_CGREAL*cyp
      czp = 2.0_CGREAL*czp
      cxs = 2.0_CGREAL*cxs
      cys = 2.0_CGREAL*cys
      czs = 2.0_CGREAL*czs
      cpw = 2.0_CGREAL*cpw      !Added by Wanh on 04/10/11
      cmxp= 2.0_CGREAL*cmxp
      cmxs= 2.0_CGREAL*cmxs
      cmyp= 2.0_CGREAL*cmyp
      cmys= 2.0_CGREAL*cmys
      cmzp= 2.0_CGREAL*cmzp
      cmzs= 2.0_CGREAL*cmzs

      bodyMarkerForce(3*mG-2)=2.0*bodyMarkerForce(3*mG-2)
      bodyMarkerForce(3*mG-1)=2.0*bodyMarkerForce(3*mG-1)
      bodyMarkerForce(3*mG)=2.0*bodyMarkerForce(3*mG)

      IF ( ndim == DIM_2D ) THEN
        cxp = cxp/zout
        cyp = cyp/zout
        cxs = cxs/zout
        cys = cys/zout
        cpw = cpw/zout      !add by yan
        cmzp = cmzp/zout
        cmzs = cmzs/zout
      ENDIF

      cx  = cxp  + cxs         ! for non-cylinders these are
      cy  = cyp  + cys         ! raw
      cz  = czp  + czs         ! forces
      cmx = cmxp + cmxs        ! and moments
      cmy = cmyp + cmys        ! need to non-dmensionalize them
      cmz = cmzp + cmzs        ! during post-processing

!---  Write File
!     IF statement is added by Wanh
      IF (boundary_motion_type(iBody)<=PRESCRIBED) THEN
          WRITE(ifuDragOut+iBody-1,121) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
                                        ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw
      ENDIF
      IF (boundary_motion_type(iBody)>PRESCRIBED .and. Converged_FSI(iBody) ) THEN
         IF (ndim == DIM_2D) THEN
!          WRITE(ifuDragOut+iBody-1,122) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
!                                        ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw &
!                                        ,cmx/cx, cmy/cy !Added on 04/01/11 for aero center.
          WRITE(ifuDragOut+iBody-1,121) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
                                        ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw

         ELSEIF (ndim == DIM_3D) THEN
!            WRITE(ifuDragOut+iBody-1,122) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
!                                         ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw &
!                                         ,cmx/cx,cmy/cy,cmz/cz !Added on 04/01/11 for aero center.
            WRITE(ifuDragOut+iBody-1,121) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
                                         ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw
            write(ifuforce_debug+ibody-1,124) ntime,niterFS,cx,cy,cz,cmx,cmy,cmz   ! Added by Geng for debug
         ENDIF
      ENDIF
   END DO ! ibody

   print *, 'output bodyMarkerForce for debug'
   write(bmfFileName,'(a8)')ntime
   open(unit=66,file='IBM_BMF.dat',status='unknown')
!   open(unit=66,file=trim(adjustl(bmfFileName))//'.IBMBMF.dat',status='unknown')
   DO i=1,nPtsBodyMarker(1)
     write(66,*) bodyMarkerForce(i*3-2), bodyMarkerForce(i*3-1), bodyMarkerForce(i*3)
   END DO
   close(66)


121 FORMAT(20(3x,1PE12.5))
122 FORMAT(22(3x,1PE12.5))
123 FORMAT(23(3x,1PE12.5))
124 format(I5,I3,6(1x,1PE12.5))
END SUBROUTINE  drag_lift_solid
!-------------------------------------------------------------------------------

SUBROUTINE drag_lift_membrane()
!
!   Compute the Lift and Drag coefficients for n-Bodies in the flow
!
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE usr_module ,ONLY : scx,scy,scz,scmx,scmy,scmz	!VEERA
    USE body_dynamics  !, ONLY : niterFS,thickoverlength ! Wanh temp

    IMPLICIT NONE

    INTEGER              :: i,j,k
    INTEGER              :: ibody, izoneMax, MG, m, zoneNum, zoneMaxNum   !zoneNum, added by Chengyu
    INTEGER              :: K_ini, K_end

    REAL(KIND=CGREAL)    :: amx,apx,acx
    REAL(KIND=CGREAL)    :: amy,apy,acy
    REAL(KIND=CGREAL)    :: amz,apz,acz

    REAL(KIND=CGREAL)    :: cxp,cyp,czp,cxs,cys,czs
    REAL(KIND=CGREAL)    :: cmxp,cmyp,cmzp,cmxs,cmys,cmzs
    REAL(KIND=CGREAL)    :: cx,cy,cz
    REAL(KIND=CGREAL)    :: cmx,cmy,cmz, cpw, distMin, dist
    REAL(KIND=CGREAL)    :: xp,yp,zp,xs,ys,zs
    REAL(KIND=CGREAL)    :: segmentLength,segmentMass,kmid,ktotal
    REAL(KIND=CGREAL)    :: chordLength,wingMass,wingMOI,kTrans,kRot,kRigid,wingAngw,xTemp1,xTemp2,wingAlpha,wingAlphaNext
    REAL(KIND=CGREAL)    :: xc_xcent_min,xc_xcent_max,yc_ycent_min,yc_ycent_max
    REAL(KIND=CGREAL)    :: xp_min,xp_max,yp_min,yp_max
    REAL(KIND=CGREAL)    :: xmid,ymid,umid,vmid

    integer          :: itt, izone

    REAL(KIND=CGREAL)    :: pum, pup

    CHARACTER*9          :: dragfile
    CHARACTER*25         :: indragfile

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: zone_cxp,zone_cxs,zone_cx, &         !added by Chengyu
                                                      zone_cyp,zone_cys,zone_cy, &
                                                      zone_czp,zone_czs,zone_cz, &
                                                      zone_cmxp,zone_cmxs,zone_cmx, &
                                                      zone_cmyp,zone_cmys,zone_cmy, &
                                                      zone_cmzp,zone_cmzs,zone_cmz, &
                                                      zone_cpw

    !WRITE (*,*) 'Calling subroutine drag_lift_membrane()!!!!!!!!!!!!!!!!!!!!!!!! ...'
    !pause
      zoneMaxNum=maxval(zoneMax)

      ALLOCATE(zone_cxp(zoneMaxNum))
      ALLOCATE(zone_cyp(zoneMaxNum))
      ALLOCATE(zone_czp(zoneMaxNum))

      ALLOCATE(zone_cmxp(zoneMaxNum))
      ALLOCATE(zone_cmyp(zoneMaxNum))
      ALLOCATE(zone_cmzp(zoneMaxNum))
      ALLOCATE(zone_cpw(zoneMaxNum))
      ALLOCATE(zone_cxs(zoneMaxNum))
      ALLOCATE(zone_cys(zoneMaxNum))
      ALLOCATE(zone_czs(zoneMaxNum))
      ALLOCATE(zone_cmxs(zoneMaxNum))
      ALLOCATE(zone_cmys(zoneMaxNum))
      ALLOCATE(zone_cmzs(zoneMaxNum))
      ALLOCATE(zone_cx(zoneMaxNum))
      ALLOCATE(zone_cy(zoneMaxNum))
      ALLOCATE(zone_cz(zoneMaxNum))

      ALLOCATE(zone_cmx(zoneMaxNum))
      ALLOCATE(zone_cmy(zoneMaxNum))
      ALLOCATE(zone_cmz(zoneMaxNum))

    IF (ndim == DIM_2D) THEN
          K_ini = 1
          K_end = 2
    ELSE
         K_ini = 2
         K_end = NZ-2
    END IF

  DO iBody = nBody_solid+1, nBody_solid+nBody_membrane
      DO i=1,nPtsBodyMarker(iBody)
        bodyMarkerForce(3*i-2)=0.0
        bodyMarkerForce(3*i-1)=0.0
        bodyMarkerForce(3*i)=0.0
      ENDDO
!---- Pressure Term and Shear Stress Term

      cxp  = zero
      cyp  = zero
      czp  = zero
      cmxp = zero
      cmyp = zero
      cmzp = zero
      cpw =  zero
      cxs  = zero
      cys  = zero
      czs  = zero
      cmxs = zero
      cmys = zero
      cmzs = zero
      zone_cxp  = zero           !added by Chengyu
      zone_cyp  = zero
      zone_czp  = zero
      zone_cmxp = zero
      zone_cmyp = zero
      zone_cmzp = zero
      zone_cpw =  zero
      zone_cxs  = zero
      zone_cys  = zero
      zone_czs  = zero
      zone_cmxs = zero
      zone_cmys = zero
      zone_cmzs = zero

      xc_xcent_min = 1e+10
      xc_xcent_max = -1e-10
      yc_ycent_min = 1e+10
      yc_ycent_max = -1e-10

      xp_min = 1e+10
      xp_max = -1e-10
      yp_min = 1e+10
      yp_max = -1e-10

      itt = 0
      pup = 0.d0
      pum = 0.d0

      DO k = K_ini, K_end
      DO j = 2, ny-2
      DO i = 2, nx-2
        distMin = 1.0E8_CGREAL
        IF ( bodyNum(i-1,j  ,k  ) == iBody .OR. &
             bodyNum(i+1,j  ,k  ) == iBody .OR. &
             bodyNum(i  ,j+1,k  ) == iBody .OR. &
             bodyNum(i  ,j-1,k  ) == iBody .OR. &
             bodyNum(i  ,j  ,k+1) == iBody .OR. &
             bodyNum(i  ,j  ,k-1) == iBody      ) THEN

             DO m = 1,nPtsBodyMarker(iBody)
               dist = (xc(i) - xBodyMarker(iBody,m))**2 &
                     +(yc(j) - yBodyMarker(iBody,m))**2 &
                     +(zc(k) - zBodyMarker(iBody,m))**2

               IF ( dist <= distMin ) THEN
                      distMin = dist
                      mG  = m
                      IF (zoneSeparate) THEN
                      zoneNum = zoneMarker(iBody,mG)       !added by Chengyu
                      ENDIF
               ENDIF
            ENDDO ! m
            zoneCheck(i,j,k)=zoneNum

            xp = -p(i,j,k)*ium(i,j,k)*dy(j)*dz(k) &
	           +p(i,j,k)*iup(i,j,k)*dy(j)*dz(k)
	        yp = -p(i,j,k)*jum(i,j,k)*dx(i)*dz(k) &
	           +p(i,j,k)*jup(i,j,k)*dx(i)*dz(k)
            zp = -p(i,j,k)*kum(i,j,k)*dx(i)*dy(j) &
	           +p(i,j,k)*kup(i,j,k)*dx(i)*dy(j)

            cxp = cxp + xp
	        cyp = cyp + yp
	        czp = czp + zp

            bodyMarkerForce(3*mG-2)=bodyMarkerForce(3*mG-2)+xp
            bodyMarkerForce(3*mG-1)=bodyMarkerForce(3*mG-1)+yp
            bodyMarkerForce(3*mG)=bodyMarkerForce(3*mG)+zp

            IF (zoneSeparate) THEN
            zone_cxp(zoneNum) = zone_cxp(zoneNum) + xp      !added by Chengyu
	        zone_cyp(zoneNum) = zone_cyp(zoneNum) + yp
	        zone_czp(zoneNum) = zone_czp(zoneNum) + zp
            ENDIF

            IF (.NOT. Prsb_MomentRef) THEN
              IF(boundary_motion_type(ibody)==BIO_FOLLOWED_DYNAMICS_COUPLED)THEN
                 cmxp = cmxp + ( yc(j) - ycent(1) )*zp
                 cmyp = cmyp - ( xc(i) - xcent(1) )*zp
                 cmzp = cmzp + ( xc(i) - xcent(1) )*yp &
                             - ( yc(j) - ycent(1) )*xp
                IF (zoneSeparate) THEN
                 zone_cmxp(zoneNum) = zone_cmxp(zoneNum) + ( yc(j) - ycent(iBody) )*zp          !added by Chengyu
                 zone_cmyp(zoneNum) = zone_cmyp(zoneNum) - ( xc(i) - xcent(iBody) )*zp
                 zone_cmzp(zoneNum) = zone_cmzp(zoneNum) + ( xc(i) - xcent(iBody) )*yp &
                                                         - ( yc(j) - ycent(iBody) )*xp
                ENDIF

              ELSE
                 cmxp = cmxp + ( yc(j) - ycent(iBody) )*zp
                 cmyp = cmyp - ( xc(i) - xcent(iBody) )*zp
                 cmzp = cmzp + ( xc(i) - xcent(iBody) )*yp &
                             - ( yc(j) - ycent(iBody) )*xp
                IF (zoneSeparate) THEN
                 zone_cmxp(zoneNum) = zone_cmxp(zoneNum) + ( yc(j) - ycent(iBody) )*zp
                 zone_cmyp(zoneNum) = zone_cmyp(zoneNum) - ( xc(i) - xcent(iBody) )*zp
                 zone_cmzp(zoneNum) = zone_cmzp(zoneNum) + ( xc(i) - xcent(iBody) )*yp &
                             - ( yc(j) - ycent(iBody) )*xp
                ENDIF

              ENDIF

              IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                if(boundary_motion_type(ibody)==BIO_FOLLOWED_DYNAMICS_COUPLED)then
                    cmxp = cmxp - ( zc(k) - zcent(1) )*yp
                    cmyp = cmyp + ( zc(k) - zcent(1) )*xp
                    IF (zoneSeparate) THEN
                    zone_cmxp(zoneNum) = zone_cmxp(zoneNum) - ( zc(k) - zcent(1) )*yp
                    zone_cmyp(zoneNum) = zone_cmyp(zoneNum) + ( zc(k) - zcent(1) )*xp
                    ENDIF
                else
                    cmxp = cmxp - ( zc(k) - zcent(iBody) )*yp
                    cmyp = cmyp + ( zc(k) - zcent(iBody) )*xp
                    IF (zoneSeparate) THEN
                    zone_cmxp(zoneNum) = zone_cmxp(zoneNum) - ( zc(k) - zcent(iBody) )*yp
                    zone_cmyp(zoneNum) = zone_cmyp(zoneNum) + ( zc(k) - zcent(iBody) )*xp
                    ENDIF
                end if
              ENDIF

            ELSEIF (Prsb_MomentRef) THEN
              cmxp = cmxp + ( yc(j) - Moment_refy )*zp
              cmyp = cmyp - ( xc(i) - Moment_refx )*zp
              cmzp = cmzp + ( xc(i) - Moment_refx )*yp &
                         - ( yc(j) - Moment_refy )*xp
              IF (zoneSeparate) THEN
              zone_cmxp(zoneNum) = zone_cmxp(zoneNum) + ( yc(j) - Moment_refy )*zp               !added by Chengyu
              zone_cmyp(zoneNum) = zone_cmyp(zoneNum) - ( xc(i) - Moment_refx )*zp
              zone_cmzp(zoneNum) = zone_cmzp(zoneNum) + ( xc(i) - Moment_refx )*yp &
                                                               - ( yc(j) - Moment_refy )*xp
              ENDIF

              IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                cmxp = cmxp - ( zc(k) - Moment_refz )*yp
                cmyp = cmyp + ( zc(k) - Moment_refz )*xp
                IF (zoneSeparate) THEN
                zone_cmxp(zoneNum) = zone_cmxp(zoneNum) - ( zc(k) - Moment_refz )*yp                !added by Chengyu
                zone_cmyp(zoneNum) = zone_cmyp(zoneNum) + ( zc(k) - Moment_refz )*xp
                ENDIF
              ENDIF

            ENDIF

            amx = dxinv(i)*ium(i,j,k)*2.0_CGREAL
            apx = dxinv(i)*iup(i,j,k)*2.0_CGREAL
            amy = dyinv(j)*jum(i,j,k)*2.0_CGREAL
            apy = dyinv(j)*jup(i,j,k)*2.0_CGREAL
            amz = dzinv(k)*kum(i,j,k)*2.0_CGREAL
            apz = dzinv(k)*kup(i,j,k)*2.0_CGREAL

            xs = reinv * dx(i)* dz(k) *                              &
                    ( amy*(    u(i,j,k) -bcyu(i,j,k) )   &
                    - apy*( bcyu(i,j,k) -   u(i,j,k) ) ) &
                    +reinv * dx(i)* dy(j) *                              &
                    ( amz*(    u(i,j,k) -bczu(i,j,k) )   &
                    - apz*( bczu(i,j,k) -   u(i,j,k) ) )
            ys = reinv * dy(j)* dz(k) *                              &
                    ( amx*(    v(i,j,k) -bcxv(i,j,k) )   &
                    - apx*( bcxv(i,j,k) -   v(i,j,k) ) ) &
                    +reinv * dx(i)* dy(j) *                              &
                    ( amz*(    v(i,j,k) -bczv(i,j,k) )   &
                    - apz*( bczv(i,j,k) -   v(i,j,k) ) )
            zs = reinv * dy(j)* dz(k) *                              &
                    ( amx*(    w(i,j,k) -bcxw(i,j,k) )   &
                    - apx*( bcxw(i,j,k) -   w(i,j,k) ) ) &
                    +reinv * dx(i)* dz(k) *                              &
                    ( amy*(    w(i,j,k) -bcyw(i,j,k) )   &
                    - apy*( bcyw(i,j,k) -   w(i,j,k) ) )

            cxs = cxs + xs
            cys = cys + ys
            czs = czs + zs

            bodyMarkerForce(3*mG-2)=xs+bodyMarkerForce(3*mG-2)
            bodyMarkerForce(3*mG-1)=ys+bodyMarkerForce(3*mG-1)
            bodyMarkerForce(3*mG)=zs+bodyMarkerForce(3*mG)
            bodyMarkerForce(3*mG-2)=2.0*bodyMarkerForce(3*mG-2)
            bodyMarkerForce(3*mG-1)=2.0*bodyMarkerForce(3*mG-1)
            bodyMarkerForce(3*mG)=2.0*bodyMarkerForce(3*mG)

            cpw = cpw + (xp+xs)*(uBodyMarker(ibody,mG)-uinit) + &
                      (yp+ys)*(vBodyMarker(ibody,mG)-vinit) + &
                      (zp+zs)*(wBodyMarker(ibody,mG)-winit)

          IF (zoneSeparate) THEN
          zone_cxs(zoneNum) = zone_cxs(zoneNum) + xs        !added by Chengyu
          zone_cys(zoneNum) = zone_cys(zoneNum) + ys
          zone_czs(zoneNum) = zone_czs(zoneNum) + zs

          zone_cpw(zoneNum) = zone_cpw(zoneNum) + (xp+xs)*uBodyMarker(ibody,mG) + &     !added by Chengyu
                                                  (yp+ys)*vBodyMarker(ibody,mG) + &
                                                  (zp+zs)*wBodyMarker(ibody,mG)
          ENDIF

          IF (.NOT. Prsb_MomentRef) THEN
            if(boundary_motion_type(ibody)==BIO_FOLLOWED_DYNAMICS_COUPLED)then
                cmxs = cmxs + ( yc(j) - ycent(1) )*zs
                cmys = cmys - ( xc(i) - xcent(1) )*zs
                cmzs = cmzs + ( xc(i) - xcent(1) )*ys &
                            - ( yc(j) - ycent(1) )*xs
                IF (zoneSeparate) THEN
                  zone_cmxs(zoneNum) = zone_cmxs(zoneNum) + ( yc(j) - ycent(iBody) )*zs         !added by Chengyu
                  zone_cmys(zoneNum) = zone_cmys(zoneNum) - ( xc(i) - xcent(iBody) )*zs
                  zone_cmzs(zoneNum) = zone_cmzs(zoneNum) + ( xc(i) - xcent(iBody) )*ys &
                                                          - ( yc(j) - ycent(iBody) )*xs
                ENDIF
            else
                 cmxs = cmxs + ( yc(j) - ycent(iBody) )*zs
                 cmys = cmys - ( xc(i) - xcent(iBody) )*zs
                 cmzs = cmzs + ( xc(i) - xcent(iBody) )*ys &
                             - ( yc(j) - ycent(iBody) )*xs
                 IF (zoneSeparate) THEN
                 zone_cmxs(zoneNum) = zone_cmxs(zoneNum) + ( yc(j) - ycent(iBody) )*zs
                 zone_cmys(zoneNum) = zone_cmys(zoneNum) - ( xc(i) - xcent(iBody) )*zs
                 zone_cmzs(zoneNum) = zone_cmzs(zoneNum) + ( xc(i) - xcent(iBody) )*ys &
                                                         - ( yc(j) - ycent(iBody) )*xs
                 ENDIF
            end if
            if(boundary_motion_type(ibody)==BIO_FOLLOWED_DYNAMICS_COUPLED)then
                IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                    cmxs = cmxs - ( zc(k) - zcent(1) )*ys
                    cmys = cmys + ( zc(k) - zcent(1) )*xs
                    IF (zoneSeparate) THEN
                    zone_cmxs(zoneNum) = zone_cmxs(zoneNum) - ( zc(k) - zcent(iBody) )*ys       !added by Chengyu
                    zone_cmys(zoneNum) = zone_cmys(zoneNum) + ( zc(k) - zcent(iBody) )*xs
                    ENDIF
                ENDIF
            else
             IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                cmxs = cmxs - ( zc(k) - zcent(iBody) )*ys
                cmys = cmys + ( zc(k) - zcent(iBody) )*xs
             ENDIF
            end if
          ELSEIF (Prsb_MomentRef) THEN
                 cmxs = cmxs + ( yc(j) - Moment_refy )*zs
                 cmys = cmys - ( xc(i) - Moment_refx )*zs
                 cmzs = cmzs + ( xc(i) - Moment_refx )*ys &
                             - ( yc(j) - Moment_refy )*xs
                IF (zoneSeparate) THEN
                 zone_cmxs(zoneNum) =  zone_cmxs(zoneNum) + ( yc(j) - Moment_refy )*zs           !added by Chengyu
                 zone_cmys(zoneNum) =  zone_cmys(zoneNum) - ( xc(i) - Moment_refx )*zs
                 zone_cmzs(zoneNum) =  zone_cmzs(zoneNum) + ( xc(i) - Moment_refx )*ys &
                                                          - ( yc(j) - Moment_refy )*xs
                ENDIF
             IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
                cmxs = cmxs - ( zc(k) - Moment_refz )*ys
                cmys = cmys + ( zc(k) - Moment_refz )*xs
                IF (zoneSeparate) THEN
                zone_cmxs(zoneNum) = zone_cmxs(zoneNum) - ( zc(k) - Moment_refz )*ys            !added by Chengyu
                zone_cmys(zoneNum) = zone_cmys(zoneNum) + ( zc(k) - Moment_refz )*xs
                ENDIF
             ENDIF
          ENDIF

          if (xc(i) - xcent(iBody) < xc_xcent_min) xc_xcent_min = xc(i) - xcent(iBody)
          if (xc(i) - xcent(iBody) > xc_xcent_max) xc_xcent_max = xc(i) - xcent(iBody)
          if (yc(j) - ycent(iBody) < yc_ycent_min) yc_ycent_min = yc(j) - ycent(iBody)
          if (yc(j) - ycent(iBody) > yc_ycent_max) yc_ycent_max = yc(j) - ycent(iBody)

          if (xp > xp_max) xp_max = xp
          if (xp < xp_min) xp_min = xp
          if (yp > yp_max) yp_max = yp
          if (yp < yp_min) yp_min = yp

        ENDIF ! bodyNum

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

      !-----------------------------------------------------------------------
      !calculate power due to inertial force, 2D membrane body   : add by yan
      !-----------------------------------------------------------------------
      if(ndim==DIM_2D)then
        chordLength=sqrt((xBodyMarker(iBody,1)-xBodyMarker(iBody,nPtsBodyMarker(iBody)/3))**twod+&
                    (yBodyMarker(iBody,1)-yBodyMarker(iBody,nPtsBodyMarker(iBody)/3))**twod)
        wingMass=chordLength*thickoverlength*density_solid(iBody)
        wingMOI=wingMass*chordLength**twod/12.0_CGREAL
        wingAlpha=acos((xBodyMarker(iBody,1)-xBodyMarker(iBody,nPtsBodyMarker(iBody)/3))/chordLength)
        xTemp1=xBodyMarker(iBody,1)+dt*uBodyMarker(iBody,1)
        xTemp2=xBodyMarker(iBody,nPtsBodyMarker(iBody)/3)+dt*uBodyMarker(iBody,nPtsBodyMarker(iBody)/3)
        wingAlphaNext=acos((xTemp1-xTemp2)/chordLength)
        wingAngw=(wingAlphaNext-wingAlpha)/dt
        kTrans=half*wingMass*((half*(uBodyMarker(iBody,1)+uBodyMarker(iBody,nPtsBodyMarker(iBody)/3)))**twod+&
                            (half*(vBodyMarker(iBody,1)+vBodyMarker(iBody,nPtsBodyMarker(iBody)/3)))**twod)
        kRot=half*wingMOI*wingAngw**twod
        kRigid=kTrans+kRot
        do i=1,nPtsBodyMarker(iBody)/3-1
            umid=(uBodyMarker(iBody,i)+uBodyMarker(iBody,i+1))*half
            vmid=(vBodyMarker(iBody,i)+vBodyMarker(iBody,i+1))*half

            segmentLength=sqrt((xBodyMarker(iBody,i)-xBodyMarker(iBody,i+1))**twod+&
                            (yBodyMarker(iBody,i)-yBodyMarker(iBody,i+1))**twod)

            segmentMass=segmentLength*thickoverlength*density_solid(ibody)

            kmid=half*segmentMass*(umid**twod+vmid**twod)

            ktotal=ktotal+kmid

        end do
        if(ntime>ntime_start+1)then
            cpwInertiaPre=(ktotal-ktotalPre)/dt*twod
            cpwRigidPre=(kRigid-kRigidPre)/dt*twod
        end if

        kRigidPre=kRigid
        ktotalPre=ktotal
      end if
      !-----------------------------------------------------------------------

      scx(ibody)  = cxp  + cxs         !VEERA ... Fed to Flow Induced Motion subroutine
      scy(ibody)  = cyp  + cys         !
      scz(ibody)  = czp  + czs         !
      scmx(ibody) = cmxp + cmxs        !
      scmy(ibody) = cmyp + cmys        !
      scmz(ibody) = cmzp + cmzs        !VEERA ... Fed to Flow Induced Motion subroutine

!--- Construct Components
      cxp = 2.0_CGREAL*cxp
      cyp = 2.0_CGREAL*cyp
      czp = 2.0_CGREAL*czp
      cxs = 2.0_CGREAL*cxs
      cys = 2.0_CGREAL*cys
      czs = 2.0_CGREAL*czs

      cpw = 2.0_CGREAL*cpw
      cpwAeroPre=cpwAero
      cpwAero=-cpw

      cmxp= 2.0_CGREAL*cmxp
      cmxs= 2.0_CGREAL*cmxs
      cmyp= 2.0_CGREAL*cmyp
      cmys= 2.0_CGREAL*cmys
      cmzp= 2.0_CGREAL*cmzp
      cmzs= 2.0_CGREAL*cmzs
      IF (zoneSeparate) THEN
          zone_cxp(:)  = 2.0_CGREAL*zone_cxp(:)               !added by Chengyu
          zone_cyp(:)  = 2.0_CGREAL*zone_cyp(:)
          zone_czp(:)  = 2.0_CGREAL*zone_czp(:)
          zone_cxs(:)  = 2.0_CGREAL*zone_cxs(:)
          zone_cys(:)  = 2.0_CGREAL*zone_cys(:)
          zone_czs(:)  = 2.0_CGREAL*zone_czs(:)
          zone_cpw(:)  = 2.0_CGREAL*zone_cpw(:)
          zone_cmxp(:) = 2.0_CGREAL*zone_cmxp(:)
          zone_cmxs(:) = 2.0_CGREAL*zone_cmxs(:)
          zone_cmyp(:) = 2.0_CGREAL*zone_cmyp(:)
          zone_cmys(:) = 2.0_CGREAL*zone_cmys(:)
          zone_cmzp(:) = 2.0_CGREAL*zone_cmzp(:)
          zone_cmzs(:) = 2.0_CGREAL*zone_cmzs(:)
      ENDIF

      IF ( ndim == DIM_2D ) THEN

        cxp = cxp/zout
        cyp = cyp/zout
        cxs = cxs/zout
        cys = cys/zout

        cpw = cpw/zout      !add by yan
        cpwAeroPre = cpwAeroPre/zout

        cmzp = cmzp/zout
        cmzs = cmzs/zout

        IF (zoneSeparate) THEN
            zone_cxp(:)  = zone_cxp(:)/zout          !added by Chengyu
            zone_cyp(:)  = zone_cyp(:)/zout
            zone_cxs(:)  = zone_cxs(:)/zout
            zone_cys(:)  = zone_cys(:)/zout
            zone_cmzp(:) = zone_cmzp(:)/zout
            zone_cmzs(:) = zone_cmzs(:)/zout
            zone_cpw(:)  = zone_cpw(:)/zout
        ENDIF
      ENDIF

      cx  = cxp  + cxs         ! for non-cylinders these are
      cy  = cyp  + cys         ! raw
      cz  = czp  + czs         ! forces
      cmx = cmxp + cmxs        ! and moments
      cmy = cmyp + cmys        ! need to non-dmensionalize them
      cmz = cmzp + cmzs        ! during post-processing
      IF (zoneSeparate) THEN
          zone_cx(:)  = zone_cxp(:)  + zone_cxs(:)         ! for non-cylinders these are                !added by Chengyu
          zone_cy(:)  = zone_cyp(:)  + zone_cys(:)         ! raw
          zone_cz(:)  = zone_czp(:)  + zone_czs(:)         ! forces
          zone_cmx(:) = zone_cmxp(:) + zone_cmxs(:)        ! and moments
          zone_cmy(:) = zone_cmyp(:) + zone_cmys(:)        ! need to non-dmensionalize them
          zone_cmz(:) = zone_cmzp(:) + zone_cmzs(:)        ! during post-processing
      ENDIF
!--- Write File
    if(ndim==DIM_2D.and.ntime>ntime_start+1)then
      write(ifuPowerMembrane2D+iBody-1,122) time-dt,cpwAeroPre,cpwInertiaPre,cpwRigidPre,cpwAeroPre+cpwInertiaPre
    end if

    IF (boundary_motion_type(iBody)<=PRESCRIBED) THEN
      WRITE(ifuDragOut+iBody-1,121) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
                                        ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw  !
    ENDIF

    IF(boundary_motion_type(iBody)>PRESCRIBED .and. Converged_FSI(iBody))THEN
       IF (ndim == DIM_2D) THEN
         WRITE(ifuDragOut+iBody-1,121) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
                                        ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw
       ELSEIF(ndim == DIM_3D) THEN
            WRITE(ifuDragOut+iBody-1,121) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
                                         ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw
            write(ifuforce_debug+ibody-1,124) ntime,niterFS,cx,cy,cz,cmx,cmy,cmz   ! Added by Geng for debug
       ENDIF
    ENDIF

    IF(zoneSeparate) THEN
      DO izoneMax=1, zoneMax(ibody)
        WRITE(ifuDragOutZone+izoneMax-1,121) time,zone_cxp(izoneMax),zone_cxs(izoneMax),zone_cx(izoneMax),zone_cyp(izoneMax),zone_cys(izoneMax),zone_cy(izoneMax),zone_czp(izoneMax),zone_czs(izoneMax),zone_cz(izoneMax) &                             !added by Chengyu
                                                      ,zone_cmxp(izoneMax),zone_cmxs(izoneMax),zone_cmx(izoneMax),zone_cmyp(izoneMax),zone_cmys(izoneMax),zone_cmy(izoneMax),zone_cmzp(izoneMax),zone_cmzs(izoneMax),zone_cmz(izoneMax),zone_cpw(izoneMax)
      ENDDO !izoneMax
    ENDIF
    IF(optimization)THEN
      if(mod(nTime,evaluationInterval)<=(evaluationInterval/2).and.mod(nTime,evaluationInterval)/=0)then
        sumV4=sumV4-cx
      else
        sumV4=sumV4+cx
      end if
      sumV7=sumV7+cy
      sumV10=sumV10+cz
      IF(mod(nTime,evaluationInterval)==0)THEN
        optCount=optCount+1
        if(optCount>=2)then
          V4Pre=V4
          V7Pre=V7
          V10Pre=V10
        end if
        V4=sumV4
        V7=sumV7
        V10=sumV10
        sumV4=zero
        sumV7=zero
        sumV10=zero
        write(4030,4040) optCount,V4/real(evaluationInterval,8),V7/real(evaluationInterval,8),V10/real(evaluationInterval,8)
        if(optCount>=2.and.(abs(V7-V7Pre)/abs(V7))<=errorPermitted)then
          optStop=.true.

          open(4040,file='opt_evaluation_output.dat',status='replace')
          write(4040,4040) optCount,V4/real(evaluationInterval,8),V7/real(evaluationInterval,8),V10/real(evaluationInterval,8)
          close(4040)
        end if
      ENDIF
    ENDIF
  END DO ! ibody
4040 format(i4,3f16.8)
121 FORMAT(20(3x,1PE12.5))
122 FORMAT(5(3x,1PE12.5))
124 format(I5,I3,6(1x,1PE12.5))

   print *, 'output bodyMarkerForce for debug'
   open(unit=66,file='IBM_BMF.dat',status='unknown')
   DO i=1,nPtsBodyMarker(1)
     write(66,*) bodyMarkerForce(i*3-2), bodyMarkerForce(i*3-1), bodyMarkerForce(i*3)
   ENDDO
   close(66)

  ! Interpolate markerPressure & markerInterpolateVelocity using markerInterpolateIndex & markerInterpolateRatio
  iBody=1
  DO i=1,nPtsBodyMarker(iBody)
    DO j=1,2
      markerPressure(i*2-2+j)=p(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3))*&
                          markerInterpolateRatio(i*6+j*3-9+1)*markerInterpolateRatio(i*6+j*3-9+2)*markerInterpolateRatio(i*6+j*3-9+3) + &
                          p(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3))*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*markerInterpolateRatio(i*6+j*3-9+2)*markerInterpolateRatio(i*6+j*3-9+3) + &
                          p(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3))*&
                          markerInterpolateRatio(i*6+j*3-9+1)*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*markerInterpolateRatio(i*6+j*3-9+3) + &
                          p(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          markerInterpolateRatio(i*6+j*3-9+1)*markerInterpolateRatio(i*6+j*3-9+2)*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          p(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3))*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*markerInterpolateRatio(i*6+j*3-9+3) + &
                          p(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          markerInterpolateRatio(i*6+j*3-9+1)*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          p(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*markerInterpolateRatio(i*6+j*3-9+2)*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          p(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*(1.0-markerInterpolateRatio(i*6+j*3-9+3))
      markerInterpolateVelocity(i*6+j*3-9+1)=u(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3))*&
                          markerInterpolateRatio(i*6+j*3-9+1)*markerInterpolateRatio(i*6+j*3-9+2)*markerInterpolateRatio(i*6+j*3-9+3) + &
                          u(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3))*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*markerInterpolateRatio(i*6+j*3-9+2)*markerInterpolateRatio(i*6+j*3-9+3) + &
                          u(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3))*&
                          markerInterpolateRatio(i*6+j*3-9+1)*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*markerInterpolateRatio(i*6+j*3-9+3) + &
                          u(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          markerInterpolateRatio(i*6+j*3-9+1)*markerInterpolateRatio(i*6+j*3-9+2)*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          u(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3))*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*markerInterpolateRatio(i*6+j*3-9+3) + &
                          u(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          markerInterpolateRatio(i*6+j*3-9+1)*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          u(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*markerInterpolateRatio(i*6+j*3-9+2)*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          u(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*(1.0-markerInterpolateRatio(i*6+j*3-9+3))
      markerInterpolateVelocity(i*6+j*3-9+2)=v(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3))*&
                          markerInterpolateRatio(i*6+j*3-9+1)*markerInterpolateRatio(i*6+j*3-9+2)*markerInterpolateRatio(i*6+j*3-9+3) + &
                          v(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3))*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*markerInterpolateRatio(i*6+j*3-9+2)*markerInterpolateRatio(i*6+j*3-9+3) + &
                          v(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3))*&
                          markerInterpolateRatio(i*6+j*3-9+1)*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*markerInterpolateRatio(i*6+j*3-9+3) + &
                          v(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          markerInterpolateRatio(i*6+j*3-9+1)*markerInterpolateRatio(i*6+j*3-9+2)*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          v(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3))*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*markerInterpolateRatio(i*6+j*3-9+3) + &
                          v(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          markerInterpolateRatio(i*6+j*3-9+1)*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          v(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*markerInterpolateRatio(i*6+j*3-9+2)*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          v(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*(1.0-markerInterpolateRatio(i*6+j*3-9+3))
      markerInterpolateVelocity(i*6+j*3-9+3)=w(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3))*&
                          markerInterpolateRatio(i*6+j*3-9+1)*markerInterpolateRatio(i*6+j*3-9+2)*markerInterpolateRatio(i*6+j*3-9+3) + &
                          w(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3))*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*markerInterpolateRatio(i*6+j*3-9+2)*markerInterpolateRatio(i*6+j*3-9+3) + &
                          w(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3))*&
                          markerInterpolateRatio(i*6+j*3-9+1)*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*markerInterpolateRatio(i*6+j*3-9+3) + &
                          w(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          markerInterpolateRatio(i*6+j*3-9+1)*markerInterpolateRatio(i*6+j*3-9+2)*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          w(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3))*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*markerInterpolateRatio(i*6+j*3-9+3) + &
                          w(markerInterpolateIndex(i*6+j*3-9+1),markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          markerInterpolateRatio(i*6+j*3-9+1)*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          w(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2),markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*markerInterpolateRatio(i*6+j*3-9+2)*(1.0-markerInterpolateRatio(i*6+j*3-9+3)) + &
                          w(markerInterpolateIndex(i*6+j*3-9+1)+1,markerInterpolateIndex(i*6+j*3-9+2)+1,markerInterpolateIndex(i*6+j*3-9+3)+1)*&
                          (1.0-markerInterpolateRatio(i*6+j*3-9+1))*(1.0-markerInterpolateRatio(i*6+j*3-9+2))*(1.0-markerInterpolateRatio(i*6+j*3-9+3))
    ENDDO
  ENDDO

END SUBROUTINE  drag_lift_membrane
!-------------------------------------------------------------------------------

   SUBROUTINE GCM_drag_lift()
!
!   Compute the Lift and Drag coefficients for n-Bodies in the flow
!
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    USE usr_module ,ONLY : scx,scy,scz,scmx,scmy,scmz	!VEERA

    IMPLICIT NONE

    INTEGER           :: i,j,k
    INTEGER           :: ibody,ifudrag
    INTEGER           :: iG, jG, kG, indZ
    INTEGER           :: nG,m ,m1, m2, mPointsOrig, nTrielemMax
    INTEGER           :: ii,jj,kk,n,iRow

    REAL(KIND=CGREAL) :: uIP, vIP, wIP,        &
                         uBI, vBI, wBI,        &
                         dUt1Dn,dUt2Dn,        &
                         alphaX,alphaY,alphaZ, &
                         rectElemCentX,rectElemCentY,rectElemCentZ,rectElemArea, &
                         dist, dist0, dist1, distMin, distMin0, distMin1, xp, yp, zp, xs, ys, zs

    REAL(KIND=CGREAL) :: p_inv_dist_int0,p_inv_dist_int1,dist_int0,dist_int1


    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: cxp,cyp,czp,     &
                                                    cxs,cys,czs,     &
                                                    cmxp,cmyp,cmzp,  &
                                                    cmxs,cmys,cmzs,  &
                                                    cx,cy,cz,        &
                                                    cmx,cmy,cmz


    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: cxpSurf, cypSurf,    &
    						    cxsSurf, cysSurf

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:) :: InterceptElemCentX,  &
						    InterceptElemCentY

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:) :: cxpTemp, cypTemp,  &
    						      cxsTemp, cysTemp

    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:) :: pTriElemCent

    INTEGER(1), DIMENSION(:,:,:), POINTER :: iblankTemp

    CHARACTER*9          :: dragfile
    CHARACTER*25         :: indragfile

    CHARACTER*16         :: Ptname1		!output file names  ..abel
    CHARACTER*16         :: Pbname1

    ifudrag = 155

! allocate local arrays

    ALLOCATE(cxp(nBody))
    ALLOCATE(cyp(nBody))
    ALLOCATE(czp(nBody))

    ALLOCATE(cxs(nBody))
    ALLOCATE(cys(nBody))
    ALLOCATE(czs(nBody))

    ALLOCATE(cxpSurf(nPtsBodyMarkerOrig(nBody)))	!... abel
    ALLOCATE(cypSurf(nPtsBodyMarkerOrig(nBody)))
    ALLOCATE(cxsSurf(nPtsBodyMarkerOrig(nBody)))	!... abel
    ALLOCATE(cysSurf(nPtsBodyMarkerOrig(nBody)))

    ALLOCATE(cxpTemp(nPtsBodyMarkerOrig(nBody),nz))	!... abel
    ALLOCATE(cypTemp(nPtsBodyMarkerOrig(nBody),nz))
    ALLOCATE(cxsTemp(nPtsBodyMarkerOrig(nBody),nz))	!... abel
    ALLOCATE(cysTemp(nPtsBodyMarkerOrig(nBody),nz))

    ALLOCATE(InterceptElemCentX(nPtsBodyMarkerOrig(nBody)))
    ALLOCATE(InterceptElemCentY(nPtsBodyMarkerOrig(nBody)))

    ALLOCATE(cmxp(nBody))
    ALLOCATE(cmyp(nBody))
    ALLOCATE(cmzp(nBody))

    ALLOCATE(cmxs(nBody))
    ALLOCATE(cmys(nBody))
    ALLOCATE(cmzs(nBody))

    ALLOCATE(cx(nBody))
    ALLOCATE(cy(nBody))
    ALLOCATE(cz(nBody))

    ALLOCATE(cmx(nBody))
    ALLOCATE(cmy(nBody))
    ALLOCATE(cmz(nBody))

    nTriElemMax = MAXVAL(totNumTriElem(:))
    ALLOCATE(pTriElemCent(nBody,nTriElemMax))

! initialize variables

    cxp = zero
    cyp = zero
    czp = zero

    cxs = zero
    cys = zero
    czs = zero

    cxpSurf = zero		!...abel
    cypSurf = zero
    cxsSurf = zero		!...abel
    cysSurf = zero

    cxpTemp = zero		!...abel
    cypTemp = zero
    cxsTemp = zero		!...abel
    cysTemp = zero

    InterceptElemCentX = zero	!...abel
    InterceptElemCentY = zero

    cmxp = zero
    cmyp = zero
    cmzp = zero

    cmxs = zero
    cmys = zero
    cmzs = zero

    cx = zero
    cy = zero
    cz = zero

    cmx = zero
    cmy = zero
    cmz = zero

    DO iBody = 1,nBody

          IF (unstruc_surface_type(iBody) == SOLID_BODY) THEN
              iblankTemp => iblank_solid
          ELSE
              iblankTemp => iblank_memb
          ENDIF

      SELECT CASE( ndim )
        CASE(DIM_2D)

          mPointsOrig = totNumTriElem(iBody)/2/(nz-1)

          DO indZ = 1,nz-1

           DO m=1,mPointsOrig   ! looping over original marker points

             m1 = 2*mPointsOrig*(indZ-1)+(2*m-1)      !counter results in odd values
	        m2 = m1+1				      !counter results in even values

  ! compute pressure at centroid of rectangle formed by these two triangular elements
             rectElemCentX = half*(triElemCentX(iBody,m1) + triElemCentX(iBody,m2))
             rectElemCentY = half*(triElemCentY(iBody,m1) + triElemCentY(iBody,m2))
             rectElemCentZ = half*(triElemCentZ(iBody,m1) + triElemCentZ(iBody,m2))

	        rectElemArea  = triElemArea(iBody,m1) + triElemArea(iBody,m2)

  ! Storing body intercept coordinates 			...abel
	        InterceptElemCentX(m) = rectElemCentX
	        InterceptElemCentY(m) = rectElemCentY

!  ! search grid to find nodes surrounding element
!             DO i = 1, nx-1
!               IF ( ( xc(i) <= rectElemCentX ) .AND. ( xc(i+1) > rectElemCentX ) ) ii = i
!             ENDDO
!             DO j = 1,ny-1
!               IF ( ( yc(j) <= rectElemCentY ) .AND. ( yc(j+1) > rectElemCentY ) ) jj = j
!             ENDDO
!             kk = indz
!
!! inverse distance sqr weighted interpolation
!             p_inv_dist_int0 = zero
!             p_inv_dist_int1 = zero
!             dist_int0       = zero
!             dist_int1       = zero
!
!             DO j=0,1
!             DO i=0,1
!
!! Pressure from the iblank=0 side
!               IF (  iblankTemp(ii+i,jj+j,kk) == 0 .OR.  &
!                    (iblankTemp(ii+i,jj+j,kk) == 1 .AND. ghostCellMark(ii+i,jj+j,kk) == 1)  ) THEN
!                  dist0     =  (xc(ii+i)-rectElemCentX)**2 &
!                              +(yc(jj+j)-rectElemCentY)**2 &
!                              +(zc(kk)  -rectElemCentZ)**2
!                  dist_int0 = dist_int0 + (oned/dist0)
!                  IF (iblankTemp(ii+i,jj+j,kk) == 1 .AND. ghostCellMark(ii+i,jj+j,kk) == 1) THEN
!                    p_inv_dist_int0 = p_inv_dist_int0  &
!                                         + (oned/dist0)*p(ii+i,jj+j,kk)
!                  ELSE
!                    p_inv_dist_int0 = p_inv_dist_int0  &
!                                         + (oned/dist0)*p(ii+i,jj+j,kk)
!                  ENDIF
!               ENDIF
!
!! Pressure from the iblank=1 side  (relevant for membrane)
!               IF (  iblankTemp(ii+i,jj+j,kk) == 1 .OR.  &
!                   (iblankTemp(ii+i,jj+j,kk) == 0 .AND. ghostCellMark(ii+i,jj+j,kk) == 1)  ) THEN
!                  dist1     =  (xc(ii+i)-rectElemCentX)**2 &
!                              +(yc(jj+j)-rectElemCentY)**2 &
!                              +(zc(kk)  -rectElemCentZ)**2
!                  dist_int1 = dist_int1 + (oned/dist1)
!                  IF (iblankTemp(ii+i,jj+j,kk) == 0 .AND. ghostCellMark(ii+i,jj+j,kk) == 1) THEN
!                    p_inv_dist_int1 = p_inv_dist_int1  &
!                                         + (oned/dist1)*p(ii+i,jj+j,kk)
!                  ELSE
!                    p_inv_dist_int1 = p_inv_dist_int1  &
!                                         + (oned/dist1)*p(ii+i,jj+j,kk)
!                  ENDIF
!               ENDIF
!
!             ENDDO
!             ENDDO
!
!             pTriElemCent0(iBody,m) = p_inv_dist_int0/dist_int0
!
!             IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
!               pTriElemCent1(iBody,m) = p_inv_dist_int1/dist_int1
!             ELSE
!               pTriElemCent1(iBody,m) = zero
!             ENDIF

! find closest ghost node
! for membranes, need closest ghost nodes on both sides
             distMin = 1.0E8_CGREAL
             DO n=1,nGhost
               i = iGhost(n)
               j = jGhost(n)
               k = kGhost(n)

               dist = (xc(i)-rectElemCentX)**2 &
                     +(yc(j)-rectElemCentY)**2 &
                     +(zc(k)-rectElemCentZ)**2

               IF ( dist <= distMin ) THEN
                  distMin = dist
                  nG  = n
                  iG  = i
                  jG  = j
                  kG  = k
	       ENDIF
             ENDDO

             IF (kG /= indZ) THEN
                PRINT*,' WARNING #################something wrong in lift-drag routine'
                PRINT*,nG,iG,jG,kG
             ENDIF

!             IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
!
!               distMin1 = 1.0E8_CGREAL
!
!               DO n=1,nGhost
!
!                 i = iGhost(n)
!                 j = jGhost(n)
!                 k = kGhost(n)
!
!                 IF (iblankTemp(i,j,k)==0) THEN
!                   dist1 = (xc(i)-rectElemCentX)**2 &
!                          +(yc(j)-rectElemCentY)**2 &
!                          +(zc(k)-rectElemCentZ)**2
!
!                   IF ( dist1 <= distMin1 ) THEN
!                      distMin1 = dist1
!                      nG1  = n
!                      iG1  = i
!                      jG1  = j
!                      kG1  = k
!	           ENDIF
!                 ENDIF
!
!               ENDDO
!
!               IF (kG1 /= indZ) THEN
!                 PRINT*,' WARNING #################something wrong in lift-drag routine'
!                 PRINT*,nG1,iG1,jG1,kG1
!               ENDIF
!
!             ENDIF

  ! search vicinity of closest ghost node to find nodes surrounding element
             DO i = iG-3, iG+3
               IF ( ( xc(i) <= rectElemCentX ) .AND. ( xc(i+1) > rectElemCentX ) ) ii = i
             ENDDO
             DO j = jG-3, jG+3
               IF ( ( yc(j) <= rectElemCentY ) .AND. ( yc(j+1) > rectElemCentY ) ) jj = j
             ENDDO
             kk = indz

  ! trilinear interpolation
	     alphaX = (rectElemCentX-xc(ii))*dxinv(ii)
	     alphaY = (rectElemCentY-yc(jj))*dyinv(jj)
	     alphaZ = zero

             pTriElemCent(iBody,m)  &
                  = p(ii,jj,kk)*(oned - alphaX)*(oned - alphaY)*(oned - alphaZ) &
                  + p(ii+1,jj,kk)*alphaX*(oned - alphaY)*(oned - alphaZ)              &
                  + p(ii,jj+1,kk)*(oned - alphaX)*alphaY*(oned - alphaZ)              &
                  + p(ii,jj,kk+1)*(oned - alphaX)*(oned - alphaY)*alphaZ              &
                  + p(ii+1,jj+1,kk)*alphaX*alphaY*(oned - alphaZ)                           &
                  + p(ii+1,jj,kk+1)*alphaX*(oned - alphaY)*alphaZ                           &
                  + p(ii,jj+1,kk+1)*(oned - alphaX)*alphaY*alphaZ                           &
                  + p(ii+1,jj+1,kk+1)*alphaX*alphaY*alphaZ

!---- Compute shear stress at element centroid

  ! find velocity at image point corresponding to ghost point

             uIP = zero
             vIP = zero
             wIP = zero

             CALL GCM_calc_BIVelocity_Unstruc( iG, jG, kG, xC(iG), yC(jG), zC(kG),  &
                                               closestElementGC(nG),                &
                                               uBodyIntercept(nG),                  &
                                               vBodyIntercept(nG),                  &
                                               wBodyIntercept(nG)                     )
             DO iRow = 1, iRowMax

               ii = iCellIndex(nG) + incI(iRow)
               jj = jCellIndex(nG) + incJ(iRow)
               kk = kCellIndex(nG) + incK(iRow)


               IF ( ii /= iG .OR. jj /= jG .OR. kk /= kG) THEN
                 uIP = uIP + coeffGCMD(iRow,nG)* u(ii,jj,kk)
                 vIP = vIP + coeffGCMD(iRow,nG)* v(ii,jj,kk)
                 wIP = wIP + coeffGCMD(iRow,nG)* w(ii,jj,kk)
               ELSE
                 uIP = uIP + coeffGCMD(iRow,nG)* uBodyIntercept(nG)
                 vIP = vIP + coeffGCMD(iRow,nG)* vBodyIntercept(nG)
                 wIP = wIP + coeffGCMD(iRow,nG)* wBodyIntercept(nG)
               ENDIF ! ii

             ENDDO ! iRow

  ! compute normal-gradient of tangential velocity  at body intercept


             dUt1Dn = (   (uIP-u(iG,jG,kG))*triElemTang1x(iBody,m1) &
                         +(vIP-v(iG,jG,kG))*triElemTang1y(iBody,m1) &
                         +(wIP-w(iG,jG,kG))*triElemTang1z(iBody,m1)   )/probeLength(nG)
             dUt2Dn = (   (uIP-u(iG,jG,kG))*triElemTang2x(iBody,m1) &
                         +(vIP-v(iG,jG,kG))*triElemTang2y(iBody,m1) &
                         +(wIP-w(iG,jG,kG))*triElemTang2z(iBody,m1)   )/probeLength(nG)

!              IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
!                ELSE
!
!               dUt1Dn1 = zero
!               dUt2Dn1 = zero
!
!             ENDIF

!---- Pressure Force Term
!   C_p = Integral (p n.ds)  (no negative sign since n acts into the body)

             xp =  pTriElemCent(iBody,m)*rectElemArea*triElemNormx(iBody,m1)
             yp =  pTriElemCent(iBody,m)*rectElemArea*triElemNormy(iBody,m1)
             zp =  pTriElemCent(iBody,m)*rectElemArea*triElemNormz(iBody,m1)

	     cxp(iBody) = cxp(iBody) + xp
             cyp(iBody) = cyp(iBody) + yp
             czp(iBody) = czp(iBody) + zp

	     cmxp(iBody) = cmxp(iBody) + ( yc(j) - ycent(iBody) )*zp
             cmyp(iBody) = cmyp(iBody) - ( xc(i) - xcent(iBody) )*zp
             cmzp(iBody) = cmzp(iBody) + ( xc(i) - xcent(iBody) )*yp - ( yc(j) - ycent(iBody) )*xp

!---- Shear Stress Term
!   C_s = Integral(tau.n ds)

             xs =  reinv*(dUt1Dn*triElemTang1x(iBody,m1) + dUt2Dn*triElemTang2x(iBody,m1))*rectElemArea
             ys =  reinv*(dUt1Dn*triElemTang1y(iBody,m1) + dUt2Dn*triElemTang2y(iBody,m1))*rectElemArea
             zs =  reinv*(dUt1Dn*triElemTang1z(iBody,m1) + dUt2Dn*triElemTang2z(iBody,m1))*rectElemArea

             cxs(iBody) = cxs(iBody) + xs
             cys(iBody) = cys(iBody) + ys
             czs(iBody) = czs(iBody) + zs

             cmxs(iBody) = cmxs(iBody) + ( yc(j) - ycent(iBody) )*zs
             cmys(iBody) = cmys(iBody) - ( xc(i) - xcent(iBody) )*zs
             cmzs(iBody) = cmzs(iBody) + ( xc(i) - xcent(iBody) )*ys - ( yc(j) - ycent(iBody) )*xs


  !~~~~~~~~~~~~~~ Storing Surface Components ~~~~~~~~~~~~~~~~~~~~~~~~~~
	    cxpTemp(m,indZ) = 2.0_CGREAL*xp/rectElemArea / (nz-1)		!remove Area, want comp. only
    	cypTemp(m,indZ) = 2.0_CGREAL*yp/rectElemArea / (nz-1)	 	! values normalized
    	cxsTemp(m,indZ) = 2.0_CGREAL*xs/rectElemArea / (nz-1)	 	!...abel
   	    cysTemp(m,indZ) = 2.0_CGREAL*ys/rectElemArea / (nz-1)

           ENDDO ! m - markerPts loop

          ENDDO  ! indZ


  !--------Summing up components along the span---------	...abel
	DO m = 1, mPointsOrig
		cxpSurf(m) = cxpTemp(m,nz-2) + cxpTemp(m,nz-1)
		cypSurf(m) = cypTemp(m,nz-2) + cypTemp(m,nz-1)
		cxsSurf(m) = cxsTemp(m,nz-2) + cxsTemp(m,nz-1)
		cysSurf(m) = cysTemp(m,nz-2) + cysTemp(m,nz-1)
	ENDDO



!******************************************************************************************
!*********** Printing terms along the surface --- modified by abel **************************

 IF ( MOD(ntime,nrestart) == 0 .OR. ntime==ntime_start+no_tsteps) THEN
 !IF ( MOD(ntime,ndump) == 0 )    THEN

    IF ( ntime >= 0 .AND. ntime .le. 9  )            &
        write(Ptname1,341)ntime
    IF ( ntime >= 10 .AND. ntime .le.99 )            &
        write(Ptname1,342)ntime
    IF ( ntime >= 100 .AND. ntime .le. 999 )         &
        write(Ptname1,343)ntime
    IF ( ntime >= 1000 .AND. ntime .le. 9999 )       &
        write(Ptname1,344)ntime
    IF ( ntime >= 10000 .AND. ntime .le. 99999 )     &
        write(Ptname1,345)ntime
    IF ( ntime >= 100000 .AND. ntime .le. 999999 )   &
        write(Ptname1,346)ntime
    IF ( ntime >= 1000000 .AND. ntime .le. 9999999 ) &
        write(Ptname1,347)ntime
341     format('Ptop.000000',i1)
342     format('Ptop.00000',i2)
343     format('Ptop.0000',i3)
344     format('Ptop.000',i4)
345     format('Ptop.00',i5)
346     format('Ptop.0',i6)
347     format('Ptop.',i7)


   IF ( ntime >= 0 .AND. ntime .le. 9  )            &
        write(Pbname1,441)ntime
   IF ( ntime >= 10 .AND. ntime .le.99 )            &
        write(Pbname1,442)ntime
   IF ( ntime >= 100 .AND. ntime .le. 999 )         &
        write(Pbname1,443)ntime
   IF ( ntime >= 1000 .AND. ntime .le. 9999 )       &
        write(Pbname1,444)ntime
   IF ( ntime >= 10000 .AND. ntime .le. 99999 )     &
        write(Pbname1,445)ntime
   IF ( ntime >= 100000 .AND. ntime .le. 999999 )   &
        write(Pbname1,446)ntime
   IF ( ntime >= 1000000 .AND. ntime .le. 9999999 ) &
        write(Pbname1,447)ntime
441     format('Pbot.000000',i1)
442     format('Pbot.00000',i2)
443     format('Pbot.0000',i3)
444     format('Pbot.000',i4)
445     format('Pbot.00',i5)
446     format('Pbot.0',i6)
447     format('Pbot.',i7)


    OPEN(UNIT=1,FILE=Ptname1,STATUS='UNKNOWN')		!top surface
    OPEN(UNIT=2,FILE=Pbname1,STATUS='UNKNOWN')		!bott surface

      WRITE(1,*)'VARIABLES = x/c y Pressure C_x_p C_y_p C_x_s C_y_s'
      WRITE(2,*)'VARIABLES = x/c y Pressure C_x_p C_y_p C_x_s C_y_s'

      WRITE(*,*)'  '
      WRITE(*,*)'****** Saving pressure distribution  *******'

   !-------Determine top & bottom surface by the sign of the triElemNormy --------------

      DO m=1,mPointsOrig

	if(triElemNormy(iBody,2*m-1) <= 0)then			!determining top surface
	   InterceptElemCentX(m) = InterceptElemCentX(m) - 4

	   WRITE(1,122) InterceptElemCentX(m), InterceptElemCentY(m), pTriElemCent(iBody,m), &
	   		cxpSurf(m), cypSurf(m),cxsSurf(m),cysSurf(m)
	endif

	if(triElemNormy(iBody,2*m-1) > 0)then 			!determinig bottom surface
	   InterceptElemCentX(m) = InterceptElemCentX(m) - 4

           WRITE(2,122) InterceptElemCentX(m), InterceptElemCentY(m), pTriElemCent(iBody,m),  &
	   		cxpSurf(m),cypSurf(m),cxsSurf(m),cysSurf(m)
	endif

      ENDDO

    CLOSE(1)
    CLOSE(2)

ENDIF

!******************************* END OF 2D CASE  ****************************************
!****************************************************************************************

        CASE(DIM_3D)

          IF(canonical_body_type(iBody) ==  ELLIPSOID) THEN
            PRINT*,'CANNOT COMPUTE SURFACE QUANTITIES FOR "ELLIPSOID"'
            PRINT*,'ABORTING RUN'
            STOP
	  ENDIF

          DO m=1,totNumTriElem(iBody)


  ! find closest ghost node
             distMin = 1.0E8_CGREAL
             DO n=1,nGhost
               i = iGhost(n)
               j = jGhost(n)
               k = kGhost(n)

               dist = (xc(i)-triElemCentX(iBody,m))**2 &
                     +(yc(j)-triElemCentY(iBody,m))**2 &
                     +(zc(k)-triElemCentZ(iBody,m))**2

               IF ( dist <= distMin ) THEN
                  distMin = dist
                  nG  = n
                  iG  = i
                  jG  = j
                  kG  = k
	       ENDIF
             ENDDO

  ! search vicinity of closest ghost node to find nodes surrounding element
             DO i = iG-3, iG+3
               IF ( ( xc(i)-triElemCentX(iBody,m) )*(xc(i+1)-triElemCentX(iBody,m) ) <= zero ) ii = i
             ENDDO
             DO j = jG+3, jG-3
               IF ( ( yc(j)-triElemCentY(iBody,m) )*(yc(j+1)-triElemCentY(iBody,m) ) <= zero ) jj = j
             ENDDO
             DO k = kG+3, kG-3
               IF ( ( zc(k)-triElemCentZ(iBody,m) )*(zc(k+1)-triElemCentZ(iBody,m) ) <= zero ) kk = k
             ENDDO

! compute pressure at element centroid
  ! trilinear interpolation
	     alphaX = (triElemCentX(iBody,m)-xc(ii))*dxinv(ii)
	     alphaY = (triElemCentY(iBody,m)-yc(jj))*dyinv(jj)
	     alphaZ = (triElemCentZ(iBody,m)-zc(kk))*dzinv(kk)

             pTriElemCent(iBody,m)  &
                  = p(ii,jj,kk)*(oned - alphaX)*(oned - alphaY)*(oned - alphaZ) &
                  + p(ii+1,jj,kk)*alphaX*(oned - alphaY)*(oned - alphaZ)             &
                  + p(ii,jj+1,kk)*(oned - alphaX)*alphaY*(oned - alphaZ)             &
                  + p(ii,jj,kk+1)*(oned - alphaX)*(oned - alphaY)*alphaZ             &
                  + p(ii+1,jj+1,kk)*alphaX*alphaY*(oned - alphaZ)                           &
                  + p(ii+1,jj,kk+1)*alphaX*(oned - alphaY)*alphaZ                           &
                  + p(ii,jj+1,kk+1)*(oned - alphaX)*alphaY*alphaZ                          &
                  + p(ii+1,jj+1,kk+1)*alphaX*alphaY*alphaZ


! compute shear stress at element centroid

  ! find velocity at image point corresponding to ghost point

             uIP = zero
             vIP = zero
             wIP = zero

             DO iRow = 1,iRowMax

               ii = iCellIndex(nG) + incI(iRow)
               jj = jCellIndex(nG) + incJ(iRow)
               kk = kCellIndex(nG) + incK(iRow)

               IF ( ii /= iG .OR. jj /= jG .OR. kk /= kG) THEN
                 uIP = uIP + coeffGCMD(iRow,nG)* u(ii,jj,kk)
                 vIP = vIP + coeffGCMD(iRow,nG)* v(ii,jj,kk)
                 wIP = wIP + coeffGCMD(iRow,nG)* w(ii,jj,kk)
               ELSE
                 uIP = uIP + coeffGCMD(iRow,nG)* uBodyIntercept(nG)
                 vIP = vIP + coeffGCMD(iRow,nG)* vBodyIntercept(nG)
                 wIP = wIP + coeffGCMD(iRow,nG)* wBodyIntercept(nG)
               ENDIF ! ii

             ENDDO ! iRow

  ! compute normal-gradient of tangential velocity  at body intercept

             dUt1Dn = (   (uIP-u(iG,jG,kG))*triElemTang1x(iBody,m) &
                         +(vIP-v(iG,jG,kG))*triElemTang1y(iBody,m) &
                         +(wIP-w(iG,jG,kG))*triElemTang1z(iBody,m)   )/probeLength(nG)
             dUt2Dn = (   (uIP-u(iG,jG,kG))*triElemTang2x(iBody,m) &
                         +(vIP-v(iG,jG,kG))*triElemTang2y(iBody,m) &
                         +(wIP-w(iG,jG,kG))*triElemTang2z(iBody,m)   )/probeLength(nG)

!---- Pressure Force Term
!   C_p = Integral (p n.ds)  (no negative sign since n acts into the body)

             xp = pTriElemCent(iBody,m)*triElemArea(iBody,m)*triElemNormx(iBody,m)
             yp = pTriElemCent(iBody,m)*triElemArea(iBody,m)*triElemNormy(iBody,m)
             zp = pTriElemCent(iBody,m)*triElemArea(iBody,m)*triElemNormz(iBody,m)

             cxp(iBody) = cxp(iBody) + xp
             cyp(iBody) = cyp(iBody) + yp
             czp(iBody) = czp(iBody) + zp

             cmxp = cmxp + ( yc(j) - ycent(iBody) )*zp - ( zc(k) - zcent(iBody) )*yp
             cmyp = cmyp - ( xc(i) - xcent(iBody) )*zp + ( zc(k) - zcent(iBody) )*xp
             cmzp = cmzp + ( xc(i) - xcent(iBody) )*yp - ( yc(j) - ycent(iBody) )*xp

!---- Shear Stress Term
!   C_s = Integral(tau.n ds)

             xs =  reinv*(dUt1Dn*triElemTang1x(iBody,m) + dUt2Dn*triElemTang2x(iBody,m))*triElemArea(iBody,m)
             ys =  reinv*(dUt1Dn*triElemTang1y(iBody,m) + dUt2Dn*triElemTang2y(iBody,m))*triElemArea(iBody,m)
             zs =  reinv*(dUt1Dn*triElemTang1z(iBody,m) + dUt2Dn*triElemTang2z(iBody,m))*triElemArea(iBody,m)

             cxs(iBody) = cxs(iBody) + xs
             cys(iBody) = cys(iBody) + ys
             czs(iBody) = czs(iBody) + zs

             cmxs = cmxs + ( yc(j) - ycent(iBody) )*zs - ( zc(k) - zcent(iBody) )*ys
             cmys = cmys - ( xc(i) - xcent(iBody) )*zs + ( zc(k) - zcent(iBody) )*xs
             cmzs = cmzs + ( xc(i) - xcent(iBody) )*ys - ( yc(j) - ycent(iBody) )*xs

        ENDDO ! m

    END SELECT

    ENDDO ! iBody

    DO iBody = 1, nBody

!--- Construct  Coefficients
!    Note    Cf = F / (0.5 rho U^2 Af )      where Af is frontal area
!    Note    Cm = M / (0.5 rho U^2 Af L )    where L is lenght scale
!    Here we assume that  U = L = 1,
!    Af = 1 x zout for 2D ;  Af = 1 for 3D
!    If this is not the case, then divide by appropriate values

	scx(iBody) = cxp(iBody)+cxs(iBody)    !VEERA ... Fed to Flow Induced Motion subroutine
        scy(iBody) = cyp(iBody)+cys(iBody)    !VEERA ... Fed to Flow Induced Motion subroutine
        scz(iBody) = czp(iBody)+czs(iBody)    !VEERA ... Fed to Flow Induced Motion subroutine

        scmx(iBody) = cmxp(iBody)+cmxs(iBody) !VEERA ... Fed to Flow Induced Motion subroutine
        scmy(iBody) = cmyp(iBody)+cmys(iBody) !VEERA ... Fed to Flow Induced Motion subroutine
        scmz(iBody) = cmzp(iBody)+cmzs(iBody) !VEERA ... Fed to Flow Induced Motion subroutine

        cxp(iBody) = 2.0_CGREAL*cxp(iBody)
        cyp(iBody) = 2.0_CGREAL*cyp(iBody)
        czp(iBody) = 2.0_CGREAL*czp(iBody)

        cxs(iBody) = 2.0_CGREAL*cxs(iBody)
        cys(iBody) = 2.0_CGREAL*cys(iBody)
        czs(iBody) = 2.0_CGREAL*czs(iBody)

        cmxp(iBody) = 2.0_CGREAL*cmxp(iBody)
        cmyp(iBody) = 2.0_CGREAL*cmyp(iBody)
        cmzp(iBody) = 2.0_CGREAL*cmzp(iBody)

        cmxs(iBody) = 2.0_CGREAL*cmxs(iBody)
        cmys(iBody) = 2.0_CGREAL*cmys(iBody)
        cmzs(iBody) = 2.0_CGREAL*cmzs(iBody)

        IF ( ndim == DIM_2D ) THEN
          cxp(iBody)  = cxp(iBody)/zout
          cyp(iBody)  = cyp(iBody)/zout
          cxs(iBody)  = cxs(iBody)/zout
          cys(iBody)  = cys(iBody)/zout
          cmzp(iBody) = cmzp(iBody)/zout
          cmzs(iBody) = cmzs(iBody)/zout
        ENDIF

        cx(iBody) = cxp(iBody)+cxs(iBody)
        cy(iBody) = cyp(iBody)+cys(iBody)
        cz(iBody) = czp(iBody)+czs(iBody)

        cmx(iBody) = cmxp(iBody)+cmxs(iBody)
        cmy(iBody) = cmyp(iBody)+cmys(iBody)
        cmz(iBody) = cmzp(iBody)+cmzs(iBody)

!--- Write File

        WRITE(ifuDragOut+iBody-1,121) time,                                    &
                                    cxp(iBody), cxs(iBody), cx(iBody),  &
                                    cyp(iBody), cys(iBody), cy(iBody),  &
                                    czp(iBody), czs(iBody), cz(iBody),  &
                                    cmxp(iBody),cmxs(iBody),cmx(iBody), &
                                    cmyp(iBody),cmys(iBody),cmy(iBody), &
                                    cmzp(iBody),cmzs(iBody),cmz(iBody), &
                                    surfArea(iBody)

    ENDDO ! iBody

! deallocate local arrays

    DEALLOCATE(cxp)
    DEALLOCATE(cyp)
    DEALLOCATE(czp)

    DEALLOCATE(cxs)
    DEALLOCATE(cys)
    DEALLOCATE(czs)

    DEALLOCATE(cxpSurf)			!... abel
    DEALLOCATE(cypSurf)
    DEALLOCATE(cxsSurf)			!... abel
    DEALLOCATE(cysSurf)

    DEALLOCATE(cxpTemp)			!... abel
    DEALLOCATE(cypTemp)
    DEALLOCATE(cxsTemp)
    DEALLOCATE(cysTemp)

    DEALLOCATE(InterceptElemCentX)	!... abel
    DEALLOCATE(InterceptElemCentY)

    DEALLOCATE(cmxp)
    DEALLOCATE(cmyp)
    DEALLOCATE(cmzp)

    DEALLOCATE(cmxs)
    DEALLOCATE(cmys)
    DEALLOCATE(cmzs)

    DEALLOCATE(cx)
    DEALLOCATE(cy)
    DEALLOCATE(cz)

    DEALLOCATE(cmx)
    DEALLOCATE(cmy)
    DEALLOCATE(cmz)

    DEALLOCATE(pTriElemCent)

121 FORMAT(20(3x,1PE12.5))
122 FORMAT(20(3x,1PE14.7))		!...abel

   END SUBROUTINE  GCM_drag_lift
!-------------------------------------------------------------------------------


   SUBROUTINE drag_lift_membrane_partial_dynamics(iBody)
!
!   Compute the Lift and Drag coefficients for n-Bodies in the flow
!
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE usr_module ,ONLY : scx,scy,scz,scmx,scmy,scmz	!VEERA
    USE fea_unstructure_surface
    USE body_dynamics

    IMPLICIT NONE

    INTEGER              :: i,j,k
    INTEGER              :: iBody, MG, m, mD, n
    INTEGER              :: K_ini, K_end
    INTEGER              :: iSection
    INTEGER              :: iG,jG,kG,nG
    INTEGER              :: ih,jh,kh

    REAL(KIND=CGREAL)    :: amx,apx,acx
    REAL(KIND=CGREAL)    :: amy,apy,acy
    REAL(KIND=CGREAL)    :: amz,apz,acz

    REAL(KIND=CGREAL)    :: cxp,cyp,czp,cxs,cys,czs
    REAL(KIND=CGREAL)    :: cmxp,cmyp,cmzp,cmxs,cmys,cmzs
    REAL(KIND=CGREAL)    :: cx,cy,cz
    REAL(KIND=CGREAL)    :: cmx,cmy,cmz, cpw, distMin, dist
    REAL(KIND=CGREAL)    :: xp,yp,zp,xs,ys,zs
    REAL(KIND=CGREAL)    :: hingex, hingey, hingez
!    REAL(KIND=CGREAL)    :: moment_refx, moment_refy, moment_refz
    REAL(KIND=CGREAL)    :: xmarker,ymarker,zmarker
    REAL(KIND=CGREAL)    :: p1, p0

    CHARACTER*9          :: dragfile
    CHARACTER*25         :: indragfile

    WRITE (*,*) 'Calling subroutine drag_lift_membrane_partial_dynamics ...'

    IF (ndim == DIM_2D) THEN
         K_ini = 1
         K_end = 2
    ELSE
         K_ini = 2
         K_end = NZ-2
    END IF

!    DO iBody = nBody_solid+1, nBody_membrane
    Do iSection = 2,nSection

!---- Pressure Term and Shear Stress Term

      cxp  = zero
      cyp  = zero
      czp  = zero
      cmxp = zero
      cmyp = zero
      cmzp = zero

      cpw =  zero

      cxs  = zero
      cys  = zero
      czs  = zero
      cmxs = zero
      cmys = zero
      cmzs = zero

      hingex = xBodyMarker(iBody,hingemarker)
      hingey = yBodyMarker(iBody,hingemarker)
      hingez = zBodyMarker(iBody,hingemarker)

      moment_refx = hingex
      moment_refy = hingey
      moment_refz = hingez

      print *, '   moment_refx, moment_refy, moment_refz:', moment_refx, moment_refy, moment_refz

      CALL LOCATION2(hingex,hingey,hingez,ih,jh,kh)

!     Only jh is used for hovering case.

      DO k = K_ini, K_end
      DO j = 2, jh
      DO i = 2, nx-2

         distMin = 1.0E8_CGREAL

         DO md = 1, SectionMarker(iBody,iSection)
            m = DynamicMarker(iBody,iSection,mD)
            xmarker = xBodyMarker(iBody,m)
            ymarker = yBodyMarker(iBody,m)
            zmarker = zBodyMarker(iBody,m)
            dist = (xc(i) - xBodyMarker(iBody,m))**2 &
                 +(yc(j) - yBodyMarker(iBody,m))**2 &
                 +(zc(k) - zBodyMarker(iBody,m))**2

            IF ( dist <= distMin ) THEN
                  distMin = dist
                  mG  = m
            ENDIF
         ENDDO ! md

         IF ( bodyNum(i-1,j  ,k  ) == iBody .OR. &
             bodyNum(i+1,j  ,k  ) == iBody .OR. &
             bodyNum(i  ,j+1,k  ) == iBody .OR. &
             bodyNum(i  ,j-1,k  ) == iBody .OR. &
             bodyNum(i  ,j  ,k+1) == iBody .OR. &
             bodyNum(i  ,j  ,k-1) == iBody      ) THEN


            xp = -p(i,j,k)*ium(i,j,k)*dy(j)*dz(k) &
	         +p(i,j,k)*iup(i,j,k)*dy(j)*dz(k)

	    yp = -p(i,j,k)*jum(i,j,k)*dx(i)*dz(k) &
	         +p(i,j,k)*jup(i,j,k)*dx(i)*dz(k)

	    zp = -p(i,j,k)*kum(i,j,k)*dx(i)*dy(j) &
	         +p(i,j,k)*kup(i,j,k)*dx(i)*dy(j)

            cxp = cxp + xp
	    cyp = cyp + yp
	    czp = czp + zp

            cmxp = cmxp + ( yc(j) - moment_refy )*zp
            cmyp = cmyp - ( xc(i) - moment_refx )*zp
            cmzp = cmzp + ( xc(i) - moment_refx )*yp &
                        - ( yc(j) - moment_refy )*xp

            IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
               cmxp = cmxp - ( zc(k) - moment_refz )*yp
               cmyp = cmyp + ( zc(k) - moment_refz )*xp
            ENDIF

            amx = dxinv(i)*ium(i,j,k)*2.0_CGREAL
            apx = dxinv(i)*iup(i,j,k)*2.0_CGREAL

            amy = dyinv(j)*jum(i,j,k)*2.0_CGREAL
            apy = dyinv(j)*jup(i,j,k)*2.0_CGREAL


            amz = dzinv(k)*kum(i,j,k)*2.0_CGREAL
            apz = dzinv(k)*kup(i,j,k)*2.0_CGREAL

            xs = reinv * dx(i)* dz(k) *                              &
                    ( amy*(    u(i,j,k) -bcyu(i,j,k) )   &
                    - apy*( bcyu(i,j,k) -   u(i,j,k) ) ) &
                    +reinv * dx(i)* dy(j) *                              &
                    ( amz*(    u(i,j,k) -bczu(i,j,k) )   &
                    - apz*( bczu(i,j,k) -   u(i,j,k) ) )

            ys = reinv * dy(j)* dz(k) *                              &
                    ( amx*(    v(i,j,k) -bcxv(i,j,k) )   &
                    - apx*( bcxv(i,j,k) -   v(i,j,k) ) ) &
                    +reinv * dx(i)* dy(j) *                              &
                    ( amz*(    v(i,j,k) -bczv(i,j,k) )   &
                    - apz*( bczv(i,j,k) -   v(i,j,k) ) )

            zs = reinv * dy(j)* dz(k) *                              &
                    ( amx*(    w(i,j,k) -bcxw(i,j,k) )   &
                    - apx*( bcxw(i,j,k) -   w(i,j,k) ) ) &
                    +reinv * dx(i)* dz(k) *                              &
                    ( amy*(    w(i,j,k) -bcyw(i,j,k) )   &
                    - apy*( bcyw(i,j,k) -   w(i,j,k) ) )

            cxs = cxs + xs
            cys = cys + ys
            czs = czs + zs

            cpw = cpw + (xp+xs)*uBodyMarker(ibody,mG) + &
                        (yp+ys)*vBodyMarker(ibody,mG) + &
                        (zp+zs)*wBodyMarker(ibody,mG)

            cmxs = cmxs + ( yc(j) - moment_refy )*zs
            cmys = cmys - ( xc(i) - moment_refx )*zs
            cmzs = cmzs + ( xc(i) - moment_refx )*ys &
                        - ( yc(j) - moment_refy )*xs

            IF ( canonical_body_type(iBody) > GENERAL_CYLINDER ) THEN
               cmxs = cmxs - ( zc(k) - moment_refz )*ys
               cmys = cmys + ( zc(k) - moment_refz )*xs
            ENDIF

         ENDIF ! bodyNum

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

      scx(iSection)  = cxp  + cxs
      scy(iSection)  = cyp  + cys
      scz(iSection)  = czp  + czs
      scmx(iSection) = cmxp + cmxs
      scmy(iSection) = cmyp + cmys
      scmz(iSection) = cmzp + cmzs

      print *, '   cxp,cxs =',cxp,cxs
      print *, '   scx(iSection)=',scx(iSection)

      print *, '   cyp,cys =',cyp,cys
      print *, '   scy(iSection)=',scy(iSection)

      print *, '   cmzp,cmzs =',cmzp,cmzs
      print *, '   scmz(iSection)=',scmz(iSection)

!--- Construct Components

      cxp = 2.0_CGREAL*cxp
      cyp = 2.0_CGREAL*cyp
      czp = 2.0_CGREAL*czp
      cxs = 2.0_CGREAL*cxs
      cys = 2.0_CGREAL*cys
      czs = 2.0_CGREAL*czs

      cpw = 2.0_CGREAL*cpw

      cmxp= 2.0_CGREAL*cmxp
      cmxs= 2.0_CGREAL*cmxs
      cmyp= 2.0_CGREAL*cmyp
      cmys= 2.0_CGREAL*cmys
      cmzp= 2.0_CGREAL*cmzp
      cmzs= 2.0_CGREAL*cmzs

      IF ( ndim == DIM_2D ) THEN

        cxp = cxp/zout
        cyp = cyp/zout
        cxs = cxs/zout
        cys = cys/zout

        cmzp = cmzp/zout
        cmzs = cmzs/zout

      ENDIF

      cx  = cxp  + cxs         ! for non-cylinders these are
      cy  = cyp  + cys         ! raw
      cz  = czp  + czs         ! forces
      cmx = cmxp + cmxs        ! and moments
      cmy = cmyp + cmys        ! need to non-dmensionalize them
      cmz = cmzp + cmzs        ! during post-processing

!--- Write File

      IF (Converged_FSI(iBody)) THEN
         WRITE(ifuDragOutHinge+iBody-1,131) time,cxp,cxs,cx,cyp,cys,cy,czp,czs,cz &
                                        ,cmxp,cmxs,cmx,cmyp,cmys,cmy,cmzp,cmzs,cmz,cpw,aero_moment,grav_moment
      ENDIF

   END DO ! ibody

131 FORMAT(22(3x,1PE12.5))

   END SUBROUTINE drag_lift_membrane_partial_dynamics



   SUBROUTINE LOCATION2(xp,yp,zp,i1,j1,k1)

!  Find out i1,i2 in which x(i1)<xp<x(i2); similarly for y and z
!  Return back i1,j1,k1

   USE flow_arrays
   USE flow_parameters
   USE pressure_arrays
   USE grid_arrays
   USE boundary_arrays
   USE unstructured_surface_arrays

   implicit none

   REAL(KIND=CGREAL) :: xp,yp,zp
   integer, parameter :: NP=8
   double precision :: xe(NP),ye(NP)

   REAL(KIND=CGREAL) :: x1,x2,y1,y2
   double precision :: xdebug,ydebug,zdebug

   INTEGER :: nsearch, maxsearch
   integer :: i1,i2,j1,j2,k1,k2
   integer :: debug

   maxsearch = 1000

   debug = 0
   if (debug ==1) then
      xdebug = 12.9263721552879
      ydebug = 15.6456500000000
      zdebug = 1.000525357730557E-003
   endif

   i1 = 2
   i2 = nx-1

   nsearch = 0
   do while ( (i2-i1)>1.1 )
      if (nsearch >maxsearch) then
         write (*,*) 'max number of searching has reached in x'
         write (*,*) 'xp = ', xp
         stop
      endif
      if (xp>=xc(i1) .and. xp<xc((i1+i2)/2)) then
         i2 = (i1+i2)/2
      else if ( xp>=xc((i1+i2)/2) .and. xp<xc(i2) ) then
         i1 = (i1+i2)/2
      endif
      nsearch = nsearch+1
   enddo

   xe(1) = xc(i1)
   xe(2) = xc(i2)
   xe(3) = xe(2)
   xe(4) = xe(1)


   if (debug ==1 .and. abs(xp-xdebug)<1e-4 .and. abs(yp-ydebug)<1e-4) then
      print *, 'i1,i2,xe(1),xe(2):',i1,i2,xe(1),xe(2)
   endif

   j1 = 2
   j2 = ny-1

   nsearch = 0
   do while ( (j2-j1)>1.1 )
      if (nsearch >maxsearch) then
         write (*,*) 'max number of searching has reached in y'
         write (*,*) 'yp = ', yp
         stop
      endif
      if (yp>=yc(j1) .and. yp<yc((j1+j2)/2)) then
         j2 = (j1+j2)/2
      else if (yp>=yc((j1+j2)/2) .and. yp<yc(j2)) then
         j1 = (j1+j2)/2
      endif
      nsearch = nsearch+1
   enddo

   ye(1) = yc(j1)
   ye(2) = ye(1)
   ye(3) = yc(j2)
   ye(4) = ye(3)

   k1 = 2
   k2 = nz-1

   nsearch = 0
   do while ( (k2-k1)>1.1 )
      if (nsearch >maxsearch) then
         write (*,*) 'max number of searching has reached in z'
         write (*,*) 'zp = ', zp
         stop
      endif
      if (zp>=zc(k1) .and. zp<zc((k1+k2)/2)) then
         k2 = (k1+k2)/2
      else if (zp>=zc((k1+k2)/2) .and. zp<zc(k2)) then
         k1 = (k1+k2)/2
      endif

      nsearch = nsearch+1

   enddo

   if (debug ==1 .and. abs(xp-xdebug)<1e-4 .and. abs(yp-ydebug)<1e-4 .and. &
       abs(zp-zdebug)<1e-4 ) then
      print *, 'i1,i2,xe(1),xe(3):',j1,j2,ye(1),ye(3)
   endif

   END SUBROUTINE LOCATION2
!===========================================================================

!===========================================================================
   subroutine moment_rot(ibody)
! added by Geng Liu

   USE global_parameters
   USE flow_parameters
   USE usr_module
   USE body_dynamics
   USE boundary_arrays

   integer :: iBody,m

   REAL(KIND=CGREAL) :: vec1(3),quat_inv(4)

   quat_inv(1)=quat_prev(1)
   do m=2,4
     quat_inv(m)=-quat_prev(m)
   enddo

   vec1(1) = scmx(ibody)
   vec1(2) = scmy(ibody)
   vec1(3) = scmz(ibody)
   call quaternion_rotation(vec1,quat_inv)
   scmx(ibody)=vec1(1)
   scmy(ibody)=vec1(2)
   scmz(ibody)=vec1(3)

   endsubroutine moment_rot
!===========================================================================

