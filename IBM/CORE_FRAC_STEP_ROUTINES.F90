!-------------------------------------------------
! Subroutines :  face_vel() ;  correct_vel()
!-------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE face_vel()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays
    USE multiuse_arrays
    USE pressure_arrays
    USE solver_arrays
    use hybrid_cell_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,iBody,n
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
    REAL(KIND=CGREAL)    :: pgxe,pgxw,pgyn,pgys,pgzf,pgzb
    REAL(KIND=CGREAL)    :: ghostY, ghostN

    real(cgreal) :: dist,distx,disty,distz,distb,ratiox,ratioy,ratioz
    real(cgreal) :: xBIbeta,yBIbeta,zBIbeta
    real(cgreal) :: xBIalpha1,yBIalpha1,zBIalpha1
    real(cgreal) :: xBIalpha2,yBIalpha2,zBIalpha2
    real(cgreal) :: xBIalpha3,yBIalpha3,zBIalpha3
    integer :: h_mark(3)

! compute face velocities

    IF (frac_step_type == VAN_KAN) THEN

      DO k = 1,nzc
      DO j = 1,nyc
      DO i = 1,nxc

        pe = ( fx(i+1)*p(i+1,j,k) + (oned-fx(i+1))*p(i,j,k)   )*(1-iup(i,j,k)) &
             + p(i,j,k)*iup(i,j,k)

        pw = ( fx(i)  *p(i,j,k)   + (oned-fx(i))  *p(i-1,j,k) )*(1-ium(i,j,k)) &
             + p(i,j,k)*ium(i,j,k)

        pn = ( fy(j+1)*p(i,j+1,k) + (oned-fy(j+1))*p(i,j,k)   )*(1-jup(i,j,k)) &
             + p(i,j,k)*jup(i,j,k)

        ps = ( fy(j)  *p(i,j,k)   + (oned-fy(j))  *p(i,j-1,k) )*(1-jum(i,j,k)) &
             + p(i,j,k)*jum(i,j,k)

        pf = ( fz(k+1)*p(i,j,k+1) + (oned-fz(k+1))*p(i,j,k)   )*(1-kup(i,j,k)) &
             + p(i,j,k)*kup(i,j,k)

        pb = ( fz(k)  *p(i,j,k)   + (oned-fz(k))  *p(i,j,k-1) )*(1-kum(i,j,k)) &
             + p(i,j,k)*kum(i,j,k)

        pgx= (pe-pw)*dxinv(i)
        pgy= (pn-ps)*dyinv(j)
        pgz= (pf-pb)*dzinv(k)

        uTilde(i,j,k) = u(i,j,k) + dt*pgx*REAL(1-iblank(i,j,k),KIND=CGREAL)
        vTilde(i,j,k) = v(i,j,k) + dt*pgy*REAL(1-iblank(i,j,k),KIND=CGREAL)
        wTilde(i,j,k) = w(i,j,k) + dt*pgz*REAL(1-iblank(i,j,k),KIND=CGREAL)

      ENDDO
      ENDDO
      ENDDO

    ELSE

      uTilde = u
      vTilde = v
      wTilde = w

    ENDIF ! frac_step_type

    DO k = 1,nzc
    DO j = 1,nyc
      DO i = 1,nxc
        face_u(i+1,j,k) = (        fx(i+1) *uTilde(i+1,j,k)                   &
                           + (oned-fx(i+1))*uTilde(i,j,k)   )*(1-iup(i,j,k))  &
                         + bcxu(i,j,k)*iup(i,j,k)

        face_u(i,j,k)   = (        fx(i) *uTilde(i,j,k)                       &
                           + (oned-fx(i))*uTilde(i-1,j,k)   )*(1-ium(i,j,k))  &
                         + bcxu(i,j,k)*ium(i,j,k)

        face_v(i,j+1,k) = (        fy(j+1) *vTilde(i,j+1,k)                   &
                           + (oned-fy(j+1))*vTilde(i,j,k)   )*(1-jup(i,j,k))  &
                         + bcyv(i,j,k)*jup(i,j,k)

        face_v(i,j,k)   = (        fy(j) *vTilde(i,j,k)                       &
                           + (oned-fy(j))*vTilde(i,j-1,k)   )*(1-jum(i,j,k))  &
                         + bcyv(i,j,k)*jum(i,j,k)

        face_w(i,j,k+1) = (        fz(k+1) *wTilde(i,j,k+1)                   &
                           + (oned-fz(k+1))*wTilde(i,j,k)   )*(1-kup(i,j,k))  &
                         + bczw(i,j,k)*kup(i,j,k)

        face_w(i,j,k)   = (        fz(k) *wTilde(i,j,k)                       &
                           + (oned-fz(k))*wTilde(i,j,k-1)   )*(1-kum(i,j,k))  &
                         + bczw(i,j,k)*kum(i,j,k)

        ENDDO
    ENDDO
    ENDDO

    IF (frac_step_type == VAN_KAN) THEN

! add pressure gradient to face velocities
      DO k = 1,nzc
      DO j = 1,nyc
        DO i = 1,nxc
            pgxw       = (p(i,j,k)  -p(i-1,j,k))*dxcinv(i)
            pgxe       = (p(i+1,j,k)-p(i,j,k)  )*dxcinv(i+1)
            face1(i)   = face_u(i,j,k)   - dt*pgxw*(1-ium(i,j,k))
            face2(i)   = face_u(i+1,j,k) - dt*pgxe*(1-iup(i,j,k))
        ENDDO
        DO i = 1,nxc
            face_u(i,j,k)     = face1(i)
            face_u(i+1,j,k)   = face2(i)
        ENDDO
        DO i = 1,nxc
            face_u(i,j,k)     = face_u(i,j,k)  *(1 - ium(i,j,k)) &
                              + bcxu(i,j,k)*ium(i,j,k)
            face_u(i+1,j,k)   = face_u(i+1,j,k)*(1 - iup(i,j,k)) &
                              + bcxu(i,j,k)*iup(i,j,k)
        ENDDO
      ENDDO
      ENDDO

      DO k = 1,nzc
      DO i = 1,nxc
        DO j = 1,nyc
            pgys       = (p(i,j,k)  -p(i,j-1,k))*dycinv(j)
            pgyn       = (p(i,j+1,k)-p(i,j,k)  )*dycinv(j+1)
            face1(j)   = face_v(i,j,k)   - dt*pgys*(1-jum(i,j,k))  
            face2(j)   = face_v(i,j+1,k) - dt*pgyn*(1-jup(i,j,k)) 
        ENDDO
        DO j = 1,nyc
            face_v(i,j,k)   = face1(j)
            face_v(i,j+1,k) = face2(j)
        ENDDO
        DO j = 1,nyc
            face_v(i,j,k)     = face_v(i,j,k)  *(1 - jum(i,j,k)) &
                              + bcyv(i,j,k)*jum(i,j,k)
            face_v(i,j+1,k)   = face_v(i,j+1,k)*(1 - jup(i,j,k)) &
                              + bcyv(i,j,k)*jup(i,j,k)
        ENDDO
      ENDDO
      ENDDO

      DO j = 1,nyc
      DO i = 1,nxc
        DO k = 1,nzc
            pgzb     = (p(i,j,k)  -p(i,j,k-1))*dzcinv(k)
            pgzf     = (p(i,j,k+1)-p(i,j,k)  )*dzcinv(k+1)
            face1(k) = face_w(i,j,k)   - dt*pgzb*(1-kum(i,j,k))  
            face2(k) = face_w(i,j,k+1) - dt*pgzf*(1-kup(i,j,k))
        ENDDO
        DO k = 1,nzc
            face_w(i,j,k)   = face1(k)
            face_w(i,j,k+1) = face2(k)
        ENDDO
        DO k = 1,nzc
            face_w(i,j,k)     = face_w(i,j,k)  *(1 - kum(i,j,k)) &
                              + bczw(i,j,k)*kum(i,j,k)
            face_w(i,j,k+1)   = face_w(i,j,k+1)*(1 - kup(i,j,k)) &
                              + bczw(i,j,k)*kup(i,j,k)
        ENDDO
      ENDDO
      ENDDO

    ENDIF ! frac_step_type

    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
       face_u(i,j,k)     = face_u(i,j,k)  *(1 - ium(i,j,k)) &
                         + bcxu(i,j,k)*ium(i,j,k)
       face_u(i+1,j,k)   = face_u(i+1,j,k)*(1 - iup(i,j,k)) &
                         + bcxu(i,j,k)*iup(i,j,k)
       face_v(i,j,k)     = face_v(i,j,k)  *(1 - jum(i,j,k)) &
                         + bcyv(i,j,k)*jum(i,j,k)
       face_v(i,j+1,k)   = face_v(i,j+1,k)*(1 - jup(i,j,k)) &
                         + bcyv(i,j,k)*jup(i,j,k)
       face_w(i,j,k)     = face_w(i,j,k)  *(1 - kum(i,j,k)) &
                         + bczw(i,j,k)*kum(i,j,k)
       face_w(i,j,k+1)   = face_w(i,j,k+1)*(1 - kup(i,j,k)) &
                         + bczw(i,j,k)*kup(i,j,k)
    ENDDO
    ENDDO
    ENDDO

! for GCM, face velocities at boundary-ghost faces are not adjusted
! for the compact/non-compact pressure gradient 

    IF ( frac_step_type == VAN_KAN .AND. &
         boundary_formulation == GCM_METHOD ) THEN 
  
       DO k = 1,nzc
       DO j = 1,nyc
       DO i = 1,nxc

          ghostY = REAL(ghostCellMark(i,j,k),KIND=CGREAL)
          ghostN = oned - ghostY

          face_u(i+1,j,k) = ghostN*face_u(i+1,j,k) &
                           +ghostY*(fx(i+1)*u(i+1,j,k) + (oned-fx(i+1))*u(i,j,k)   )

          face_u(i,j,k)   = ghostN*face_u(i,j,k)  &
                           +ghostY*(fx(i) *u(i,j,k)    + (oned-fx(i))*u(i-1,j,k)   )

          face_v(i,j+1,k) = ghostN*face_v(i,j+1,k) &
                           +ghostY*(fy(j+1)*v(i,j+1,k) + (oned-fy(j+1))*v(i,j,k)   )

          face_v(i,j,k)   = ghostN*face_v(i,j,k)  &
                           +ghostY*(fy(j) *v(i,j,k)    + (oned-fy(j))*v(i,j-1,k)   )

          face_w(i,j,k+1) = ghostN*face_w(i,j,k+1) &
                           +ghostY*(fz(k+1)*w(i,j,k+1) + (oned-fz(k+1))*w(i,j,k)   )

          face_w(i,j,k)   = ghostN*face_w(i,j,k)  &
                           +ghostY*(fz(k) *w(i,j,k)    + (oned-fz(k))*w(i,j,k-1)   )


       ENDDO
       ENDDO
       ENDDO

    ENDIF

! Zeroing out face velocities for solid cells.
! Also for "dead" faces of ghost cells
    IF (boundary_formulation == GCM_METHOD) CALL GCM_set_face_vel_body()

   END SUBROUTINE face_vel   
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE correct_vel()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays

    IMPLICIT NONE
    
    INTEGER :: i,j,k,nmax
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
    REAL(KIND=CGREAL)    :: pgxe,pgxw,pgyn,pgys,pgzf,pgzb
    REAL(KIND=CGREAL)    :: AAe, AAw, AAn, AAs, AAf, AAb, AAg, Ag
    REAL(KIND=CGREAL)    :: CVF

    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
       CALL enforce_p_periodic(pPrime)
    END IF

! correct nodal velocities
! need to account for boundaries...
    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc

       IF (bcx1 .EQ. BC_TYPE_PERIODIC .AND. &
           bcx2 .EQ. BC_TYPE_PERIODIC ) THEN

        pe = fx(i+1)*pPrime(i+1,j,k) + (oned-fx(i+1))*pPrime(i,j,k)
 
        pw = fx(i)  *pPrime(i,j,k)   + (oned-fx(i))  *pPrime(i-1,j,k)

       ELSE
        pe = ( fx(i+1)*pPrime(i+1,j,k) + (oned-fx(i+1))*pPrime(i,j,k)   )  &
           *(1-iup(i,j,k)) &
           + pPrime(i,j,k)*iup(i,j,k)

        pw = ( fx(i)  *pPrime(i,j,k)   + (oned-fx(i))  *pPrime(i-1,j,k) )  &
           *(1-ium(i,j,k)) &
           + pPrime(i,j,k)*ium(i,j,k)

       END IF

     IF (bcy1 .EQ. BC_TYPE_PERIODIC .AND. &
          bcy2 .EQ. BC_TYPE_PERIODIC) THEN

      pn = fy(j+1)*pPrime(i,j+1,k) + (oned-fy(j+1))*pPrime(i,j,k)
 
      ps = fy(j)  *pPrime(i,j,k)   + (oned-fy(j))  *pPrime(i,j-1,k)

     ELSE
      pn = ( fy(j+1)*pPrime(i,j+1,k) + (oned-fy(j+1))*pPrime(i,j,k)   )  &
          *(1-jup(i,j,k)) &
           + pPrime(i,j,k)*jup(i,j,k)

      ps = ( fy(j)  *pPrime(i,j,k)   + (oned-fy(j))  *pPrime(i,j-1,k) )  &
          *(1-jum(i,j,k)) &
           + pPrime(i,j,k)*jum(i,j,k)

     END IF


     IF (bcz1 .EQ. BC_TYPE_PERIODIC .AND. &
         bcz2 .EQ. BC_TYPE_PERIODIC ) THEN

      pf =  fz(k+1)*pPrime(i,j,k+1) + (oned-fz(k+1))*pPrime(i,j,k)
      pb =  fz(k)  *pPrime(i,j,k)   + (oned-fz(k))  *pPrime(i,j,k-1)

     ELSE
      pf = ( fz(k+1)*pPrime(i,j,k+1) + (oned-fz(k+1))*pPrime(i,j,k)   )  &
          *(1-kup(i,j,k)) &
           + pPrime(i,j,k)*kup(i,j,k)

      pb = ( fz(k)  *pPrime(i,j,k)   + (oned-fz(k))  *pPrime(i,j,k-1) )  &
          *(1-kum(i,j,k)) &
           + pPrime(i,j,k)*kum(i,j,k)

     END IF

      IF (cure_pressure_oscillations .AND. iblank(i,j,k)==0 .AND. ivc(i,j,k)>0) THEN
        AAe=face(cell(ivc(i,j,k))%F_ip)%a
        AAw=face(cell(ivc(i,j,k))%F_im)%a
        AAn=face(cell(ivc(i,j,k))%F_jp)%a
        AAs=face(cell(ivc(i,j,k))%F_jm)%a
        AAf=face(cell(ivc(i,j,k))%F_kp)%a
        AAb=face(cell(ivc(i,j,k))%F_km)%a
        AAg=face(cell(ivc(i,j,k))%F_slice)%a

        pe = pe * AAe
        pw = pw * AAw
        pn = pn * AAn
        ps = ps * AAs
        pf = pf * AAf
        pb = pb * AAb
      ENDIF

      pgx= (pe-pw)*dxinv(i)
      pgy= (pn-ps)*dyinv(j)
      pgz= (pf-pb)*dzinv(k)

      IF (cure_pressure_oscillations .AND. iblank(i,j,k)==0 .AND. ivc(i,j,k)>0) THEN
        CVF=cell(ivc(i,j,k))%volumn
        pgx=pgx+pPrime(i,j,k)*cell(ivc(i,j,k))%slice_normal(1)*AAg
        pgy=pgy+pPrime(i,j,k)*cell(ivc(i,j,k))%slice_normal(2)*AAg
        pgz=pgz+pPrime(i,j,k)*cell(ivc(i,j,k))%slice_normal(3)*AAg
        u(i,j,k) = u(i,j,k) - dt*pgx*REAL(1-iblank(i,j,k),KIND=CGREAL)/CVF
        v(i,j,k) = v(i,j,k) - dt*pgy*REAL(1-iblank(i,j,k),KIND=CGREAL)/CVF
        w(i,j,k) = w(i,j,k) - dt*pgz*REAL(1-iblank(i,j,k),KIND=CGREAL)/CVF
      ELSE
        u(i,j,k) = u(i,j,k) - dt*pgx*REAL(1-iblank(i,j,k),KIND=CGREAL)
        v(i,j,k) = v(i,j,k) - dt*pgy*REAL(1-iblank(i,j,k),KIND=CGREAL)
        w(i,j,k) = w(i,j,k) - dt*pgz*REAL(1-iblank(i,j,k),KIND=CGREAL)
      ENDIF
    ENDDO
    ENDDO
    ENDDO

    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. & 
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. & 
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
       CALL enforce_u_periodic
    END IF

! update value of velocity at ghost points through interpolation
!   WRITE(*,*)'entering GCM_vel_set_bc_internal'
    IF (boundary_formulation == GCM_METHOD) CALL GCM_vel_set_bc_internal()

! correct face velocities
    DO k = 1,nzc
    DO j = 1,nyc
      DO i = 1,nxc
        pgxw       = (pPrime(i,j,k)  -pPrime(i-1,j,k))*dxcinv(i)
        pgxe       = (pPrime(i+1,j,k)-pPrime(i,j,k)  )*dxcinv(i+1)
        face1(i)   = face_u(i,j,k)   - dt*pgxw*(1-ium(i,j,k))
        face2(i)   = face_u(i+1,j,k) - dt*pgxe*(1-iup(i,j,k))
      ENDDO
      DO i = 1,nxc
          face_u(i,j,k)     = face1(i)
          face_u(i+1,j,k)   = face2(i)
      ENDDO
      DO i = 1,nxc
          face_u(i,j,k)     = face_u(i,j,k)  *(1 - ium(i,j,k)) &
                            + bcxu(i,j,k)*ium(i,j,k)
          face_u(i+1,j,k)   = face_u(i+1,j,k)*(1 - iup(i,j,k)) &
                            + bcxu(i,j,k)*iup(i,j,k)
      ENDDO
    ENDDO
    ENDDO

    DO k = 1,nzc
    DO i = 1,nxc
      DO j = 1,nyc
        pgys       = (pPrime(i,j,k)  -pPrime(i,j-1,k))*dycinv(j)
        pgyn       = (pPrime(i,j+1,k)-pPrime(i,j,k)  )*dycinv(j+1)
        face1(j)   = face_v(i,j,k)   - dt*pgys*(1-jum(i,j,k))  
        face2(j)   = face_v(i,j+1,k) - dt*pgyn*(1-jup(i,j,k)) 
      ENDDO
      DO j = 1,nyc
          face_v(i,j,k)   = face1(j)
          face_v(i,j+1,k) = face2(j)
      ENDDO
      DO j = 1,nyc
          face_v(i,j,k)     = face_v(i,j,k)  *(1 - jum(i,j,k)) &
                            + bcyv(i,j,k)*jum(i,j,k)
          face_v(i,j+1,k)   = face_v(i,j+1,k)*(1 - jup(i,j,k)) &
                            + bcyv(i,j,k)*jup(i,j,k)
      ENDDO
    ENDDO
    ENDDO

    DO j = 1,nyc
    DO i = 1,nxc
      DO k = 1,nzc
        pgzb     = (pPrime(i,j,k)  -pPrime(i,j,k-1))*dzcinv(k)
        pgzf     = (pPrime(i,j,k+1)-pPrime(i,j,k)  )*dzcinv(k+1)
        face1(k) = face_w(i,j,k)   - dt*pgzb*(1-kum(i,j,k))  
        face2(k) = face_w(i,j,k+1) - dt*pgzf*(1-kup(i,j,k))
      ENDDO
      DO k = 1,nzc
          face_w(i,j,k)   = face1(k)
          face_w(i,j,k+1) = face2(k)
      ENDDO
      DO k = 1,nzc
          face_w(i,j,k)     = face_w(i,j,k)  *(1 - kum(i,j,k)) &
                            + bczw(i,j,k)*kum(i,j,k)
          face_w(i,j,k+1)   = face_w(i,j,k+1)*(1 - kup(i,j,k)) &
                            + bczw(i,j,k)*kup(i,j,k)
      ENDDO
    ENDDO
    ENDDO


! set face vel. to zero for solid cells.
! also for ghost-ghost cell faces

    IF (boundary_formulation == GCM_METHOD) CALL GCM_set_face_vel_body()

   END SUBROUTINE correct_vel
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE correct_face_vel()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays

    IMPLICIT NONE
    
    INTEGER :: i,j,k,nmax
    INTEGER          :: closest_marker,body_dist_min
    REAL(KIND=CGREAL), DIMENSION(3)  :: Ub, vec
	  INTEGER :: fc

! correct nodal velocities
! need to account for boundaries...
    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      IF (ivc(i,j,k)>0) THEN
		  call find_closest_element(face(cell(ivc(i,j,k))%F_slice)%centroid, body_dist_min, closest_marker)

		  Ub(1)=uBodyMarker(body_dist_min,closest_marker)
		  Ub(2)=vBodyMarker(body_dist_min,closest_marker)
		  Ub(3)=wBodyMarker(body_dist_min,closest_marker)
!!!******************
		  fc=cell(ivc(i,j,k))%F_ip
		  IF (fc>0) THEN
!          IF (iup(i,j,k)>0) THEN
			  vec=face(fc)%centroid
  		  face(fc)%vel=taylor_series_matching(i,j,k,Ub(1), vec, u)
!         ELSE 
!            face(fc)%vel=face_u(i+1,j,k)
!		  ENDIF
      ENDIF
!!!******************
		  fc=cell(ivc(i,j,k))%F_im
      IF (fc>0) THEN
!          IF (ium(i,j,k)>0) THEN
			  vec=face(fc)%centroid
			  face(fc)%vel=taylor_series_matching(i,j,k,Ub(1), vec, u)
!		  ELSE
!			face(fc)%vel=face_u(i,j,k)
!		  ENDIF
      ENDIF
!!!******************
		  fc=cell(ivc(i,j,k))%F_jp
      IF (fc>0) THEN
!          IF (jup(i,j,k)>0) THEN
			  vec=face(fc)%centroid
			  face(fc)%vel=taylor_series_matching(i,j,k,Ub(2), vec, v)
!		  ELSE
!			face(fc)%vel=face_v(i,j+1,k)
!		  ENDIF
      ENDIF
!!!******************
		  fc=cell(ivc(i,j,k))%F_jm
      IF (fc>0) THEN
!          IF (jum(i,j,k)>0) THEN
			  vec=face(fc)%centroid
			  face(fc)%vel=taylor_series_matching(i,j,k,Ub(2), vec, v)
!          ELSE
!	        face(fc)%vel=face_v(i,j,k)
!	      ENDIF
  		ENDIF
!!!******************

		  IF (ndim==DIM_3D) THEN
		    fc=cell(ivc(i,j,k))%F_kp
        IF (fc>0) THEN
          IF (kup(i,j,k)>0) THEN
			      vec=face(fc)%centroid
			      face(fc)%vel=taylor_series_matching(i,j,k,Ub(3), vec, w)
		      ELSE
			      face(fc)%vel=face_w(i,j,k+1)
		      ENDIF
  		  ENDIF
!!!******************
		    fc=cell(ivc(i,j,k))%F_km
        IF (fc>0) THEN
          IF (kum(i,j,k)>0) THEN
			      vec=face(fc)%centroid
			      face(fc)%vel=taylor_series_matching(i,j,k,Ub(3), vec, w)
		      ELSE
			      face(fc)%vel=face_w(i,j,k)
		      ENDIF
		    ENDIF
		  ENDIF
!!!******************
				  
      ENDIF
    ENDDO
    ENDDO
    ENDDO


! set face vel. to zero for solid cells.
! also for ghost-ghost cell faces

CONTAINS
!-------------------------------------------------------------------------------

  FUNCTION taylor_series_matching(i, j, k, Ub, Vec,uv)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
    USE body_dynamics

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:,0:,0:), INTENT(IN)  :: uv
    INTEGER          :: m, body_dist_min1, body_dist_min2, body_dist_min3
    INTEGER          :: i, j, k
    INTEGER          :: i1, j1, k1
    INTEGER          :: i2, j2, k2
    INTEGER          :: i3, j3, k3
    INTEGER          :: i4, j4, k4
    INTEGER, DIMENSION(1) :: iclose
    INTEGER, PARAMETER :: wo=2, ws=wo*2+1, ws2=ws*ws, ws3=ws*ws*ws
    REAL(KIND=CGREAL) :: dist_min,dist_min_body,taylor_series_matching
    REAL(KIND=CGREAL) :: distX, distY, distZ
    REAL(KIND=CGREAL) :: dist_sqr(ws3)
    REAL(KIND=CGREAL) :: dx1, dx2, dx3, dx4
    REAL(KIND=CGREAL) :: dy1, dy2, dy3, dy4
    REAL(KIND=CGREAL) :: dz1, dz2, dz3, dz4
    REAL(KIND=CGREAL) :: x1, x2, x3, x4
    REAL(KIND=CGREAL) :: y1, y2, y3, y4
    REAL(KIND=CGREAL) :: z1, z2, z3, z4
    REAL(KIND=CGREAL) :: c1, c2, c3, c4
    REAL(KIND=CGREAL) :: Ub
    REAL(KIND=CGREAL), DIMENSION(3)  :: Vec

    m=0
    dist_sqr=1.0E10_CGREAL
	  SELECT CASE(ndim)
	  CASE(DIM_2D)
	    DO j1=j-wo, j+wo
	    DO i1=i-wo, i+wo
		  m=m+1
		  IF (iblank(I1,J1,K)==0 .and. ivc(i1,j1,k)==0) THEN
		    distX = (vec(1)-xc(i1))
		    distY = (vec(2)-yc(j1))
		    dist_sqr(m) = distX**2 + distY**2
		  ENDIF
	    ENDDO
	    ENDDO
	  CASE(DIM_3D)
	    DO k1=k-wo, k+wo
	    DO j1=j-wo, j+wo
	    DO i1=i-wo, i+wo
		  m=m+1
		  IF (iblank(I1,J1,K1)==0 .and. ivc(i1,j1,k1)==0) THEN
		    distX = (vec(1)-xc(i1))
		    distY = (vec(2)-yc(j1))
		    distZ = (vec(3)-zc(k1))
		    dist_sqr(m) = distX**2 + distY**2 + distZ**2
		  ENDIF
	    ENDDO
	    ENDDO
	    ENDDO
	  END SELECT ! canonical_body_type
    iclose            = MINLOC(dist_sqr)
    body_dist_min1 = iclose(1)
    dist_sqr(body_dist_min1)=1.0E10_CGREAL
    iclose            = MINLOC(dist_sqr)
    body_dist_min2 = iclose(1)
	  dist_sqr(body_dist_min2)=1.0E10_CGREAL
	  iclose            = MINLOC(dist_sqr)
	  body_dist_min3 = iclose(1)
    IF (ndim==DIM_2D) THEN
	  j2=j+(body_dist_min1-1)/ws-wo
	  i2=i+mod(body_dist_min1-1,ws)-wo

	  j3=j+(body_dist_min2-1)/ws-wo
	  i3=i+mod(body_dist_min2-1,ws)-wo

	  j4=j+(body_dist_min3-1)/ws-wo
	  i4=i+mod(body_dist_min3-1,ws)-wo

	  dx1=face(cell(ivc(i,j,k))%F_slice)%centroid(1)-vec(1)
	  dy1=face(cell(ivc(i,j,k))%F_slice)%centroid(2)-vec(2)
	  dx2=xc(i2)-vec(1)
	  dy2=yc(j2)-vec(2)
	  dx3=xc(i3)-vec(1)
	  dy3=yc(j3)-vec(2)
	  
	  IF ((dx2*dy1 - dx3*dy1 - dx1*dy2 + dx3*dy2 + dx1*dy3 - dx2*dy3)==ZERO) THEN
	  dx3=xc(i4)-vec(1)
	  dy3=yc(j4)-vec(2)
	  ENDIF
c1=(dx3*dy2 - dx2*dy3)/(dx2*dy1 - dx3*dy1 - dx1*dy2 + dx3*dy2 + dx1*dy3 - dx2*dy3)

c2=(-dx3*dy1 + dx1*dy3)/(dx2*dy1 - dx3*dy1 - dx1*dy2 + dx3*dy2 + dx1*dy3 - dx2*dy3)

c3=(-dx2*dy1 + dx1*dy2)/(-dx2*dy1 + dx3*dy1 + dx1*dy2 - dx3*dy2 - dx1*dy3 + dx2*dy3)
     
    taylor_series_matching=c1*ub+c2*uv(i2,j2,k)+c3*uv(i3,j3,k)
    ELSE

    
	  k2=k+(body_dist_min1-1)/ws2-wo
	  j2=j+mod(body_dist_min1-1,ws2)/ws-wo
	  i2=i+mod(mod(body_dist_min1-1,ws2),ws)-wo

	  k3=k+(body_dist_min2-1)/ws2-wo
	  j3=j+mod(body_dist_min2-1,ws2)/ws-wo
	  i3=i+mod(mod(body_dist_min2-1,ws2),ws)-wo

	  k4=k+(body_dist_min3-1)/ws2-wo
	  j4=j+mod(body_dist_min3-1,ws2)/ws-wo
	  i4=i+mod(mod(body_dist_min3-1,ws2),ws)-wo
	
	  dx1=face(cell(ivc(i,j,k))%F_slice)%centroid(1)-vec(1)
	  dy1=face(cell(ivc(i,j,k))%F_slice)%centroid(2)-vec(2)
	  dz1=face(cell(ivc(i,j,k))%F_slice)%centroid(3)-vec(3)
	  dx2=xc(i2)-vec(1)
	  dy2=yc(j2)-vec(2)
	  dz2=zc(k2)-vec(3)
	  dx3=xc(i3)-vec(1)
	  dy3=yc(j3)-vec(2)
	  dz3=zc(k3)-vec(3)
	  dx4=xc(i4)-vec(1)
	  dy4=yc(j4)-vec(2)
	  dz4=zc(k4)-vec(3)
	  
c1=(-dx4*dy3*dz2 + dx3*dy4*dz2 + dx4*dy2*dz3 - dx2*dy4*dz3 - &
   dx3*dy2*dz4 + dx2*dy3*dz4)/(dx3*dy2*dz1 - dx4*dy2*dz1 - &
   dx2*dy3*dz1 + dx4*dy3*dz1 + dx2*dy4*dz1 - dx3*dy4*dz1 - &
   dx3*dy1*dz2 + dx4*dy1*dz2 + dx1*dy3*dz2 - dx4*dy3*dz2 - &
   dx1*dy4*dz2 + dx3*dy4*dz2 + dx2*dy1*dz3 - dx4*dy1*dz3 - &
   dx1*dy2*dz3 + dx4*dy2*dz3 + dx1*dy4*dz3 - dx2*dy4*dz3 - &
   dx2*dy1*dz4 + dx3*dy1*dz4 + dx1*dy2*dz4 - dx3*dy2*dz4 - &
   dx1*dy3*dz4 + dx2*dy3*dz4)
      
c2=(-dx4*dy3*dz1 + dx3*dy4*dz1 + dx4*dy1*dz3 - dx1*dy4*dz3 - &
   dx3*dy1*dz4 + dx1*dy3*dz4)/(-dx3*dy2*dz1 + dx4*dy2*dz1 + & 
   dx2*dy3*dz1 - dx4*dy3*dz1 - dx2*dy4*dz1 + dx3*dy4*dz1 + &
   dx3*dy1*dz2 - dx4*dy1*dz2 - dx1*dy3*dz2 + dx4*dy3*dz2 + &
   dx1*dy4*dz2 - dx3*dy4*dz2 - dx2*dy1*dz3 + dx4*dy1*dz3 + &
   dx1*dy2*dz3 - dx4*dy2*dz3 - dx1*dy4*dz3 + dx2*dy4*dz3 + &
   dx2*dy1*dz4 - dx3*dy1*dz4 - dx1*dy2*dz4 + dx3*dy2*dz4 + &
   dx1*dy3*dz4 - dx2*dy3*dz4)
   
c3=(-dx4*dy2*dz1 + dx2*dy4*dz1 + dx4*dy1*dz2 - dx1*dy4*dz2 - &
   dx2*dy1*dz4 + dx1*dy2*dz4)/(dx3*dy2*dz1 - dx4*dy2*dz1 - &
   dx2*dy3*dz1 + dx4*dy3*dz1 + dx2*dy4*dz1 - dx3*dy4*dz1 - &
   dx3*dy1*dz2 + dx4*dy1*dz2 + dx1*dy3*dz2 - dx4*dy3*dz2 - &
   dx1*dy4*dz2 + dx3*dy4*dz2 + dx2*dy1*dz3 - dx4*dy1*dz3 - &
   dx1*dy2*dz3 + dx4*dy2*dz3 + dx1*dy4*dz3 - dx2*dy4*dz3 - &
   dx2*dy1*dz4 + dx3*dy1*dz4 + dx1*dy2*dz4 - dx3*dy2*dz4 - &
   dx1*dy3*dz4 + dx2*dy3*dz4)

c4=(-dx3*dy2*dz1 + dx2*dy3*dz1 + dx3*dy1*dz2 - dx1*dy3*dz2 - &
   dx2*dy1*dz3 + dx1*dy2*dz3)/(-dx3*dy2*dz1 + dx4*dy2*dz1 + & 
   dx2*dy3*dz1 - dx4*dy3*dz1 - dx2*dy4*dz1 + dx3*dy4*dz1 + & 
   dx3*dy1*dz2 - dx4*dy1*dz2 - dx1*dy3*dz2 + dx4*dy3*dz2 + &
   dx1*dy4*dz2 - dx3*dy4*dz2 - dx2*dy1*dz3 + dx4*dy1*dz3 + &
   dx1*dy2*dz3 - dx4*dy2*dz3 - dx1*dy4*dz3 + dx2*dy4*dz3 + &
   dx2*dy1*dz4 - dx3*dy1*dz4 - dx1*dy2*dz4 + dx3*dy2*dz4 + &
   dx1*dy3*dz4 - dx2*dy3*dz4)
   
    taylor_series_matching=c1*ub+c2*uv(i2,j2,k2)+c3*uv(i3,j3,k3)+c4*uv(i4,j4,k4)
	ENDIF
   

  END FUNCTION taylor_series_matching
!-------------------------------------------------------------------------------
   END SUBROUTINE correct_face_vel
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE correct_cell_vel()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays

    IMPLICIT NONE
    
    INTEGER :: i,j,k,nmax
    INTEGER          :: closest_marker,body_dist_min
    REAL(KIND=CGREAL), DIMENSION(3)  :: Ub, vec
	  INTEGER :: fc

! correct nodal velocities
! need to account for boundaries...
    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      IF (ivc(i,j,k)>0) THEN
		    call find_closest_element(face(cell(ivc(i,j,k))%F_slice)%centroid, body_dist_min, closest_marker)

		    Ub(1)=uBodyMarker(body_dist_min,closest_marker)
		    Ub(2)=vBodyMarker(body_dist_min,closest_marker)
		    Ub(3)=wBodyMarker(body_dist_min,closest_marker)

		    vec=(/xc(i),yc(j),zc(k)/)
        call taylor_series_matching_cell(i,j,k,Ub, vec)
				  
      ENDIF
    ENDDO
    ENDDO
    ENDDO


! set face vel. to zero for solid cells.
! also for ghost-ghost cell faces

    CALL face_vel()

CONTAINS
!-------------------------------------------------------------------------------

  SUBROUTINE taylor_series_matching_cell(i, j, k, Ub, Vec)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
    USE body_dynamics

    IMPLICIT NONE

    INTEGER          :: m, body_dist_min1, body_dist_min2, body_dist_min3
    INTEGER          :: i, j, k
    INTEGER          :: i1, j1, k1
    INTEGER          :: i2, j2, k2
    INTEGER          :: i3, j3, k3
    INTEGER          :: i4, j4, k4
    INTEGER, DIMENSION(1) :: iclose
    INTEGER, PARAMETER :: wo=2, ws=wo*2+1, ws2=ws*ws, ws3=ws*ws*ws
    REAL(KIND=CGREAL) :: dist_min,dist_min_body
    REAL(KIND=CGREAL) :: distX, distY, distZ
    REAL(KIND=CGREAL) :: dist_sqr(ws3)
    REAL(KIND=CGREAL) :: dx1, dx2, dx3, dx4
    REAL(KIND=CGREAL) :: dy1, dy2, dy3, dy4
    REAL(KIND=CGREAL) :: dz1, dz2, dz3, dz4
    REAL(KIND=CGREAL) :: x1, x2, x3, x4
    REAL(KIND=CGREAL) :: y1, y2, y3, y4
    REAL(KIND=CGREAL) :: z1, z2, z3, z4
    REAL(KIND=CGREAL) :: c1, c2, c3, c4

    REAL(KIND=CGREAL), DIMENSION(3)  :: Ub, vec

    m=0
    dist_sqr=1.0E10_CGREAL
	  SELECT CASE(ndim)
	  CASE(DIM_2D)
	    DO j1=j-wo, j+wo
	    DO i1=i-wo, i+wo
		  m=m+1
		  IF (iblank(I1,J1,K)==0 .and. ivc(i1,j1,k)==0) THEN
		    distX = (vec(1)-xc(i1))
		    distY = (vec(2)-yc(j1))
		    dist_sqr(m) = distX**2 + distY**2
		  ENDIF
	    ENDDO
	    ENDDO
	  CASE(DIM_3D)
	    DO k1=k-wo, k+wo
	    DO j1=j-wo, j+wo
	    DO i1=i-wo, i+wo
		  m=m+1
		  IF (iblank(I1,J1,K1)==0 .and. ivc(i1,j1,k1)==0) THEN
		    distX = (vec(1)-xc(i1))
		    distY = (vec(2)-yc(j1))
		    distZ = (vec(3)-zc(k1))
		    dist_sqr(m) = distX**2 + distY**2 + distZ**2
		  ENDIF
	    ENDDO
	    ENDDO
	    ENDDO
	  END SELECT ! canonical_body_type
    iclose            = MINLOC(dist_sqr)
    body_dist_min1 = iclose(1)
    dist_sqr(body_dist_min1)=1.0E10_CGREAL
    iclose            = MINLOC(dist_sqr)
    body_dist_min2 = iclose(1)
	  dist_sqr(body_dist_min2)=1.0E10_CGREAL
	  iclose            = MINLOC(dist_sqr)
	  body_dist_min3 = iclose(1)
    IF (ndim==DIM_2D) THEN
	    j2=j+(body_dist_min1-1)/ws-wo
	    i2=i+mod(body_dist_min1-1,ws)-wo

	    j3=j+(body_dist_min2-1)/ws-wo
	    i3=i+mod(body_dist_min2-1,ws)-wo

	    j4=j+(body_dist_min3-1)/ws-wo
	    i4=i+mod(body_dist_min3-1,ws)-wo

	    dx1=face(cell(ivc(i,j,k))%F_slice)%centroid(1)-vec(1)
	    dy1=face(cell(ivc(i,j,k))%F_slice)%centroid(2)-vec(2)
	    dx2=xc(i2)-vec(1)
	    dy2=yc(j2)-vec(2)
	    dx3=xc(i3)-vec(1)
	    dy3=yc(j3)-vec(2)
	  
	    IF ((dx2*dy1 - dx3*dy1 - dx1*dy2 + dx3*dy2 + dx1*dy3 - dx2*dy3)==ZERO) THEN
	      dx3=xc(i4)-vec(1)
	      dy3=yc(j4)-vec(2)
	    ENDIF
      c1=(dx3*dy2 - dx2*dy3)/(dx2*dy1 - dx3*dy1 - dx1*dy2 + dx3*dy2 + dx1*dy3 - dx2*dy3)

      c2=(-dx3*dy1 + dx1*dy3)/(dx2*dy1 - dx3*dy1 - dx1*dy2 + dx3*dy2 + dx1*dy3 - dx2*dy3)

      c3=(-dx2*dy1 + dx1*dy2)/(-dx2*dy1 + dx3*dy1 + dx1*dy2 - dx3*dy2 - dx1*dy3 + dx2*dy3)
     
      u(i,j,k)=c1*ub(1)+c2*u(i2,j2,k)+c3*u(i3,j3,k)
      v(i,j,k)=c1*ub(2)+c2*v(i2,j2,k)+c3*v(i3,j3,k)
    ELSE

    
	    k2=k+(body_dist_min1-1)/ws2-wo
	    j2=j+mod(body_dist_min1-1,ws2)/ws-wo
	    i2=i+mod(mod(body_dist_min1-1,ws2),ws)-wo

	    k3=k+(body_dist_min2-1)/ws2-wo
	    j3=j+mod(body_dist_min2-1,ws2)/ws-wo
	    i3=i+mod(mod(body_dist_min2-1,ws2),ws)-wo

	    k4=k+(body_dist_min3-1)/ws2-wo
	    j4=j+mod(body_dist_min3-1,ws2)/ws-wo
	    i4=i+mod(mod(body_dist_min3-1,ws2),ws)-wo
	
	  dx1=face(cell(ivc(i,j,k))%F_slice)%centroid(1)-vec(1)
	  dy1=face(cell(ivc(i,j,k))%F_slice)%centroid(2)-vec(2)
	  dz1=face(cell(ivc(i,j,k))%F_slice)%centroid(3)-vec(3)
	  dx2=xc(i2)-vec(1)
	  dy2=yc(j2)-vec(2)
	  dz2=zc(k2)-vec(3)
	  dx3=xc(i3)-vec(1)
	  dy3=yc(j3)-vec(2)
	  dz3=zc(k3)-vec(3)
	  dx4=xc(i4)-vec(1)
	  dy4=yc(j4)-vec(2)
	  dz4=zc(k4)-vec(3)
	  
c1=(-dx4*dy3*dz2 + dx3*dy4*dz2 + dx4*dy2*dz3 - dx2*dy4*dz3 - &
     dx3*dy2*dz4 + dx2*dy3*dz4)/(dx3*dy2*dz1 - dx4*dy2*dz1 - &
     dx2*dy3*dz1 + dx4*dy3*dz1 + dx2*dy4*dz1 - dx3*dy4*dz1 - &
     dx3*dy1*dz2 + dx4*dy1*dz2 + dx1*dy3*dz2 - dx4*dy3*dz2 - &
     dx1*dy4*dz2 + dx3*dy4*dz2 + dx2*dy1*dz3 - dx4*dy1*dz3 - &
     dx1*dy2*dz3 + dx4*dy2*dz3 + dx1*dy4*dz3 - dx2*dy4*dz3 - &
     dx2*dy1*dz4 + dx3*dy1*dz4 + dx1*dy2*dz4 - dx3*dy2*dz4 - &
     dx1*dy3*dz4 + dx2*dy3*dz4)
      
c2=(-dx4*dy3*dz1 + dx3*dy4*dz1 + dx4*dy1*dz3 - dx1*dy4*dz3 - &
     dx3*dy1*dz4 + dx1*dy3*dz4)/(-dx3*dy2*dz1 + dx4*dy2*dz1 + & 
     dx2*dy3*dz1 - dx4*dy3*dz1 - dx2*dy4*dz1 + dx3*dy4*dz1 + &
     dx3*dy1*dz2 - dx4*dy1*dz2 - dx1*dy3*dz2 + dx4*dy3*dz2 + &
     dx1*dy4*dz2 - dx3*dy4*dz2 - dx2*dy1*dz3 + dx4*dy1*dz3 + &
     dx1*dy2*dz3 - dx4*dy2*dz3 - dx1*dy4*dz3 + dx2*dy4*dz3 + &
     dx2*dy1*dz4 - dx3*dy1*dz4 - dx1*dy2*dz4 + dx3*dy2*dz4 + &
     dx1*dy3*dz4 - dx2*dy3*dz4)
   
c3=(-dx4*dy2*dz1 + dx2*dy4*dz1 + dx4*dy1*dz2 - dx1*dy4*dz2 - &
     dx2*dy1*dz4 + dx1*dy2*dz4)/(dx3*dy2*dz1 - dx4*dy2*dz1 - &
     dx2*dy3*dz1 + dx4*dy3*dz1 + dx2*dy4*dz1 - dx3*dy4*dz1 - &
     dx3*dy1*dz2 + dx4*dy1*dz2 + dx1*dy3*dz2 - dx4*dy3*dz2 - &
     dx1*dy4*dz2 + dx3*dy4*dz2 + dx2*dy1*dz3 - dx4*dy1*dz3 - &
     dx1*dy2*dz3 + dx4*dy2*dz3 + dx1*dy4*dz3 - dx2*dy4*dz3 - &
     dx2*dy1*dz4 + dx3*dy1*dz4 + dx1*dy2*dz4 - dx3*dy2*dz4 - &
     dx1*dy3*dz4 + dx2*dy3*dz4)

c4=(-dx3*dy2*dz1 + dx2*dy3*dz1 + dx3*dy1*dz2 - dx1*dy3*dz2 - &
     dx2*dy1*dz3 + dx1*dy2*dz3)/(-dx3*dy2*dz1 + dx4*dy2*dz1 + & 
     dx2*dy3*dz1 - dx4*dy3*dz1 - dx2*dy4*dz1 + dx3*dy4*dz1 + & 
     dx3*dy1*dz2 - dx4*dy1*dz2 - dx1*dy3*dz2 + dx4*dy3*dz2 + &
     dx1*dy4*dz2 - dx3*dy4*dz2 - dx2*dy1*dz3 + dx4*dy1*dz3 + &
     dx1*dy2*dz3 - dx4*dy2*dz3 - dx1*dy4*dz3 + dx2*dy4*dz3 + &
     dx2*dy1*dz4 - dx3*dy1*dz4 - dx1*dy2*dz4 + dx3*dy2*dz4 + &
     dx1*dy3*dz4 - dx2*dy3*dz4)
   
    u(i,j,k)=c1*ub(1)+c2*u(i2,j2,k2)+c3*u(i3,j3,k3)+c4*u(i4,j4,k4)
    v(i,j,k)=c1*ub(2)+c2*v(i2,j2,k2)+c3*v(i3,j3,k3)+c4*v(i4,j4,k4)
    w(i,j,k)=c1*ub(3)+c2*w(i2,j2,k2)+c3*w(i3,j3,k3)+c4*w(i4,j4,k4)
	ENDIF
   

  END SUBROUTINE taylor_series_matching_cell
!-------------------------------------------------------------------------------

   END SUBROUTINE correct_cell_vel
!-------------------------------------------------------------------------------

   SUBROUTINE update_pressure()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays

    IMPLICIT NONE


    IF (frac_step_type == VAN_KAN) THEN
      p = p + pPrime
    ELSE
      p = pPrime
    ENDIF

    IF (boundary_formulation == GCM_METHOD) CALL GCM_p_set_bc_internal(p,nx+1,ny+1,nz+1)

  END SUBROUTINE update_pressure

  !-----------------------------------------------------------------------------

subroutine faceVelCorret()
USE global_parameters
USE flow_parameters
USE flow_arrays
USE boundary_arrays
USE grid_arrays
USE GCM_arrays
USE multiuse_arrays
USE pressure_arrays
USE solver_arrays
use hybrid_cell_arrays

IMPLICIT NONE

INTEGER :: i,j,k,iBody,n


real(cgreal) :: dist,distx,disty,distz,distb,ratiox,ratioy,ratioz
real(cgreal) :: xBIbeta,yBIbeta,zBIbeta
real(cgreal) :: xBIalpha1,yBIalpha1,zBIalpha1
real(cgreal) :: xBIalpha2,yBIalpha2,zBIalpha2
real(cgreal) :: xBIalpha3,yBIalpha3,zBIalpha3
integer :: h_mark(3)

integer :: testx(0:nx+1,0:ny+1,0:nz+1)
integer :: testy(0:nx+1,0:ny+1,0:nz+1)

testx=0
testy=0

do k=1,nzc
do j=1,nyc
do i=1,nxc
    if(iup(i,j,k)==1) testx(i,j,k)=-1
    if(ium(i,j,k)==1) testx(i,j,k)=1
    if(jup(i,j,k)==1) testy(i,j,k)=-1
    if(jum(i,j,k)==1) testy(i,j,k)=1




end do
end do
end do

!if(ntime==2)then
!    write(568,*) 'VARIABLES="X","Y","hybridMarkMemb","hybrid_mark","testx","testy"'
!    write(568,*) 'ZONE F=POINT, I=',nxc,' , J=',nyc
!    do j=1,nyc
!    do i=1,nxc
!        write(568,*) xc(i),yc(j),hybridMarkMemb(i,j,1),hybrid_mark(i,j,1),testx(i,j,1),testy(i,j,1)
!    end do
!    end do
!end if


    do n=1,nhybrid
    
    xBIbeta=zero
    yBIbeta=zero
    zBIbeta=zero

    xBIalpha1=zero
    yBIalpha1=zero
    zBIalpha1=zero
    xBIalpha2=zero
    yBIalpha2=zero
    zBIalpha2=zero
    xBIalpha3=zero
    yBIalpha3=zero
    zBIalpha3=zero

    dist=zero
    distx=zero
    disty=zero
    distz=zero
    distb=zero

    i=ihybrid(n)
    j=jhybrid(n)
    k=khybrid(n)
    h_mark=(/i,j,k/)
    !if(n==100) write(*,*) i,j,k
    CALL hybrid_Calc_BodyIntercept_Unstruc( i, j, k, xc(i), yc(j), zc(k),    &
                                        xBIbeta, yBIbeta, zBIbeta, closestElementHC(n) )
    distb=sqrt((xBIbeta-xc(i))**2+(yBIbeta-yc(j))**2+(zBIbeta-zc(k))**2)
    if(iup(i,j,k)==1) then
        CALL hybrid_Calc_BodyIntercept_Unstruc( i+1, j, k, xc(i+1), yc(j), zc(k),    &
                                    xBIalpha1, yBIalpha1, zBIalpha1, closestElementR1(n) )
        distx=sqrt((xBIalpha1-xc(i+1))**2+(yBIalpha1-yc(j))**2+(zBIalpha1-zc(k))**2)
    end if
    if(ium(i,j,k)==1) then
        CALL hybrid_Calc_BodyIntercept_Unstruc( i-1, j, k, xc(i-1), yc(j), zc(k),    &
                                    xBIalpha1, yBIalpha1, zBIalpha1, closestElementR1(n) )
        distx=sqrt((xBIalpha1-xc(i-1))**2+(yBIalpha1-yc(j))**2+(zBIalpha1-zc(k))**2)
    end if
    if(jup(i,j,k)==1) then
        CALL hybrid_Calc_BodyIntercept_Unstruc( i, j+1, k, xc(i), yc(j+1), zc(k),    &
                                    xBIalpha2, yBIalpha2, zBIalpha2, closestElementR1(n) )
        disty=sqrt((xBIalpha2-xc(i))**2+(yBIalpha2-yc(j+1))**2+(zBIalpha2-zc(k))**2)
    end if
    if(jum(i,j,k)==1) then
        CALL hybrid_Calc_BodyIntercept_Unstruc( i, j-1, k, xc(i), yc(j-1), zc(k),    &
                                    xBIalpha2, yBIalpha2, zBIalpha2, closestElementR1(n) )
        disty=sqrt((xBIalpha2-xc(i))**2+(yBIalpha2-yc(j-1))**2+(zBIalpha2-zc(k))**2)
    end if
    if(ndim==dim_3d)then
        if(kup(i,j,k)==1) then
            CALL hybrid_Calc_BodyIntercept_Unstruc( i, j, k+1, xc(i), yc(j), zc(k+1),    &
                                        xBIalpha3, yBIalpha3, zBIalpha3, closestElementR1(n) )
            distz=sqrt((xBIalpha3-xc(i))**2+(yBIalpha3-yc(j))**2+(zBIalpha3-zc(k+1))**2)
        end if
        if(kum(i,j,k)==1) then
            CALL hybrid_Calc_BodyIntercept_Unstruc( i, j, k-1, xc(i), yc(j), zc(k-1),    &
                                        xBIalpha3, yBIalpha3, zBIalpha3, closestElementR1(n) )
            distz=sqrt((xBIalpha3-xc(i))**2+(yBIalpha3-yc(j))**2+(zBIalpha3-zc(k-1))**2)
        end if
    end if


    ratiox=distx/(distx+distb)
    ratioy=disty/(disty+distb)

    if(ndim==dim_3d)then
        ratioz=distz/(distz+distb)
    end if
    
    if(i==8.and.j==23.and.k==2)then
        write(570,*) ntime,ratiox,ratioy
    end if



    
        if((ium(i,j,k)==1.or.iup(i,j,k)==1).and.(jum(i,j,k)==0.and.jup(i,j,k)==0))then
            if(ium(i,j,k)==1)then
!                write(*,*) '------------------'
!                write(*,*) face_u(i,j,k)
!                write(*,*) face_u(i+1,j,k)
                if(ratiox<=half)then
                    face_u(i,j,k)=half*(u(i-1,j,k)+u(i,j,k))*(oned-twod*ratiox)+bcxu(i,j,k)*(twod*ratiox)
                    !face_u(i+1,j,k)=bcxu(i,j,k)*ratioX+face_u(i+1,j,k)*(oned-ratioX)
                else
                    face_u(i,j,k)=half*(u(i-1,j,k)+u(i,j,k))*(twod*ratiox-oned)+bcxu(i,j,k)*(twod-twod*ratiox)
                    !face_u(i+1,j,k)=bcxu(i,j,k)*ratioX+face_u(i+1,j,k)*(oned-ratioX)
                end if
!                write(*,*) face_u(i,j,k)
!                write(*,*) face_u(i+1,j,k)
!                face_v(i,j+1,k)=bcxv(i,j,k)*ratioX+face_v(i,j+1,k)*(oned-ratioX)
!                face_v(i,j,k)=bcxv(i,j,k)*ratioX+face_v(i,j,k)*(oned-ratioX)
            else
                if(ratiox<=half)then
                    face_u(i+1,j,k)=half*(u(i+1,j,k)+u(i,j,k))*(oned-twod*ratiox)+bcxu(i,j,k)*(twod*ratiox)
                    !face_u(i,j,k)=bcxu(i,j,k)*ratioX+face_u(i,j,k)*(oned-ratioX)
                else
                    face_u(i+1,j,k)=half*(u(i+1,j,k)+u(i,j,k))*(twod*ratiox-oned)+bcxu(i,j,k)*(twod-twod*ratiox)
                    !face_u(i,j,k)=bcxu(i,j,k)*ratioX+face_u(i,j,k)*(oned-ratioX)
                end if
!                face_v(i,j+1,k)=bcxv(i,j,k)*ratioX+face_v(i,j+1,k)*(oned-ratioX)
!                face_v(i,j,k)=bcxv(i,j,k)*ratioX+face_v(i,j,k)*(oned-ratioX)
            end if
        end if
        
        if((jum(i,j,k)==1.or.jup(i,j,k)==1).and.(ium(i,j,k)==0.and.iup(i,j,k)==0))then
            if(jum(i,j,k)==1)then
                !face_v(i,j,k)=face_v(i,j,k)*ratioY
                !face_v(i,j+1,k)=bcyv(i,j,k)*ratioY+face_v(i,j+1,k)*(oned-ratioY)
!                face_u(i+1,j,k)=bcyu(i,j,k)*ratioY+face_u(i+1,j,k)*(oned-ratioY)
!                face_u(i,j,k)=bcyu(i,j,k)*ratioY+face_u(i,j,k)*(oned-ratioY)
            else
                !face_v(i,j+1,k)=face_v(i,j+1,k)*ratioY
                !face_v(i,j,k)=bcyv(i,j,k)*ratioY+face_v(i,j,k)*(oned-ratioY)
!                face_u(i+1,j,k)=bcyu(i,j,k)*ratioY+face_u(i+1,j,k)*(oned-ratioY)
!                face_u(i,j,k)=bcyu(i,j,k)*ratioY+face_u(i,j,k)*(oned-ratioY)
            end if
        end if

        if((iup(i,j,k)==1.or.ium(i,j,k)==1).and.(jup(i,j,k)==1.or.jum(i,j,k)==1))then
            if(iup(i,j,k)==1.and.jup(i,j,k)==1)then
                face_u(i,j,k)=bcxu(i,j,k)*ratioX+face_u(i,j,k)*(oned-ratioX)
                face_v(i,j,k)=bcyv(i,j,k)*ratioY+face_v(i,j,k)*(oned-ratioY)
            else if(iup(i,j,k)==1.and.jum(i,j,k)==1)then
                face_u(i,j,k)=bcxu(i,j,k)*ratioX+face_u(i,j,k)*(oned-ratioX)
                face_v(i,j+1,k)=bcyv(i,j,k)*ratioY+face_v(i,j+1,k)*(oned-ratioY)
            else if(ium(i,j,k)==1.and.jup(i,j,k)==1)then
                face_u(i+1,j,k)=bcxu(i,j,k)*ratioX+face_u(i+1,j,k)*(oned-ratioX)
                face_v(i,j,k)=bcyv(i,j,k)*ratioY+face_v(i,j,k)*(oned-ratioY)
            else if(ium(i,j,k)==1.and.jum(i,j,k)==1)then
                face_u(i+1,j,k)=bcxu(i,j,k)*ratioX+face_u(i+1,j,k)*(oned-ratioX)
                face_v(i,j+1,k)=bcyv(i,j,k)*ratioY+face_v(i,j+1,k)*(oned-ratioY)
            end if
        end if

        if(ndim==dim_3D)then
            
        end if

end do






end subroutine faceVelCorret