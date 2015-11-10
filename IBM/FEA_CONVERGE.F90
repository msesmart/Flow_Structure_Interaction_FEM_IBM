   SUBROUTINE fea_converge(iBody)

   USE fea_unstructure_surface
   USE unstructured_surface_arrays
   USE flow_parameters
   USE boundary_arrays, ONLY : uBodyMarker,vBodyMarker,wBodyMarker
   USE usr_module, ONLY : uBodyMarker_iter,vBodyMarker_iter,wBodyMarker_iter, &
                          uBodyMarker_worelax,vBodyMarker_worelax,wBodyMarker_worelax

   IMPLICIT NONE

   INTEGER :: iBody
   INTEGER :: MAX_ITERATION
   INTEGER :: i,j
   REAL(KIND=CGREAL),DIMENSION(nodeDoF) :: maxvel_err
   REAL(KIND=CGREAL),DIMENSION(nodeDoF) :: tmp

   MAX_ITERATION = 50

   maxvel_err(:) = 0
   tmp(:) = 0

   print *, 'In converge max/min struc_vel_iter(,2) =', MAXVAL(struc_vel_iter(:,2)), MINVAL(struc_vel_iter(:,2))
   print *, 'In converge max/min struc_vel(,2)      =', MAXVAL(struc_vel(:,2)),      MINVAL(struc_vel(:,2))

   DO i=1,nPtsBodyMarker(iBody)
      DO j=1,nodeDoF
         tmp(j) = dabs( (struc_vel(i,j)-struc_vel_iter(i,j))/struc_vel(i,j) )
         IF (tmp(j)>maxvel_err(j)) maxvel_err(j) = tmp(j)
      ENDDO
   ENDDO

!  IF (MAXVAL(maxvel_err) < FSIConvg_Criteria) THEN 
   IF (maxvel_err(2) < FSIConvg_Criteria) THEN 
      FSI_CONVERGE = .TRUE.
      FSI_ITERATION = FSI_ITERATION+1

      WRITE (*,*) 'FSI iteration number is:', FSI_ITERATION
      WRITE (*,*) 'Max FSI vel_err of y-direction is: ', maxvel_err(2)

   ELSE
      FSI_ITERATION = FSI_ITERATION+1
      WRITE (*,*) 'FSI iteration number is:', FSI_ITERATION
      WRITE (*,*) 'Max FSI vel_err of y-direction is: ', maxvel_err(2)
      IF (FSI_ITERATION > MAX_ITERATION) THEN
         WRITE (*,*) 'FSI iteration number reached Max number:', MAX_ITERATION
         STOP
      ENDIF

   ENDIF

   DO i=1,nPtsBodyMarker(iBody)
      uBodyMarker(iBody,i) = struc_vel(i,1)
      vBodyMarker(iBody,i) = struc_vel(i,2)
      wBodyMarker(iBody,i) = struc_vel(i,3)
   ENDDO

   END SUBROUTINE fea_converge



!-------------------------------------------------------------------------------------------------------
   SUBROUTINE fea_converge_static(iBody)

   USE fea_unstructure_surface
   USE unstructured_surface_arrays
   USE flow_parameters
   USE boundary_arrays, ONLY : uBodyMarker,vBodyMarker,wBodyMarker
   USE usr_module, ONLY : uBodyMarker_iter,vBodyMarker_iter,wBodyMarker_iter, &
                          uBodyMarker_worelax,vBodyMarker_worelax,wBodyMarker_worelax

   IMPLICIT NONE

   INTEGER :: iBody
   INTEGER :: MAX_ITERATION
   INTEGER :: i,j
   REAL(KIND=CGREAL),DIMENSION(nodeDoF) :: maxdisp_err
   REAL(KIND=CGREAL),DIMENSION(nodeDoF) :: tmp

   MAX_ITERATION = 50

   maxdisp_err(:) = 0
   tmp(:) = 0

   struc_vel_iter(:,:) = struc_vel(:,:)         !For only one body

   print *, ' '
   print *, 'In fea_converge_static:'
   WRITE (*,*) 'FSI iteration number is:', FSI_ITERATION, ' at Ntime', Ntime
   print *, 'In converge max/min struc_disp_iter(,2) =', MAXVAL(struc_disp_iter(:,2)), MINVAL(struc_disp_iter(:,2))
   print *, 'In converge max/min struc_disp(,2)      =', MAXVAL(struc_disp(:,2)),      MINVAL(struc_disp(:,2))

   DO i=1,nPtsBodyMarker(iBody)
      struc_vel_worelax(i,1) = (struc_disp_worelax(i,1)-struc_olddisp(i,1))/dt
      struc_vel_worelax(i,2) = (struc_disp_worelax(i,2)-struc_olddisp(i,2))/dt
      struc_vel_worelax(i,3) = (struc_disp_worelax(i,3)-struc_olddisp(i,3))/dt
   ENDDO

   print *, 'struc_disp_worelax 1,250:',   struc_disp_worelax(1,2),struc_disp_worelax(250,2)
   print *, 'struc_disp_worelax 500,750:', struc_disp_worelax(500,2),struc_disp_worelax(750,2)
!   print *, 'struc_vel_worelax 1,250:', struc_vel_worelax(1,2),struc_vel_worelax(250,2)
!   print *, 'struc_vel_worelax 500,750:', struc_vel_worelax(500,2),struc_vel_worelax(750,2)

!   CALL fea_underrelation_vel(iBody)
   CALL fea_underrelation_disp(iBody)

   print *, 'Max/Min struc_disp(:,2):', Maxval(struc_disp(:,2)),Minval(struc_disp(:,2))
   print *, 'after relax, struc_disp at 1,250:',   struc_disp(1,2),struc_disp(250,2)
   print *, 'after relax, struc_disp at 500,750:', struc_disp(500,2),struc_disp(750,2)

!   DO i=1,nPtsBodyMarker(iBody)
!      DO j=1,nodeDoF
!         tmp(j) = dabs( (struc_disp(i,j)-struc_disp_iter(i,j))/struc_disp(i,j) )
!         IF (tmp(j)>maxdisp_err(j)) maxdisp_err(j) = tmp(j)
!      ENDDO
!   ENDDO

!   IF (maxdisp_err(2) < FSIConvg_Criteria) THEN 
!      FSI_CONVERGE = .TRUE.
!      WRITE (*,*) 'Max FSI disp_err of y-direction is: ', maxdisp_err(2)
!
!   ELSE
!      WRITE (*,*) 'Max FSI disp_err of y-direction is: ', maxdisp_err(2)
!      IF (FSI_ITERATION > MAX_ITERATION) THEN
!         WRITE (*,*) 'FSI iteration number reached Max number:', MAX_ITERATION
!         STOP
!      ENDIF
!
!   ENDIF

   DO i=1,nPtsBodyMarker(iBody)
      tmp(2) = dabs( (struc_disp(i,2)-struc_disp_iter(i,2))/struc_disp(i,2) )
      IF (tmp(2)>maxdisp_err(2)) maxdisp_err(2) = tmp(2)
   ENDDO

   IF (maxdisp_err(2) < FSIConvg_Criteria) THEN 
      FSI_CONVERGE = .TRUE.
      WRITE (*,*) 'Max FSI stru_disp err in y-direction is: ', maxdisp_err(2)

   ELSE
      FSI_ITERATION = FSI_ITERATION+1
      WRITE (*,*) 'Max FSI stru_disp err in y-direction is: ', maxdisp_err(2)

      DO i=1,nPtsBodyMarker(iBody)
!         uBodyMarker(iBody,i) = uBodyMarker(iBody,1)+struc_vel(i,1)
!         vBodyMarker(iBody,i) = vBodyMarker(iBody,1)+struc_vel(i,2)
!         wBodyMarker(iBody,i) = wBodyMarker(iBody,1)+struc_vel(i,3)
          struc_vel(i,1) = (struc_disp(i,1)-struc_olddisp(i,1))/dt
          struc_vel(i,2) = (struc_disp(i,2)-struc_olddisp(i,2))/dt
          struc_vel(i,3) = (struc_disp(i,3)-struc_olddisp(i,3))/dt
      ENDDO !  i-loop

      IF (FSI_ITERATION > MAX_ITERATION) THEN
         WRITE (*,*) 'FSI iteration number reached Max number:', MAX_ITERATION
         STOP
      ENDIF

   ENDIF

   END SUBROUTINE fea_converge_static


   SUBROUTINE fea_underrelation_vel(iBody)
   USE fea_unstructure_surface
   USE boundary_arrays
   USE flow_parameters
   USE usr_module, ONLY : uBodyMarker_iter,vBodyMarker_iter,wBodyMarker_iter, &
                          uBodyMarker_worelax,vBodyMarker_worelax,wBodyMarker_worelax

   IMPLICIT NONE

   INTEGER :: iBody
   INTEGER :: i

   REAL(KIND=CGREAL) :: dAitkenu,dAitkenv,dAitkenw
   REAL(KIND=CGREAL) :: Aitken_nominator, Aitken_denominator

   REAL(KIND=CGREAL), PARAMETER  :: NUM_ZERO = 1E-30

   IF (AITKEN) THEN

      IF (FSI_ITERATION < 1e-1) THEN
         Lambdau_Aitken(iBody) = Lambdau_Aitken_Init
         Lambdav_Aitken(iBody) = Lambdav_Aitken_Init
         Lambdaw_Aitken(iBody) = Lambdaw_Aitken_Init

         DO i = 1,nPtsBodyMarker(iBody)
            deltau_Aitken = struc_vel(i,1) - struc_vel_worelax(i,1)
            deltav_Aitken = struc_vel(i,2) - struc_vel_worelax(i,2)
            deltaw_Aitken = struc_vel(i,3) - struc_vel_worelax(i,3)

            struc_vel(i,1) = Lambdau_Aitken(iBody)*struc_vel(i,1) +  &
                             (1-Lambdau_Aitken(iBody))*struc_vel_worelax(i,1)
            struc_vel(i,2) = Lambdav_Aitken(iBody)*struc_vel(i,2) +  &
                             (1-Lambdaw_Aitken(iBody))*struc_vel_worelax(i,2)
            struc_vel(i,3) = Lambdav_Aitken(iBody)*struc_vel(i,3) +  &
                             (1-Lambdaw_Aitken(iBody))*struc_vel_worelax(i,3)
         ENDDO
      ELSE
         ! Relax for u velocity induced by struture.
         Aitken_denominator = 0
         Aitken_nominator   = 0
         deltau_Aitken_prev(:) = deltau_Aitken(:)

         DO i = 1,nPtsBodyMarker(iBody)
            deltau_Aitken(i) = struc_vel_iter(i,1) - struc_vel_worelax(i,1)
            dAitkenu = deltau_Aitken_prev(i) - deltau_Aitken(i)
            Aitken_denominator = Aitken_denominator + dAitkenu*dAitkenu 
            Aitken_nominator   = Aitken_nominator   + dAitkenu*deltau_Aitken(i)
         ENDDO ! i-loop

         IF (abs(Aitken_denominator)>NUM_ZERO) THEN
            Lambdau_Aitken(iBody) = Lambdau_Aitken(iBody) + (Lambdau_Aitken(iBody)-1)*   &
                                    Aitken_nominator/Aitken_denominator
         ENDIF

!         WRITE(*,*) 'UnderRelaxation factor for struc_vel(1):', 1-Lambdau_Aitken(iBody)

         DO i = 1,nPtsBodyMarker(iBody)
            struc_vel(i,1) = Lambdau_Aitken(iBody)*struc_vel_iter(i,1) +  &
                             (1-Lambdau_Aitken(iBody))*struc_vel_worelax(i,1)
         ENDDO ! i-loop

         ! Relax for v velocity induced by struture.
         Aitken_denominator = 0
         Aitken_nominator   = 0
         deltav_Aitken_prev(:) = deltav_Aitken(:)

         print *,'max/min deltav_Aitken_prev:', maxval(deltav_Aitken_prev),minval(deltav_Aitken_prev)
         print *,'max/min struc_vel_worelax:', maxval(struc_vel_worelax(:,2)),minval(struc_vel_worelax(:,2))
         print *,'max/min struc_vel_iter:', maxval(struc_vel_iter(:,2)),minval(struc_vel_iter(:,2))

         DO i = 1,nPtsBodyMarker(iBody)
            deltav_Aitken(i) = struc_vel_iter(i,2) - struc_vel_worelax(i,2)
           
            dAitkenv = deltav_Aitken_prev(i) - deltav_Aitken(i)

            Aitken_denominator = Aitken_denominator + dAitkenv*dAitkenv
            Aitken_nominator   = Aitken_nominator   + dAitkenv*deltav_Aitken(i)
         ENDDO ! i-loop

         print *,'max/min deltav_Aitken:', maxval(deltav_Aitken),minval(deltav_Aitken)
         print *, 'Aitken_denominator for v:', Aitken_denominator
         print *, 'Aitken_nominator for v:', Aitken_nominator

         IF (abs(Aitken_denominator)>NUM_ZERO) THEN
            Lambdav_Aitken(iBody) = Lambdav_Aitken(iBody) + (Lambdav_Aitken(iBody)-1)*   &
                                    Aitken_nominator/Aitken_denominator
         ENDIF

         WRITE(*,*) 'UnderRelaxation factor for struc_vel(2):', 1-Lambdav_Aitken(iBody)

         DO i = 1,nPtsBodyMarker(iBody)
            struc_vel(i,2) = Lambdav_Aitken(iBody)*struc_vel_iter(i,2) +  &
                             (1-Lambdav_Aitken(iBody))*struc_vel_worelax(i,2)
         ENDDO ! i-loop

         ! Relax for w velocity induced by struture.
         Aitken_denominator = 0
         Aitken_nominator   = 0
         deltaw_Aitken_prev(:) = deltaw_Aitken(:)

         DO i = 1,nPtsBodyMarker(iBody)
            deltaw_Aitken(i) = struc_vel_iter(i,3) - struc_vel_worelax(i,3)
            dAitkenw = deltaw_Aitken_prev(i) - deltaw_Aitken(i)
            Aitken_denominator = Aitken_denominator + dAitkenw*dAitkenw
            Aitken_nominator   = Aitken_nominator   + dAitkenw*deltaw_Aitken(i)
         ENDDO ! i-loop

         IF (abs(Aitken_denominator)>NUM_ZERO) THEN
            Lambdaw_Aitken(iBody) = Lambdaw_Aitken(iBody) + (Lambdaw_Aitken(iBody)-1)*   &
                                    Aitken_nominator/Aitken_denominator
         ENDIF

!         WRITE(*,*) 'UnderRelaxation factor for struc_vel(3):', 1-Lambdaw_Aitken(iBody)

         DO i = 1,nPtsBodyMarker(iBody)
            struc_vel(i,3) = Lambdaw_Aitken(iBody)*struc_vel_iter(i,3) +  &
                             (1-Lambdaw_Aitken(iBody))*struc_vel_worelax(i,3)
         ENDDO ! i-loop


      ENDIF !FSI_ITERATION

   ELSE  ! AITKEN = False

   ENDIF ! AITKEN


   END SUBROUTINE fea_underrelation_vel




   SUBROUTINE fea_underrelation_disp(iBody)
   USE fea_unstructure_surface
   USE boundary_arrays
   USE flow_parameters
   USE usr_module, ONLY : uBodyMarker_iter,vBodyMarker_iter,wBodyMarker_iter, &
                          uBodyMarker_worelax,vBodyMarker_worelax,wBodyMarker_worelax

   IMPLICIT NONE

   INTEGER :: iBody
   INTEGER :: i

   REAL(KIND=CGREAL) :: dAitkenu,dAitkenv,dAitkenw
   REAL(KIND=CGREAL) :: Aitken_nominator, Aitken_denominator

   REAL(KIND=CGREAL), PARAMETER  :: NUM_ZERO = 1E-30

   IF (AITKEN) THEN

      IF (FSI_ITERATION < 1e-1) THEN
         Lambdau_Aitken(iBody) = Lambdau_Aitken_Init
         Lambdav_Aitken(iBody) = Lambdav_Aitken_Init
         Lambdaw_Aitken(iBody) = Lambdaw_Aitken_Init

         DO i = 1,nPtsBodyMarker(iBody)
            deltau_Aitken = struc_disp(i,1) - struc_disp_worelax(i,1)
            deltav_Aitken = struc_disp(i,2) - struc_disp_worelax(i,2)
            deltaw_Aitken = struc_disp(i,3) - struc_disp_worelax(i,3)

            struc_disp(i,1) = Lambdau_Aitken(iBody)*struc_disp(i,1) +  &
                             (1-Lambdau_Aitken(iBody))*struc_disp_worelax(i,1)
            struc_disp(i,2) = Lambdav_Aitken(iBody)*struc_disp(i,2) +  &
                             (1-Lambdav_Aitken(iBody))*struc_disp_worelax(i,2)
            struc_disp(i,3) = Lambdaw_Aitken(iBody)*struc_disp(i,3) +  &
                             (1-Lambdaw_Aitken(iBody))*struc_disp_worelax(i,3)
         ENDDO
      ELSE
         ! Relax for disp_x from struture.
         Aitken_denominator = 0
         Aitken_nominator   = 0
         deltau_Aitken_prev(:) = deltau_Aitken(:)

         DO i = 1,nPtsBodyMarker(iBody)
            deltau_Aitken(i) = struc_disp_iter(i,1) - struc_disp_worelax(i,1)
            dAitkenu = deltau_Aitken_prev(i) - deltau_Aitken(i)
            Aitken_denominator = Aitken_denominator + dAitkenu*dAitkenu 
            Aitken_nominator   = Aitken_nominator   + dAitkenu*deltau_Aitken(i)
         ENDDO ! i-loop

         IF (abs(Aitken_denominator)>NUM_ZERO) THEN
            Lambdau_Aitken(iBody) = Lambdau_Aitken(iBody) + (Lambdau_Aitken(iBody)-1)*   &
                                    Aitken_nominator/Aitken_denominator
         ENDIF

!         WRITE(*,*) 'UnderRelaxation factor for struc_disp(1):', 1-Lambdau_Aitken(iBody)

         DO i = 1,nPtsBodyMarker(iBody)
            struc_disp(i,1) = Lambdau_Aitken(iBody)*struc_disp_iter(i,1) +  &
                             (1-Lambdau_Aitken(iBody))*struc_disp_worelax(i,1)
         ENDDO ! i-loop

         ! Relax for disp_y from struture.
         Aitken_denominator = 0
         Aitken_nominator   = 0
         deltav_Aitken_prev(:) = deltav_Aitken(:)

         print *,'max/min deltav_Aitken_prev:', maxval(deltav_Aitken_prev),minval(deltav_Aitken_prev)
         print *,'max/min struc_disp_worelax:', maxval(struc_disp_worelax(:,2)),minval(struc_disp_worelax(:,2))
         print *,'max/min struc_disp_iter:', maxval(struc_disp_iter(:,2)),minval(struc_disp_iter(:,2))

         DO i = 1,nPtsBodyMarker(iBody)
            deltav_Aitken(i) = struc_disp_iter(i,2) - struc_disp_worelax(i,2)
           
            dAitkenv = deltav_Aitken_prev(i) - deltav_Aitken(i)

            Aitken_denominator = Aitken_denominator + dAitkenv*dAitkenv
            Aitken_nominator   = Aitken_nominator   + dAitkenv*deltav_Aitken(i)
         ENDDO ! i-loop

         print *,'max/min deltav_Aitken:', maxval(deltav_Aitken),minval(deltav_Aitken)
         print *, 'Aitken_denominator for v:', Aitken_denominator
         print *, 'Aitken_nominator for v:', Aitken_nominator

         IF (abs(Aitken_denominator)>NUM_ZERO) THEN
            Lambdav_Aitken(iBody) = Lambdav_Aitken(iBody) + (Lambdav_Aitken(iBody)-1)*   &
                                    Aitken_nominator/Aitken_denominator
         ENDIF

         WRITE(*,*) 'UnderRelaxation factor for struc_disp(2):', 1-Lambdav_Aitken(iBody)

         DO i = 1,nPtsBodyMarker(iBody)
            struc_disp(i,2) = Lambdav_Aitken(iBody)*struc_disp_iter(i,2) +  &
                             (1-Lambdav_Aitken(iBody))*struc_disp_worelax(i,2)
         ENDDO ! i-loop

         ! Relax for disp_z from struture.
         Aitken_denominator = 0
         Aitken_nominator   = 0
         deltaw_Aitken_prev(:) = deltaw_Aitken(:)

         DO i = 1,nPtsBodyMarker(iBody)
            deltaw_Aitken(i) = struc_disp_iter(i,3) - struc_disp_worelax(i,3)
            dAitkenw = deltaw_Aitken_prev(i) - deltaw_Aitken(i)
            Aitken_denominator = Aitken_denominator + dAitkenw*dAitkenw
            Aitken_nominator   = Aitken_nominator   + dAitkenw*deltaw_Aitken(i)
         ENDDO ! i-loop

         IF (abs(Aitken_denominator)>NUM_ZERO) THEN
            Lambdaw_Aitken(iBody) = Lambdaw_Aitken(iBody) + (Lambdaw_Aitken(iBody)-1)*   &
                                    Aitken_nominator/Aitken_denominator
         ENDIF

!         WRITE(*,*) 'UnderRelaxation factor for struc_disp(3):', 1-Lambdaw_Aitken(iBody)

         DO i = 1,nPtsBodyMarker(iBody)
            struc_disp(i,3) = Lambdaw_Aitken(iBody)*struc_disp_iter(i,3) +  &
                             (1-Lambdaw_Aitken(iBody))*struc_disp_worelax(i,3)
         ENDDO ! i-loop

         DO i = 1,nPtsBodyMarker(iBody)
            struc_disp(i,4) = struc_disp_worelax(i,4)
            struc_disp(i,5) = struc_disp_worelax(i,5)
            struc_disp(i,6) = struc_disp_worelax(i,6)
         ENDDO ! i-loop

      ENDIF !FSI_ITERATION

   ELSE  ! AITKEN = False

   ENDIF ! AITKEN

   END SUBROUTINE fea_underrelation_disp
