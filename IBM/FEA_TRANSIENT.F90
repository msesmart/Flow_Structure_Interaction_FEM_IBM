!  TRNSIENT analysis  by NEWMARK time integration
   SUBROUTINE FEA_TRANSIENT(iBody)

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface
   USE flow_parameters, ONLY : nPtsMax, dt, nPtsBodyMarker
   
   IMPLICIT NONE
   
   INTEGER :: iBody
   INTEGER :: node,jdof
   INTEGER :: idis

   DOUBLE PRECISION :: GAMMA,BETA

   PARAMETER(GAMMA=0.5e0,BETA=0.25e0)  !Book of Doyle (on Sta & Dyn), pp354

   REAL(KIND=CGREAL) :: a0,a1,a2,a3,a4,a5,a6,a7
   REAL(KIND=CGREAL) :: LOADEFF(NEQ)
   REAL(KIND=CGREAL) :: STFTMP(MAXSTIFF),STFRHS(MAXSTIFF)
   REAL(KIND=CGREAL) :: MSSTMP(MAXSTIFF)

   REAL(KIND=CGREAL) :: dampK(nPtsMax*nodeDoF)
   REAL(KIND=CGREAL) :: olddis(nPtsMax*nodeDoF)
   REAL(KIND=CGREAL) :: oldVEL,oldACC
   REAL(KIND=CGREAL) :: aterm, atermc
   REAL(KIND=CGREAL) :: dmm, dkk
   INTEGER :: I,J,ILOC,IERROR
   
   dmm = 0.0
   dkk = 0.0

!  Set integeration constants for Newmark method

   a0 = 1.0e0/(BETA*dt*dt)
   a1 = GAMMA/(BETA*dt)
   a2 = 1.0e0/(BETA*dt)
   a3 = 1.0e0/(BETA*2.0e0) - 1.0e0
   a4 = GAMMA/BETA- 1.0e0
   a5 = (GAMMA/BETA- 2.0e0)*0.5e0*dt
   a6 = (1.0e0 - GAMMA)*dt
   a7 = GAMMA*dt

!  Form effective stiffness matrix by adding inertia & damping 
!  only proportional damping

   DO i= 1, NEQ  
      iloc=nloc(i)
      DO j= 1, iprof(i)
         STFTMP(iloc+j-1) = stf(iloc+j-1) + a0*MSS(iloc+j-1)  &
                       + a1*(dkk*stf(iloc+j-1)+dmm*MSS(iloc+j-1))

         MSSTMP(iloc+j-1) = MSS(iloc+j-1)

         STFRHS(iloc+j-1) = stf(iloc+j-1)
      ENDDO
   ENDDO

   CALL uduCOL(STFTMP,ierror)
   if (ierror.eq.0) then
      write(*,*)'ERROR: zero diagonal term in transient after udu decompose.'
      stop
   endif

   olddis(:) = DISP(:)

!  Form effective load vector by adding  + MSS x ACCeleration
   DO i= 1, NEQ  
      aterm  = a0*DISP(i) + a2*VEL(i) + a3*ACC(i)   !RHS coeff. of M
      atermc = a1*DISP(i) + a4*VEL(i) + a5*ACC(i)   !RHS coeff. of C
      DISP(i) = aterm + dmm*atermc
      dampK(i) = dkk*atermc
   ENDDO

   LOADEFF(:) = 0

   CALL AxBCOL(MSSTMP,DISP,LOADEFF) !doing LOADEFF = DISP*MSS

   CALL AxBCOL(STFRHS,dampK,LOADEFF) !doing LOADEFF = LOADEFF++DISP*STF

   DO I=1,NEQ
      LOADEFF(I) = LOADEFF(I)+LOAD(I)
   ENDDO


!  Solve for new DISPlacements: UDU already obtained, do back-substitution
   CALL bakCOL(STFTMP,maxstiff,LOADEFF,ierror) !wk is calculated.

!  Obtain new VELocities, ACCelerations
   DO i=1,NEQ  
      oldVEL = VEL(i)
      oldACC = ACC(i)
      DISP(i)= wk(i)
      ACC(i) = a0*(DISP(i) - olddis(i)) -a2*oldVEL - a3*oldACC
      VEL(i) = oldVEL + a6*oldACC + a7*ACC(i)
   ENDDO

   write(*,22)'wk,load in transient:', 7, wk(7),load(7)

 22 format(1x,7(g12.5,1x))

   struc_vel_iter(:,:) = struc_vel(:,:)
   print *,'vel_iter, disp_iter:'
   print *, 'max/min struc_vel(,2) at 1 =', MAXVAL(struc_vel(:,2)),MINVAL(struc_vel(:,2))
   print *, 'max/min struc_disp(,2) at 1 =', MAXVAL(struc_disp(:,2)),MINVAL(struc_disp(:,2))

   DO i=1,nPtsBodyMarker(iBody)*nodeDoF
      node=(i+nodeDoF-1)/nodeDoF
      jdof=i-(node-1)*nodeDoF  !jdof 1~3 translation; 4~5 rotation
      if (jbc(i) == 0) then
         struc_disp(node,jdof) = 0.0
         struc_vel(node,jdof) = 0.0    !may be problem
      else
         struc_disp(node,jdof) = DISP(jbc(i))
         struc_vel(node,jdof)  = VEL(jbc(i))
      endif
   ENDDO

   print *, 'max/min struc_vel(,2) at 1 =', MAXVAL(struc_vel(:,2)),MINVAL(struc_vel(:,2))
   print *, 'max/min struc_disp(,2) at 1 =', MAXVAL(struc_disp(:,2)),MINVAL(struc_disp(:,2))

   idis = 179

   OPEN (idis,FILE='dis.out')
   WRITE (idis, 1100)
 1100  FORMAT(1x,a50,/, &
      '  node    x-disp       y-disp       z-disp ',  &
      '      x-rot        y-rot        z-rot' )

   DO i=1, nPtsBodyMarker(iBody)
      write(idis,1110) i, (struc_disp(i,j),j=1,6)
   ENDDO

 1110  FORMAT(1x,i5,2x,6(g12.6,1x))


   END SUBROUTINE FEA_TRANSIENT


!  ---------------------------------------------------------------------------------
   SUBROUTINE AxBCOL(matrix,vecin,vecout)
!  AxB product for COLumn storage

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE flow_parameters, ONLY : nPtsMax
   USE fea_unstructure_surface

   IMPLICIT NONE
   
   REAL(KIND=CGREAL) :: matrix(MAXSTIFF),vecin(nPtsMax*nodeDoF)
   REAL(KIND=CGREAL) :: vecout(NEQ)

   REAL(KIND=CGREAL) :: val,valmat

   INTEGER :: i,j,is,jlim,io,iloc

   do 10 i=1,NEQ
            jlim=max(1,(i-iprof(i)+1))
            do 20 j=jlim,i
               is=i
               io=i-j+1
               if (io .gt. iprof(is)) goto 20
               iloc = nloc(is) + io -1
               valmat=matrix(iloc)
               val = vecin(j)
               vecout(i)=vecout(i) + val*valmat
 20         continue

            jlim=min(iprof2(i),(neq-i+1))
            do 30 j=2,jlim
               is=i+j-1
               io=j
               if (io .gt. iprof(is)) goto 30
               iloc = nloc(is) + io -1
               valmat=matrix(iloc)
               val = vecin(i+j-1)
               vecout(i)=vecout(i) + val*valmat     
 30         continue
 10 continue

   END SUBROUTINE AxBCOL



