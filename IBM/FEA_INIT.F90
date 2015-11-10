 SUBROUTINE fea_initial

   USE fea_unstructure_surface
   USE unstructured_surface_arrays
   USE flow_parameters
   USE boundary_arrays, ONLY : uBodyMarker,vBodyMarker,wBodyMarker,xBodyMarker,yBodyMarker,zBodyMarker
   
   IMPLICIT NONE

   REAL(KIND=CGREAL) :: TEMP(10)

   INTEGER :: nTriElemMax
   INTEGER :: XDOF,YDOF,ZDOF,XROT,YROT,ZROT  ! 0 or 1

   INTEGER :: ibody,iNode,iBC
   INTEGER :: I,J,ind,itmp
     
   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: LOADold   
 
        print *, 'in fea_init'  
   ALLOCATE(pBodyMarker(nBody,nPtsMax))
   ALLOCATE(pxBodyMarker(nBody,nPtsMax))
   ALLOCATE(pyBodyMarker(nBody,nPtsMax))
   ALLOCATE(pzBodyMarker(nBody,nPtsMax))

!   ALLOCATE(nPrescribedMarker(nBody))
!   ALLOCATE(PrescribedMarker(nBody,nPtsMax/10))

   ALLOCATE(PROPERTY(2,10))        !At most 2 kinds of materials
   ALLOCATE(JBC(nPtsMax*nodeDoF))  !Save the DoF at each node

   nTriElemMax = MAXVAL(totNumTriElem(:))
   ALLOCATE(triElemP(nBody,nTriElemMax))

   ALLOCATE(LOAD(nPtsMax*nodeDoF))
   ALLOCATE(dLOAD(nPtsMax*nodeDoF))
   ALLOCATE(LOADOLD(nPtsMax*nodeDoF))

   ALLOCATE(WK(nPtsMax*nodeDoF))
   ALLOCATE(DISP(nPtsMax*nodeDoF))
   ALLOCATE(VEL(nPtsMax*nodeDoF))
   ALLOCATE(ACC(nPtsMax*nodeDoF))

   ALLOCATE(struc_disp(nPtsMax,nodeDoF))
   ALLOCATE(struc_vel(nPtsMax,nodeDoF))
   ALLOCATE(struc_disp_iter(nPtsMax,nodeDoF))
   ALLOCATE(struc_vel_iter(nPtsMax,nodeDoF))
   ALLOCATE(struc_olddisp(nPtsMax,nodeDoF))

   ALLOCATE(struc_disp_worelax(nPtsMax,nodeDoF))
   ALLOCATE(struc_vel_worelax(nPtsMax,nodeDoF))

   ALLOCATE(xBodyMarker0(nBody,nPtsMax))
   ALLOCATE(yBodyMarker0(nBody,nPtsMax))
   ALLOCATE(zBodyMarker0(nBody,nPtsMax))

   nTriElemMax = MAXVAL(totNumTriElem(:))
   
   ALLOCATE(ELMatType(nTriElemMax))

   ALLOCATE(IPROF(nPtsMax*nodeDoF))
   ALLOCATE(IPROF2(nPtsMax*nodeDoF))
   ALLOCATE(NLOC(nPtsMax*nodeDoF))
   ALLOCATE(STF((nPtsMax*nodeDoF)**2))
   ALLOCATE(MSS((nPtsMax*nodeDoF)**2))
  
   ALLOCATE(yBodyMarkerOld(nBody,nPtsMax))
  
   xBodyMarker0(:,:) = xBodyMarker(:,:) 
   yBodyMarker0(:,:) = yBodyMarker(:,:) 
   zBodyMarker0(:,:) = zBodyMarker(:,:) 

   uBodyMarker(:,:) = 0
   vBodyMarker(:,:) = 0
   wBodyMarker(:,:) = 0
   
   struc_disp(:,:) = 0
   struc_olddisp(:,:) = 0
   struc_vel(:,:) = 0
   struc_disp_iter(:,:) = 0
   struc_vel_iter(:,:) = 0
   
   yBodyMarkerOld(:,:) = yBodyMarker(:,:)
   
!  ------------------------- Read in fea_input.dat -----------------------------!
   READ(ifuFEA_Input,*) !About the input file.
   READ(ifuFEA_Input,*) TRNS
   READ(ifuFEA_Input,*) ELETYPE
   READ(ifuFEA_Input,*) NMATP !Total types of material
   READ(ifuFEA_Input,*) !Materail_type  E G  Thickness  Density  Plane_stress/strain
   DO I=1,NMATP
      READ (ifuFEA_Input,*) itmp, (PROPERTY(I,J),J=1,7)
   ENDDO
 
   massratio = PROPERTY(1,7)/PROPERTY(1,4)*PROPERTY(1,6)/PROPERTY(1,3)  !only massratio of material 1 is used.
   write (*,*) 'massratio is:', massratio

!  --------------------- Initialization of boundary condition -------------------------!
!  0=fixed;    1=free (default)
   DO I = 1,nPtsMax*nodeDoF
      JBC(I) = 1
   ENDDO

   READ (ifuFEA_Bound_Cond,*)
   READ (ifuFEA_Bound_Cond,*) NBC  ! Number of nodes given explicit DoF constraint.
   DO iBC = 1, NBC
      READ (ifuFEA_Bound_Cond,*) iNode,XDOF,YDOF,ZDOF,XROT,YROT,ZROT
      !print *, iNode,XDOF,YDOF,ZDOF,XROT,YROT,ZROT
      ind = iNode*6
      JBC(ind-5) = XDOF
      JBC(ind-4) = YDOF
      JBC(ind-3) = ZDOF
      JBC(ind-2) = XROT
      JBC(ind-1) = YROT
      JBC(ind)   = ZROT
   ENDDO

   iBody = 1 ! Wanh set iBody temporarily
   nPrescribedMarker(iBody) = NBC  !set prescribed marker same as points with
                                   !constrains

   CALL DATACNV(iBody)

!  Compute locations of cols in vector
   NLOC(1)=1
   DO I=1,NEQ-1
      NLOC(i+1) = NLOC(i) + IPROF(i)
   ENDDO

   MAXSTIFF = NLOC(NEQ)+IPROF(NEQ)



!   CALL fea_readmesh
!   CALL fea_readconload
!   CALL fea_readcorrespondancetable

!   DO ibody = 1, fea_nbody
!     CALL fea_assemble(ibody)
!    CALL fea_dynamicsetup(ibody)
!    CALL fea_contactsetup(ibody)
!   END DO !ibody

!   IF(nread == 1) CALL fea_readrestart()
!   CALL fea_read_probe_inputs()
!   CALL fea_open_probe_files()

 END SUBROUTINE fea_initial


!-----------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------!
   SUBROUTINE DATACNV(iBody)

!  DATA CoNVerter

   USE global_parameters 
   USE flow_parameters, ONLY : nPtsBodyMarker
   USE boundary_arrays    
   USE unstructured_surface_arrays
   USE fea_unstructure_surface
     
   INTEGER :: iBody

   INTEGER :: npk(totNumTriElem(iBody))
   INTEGER :: I


!  Impose conditions for special global shapes
!  CALL SPCLBC(1,iBody)
!   CALL SPCLBC(2,iBody)  !original, force applied in z-drection
   CALL SPCLBC(4,iBody)   !force applied in y-direction

!  Assigning the equation numbers in JBC for non-zero DoFs 
!  From 1 up; only non-zero given a number

   NEQ=0

   DO I=1,nPtsBodyMarker(iBody)*nodeDoF
      IF (JBC(I) > 0) THEN
         NEQ = NEQ + 1
         JBC(I) = NEQ  
      ELSE
         JBC(I) = 0
      ENDIF
   ENDDO
!  Now JBC saves the numbering of equations (corresponding to body markers)
!  JBC(I) = 0, means certain DoF at certain marker is fixed.
!  NEQ is the total number of equations (unknown displacements)


!   DO I= 1, nPtsBodyMarker(iBody)
!      ind = (I-1)*6
!      if ( JBC(ind+1).ne.0 )   LOAD(JBC(ind+1)) = pxBodyMarker(iBody,I) 
!      if ( JBC(ind+2).ne.0 )   LOAD(JBC(ind+2)) = pyBodyMarker(iBody,I)
!      if ( JBC(ind+3).ne.0 )   LOAD(JBC(ind+3)) = pzBodyMarker(iBody,I)
!      if ( JBC(ind+4).ne.0 )   LOAD(JBC(ind+4)) = 0  !No external moment applied
!      if ( JBC(ind+5).ne.0 )   LOAD(JBC(ind+5)) = 0
!      if ( JBC(ind+6).ne.0 )   LOAD(JBC(ind+6)) = 0
!      if (i == 2) then
!         print *, 'i and LOAD(JBC(ind+3)):', i,LOAD(JBC(ind+1)),LOAD(JBC(ind+3))
!         print *, 'jbc(ind+1),2,3:', JBC(ind+1),JBC(ind+2),JBC(ind+3)
!         print *, 'pzbodymark', pzBodyMarker(iBody,I)
!      endif
!   ENDDO

   print *,' max load:', maxval(LOAD),maxloc(LOAD)
        
!  Calcualte IBAND
   CALL MAXBND(iBody)
   write(*,*)'@@ from datain,  half band =',iband 

!  Compute storage
   iloc=0
   do i=1,neq
      iloc=iloc+IPROF(i)
   enddo
   ipercent=(iloc*100)/(NEQ*IBAND)
!  write(ilog,1022) neq,iband,neq*iband,iloc,ipercent

!  STORE results for future reference, cal other profile
   do i=1,neq
      ibandh=1
      iend=i+IBAND-1
      if (iend .gt. neq) iend=neq
      do j=i+1,iend
                ibandv=IPROF(j)
                ji1=j-i+1
                if (ibandv .ge. ji1) then
                    ibandh = ji1
                endif
      enddo
      iprof2(i)=ibandh
   enddo

!  Compute total mass
   DO i= 1,totNumTriElem(iBody)
      i1 = triElemNeig(iBody,1,i)
      j1 = triElemNeig(iBody,2,i)
      k1 = triElemNeig(iBody,3,i)

      x1 = xBodyMarker(iBody,i1)
      x2 = xBodyMarker(iBody,j1)
      x3 = xBodyMarker(iBody,k1)

      y1 = yBodyMarker(iBody,i1)
      y2 = yBodyMarker(iBody,j1)
      y3 = yBodyMarker(iBody,k1)

      z1 = zBodyMarker(iBody,i1)
      z2 = zBodyMarker(iBody,j1)
      z3 = zBodyMarker(iBody,k1)

!      rh0 = PROPERTY(ELMatType(i),4)
      rh0 = PROPERTY(1,4)   !set material type to be 1 temporarily \\Wanh

      IF (ELETYPE .EQ. 3) THEN
!         a0 = PROPERTY(ELMatType(i),3)
         a0 = PROPERTY(1,3) !set material type to be 1 temporarily \\Wanh
         zlen = sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
         vol = zlen*a0
      ELSEIF (ELETYPE .EQ. 4) THEN
!         tt0 = PROPERTY(ELMatType(i),3)
         tt0 = PROPERTY(1,3)  !set material type to be 1 temporarily \\Wanh
!        Determine vector area
         axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
         ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
         azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
         area=sqrt( axy*axy + ayz*ayz + azx*azx)
         vol = area*tt0
      ENDIF

      zmass=vol*rh0
   ENDDO

!!  Compute resultants
!   xres=0.0
!   yres=0.0
!   zres=0.0
!   xmres=0.0
!   ymres=0.0
!   zmres=0.0
!
!   do i=1,nPtsBodyMarker(iBody)
!      xres=xres+pxBodyMarker(iBody,i)
!      yres=yres+pyBodyMarker(iBody,i)
!      zres=zres+pzBodyMarker(iBody,i)
!!      xmres=xmres+xmom(i)
!!      ymres=ymres+ymom(i) 
!!      zmres=zmres+zmom(i)  
!   enddo
!
!   write(*,*) '@@RESULTants:'
!   write(*,*) '@@force: ', xres,yres,zres
!   write(*,*) '@@moment:', xmres,ymres,zmres


   END SUBROUTINE DATACNV

!-----------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------!
   SUBROUTINE MAXBND(iBody)
!  MAX BaNDwidth calculation

   USE global_parameters
   USE flow_parameters, ONLY : nPtsBodyMarker
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE
   INTEGER :: iBody
   INTEGER :: IHFBND
   INTEGER :: i,n,idof,jdof,kdof,ieqn1,ieqn2,ieq,J,JBAND,IEQ2
   INTEGER :: ipv(18)

!  nodeDoF is defined as 6 in AMODULES.f90
   DO n=1,nPtsBodyMarker(iBody)*nodeDoF
      IPROF(n)=0
   ENDDO

   IHFBND=0

   DO 100 n=1,totNumTriElem(iBody)

      idof=(triElemNeig(iBody,1,n)-1)*nodeDoF
      jdof=(triElemNeig(iBody,2,n)-1)*nodeDoF
      kdof=(triElemNeig(iBody,3,n)-1)*nodeDoF

      DO i=1,nodeDoF
         ipv(i  )=idof+i
         ipv(i+nodeDoF)=jdof+i
         ipv(i+nodeDoF*2)=kdof+i
      ENDDO

      DO I=1,18
         IEQN1=JBC(IPV(I))
         IF (IEQN1.GT.0) THEN
            DO J=I,18
               IEQN2=JBC(IPV(J))
               IF (IEQN2.GT.0) THEN
                  IHFBND = MAX0(IHFBND,IABS(IEQN1-IEQN2))
                  JBAND=ABS(IEQN1-IEQN2)+1
                  IEQ=MAX(IEQN1,IEQN2)
                  IF (JBAND .GT. IPROF(IEQ)) THEN
                     IPROF(IEQ)=JBAND
                  ENDIF
                  IEQ2=MIN(IEQN1,IEQN2)
!                 if (jband .gt. iprof2(ieq2)) then
!                    iprof2(ieq2)=jband
!                 endif
               ENDIF
            ENDDO ! END LOOP OF J
         ENDIF
      ENDDO

 100 CONTINUE

   IBAND=IHFBND+1

   END SUBROUTINE MAXBND 


!-----------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------!
   SUBROUTINE SPCLBC(IGLOBAL,iBody)
!  Impose SPeCiaL global BCs

   USE global_parameters
   USE flow_parameters, ONLY : nPtsBodyMarker
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   INTEGER :: IGOLBAL,iBody

   INTEGER :: I

   IF (iglobal .eq. 1) THEN
!     flat 2-D in-plane with rot plate 
!     dof: u,v,rotZ 
      DO i=1,nPtsBodyMarker(iBody)
         jbc(i*6-1) = 0
         jbc(i*6-2) = 0
         jbc(i*6-3) = 0
      ENDDO
   ELSEIF (iglobal .eq. 2) THEN
!     flat 2-D flex plate 
!     dof: w,rotX,rotY
      DO i=1,nPtsBodyMarker(iBody)
         jbc(i*6-0) = 0  ! rotz = 0
         jbc(i*6-4) = 0  ! v=0
         jbc(i*6-5) = 0  ! u=0
      ENDDO     
   ELSEIF (iglobal .eq. 3) THEN
!     flat 2-D flex plate with in-plane loading
      DO i=1,nPtsBodyMarker(iBody)
         jbc(i*6-0) = 0  !rotz = 0
      ENDDO
   ELSEIF (iglobal .eq. 4) THEN
!     Flex plate with in-plane loading, force is applied in y-direction
      DO i=1,nPtsBodyMarker(iBody)
         jbc(i*6-1) = 0  ! roty = 0
         jbc(i*6-3) = 0  ! w=0
         jbc(i*6-5) = 0  ! u=0
      ENDDO
   ENDIF

   END SUBROUTINE SPCLBC
