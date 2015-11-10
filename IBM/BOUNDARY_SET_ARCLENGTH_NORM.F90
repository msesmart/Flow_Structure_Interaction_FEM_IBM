!------------------------------------------------------------------------------
   SUBROUTINE calculate_arclength_norm_ds()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... Loop variables
    INTEGER           :: iBody,m, mMax, mMin, cNode, cElement
    REAL(KIND=CGREAL) :: tangX, tangY, distNode, dMinUnstruc
    REAL(KIND=CGREAL) :: normMag, tangMag
    REAL(KIND=CGREAL) :: vectorx,vectory,vectorz,vectProduct
    REAL(KIND=CGREAL) :: triElemVectorX(2),triElemVectorY(2),triElemVectorZ(2)
    REAL(KIND=CGREAL) :: dxMin,dyMin,dzMin,cellAreaMin,triElemAreaAv,bodyResolutionNormalized


!... Local variables

! check direction of surface normal and adjust ordering of triangle if normal does not point into body
! this is done only during fresh start

    IF (ntime == 0 .AND. nread == 0) THEN
      DO iBody = 1, nBody
        SELECT CASE( canonical_body_type(iBody) )
          CASE( ELLIPSOID ) 
           IF (boundary_formulation == GCM_METHOD) THEN
             PRINT*,'CANNOT COMPUTE SURFACE QUANTITIES FOR "ELLIPSOID"'
             PRINT*,'ABORTING RUN' 
             STOP
           ENDIF
          CASE(ELLIPTIC_CYLINDER,GENERAL_CYLINDER,UNSTRUCTURED_SURFACE)
            normDirFlag = oned
            PRINT*,'CALL determine_norm_dir_unstruc(iBody)'
            CALL determine_norm_dir_unstruc(iBody)
        END SELECT ! canonical_body_type
       ENDDO ! iBody
    ENDIF
 
    PRINT*,'Computing Surface Quantities for Unstructured surface'
    sBodyMarker = zero
    surfArea    = zero
    
! Select appropriate body type

!    PRINT*,'SETTING UP CANONICAL BODIES: Arc Length '

    DO iBody = 1, nBody

      SELECT CASE( canonical_body_type(iBody) )

        CASE( ELLIPSOID ) 
  
         IF (boundary_formulation == GCM_METHOD) THEN
           PRINT*,'CANNOT COMPUTE SURFACE QUANTITIES FOR "ELLIPSOID"'
           PRINT*,'ABORTING RUN' 
           STOP
	 ENDIF

        CASE(ELLIPTIC_CYLINDER,GENERAL_CYLINDER,UNSTRUCTURED_SURFACE)
           

!--Sweep surface of each triangular element
          surfArea(iBody)  = zero
                                    
          DO m=1,totNumTriElem(iBody)
                                    
!--First vector of each element
           triElemVectorx(1)= xBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                             -xBodyMarker(iBody,triElemNeig(iBody,1,m))
           triElemVectory(1)= yBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                             -yBodyMarker(iBody,triElemNeig(iBody,1,m))
           triElemVectorz(1)= zBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                             -zBodyMarker(iBody,triElemNeig(iBody,1,m))

!--Second vector of each element
           triElemVectorx(2)= xBodyMarker(iBody,triElemNeig(iBody,3,m)) &
                             -xBodyMarker(iBody,triElemNeig(iBody,2,m))
           triElemVectory(2)= yBodyMarker(iBody,triElemNeig(iBody,3,m)) &
                             -yBodyMarker(iBody,triElemNeig(iBody,2,m))
           triElemVectorz(2)= zBodyMarker(iBody,triElemNeig(iBody,3,m)) &
                             -zBodyMarker(iBody,triElemNeig(iBody,2,m))

!--Normal vector of the element
           triElemNormx(iBody,m)= (  triElemVectory(1)*triElemVectorz(2)  &
                                   - triElemVectorz(1)*triElemVectory(2)  )
           triElemNormy(iBody,m)= (  triElemVectorz(1)*triElemVectorx(2)  &
                                   - triElemVectorx(1)*triElemVectorz(2)  )
           triElemNormz(iBody,m)= (  triElemVectorx(1)*triElemVectory(2)  &
                                   - triElemVectory(1)*triElemVectorx(2)  )
!           IF (body_dim(iBody) == BODY_DIM2) THEN
!
!              IF( ABS(triElemNormz(iBody,m)) > 1.0E-6_CGREAL ) THEN
!                IF( ABS(triElemNormz(iBody,m)) > 1.0E-4_CGREAL ) THEN
!                 PRINT*,'Body:Element = ',ibody,m
!                 PRINT*, triElemNormx(iBody,m),triElemNormy(iBody,m),triElemNormz(iBody,m)
!                 PRINT*,'2D Body Element has normal in the z-direction'
!                 PRINT*,'Aborting'
!                ELSE
!                 PRINT*,'Body:Element = ',ibody,m
!                 PRINT*, triElemNormx(iBody,m),triElemNormy(iBody,m),triElemNormz(iBody,m)
!                 PRINT*,'Zeroing out the small z-normal'
!                 triElemNormz(iBody,m) = zero
!                ENDIF
!              ENDIF
! 
!	   ENDIF
                 
           normMag  = SQRT(  triElemNormx(iBody,m)**2   &
                           + triElemNormy(iBody,m)**2   &
                           + triElemNormz(iBody,m)**2     )

! Area of element
           triElemArea(iBody,m)  = half*normMag
           surfArea(iBody)       = surfArea(iBody) + triElemArea(iBody,m)

! Unit Normal vector
           triElemNormx(iBody,m) = triElemNormx(iBody,m)/normMag
           triElemNormy(iBody,m) = triElemNormy(iBody,m)/normMag
           triElemNormz(iBody,m) = triElemNormz(iBody,m)/normMag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	  write(876,123)m, triElemNormx(iBody,m),triElemNormy(iBody,m),triElemNormz(iBody,m)
!123	format(i4,3(2x,e14.7))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit Tangents
! Tangent-2 defined parallel to vector from vertex-1 to vertex-2

           triElemTang2x(iBody,m)  = xBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                                    -xBodyMarker(iBody,triElemNeig(iBody,1,m))
           triElemTang2y(iBody,m)  = yBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                                    -yBodyMarker(iBody,triElemNeig(iBody,1,m))
           triElemTang2z(iBody,m)  = zBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                                    -zBodyMarker(iBody,triElemNeig(iBody,1,m))

           tangMag     = SQRT(  triElemTang2x(iBody,m)**2   &
                              + triElemTang2y(iBody,m)**2   &
                              + triElemTang2z(iBody,m)**2 )

           triElemTang2x(iBody,m)  = triElemTang2x(iBody,m)/tangMag
           triElemTang2y(iBody,m)  = triElemTang2y(iBody,m)/tangMag
           triElemTang2z(iBody,m)  = triElemTang2z(iBody,m)/tangMag

! t1 = t2 x n

           triElemTang1x(iBody,m)  =  triElemTang2y(iBody,m)*triElemNormz(iBody,m)  &
                                    - triElemTang2z(iBody,m)*triElemNormy(iBody,m)
           triElemTang1y(iBody,m)  =- triElemTang2x(iBody,m)*triElemNormz(iBody,m)  &
                                    + triElemTang2z(iBody,m)*triElemNormx(iBody,m)
           triElemTang1z(iBody,m)  =  triElemTang2x(iBody,m)*triElemNormy(iBody,m)  &
                                    - triElemTang2y(iBody,m)*triElemNormx(iBody,m)

!--Centroid of the each element
           triElemCentx(iBody,m)=(xBodyMarker(iBody,triElemNeig(iBody,1,m)) &
                                 +xBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                                 +xBodyMarker(iBody,triElemNeig(iBody,3,m)))/3.0_CGREAL
           triElemCenty(iBody,m)=(yBodyMarker(iBody,triElemNeig(iBody,1,m)) &
                                 +yBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                                 +yBodyMarker(iBody,triElemNeig(iBody,3,m)))/3.0_CGREAL
           triElemCentz(iBody,m)=(zBodyMarker(iBody,triElemNeig(iBody,1,m)) &
                                 +zBodyMarker(iBody,triElemNeig(iBody,2,m)) &
                                 +zBodyMarker(iBody,triElemNeig(iBody,3,m)))/3.0_CGREAL

          ENDDO ! m

          dxMin                    = MINVAL(dx(1:nx-1))
          dyMin                    = MINVAL(dy(1:ny-1))
          dzMin                    = MINVAL(dz(1:nz-1))
          cellAreaMin              = (dxMin*dyMIn*dzMin)**(2.0_CGREAL/3.0_CGREAL)
          triElemAreaAv            = surfArea(iBody)/totNumTriElem(iBody)
          bodyResolutionNormalized = cellAreaMin/triElemAreaAv
          print*,'Min Cell Area                             = ',cellAreaMin
          print*,'Surface area of body',ibody,            ' = ',surfArea(iBody)
          print*,'Average Element Area for',ibody,        ' = ',triElemAreaAv
          print*,'Norm. Surface resolution of body',ibody,' = ',bodyResolutionNormalized

      END SELECT ! canonical_body_type

    ENDDO ! iBody


   END SUBROUTINE calculate_arclength_norm_ds
!
!
!------------------------------------------------------------------------------
! Subroutine determines if for a given unstructured surface mesh,
! the norm points in or out and also reorders vertex numbers

! Note: our convention is that positive normal vector points into body
! Per convention of solid geometry, if three vertices (1,2,3) of a triangle are
! ordered clockwise when viewed from one side, then (p2-p1)x(p3-p1) produces a vector
! that points out towards the direction from where it is being viewed. If we force
! the normal to point inwards, then we should also change the vertex ordering
!
!------------------------------------------------------------------------------

   SUBROUTINE determine_norm_dir_unstruc(iBody)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER,INTENT(IN)::iBody

!... Loop variables
    INTEGER           :: m,cNode, cElement

!... Local variables
    INTEGER           :: node2,node3
    REAL(KIND=CGREAL) :: distNode, dMinUnstruc
    REAL(KIND=CGREAL) :: vectorx,vectory,vectorz,vectProduct
    REAL(KIND=CGREAL) :: triElemVectorX(2),triElemVectorY(2),triElemVectorZ(2)
    REAL(KIND=CGREAL) :: cTriElemNormx,cTriElemNormy,cTriElemNormz
    REAL(KIND=CGREAL) :: cTriElemCentx,cTriElemCenty,cTriElemCentz


    cNode = 0

!---Find node closest to outside point
    dMinUnstruc = 1.0E8_CGREAL
    DO m=1,nPtsBodyMarker(iBody)
       distNode= SQRT( (xBodyMarker(iBody,m)-pointOutsideBodyX(iBody))**2  &
                      +(yBodyMarker(iBody,m)-pointOutsideBodyY(iBody))**2  &
                      +(zBodyMarker(iBody,m)-pointOutsideBodyZ(iBody))**2   )
       IF(distNode <= dMinUnstruc) THEN
         dMinUnstruc  = distNode
         cNode        = m
       ENDIF
    ENDDO

!---Find element corresponing to closest node
    DO m=1,totNumTriElem(iBody)
       IF ( triElemNeig(iBody,1,m) == cNode .OR. &
            triElemNeig(iBody,2,m) == cNode .OR. &
            triElemNeig(iBody,3,m) == cNode       ) cElement = m
    ENDDO

!RRRRRRRRRRRRRRRRRRRRRRRRRRR
    print*,'Closest Node to outside point            = ',cNode
    print*,'Element Corresponding to closest node is = ',cElement
!RRRRRRRRRRRRRRRRRRRRRRRRRRR


!     1 *-------------* 3
!        \           /
!         \         /
!          \       /
!           \     /
!            \   /
!             \ /
!            2 *
!
!--Sweep surface of CLOSEST element to determine normal direction of triangular elements

    triElemVectorx(1)= xBodyMarker(iBody,triElemNeig(iBody,2,cElement)) &
                     - xBodyMarker(iBody,triElemNeig(iBody,1,cElement))
    triElemVectory(1)= yBodyMarker(iBody,triElemNeig(iBody,2,cElement)) &
                     - yBodyMarker(iBody,triElemNeig(iBody,1,cElement))
    triElemVectorz(1)= zBodyMarker(iBody,triElemNeig(iBody,2,cElement)) &
                     - zBodyMarker(iBody,triElemNeig(iBody,1,cElement))

    triElemVectorx(2)= xBodyMarker(iBody,triElemNeig(iBody,3,cElement)) &
                     - xBodyMarker(iBody,triElemNeig(iBody,2,cElement))
    triElemVectory(2)= yBodyMarker(iBody,triElemNeig(iBody,3,cElement)) &
                     - yBodyMarker(iBody,triElemNeig(iBody,2,cElement))
    triElemVectorz(2)= zBodyMarker(iBody,triElemNeig(iBody,3,cElement)) &
                     - zBodyMarker(iBody,triElemNeig(iBody,2,cElement))

!-- Normal of closest element
    cTriElemNormx = ( triElemVectory(1)*triElemVectorz(2)  &
                    - triElemVectorz(1)*triElemVectory(2) )
    cTriElemNormy = ( triElemVectorz(1)*triElemVectorx(2)  &
                    - triElemVectorx(1)*triElemVectorz(2) )
    cTriElemNormz = ( triElemVectorx(1)*triElemVectory(2)  &
                    - triElemVectory(1)*triElemVectorx(2) )

!--Centroid of the closest element
    cTriElemCentx = (xBodyMarker(iBody,triElemNeig(iBody,1,cElement)) &
                   + xBodyMarker(iBody,triElemNeig(iBody,2,cElement)) &
                   + xBodyMarker(iBody,triElemNeig(iBody,3,cElement)))/3.0_CGREAL
    cTriElemCenty = (yBodyMarker(iBody,triElemNeig(iBody,1,cElement)) &
                   + yBodyMarker(iBody,triElemNeig(iBody,2,cElement)) &
                   + yBodyMarker(iBody,triElemNeig(iBody,3,cElement)))/3.0_CGREAL
    cTriElemCentz = (zBodyMarker(iBody,triElemNeig(iBody,1,cElement)) &
                   + zBodyMarker(iBody,triElemNeig(iBody,2,cElement)) &
                   + zBodyMarker(iBody,triElemNeig(iBody,3,cElement)))/3.0_CGREAL

    vectorx = cTriElemCentX - pointOutsideBodyX(iBody)
    vectory = cTriElemCentY - pointOutsideBodyY(iBody)
    vectorz = cTriElemCentZ - pointOutsideBodyZ(iBody)

    vectProduct = vectorx*cTriElemNormx  &
                + vectory*cTriElemNormy  &
                + vectorz*cTriElemNormz

!RRRRRRRRRRRRRRRRRRRRRRRRRRRRr
    print*,'vectProduct = ',vectproduct
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRr

    normDirFlag = oned
    IF (vectProduct < zero)THEN                 
        normDirFlag =-oned
        print*,'Reordering triangle vertices since normal points out'
        DO m=1,totNumTriElem(iBody)         ! changing vertex ordering
            node2 = triElemNeig(iBody,2,m)
            node3 = triElemNeig(iBody,3,m)
            triElemNeig(iBody,2,m) = node3
            triElemNeig(iBody,3,m) = node2
        ENDDO
    ENDIF


   END SUBROUTINE determine_norm_dir_unstruc
!------------------------------------------------------------------------------

