   SUBROUTINE fea_getpressure
!  Extract pressure from fluid to solid.
!  iblankm is used in this version (Get-pressure from Dr. Dong) 

   USE flow_arrays
   USE flow_parameters
   USE pressure_arrays
   USE grid_arrays
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface
   USE GCM_arrays

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: triElemPp, triElemPm
   REAL(KIND=CGREAL) :: pmin, pmax, pdiffmin, pdiffmax, areamin, areamax
   REAL(KIND=CGREAL) :: fluidpmax, fluidpmin
   REAL(KIND=CGREAL) :: distMin,dist
   REAL(KIND=CGREAL) :: alphaX,alphaY,alphaZ

   REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:) :: pBodyMarker_1, pBodyMarker_2

   INTEGER :: Ipmin,Jpmin,Kpmin,Ipmax,Jpmax,Kpmax
   INTEGER :: ii,jj,kk

   REAL(KIND=CGREAL)    :: xp,yp,zp,pc_x,pc_y,pc_z,costheta
   REAL(KIND=CGREAL)    :: totlift,totx,totz,tot,totarea
   INTEGER :: I,J,K, iBody, iB,m,ie,im,ind,n
   INTEGER :: Ic,Jc,Kc

   INTEGER :: nG,iG,jG,kG

   CHARACTER*100 :: tmpchar
   CHARACTER*1 :: charid

   LOGICAL :: debug_getpressure

   debug_getpressure = .false.
!   debug_getpressure = .true.

   nPtsMax = MAXVAL(nPtsBodyMarker(:))

   ALLOCATE(pBodyMarker_1(nBody, nPtsMax), pBodyMarker_2(nBody, nPtsMax))

   DO iBody = 1,nBody

     DO m = 1,nPtsBodyMarker(iBody)

     ! find closest ghost node
        distMin = 1.0E8_CGREAL
        DO n=1,nGhost
           i = iGhostP(n)
           j = jGhostP(n)
           k = kGhostP(n)

           IF (iblank_memb(i,j,k) == 1) THEN
!           IF (ghostCellMemb(i,j,k) == 1) THEN

                 dist = (xc(i)-xBodyMarker(iBody,m))**2 &
                     +(yc(j)-yBodyMarker(iBody,m))**2 &
                     +(zc(k)-zBodyMarker(iBody,m))**2

                 IF ( dist <= distMin ) THEN
                  distMin = dist
                  nG  = n
                  iG  = i
                  jG  = j
                  kG  = k
                 ENDIF

           ENDIF
        ENDDO

        pBodyMarker_1(ibody,m) = p(iG, jG, kG)


  ! find closest ghost node
        distMin = 1.0E8_CGREAL
        DO n=1,nGhost
           i = iGhostP(n)
           j = jGhostP(n)
           k = kGhostP(n)

           IF (iblank_memb(i,j,k) == 0) THEN
!           IF (ghostCellMark_memb(i,j,k) == 0) THEN

              dist = (xc(i)-xBodyMarker(iBody,m))**2 &
                    +(yc(j)-yBodyMarker(iBody,m))**2 &
                    +(zc(k)-zBodyMarker(iBody,m))**2

              IF ( dist <= distMin ) THEN
                  distMin = dist
                  nG  = n
                  iG  = i
                  jG  = j
                  kG  = k
              ENDIF

           END IF

        ENDDO

       pBodyMarker_2(ibody,m) = p(iG, jG, kG)

     END DO ! end for m


!      Do m=1,totNumTriElem(iBody)

!         xp = triElemCentx(iBody,m)
!         yp = triElemCenty(iBody,m)
!         zp = triElemCentz(iBody,m)
        
!         CALL LOCATION(xp,yp,zp,Ic,Jc,Kc)

!         fluidpmin = 1E10_CGREAL
!         fluidpmax = -1E10_CGREAL

!         Do k= Kc-1,Kc+1
!         Do j= Jc-1,Jc+1
!         Do i= Ic-1,Ic+1

!            IF (fluidpmax < pPrime(i,j,k)) then
!               fluidpmax = pPrime(i,j,k)
!               Ipmax = i
!               Jpmax = j
!               Kpmax = k
!            ENDIF
!            IF (fluidpmin > pPrime(i,j,k)) then
!               fluidpmin = pPrime(i,j,k)
!               Ipmin = i
!               Jpmin = j
!               Kpmin = k
!            ENDIF 

!         ENDDO ! end i loop
!         ENDDO ! end j loop
!         ENDDO ! end k loop

!      Enddo ! end m loop

!     print *, 'Calling Filter ...'
!     After calling FILTER once, the pressure at center is distributed to apex
!     and saved in pBodyMarker.
!      call FILTER_pxyz (iBody,10)


      DO I= 1, nPtsBodyMarker(iBody)

         pyBodyMarker(iBody,I) = pBodyMarker_2(iBody,I)-pBodyMarker_1(iBody,I)

         ind = (I-1)*6
         if ( JBC(ind+1).ne.0 )   LOAD(JBC(ind+1)) = pxBodyMarker(iBody,I)*massratio
         if ( JBC(ind+2).ne.0 )   LOAD(JBC(ind+2)) = pyBodyMarker(iBody,I)*massratio 
         if ( JBC(ind+3).ne.0 )   LOAD(JBC(ind+3)) = pzBodyMarker(iBody,I)*massratio 
         if ( JBC(ind+4).ne.0 )   LOAD(JBC(ind+4)) = 0  !No external moment applied
         if ( JBC(ind+5).ne.0 )   LOAD(JBC(ind+5)) = 0
         if ( JBC(ind+6).ne.0 )   LOAD(JBC(ind+6)) = 0
!        for debug only
         if ( JBC(ind+1).ne.0 )   LOAD(JBC(ind+1)) = 0
         if ( JBC(ind+3).ne.0 )   LOAD(JBC(ind+3)) = 0
      ENDDO     

      print *, 'In getpressure, max/min pybodymarker:', maxval(pyBodyMarker), minval(pyBodyMarker)

!      IF (debug_getpressure) THEN
!          write (charid,"(I1.1)") iBody
!!          OPEN (21,FILE='wing'//trim(charid)//'pre.dat')
!          OPEN (21,FILE='platepressure.dat')
!!          write (21,*) ' VARIABLES= "X","Y","Z","P","Pm","Pp","Norm"'
!          write (21,*) ' VARIABLES= "X","Y","Z","P","Norm"'
!          write (21,*) ' ZONE T="wing pre " N=',nPtsBodyMarker(iBody), 'E=',totNumTriElem(iBody), & 
!!                       'DATAPACKING=BLOCK, VARLOCATION=([4,5,6,7]=CELLCENTERED), ZONETYPE=FEtriangle'
!                       'DATAPACKING=BLOCK, VARLOCATION=([4,5]=CELLCENTERED), ZONETYPE=FEtriangle'
!
!          write(21,100) (xBodyMarker(iBody,im),im=1,nPtsBodyMarker(iBody))
!          write(21,100) (yBodyMarker(iBody,im),im=1,nPtsBodyMarker(iBody))
!          write(21,100) (zBodyMarker(iBody,im),im=1,nPtsBodyMarker(iBody))
!
!          write(21,100) (triElemP(iBody,ie)*triElemNormy(iBody,ie),ie=1,totNumTriElem(iBody))
!!          write(21,100) (triElemPm(iBody,ie),ie=1,totNumTriElem(iBody))
!!          write(21,100) (triElemPp(iBody,ie),ie=1,totNumTriElem(iBody))
!          write(21,100) (triElemNormy(iBody,ie),ie=1,totNumTriElem(iBody))
!     
!          do ie=1,totNumTriElem(iBody)
!             write (21,*) triElemNeig(iBody,1,ie),triElemNeig(iBody,2,ie),triElemNeig(iBody,3,ie)
!          end do
!
!          close(21)
!
!        ENDIF  ! End of debug_getpressure
        
        IF (debug_getpressure) THEN
           print *, 'In getpressure, debug is on'
           write (tmpchar,"('pressure_platechord.',I7.7,'.',I3.3)") ntime,FSI_ITERATION
           open (22,file=trim(tmpchar))
           DO i=1,nPtsBodyMarker(iBody)
              write (22,122) xBodyMarker(iBody,i), pBodyMarker_2(iBody,i)-pBodyMarker_1(iBody,I),  &
                pBodyMarker_2(iBody,I),pBodyMarker_1(iBody,I)
           ENDDO
           CLOSE (22)
        ENDIF

   ENDDO ! loop of iBody

 100 format (20E18.6)
 122 format (4E18.6)

   DEALLOCATE(pBodyMarker_1, pBodyMarker_2)

   END SUBROUTINE fea_getpressure



   SUBROUTINE LOCATION(xp,yp,zp,i1,j1,k1)

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
   double precision :: xe(NP),ye(NP),ue(NP),ve(NP)
        
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

   i1 = 1
   i2 = nx

   nsearch = 0 
   do while ( (i2-i1)>1.1 ) 
      if (nsearch >maxsearch) then 
         write (*,*) 'max number of searching has reached in x'
         write (*,*) 'xp = ', xp
         stop
      endif
      if (xp>=x(i1) .and. xp<x((i1+i2)/2)) then
         i2 = (i1+i2)/2
      else if ( xp>=x((i1+i2)/2) .and. xp<x(i2) ) then
         i1 = (i1+i2)/2
      endif
      nsearch = nsearch+1
   enddo

   xe(1) = x(i1)
   xe(2) = x(i2)
   xe(3) = xe(2)
   xe(4) = xe(1)


   if (debug ==1 .and. abs(xp-xdebug)<1e-4 .and. abs(yp-ydebug)<1e-4) then
      print *, 'i1,i2,xe(1),xe(2):',i1,i2,xe(1),xe(2)
   endif

   j1 = 1
   j2 = ny

   nsearch = 0 
   do while ( (j2-j1)>1.1 ) 
      if (nsearch >maxsearch) then 
         write (*,*) 'max number of searching has reached in y'
         write (*,*) 'yp = ', yp
         stop
      endif
      if (yp>=y(j1) .and. yp<y((j1+j2)/2)) then
         j2 = (j1+j2)/2
      else if (yp>=y((j1+j2)/2) .and. yp<y(j2)) then
         j1 = (j1+j2)/2
      endif
      nsearch = nsearch+1
   enddo

   ye(1) = y(j1)
   ye(2) = ye(1)
   ye(3) = y(j2)
   ye(4) = ye(3)

   k1 = 1
   k2 = nz

   if (debug == 1) then 
        print *, 'before zp'
        print *, 'xp=', xp
        print *, 'yp=', yp
        print *, 'zp=', zp
   endif

   nsearch = 0
   do while ( (k2-k1)>1.1 )  
      if (nsearch >maxsearch) then 
         write (*,*) 'max number of searching has reached in z'
         write (*,*) 'zp = ', zp
         stop
      endif
      if (zp>=z(k1) .and. zp<z((k1+k2)/2)) then
         k2 = (k1+k2)/2
      else if (zp>=z((k1+k2)/2) .and. zp<z(k2)) then
         k1 = (k1+k2)/2
      endif

      nsearch = nsearch+1

   enddo

   if (debug ==1 .and. abs(xp-xdebug)<1e-4 .and. abs(yp-ydebug)<1e-4 .and. &
       abs(zp-zdebug)<1e-4 ) then
      print *, 'i1,i2,xe(1),xe(3):',j1,j2,ye(1),ye(3)
!      print *, 'k1,k2,ze(1),ze(3):',k1,k2,ze(1),ze(3)
   endif
        
   END SUBROUTINE LOCATION



   SUBROUTINE FILTER(iB,Nfilter)
!  Smooth out the pressure distribution on the wing surface

   USE flow_arrays
   USE flow_parameters
   USE pressure_arrays
   USE grid_arrays
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE

   REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: marker_in_element

   INTEGER iB,Nfilter
   INTEGER iM,IE,I,ifilter
   INTEGER MAXmarker_in_ele, MINmarker_in_ele

   ALLOCATE (marker_in_element(nBody,nPtsMax,2))


!  marker_in_element(:,:,1) saves the number of element associated with 
!                           that marker.
!  marker_in_element(:,:,2) saves the total area of elements associated with 
!                           that marker.
   marker_in_element(:,:,1:2) = 0

   MAXmarker_in_ele = 0
   Minmarker_in_ele = 100 

   DO IE=1,totNumTriElem(iB) 
      DO iM=1,3
         i = triElemNeig(iB,iM,IE)
         marker_in_element(iB,i,1) = marker_in_element(iB,i,1)+1
         marker_in_element(iB,i,2) = marker_in_element(iB,i,2)+triElemArea(iB,IE)
      ENDDO
!     MAXmarker_in_ele and MINmarker_in_ele just saves the Max/Min number of 
!     elements associated with certain nodes. Not used though.
      if (marker_in_element(iB,i,1)>MAXmarker_in_ele) MAXmarker_in_ele=marker_in_element(iB,i,1)
      if (marker_in_element(iB,i,1)<MINmarker_in_ele) MINmarker_in_ele=marker_in_element(iB,i,1)
   ENDDO

!   print *, 'MAXmarker_in_ele is:', MAXmarker_in_ele
!   print *, 'MINmarker_in_ele is:', MINmarker_in_ele

   DO ifilter=1,Nfilter

!      IF (ifilter == 1) THEN
         DO IE=1,totNumTriElem(iB) 
            DO iM=1,3
               i = triElemNeig(iB,iM,IE)
               pBodyMarker(iB,i) = 0.0
            ENDDO
         ENDDO
!      ENDIF

      DO IE=1,totNumTriElem(iB)
         DO iM=1,3
            i = triElemNeig(iB,iM,IE)
            pBodyMarker(iB,i) = pBodyMarker(iB,i)+ &
                                triElemP(iB,IE)*triElemArea(iB,IE)/3.0_CGREAL
         ENDDO
!        pBodyMarker(iB,i) is force at this step.
      ENDDO  !end loop of IE=1,totNumTriElem(iB)

      DO i = 1,nPtsBodyMarker(iB)
         pBodyMarker(iB,i) = pBodyMarker(iB,i)/marker_in_element(iB,i,2)
!        pBodyMarker(iB,i) is pressure at this step.
      ENDDO

      DO IE=1,totNumTriElem(iB)
!        Note there is no need to divide 3
         triElemP(iB,IE) = ( pBodyMarker(iB,triElemNeig(iB,1,IE)) + &
                             pBodyMarker(iB,triElemNeig(iB,2,IE)) + &
                             pBodyMarker(iB,triElemNeig(iB,3,IE)) )
      ENDDO

   ENDDO ! end loop of ifilter

   END SUBROUTINE FILTER


   SUBROUTINE FILTER_pxyz(iB,Nfilter)
!  Smooth out the pressure distribution on the wing surface

   USE flow_arrays
   USE flow_parameters
   USE pressure_arrays
   USE grid_arrays
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE

   REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: marker_in_element
   REAL(KIND=CGREAL) :: elemP

   INTEGER iB,Nfilter
   INTEGER iM,IE,I,ifilter
   INTEGER MAXmarker_in_ele, MINmarker_in_ele

   ALLOCATE (marker_in_element(nBody,nPtsMax,2))


!  marker_in_element(:,:,1) saves the number of element associated with 
!                           that marker.
!  marker_in_element(:,:,2) saves the total area of elements associated with 
!                           that marker.
   marker_in_element(:,:,1:2) = 0

   MAXmarker_in_ele = 0
   Minmarker_in_ele = 100 

   DO IE=1,totNumTriElem(iB) 
      DO iM=1,3
         i = triElemNeig(iB,iM,IE)
         marker_in_element(iB,i,1) = marker_in_element(iB,i,1)+1
         marker_in_element(iB,i,2) = marker_in_element(iB,i,2)+triElemArea(iB,IE)
      ENDDO
!     MAXmarker_in_ele and MINmarker_in_ele just saves the Max/Min number of 
!     elements associated with certain nodes. Not used though.
      if (marker_in_element(iB,i,1)>MAXmarker_in_ele) MAXmarker_in_ele=marker_in_element(iB,i,1)
      if (marker_in_element(iB,i,1)<MINmarker_in_ele) MINmarker_in_ele=marker_in_element(iB,i,1)
   ENDDO

!   print *, 'MAXmarker_in_ele is:', MAXmarker_in_ele
!   print *, 'MINmarker_in_ele is:', MINmarker_in_ele

   DO ifilter=1,Nfilter

!      IF (ifilter == 1) THEN
         DO IE=1,totNumTriElem(iB) 
            DO iM=1,3
               i = triElemNeig(iB,iM,IE)
               pBodyMarker(iB,i) = 0.0
               pxBodyMarker(iB,i) = 0.0
               pyBodyMarker(iB,i) = 0.0
               pzBodyMarker(iB,i) = 0.0
            ENDDO
         ENDDO
!      ENDIF

      DO IE=1,totNumTriElem(iB)
         DO iM=1,3
            i = triElemNeig(iB,iM,IE)
            elemP = triElemP(iB,IE)*triElemArea(iB,IE)/3.0_CGREAL
            pBodyMarker(iB,i) = pBodyMarker(iB,i)+elemP
                                
            pxBodyMarker(iB,i) = pxBodyMarker(iB,i)+   &
                                 triElemNormx(iB,IE)*elemP
            pyBodyMarker(iB,i) = pyBodyMarker(iB,i)+   &
                                 triElemNormy(iB,IE)*elemP 
            pzBodyMarker(iB,i) = pzBodyMarker(iB,i)+   &
                                 triElemNormz(iB,IE)*elemP
         ENDDO
!        pBodyMarker(iB,i) is force at this step.
      ENDDO  !end loop of IE=1,totNumTriElem(iB)

      DO i = 1,nPtsBodyMarker(iB)
         pBodyMarker(iB,i) = pBodyMarker(iB,i)/marker_in_element(iB,i,2)
!        pBodyMarker(iB,i) is pressure at this step.
         pxBodyMarker(iB,i) = pxBodyMarker(iB,i)/marker_in_element(iB,i,2)
         pyBodyMarker(iB,i) = pyBodyMarker(iB,i)/marker_in_element(iB,i,2)
         pzBodyMarker(iB,i) = pzBodyMarker(iB,i)/marker_in_element(iB,i,2)
      ENDDO

      DO IE=1,totNumTriElem(iB)
!        Note there is no need to divide 3
         triElemP(iB,IE) = ( pBodyMarker(iB,triElemNeig(iB,1,IE)) + &
                             pBodyMarker(iB,triElemNeig(iB,2,IE)) + &
                             pBodyMarker(iB,triElemNeig(iB,3,IE)) )
      ENDDO

   ENDDO ! end loop of ifilter

   END SUBROUTINE FILTER_pxyz


