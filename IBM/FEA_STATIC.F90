   SUBROUTINE FEA_STATIC(iBody)
!  STATIC solution
   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface
   USE flow_parameters, ONLY : nPtsBodyMarker

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: STFTMP(MAXSTIFF)
   REAL(KIND=CGREAL) :: LOADEFF(NEQ)

   INTEGER :: iBody
   INTEGER :: I,ierror


!  Decompose effective stiffness matrix

   DO I=1,maxstiff
      STFTMP(I)=STF(I)  
   ENDDO

   ierror=0
   CALL  uduCOL(STFTMP,ierror)

   IF (ierror .eq. 0) THEN
      write(*,*)'ERROR: zero diagonal term'
      RETURN
   ENDIF

   DO I = 1,NEQ
      LOADEFF(I) = LOAD(I)
   ENDDO

   WRITE (*,*) 'max/min(load) in FEA_STATIC is:', MAXVAL(LOAD),MINVAL(LOAD)   
                
!  SOLVE for new displacements 
   CALL bakCOL(STFTMP,maxstiff,LOADEFF,ierror)


! 22 format(1x,7(g12.5,1x))

   CALL STATOUT(iBody)

   END SUBROUTINE FEA_STATIC


   SUBROUTINE STATOUT(iBody)
!  STATic OUTput routine for displacements and nodal loads

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface
   USE flow_parameters, ONLY : nPtsBodyMarker

   IMPLICIT NONE
   INTEGER :: iBody

!   REAL(KIND=CGREAL) :: dispful(nPtsBodyMarker(iBody),nodeDoF)

   INTEGER :: i,j,node,jdof,ieqnum,idis

   idis = 179

   OPEN (idis,FILE='dis.out')
   WRITE (idis, 1100)
 1100  FORMAT(1x,a50,/, &
      '  node    x-disp       y-disp       z-disp ',  &
      '      x-rot        y-rot        z-rot' )

!  Fill in  the full displacement vector
   DO i=1, nPtsBodyMarker(iBody)*nodeDoF
      node=(i+5)/6       ! node is marker index
      jdof=i-(node-1)*6  ! jdof is always from 1 to 6.
      ieqnum = jbc(i)    ! index of equation number.

      if (ieqnum .gt. 0) then
         struc_disp_worelax(node,jdof) = wk(ieqnum)
      else
         struc_disp_worelax(node,jdof) = 0.0e0
      endif
   ENDDO

   DO i=1, nPtsBodyMarker(iBody)
      write(idis,1110) i, (struc_disp_worelax(i,j),j=1,6) 
   ENDDO

 1110  FORMAT(1x,i5,2x,6(g12.6,1x))


   END SUBROUTINE STATOUT





