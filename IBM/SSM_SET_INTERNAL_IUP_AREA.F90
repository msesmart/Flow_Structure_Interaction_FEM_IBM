!------------------------------------------------------------------------------

!!------------------------------------------------------------
!  SUBROUTINE SSM_set_internal_iup_membrane()
!!
!! Give a set of marker points, this subroutine computes
!!   1) normal intercepts and associated geometrical info. from ghost nodes to body
!!   2) Image point location
!!   3) Weights in stencil for computing values at the image points
!!
! 
!    USE global_parameters
!    USE flow_parameters
!    USE flow_arrays
!    USE grid_arrays
!    USE boundary_arrays
!    USE gcm_arrays
!    USE unstructured_surface_arrays
! 
!    IMPLICIT NONE
! 
!!... loop variables
! 
!    INTEGER :: i,iBody,iRow,j,k,m,n
! 
!!... local variables
!    INTEGER(1) :: boundCell(0:nx+1,0:ny+1,0:nz+1)
!
!!    INTEGER :: iG, jG, kG, nbdr, iCIndx, jCIndx, kCIndx , iCIndxS, jCIndxS
! 
!    INTEGER :: iG, jG, kG, nbdr, nbdrG, iCIndx, jCIndx, kCIndx , iCIndxS, jCIndxS  ! new statement
!    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
!    INTEGER :: iM, iP, jM, jP, kM, kP
!!    INTEGER :: iRange, jRange, kRange
!! new statements
!    INTEGER :: iRange, jRange, kRange, iC, jC, kC
!    INTEGER :: checkGhostFlag
! 
!    INTEGER ::  nGhost_mem_1side,iitemp,nVertex
! 
!    REAL(KIND=CGREAL) :: cosTheta,dsIntercept,sinTheta,               &
!                         xBI,xBIN,xGC,xIP,xIPS,yBI,yBIN,yGC,yIP,yIPS,zBI,zBIN,zGC,zIP, &
!                         minProbeLengthShear,slopeX, slopeY, slopeZ, maxDelta
! 
!!*****************************************************************************************
! 
!! Initialize arrays 
!    
!      write(*,*) 'This process is for membrane bodies !!'
! 
!      iblank = iblank_memb   ! for detecting the ghostcells
!
!      ghostCellMark  = 0
! 
!! Find all cells that contain a IB node
! 
!      boundCell = 0
! 
!    DO iBody = nBody_solid+1, nBody
!
!
!      ic = -1
!      jc = -1
!      kc = -1
!
!     DO m = 1,nPtsBodyMarker(iBody)
! 
!       DO i = 1,nx
!         IF ( ( xc(i) <= xBodyMarker(iBody,m)  .AND. xc(i+1) > xBodyMarker(iBody,m) ) ) THEN
!           iC = i
!         ENDIF ! xc
!       ENDDO ! i
!
!       DO j = 1,ny
!         IF ( ( yc(j) <= yBodyMarker(iBody,m)  .AND. yc(j+1) > yBodyMarker(iBody,m) ) ) THEN
!           jC = j
!         ENDIF ! yc   
!       ENDDO ! j
!   
!       DO k = 1,nz
!         IF ( ( zc(k) <= zBodyMarker(iBody,m)  .AND. zc(k+1) > zBodyMarker(iBody,m) ) ) THEN
!           kC = k
!         ENDIF ! zc 
!       ENDDO ! k 
!
! 
!       boundCell(iC,jC,kC)       = 1
!       boundCell(iC+1,jC,kC)     = 1 
!       boundCell(iC,jC+1,kC)     = 1
!       boundCell(iC,jC,kC+1)     = 1
!       boundCell(iC+1,jC+1,kC)   = 1
!       boundCell(iC,jC+1,kC+1)   = 1
!       boundCell(iC+1,jC,kC+1)   = 1
!       boundCell(iC+1,jC+1,kC+1) = 1 
!            
!!write(897,*)m,ic,jc,kc
!
!     ENDDO
!
!   ENDDO ! end iBody
!
!! Above parts are for detecting boundcells.
!!
!! Mark boundary cells (ie. cells inside body which have
!! at least one neighbor in fluid)
!
!!      CALL GCM_AllocateGhostCellArrays()  !Added on 05/21/10 to allocate iGhost
!
!      nbdr = 0
!
!      DO k = 1, nzc
!      DO j = 1, nyc
!      DO i = 1, nxc
! 
!        iM = MAX(i-1,1)
!        iP = MIN(i+1,nxc)
!        jM = MAX(j-1,1)
!        jP = MIN(j+1,nyc)
!        kM = MAX(k-1,1)
!        kP = MIN(k+1,nzc)
!      
!!  new replacement
! 
!        IF (boundCell(i,j,k) == 1) THEN 
! 
!          IF ( ( iblank(i,j,k)+iblank(iM,j,k) == 1 .AND. boundCell(iM,j,k) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(iP,j,k) == 1 .AND. boundCell(iP,j,k) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(i,jM,k) == 1 .AND. boundCell(i,jM,k) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(i,jP,k) == 1 .AND. boundCell(i,jP,k) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(i,j,kM) == 1 .AND. boundCell(i,j,kM) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(i,j,kP) == 1 .AND. boundCell(i,j,kP) == 1 )      )  THEN
!              
!              ghostCellMark(i,j,k) = 1
!              nbdr = nbdr + 1
!              iGhostP(nbdr) = i
!              jGhostP(nbdr) = j
!              kGhostP(nbdr) = k
!          ENDIF
! 
!        ENDIF
!
!      ENDDO ! i
!      ENDDO ! j
!      ENDDO ! k
!
!    nGhost = nbdr
!
!    CALL set_internal_iup_GCM()
!       
!    ghostCellMemb =  ghostCellMark
!
!    iblank = 0
!
!    ghostCellMark = 0
!
!  END SUBROUTINE SSM_set_internal_iup_membrane


  
!  SUBROUTINE set_internal_iup()
!
!    USE global_parameters
!    USE flow_parameters
!    USE grid_arrays
!    USE boundary_arrays
!
!    IMPLICIT NONE
!
!    INTEGER           :: i,j,k
!
!    iblank = iblank_solid
!
!    DO k=1,nzc
!    DO j=1,nyc
!    DO i=1,nxc
!      IF (iblank(i,j,k) == 0 ) THEN
!        IF (iblank(i+1,j,k) == 1) iup(i,j,k)=oned 
!        IF (iblank(i-1,j,k) == 1) ium(i,j,k)=oned 
!        IF (iblank(i,j+1,k) == 1) jup(i,j,k)=oned 
!        IF (iblank(i,j-1,k) == 1) jum(i,j,k)=oned 
!        IF (iblank(i,j,k+1) == 1) kup(i,j,k)=oned 
!        IF (iblank(i,j,k-1) == 1) kum(i,j,k)=oned
!
!
!        IF (iblank(i+1,j,k) == 1) iupp(i,j,k)=oned    !
!        IF (iblank(i-1,j,k) == 1) iumm(i,j,k)=oned    ! Added by Rupesh
!        IF (iblank(i,j+1,k) == 1) jupp(i,j,k)=oned    ! Flags the gridline
!        IF (iblank(i,j-1,k) == 1) jumm(i,j,k)=oned    ! making the boundary
!        IF (iblank(i,j,k+1) == 1) kupp(i,j,k)=oned    ! just as iup, etc do.
!        IF (iblank(i,j,k-1) == 1) kumm(i,j,k)=oned    !
! 
!        IF (iblank(i+1,j,k) == 1) iupp(i-1,j,k)=oned  !
!        IF (iblank(i-1,j,k) == 1) iumm(i+1,j,k)=oned  ! Added by Rupesh
!        IF (iblank(i,j+1,k) == 1) jupp(i,j-1,k)=oned  ! Flags the gridline next
!        IF (iblank(i,j-1,k) == 1) jumm(i,j+1,k)=oned  ! (in the  outward normal dir)
!        IF (iblank(i,j,k+1) == 1) kupp(i,j,k-1)=oned  ! to the already flagged gridline
!        IF (iblank(i,j,k-1) == 1) kumm(i,j,k+1)=oned  ! making the boundary.
!
!      ENDIF
!    ENDDO
!    ENDDO
!    ENDDO
!
!   END SUBROUTINE set_internal_iup
!------------------------------------------------------------------------------

!------------------------------------------------------------
!  SUBROUTINE SSM_set_internal_iup_solid()
!!
!! Give a set of marker points, this subroutine computes
!!   1) normal intercepts and associated geometrical info. from ghost nodes to body
!!   2) Image point location
!!   3) Weights in stencil for computing values at the image points
!!
! 
!    USE global_parameters
!    USE flow_parameters
!    USE flow_arrays
!    USE grid_arrays
!    USE boundary_arrays
!    USE gcm_arrays
!    USE unstructured_surface_arrays
! 
!    IMPLICIT NONE
! 
!!... loop variables
! 
!    INTEGER :: i,iBody,iRow,j,k,m,n
! 
!!... local variables
!    INTEGER(1) :: boundCell(0:nx+1,0:ny+1,0:nz+1)
!    
!!    INTEGER :: iG, jG, kG, nbdr, iCIndx, jCIndx, kCIndx , iCIndxS, jCIndxS
! 
!    INTEGER :: iG, jG, kG, nbdr, nbdrG, iCIndx, jCIndx, kCIndx , iCIndxS, jCIndxS  ! new statement
!    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
!    INTEGER :: iM, iP, jM, jP, kM, kP
!!    INTEGER :: iRange, jRange, kRange
!! new statements
!    INTEGER :: iRange, jRange, kRange, iC, jC, kC
!    INTEGER :: checkGhostFlag
! 
!    INTEGER ::  nGhost_mem_1side,iitemp,nVertex
! 
!    REAL(KIND=CGREAL) :: cosTheta,dsIntercept,sinTheta,               &
!                         xBI,xBIN,xGC,xIP,xIPS,yBI,yBIN,yGC,yIP,yIPS,zBI,zBIN,zGC,zIP, &
!                         minProbeLengthShear,slopeX, slopeY, slopeZ, maxDelta
! 
!!*****************************************************************************************
! 
!! Initialize arrays 
!    
!      write(*,*) 'This process is for solid bodies !!'
! 
!      iblank = iblank_solid   ! for detecting the ghostcells
!
!      ghostCellMark  = 0
! 
!! Find all cells that contain a IB node
! 
!      boundCell = 0
! 
!    DO iBody = 1, nBody_solid
!
!      ic = 1
!      jc = 1
!      kc = 1
!
!     DO m = 1,nPtsBodyMarker(iBody)
! 
!       DO i = 1,nx
!         IF ( ( xc(i) <= xBodyMarker(iBody,m)  .AND. xc(i+1) > xBodyMarker(iBody,m) ) ) THEN
!           iC = i
!         ENDIF ! xc
!       ENDDO ! i
!
!       DO j = 1,ny
!         IF ( ( yc(j) <= yBodyMarker(iBody,m)  .AND. yc(j+1) > yBodyMarker(iBody,m) ) ) THEN
!           jC = j
!         ENDIF ! yc   
!       ENDDO ! j
!   
!       DO k = 1,nz
!         IF ( ( zc(k) <= zBodyMarker(iBody,m)  .AND. zc(k+1) > zBodyMarker(iBody,m) ) ) THEN
!           kC = k
!         ENDIF ! zc 
!       ENDDO ! k 
! 
!       boundCell(iC,jC,kC)       = 1
!       boundCell(iC+1,jC,kC)     = 1 
!       boundCell(iC,jC+1,kC)     = 1
!       boundCell(iC,jC,kC+1)     = 1
!       boundCell(iC+1,jC+1,kC)   = 1
!       boundCell(iC,jC+1,kC+1)   = 1
!       boundCell(iC+1,jC,kC+1)   = 1
!       boundCell(iC+1,jC+1,kC+1) = 1 
!            
!     ENDDO
!
!   ENDDO ! end iBody
!
!! Above parts are for detecting boundcells.
!!
!! Mark boundary cells (ie. cells inside body which have
!! at least one neighbor in fluid)
! 
!      DO k = 1, nzc
!      DO j = 1, nyc
!      DO i = 1, nxc
! 
!        iM = MAX(i-1,1)
!        iP = MIN(i+1,nxc)
!        jM = MAX(j-1,1)
!        jP = MIN(j+1,nyc)
!        kM = MAX(k-1,1)
!        kP = MIN(k+1,nzc)
!      
!!  new replacement
! 
!        IF (boundCell(i,j,k) == 1) THEN 
! 
!          IF ( ( iblank(i,j,k)+iblank(iM,j,k) == 1 .AND. boundCell(iM,j,k) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(iP,j,k) == 1 .AND. boundCell(iP,j,k) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(i,jM,k) == 1 .AND. boundCell(i,jM,k) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(i,jP,k) == 1 .AND. boundCell(i,jP,k) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(i,j,kM) == 1 .AND. boundCell(i,j,kM) == 1 )  .OR.  &
!               ( iblank(i,j,k)+iblank(i,j,kP) == 1 .AND. boundCell(i,j,kP) == 1 )      )  THEN
!           IF (iblank(i,j,k) == 1) THEN
!             ghostCellMark(i,j,k) = 1
!             nbdr = nbdr + 1
!           ENDIF 
!             IF (iblank(i,j,k) == 0) THEN
!              ghostCellMark(i,j,k) = 1
!              nbdr = nbdr + 1
!             ENDIF
!          ENDIF
! 
!        ENDIF
! 
!      ENDDO ! i
!      ENDDO ! j
!      ENDDO ! k
!
!!end here 
!       
!    ghostCellSolid =  ghostCellMark
!
!    DO k=1,nzc
!    DO j=1,nyc
!    DO i=1,nxc
!      IF (iblank(i,j,k) == 0 ) THEN
!        IF (iblank(i+1,j,k) == 1) iup(i,j,k)=oned
!        IF (iblank(i-1,j,k) == 1) ium(i,j,k)=oned
!        IF (iblank(i,j+1,k) == 1) jup(i,j,k)=oned
!        IF (iblank(i,j-1,k) == 1) jum(i,j,k)=oned
!        IF (iblank(i,j,k+1) == 1) kup(i,j,k)=oned
!        IF (iblank(i,j,k-1) == 1) kum(i,j,k)=oned
! 
! 
!        IF (iblank(i+1,j,k) == 1) iupp(i,j,k)=oned    !
!        IF (iblank(i-1,j,k) == 1) iumm(i,j,k)=oned    ! Added by Rupesh
!        IF (iblank(i,j+1,k) == 1) jupp(i,j,k)=oned    ! Flags the gridline
!        IF (iblank(i,j-1,k) == 1) jumm(i,j,k)=oned    ! making the boundary
!        IF (iblank(i,j,k+1) == 1) kupp(i,j,k)=oned    ! just as iup, etc do.
!        IF (iblank(i,j,k-1) == 1) kumm(i,j,k)=oned    !
! 
!        IF (iblank(i+1,j,k) == 1) iupp(i-1,j,k)=oned  !
!        IF (iblank(i-1,j,k) == 1) iumm(i+1,j,k)=oned  ! Added by Rupesh
!        IF (iblank(i,j+1,k) == 1) jupp(i,j-1,k)=oned  ! Flags the gridline next
!        IF (iblank(i,j-1,k) == 1) jumm(i,j+1,k)=oned  ! (in the  outward normal dir)
!        IF (iblank(i,j,k+1) == 1) kupp(i,j,k-1)=oned  ! to the already flagged gridline
!        IF (iblank(i,j,k-1) == 1) kumm(i,j,k+1)=oned  ! making the boundary.
! 
!      ENDIF
!    ENDDO
!    ENDDO
!    ENDDO
!
!  END SUBROUTINE SSM_set_internal_iup_solid

   SUBROUTINE set_internal_iup_GCM()
    
    USE global_parameters
    USE flow_parameters 
    USE grid_arrays 
    USE boundary_arrays
     
    IMPLICIT NONE
     
    INTEGER           :: i,j,k
     
    DO k=1,nzc
    DO j=1,nyc
    DO i=1,nxc
      IF (ghostCellMark(i,j,k) == 1 ) THEN
        IF (iblank(i+1,j,k) /= iblank(i,j,k) .AND. ghostCellMark(i+1,j,k) == 1) iup(i,j,k)=1
        IF (iblank(i-1,j,k) /= iblank(i,j,k) .AND. ghostCellMark(i-1,j,k) == 1) ium(i,j,k)=1
        IF (iblank(i,j+1,k) /= iblank(i,j,k) .AND. ghostCellMark(i,j+1,k) == 1) jup(i,j,k)=1
        IF (iblank(i,j-1,k) /= iblank(i,j,k) .AND. ghostCellMark(i,j-1,k) == 1) jum(i,j,k)=1
        IF (iblank(i,j,k+1) /= iblank(i,j,k) .AND. ghostCellMark(i,j,k+1) == 1) kup(i,j,k)=1
        IF (iblank(i,j,k-1) /= iblank(i,j,k) .AND. ghostCellMark(i,j,k-1) == 1) kum(i,j,k)=1
      ENDIF
    ENDDO
    ENDDO
    ENDDO
 
   END SUBROUTINE set_internal_iup_GCM



!------------------------------------------------------------------------------

   SUBROUTINE SSM_set_internal_area()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER           :: i,j,k

!...Compute surface area of immersed boundary

      areax1=zero
      areax2=zero
      areay1=zero
      areay2=zero
      areaz1=zero
      areaz2=zero

      DO k = 1, nzc
      DO j = 1, nyc
      DO i = 1, nxc

        areax1 = areax1 + iup(i,j,k)*REAL(iblank(i+1,j,k),KIND=CGREAL)*dy(j)*dz(k)
        areax2 = areax2 + ium(i,j,k)*REAL(iblank(i-1,j,k),KIND=CGREAL)*dy(j)*dz(k)

        areay1 = areay1 + jup(i,j,k)*REAL(iblank(i,j+1,k),KIND=CGREAL)*dx(i)*dz(k)
        areay2 = areay2 + jum(i,j,k)*REAL(iblank(i,j-1,k),KIND=CGREAL)*dx(i)*dz(k) 

        areaz1 = areaz1 + kup(i,j,k)*REAL(iblank(i,j,k+1),KIND=CGREAL)*dx(i)*dy(j)
        areaz2 = areaz2 + kum(i,j,k)*REAL(iblank(i,j,k-1),KIND=CGREAL)*dx(i)*dy(j) 

      END DO
      END DO
      END DO

   END SUBROUTINE SSM_set_internal_area
!------------------------------------------------------------------------------


SUBROUTINE SSM_set_internal_iup_solid()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER           :: i,j,k
    logical :: keepGo


!888 continue


    DO k=1,nzc
    DO j=1,nyc
    DO i=1,nxc
      IF (iblank(i,j,k) == 0) THEN
        IF (ghostCellSolid(i+1,j,k) == -1) iup(i,j,k)=1
        IF (ghostCellSolid(i-1,j,k) == -1) ium(i,j,k)=1
        IF (ghostCellSolid(i,j+1,k) == -1) jup(i,j,k)=1
        IF (ghostCellSolid(i,j-1,k) == -1) jum(i,j,k)=1
        IF (ghostCellSolid(i,j,k+1) == -1) kup(i,j,k)=1
        IF (ghostCellSolid(i,j,k-1) == -1) kum(i,j,k)=1
      END IF
    ENDDO
    ENDDO
    ENDDO

!    keepGo=.true.
!    do k=1,nzc
!    do j=1,nyc
!    do i=1,nxc
!        if(iblank(i,j,k)==0)then
!        if((iup(i,j,k)==1.and.ium(i,j,k)==1).or.&
!            (jup(i,j,k)==1.and.jum(i,j,k)==1).or.&
!            (kup(i,j,k)==1.and.kum(i,j,k)==1))then
!            iup(i,j,k)=0
!            ium(i+1,j,k)=0
!            ium(i,j,k)=0
!            iup(i-1,j,k)=0
!            jup(i,j,k)=0
!            jum(i,j+1,k)=0
!            jum(i,j,k)=0
!            jup(i,j-1,k)=0
!            kup(i,j,k)=0
!            kum(i,j,k+1)=0
!            kum(i,j,k)=0
!            kup(i,j,k-1)=0
!
!            iup(i,j-1,k)=0
!            ium(i+1,j-1,k)=0
!            ium(i,j-1,k)=0
!            iup(i-1,j-1,k)=0
!            jup(i,j-1,k)=0
!            jum(i,j,k)=0
!            jum(i,j-1,k)=0
!            jup(i,j-2,k)=0
!            kup(i,j-1,k)=0
!            kum(i,j-1,k+1)=0
!            kum(i,j-1,k)=0
!            kup(i,j-1,k-1)=0
!
!            if(iup(i,j,k)==1.and.ium(i,j,k)==1)then
!                iblank_memb(i-1:i+1,j,k)=0
!                ghostCellMemb(i-1:i+1,j,k)=0
!            end if
!            if(jup(i,j,k)==1.and.jum(i,j,k)==1)then
!                iblank_memb(i,j-1:j+1,k)=0
!                ghostCellMemb(i,j-1:j+1,k)=0
!            end if
!            if(kup(i,j,k)==1.and.kum(i,j,k)==1)then
!                iblank_memb(i,j,k-1:k+1)=0
!                ghostCellMemb(i,j,k-1:k+1)=0
!            end if
!
!            keepGo=.false.
!        end if
!        end if
!    end do
!    end do
!    end do
!    
!
!    if(.not.keepGo)go to 888

    if(ntime==0)then
    write(5679,*)  'VARIABLES="X","Y","iblank","ghostCellSolid","iup","ium","jup","jum"'
    write(5679,*)  'ZONE F=POINT, I=128 , J=160'
    do k=1,1
    do j=1,nyc
    do i=1,nxc
      write(5679,5679) xc(i),yc(j),iblank(i,j,k),ghostCellSolid(i,j,k),iup(i,j,k),ium(i,j,k),jup(i,j,k),jum(i,j,k)
    end do
    end do
    end do
    close(5679)
    end if


    IF (Hybrid) THEN
      DO k=1,nzc
      DO j=1,nyc
      DO i=1,nxc
        IF (iblank(i,j,k) == 0) THEN
          IF (ghostCellSolid(i+1,j,k) == -1) iupp(i,j,k)=1    !
          IF (ghostCellSolid(i-1,j,k) == -1) iumm(i,j,k)=1    ! Added by Rupesh
          IF (ghostCellSolid(i,j+1,k) == -1) jupp(i,j,k)=1    ! Flags the gridline making 
          IF (ghostCellSolid(i,j-1,k) == -1) jumm(i,j,k)=1    ! the boundary just as iup, etc do.
          IF (ghostCellSolid(i,j,k+1) == -1) kupp(i,j,k)=1    ! 
          IF (ghostCellSolid(i,j,k-1) == -1) kumm(i,j,k)=1    !

          IF (ghostCellSolid(i+1,j,k) == -1) iupp(i-1,j,k)=1  ! 
          IF (ghostCellSolid(i-1,j,k) == -1) iumm(i+1,j,k)=1  ! Added by Rupesh
          IF (ghostCellSolid(i,j+1,k) == -1) jupp(i,j-1,k)=1  ! Flags the gridline next (in the 
          IF (ghostCellSolid(i,j-1,k) == -1) jumm(i,j+1,k)=1  ! outward normal dir) to the already 
          IF (ghostCellSolid(i,j,k+1) == -1) kupp(i,j,k-1)=1  ! flagged gridline making the boundary. 
          IF (ghostCellSolid(i,j,k-1) == -1) kumm(i,j,k+1)=1  !
        END IF
      ENDDO
      ENDDO
      ENDDO

    ENDIF

    5679 format(2f16.8,6i4)
END SUBROUTINE SSM_set_internal_iup_solid
!---------------------------------------------------------------------



SUBROUTINE SSM_set_internal_iup_membrane()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER           :: i,j,k


    DO k=1,nzc
    DO j=1,nyc
    DO i=1,nxc
      IF (iblank_memb(i,j,k) == 1) THEN
        IF ( iblank_memb(i+1,j,k) == 1 &
             .AND. ghostCellMark(i,j,k) * ghostCellMark(i+1,j,k) == -1) iup(i,j,k)=1
        IF ( iblank_memb(i,j+1,k) == 1  &
             .AND. ghostCellMark(i,j,k) * ghostCellMark(i,j+1,k) == -1) jup(i,j,k)=1
        IF ( iblank_memb(i,j,k+1) == 1 &
             .AND. ghostCellMark(i,j,k) * ghostCellMark(i,j,k+1) == -1) kup(i,j,k)=1
        IF ( iblank_memb(i-1,j,k) == 1 &
             .AND. ghostCellMark(i,j,k) * ghostCellMark(i-1,j,k) == -1) ium(i,j,k)=1
        IF ( iblank_memb(i,j-1,k) == 1 &
             .AND. ghostCellMark(i,j,k) * ghostCellMark(i,j-1,k) == -1) jum(i,j,k)=1
        IF ( iblank_memb(i,j,k-1) == 1 &
             .AND. ghostCellMark(i,j,k) * ghostCellMark(i,j,k-1) == -1) kum(i,j,k)=1

!
!!!!!!!!!!!!THIS IS NOT SET UP FOR MEMBRANE
!        IF (ghostCellMark(i+1,j,k) == 1 .AND. iblank_memb(i,j,k) == 0) iupp(i,j,k)=1    !
!        IF (ghostCellMark(i-1,j,k) == 1 .AND. iblank_memb(i,j,k) == 0) iumm(i,j,k)=1    ! Added by Rupesh
!        IF (ghostCellMark(i,j+1,k) == 1 .AND. iblank_memb(i,j,k) == 0) jupp(i,j,k)=1    ! Flags the gridline making 
!        IF (ghostCellMark(i,j-1,k) == 1 .AND. iblank_memb(i,j,k) == 0) jumm(i,j,k)=1    ! the boundary just as iup, etc do.
!        IF (ghostCellMark(i,j,k+1) == 1 .AND. iblank_memb(i,j,k) == 0) kupp(i,j,k)=1    ! 
!        IF (ghostCellMark(i,j,k-1) == 1 .AND. iblank_memb(i,j,k) == 0) kumm(i,j,k)=1    !
!
!        IF (ghostCellMark(i+1,j,k) == 1 .AND. iblank_memb(i,j,k) == 0) iupp(i-1,j,k)=1  ! 
!        IF (ghostCellMark(i-1,j,k) == 1 .AND. iblank_memb(i,j,k) == 0) iumm(i+1,j,k)=1  ! Added by Rupesh
!        IF (ghostCellMark(i,j+1,k) == 1 .AND. iblank_memb(i,j,k) == 0) jupp(i,j-1,k)=1  ! Flags the gridline next (in the 
!        IF (ghostCellMark(i,j-1,k) == 1 .AND. iblank_memb(i,j,k) == 0) jumm(i,j+1,k)=1  ! outward normal dir) to the already 
!        IF (ghostCellMark(i,j,k+1) == 1 .AND. iblank_memb(i,j,k) == 0) kupp(i,j,k-1)=1  ! flagged gridline making the boundary. 
!        IF (ghostCellMark(i,j,k-1) == 1 .AND. iblank_memb(i,j,k) == 0) kumm(i,j,k+1)=1  !
	  END IF
    ENDDO
    ENDDO
    ENDDO

    if(ntime==0)then
    write(5678,*)  'VARIABLES="X","Y","iblank_memb","ghostCellMark","iup","ium","jup","jum"'
    write(5678,*)  'ZONE F=POINT, I=128 , J=160'
    do k=1,1
    do j=1,nyc
    do i=1,nxc
      write(5678,5678) xc(i),yc(j),iblank_memb(i,j,k),ghostCellMark(i,j,k),iup(i,j,k),ium(i,j,k),jup(i,j,k),jum(i,j,k)
    end do
    end do
    end do
    close(5678)
    end if


    5678 format(2f16.8,6i4)

END SUBROUTINE SSM_set_internal_iup_membrane
!---------------------------------------------------------------------
