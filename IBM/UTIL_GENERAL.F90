!------------------------------------------    
!  SUBROUTINE write_monitor() 
!  SUBROUTINE write_dump()
!  SUBROUTINE write_minmaxvals
!  SUBROUTINE vorticity()
!  SUBROUTINE divergence()
!------------------------------------------    



!------------------------------------------    
!------------------------------------------    
    SUBROUTINE write_monitor() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE nlold_arrays
    USE pressure_arrays

    IMPLICIT NONE
    INTEGER                 :: i,j,k
    INTEGER, DIMENSION(3)   :: maxcflloc,minddfloc
    REAL(KIND=CGREAL)       :: maxcfl,minddf,rnDim

! cfl = dt*{ |u|/dx  + |v|/dy + |w|/dz }
! diagonal dominance factor of advection diffusion equation
!  = (1 + rx + ry + rz)/(rx + ry + rz)

    nlv = zero
    nlw = zero

    rnDim = REAL((nDim-DIM_2D),KIND=CGREAL)

    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      nlv(i,j,k) = dt*( abs(u(i,j,k))*dxinv(i)   &
                       +abs(v(i,j,k))*dyinv(j)   &
                       +abs(w(i,j,k))*dzinv(k) ) &
                 *REAL(1-iblank(i,j,k),KIND=CGREAL)

      nlw(i,j,k) =  dt*reinv* ( dxinv(i)**2       &
                               +dyinv(j)**2       &
                               +rnDim*dzinv(k)**2 )
      nlw(i,j,k) = oned + oned/nlw(i,j,k)
    ENDDO
    ENDDO
    ENDDO

    maxcfl    = MAXVAL(nlv(1:nxc,1:nyc,1:nzc))
    maxcflloc = MAXLOC(nlv(1:nxc,1:nyc,1:nzc))
    WRITE(STDOUT,'(A,E15.7,A,3(2X,I4))')'Max CFL              = ',maxcfl,' at ',maxcflloc

    minddf    = MINVAL(nlw(1:nxc,1:nyc,1:nzc))
    minddfloc = MINLOC(nlw(1:nxc,1:nyc,1:nzc))
    WRITE(STDOUT,'(A,E15.7,A,3(2X,I4))')'Min Diagonal Dom. Fac AD = ',minddf,' at ',minddfloc
!   IF ( minddf .LT. 20.0 ) THEN
!     print*,'Warning----advection diffusion operator has weak diagonal dominance'
!     print*,'If you observe problems with convergence of advection-diffusion equation then might need to reduce dt'
!   ENDIF

    CALL divergence()

    CALL write_minmaxvals()

    END SUBROUTINE write_monitor 
!------------------------------------------
!    
!------------------------------------------    
    SUBROUTINE write_minmaxvals() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE nlold_arrays
    USE pressure_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE turb_parameters

    IMPLICIT NONE
    INTEGER                 :: i,j,k

    PRINT*, 'Min-Max of U = ',minval(u(1:nxc,1:nyc,1:nzc)),maxval(u(1:nxc,1:nyc,1:nzc))
    PRINT*, 'Min-Max of V = ',minval(v(1:nxc,1:nyc,1:nzc)),maxval(v(1:nxc,1:nyc,1:nzc))
    PRINT*, 'Min-Max of W = ',minval(w(1:nxc,1:nyc,1:nzc)),maxval(w(1:nxc,1:nyc,1:nzc))
    PRINT*, 'Min-Max of P = ',minval(p(1:nxc,1:nyc,1:nzc)),maxval(p(1:nxc,1:nyc,1:nzc))

    PRINT*, 'Min-Max of FACE_U = ',minval(face_u(1:nxc,1:nyc,1:nzc)),maxval(face_u(1:nxc,1:nyc,1:nzc))
    PRINT*, 'Min-Max of FACE_V = ',minval(face_v(1:nxc,1:nyc,1:nzc)),maxval(face_v(1:nxc,1:nyc,1:nzc))
    PRINT*, 'Min-Max of FACE_W = ',minval(face_w(1:nxc,1:nyc,1:nzc)),maxval(face_w(1:nxc,1:nyc,1:nzc))

    PRINT*, 'Min-Max of NLU = ',minval(nlu(1:nxc,1:nyc,1:nzc)),maxval(nlu(1:nxc,1:nyc,1:nzc))
    PRINT*, 'Min-Max of NLV = ',minval(nlv(1:nxc,1:nyc,1:nzc)),maxval(nlv(1:nxc,1:nyc,1:nzc))
    PRINT*, 'Min-Max of NLW = ',minval(nlw(1:nxc,1:nyc,1:nzc)),maxval(nlw(1:nxc,1:nyc,1:nzc))
    
    PRINT*, 'Min-Max of U(0: = ',minval(u(0:nx+1,0:ny+1,0:nz+1)),maxval(u(0:nx+1,0:ny+1,0:nz+1))
    PRINT*, 'Min-Max of V(0: = ',minval(v(0:nx+1,0:ny+1,0:nz+1)),maxval(v(0:nx+1,0:ny+1,0:nz+1))
    PRINT*, 'Min-Max of W(0: = ',minval(w(0:nx+1,0:ny+1,0:nz+1)),maxval(w(0:nx+1,0:ny+1,0:nz+1))
    
    PRINT*, 'Min-Max of NLU(0: = ',minval(nlu(0:nx+1,0:ny+1,0:nz+1)),maxval(nlu(0:nx+1,0:ny+1,0:nz+1))
    PRINT*, 'Min-Max of NLV(0: = ',minval(nlv(0:nx+1,0:ny+1,0:nz+1)),maxval(nlv(0:nx+1,0:ny+1,0:nz+1))
    PRINT*, 'Min-Max of NLW(0: = ',minval(nlw(0:nx+1,0:ny+1,0:nz+1)),maxval(nlw(0:nx+1,0:ny+1,0:nz+1))

    IF ( turbActive /= ACTIVE ) THEN
      PRINT*, 'Min-Max of viscTot = ',&
      MINVAL(viscTot(1:nxc,1:nyc,1:nzc)),MAXVAL(viscTot(1:nxc,1:nyc,1:nz- 1))
    ENDIF ! turbModel
    
!    PRINT*,'Planar Min-Max Velocity Field Values'
!    DO k = 1, nzc
!      PRINT*, '  Min-Max of U-kplane(',k,')=',MINVAL(u(1:nxc,1:nyc,k)),&
!                                              MAXVAL(u(1:nxc,1:nyc,k))
!    ENDDO
!    DO k = 1, nzc
!      PRINT*, '  Min-Max of V-kplane(',k,')=',MINVAL(v(1:nxc,1:nyc,k)),&
!                                              MAXVAL(v(1:nxc,1:nyc,k))
!    ENDDO
!    DO k = 1, nzc
!      PRINT*, '  Min-Max of W-kplane(',k,')=',MINVAL(w(1:nxc,1:nyc,k)),&
!                                              MAXVAL(w(1:nxc,1:nyc,k))
!    ENDDO
!    DO k = 1, nzc
!      PRINT*, '  Min-Max of P-kplane(',k,')=',MINVAL(p(1:nxc,1:nyc,k)),&
!                                              MAXVAL(p(1:nxc,1:nyc,k))
!    ENDDO

    END SUBROUTINE write_minmaxvals 
!------------------------------------------
!
!------------------------------------------    
    SUBROUTINE write_dump() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE nlold_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays


    IMPLICIT NONE
   
    REAL            :: fsmach,alp,density,xhat,yhat
    INTEGER         :: i,j,k,iBody,n,m,nCylinder, ii, jj, kk
    CHARACTER*13    :: fname1
    CHARACTER*20    :: fname2
    CHARACTER*256	  :: tmp, tmp1
	
    PRINT*,'Writing out dump file'
    if(.not.dryRun)then
        IF ( ntime >= 0 .AND. ntime .le. 9  )            &
            write(fname1,341)ntime
        IF ( ntime >= 10 .AND. ntime .le.99 )            &
            write(fname1,342)ntime
        IF ( ntime >= 100 .AND. ntime .le. 999 )         &
            write(fname1,343)ntime
        IF ( ntime >= 1000 .AND. ntime .le. 9999 )       &
            write(fname1,344)ntime
        IF ( ntime >= 10000 .AND. ntime .le. 99999 )     &
            write(fname1,345)ntime
        IF ( ntime >= 100000 .AND. ntime .le. 999999 )   &
            write(fname1,346)ntime
        IF ( ntime >= 1000000 .AND. ntime .le. 9999999 ) &
            write(fname1,347)ntime
    341     format('q.000000',i1)
    342     format('q.00000',i2)
    343     format('q.0000',i3)
    344     format('q.000',i4)
    345     format('q.00',i5)
    346     format('q.0',i6)
    347     format('q.',i7)
	
	    IF (BinaryOutput==0) THEN
	      OPEN(UNIT=70,FILE=fname1,STATUS='UNKNOWN')
	    ELSE
	      OPEN(UNIT=70,FILE=fname1,form='unformatted')
	    END IF
	
        call vorticity()

        IF (format_dump == TECPLOT) THEN


        IF (nDim == DIM_2D) THEN
	      IF (AdditionalOutput == 0) THEN
		      write(70,*)'VARIABLES="X","Y","U","V","P","OZ","IBK_memb"'
		      write(70,*)'ZONE F=POINT, I=',nxc,', J=',nyc
		      k = 1
		      do j=1,nyc       ! j has been changed from 1,nyc to 2,ny-2 to output jum, jup
		      do i=1,nxc       ! i has been changed from 1,nyc to 2,ny-2 to output ium, iup
        
	            write(70,121)xc(i),yc(j),u(i,j,k),v(i,j,k),p(i,j,k) &
    !                    ,nlw(i,j,k), iblank(i,j,k),ghostCellMemb(i,j,k)+ & changed 05/27/10
    !                    ,nlw(i,j,k), iblank_memb(i,j,k),ghostCellMemb(i,j,k)+ &
    !                        ghostCellSolid(i,j,k)
    !                    ,nlw(i,j,k), iblank_memb(i,j,k)
                        ,nlw(i,j,k), iblank(i,j,k)
    121 format(6(2x,e14.7),1(2x,i2)) 
		      enddo
		      enddo
          ELSE
		    write(tmp,*) 'VARIABLES="X","Y","U","V","P","OZ"'
		    do ii=1, 14
		      jj=iand(output_para(ii),AdditionalOutput)
		      if (jj>0) then
		      select case (jj)
		      case (OUTPUT_IBLANK)
			    tmp=TRIM(tmp)//',"IBK"'
		      case (OUTPUT_FRESH_CELL)
			    tmp=TRIM(tmp)//',"FRESH"'
		      case (OUTPUT_GHOSTCELLMARK)
			    tmp=TRIM(tmp)//',"ghost"'
		      case (OUTPUT_BODYNUM)
			    tmp=TRIM(tmp)//',"BODYNUM"'
		      case (OUTPUT_IUM)
			    tmp=TRIM(tmp)//',"IUM"'
		      case (OUTPUT_IUP)
			    tmp=TRIM(tmp)//',"IUP"'
		      case (OUTPUT_JUM)
			    tmp=TRIM(tmp)//',"JUM"'
		      case (OUTPUT_JUP)
			    tmp=TRIM(tmp)//',"JUP"'
		      case (OUTPUT_KUM)
			    tmp=TRIM(tmp)//',"KUM"'
		      case (OUTPUT_KUP)
			    tmp=TRIM(tmp)//',"KUP"'
		      case (OUTPUT_GHOSTCELLSOLID)
			    tmp=TRIM(tmp)//',"GC_Solid"'
		      case (OUTPUT_GHOSTCELLMEMB)
			    tmp=TRIM(tmp)//',"GC_Memb"'
		      case (OUTPUT_IBLANK_SOLID)
			    tmp=TRIM(tmp)//',"IBK_Solid"'
		      case (OUTPUT_IBLANK_MEMB)
			    tmp=TRIM(tmp)//',"IBK_Memb"'
    !		  case ()
		      end select
		      end if
		    end do
		    write(70,'(A)') TRIM(tmp)
		    write(70,*)'ZONE F=POINT, I=',nxc,', J=',nyc
		
		    k = 1
		    do j=1,nyc       ! j has been changed from 1,nyc to 2,ny-2 to output jum, jup
		    do i=1,nxc       ! i has been changed from 1,nyc to 2,ny-2 to output ium, iup
        
            write(tmp,122)xc(i),yc(j),u(i,j,k),v(i,j,k),p(i,j,k),nlw(i,j,k)
    122 format(6(2x,e14.7)) 
		    do ii=1, 14
		      jj=iand(output_para(ii),AdditionalOutput)
		      if (jj>0) then
		      select case (jj)
		      case (OUTPUT_IBLANK)
			    write(tmp1,'(2x,i2)') iblank(i,j,k)
		      case (OUTPUT_FRESH_CELL)
			    write(tmp1,'(2x,i2)') fresh_cell(i,j,k)
		      case (OUTPUT_GHOSTCELLMARK)
			    write(tmp1,'(2x,i2)') ghostCellMark(i,j,k)
		      case (OUTPUT_BODYNUM)
			    write(tmp1,'(2x,i2)') BODYNUM(i,j,k)
		      case (OUTPUT_IUM)
			    write(tmp1,'(2x,i2)') int(ium(i,j,k))
		      case (OUTPUT_IUP)
			    write(tmp1,'(2x,i2)') int(iup(i,j,k))
		      case (OUTPUT_JUM)
			    write(tmp1,'(2x,i2)') int(jum(i,j,k))
		      case (OUTPUT_JUP)
			    write(tmp1,'(2x,i2)') int(jup(i,j,k))
		      case (OUTPUT_KUM)
			    write(tmp1,'(2x,i2)') int(kum(i,j,k))
		      case (OUTPUT_KUP)
			    write(tmp1,'(2x,i2)') int(kup(i,j,k))
		      case (OUTPUT_GHOSTCELLSOLID)
			    write(tmp1,'(2x,i2)') ghostCellSolid(i,j,k)
		      case (OUTPUT_GHOSTCELLMEMB)
			    write(tmp1,'(2x,i2)') ghostCellMemb(i,j,k)
		      case (OUTPUT_IBLANK_SOLID)
			    write(tmp1,'(2x,i2)') iblank_solid(i,j,k)
		      case (OUTPUT_IBLANK_MEMB)
			    write(tmp1,'(2x,i2)') iblank_memb(i,j,k)
    !		  case ()
		      end select
		      tmp=TRIM(tmp)//TRIM(tmp1)
		      end if
		    end do
		    write(70,'(A)') TRIM(tmp)
            enddo
	        enddo
	      END IF
        ELSE
	      IF (AdditionalOutput == 0) THEN
		    IF (BinaryOutput==1) THEN
		      write(70)'VARIABLES="X","Y","Z","U","V","W","P"'
		      write(70)'ZONE F=POINT, I=',nxc,', J=',nyc,' K=',nzc
		      do k=1,nzc
		      do j=1,nyc
		      do i=1,nxc
			     write(70)xc(i),yc(j),zc(k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)
		      enddo
		      enddo
		      enddo
		    ELSE
		      write(70,*)'VARIABLES="X","Y","Z","U","V","W","P"'
		      write(70,*)'ZONE F=POINT, I=',nxc,', J=',nyc,' K=',nzc
		      do k=1,nzc
		      do j=1,nyc
		      do i=1,nxc
			     write(70,123)xc(i),yc(j),zc(k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)
    123 format(7(2x,e14.7))
		      enddo
		      enddo
		      enddo
		    END IF
          ELSE
		    write(tmp,*) 'VARIABLES="X","Y","Z","U","V","W","P"'
		    do ii=1, 14
		      jj=iand(output_para(ii),AdditionalOutput)
		      if (jj>0) then
		      select case (jj)
		      case (OUTPUT_IBLANK)
			    tmp=TRIM(tmp)//',"IBK"'
		      case (OUTPUT_FRESH_CELL)
			    tmp=TRIM(tmp)//',"FRESH"'
		      case (OUTPUT_GHOSTCELLMARK)
			    tmp=TRIM(tmp)//',"ghost"'
		      case (OUTPUT_BODYNUM)
			    tmp=TRIM(tmp)//',"BODYNUM"'
		      case (OUTPUT_IUM)
			    tmp=TRIM(tmp)//',"IUM"'
		      case (OUTPUT_IUP)
			    tmp=TRIM(tmp)//',"IUP"'
		      case (OUTPUT_JUM)
			    tmp=TRIM(tmp)//',"JUM"'
		      case (OUTPUT_JUP)
			    tmp=TRIM(tmp)//',"JUP"'
		      case (OUTPUT_KUM)
			    tmp=TRIM(tmp)//',"KUM"'
		      case (OUTPUT_KUP)
			    tmp=TRIM(tmp)//',"KUP"'
		      case (OUTPUT_GHOSTCELLSOLID)
			    tmp=TRIM(tmp)//',"GC_Solid"'
		      case (OUTPUT_GHOSTCELLMEMB)
			    tmp=TRIM(tmp)//',"GC_Memb"'
		      case (OUTPUT_IBLANK_SOLID)
			    tmp=TRIM(tmp)//',"IBK_Solid"'
		      case (OUTPUT_IBLANK_MEMB)
			    tmp=TRIM(tmp)//',"IBK_Memb"'
    !		  case ()
		      end select
		      end if
		    end do
		    write(70,'(A)') TRIM(tmp)
		    write(70,*)'ZONE F=POINT, I=',nxc,', J=',nyc,' K=',nzc
		    do k=1,nzc
		    do j=1,nyc
		    do i=1,nxc
		      write(tmp,123)xc(i),yc(j),zc(k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)
		      do ii=1, 14
			    jj=iand(output_para(ii),AdditionalOutput)
			    if (jj>0) then
			    select case (jj)
			    case (OUTPUT_IBLANK)
			      write(tmp1,'(2x,i2)') iblank(i,j,k)
			    case (OUTPUT_FRESH_CELL)
			      write(tmp1,'(2x,i2)') fresh_cell(i,j,k)
			    case (OUTPUT_GHOSTCELLMARK)
			      write(tmp1,'(2x,i2)') ghostCellMark(i,j,k)
			    case (OUTPUT_BODYNUM)
			      write(tmp1,'(2x,i2)') BODYNUM(i,j,k)
			    case (OUTPUT_IUM)
			      write(tmp1,'(2x,i2)') int(ium(i,j,k))
			    case (OUTPUT_IUP)
			      write(tmp1,'(2x,i2)') int(iup(i,j,k))
			    case (OUTPUT_JUM)
			      write(tmp1,'(2x,i2)') int(jum(i,j,k))
			    case (OUTPUT_JUP)
			      write(tmp1,'(2x,i2)') int(jup(i,j,k))
			    case (OUTPUT_KUM)
			      write(tmp1,'(2x,i2)') int(kum(i,j,k))
			    case (OUTPUT_KUP)
			      write(tmp1,'(2x,i2)') int(kup(i,j,k))
			    case (OUTPUT_GHOSTCELLSOLID)
			      write(tmp1,'(2x,i2)') ghostCellSolid(i,j,k)
			    case (OUTPUT_GHOSTCELLMEMB)
			      write(tmp1,'(2x,i2)') ghostCellMemb(i,j,k)
			    case (OUTPUT_IBLANK_SOLID)
			      write(tmp1,'(2x,i2)') iblank_solid(i,j,k)
			    case (OUTPUT_IBLANK_MEMB)
			      write(tmp1,'(2x,i2)') iblank_memb(i,j,k)
      !		  case ()
			    end select
			    tmp=TRIM(tmp)//TRIM(tmp1)
			    end if
		      end do
		      write(70,'(A)') TRIM(tmp)
		    enddo
		    enddo
		    enddo

    !124 format(7(2x,e14.7),3(2x,i2))
          END IF
        ENDIF

    !     write(70,*)'VARIABLES="X","Y","U","V","W","P","RHS","IBLANK","FRESH"'
    !     write(70,*)'ZONE F=POINT, I=',nxc,', J=',nyc
    !     k=nz/2
    !     do j=1,nyc
    !     do i=1,nxc
    !        write(70,127)xc(i),yc(j),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k) &
    !                   ,nlu(i,j,k),iblank(i,j,k),fresh_cell(i,j,k)
    !     enddo
    !     enddo
    127 FORMAT(7(2x,e14.7),2(2x,i4))

      
    !     IF (nDim == DIM_2D .AND. internal_boundary_present == INTR_BOUND_PRESENT .AND. body_type == CANONICAL) THEN
    !
    !        nCylinder = 0
    !        DO n=1,nBody
    ! 	  IF (canonical_body_type(n) ==  ELLIPTIC_CYLINDER .OR. &
    !              canonical_body_type(n) ==  GENERAL_CYLINDER       ) nCylinder = nCylinder + 1 
    !        ENDDO
    !
    !        write(70,*) 'GEOMETRY T=LINE,C=BLACK,CS=GRID'
    !        write(70,*) nCylinder
    !
    !        DO n=1,nBody
    !          SELECT CASE (canonical_body_type(n))
    !            CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)
    !               write(70,*) nPtsBodyMarkerOrig(n)
    !               DO m=1,nPtsBodyMarkerOrig(n)
    !                  WRITE(70,*) xBodyMarker(n,m), yBodyMarker(n,m)
    !               ENDDO
    !          END SELECT
    !        ENDDO
    !
    !     ENDIF

        ENDIF ! tecplot
    end if  !dryRun

!   Marker points viz file
    IF ( boundary_motion == MOVING_BOUNDARY ) THEN
      IF ( ntime >= 0 .AND. ntime .le. 9  )            &
        write(fname2,441)ntime
      IF ( ntime >= 10 .AND. ntime .le.99 )            &
        write(fname2,442)ntime
      IF ( ntime >= 100 .AND. ntime .le. 999 )         &
        write(fname2,443)ntime
      IF ( ntime >= 1000 .AND. ntime .le. 9999 )       &
        write(fname2,444)ntime
      IF ( ntime >= 10000 .AND. ntime .le. 99999 )     &
        write(fname2,445)ntime
      IF ( ntime >= 100000 .AND. ntime .le. 999999 )   &
        write(fname2,446)ntime
      IF ( ntime >= 1000000 .AND. ntime .le. 9999999 ) &
        write(fname2,447)ntime
441     format('marker.000000',i1)
442     format('marker.00000',i2)
443     format('marker.0000',i3)
444     format('marker.000',i4)
445     format('marker.00',i5)
446     format('marker.0',i6)
447     format('marker.',i7)

      OPEN(UNIT=234,FILE=fname2,STATUS='UNKNOWN')

      WRITE(234,*)'TITLE="3D TRIANGULAR SURFACE DATA"'
      WRITE(234,*)'VARIABLES= "X","Y","Z"'
 
     DO iBody = 1, nBody
      WRITE(234,*)'ZONE T="unstruc"','N=',nPtsBodyMarker(iBody),'E=',totNumTriElem(iBody),'F=FEPOINT  ET=TRIANGLE'
       DO i=1,nPtsBodyMarker(iBody)
         write(234,*)xBodyMarker(iBody,i),yBodyMarker(iBody,i),zBodyMarker(iBody,i)
       ENDDO ! i
 
      DO  j=1,totNumTriElem(iBody)
        WRITE(234,*) triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
      ENDDO
 
    ENDDO
      CLOSE(234)

     IF (nDim == DIM_3D.AND.boundary_formulation==GCM_method) THEN
         CALL OUTPUT_surfacePressure
     END IF

    END IF ! boundary_motion
	
    IF (format_dump == FIELDVIEW) THEN
      fsmach = 0.0
      alp    = 0.0 
      density= 1.0
      WRITE(70,*)nx,ny,nz
      WRITE(70,*)fsmach,alp,re,time
      WRITE(70,*)(((density,i=1,nx),j=1,ny),k=1,nz ),   &
                 (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz),   & 
                 (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz),   & 
                 (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz),   & 
                 (((p(i,j,k),i=1,nx),j=1,ny),k=1,nz) 
    ENDIF


    CLOSE(70)

!    OPEN(UNIT=71,FILE='face_vel.dat',STATUS='UNKNOWN')
!
!    write(71,*)'VARIABLES="X","Y","Z","FACEU","FACEV","FACEW","IBLANK"'
!    write(71,*)'ZONE F=POINT, I=',nx,', J=',ny,' K=',nz
!    do k=1,nz
!    do j=1,ny
!    do i=1,nx
!       write(71,124)x(i),y(j),z(k),face_u(i,j,k),face_v(i,j,k),face_w(i,j,k), &
!         iblank(i,j,k)
!    enddo
!    enddo
!    enddo
!
!    CLOSE(71)

129 format(7(2x,e14.7),2(2x,i2))

    END SUBROUTINE write_dump


!------------------------------------------------    
! three components of vorticity in nlu, nlv, nlw
!------------------------------------------------    
    SUBROUTINE vorticity()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays

    IMPLICIT NONE
    
    INTEGER           :: i,j,k
    REAL(KIND=CGREAL) :: un,us,uf,ub
    REAL(KIND=CGREAL) :: ve,vw,vf,vb
    REAL(KIND=CGREAL) :: wn,ws,we,ww

    nlu(0:,0:,0:) = zero
    nlv(0:,0:,0:) = zero
    nlw(0:,0:,0:) = zero
    
    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      wn = ( fy(j+1)*w(i,j+1,k) + (oned-fy(j+1))*w(i,j,k)   )*(1-jup(i,j,k))   &
           + bcyw(i,j,k)*jup(i,j,k)

      ws = ( fy(j)  *w(i,j,k)   + (oned-fy(j)  )*w(i,j-1,k) )*(1-jum(i,j,k))   &
           + bcyw(i,j,k)*jum(i,j,k)

      vf = ( fz(k+1)*v(i,j,k+1) + (oned-fz(k+1))*v(i,j,k)   )*(1-kup(i,j,k))   &
           + bczv(i,j,k)*kup(i,j,k)

      vb = ( fz(k)  *v(i,j,k)   + (oned-fz(k)  )*v(i,j,k-1) )*(1-kum(i,j,k))   &
           + bczv(i,j,k)*kum(i,j,k)

      nlu(i,j,k) =(  ( wn - ws )*dyinv(j)       &
                    -( vf - vb )*dzinv(k) )*(1-iblank(i,j,k))
    ENDDO
    ENDDO
    ENDDO

    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      we = ( fx(i+1)*w(i+1,j,k) + (oned-fx(i+1))*w(i,j,k)   )*(1-iup(i,j,k))   &
           + bcxw(i,j,k)*iup(i,j,k)

      ww = ( fx(i)  *w(i,j,k)   + (oned-fx(i)  )*w(i-1,j,k) )*(1-ium(i,j,k))   &
           + bcxw(i,j,k)*ium(i,j,k)

      uf = ( fz(k+1)*u(i,j,k+1) + (oned-fz(k+1))*u(i,j,k)   )*(1-kup(i,j,k))   &
           + bczu(i,j,k)*kup(i,j,k)

      ub = ( fz(k)  *u(i,j,k)   + (oned-fz(k)  )*u(i,j,k-1) )*(1-kum(i,j,k))   &
           + bczu(i,j,k)*kum(i,j,k)

      nlv(i,j,k) =(  ( uf - ub )*dzinv(k)  &
                    -( we - ww )*dxinv(i) )*(1-iblank(i,j,k))
    ENDDO
    ENDDO
    ENDDO

    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      ve = ( fx(i+1)*v(i+1,j,k) + (oned-fx(i+1))*v(i,j,k)   )*(1-iup(i,j,k))   &
           + bcxv(i,j,k)*iup(i,j,k)

      vw = ( fx(i)  *v(i,j,k)   + (oned-fx(i)  )*v(i-1,j,k) )*(1-ium(i,j,k))   &
           + bcxv(i,j,k)*ium(i,j,k)

      un = ( fy(j+1)*u(i,j+1,k) + (oned-fy(j+1))*u(i,j,k)   )*(1-jup(i,j,k))   &
           + bcyu(i,j,k)*jup(i,j,k)

      us = ( fy(j)  *u(i,j,k)   + (oned-fy(j)  )*u(i,j-1,k) )*(1-jum(i,j,k))   &
           + bcyu(i,j,k)*jum(i,j,k)

      nlw(i,j,k) =(  ( ve - vw )*dxinv(i)  & 
                    -( un - us )*dyinv(j) )*(1-iblank(i,j,k))
      ! Changed by Wanh on 08/03/11
!      nlw(i,j,k) =(  ( ve - vw )*dxinv(i)*(1+ium(i,j,k)+iup(i,j,k))  &
!                    -( un - us )*dyinv(j)*(1+jum(i,j,k)+jup(i,j,k)) )*(oned-iblank(i,j,k))
    ENDDO
    ENDDO
    ENDDO

    END SUBROUTINE vorticity

!-------------------------------------------------------------------------------
   SUBROUTINE divergence()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    INTEGER              :: i,j,k
    REAL(KIND=CGREAL)    :: sum,divmax,divmin
    INTEGER , DIMENSION(1:3):: maxl,minl

! div = [ d(U)/dx + d(V)/dy + d(W)/dz ] 

    sum           = zero
    nlu(0:,0:,0:) = zero

    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      nlu(i,j,k)    = ( ( face_u(i+1,j,k) - face_u(i,j,k) )*dxinv(i)  &
                       +( face_v(i,j+1,k) - face_v(i,j,k) )*dyinv(j)  &
                       +( face_w(i,j,k+1) - face_w(i,j,k) )*dzinv(k) )
      nlu(i,j,k)    = nlu(i,j,k)*REAL(1-iblank(i,j,k),KIND=CGREAL)
      sum = sum + nlu(i,j,k)
    ENDDO
    ENDDO
    ENDDO

    nlu    = ABS(nlu)
    divmax = MAXVAL(nlu(1:nxc,1:nyc,1:nzc))
    maxl   = MAXLOC(nlu(1:nxc,1:nyc,1:nzc))
    WRITE(STDOUT,'(A,E15.7,A,3(2X,I4))') 'Max Divergence       = ',divmax,' at ',maxl

   END SUBROUTINE divergence


!-------------------------------------------------------------------------

   SUBROUTINE OUTPUT_surfacePressure() 
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

    IMPLICIT NONE

    INTEGER           :: i,j,k
    INTEGER           :: ibody,ifudrag
    INTEGER           :: iG, jG, kG, indZ
    INTEGER           :: nG,m ,m1, m2, mPointsOrig
    INTEGER           :: ii,jj,kk,n,iRow

    REAL(KIND=CGREAL) :: uIP, vIP, wIP,        &
                         uBI, vBI, wBI,        &
                         dUt1Dn,dUt2Dn,        &
                         alphaX,alphaY,alphaZ, &
                         dist, distMin, xp, yp, zp, xs, ys, zs
                              
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:) :: pBodyMarker_1, pBodyMarker_2

    CHARACTER*20              :: fname2
    
! allocate local arrays 

    write(*,*)'nBody, nPtsMax=', nBody, nPtsMax

    ALLOCATE(pBodyMarker_1(nBody, nPtsMax), pBodyMarker_2(nBody, nPtsMax))
               
! initialize variables

    DO iBody = 1, nBody

     DO m = 1,nPtsBodyMarker(iBody)

  ! find closest ghost node
             distMin = 1.0E8_CGREAL
             DO n=1,nGhost
               i = iGhost(n)
               j = jGhost(n)
               k = kGhost(n)
 
               if (iblank(i,j,k) == 0) Then
 
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
               end if
             ENDDO

  ! search vicinity of closest ghost node to find nodes surrounding element
             DO i = iG-3, iG+3
               IF ( ( xc(i)-xBodyMarker(ibody,m) )*(xc(i+1)-xBodyMarker(iBody,m) ) <= 0.0_CGREAL &
                                                .AND.   iblank(i,j,k) == 0) ii = i
             ENDDO
             DO j = jG-3, jG+3
               IF ( ( yc(j)-yBodyMarker(ibody,m) )*(yc(j+1)-yBodyMarker(iBody,m) ) <= 0.0_CGREAL &
                                                .AND.   iblank(i,j,k) == 0) jj = j
             ENDDO
             DO k = kG-3, kG+3
               IF ( ( zc(k)-zBodyMarker(ibody,m) )*(zc(k+1)-zBodyMarker(iBody,m) ) <= 0.0_CGREAL &
                                                .AND.   iblank(i,j,k) == 0) kk = k
             ENDDO
 
  ! trilinear interpolation
             alphaX = (xBodyMarker(iBody,m)-xc(ii))*dxcinv(ii+1)
             alphaY = (yBodyMarker(iBody,m)-yc(jj))*dycinv(jj+1)
             alphaZ = (zBodyMarker(iBody,m)-zc(kk))*dzcinv(kk+1)
             pBodyMarker_1(iBody,m)  &
                  = p(ii,jj,kk)*(1.0_CGREAL - alphaX)*(1.0_CGREAL - alphaY)*(1.0_CGREAL - alphaZ) &
              + p(ii+1,jj,kk)*alphaX*(1.0_CGREAL - alphaY)*(1.0_CGREAL - alphaZ)*(1.0_CGREAL-iblank(ii+1,jj,kk))    &
              + p(ii,jj+1,kk)*(1.0_CGREAL - alphaX)*alphaY*(1.0_CGREAL - alphaZ)*(1.0_CGREAL-iblank(ii,jj+1,kk))    &
             + p(ii,jj,kk+1)*(1.0_CGREAL - alphaX)*(1.0_CGREAL - alphaY)*alphaZ*(1.0_CGREAL-iblank(ii,jj,kk+1))  &   
              + p(ii+1,jj+1,kk)*alphaX*alphaY*(1.0_CGREAL - alphaZ)*(1.0_CGREAL-iblank(ii+1,jj+1,kk))       &
             + p(ii+1,jj,kk+1)*alphaX*(1.0_CGREAL - alphaY)*alphaZ*(1.0_CGREAL-iblank(ii+1,jj,kk+1))       & 
             + p(ii,jj+1,kk+1)*(1.0_CGREAL - alphaX)*alphaY*alphaZ*(1.0_CGREAL-iblank(ii,jj+1,kk+1))        & 
             + p(ii+1,jj+1,kk+1)*alphaX*alphaY*alphaZ*(1.0_CGREAL-iblank(ii+1,jj+1,kk+1))

       END DO ! end for m

    END DO ! end for ibody

!   Marker points viz file
      IF ( ntime >= 0 .AND. ntime .le. 9  )            &
        write(fname2,441)ntime
      IF ( ntime >= 10 .AND. ntime .le.99 )            &
        write(fname2,442)ntime
      IF ( ntime >= 100 .AND. ntime .le. 999 )         &
        write(fname2,443)ntime
      IF ( ntime >= 1000 .AND. ntime .le. 9999 )       &
        write(fname2,444)ntime
      IF ( ntime >= 10000 .AND. ntime .le. 99999 )     &
        write(fname2,445)ntime
      IF ( ntime >= 100000 .AND. ntime .le. 999999 )   &
        write(fname2,446)ntime
      IF ( ntime >= 1000000 .AND. ntime .le. 9999999 ) &
        write(fname2,447)ntime
441     format('markerP.000000',i1)
442     format('markerP.00000',i2)
443     format('markerP.0000',i3)
444     format('markerP.000',i4)
445     format('markerP.00',i5)
446     format('markerP.0',i6)
447     format('markerP.',i7)
 
      OPEN(UNIT=234,FILE=fname2,STATUS='UNKNOWN')
     DO iBody = 1, nBody
       WRITE(234,*)'ZONE T="unstruc"','N=',nPtsBodyMarker(iBody),'E=',totNumTriElem(iBody),'F=FEPOINT  ET=TRIANGLE'
       DO i=1,nPtsBodyMarker(iBody)
         write(234,*)xBodyMarker(iBody,i),yBodyMarker(iBody,i),zBodyMarker(iBody,i), pBodyMarker_1(ibody,i)
                      
       ENDDO ! i

       DO  j=1,totNumTriElem(iBody)
         WRITE(234,*) triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
       ENDDO

     ENDDO ! iBody
      CLOSE(234)

    DEALLOCATE(pBodyMarker_1, pBodyMarker_2)
    
   END SUBROUTINE OUTPUT_surfacePressure 
!-------------------------------------------------------------------------------

SUBROUTINE write_dump_debug(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  REAL(KIND=CGREAL) :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
write(70,*)'VARIABLES="X","Y","Z","U"'
write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc,', K=',nzc
!  WRITE(fname1,"(A,'.001.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=71,FILE=fname1,FORM='FORMATTED')
!write(71,*)'VARIABLES="X","Y","U"'
!write(71,*) 'ZONE F=POINT, I=',256,', J=',512
!  WRITE(fname1,"(A,'.002.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=72,FILE=fname1,FORM='FORMATTED')
!write(72,*)'VARIABLES="X","Y","U"'
!write(72,*) 'ZONE F=POINT, I=',256,', J=',512
!  WRITE(fname1,"(A,'.003.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=73,FILE=fname1,FORM='FORMATTED')
!write(73,*)'VARIABLES="X","Y","U"'
!write(73,*) 'ZONE F=POINT, I=',256,', J=',512

  do k=1,nzc
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,1PE22.15))') xc(i), yc(j), zc(k), var(i,j,k)
  END DO
  END DO
  enddo

!  write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
!  k=2
!  DO j=1,nyc
!  DO i=1,nxc
!      WRITE(70,'(2(3X,1PE12.5),(3X,1PE22.15))') xc(i), yc(j), var(i,j,k)
!  END DO
!  END DO
!
!  CLOSE(70)
!  CLOSE(71)
!  CLOSE(72)
!  CLOSE(73)

  RETURN
END SUBROUTINE

SUBROUTINE write_dump_debug2d(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  REAL(KIND=CGREAL) :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
write(70,*)'VARIABLES="X","Y","U"'
write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
!  WRITE(fname1,"(A,'.001.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=71,FILE=fname1,FORM='FORMATTED')
!write(71,*)'VARIABLES="X","Y","U"'
!write(71,*) 'ZONE F=POINT, I=',256,', J=',512
!  WRITE(fname1,"(A,'.002.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=72,FILE=fname1,FORM='FORMATTED')
!write(72,*)'VARIABLES="X","Y","U"'
!write(72,*) 'ZONE F=POINT, I=',256,', J=',512
!  WRITE(fname1,"(A,'.003.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
!  OPEN(UNIT=73,FILE=fname1,FORM='FORMATTED')
!write(73,*)'VARIABLES="X","Y","U"'
!write(73,*) 'ZONE F=POINT, I=',256,', J=',512

  k=2
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,1PE22.15))') xc(i), yc(j), var(i,j,k)
  END DO
  END DO

  CLOSE(70)
!  CLOSE(71)
!  CLOSE(72)
!  CLOSE(73)

  RETURN
END SUBROUTINE

!---------------------------------------------------------------------
SUBROUTINE write_dump_debug_i(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  INTEGER(1) :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  !print *,'output variable ', trim(vname)
  !WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  !OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  !write(70,*)'VARIABLES="Y","Z","U"'
  !
  !write(70,*) 'ZONE F=POINT, I=',nyc,', J=',nzc
  !do k=1,nzc
  !DO j=1,nyc
  !DO i=155,155
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') yc(j), zc(k), var(i,j,k)
  !END DO
  !END DO
  !end do
  
  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  write(70,*)'VARIABLES="X","Y","Z","U"'
  
  write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc,', K=',nzc
  do k=1,nzc
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(3(3X,1PE12.5),(3X,I8))') xc(i), yc(j), zc(k), var(i,j,k)
  END DO
  END DO
  end do
  
  !k=2
  !write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  !DO j=1,nyc
  !DO i=1,nxc
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,k)
  !END DO
  !END DO
  CLOSE(70)

  RETURN
    END SUBROUTINE
!---------------------------------------------------------------------
    
SUBROUTINE write_dump_debug_f(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  real*8 :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  !print *,'output variable ', trim(vname)
  !WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  !OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  !write(70,*)'VARIABLES="Y","Z","U"'
  !
  !write(70,*) 'ZONE F=POINT, I=',nyc,', J=',nzc
  !do k=1,nzc
  !DO j=1,nyc
  !DO i=155,155
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') yc(j), zc(k), var(i,j,k)
  !END DO
  !END DO
  !end do
  
  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  write(70,*)'VARIABLES="X","Y","Z","U"'
  
  write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc,', K=',nzc
  do k=1,nzc
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(4(3X,1PE12.5))') xc(i), yc(j), zc(k), var(i,j,k)
  END DO
  END DO
  end do
  
  !k=2
  !write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  !DO j=1,nyc
  !DO i=1,nxc
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,k)
  !END DO
  !END DO
  CLOSE(70)

  RETURN
    END SUBROUTINE
    
SUBROUTINE write_dump_debug_f_2D(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  real*8 :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  !print *,'output variable ', trim(vname)
  !WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  !OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  !write(70,*)'VARIABLES="Y","Z","U"'
  !
  !write(70,*) 'ZONE F=POINT, I=',nyc,', J=',nzc
  !do k=1,nzc
  !DO j=1,nyc
  !DO i=155,155
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') yc(j), zc(k), var(i,j,k)
  !END DO
  !END DO
  !end do
  
  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  write(70,*)'VARIABLES="X","Y","U"'
  
  write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(3(3X,1PE12.5))') xc(i), yc(j), var(i,j,1)
  END DO
  END DO
  
  !k=2
  !write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  !DO j=1,nyc
  !DO i=1,nxc
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,k)
  !END DO
  !END DO
  CLOSE(70)

  RETURN
END SUBROUTINE
    
SUBROUTINE write_dump_debug_i_2D(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  INTEGER(1) :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  !print *,'output variable ', trim(vname)
  !WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  !OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  !write(70,*)'VARIABLES="Y","Z","U"'
  !
  !write(70,*) 'ZONE F=POINT, I=',nyc,', J=',nzc
  !do k=1,nzc
  !DO j=1,nyc
  !DO i=155,155
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') yc(j), zc(k), var(i,j,k)
  !END DO
  !END DO
  !end do
  
  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  write(70,*)'VARIABLES="X","Y","U"'
  
  write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,1)
  END DO
  END DO
  
  !k=2
  !write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  !DO j=1,nyc
  !DO i=1,nxc
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,k)
  !END DO
  !END DO
  CLOSE(70)

  RETURN
END SUBROUTINE
    
SUBROUTINE write_dump_debug_i4_2D(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  INTEGER :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  !print *,'output variable ', trim(vname)
  !WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  !OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  !write(70,*)'VARIABLES="Y","Z","U"'
  !
  !write(70,*) 'ZONE F=POINT, I=',nyc,', J=',nzc
  !do k=1,nzc
  !DO j=1,nyc
  !DO i=155,155
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') yc(j), zc(k), var(i,j,k)
  !END DO
  !END DO
  !end do
  
  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  write(70,*)'VARIABLES="X","Y","U"'
  
  write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,1)
  END DO
  END DO
  
  !k=2
  !write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  !DO j=1,nyc
  !DO i=1,nxc
  !    WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,k)
  !END DO
  !END DO
  CLOSE(70)

  RETURN
END SUBROUTINE
    
SUBROUTINE write_dump_debug_i4(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  INTEGER :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  write(70,*)'VARIABLES="X","Y","U"'

  write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  k=1
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,k)
  END DO
  END DO
  k=2
  write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,I8))') xc(i), yc(j), var(i,j,k)
  END DO
  END DO
  CLOSE(70)

  RETURN
END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE write_dump_debug_i_2(vname, dbg,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  INTEGER(1) :: VAR(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
write(70,*)'VARIABLES="X","Y","Z","U"'
write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc,', K=',nzc

  DO k=1,nzc
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(3(3X,1PE12.5),(3X,I8))') xc(i), yc(j), zc(k), var(i,j,k)
  END DO
  END DO
  END DO

  CLOSE(70)

  RETURN
END SUBROUTINE


!---------------------------------------------------------------------
SUBROUTINE write_dump_debug_body(vname, dbg,iBody,var)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg, iBody
  REAL(KIND=CGREAL) :: VAR(nBody,nPtsMax)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  WRITE(70,*) 'VARIABLES="index","velocity"'

  DO i=1,nPtsBodyMarker(iBody)
!      WRITE(70,'(1x,i8), (1(3X,1PE12.5))') i, var(iBody,i)
!      WRITE(70,'(1X,i8), (2(3X,1PE12.5))') xc(i), yc(j), var(i,j,k)
     WRITE(70,'(1(3X,1PE12.5),(3X,I8))') var(iBody,i), i
  END DO

  CLOSE(70)

  RETURN
END SUBROUTINE

!======================================================================
SUBROUTINE write_dump_debug_body_vel(iBody)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays
  use body_dynamics
  use unstructured_surface_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg, iBody
  REAL(KIND=CGREAL) :: VAR(nBody,nPtsMax)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A8,I7.7,A1,I2.2,'.dat')") 'BodyVel_',ntime,'_',niterFS
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  WRITE(70,*) 'VARIABLES=x,y,z,u,v,w'
  write(70,*) 'ZONE T="unstruc"N=', nPtsBodyMarker(iBody)  ,' E=', totNumTriElem(iBody) ,' F=FEPOINT ET=TRIANGLE'

  DO i=1,nPtsBodyMarker(iBody)
!      WRITE(70,'(1x,i8), (1(3X,1PE12.5))') i, var(iBody,i)
!      WRITE(70,'(1X,i8), (2(3X,1PE12.5))') xc(i), yc(j), var(i,j,k)
     WRITE(70,'(6(1X,1PE12.5))') xbodymarker(iBody,i),ybodymarker(iBody,i),zBodyMarker(iBody,i), &
                                 ubodymarker(iBody,i),vbodymarker(iBody,i),wBodyMarker(iBody,i) 
  END DO

  DO  j=1,totNumTriElem(iBody)
      WRITE(70,*) triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
  ENDDO

  CLOSE(70)

  RETURN
END SUBROUTINE


SUBROUTINE write_dump_debug_mg1(vname, dbg,var,ilb,iub,jlb,jub,klb,kub)

  USE global_parameters
  USE flow_parameters
  USE grid_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, dbg,ilb,iub,jlb,jub,klb,kub
  REAL(KIND=CGREAL) :: VAR(ilb:iub,jlb:jub,klb:kub)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
  WRITE(70,*)'VARIABLES="X","Y","U"'
  WRITE(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc

  k=1
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,1PE19.11))') xc(i), yc(j), var(i,j,k)
  END DO
  END DO

  CLOSE(70)
!  CLOSE(71)
!  CLOSE(72)
!  CLOSE(73)

  RETURN
END SUBROUTINE
!=============================================
SUBROUTINE write_marker_vel(n1,n2)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE nlold_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    use fea_unstructure_surface

    integer :: ibody,n1,n2
    character fnametemp*100

    write(fnametemp,10)'marker_vel.',n1,'.',n2,'.dat'
10  format(A11,I7.7,A1,I3.3,A4)

      OPEN(UNIT=234,FILE=trim(fnametemp))

      WRITE(234,*)'TITLE="3D TRIANGULAR SURFACE DATA"'
      WRITE(234,*)'VARIABLES= "X","Y","Z","u","v","w"'
 
     DO iBody = 1, nBody
      WRITE(234,*)'ZONE T="unstruc"','N=',nPtsBodyMarker(iBody),'E=',totNumTriElem(iBody),'F=FEPOINT  ET=TRIANGLE'
       DO i=1,nPtsBodyMarker(iBody)
         write(234,*)xBodyMarker(iBody,i),yBodyMarker(iBody,i),zBodyMarker(iBody,i),&
                      uBodyMarker(iBody,i),vBodyMarker(iBody,i),wBodyMarker(iBody,i)
       ENDDO ! i
 
      DO  j=1,totNumTriElem(iBody)
        WRITE(234,*) triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
      ENDDO
 
    ENDDO
      CLOSE(234)


END SUBROUTINE

!===================================================
SUBROUTINE write_dump_debug2d_3vars(vname,dbg,var1,var2,var3)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  REAL(KIND=CGREAL) :: VAR1(0:nx+1,0:ny+1,0:nz+1),VAR2(0:nx+1,0:ny+1,0:nz+1),VAR3(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
write(70,*)'VARIABLES="X","Y","U","V","P"'
write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc


  k=2
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,1PE22.15))') xc(i), yc(j), var1(i,j,k),var2(i,j,k),var3(i,j,k)
  END DO
  END DO

  CLOSE(70)


  RETURN
END SUBROUTINE write_dump_debug2d_3vars

!===================================================
SUBROUTINE write_dump_debug2d_2vars(vname,dbg,var1,var2)

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE grid_arrays
  USE pressure_arrays
  USE boundary_arrays
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i, j, k, iG, jG, dbg
  REAL(KIND=CGREAL) :: VAR1(0:nx+1,0:ny+1,0:nz+1),VAR2(0:nx+1,0:ny+1,0:nz+1)

  CHARACTER*30 :: fname1
  CHARACTER*4 :: vname

  print *,'output variable ', trim(vname)
  WRITE(fname1,"(A,'.',I7.7,'.',I2.2,'.dat')") trim(vname),ntime,dbg
  OPEN(UNIT=70,FILE=fname1,FORM='FORMATTED')
write(70,*)'VARIABLES="X","Y","U","V"'
write(70,*) 'ZONE F=POINT, I=',nxc,', J=',nyc


  k=2
  DO j=1,nyc
  DO i=1,nxc
      WRITE(70,'(2(3X,1PE12.5),(3X,1PE22.15))') xc(i), yc(j), var1(i,j,k),var2(i,j,k)
  END DO
  END DO

  CLOSE(70)


  RETURN
END SUBROUTINE write_dump_debug2d_2vars