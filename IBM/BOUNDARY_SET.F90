!---------------------------------------------
!   SUBROUTINE set_boundary()
!   SUBROUTINE set_internal_boundary()
!   SUBROUTINE set_iblank_canonical_body()
!---------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE set_boundary()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,m, ibody
    REAL(KIND=CGREAL) :: rbcx1,rbcx2,rbcy1,rbcy2,rbcz1,rbcz2
    REAL(KIND=CGREAL) :: facScalInv

    iup = 0
    jup = 0
    kup = 0
    ium = 0
    jum = 0
    kum = 0

    iMarkp = 0
    jMarkp = 0
    kMarkp = 0
    iMarkm = 0
    jMarkm = 0
    kMarkm = 0

! 	set outer boundary
!		------------------
    WRITE(STDOUT,'(5X,A)') 'Entering set_outer_iup'
    CALL set_outer_iup()

!   sets internal boundary
!   ----------------------

    IF ( boundary_formulation /= NO_INTERNAL_BOUNDARY ) THEN
      WRITE(STDOUT,'(5X,A)') 'Entering set_internal_boundary....'
      CALL set_internal_boundary()
    ENDIF ! boundary_formulation

    WRITE(STDOUT,'(5X,A)') 'set_outflow_area'
    CALL set_outflow_area()

   END SUBROUTINE set_boundary
!------------------------------------------------------------------------------

SUBROUTINE set_outer_iup()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k

! outer boundary

    DO j=1,ny
    DO k=1,nz
      ium(1,j,k)    = 1
      iup(nxc,j,k) = 1
    ENDDO
    ENDDO

    DO i=1,nx
    DO k=1,nz
      jum(i,1,k)    = 1
      jup(i,nyc,k) = 1
    ENDDO
    ENDDO

    IF (nDim == DIM_3D) THEN
      DO i=1,nx
      DO j=1,ny
        kum(i,j,1)    = 1
        kup(i,j,nzc) = 1
      ENDDO
      ENDDO
    END IF

    IF (Hybrid) THEN
      iupp = 0    !
      jupp = 0    !
      kupp = 0    !   Added by Rupesh (used for 2nd upwinding)
      iumm = 0    !
      jumm = 0    !
      kumm = 0    !

      DO j=1,ny
      DO k=1,nz
        iumm(1,j,k)    = 1  !
        iupp(nxc,j,k) = 1  !   Added by Rupesh (used for 2nd upwinding)
        iumm(2,j,k)    = 1  !
        iupp(nx-2,j,k) = 1  !
      ENDDO
      ENDDO

      DO i=1,nx
      DO k=1,nz
        jumm(i,1,k)    = 1  !
        jupp(i,nyc,k) = 1  !   Added by Rupesh (used for 2nd upwinding)
        jumm(i,2,k)    = 1  !
        jupp(i,ny-2,k) = 1  !

      ENDDO
      ENDDO

      IF (nDim == DIM_3D) THEN
        DO i=1,nx
        DO j=1,ny
          kumm(i,j,1)    = 1  !
          kupp(i,j,nzc) = 1  !   Added by Rupesh (used for 2nd upwinding)
          kumm(i,j,2)    = 1  !
          kupp(i,j,nz-2) = 1  !
        ENDDO
        ENDDO
      END IF

    ENDIF

END SUBROUTINE set_outer_iup
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE set_outflow_area()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k
    REAL(KIND=CGREAL) :: rbcx1,rbcx2,rbcy1,rbcy2,rbcz1,rbcz2
    REAL(KIND=CGREAL) :: facScalInv

! compute areas of external boundaries

    area_left  = zero
    area_right = zero
    area_bot   = zero
    area_top   = zero
    area_back  = zero
    area_front = zero

    DO j=1,nyc
    DO k=1,nzc
      area_left  = area_left  + dy(j)*dz(k)*(1-iblank(1,   j,k) )
      area_right = area_right + dy(j)*dz(k)*(1-iblank(nxc,j,k) )
    ENDDO
    ENDDO

    DO i=1,nxc
    DO k=1,nzc
      area_bot   = area_bot   + dx(i)*dz(k)*(1-iblank(i,1,   k) )
      area_top   = area_top   + dx(i)*dz(k)*(1-iblank(i,nyc,k) )
    ENDDO
    ENDDO

    DO j=1,nyc
    DO i=1,nxc
      area_back  = area_back  + dx(i)*dy(j)*(1-iblank(i,j,1   ) )
      area_front = area_front + dx(i)*dy(j)*(1-iblank(i,j,nzc) )
    ENDDO
    ENDDO

    rbcx1 = REAL(bcx1,KIND=CGREAL)
    rbcx2 = REAL(bcx2,KIND=CGREAL)
    rbcy1 = REAL(bcy1,KIND=CGREAL)
    rbcy2 = REAL(bcy2,KIND=CGREAL)
    rbcz1 = REAL(bcz1,KIND=CGREAL)
    rbcz2 = REAL(bcz2,KIND=CGREAL)

    facScalInv = oned/ &
                (oned)/(-oned)/(-2.0_CGREAL)/&
                (-3.0_CGREAL)/(-4.0_CGREAL)/(-5.0_CGREAL)

! calculate outflow area without using IF statements
    outflow_area = (rbcx1-oned)*(rbcx1-3.0_CGREAL)* &
                   (rbcx1-4.0_CGREAL)*(rbcx1-5.0_CGREAL)* &
                   (rbcx1-6.0_CGREAL)*(rbcx1-7.0_CGREAL)* &
                    facScalInv * area_left                &

                  +(rbcx2-oned)*(rbcx2-3.0_CGREAL)* &
                   (rbcx2-4.0_CGREAL)*(rbcx2-5.0_CGREAL)* &
                   (rbcx2-6.0_CGREAL)*(rbcx2-7.0_CGREAL)* &
                    facScalInv * area_right               &

                  +(rbcy1-oned)*(rbcy1-3.0_CGREAL)* &
                   (rbcy1-4.0_CGREAL)*(rbcy1-5.0_CGREAL)* &
                   (rbcy1-6.0_CGREAL)*(rbcy1-7.0_CGREAL)* &
                    facScalInv * area_bot                 &

                  +(rbcy2-oned)*(rbcy2-3.0_CGREAL)* &
                   (rbcy2-4.0_CGREAL)*(rbcy2-5.0_CGREAL)* &
                   (rbcy2-6.0_CGREAL)*(rbcy2-7.0_CGREAL)* &
                    facScalInv * area_top                 &

                  +(rbcz1-oned)*(rbcz1-3.0_CGREAL)* &
                   (rbcz1-4.0_CGREAL)*(rbcz1-5.0_CGREAL)* &
                   (rbcz1-6.0_CGREAL)*(rbcz1-7.0_CGREAL)* &
                    facScalInv * area_back                &

                  +(rbcz2-oned)*(rbcz2-3.0_CGREAL)* &
                   (rbcz2-4.0_CGREAL)*(rbcz2-5.0_CGREAL)* &
                   (rbcz2-6.0_CGREAL)*(rbcz2-7.0_CGREAL)* &
                    facScalInv * area_front

END SUBROUTINE set_outflow_area
!---------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE set_initial_iblank()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k
    INTEGER :: ifuGrid

! by pass setting iblank only at initialization stage after restart

    IF(readIblankFlag) GOTO 1000

! Read iblank for general body

    IF ( body_type == GENERAL) CALL set_iblank_general_body()

! initialize iblank for canonical bodies

!    IF (body_type == CANONICAL) CALL set_iblank_canonical_body()

! write file

1000 CONTINUE

    IF (format_dump == FIELDVIEW) THEN
      ifuGrid = 265
      OPEN(UNIT=ifuGrid,FILE='grid_fieldview.dat')
      WRITE(ifuGrid,*)nx,ny,nz
      WRITE(ifuGrid,*)(((x(i),i=1,nx),j=1,ny),k=1,nz)    &
                     ,(((y(j),i=1,nx),j=1,ny),k=1,nz)    &
                    ,(((z(k),i=1,nx),j=1,ny),k=1,nz)     &
                    ,(((1-iblank(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    ENDIF

   END SUBROUTINE set_initial_iblank
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE set_internal_boundary()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,iBody,iErr
    INTEGER(1) :: iblankOld(0:nx+1,0:ny+1,0:nz+1)
! set iblank at current time step

    ghostCellMark = 0
    fresh_cell = 0
    num_iblank = 0

    conflictCell=0
    conflictBCi=0
    conflictBCj=0
    conflictBCk=0
    WRITE(*,*)  'start set_internal_boundary ...'
    !body_type=CANONICAL
    SELECT CASE(body_type)
    CASE(CANONICAL)

      IF(nBody_membrane /=0) THEN
        WRITE(*,*) 'Setting up for membranes...'

        iblankOld=iblank_memb
        iblank_memb=0
! ghostCellMemb stores the info of ghost cell for membrane in last time step

        DO iBody=nBody_solid+1, nBody

          iblank = iblankOld

          CALL set_iblank_canonical_body_fast(iBody,iBody)

!!! The fresh cell for membrance is set in identify_ghostcells_membrane
          CALL identify_ghostcells_membrane(iBody)

          CALL find_IsolatedIblank(1)
          CALL find_iblankHoles(1)
!
!        CALL write_dump_debug_i('GCM ',iBody,ghostCellMark)
!        CALL write_dump_debug_i('ibme',iBody,iblank_memb)

        ENDDO

        IF (boundary_formulation == SSM_METHOD) THEN
          CALL SSM_set_internal_iup_membrane()
		    END IF

        ghostCellMemb = ghostCellMark

        iblank=0

!        ghostCellMark = 0

      END IF ! membrane

      IF(nBody_solid /=0) THEN
        WRITE(STDOUT,'(5X,A)') 'Setting up for solid bodies...'
        WRITE(*,*) 'Setting up for solid bodies...'

        iblank=iblank_solid
        if(channel_flow)then

            CALL set_iblank_canonical_body_fast(1, 1)

            CALL find_IsolatedIblank(1)
            CALL find_iblankHoles(1)

            do k=1,nzc
            do j=1,nyc
            do i=1,nxc
                if(iblank(i,j,k)==0)then
                    iblank(i,j,k)=1
                else if(iblank(i,j,k)==1)then
                    iblank(i,j,k)=0
                end if

            end do
            end do
            end do
            iblank_solid = iblank
            if(nBody_Solid>1)then
                CALL set_iblank_canonical_body_fast(2, nBody_Solid)
                CALL find_IsolatedIblank(1)
                CALL find_iblankHoles(1)
                iblank=iblank+iblank_solid
            end if
        else

            CALL set_iblank_canonical_body_fast(1, nBody_Solid)

            CALL find_IsolatedIblank(1)
            CALL find_iblankHoles(1)
        end if

        WRITE(STDOUT,'(5X,A)') 'call fill_cell...'
        call fill_cell()

        iblank_solid = iblank

		    IF (boundary_formulation == SSM_METHOD) THEN
              CALL identify_ghostcells_solid()
		      CALL SSM_set_internal_iup_solid()
		      CALL SSM_set_internal_area()
		    END IF

    		IF (cure_pressure_oscillations)  CALL calc_area_n_volumn()

      END IF ! solid_body

      WRITE(*,*) 'call identify_conflict_cell...'
      call identify_conflict_cell

      !CALL write_dump_debug_i4_2D('bdNm',1,bodyNum)                !yan_dbg
      !CALL write_dump_debug_i4_2D('conf',1,conflictCell)           !yan_dbg
      !CALL write_dump_debug_i_2D('iblk',1,iblank)                  !yan_dbg
      !CALL write_dump_debug_i4_2D('ibci',1,conflictBCi)            !yan_dbg
      !CALL write_dump_debug_i4_2D('ibcj',1,conflictBCj)            !yan_dbg
      if(ibkOut)then
          if(mod(nTime,ndump)==0)then
            if(ndim==dim_3D) CALL write_dump_debug_i('iblk',1,iblank)            !yan_dbg
            if(ndim==dim_2D) CALL write_dump_debug_i_2D('iblk',1,iblank)
          end if
      end if
!      CALL identify_freshcells()

      iblank = iblank_solid

      ghostCellMark  = 0
      IF ( boundary_formulation == GCM_METHOD ) THEN
        CALL GCM_set_internal_boundary()
      END IF ! boundary_formulation

    CASE(GENERAL)
      CALL set_iblank_general_body()

    END SELECT ! body_type
   END SUBROUTINE set_internal_boundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE set_iblank_general_body()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k

      DO k=1,nzc
      DO j=1,nyc
      DO i=1,nxc
        READ(52,*)iblank(i,j,k)
      ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE set_iblank_general_body
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------

!-------------------------------------------------------
   SUBROUTINE find_isolatedIblank(fix)
!-------------------------------------------------------
! Checking for isolated iblanked cells and eliminating  then
! such isolated cells can be formed due to roundoff error

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: fix

    INTEGER              :: i,j,k,scanRange,sumIblank
    INTEGER              :: ii,jj,kk
    INTEGER              :: iMin,iMax,jMin,jMax,kMin,kMax

    DO k=1,nzc
    DO j=1,nyc
    DO i=1,nxc

      IF (iBlank(i,j,k) == 1) THEN
        scanRange = 1                 ! how far to search to detect boundary
        iMin = MAX(1   ,i-scanRange)
        iMax = MIN(nxc,i+scanRange)
        jMin = MAX(1   ,j-scanRange)
        jMax = MIN(nyc,j+scanRange)
        kMin = MAX(1   ,k-scanRange)
        kMax = MIN(nzc,k+scanRange)

        sumIblank = 0
        DO kk = kMin,kMax
        DO jj = jMin,jMax
        DO ii = iMin,iMax
          sumIblank = sumIblank + iBlank(ii,jj,kk)
        ENDDO
        ENDDO
        ENDDO

!        IF (body_dim(bodyNum(i,j,k)) == BODY_DIM2) sumIblank = sumIblank/2

        IF (sumIblank <= 1) THEN
          PRINT*,'Found Isolated Iblanked Cell at ',i,j,k
          IF (fix == 1) THEN
            PRINT*,'Iblanking removed for this cell'
            iBlank(i,j,k)  = 0
            bodyNum(i,j,k) = 0
          ENDIF
        ENDIF
      ENDIF

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

   END SUBROUTINE find_isolatedIblank

!------------------------------------------------------------------------------
   SUBROUTINE find_iblankHoles(fix)
!------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: fix

    INTEGER              :: i,j,k

    logical :: hole

! finding holes

    hole=.true.
    IF ( boundary_formulation == SSM_METHOD ) THEN

      do while(hole)
      DO k=1,nzc
      DO j=1,nyc
      DO i=1,nxc

        IF (iBlank(i,j,k) == 0 ) THEN

          IF ( (iblank(i+1,j,k) == 1 .OR. i == nxc) .AND. (iblank(i-1,j,k) == 1 .OR. i == 1) ) THEN
           PRINT*,'flow hole-1 at',i,j,k
           IF (fix == 1) THEN
             PRINT*,'closing hole'
             iblank(i,j,k) = 1
             bodyNum(i,j,k)= bodynum(i,j,k)
             IF (fresh_cell(i,j,k) == 1) THEN
               fresh_cell(i,j,k) = 0
!              num_fresh = num_fresh-1
             ENDIF
           ENDIF
          ENDIF

          IF ( (iblank(i,j+1,k) == 1 .OR. j == nyc) .AND. (iblank(i,j-1,k) == 1 .OR. j == 1) ) THEN
           PRINT*,'flow hole-2 at',i,j,k
           IF (fix == 1) THEN
             PRINT*,'closing hole'
             iblank(i,j,k) = 1
             bodyNum(i,j,k)= bodynum(i,j,k)
             IF (fresh_cell(i,j,k) == 1) THEN
               fresh_cell(i,j,k) = 0
!              num_fresh = num_fresh-1
             ENDIF
           ENDIF
          ENDIF

          IF ( (iblank(i,j,k+1) == 1 .OR. k == nzc) .AND. (iblank(i,j,k-1) == 1 .OR. k == 1) ) THEN
           PRINT*,'flow hole-3 at',i,j,k
           IF (fix == 1) THEN
             PRINT*,'closing hole'
             iblank(i,j,k) = 1
             bodyNum(i,j,k)= bodynum(i,j,k)
             IF (fresh_cell(i,j,k) == 1) THEN
               fresh_cell(i,j,k) = 0
!              num_fresh = num_fresh-1
             ENDIF
           ENDIF
          ENDIF

        ENDIF

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k
      hole=.false.
          DO k=1,nzc
          DO j=1,nyc
          DO i=1,nxc
              IF (iBlank(i,j,k) == 0 ) THEN
                    IF ( ((iblank(i+1,j,k) == 1 .OR. i == nxc) .AND. (iblank(i-1,j,k) == 1 .OR. i == 1)) .or.&
                         ((iblank(i,j+1,k) == 1 .OR. j == nyc) .AND. (iblank(i,j-1,k) == 1 .OR. j == 1)) .or.&
                         ((iblank(i,j,k+1) == 1 .OR. k == nzc) .AND. (iblank(i,j,k-1) == 1 .OR. k == 1))) THEN

                         hole=.true.
                         goto 888
                    end if


              end if

          end do
          end do
          end do

888 continue
    end do

    ENDIF

   END SUBROUTINE find_iblankHoles

!------------------------------------------------------------------------------
