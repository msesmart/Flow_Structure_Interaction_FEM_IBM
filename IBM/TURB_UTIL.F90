!******************************************************************************
!
! Purpose: utilities to write Tecplot file for turbulence models
!
! Description: none.
!
! Input: viscTurb = turbulent viscosity
!
! Output: tecplot files, min-max values.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: Exp $
!
! Copyright: (c) 2004 by the George Washington University
!
!******************************************************************************

!------------------------------------------------------------------------------   
    SUBROUTINE TURB_write_monitor() 

!==============================================================================
!  Purpose: Write min-max values to standard output
!==============================================================================

    USE global_parameters
    USE turb_global_parameters
    USE flow_parameters
    USE flow_arrays
    USE turb_parameters
    USE turb_arrays
    
    IMPLICIT NONE

    WRITE(STDOUT,'(A,2(2X,1PE15.7))') 'Min-Max of viscTot = ',&
       MINVAL(viscTot(1:nx-1,1:ny-1,1:nz-1)),MAXVAL(viscTot(1:nx-1,1:ny-1,1:nz-1))

    WRITE(STDOUT,'(A,2(2X,1PE15.7))') 'Min-Max of S_ii     = ',&
      MINVAL(sij(S11,1:nx-1,1:ny-1,1:nz-1)&
            +sij(S22,1:nx-1,1:ny-1,1:nz-1)&
            +sij(S33,1:nx-1,1:ny-1,1:nz-1)),&
      MAXVAL(sij(S11,1:nx-1,1:ny-1,1:nz-1)&
            +sij(S22,1:nx-1,1:ny-1,1:nz-1)&
            +sij(S33,1:nx-1,1:ny-1,1:nz-1)) 

    END SUBROUTINE TURB_write_monitor
!------------------------------------------------------------------------------ 
 
    SUBROUTINE TURB_write_dump() 

!==============================================================================
!  Purpose: Write dump files of the turbulent viscosity 
!          in tecplot or fieldview fornat
!==============================================================================

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE turb_parameters
    USE turb_arrays
    
    IMPLICIT NONE

!... Loop variables
    INTEGER      :: i,j,k
    
!... Local variables

    CHARACTER*13 :: fname1   
    REAL         :: fsmach,alp,density 

!******************************************************************************

    PRINT*,'Writing out dump turbulence file'
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
341     format('nut.000000',i1)
342     format('nut.00000',i2)
343     format('nut.0000',i3)
344     format('nut.000',i4)
345     format('nut.00',i5)
346     format('nut.0',i6)
347     format('nut.',i7)

    OPEN(UNIT=70,FILE=fname1,STATUS='UNKNOWN')

    SELECT CASE (format_dump)
      CASE(TECPLOT)
        WRITE(70,*)'VARIABLES="X","Y","Z","NUTOT","IBLANK","FRESH"'
        WRITE(70,*)'ZONE F=POINT, I=',nx-1,', J=',ny-1,' K=',nz-1
        DO k=1,nz-1
        DO j=1,ny-1
        DO i=1,nx-1
         write(70,123)xc(i),yc(j),zc(k),viscTot(i,j,k),iblank(i,j,k),fresh_cell(i,j,k) 
        ENDDO
        ENDDO
        ENDDO
      
      CASE(FIELDVIEW)
        fsmach = 0.0
        alp    = 0.0 
        density= 1.0 
        WRITE(70,*)nx,ny,nz
        WRITE(70,*)fsmach,alp,re,time
        WRITE(70,*)(((density        ,i=1,nx),j=1,ny),k=1,nz),   &
                   (((viscTot(i,j,k),i=1,nx),j=1,ny),k=1,nz)
 
    END SELECT ! format_dump 

    CLOSE(70)

123 FORMAT(4(2x,e14.7),2(2x,i2))

    END SUBROUTINE TURB_write_dump
!------------------------------------------------------------------------------ 
