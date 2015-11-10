!---------------------------------------
!  SUBROUTINE open_probe_files()
!  SUBROUTINE read_probe_inputs()
!  SUBROUTINE write_probe_files()
!---------------------------------------
!
!---------------------------------------------

   SUBROUTINE open_probe_files()
    
    USE probe_parameters
    USE flow_parameters
    USE grid_arrays
    
    IMPLICIT NONE

    CHARACTER*9          :: probeFile
    CHARACTER*25         :: inProbeFile

    INTEGER :: m

!   Open files for Probe
 
    DO m = 1, nProbe
      probeFile = TRIM("probe_out")
      WRITE(inProbeFile,101) probeFile,m
      OPEN(UNIT=ifuProbeOut+m-1,FILE=inProbeFile,FORM='formatted',ACTION="WRITE")

!   Write Header

      WRITE(ifuProbeOut+m-1,1000) &
            iProbe(m),jProbe(m), kProbe(m), &
             x(iProbe(m)),y(jProbe(m)),z(kProbe(m))
    END DO ! m

!   formats

101  FORMAT(a,'_',i3.3,'.dat')
1000 FORMAT('# probe data (time, u, v, w, p)',/, &
            '# icell ',I5,', jcell ',I5,', kcell ',I5,/, &
            '# x=',E13.5,', y=',E13.5,', z=',E13.5)

   END SUBROUTINE open_probe_files
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE read_probe_inputs()

    USE flow_parameters
    USE probe_parameters

    IMPLICIT NONE

    INTEGER :: m

    OPEN(ifuProbeIn,FILE='probe_in.dat',STATUS='UNKNOWN')
    READ(ifuProbeIn,*)nProbe
    PRINT*,'   nProbe = ',nProbe

    ALLOCATE(iProbe(nProbe))
    ALLOCATE(jProbe(nProbe))
    ALLOCATE(kProbe(nProbe))
    
    PRINT*,'Reading probe_in.dat File'
    DO m= 1, nProbe
      READ(ifuProbeIn,*)iProbe(m),jProbe(m),kProbe(m)   
    ENDDO ! m

   END SUBROUTINE read_probe_inputs
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE write_probe_files()

    USE global_parameters
    USE flow_parameters
    USE probe_parameters
    USE flow_arrays
    USE pressure_arrays

    IMPLICIT NONE
    
    INTEGER :: i,j,k,m
    REAL(KIND=CGREAL) :: uProbe,vProbe,wProbe,pProbe
    
    DO m= 1, nProbe
      i = iProbe(m)
      j = jProbe(m)
      k = kProbe(m)
      
      uProbe = u(i,j,k)
      vProbe = v(i,j,k)
      wProbe = w(i,j,k)
      pProbe = p(i,j,k)
      
      WRITE(ifuProbeOut+m-1,1001) time,uProbe,vProbe,wProbe,pProbe
    ENDDO ! m  

1001 FORMAT(1PE14.7,4E18.7)

   END SUBROUTINE write_probe_files
!-------------------------------------------   
