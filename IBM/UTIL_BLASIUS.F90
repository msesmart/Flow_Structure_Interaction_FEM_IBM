        SUBROUTINE BLASIUS_VELOCITY

        USE flow_parameters,ONLY : nx,ny,nz
        USE grid_arrays, ONLY : x,y
        USE boundary_arrays,ONLY : iblank
        USE blasius_profile

        IMPLICIT NONE
!
!... Local Variables
!
        INTEGER :: i,j,k,j_surf,j_slot
        REAL(KIND=CGREAL) :: d_start,d_end

!... Allocation for Blasius Velocity Profiles

       ALLOCATE (eta(0:ny+1))
       ALLOCATE (u_blasius(0:ny+1))


        OPEN (9,file='vprof.dat')
        READ(9,*) ddratio             ! delta-d ratio
        READ(9,*) uinf                ! free-stream Vel
        CLOSE(9)
!
!... Determine cavity height and slot height using iblank file
!
        i = 5
        k = IDNINT(REAL(nz,KIND=CGREAL)/2.0_CGREAL)
         DO j = ny,1,-1
           IF ( (iblank(i,j,k) .EQ. 0) .AND. (iblank(i,j-1,k) .EQ. 1) ) THEN
             j_surf = j
           ENDIF
           IF ( (iblank(i,j,k) .EQ. 1) .AND. (iblank(i,j-1,k) .EQ. 0) ) THEN
             j_slot = j
           ENDIF
         ENDDO
        cavity_H = y(j_slot) - y(1)
        slot_H   = y(j_surf) - y(j_slot)
!
!... Determine Slot width
!
        k = IDNINT(REAL(nz,KIND=CGREAL)/2.0_CGREAL)
        DO i = 1,nx-2
         IF ( (iblank(i,j_surf-1,k) .EQ. 1) .AND.             &
              (iblank(i+1,j_surf-1,k)) .EQ. 0 ) THEN 
            d_start = x(i)
         ENDIF
         IF ( (iblank(i,j_surf-1,k) .EQ. 0) .AND.             &
              (iblank(i+1,j_surf-1,k)) .EQ. 1 ) THEN
            d_end = x(i)
         ENDIF
        ENDDO

        d = d_end - d_start
!
!...  Determine delta (BL thickness) based on ddratio
!
        delta = d*ddratio

        DO l = 1,ny
         IF ( (y(l)-(cavity_H+slot_H)) .LT. 1.0E-06 ) THEN
          i_start = l
         ENDIF
        ENDDO

        OPEN(12,file='blasius_u.dat')
        DO l =1,ny
         IF ( y(l) .LT. cavity_H+slot_H) THEN
          u_blasius(l) = zero
         ELSE
          IF ( y(l) .GE. (cavity_H+slot_H) .AND. y(l) .LT.       &
              (cavity_H+slot_H+delta)) THEN
           eta(l) = ( y(l)-y(i_start) )/delta
           u_blasius(l) = uinf*( 3.0_CGREAL*eta(l)/2.0_CGREAL -  &
                                 eta(l)**3/2.0_CGREAL )
          ELSE
           u_blasius(l) = uinf
          ENDIF
         ENDIF
         WRITE(12,*) y(l),u_blasius(l)
        ENDDO
        CLOSE(12)

!        DEALLOCATE (eta)
!        DEALLOCATE (u_blasius)

        END SUBROUTINE BLASIUS_VELOCITY
