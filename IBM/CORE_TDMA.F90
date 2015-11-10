!--------------------------------------------------------
!  SOLVE THE TRIDIAGONAL MATRIX 
!    [b1  c1           ]
!    [a2  b2  c2       ]
!    [...              ]
!    [...              ]
!    [           aN  bN]
!--------------------------------------------------------

   SUBROUTINE tdma(a,b,c,r,u,n1,n2)

    USE global_parameters    
    IMPLICIT NONE

    INTEGER                           , INTENT (IN)  :: n1,n2
    REAL(KIND=CGREAL), DIMENSION(0:n2), INTENT (IN)  :: a,b,c,r
    REAL(KIND=CGREAL), DIMENSION(0:n2), INTENT (OUT) :: u

    REAL(KIND=CGREAL)                             :: binv
    REAL(KIND=CGREAL), DIMENSION(0:n2)            :: gam
    INTEGER                          :: j

    binv  = oned/b(n1)
    u(n1) = r(n1)*binv
    DO j=n1+1,n2
      gam(j)  = c(j-1)*binv
      binv    = oned/(b(j)-a(j)*gam(j))
      u(j)    = (r(j)-a(j)*u(j-1))*binv
    ENDDO
    DO j=n2-1,n1,-1
      u(j) = u(j) - gam(j+1)*u(j+1)
    ENDDO

   END SUBROUTINE tdma

!******************************************************************
!...Subroutine to solve Ax=b where A is a cyclic tridiagonal matrix
!    [b c   Be]
!    [a b c   ]
!    [  a b c ]
!    [Al  a b ]
!******************************************************************
    SUBROUTINE cyclic_tdma(a,b,c,r,u,alp,bet,n1,n2)

      USE global_parameters    
      IMPLICIT NONE

      INTEGER                           , INTENT (IN)  :: n1,n2
      REAL(KIND=CGREAL), DIMENSION(0:n2), INTENT (IN)  :: a,b,c,r
      REAL(KIND=CGREAL)                 , INTENT (IN)  :: alp,bet
      REAL(KIND=CGREAL), DIMENSION(0:n2), INTENT (OUT) :: u

      REAL(KIND=CGREAL)                                :: gam,fact
      REAL(KIND=CGREAL), DIMENSION(0:n2)               :: x,y,z,bb
      INTEGER                                          :: j

      gam      = -b(n1)
      bb(n1)   =  b(n1)-gam
      bb(n2)   =  b(n2)-alp*bet/gam
 
      DO J=n1+1,n2-1
        bb(j) = b(j)
      ENDDO
 
      CALL tdma(a,bb,c,r,x,n1,n2) 

      y(n1) = gam
      y(n2) = alp
 
      DO j=n1+1,n2-1
        y(j) = zero
      ENDDO
 
      CALL tdma(a,bb,c,y,z,n1,n2) 

      fact =         ( x(n1) + bet*x(n2)/gam ) &
               /( oned + z(n1) + bet*z(n2)/gam )
 
      DO j=n1,n2
         x(j) = x(j) - fact*z(j)
         u(j) = x(j)
      ENDDO
 
   END SUBROUTINE cyclic_tdma
     
