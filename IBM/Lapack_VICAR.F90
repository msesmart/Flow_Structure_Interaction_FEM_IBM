      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!         
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992  
!         
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N 
!     ..  
!     .. Array Arguments .. 
      INTEGER            IPIV( * )  
      DOUBLE PRECISION   A( LDA, * )
!     ..
!      
!  Purpose 
!  ======= 
!      
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the 
!          matrix was interchanged with row IPIV(i). 
!   
!  INFO    (output) INTEGER
!          = 0: successful exit 
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!           
!  =====================================================================
!   
!     .. Parameters .. 
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 ) 
!     ..    
!     .. Local Scalars .. 
      INTEGER            J, JP  
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX 
      EXTERNAL           IDAMAX
!     ..    
!     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      if ( M < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!      
      if ( M == 0 .OR. N.EQ.0 ) & 
         RETURN
!      
      DO 10 J = 1, MIN( M, N ) 
!      
!        Find pivot and test for singularity.
!
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         if ( A( JP, J ).NE.ZERO ) then
!      
!           Apply the interchange to columns 1:N.
!      
            if ( JP.NE.J ) &
               CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!         
!           Compute elements J+1:M of J-th column.
!      
            if ( J < M ) & 
               CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!      
         else if ( INFO == 0 ) then
!      
            INFO = J
         end if
!
         if ( J < MIN( M, N ) ) then
!
!           Update trailing submatrix.
!
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, &
                       A( J+1, J+1 ), LDA )
         end if
   10 CONTINUE
      RETURN
!
!     End of DGETF2
!
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..  
!     .. Array Arguments ..
      INTEGER            IPIV( * ) 
      DOUBLE PRECISION   A( LDA, * )
!     ..     
!
!  Purpose   
!  =======              
!         
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!   
!  This is the right-looking Level 3 BLAS version of the algorithm.
!      
!  Arguments 
!  =========
!      
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!           
!  =====================================================================
!   
!     .. Parameters .. 
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 ) 
!     ..    
!     .. Local Scalars .. 
      INTEGER            I, IINFO, J, JB, NB 
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA 
!     ..
!     .. External Functions .. 
      INTEGER            ILAENV 
      EXTERNAL           ILAENV 
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..    
!     .. Executable Statements .. 
!                 
!     Test the input parameters. 
!
      INFO = 0
      if ( M < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      if ( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) then
!
!        Use unblocked code.
!      
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE 
!      
!        Use blocked code.
!      
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!      
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!      
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!      
!           Adjust INFO and the pivot indices.
!      
            if ( INFO == 0 .AND. IINFO.GT.0 ) &
               INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 ) 
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE 
!      
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            if ( J+JB.LE.N ) then
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                            IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                           N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
                           LDA )
               if ( J+JB.LE.M ) then
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                              N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                              A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                              LDA )
               end if
            end if
   20    CONTINUE 
      end if
      RETURN 
!
!     End of DGETRF 
!
      END       
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) -- 
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University 
!     June 30, 1999         
!                           
!     .. Scalar Arguments .. 
      INTEGER            INFO, LDA, LWORK, N
!     ..           
!     .. Array Arguments ..
      INTEGER            IPIV( * ) 
      DOUBLE PRECISION   A( LDA, * ), WORK( * ) 
!     ..                       
!                              
!  Purpose      
!  =======
!
!  DGETRI computes the inverse of a matrix using the LU factorization
!  computed by DGETRF.
!
!  This method inverts U and then computes inv(A) by solving the system
!  inv(A)*L = inv(U) for inv(A).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the factors L and U from the factorization
!          A = P*L*U as computed by DGETRF.
!          On exit, if INFO = 0, the inverse of the original matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER 
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimal performance LWORK >= N*NB, where NB is
!          the optimal blocksize returned by ILAENV.
!   
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!   
!  INFO    (output) INTEGER 
!          = 0:  successful exit  
!          < 0:  if INFO = -i, the i-th argument had an illegal value 
!          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!                singular and its inverse could not be computed.
!           
!  =====================================================================
!   
!     .. Parameters .. 
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB, &
                         NBMIN, NN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB 
      WORK( 1 ) = LWKOPT  
      LQUERY = ( LWORK == -1 )
      if ( N < 0 ) then 
         INFO = -1        
      else if ( LDA < MAX( 1, N ) ) then 
         INFO = -3        
      else if ( LWORK < MAX( 1, N ) .AND. .NOT.LQUERY ) then
         INFO = -6 
      end if 
      if ( INFO.NE.0 ) then 
         CALL XERBLA( 'DGETRI', -INFO )
         RETURN 
      else if ( LQUERY ) then 
         RETURN
      end if 
!      
!     Quick return if possible
!      
      if ( N == 0 ) &
         RETURN 
!
!     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
!     and the inverse is not computed. 
!
      CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      if ( INFO.GT.0 ) &
         RETURN
!
      NBMIN = 2
      LDWORK = N
      if ( NB.GT.1 .AND. NB < N ) then
         IWS = MAX( LDWORK*NB, 1 )
         if ( LWORK < IWS ) then
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
         end if
      ELSE
         IWS = N
      end if
!
!     Solve the equation inv(A)*L = inv(U) for inv(A).
!
      if ( NB < NBMIN .OR. NB.GE.N ) then
!
!        Use unblocked code.
!
         DO 20 J = N, 1, -1
!
!           Copy current column of L to WORK and replace with zeros.
!      
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
!      
!           Compute current column of inv(A).
!         
            if ( J < N ) & 
               CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ), &
                           LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE 
!      
!        Use blocked code.
!      
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB 
            JB = MIN( NB, N-J+1 )
!         
!           Copy current block column of L to WORK and replace with
!           zeros.  
!
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
!
!           Compute current block column of inv(A).
!
            if ( J+JB.LE.N ) &
               CALL DGEMM( 'No transpose', 'No transpose', N, JB, &
                           N-J-JB+1, -ONE, A( 1, J+JB ), LDA, &
                           WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB, &
                        ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      end if
!
!     Apply column interchanges.
!
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         if ( JP.NE.J ) &
            CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE  
!               
      WORK( 1 ) = IWS 
      RETURN       
!   
!     End of DGETRI 
!
      END    

      subroutine  dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1 
      do 50 i = mp1,n,3  
        dtemp = dx(i)
        dx(i) = dy(i) 
        dy(i) = dtemp 
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue 
      return
      end
      SUBROUTINE DSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, &
                         BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     .. 
!
!  Purpose
!  =======
!
!  DSYMM  performs one of the matrix-matrix operations
!
!     C := alpha*A*B + beta*C,
!
!  or
!
!     C := alpha*B*A + beta*C,
!
!  where alpha and beta are scalars,  A is a symmetric matrix and  B and
!  C are  m by n matrices.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE  specifies whether  the  symmetric matrix  A
!           appears on the  left or right  in the  operation as follows:
!
!              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
!
!              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
!   
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1. 
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of  the  symmetric  matrix   A  is  to  be
!           referenced as follows:
!   
!              UPLO = 'U' or 'u'   Only the upper triangular part of the
!                                  symmetric matrix is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of the
!                                  symmetric matrix is to be referenced.
!
!           Unchanged on exit.
!   
!  M      - INTEGER.
!           On entry,  M  specifies the number of rows of the matrix  C.
!           M  must be at least zero. 
!           Unchanged on exit. 
!
!  N      - INTEGER. 
!           On entry, N specifies the number of columns of the matrix C.
!           N  must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
!           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
!           the array  A  must contain the  symmetric matrix,  such that
!           when  UPLO = 'U' or 'u', the leading m by m upper triangular
!           part of the array  A  must contain the upper triangular part
!           of the  symmetric matrix and the  strictly  lower triangular
!           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!           the leading  m by m  lower triangular part  of the  array  A
!           must  contain  the  lower triangular part  of the  symmetric
!           matrix and the  strictly upper triangular part of  A  is not
!           referenced.
!           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
!           the array  A  must contain the  symmetric matrix,  such that
!           when  UPLO = 'U' or 'u', the leading n by n upper triangular
!           part of the array  A  must contain the upper triangular part
!           of the  symmetric matrix and the  strictly  lower triangular
!           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!           the leading  n by n  lower triangular part  of the  array  A
!           must  contain  the  lower triangular part  of the  symmetric
!           matrix and the  strictly upper triangular part of  A  is not
!           referenced. 
!           Unchanged on exit.
!
!  LDA    - INTEGER. 
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, n ).  
!           Unchanged on exit.  
!            
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ). 
!           Before entry, the leading  m by n part of the array  B  must
!           contain the matrix B. 
!           Unchanged on exit. 
!            
!  LDB    - INTEGER. 
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ). 
!           Unchanged on exit. 
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n updated
!           matrix.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd. 
!            
!            
!     .. External Functions ..
      LOGICAL            LSAME 
      EXTERNAL           LSAME 
!     .. External Subroutines .. 
      EXTERNAL           XERBLA 
!     .. Intrinsic Functions ..  
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, K, NROWA  
      DOUBLE PRECISION   TEMP1, TEMP2 
!     .. Parameters .. 
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Set NROWA as the number of rows of A.
!      
      if ( LSAME( SIDE, 'L' ) )THEN
         NROWA = M
      ELSE
         NROWA = N
      end if
      UPPER = LSAME( UPLO, 'U' )
!
!     Test the input parameters.
!
      INFO = 0
      if (      ( .NOT.LSAME( SIDE, 'L' ) ).AND. &
               ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
         INFO = 1
      else if ( ( .NOT.UPPER              ).AND. &
               ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
         INFO = 2
      else if ( M   < 0               )THEN
         INFO = 3
      else if ( N   < 0               )THEN
         INFO = 4
      else if ( LDA < MAX( 1, NROWA ) )THEN
         INFO = 7
      else if ( LDB < MAX( 1, M     ) )THEN
         INFO = 9
      else if ( LDC < MAX( 1, M     ) )THEN
         INFO = 12
      end if
      if ( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYMM ', INFO )
         RETURN 
      end if
!      
!     Quick return if possible.
!      
      if ( ( M == 0 ).OR.( N.EQ.0 ).OR. &  
          ( ( ALPHA == ZERO ).AND.( BETA.EQ.ONE ) ) ) &
         RETURN  
!      
!     And when  alpha.eq.zero. 
!         
      if ( ALPHA == ZERO )THEN         
         if ( BETA == ZERO )THEN
            DO 20, J = 1, N            
               DO 10, I = 1, M
                  C( I, J ) = ZERO 
   10          CONTINUE
   20       CONTINUE 
         ELSE  
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         end if
         RETURN
      end if
!
!     Start the operations.
!
      if ( LSAME( SIDE, 'L' ) )THEN
!
!        Form  C := alpha*A*B + beta*C.
!
         if ( UPPER )THEN
            DO 70, J = 1, N
               DO 60, I = 1, M
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 50, K = 1, I - 1
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   50             CONTINUE
                  if ( BETA == ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE 
                     C( I, J ) = BETA *C( I, J ) + &
                                 TEMP1*A( I, I ) + ALPHA*TEMP2
                  end if
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 100, J = 1, N
               DO 90, I = M, 1, -1
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 80, K = I + 1, M 
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   80             CONTINUE 
                  if ( BETA == ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE 
                     C( I, J ) = BETA *C( I, J ) + &
                                 TEMP1*A( I, I ) + ALPHA*TEMP2  
                  end if 
   90          CONTINUE 
  100       CONTINUE
         end if
      ELSE
!
!        Form  C := alpha*B*A + beta*C.
!
         DO 170, J = 1, N
            TEMP1 = ALPHA*A( J, J )
            if ( BETA == ZERO )THEN
               DO 110, I = 1, M
                  C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
            ELSE
               DO 120, I = 1, M
                  C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
            end if
            DO 140, K = 1, J - 1
               if ( UPPER )THEN
                  TEMP1 = ALPHA*A( K, J )
               ELSE
                  TEMP1 = ALPHA*A( J, K )
               end if
               DO 130, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
            DO 160, K = J + 1, N
               if ( UPPER )THEN  
                  TEMP1 = ALPHA*A( J, K )
               ELSE 
                  TEMP1 = ALPHA*A( K, J )
               end if  
               DO 150, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE  
      end if       
!  
      RETURN 
!            
!     End of DSYMM . 
!                  
      END       

   SUBROUTINE DGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
                         &BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1. 
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'. 
!   
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!   
!           Unchanged on exit.
!
!  M      - INTEGER. 
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!               
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero. 
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!   
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. 
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B. 
!           Unchanged on exit. 
!            
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ). 
!           Unchanged on exit. 
!            
!  BETA   - DOUBLE PRECISION. 
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit. 
!            
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). 
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA 
!     .. Intrinsic Functions ..  
      INTRINSIC          MAX 
!     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NROWA, NROWB 
!     INTEGER            NCOLA 
      DOUBLE PRECISION   TEMP
!     .. Parameters .. 
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..  
!     .. Executable Statements ..
!   
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!      
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      if ( NOTA )THEN 
         NROWA = M        
!        NCOLA = K        
      ELSE 
         NROWA = K
!        NCOLA = M
      end if
      if ( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      end if
!
!     Test the input parameters.
!
      INFO = 0
      if (      ( .NOT.NOTA                 ).AND. &
               ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      else if ( ( .NOT.NOTB                 ).AND. &
               ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
               ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      else if ( M   < 0               )THEN
         INFO = 3
      else if ( N   < 0               )THEN
         INFO = 4
      else if ( K   < 0               )THEN
         INFO = 5 
      else if ( LDA < MAX( 1, NROWA ) )THEN
         INFO = 8 
      else if ( LDB < MAX( 1, NROWB ) )THEN
         INFO = 10
      else if ( LDC < MAX( 1, M     ) )THEN
         INFO = 13
      end if
      if ( INFO.NE.0 )THEN 
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      end if     
!               
!     Quick return if possible. 
!         
      if ( ( M == 0 ).OR.( N.EQ.0 ).OR. &    
          ( ( ( ALPHA == ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )&
         RETURN 
!         
!     And if  alpha.eq.zero.           
!         
      if ( ALPHA == ZERO )THEN         
         if ( BETA == ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         end if
         RETURN
      end if
!
!     Start the operations.
!
      if ( NOTB )THEN
         if ( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               if ( BETA == ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE 
               else if ( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE 
               end if  
               DO 80, L = 1, K 
                  if ( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  end if
   80          CONTINUE 
   90       CONTINUE
         ELSE 
!         
!           Form  C := alpha*A'*B + beta*C
!            
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO 
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  if ( BETA == ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  end if
  110          CONTINUE
  120       CONTINUE
         end if
      ELSE
         if ( NOTA )THEN
!
!  Form  C := alpha*A*B' + beta*C
!
            DO J = 1, N
               if ( BETA == ZERO )THEN
                  DO I = 1, M
                     C( I, J ) = ZERO
                  end do
               else if ( BETA.NE.ONE )THEN
                  DO I = 1, M
                     C( I, J ) = BETA*C( I, J )
                  end do 
               end if 
               DO L = 1, K
                  if ( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
                     end do
                  end if
               end do
            end do
         ELSE
!         
!  Form  C := alpha*A'*B' + beta*C
!   
            DO J = 1, N
               DO I = 1, M
                  TEMP = ZERO 
                  DO L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
                  end do
                  if ( BETA == ZERO )THEN 
                     C( I, J ) = ALPHA*TEMP
                  ELSE 
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  end if
               end do
            end do
         end if
      end if
 
      RETURN
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT(' ** On entry to ',A6,' parameter number ', I2, ' had',&
            &'an illegal value' )
!
!     End of XERBLA
!
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!   
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!            
!  INCY   - INTEGER. 
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!            
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!            
!  LDA    - INTEGER. 
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ). 
!           Unchanged on exit. 
!            
!            
!  Level 2 Blas routine.
!   
!  -- Written on 22-October-1986. 
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M < 0 )THEN
         INFO = 1
      else if ( N < 0 )THEN
         INFO = 2
      else if ( INCX == 0 )THEN
         INFO = 5
      else if ( INCY == 0 )THEN
         INFO = 7
      else if ( LDA < MAX( 1, M ) )THEN
         INFO = 9 
      end if 
      if ( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN 
      end if  
!      
!     Quick return if possible.
!      
      if ( ( M == 0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
         RETURN 
!      
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!      
      if ( INCY.GT.0 )THEN 
         JY = 1
      ELSE  
         JY = 1 - ( N - 1 )*INCY
      end if  
      if ( INCX == 1 )THEN
         DO 20, J = 1, N
            if ( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M 
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE  
            end if
            JY = JY + INCY  
   20    CONTINUE
      ELSE
         if ( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         end if
         DO 40, J = 1, N
            if ( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            end if
            JY = JY + INCY
   40    CONTINUE
      end if
!
      RETURN
!
!     End of DGER  .
!
      END

      LOGICAL          FUNCTION LSAME( CA, CB )
!   
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     .. 
!
!  Purpose
!  ======= 
!      
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.  
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA == CB
      if ( LSAME ) &
         RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA ) 
      INTB = ICHAR( CB ) 
!           
      if ( ZCODE == 90 .OR. ZCODE.EQ.122 ) then
!  
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'. 
!      
         if ( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         if ( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!      
      else if ( ZCODE == 233 .OR. ZCODE.EQ.169 ) then
!      
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'. 
!
         if ( INTA.GE.129 .AND. INTA.LE.137 .OR. &
             INTA.GE.145 .AND. INTA.LE.153 .OR. &
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         if ( INTB.GE.129 .AND. INTB.LE.137 .OR. &
             INTB.GE.145 .AND. INTB.LE.153 .OR. &  
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!      
      else if ( ZCODE == 218 .OR. ZCODE.EQ.250 ) then
!      
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         if ( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         if ( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      end if
      LSAME = INTA == INTB
!
!     RETURN
!
!     End of LSAME
!
      END
      double precision function dasum(n,dx,incx)
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
!
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!         
      nincx = n*incx  
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i)) 
   10 continue 
      dasum = dtemp
      return 
!
!        code for increment equal to 1
!
!      
!        clean-up loop
!      
   20 m = mod(n,6) 
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m  
        dtemp = dtemp + dabs(dx(i)) 
   30 continue  
      if( n .lt. 6 ) go to 60 
   40 mp1 = m + 1
      do 50 i = mp1,n,6 
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))&
        + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue 
   60 dasum = dtemp
      return 
      end 
      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRTRI computes the inverse of a real upper or lower triangular
!  matrix A.
!
!  This is the Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  A is upper triangular;
!          = 'L':  A is lower triangular.
!
!  DIAG    (input) CHARACTER*1
!          = 'N':  A is non-unit triangular;
!          = 'U':  A is unit triangular.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of the array A contains
!          the upper triangular matrix, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of the array A contains
!          the lower triangular matrix, and the strictly upper
!          triangular part of A is not referenced.  If DIAG = 'U', the
!          diagonal elements of A are also not referenced and are
!          assumed to be 1.
!          On exit, the (triangular) inverse of the original matrix, in
!          the same storage format.
!
!  LDA     (input) INTEGER 
!          The leading dimension of the array A.  LDA >= max(1,N).
!   
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!               matrix is singular and its inverse can not be computed.
!
!  =====================================================================
!           
!     .. Parameters .. 
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..    
!     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER 
      INTEGER            J, JB, NB, NN 
!     ..    
!     .. External Functions .. 
      LOGICAL            LSAME 
      INTEGER            ILAENV  
      EXTERNAL           LSAME, ILAENV 
!     ..    
!     .. External Subroutines ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -5
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DTRTRI', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) & 
         RETURN           
!      
!     Check for singularity if non-unit.
!      
      if ( NOUNIT ) then
         DO 10 INFO = 1, N 
            if ( A( INFO, INFO ) == ZERO ) &
               RETURN 
   10    CONTINUE
         INFO = 0
      end if 
!      
!     Determine the block size for this environment. 
!         
      NB = ILAENV( 1, 'DTRTRI', UPLO,  DIAG, N, -1, -1, -1 ) 
!      NB = ILAENV( 1, 'DTRTRI', UPLO//DIAG, N, -1, -1, -1 ) !Wanh changed "UPLO,DIAG" to "UPLO//DIAG"

      if ( NB.LE.1 .OR. NB.GE.N ) then
!      
!        Use unblocked code
!      
         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE 
!      
!        Use blocked code 
!         
         if ( UPPER ) then
!
!           Compute inverse of upper triangular matrix
!
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
!
!              Compute rows 1:j-1 of current block column
!
               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1, &
                           JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG,J-1,&
                           JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )

!
!              Compute inverse of current diagonal block
!
               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
!
!           Compute inverse of lower triangular matrix
!
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               if ( J+JB.LE.N ) then
!
!                 Compute rows j+jb:n of current block column
!
                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG, &
                              N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, &
                              A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG, &
                              N-J-JB+1, JB, -ONE, A( J, J ), LDA, &
                              A( J+JB, J ), LDA )
               end if
!               
!              Compute inverse of current diagonal block
!               
               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE 
         end if             
      end if
!               
      RETURN
!               
!     End of DTRTRI 
!         
      END

      SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&
                         B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRMM  performs one of the matrix-matrix operations
!
!     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!            
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!            
!           Unchanged on exit. 
!
!  M      - INTEGER. 
!           On entry, M specifies the number of rows of B. M must be at
!           least zero. 
!           Unchanged on exit.
!            
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero. 
!           Unchanged on exit. 
!
!  ALPHA  - DOUBLE PRECISION. 
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.  
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!            
!  -- Written on 8-February-1989. 
!     Jack Dongarra, Argonne National Laboratory.  
!     Iain Duff, AERE Harwell.  
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd. 
!            
!            
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME 
!     .. External Subroutines ..  
      EXTERNAL           XERBLA 
!     .. Intrinsic Functions ..  
      INTRINSIC          MAX 
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER 
      INTEGER            I, INFO, J, K, NROWA 
      DOUBLE PRECISION   TEMP 
!     .. Parameters .. 
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..     
!     .. Executable Statements .. 
!            
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      if ( LSIDE )THEN 
         NROWA = M
      ELSE
         NROWA = N
      end if
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      if (      ( .NOT.LSIDE                ).AND. &
               ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      else if ( ( .NOT.UPPER                ).AND. &
               ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      else if ( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      else if ( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
               ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      else if ( M   < 0               )THEN
         INFO = 5
      else if ( N   < 0               )THEN
         INFO = 6
      else if ( LDA < MAX( 1, NROWA ) )THEN
         INFO = 9
      else if ( LDB < MAX( 1, M     ) )THEN
         INFO = 11
      end if 
      if ( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      end if  
!      
!     Quick return if possible.
!      
      if ( N == 0 ) & 
         RETURN 
!         
!     And when  alpha.eq.zero.               
!               
      if ( ALPHA == ZERO )THEN
         DO 20, J = 1, N 
            DO 10, I = 1, M 
               B( I, J ) = ZERO 
   10       CONTINUE
   20    CONTINUE  
         RETURN 
      end if 
!      
!     Start the operations.
!      
      if ( LSIDE )THEN
         if ( LSAME( TRANSA, 'N' ) )THEN 
!         
!           Form  B := alpha*A*B.    
!         
            if ( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     if ( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        if ( NOUNIT ) &
                           TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     end if
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     if ( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        if ( NOUNIT ) &
                           B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     end if
   70             CONTINUE
   80          CONTINUE
            end if
         ELSE 
!               
!           Form  B := alpha*A'*B.
!                     
            if ( UPPER )THEN  
               DO 110, J = 1, N 
                  DO 100, I = M, 1, -1  
                     TEMP = B( I, J )
                     if ( NOUNIT ) & 
                        TEMP = TEMP*A( I, I ) 
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE  
            ELSE   
               DO 140, J = 1, N  
                  DO 130, I = 1, M 
                     TEMP = B( I, J ) 
                     if ( NOUNIT ) & 
                        TEMP = TEMP*A( I, I ) 
                     DO 120, K = I + 1, M  
                        TEMP = TEMP + A( K, I )*B( K, J ) 
  120                CONTINUE 
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            end if
         end if
      ELSE
         if ( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*A.
!
            if ( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  if ( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     if ( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     end if
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  if ( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N 
                     if ( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE  
                     end if 
  210             CONTINUE 
  220          CONTINUE 
            end if 
         ELSE         
!  
!           Form  B := alpha*B*A'.  
!                     
            if ( UPPER )THEN  
               DO 260, K = 1, N  
                  DO 240, J = 1, K - 1  
                     if ( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     end if 
  240             CONTINUE 
                  TEMP = ALPHA 
                  if ( NOUNIT ) & 
                     TEMP = TEMP*A( K, K )
                  if ( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  end if
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     if ( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     end if
  280             CONTINUE
                  TEMP = ALPHA
                  if ( NOUNIT ) &
                     TEMP = TEMP*A( K, K )
                  if ( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  end if
  300          CONTINUE
            end if
         end if
      end if
!
      RETURN       
!                     
!     End of DTRMM .     
!  
      END          
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, &
                       N4 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case. 
!           
!  OPTS    (input) CHARACTER*(*) 
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.  
!           
!  N1      (input) INTEGER 
!  N2      (input) INTEGER  
!  N3      (input) INTEGER  
!  N4      (input) INTEGER 
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.  
!           
! (ILAENV) (output) INTEGER 
!          >= 0: the value of the parameter specified by ISPEC 
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!           
!  Further Details 
!  =============== 
!           
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines: 
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      if ( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
!     ..
!     .. Executable Statements ..  
!   
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000, &
              1100 ) ISPEC 
!       
!     Invalid value for ISPEC 
!       
      ILAENV = -1 
      RETURN
!       
  100 CONTINUE 
!
!     Convert NAME to upper case if the first character is lower case. 
!
      ILAENV = 1 
      SUBNAM = NAME       
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )   
      if ( IZ == 90 .OR. IZ.EQ.122 ) then
!      
!        ASCII character set 
!      
         if ( IC.GE.97 .AND. IC.LE.122 ) then
            SUBNAM( 1:1 ) = CHAR( IC-32 ) 
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if ( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         end if
!
      else if ( IZ == 233 .OR. IZ.EQ.169 ) then
!
!        EBCDIC character set
!
         if ( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) then
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if ( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) ) &
                  SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         end if
!
      else if ( IZ == 218 .OR. IZ.EQ.250 ) then
!
!        Prime machines:  ASCII+128
!
         if ( IC.GE.225 .AND. IC.LE.250 ) then
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               if ( IC.GE.225 .AND. IC.LE.250 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         end if
      end if 
!
      C1 = SUBNAM( 1:1 ) 
      SNAME = C1 == 'S' .OR. C1.EQ.'D'
      CNAME = C1 == 'C' .OR. C1.EQ.'Z' 
      if ( .NOT.( CNAME .OR. SNAME ) ) &  
         RETURN 
      C2 = SUBNAM( 2:3 ) 
      C3 = SUBNAM( 4:6 )  
      C4 = C3( 2:3 ) 
!               
      GO TO ( 110, 200, 300 ) ISPEC 
!                   
  110 CONTINUE     
!   
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision. 
!
      NB = 1 
!            
      if ( C2 == 'GE' ) then
         if ( C3 == 'TRF' ) then 
            if ( SNAME ) then  
               NB = 64 
            ELSE
               NB = 64
            end if
         else if ( C3 == 'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3 == 'QLF' ) then
            if ( SNAME ) then
               NB = 32
            ELSE
               NB = 32
            end if
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then
               NB = 32
            ELSE
               NB = 32
            end if
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NB = 32
            ELSE
               NB = 32
            end if
         else if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         end if
      else if ( C2 == 'PO' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NB = 64 
            ELSE   
               NB = 64  
            end if 
         end if 
      else if ( C2 == 'SY' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then 
               NB = 64  
            ELSE 
               NB = 64
            end if 
         else if ( SNAME .AND. C3 == 'TRD' ) then
            NB = 32 
         else if ( SNAME .AND. C3 == 'GST' ) then
            NB = 64  
         end if 
      else if ( CNAME .AND. C2 == 'HE' ) then
         if ( C3 == 'TRF' ) then
            NB = 64 
         else if ( C3 == 'TRD' ) then
            NB = 32  
         else if ( C3 == 'GST' ) then
            NB = 64  
         end if  
      else if ( SNAME .AND. C2 == 'OR' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NB = 32
            end if
         else if ( C3( 1:1 ) == 'M' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NB = 32
            end if
         end if
      else if ( CNAME .AND. C2 == 'UN' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NB = 32
            end if
         else if ( C3( 1:1 ) == 'M' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NB = 32
            end if
         end if
      else if ( C2 == 'GB' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then 
               if ( N4.LE.64 ) then 
                  NB = 1 
               ELSE 
                  NB = 32
               end if
            ELSE  
               if ( N4.LE.64 ) then 
                  NB = 1 
               ELSE 
                  NB = 32
               end if
            end if
         end if  
      else if ( C2 == 'PB' ) then 
         if ( C3 == 'TRF' ) then  
            if ( SNAME ) then 
               if ( N2.LE.64 ) then
                  NB = 1
               ELSE
                  NB = 32 
               end if 
            ELSE 
               if ( N2.LE.64 ) then
                  NB = 1
               ELSE
                  NB = 32
               end if  
            end if 
         end if
      else if ( C2 == 'TR' ) then
         if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         end if
      else if ( C2 == 'LA' ) then
         if ( C3 == 'UUM' ) then
            if ( SNAME ) then
               NB = 64
            ELSE
               NB = 64
            end if
         end if
      else if ( SNAME .AND. C2 == 'ST' ) then
         if ( C3 == 'EBZ' ) then
            NB = 1
         end if
      end if
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      if ( C2 == 'GE' ) then  
         if ( C3 == 'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3 == 'QLF' ) then
            if ( SNAME ) then
               NBMIN = 2
            ELSE 
               NBMIN = 2
            end if
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then 
               NBMIN = 2  
            ELSE 
               NBMIN = 2
            end if 
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NBMIN = 2 
            ELSE  
               NBMIN = 2
            end if
         else if ( C3 == 'TRI' ) then
            if ( SNAME ) then
               NBMIN = 2
            ELSE
               NBMIN = 2
            end if
         end if 
      else if ( C2 == 'SY' ) then
         if ( C3 == 'TRF' ) then
            if ( SNAME ) then
               NBMIN = 8
            ELSE
               NBMIN = 8
            end if
         else if ( SNAME .AND. C3 == 'TRD' ) then
            NBMIN = 2
         end if
      else if ( CNAME .AND. C2 == 'HE' ) then
         if ( C3 == 'TRD' ) then
            NBMIN = 2
         end if
      else if ( SNAME .AND. C2 == 'OR' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NBMIN = 2
            end if
         else if ( C3( 1:1 ) == 'M' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NBMIN = 2
            end if
         end if
      else if ( CNAME .AND. C2 == 'UN' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NBMIN = 2
            end if 
         else if ( C3( 1:1 ) == 'M' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NBMIN = 2 
            end if 
         end if 
      end if  
      ILAENV = NBMIN 
      RETURN 
!            
  300 CONTINUE   
!                
!     ISPEC = 3:  crossover point
!            
      NX = 0 
      if ( C2 == 'GE' ) then 
         if ( C3 == 'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3 == 'QLF' ) then 
            if ( SNAME ) then
               NX = 128
            ELSE
               NX = 128 
            end if 
         else if ( C3 == 'HRD' ) then
            if ( SNAME ) then
               NX = 128
            ELSE
               NX = 128
            end if
         else if ( C3 == 'BRD' ) then
            if ( SNAME ) then
               NX = 128
            ELSE
               NX = 128
            end if
         end if
      else if ( C2 == 'SY' ) then
         if ( SNAME .AND. C3 == 'TRD' ) then
            NX = 32
         end if
      else if ( CNAME .AND. C2 == 'HE' ) then
         if ( C3 == 'TRD' ) then
            NX = 32
         end if
      else if ( SNAME .AND. C2 == 'OR' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NX = 128
            end if
         end if
      else if ( CNAME .AND. C2 == 'UN' ) then
         if ( C3( 1:1 ) == 'G' ) then
            if ( C4 == 'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4 == 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4 == 'BR' ) then
               NX = 128
            end if  
         end if 
      end if    
      ILAENV = NX
      RETURN    
!            
  400 CONTINUE 
!      
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!            
      ILAENV = 6
      RETURN 
!         
  500 CONTINUE  
!         
!     ISPEC = 5:  minimum column dimension (not used)
!         
      ILAENV = 2  
      RETURN     
!                
  600 CONTINUE  
!            
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  900 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
 1000 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!      
!     ILAENV = 0
      ILAENV = 1
      if ( ILAENV == 1 ) then
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      end if 
      RETURN
!      
 1100 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0   
      ILAENV = 1
      if ( ILAENV == 1 ) then
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      end if
      RETURN 
!
!     End of ILAENV 
!                  
      END          
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!      
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!      
!  RETURN VALUE:  INTEGER 
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                         NEGZRO, NEWZRO, POSINF
!     .. 
!     .. Executable Statements .. 
      IEEECK = 1
!   
      POSINF = ONE / ZERO
      if ( POSINF.LE.ONE ) then
         IEEECK = 0 
         RETURN 
      end if 
!           
      NEGINF = -ONE / ZERO 
      if ( NEGINF.GE.ZERO ) then
         IEEECK = 0 
         RETURN  
      end if 
!           
      NEGZRO = ONE / ( NEGINF+ONE )
      if ( NEGZRO.NE.ZERO ) then
         IEEECK = 0 
         RETURN  
      end if
!
      NEGINF = ONE / NEGZRO
      if ( NEGINF.GE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      NEWZRO = NEGZRO + ZERO
      if ( NEWZRO.NE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      POSINF = ONE / NEWZRO
      if ( POSINF.LE.ONE ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGINF = NEGINF*POSINF
      if ( NEGINF.GE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      POSINF = POSINF*POSINF
      if ( POSINF.LE.ONE ) then
         IEEECK = 0
         RETURN
      end if
!      
!      
!         
!         
!     Return if we were only asked to check infinity arithmetic
!
      if ( ISPEC == 0 ) & 
         RETURN 
!         
      NAN1 = POSINF + NEGINF
!      
      NAN2 = POSINF / NEGINF
!      
      NAN3 = POSINF / POSINF 
!         
      NAN4 = POSINF*ZERO
!      
      NAN5 = NEGINF*NEGZRO
!      
      NAN6 = NAN5*0.0 
!         
      if ( NAN1 == NAN1 ) then
         IEEECK = 0
         RETURN
      end if  
!      
      if ( NAN2 == NAN2 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN3 == NAN3 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN4 == NAN4 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN5 == NAN5 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN6 == NAN6 ) then
         IEEECK = 0
         RETURN
      end if
!
      RETURN
      END

      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!      
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!      
!     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRTI2 computes the inverse of a real upper or lower triangular
!  matrix.
!
!  This is the Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular matrix A.  If UPLO = 'U', the
!          leading n by n upper triangular part of the array A contains
!          the upper triangular matrix, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of the array A contains
!          the lower triangular matrix, and the strictly upper
!          triangular part of A is not referenced.  If DIAG = 'U', the
!          diagonal elements of A are also not referenced and are 
!          assumed to be 1.
!
!          On exit, the (triangular) inverse of the original matrix, in
!          the same storage format.
!   
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!   
!  INFO    (output) INTEGER 
!          = 0: successful exit 
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  =====================================================================
!           
!     .. Parameters .. 
      DOUBLE PRECISION   ONE 
      PARAMETER          ( ONE = 1.0D+0 )
!     ..    
!     .. Local Scalars .. 
      LOGICAL            NOUNIT, UPPER
      INTEGER            J 
      DOUBLE PRECISION   AJJ 
!     ..    
!     .. External Functions .. 
      LOGICAL            LSAME 
      EXTERNAL           LSAME 
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DTRMV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) then
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -5
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      end if
!
      if ( UPPER ) then
!      
!        Compute inverse of upper triangular matrix.
!      
         DO 10 J = 1, N
            if ( NOUNIT ) then 
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE 
               AJJ = -ONE
            end if 
!
!           Compute elements 1:j-1 of j-th column.
!      
            CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA, &
                        A( 1, J ), 1 ) 
            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE 
      ELSE 
!      
!        Compute inverse of lower triangular matrix.
!      
         DO 20 J = N, 1, -1
            if ( NOUNIT ) then
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )  
            ELSE
               AJJ = -ONE
            end if
            if ( J < N ) then
!
!              Compute elements j+1:n of j-th column.
!
               CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J, &
                           A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
            end if
   20    CONTINUE
      end if
!
      RETURN
!
!     End of DTRTI2
!
      END
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A,LDA,&
                         B, LDB )
!     .. Scalar Arguments .. 
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA 
!     .. Array Arguments .. 
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!      
!  Purpose 
!  ======= 
!      
!  DTRSM  solves one of the matrix equations
!      
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!   
!  TRANSA - CHARACTER*1. 
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!   
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!   
!           Unchanged on exit. 
!            
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!               
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!            
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!            
!           Unchanged on exit. 
!
!  M      - INTEGER. 
!           On entry, M specifies the number of rows of B. M must be at
!           least zero. 
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit. 
!            
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.  
!            
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ). 
!           Unchanged on exit. 
!            
!            
!  Level 3 Blas routine. 
!            
!            
!  -- Written on 8-February-1989.  
!     Jack Dongarra, Argonne National Laboratory. 
!     Iain Duff, AERE Harwell. 
!     Jeremy Du Croz, Numerical Algorithms Group Ltd. 
!     Sven Hammarling, Numerical Algorithms Group Ltd. 
!            
!
!     .. External Functions ..
      LOGICAL            LSAME 
      EXTERNAL           LSAME 
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      if ( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      end if
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      if (      ( .NOT.LSIDE                ).AND. &
               ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      else if ( ( .NOT.UPPER                ).AND. &
               ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2 
      else if ( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3         
      else if ( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
               ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4 
      else if ( M   < 0               )THEN  
         INFO = 5
      else if ( N   < 0               )THEN
         INFO = 6
      else if ( LDA < MAX( 1, NROWA ) )THEN
         INFO = 9
      else if ( LDB < MAX( 1, M     ) )THEN
         INFO = 11 
      end if 
      if ( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      end if  
!      
!     Quick return if possible.
!      
      if ( N == 0 ) & 
         RETURN 
!         
!     And when  alpha.eq.zero.
!
      if ( ALPHA == ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      end if
!
!     Start the operations.
!
      if ( LSIDE )THEN
         if ( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            if ( UPPER )THEN
               DO 60, J = 1, N
                  if ( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  end if
                  DO 50, K = M, 1, -1
                     if ( B( K, J ).NE.ZERO )THEN
                        if ( NOUNIT ) &
                           B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     end if
   50             CONTINUE 
   60          CONTINUE 
            ELSE 
               DO 100, J = 1, N
                  if ( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  end if
                  DO 90 K = 1, M
                     if ( B( K, J ).NE.ZERO )THEN
                        if ( NOUNIT ) &
                           B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     end if 
   90             CONTINUE 
  100          CONTINUE  
            end if    
         ELSE      
!                  
!           Form  B := alpha*inv( A' )*B. 
!                        
            if ( UPPER )THEN 
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     if ( NOUNIT ) &
                        TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     if ( NOUNIT ) &
                        TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            end if
         end if
      ELSE
         if ( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!               
            if ( UPPER )THEN 
               DO 210, J = 1, N 
                  if ( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M  
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE 
                  end if 
                  DO 190, K = 1, J - 1
                     if ( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     end if  
  190             CONTINUE 
                  if ( NOUNIT )THEN 
                     TEMP = ONE/A( J, J ) 
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE 
                  end if 
  210          CONTINUE 
            ELSE 
               DO 260, J = N, 1, -1
                  if ( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  end if 
                  DO 240, K = J + 1, N
                     if ( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     end if
  240             CONTINUE
                  if ( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  end if
  260          CONTINUE
            end if
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            if ( UPPER )THEN
               DO 310, K = N, 1, -1
                  if ( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  end if
                  DO 290, J = 1, K - 1
                     if ( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M 
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE  
                     end if 
  290             CONTINUE 
                  if ( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE 
                  end if 
  310          CONTINUE 
            ELSE   
               DO 360, K = 1, N
                  if ( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  end if 
                  DO 340, J = K + 1, N
                     if ( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K ) 
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     end if
  340             CONTINUE 
                  if ( ALPHA.NE.ONE )THEN 
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  end if
  360          CONTINUE
            end if
         end if
      end if
!
      RETURN
!
!     End of DTRSM .
!
      END
      SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!     .. Scalar Arguments .. 
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments .. 
      DOUBLE PRECISION   A( LDA, * ), X( * )
!     .. 
!      
!  Purpose
!  =======
!   
!  DTRMV  performs one of the matrix-vector operations
!   
!     x := A*x,   or   x := A'*x,
!      
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows: 
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!   
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!            
!           Unchanged on exit. 
!
!  N      - INTEGER. 
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.  
!           Unchanged on exit.
!            
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER. 
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP 
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
!     .. External Functions ..
      LOGICAL            LSAME 
      EXTERNAL           LSAME 
!     .. External Subroutines .. 
      EXTERNAL           XERBLA 
!     .. Intrinsic Functions .. 
      INTRINSIC          MAX
!     ..   
!     .. Executable Statements .. 
!            
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND. &
               .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1 
      else if ( .NOT.LSAME( TRANS, 'N' ).AND. &
               .NOT.LSAME( TRANS, 'T' ).AND. &
               .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2 
      else if ( .NOT.LSAME( DIAG , 'U' ).AND. &
               .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3 
      else if ( N < 0 )THEN 
         INFO = 4       
      else if ( LDA < MAX( 1, N ) )THEN
         INFO = 6
      else if ( INCX == 0 )THEN
         INFO = 8
      end if
      if ( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMV ', INFO )
         RETURN
      end if
!
!     Quick return if possible.
!
      if ( N == 0 ) &
         RETURN
!
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      if ( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      else if ( INCX.NE.1 )THEN
         KX = 1
      end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      if ( LSAME( TRANS, 'N' ) )THEN
!         
!        Form  x := A*x. 
!         
         if ( LSAME( UPLO, 'U' ) )THEN
            if ( INCX == 1 )THEN
               DO 20, J = 1, N 
                  if ( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     if ( NOUNIT ) &
                        X( J ) = X( J )*A( J, J )
                  end if
   20          CONTINUE 
            ELSE
               JX = KX  
               DO 40, J = 1, N  
                  if ( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX 
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     if ( NOUNIT ) & 
                        X( JX ) = X( JX )*A( J, J ) 
                  end if
                  JX = JX + INCX 
   40          CONTINUE
            end if
         ELSE
            if ( INCX == 1 )THEN
               DO 60, J = N, 1, -1
                  if ( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     if ( NOUNIT ) &
                        X( J ) = X( J )*A( J, J )
                  end if
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  if ( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     if ( NOUNIT ) &
                        X( JX ) = X( JX )*A( J, J )
                  end if
                  JX = JX - INCX
   80          CONTINUE
            end if
         end if
      ELSE   
!               
!        Form  x := A'*x.  
!                     
         if ( LSAME( UPLO, 'U' ) )THEN  
            if ( INCX == 1 )THEN  
               DO 100, J = N, 1, -1
                  TEMP = X( J ) 
                  if ( NOUNIT ) & 
                     TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE  
                  X( J ) = TEMP
  100          CONTINUE  
            ELSE   
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX ) 
                  IX   = JX 
                  if ( NOUNIT ) &  
                     TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX  
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE 
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            end if
         ELSE
            if ( INCX == 1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  if ( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  if ( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            end if 
         end if    
      end if    
!            
      RETURN 
!            
!     End of DTRMV . 
!                  
      END   

      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
     
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!      
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999 
!      
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     .. 
!     .. Array Arguments ..
      INTEGER            IPIV( * ) 
      DOUBLE PRECISION   A( LDA, * )
!     .. 
!
!  Purpose 
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!   
!  INCX    (input) INTEGER 
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!   
!  Further Details
!  ===============  
!           
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!           
! =====================================================================
!           
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP 
!     ..
!     .. Executable Statements ..
!           
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      if ( INCX.GT.0 ) then
         IX0 = K1 
         I1 = K1 
         I2 = K2
         INC = 1 
      else if ( INCX < 0 ) then 
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      end if
!
      N32 = ( N / 32 )*32
      if ( N32.NE.0 ) then
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               if ( IP.NE.I ) then
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               end if
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      end if
      if ( N32.NE.N ) then
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX ) 
            if ( IP.NE.I ) then
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            end if
            IX = IX + INCX
   50    CONTINUE 
      end if 
!            
      RETURN 
!               
!     End of DLASWP  
!                  
      END             
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
                         BETA, Y, INCY )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha. 
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!            
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit. 
!
!  X      - DOUBLE PRECISION array of DIMENSION at least 
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least  
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.  
!           Unchanged on exit.  
!            
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.  
!           Unchanged on exit.  
!            
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 ) 
!     .. Local Scalars .. 
      DOUBLE PRECISION   TEMP 
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
!     .. External Functions .. 
      LOGICAL            LSAME 
      EXTERNAL           LSAME
!     .. External Subroutines ..  
      EXTERNAL           XERBLA 
!     .. Intrinsic Functions .. 
      INTRINSIC          MAX 
!     ..
!     .. Executable Statements ..
!            
!     Test the input parameters. 
!            
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND. &
               .NOT.LSAME( TRANS, 'T' ).AND. &
               .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1 
      else if ( M < 0 )THEN 
         INFO = 2 
      else if ( N < 0 )THEN 
         INFO = 3 
      else if ( LDA < MAX( 1, M ) )THEN
         INFO = 6
      else if ( INCX == 0 )THEN
         INFO = 8 
      else if ( INCY == 0 )THEN
         INFO = 11
      end if
      if ( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      end if
!
!     Quick return if possible.
!
      if ( ( M == 0 ).OR.( N.EQ.0 ).OR. &
          ( ( ALPHA == ZERO ).AND.( BETA.EQ.ONE ) ) ) &
         RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      if ( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      end if
      if ( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      end if
      if ( INCY.GT.0 )THEN 
         KY = 1  
      ELSE 
         KY = 1 - ( LENY - 1 )*INCY
      end if 
!         
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!      
!     First form  y := beta*y.
!      
      if ( BETA.NE.ONE )THEN 
         if ( INCY == 1 )THEN
            if ( BETA == ZERO )THEN
               DO 10, I = 1, LENY 
                  Y( I ) = ZERO   
   10          CONTINUE
            ELSE 
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            end if
         ELSE  
            IY = KY
            if ( BETA == ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            end if
         end if
      end if
      if ( ALPHA == ZERO ) &
         RETURN
      if ( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         if ( INCY == 1 )THEN
            DO 60, J = 1, N
               if ( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               end if
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               if ( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY 
   70             CONTINUE
               end if
               JX = JX + INCX
   80       CONTINUE
         end if 
      ELSE 
!      
!        Form  y := alpha*A'*x + y.
!         
         JY = KY
         if ( INCX == 1 )THEN
            DO 100, J = 1, N 
               TEMP = ZERO 
               DO 90, I = 1, M 
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE 
               Y( JY ) = Y( JY ) + ALPHA*TEMP 
               JY      = JY      + INCY
  100       CONTINUE 
         ELSE   
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX 
               DO 110, I = 1, M 
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         end if
      end if
!
      RETURN
!
!     End of DGEMV .
!
      END
      integer function idamax(n,dx,incx)  
!      
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78. 
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!      
      double precision dx(*),dmax
      integer i,incx,ix,n
!   
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return 
      if(incx.eq.1)go to 20
!   
!        code for increment not equal to 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end

