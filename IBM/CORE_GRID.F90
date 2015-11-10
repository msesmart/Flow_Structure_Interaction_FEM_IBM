!---------------------------------------------
!   SUBROUTINE make_grid()
!   SUBROUTINE metrics()
!---------------------------------------------



!---------------------------------------------
! grid convention
!
!      ny+1---------------------------------------
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!       ny +==++==+==+==+==+===+==+==+==+==+==++==+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!     ^    +--++--+--+--+--+---+--+--+--+--+--++--+
!   dy|    |  ||  |  |  |  | * |  |  |  |  |  ||  |
!     -  j +--++--+--+--+--+---+--+--+--+--+--++--+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        2 +--++--+--+--+--+---+--+--+--+--+--++--+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        1 +==++==+==+==+==+===+==+==+==+==+==++==+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        0 ---------------------------------------
!          0  1   2        i                  nx  nx+1
!                          <--->
!                           dx
!
!---------------------------------------------
   SUBROUTINE make_grid()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL)    :: delta
    INTEGER :: i,junk
 
    IF (xgrid_unif == UNIFORM_GRID) then
      delta  = xout/REAL(nx-1,KIND=CGREAL)
      DO i=1,nx
        x(i) = REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=10,FILE='xgrid.dat')
      DO i=1,nx
        read(10,*)junk,x(i)
      ENDDO
      CLOSE(10)
    ENDIF
    x(0)    = -x(2)
    x(nx+1) = x(nx) + x(nx)-x(nx-1)
  
    IF (ygrid_unif == UNIFORM_GRID) then
      delta  = yout/REAL(ny-1,KIND=CGREAL)
      DO i=1,ny
        y(i) = REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=11,FILE='ygrid.dat')
      DO i=1,ny
        read(11,*)junk,y(i)
      ENDDO
      CLOSE(11)
    ENDIF
    y(0)    = -y(2)
    y(ny+1) = y(ny) + y(ny)-y(ny-1)
  
    IF (zgrid_unif == UNIFORM_GRID) then
      delta  = zout/REAL(nz-1,KIND=CGREAL)
      DO i=1,nz
        z(i) = REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=12,FILE='zgrid.dat')
      DO i=1,nz
        read(12,*)junk,z(i)
      ENDDO
      CLOSE(12)
    ENDIF
    z(0)    = -z(2)
    z(nz+1) = z(nz) + z(nz)-z(nz-1)

! write out grid file for plotting
!    OPEN(UNIT=59,FILE='grid_plot.dat')
!
!
!    WRITE(59,*)'VARIABLES="X","Y","Z"'
!    WRITE(59,*)'ZONE F=POINT, I=',nx,', J=',ny,' K=',nz
!    DO k=1,nz
!    DO j=1,ny
!    DO i=1,nx
!       WRITE(59,123)x(i),y(j),z(k)
!    ENDDO
!    ENDDO
!    ENDDO
!
123 FORMAT(3(2x,e14.7))

    dx = zero
    dy = zero
    dz = zero
    dxc= zero
    dyc= zero
    dzc= zero
 
    dxinv = zero
    dyinv = zero
    dzinv = zero
    dxcinv= zero
    dycinv= zero
    dzcinv= zero
 
    DO i=0,nx
      xc(i)    = half*(x(i+1)+x(i)) 
      dx(i)    =             x(i+1)-x(i)  
      dxinv(i) = oned/dx(i)
    ENDDO
    xc(nx+1)= x(nx+1)
   
    DO i=0,ny
      yc(i)    = half*(y(i+1)+y(i)) 
      dy(i)    =             y(i+1)-y(i)  
      dyinv(i) = oned/dy(i)
    ENDDO
    yc(ny+1)= y(ny+1)
   
    DO i=0,nz
      zc(i)    = half*(z(i+1)+z(i)) 
      dz(i)    =             z(i+1)-z(i)  
      dzinv(i) = oned/dz(i)
    ENDDO
    zc(nz+1)= z(nz+1)

    DO i=1,nx
      dxc(i)    = xc(i)-xc(i-1)  
      dxcinv(i) = oned/dxc(i)
    ENDDO
   
    DO i=1,ny
      dyc(i)    = yc(i)-yc(i-1)  
      dycinv(i) = oned/dyc(i)
    ENDDO
   
    DO i=1,nz
      dzc(i)    = zc(i)-zc(i-1)  
      dzcinv(i) = oned/dzc(i)
    ENDDO

   END SUBROUTINE make_grid  
!------------------------------------------- 
    
!-------------------------------------------     
   SUBROUTINE metrics()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k

!            |--*--|--*--|--*--|-----|
!              i-1 w  i  e i+1 
!
!        g(w) = fx(i)g(i) + (1-f(i))*g(i-1)
!
! Interpolation metrics
    DO i=1,nx
      fx(i) = ( x(i) - xc(i-1) )/( xc(i) - xc(i-1) )
    ENDDO

    DO i=1,ny
      fy(i) = ( y(i) - yc(i-1) )/( yc(i) - yc(i-1) )
    ENDDO

    DO i=1,nz
      fz(i) = ( z(i) - zc(i-1) )/( zc(i) - zc(i-1) )
    ENDDO

   END SUBROUTINE metrics
!------------------------------------------------------------------------------
