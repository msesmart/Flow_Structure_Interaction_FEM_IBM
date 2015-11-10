   SUBROUTINE FEA_STIFF(iBody)
!  FORM STIFF MATRIX

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface
!   USE flow_parameters, ONLY : nPtsMax

   IMPLICIT NONE
   
   INTEGER :: iBody,iEle,iNode 
   INTEGER :: m1,m2,m3,ITMP
   INTEGER :: I,J,ILOC,ISTF

   REAL(KIND=CGREAL) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
   REAL(KIND=CGREAL) :: axy,ayz,azx,area
   REAL(KIND=CGREAL) :: xll,xmm,xnn
   REAL(KIND=CGREAL) :: xb1,yb1,zb1,xb2,yb2,zb2,xb3,yb3,zb3,dxb,dyb,dzb
   REAL(KIND=CGREAL), DIMENSION(:,:) :: ekb(18,18), EK(18,18)
   REAL(KIND=CGREAL), DIMENSION(:,:) :: ek9(9,9)
   REAL(KIND=CGREAL) :: e0,g0,t0,znu,zip0,zia0,zib0,r0,pL0,alpha0,beta0

   LOGICAL :: OUTSTIFF

   OUTSTIFF = .false.

   STF(:) = 0

   DO iEle = 1,totNumTriElem(iBody) 

      m1 = triElemNeig(iBody,1,iEle)
      m2 = triElemNeig(iBody,2,iEle)
      m3 = triElemNeig(iBody,3,iEle)

      x1 = xBodyMarker(iBody,m1) 
      y1 = yBodyMarker(iBody,m1) 
      z1 = zBodyMarker(iBody,m1) 

      x2 = xBodyMarker(iBody,m2) 
      y2 = yBodyMarker(iBody,m2) 
      z2 = zBodyMarker(iBody,m2) 

      x3 = xBodyMarker(iBody,m3) 
      y3 = yBodyMarker(iBody,m3) 
      z3 = zBodyMarker(iBody,m3) 

!set material type to be 1 temporarily \\Wanh
!      e0 = PROPERTY(ELMatType(iEle),1)          ! E
!      g0 = PROPERTY(ELMatType(iEle),2)          ! G
!      t0 = PROPERTY(ELMatType(iEle),3)          ! Thickness of plate
!      r0 = PROPERTY(ELMatType(iEle),4)          ! Density of material
!      pL0 = PROPERTY(ELMatType(iEle),5)         ! pL0>0 plane stress

      e0 = PROPERTY(1,1)
      g0 = PROPERTY(1,2)
      t0 = PROPERTY(1,3)
      r0 = PROPERTY(1,4)      
      pL0 = PROPERTY(1,5)      
      
      znu = e0/(2.0*g0) - 1.0                   ! Possion's ratio
      zip0 = t0*t0*t0/(1-znu*znu)/12.0          ! Bending second moment of area for uniform plate
!      zia0 = PROPERTY(ELMatType(iEle),7)        ! In-plane drilling parameter alpha
!      zib0 = PROPERTY(ELMatType(iEle),8)        ! In-plane drilling parameter beta
      zia0 = 1.5
      zib0 = 0.5

!     Determine vector area
      axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
      ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
      azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
      area=sqrt( axy*axy + ayz*ayz + azx*azx)
      xll=ayz/area
      xmm=azx/area
      xnn=axy/area

!     Transform element global coords to local X-Y
!     xb1,yb1,zb1 are saved for DKT
      xb1=x1
      yb1=y1
      zb1=z1

      CALL ROTVEC(xll,xmm,xnn,(x2-x1),(y2-y1),(z2-z1),dxb,dyb,dzb) 
      xb2=x1+dxb
      yb2=y1+dyb
      zb2=z1+dzb

      CALL ROTVEC(xll,xmm,xnn,(x3-x1),(y3-y1),(z3-z1),dxb,dyb,dzb)
      xb3=x1+dxb
      yb3=y1+dyb
      zb3=z1+dzb

!     Calculate the LOCAL stiffness matrix
      ekb(:,:) = 0

      alpha0 = zia0
      beta0 = zib0

      CALL elmstfMRT(e0,g0,t0,pL0,alpha0,beta0,area,xb1,xb2,xb3,yb1,yb2,yb3,ekb,ek9)

!     Flexure
      CALL ELMSTFDKT(e0,g0,t0,zip0,area,xb1,xb2,xb3,yb1,yb2,yb3,ekb,ek9)

!     Membrane
!      CALL elmstfCST(e0,g0,t0,pL0,area,xb1,xb2,xb3,yb1,yb2,yb3,ekb)

!     Rotate to GLOBAL and assemble
 59   FORMAT(1x,18(g11.4))
      CALL ROTATE(xll,xmm,xnn,ekb,EK)

      CALL ASSEMBCOL(STF,EK,m1,m2,m3,iBody)

   ENDDO ! end loop of iEle


!  STORE stiffness matrix
   IF (OUTSTIFF) THEN
!      OPEN (21,FILE='stf.out')
!      WRITE(21,'(a)') 'STIFFNESS: COLumn form'
!      WRITE(21,*) 'STIFFNESS: COLumn form'
      DO I=1,NEQ
         ILOC=NLOC(I)
!         WRITE(21,22) (STF(ILOC+J-1), J=1,IPROF(I))
         if (i.eq.Neq) print *, 'iloc,iprof(i)=',iloc,iprof(i)

         WRITE(*,22) (STF(ILOC+J-1), J=1,IPROF(I))
 22      FORMAT(1X,6(E13.6))
      ENDDO
!      CLOSE (21)
   ENDIF

!   CALL STOREcol(STF,ISTF,iBody)

   END SUBROUTINE FEA_STIFF



   SUBROUTINE ROTVEC(xll,xmm,xnn,dx,dy,dz,dxb,dyb,dzb)
!  ROTate VECtor, makes a 3-D coordinate transformations.

   USE global_parameters

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: xll,xmm,xnn  !Input
   REAL(KIND=CGREAL) :: dx,dy,dz     !Input
   REAL(KIND=CGREAL) :: dxb,dyb,dzb  !Output

   REAL(KIND=CGREAL), DIMENSION(3,3) :: r
   REAL(KIND=CGREAL) :: ddd

   ddd=sqrt(1-xnn**2)   ! This is equivalent to d=sqrt(l_x^2+m_x^2)

   IF (ABS(xnn) .gt. 0.9999) THEN
      r(1,1) = +xnn
      r(1,2) = +0.0
      r(1,3) = -0.0
      r(2,1) = -0.0
      r(2,2) =  1.0
      r(2,3) =  0.0
      r(3,1) =  0.0
      r(3,2) =  0.0
      r(3,3) =  xnn
      dxb = r(1,1)*dx
      dyb = r(2,2)*dy
      dzb = r(3,3)*dz
      RETURN

   ELSEIF (ABS(xll) .gt. 0.9999) THEN
      r(1,1) = +0.0
      r(1,2) = +0.0
      r(1,3) = -1.0
      r(2,1) = -0.0
      r(2,2) = +xll
      r(2,3) =  0.0
      r(3,1) = +xll
      r(3,2) =  0.0
      r(3,3) =  0.0
      dxb = r(1,3)*dz
      dyb = r(2,2)*dy
      dzb = r(3,1)*dx
      RETURN

   ELSE
      r(1,1)  = +xll*xnn/ddd
      r(1,2)  = +xmm*xnn/ddd
      r(1,3)  = -ddd
      r(2,1)  =  -xmm/ddd
      r(2,2)  =  +xll/ddd
      r(2,3)  =  0.0
      r(3,1)  =  xll
      r(3,2)  =  xmm
      r(3,3)  =  xnn
   ENDIF

!  [dxb] = [R][dx]
   dxb = r(1,1)*dx + r(1,2)*dy + r(1,3)*dz
   dyb = r(2,1)*dx + r(2,2)*dy + r(2,3)*dz
   dzb = r(3,1)*dx + r(3,2)*dy + r(3,3)*dz

!  Return dxb,dyb,dzb

   END SUBROUTINE ROTVEC




!----------------------------------------------------------------
!  ELeMent STiFfness for Discrete Kirchhoff Triangle 
!  OUTPUT: ek, ekb

   SUBROUTINE ELMSTFDKT(e0,g0,tt0,zip0,a,x1,x2,x3,y1,y2,y3,ekdkt,ekbdkt)
   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: e0,g0,tt0,zip0,a
   REAL(KIND=CGREAL) :: x1,x2,x3,y1,y2,y3
   REAL(KIND=CGREAL) :: ekdkt(18,18),ekbdkt(9,9)
   REAL(KIND=CGREAL) :: d(3,3),dd(9,9),qq(9,9),pp(3,3),pt(2,3),rs(2,3),q(3)
   REAL(KIND=CGREAL) :: gg(10,9),b(3),c(3),als(3),px(3,3)
   REAL(KIND=CGREAL) :: p1,p2,p3,r1,r2,r3
   REAL(KIND=CGREAL) :: sum,dc,det
   REAL(KIND=CGREAL) :: znu


   INTEGER :: kod(2,9),inew(9)
   INTEGER :: I,J,K,II,JJ,L,K1,K2

   data kod/1,1,2,3,3,2,4,4,5,6,6,5,7,7,8,9,9,8/
   data pp /12.0,4.0,4.0,4.0,2.0,1.0,4.0,1.0,2.0/

!  znu is the Possion's ratio, transferred from: G=E/(1+2*znu)
   znu = e0/(2.0*g0) - 1.0

   dc=(e0*tt0**3)/12.0/(1-znu*znu)  !Coefficient between M(Moment) and kapa (Curvature)
   dc=(e0*zip0)/(1-znu*znu) !zip0 is "bending second moment, Ip=h^3/12/(1-nu*nu)"

!  D matrix (relation between stree and strain.
   d(1,1) = dc
   d(1,2) =dc*znu
   d(1,3) = 0.0
   d(2,1) =dc*znu
   d(2,2) = dc
   d(2,3) = 0.0
   d(3,1) = 0.0
   d(3,2) = 0.0
   d(3,3) = dc*(1-znu)/2.0

   b(1) = y2-y3
   b(2) = y3-y1
   b(3) = y1-y2
   c(1) = x3-x2
   c(2) = x1-x3
   c(3) = x2-x1
   det = 24.0*(b(1)*c(2)-b(2)*c(1))

   do i=1,3
      do j=1,3
         px(i,j) = pp(i,j)/det
      enddo
   enddo 

   do 25 i=1,3
      do 25 j=1,3
         do 25 k1=1,3
            ii=(i-1)*3+k1
            do 25 k2=1,3
               jj=(j-1)*3+k2
               dd(ii,jj)=d(i,j)*px(k1,k2)
 25      continue

         do 30 i=1,3
            als(i) = b(i)*b(i)+c(i)*c(i)
            pt(1,i) = 6.0*c(i)/als(i)
            pt(2,i) = 6.0*b(i)/als(i)
            rs(1,i) = 3.0*c(i)*c(i)/als(i)
            rs(2,i) = 3.0*b(i)*b(i)/als(i)
            q(i)    = 3.0*b(i)*c(i)/als(i)
 30      continue

         do 720 i=1,10
            do 720 j=1,9
               gg(i,j) = 0.0
 720     continue

         do 730 i=1,2
            ii=(i-1)*5
            p1=pt(i,1)
            p2=pt(i,2)
            p3=pt(i,3)
            r1=rs(i,1)
            r2=rs(i,2)
            r3=rs(i,3)
            gg(ii+1,kod(i,1)) =  p3
            gg(ii+2,kod(i,1)) = -p2
            gg(ii+3,kod(i,1)) = -p3
            gg(ii+4,kod(i,1)) =  p2-p3
            gg(ii+5,kod(i,1)) =  p2

            gg(ii+1,kod(i,2)) = -q(3)
            gg(ii+2,kod(i,2)) = -q(2)
            gg(ii+3,kod(i,2)) =  q(3)
            gg(ii+4,kod(i,2)) =  q(2)+q(3)
            gg(ii+5,kod(i,2)) =  q(2)

            gg(ii+1,kod(i,3)) = -1.0-r3
            gg(ii+2,kod(i,3)) = -1.0-r2
            gg(ii+3,kod(i,3)) =  r3
            gg(ii+4,kod(i,3)) =  r2+r3
            gg(ii+5,kod(i,3)) =  r2

            gg(ii+1,kod(i,4)) = -p3
            gg(ii+3,kod(i,4)) =  p3
            gg(ii+4,kod(i,4)) = p1+p3

            gg(ii+1,kod(i,5)) = -q(3)
            gg(ii+3,kod(i,5)) =  q(3)
            gg(ii+4,kod(i,5)) =  q(3)-q(1)

            gg(ii+1,kod(i,6)) =  1.0-r3
            gg(ii+3,kod(i,6)) =  r3
            gg(ii+4,kod(i,6)) =  r3-r1

            gg(ii+2,kod(i,7)) =  p2
            gg(ii+4,kod(i,7)) = -p1-p2
            gg(ii+5,kod(i,7)) = -p2

            gg(ii+2,kod(i,8)) = -q(2)
            gg(ii+4,kod(i,8)) =  q(2)-q(1)
            gg(ii+5,kod(i,8)) =  q(2)

            gg(ii+2,kod(i,9)) =  1.0-r2
            gg(ii+4,kod(i,9)) =  r2-r1
            gg(ii+5,kod(i,9)) =  r2
 730     continue

         do 850 i=1,9
            qq(1,i) =     b(2)*gg(1,i) +     b(3)*gg(2,i)
            qq(2,i) = 2.0*b(2)*gg(3,i) +     b(3)*gg(4,i)
            qq(3,i) =     b(2)*gg(4,i) + 2.0*b(3)*gg(5,i)
            qq(4,i) =    -c(2)*gg(6,i) -     c(3)*gg(7,i)
            qq(5,i) =-2.0*c(2)*gg(8,i) -     c(3)*gg(9,i)
            qq(6,i) =-    c(2)*gg(9,i) - 2.0*c(3)*gg(10,i)
            qq(7,i) =     c(2)*gg(1,i) +     c(3)*gg(2 ,i) &
                         -b(2)*gg(6,i) -     b(3)*gg(7,i)
            qq(8,i) = 2.0*c(2)*gg(3,i) +     c(3)*gg(4 ,i) &
                    - 2.0*b(2)*gg(8,i) -     b(3)*gg(9,i)
            qq(9,i) =     c(2)*gg(4,i) + 2.0*c(3)*gg(5 ,i) &
                    -     b(2)*gg(9,i) - 2.0*b(3)*gg(10,i)
 850     continue

         do 855 i=1,9
            do 855 j=1,9
               gg(i,j) = 0.0
               do 855 k=1,9
                  gg(i,j) = gg(i,j) + dd(i,k)*qq(k,j)
 855      continue
!
         do 960 L=1,9
            do 960 j=L,9
               sum = 0.0
               do 900 k=1,9
                  sum = sum  + qq(k,L)*gg(k,j)
 900           continue

               ekbdkt(L,j)=sum
               ekbdkt(j,L)=sum
 960      continue
!
! ?????????????????????????????????????????
!  assign to full element matrix
         inew(1)=3
         inew(2)=4
         inew(3)=5
         inew(4)=9
         inew(5)=10
         inew(6)=11
         inew(7)=15
         inew(8)=16
         inew(9)=17
   do i=1,9
      ii=inew(i)
      do j=1,9
         jj=inew(j)
         ekdkt(ii,jj) = ekbdkt(i,j)
      enddo
   enddo
!        print *, 'ekdkt(1,1) before exit dkt is:',ekdkt(1,1)

   END SUBROUTINE ELMSTFDKT




   SUBROUTINE elmstfCST(e0,g0,tt0,pL0,a,x1,x2,x3,y1,y2,y3,ek)
   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE
   
   INTEGER :: inew(9), i,j,k, ii,jj
   REAL(KIND=CGREAL) :: e0,g0,tt0,pL0,a,x1,x2,x3,y1,y2,y3,a2,ek(18,18)
   REAL(KIND=CGREAL) :: ekb(6,6),b(3,6),c(3,3),etemp(3,6)
   REAL(KIND=CGREAL) :: nv, xkap
   REAL(KIND=CGREAL) :: att0
   
   ekb(:,:)=0.
   etemp(:,:)=0.

   nv = (e0/g0/2.0)-1.0   !Possion's ratio
!  plane strain
   if (pL0 .lt. 0) xkap=3.-4.*nv
!  plane stress
   if (pL0 .gt. 0) xkap=(3.0-nv)/(1.+nv)

!  {s}=[C]{e}
   c(1,1)=g0/(xkap-1.)*(xkap+1.)
   c(1,2)=g0/(xkap-1.)*(3.-xkap)
   c(1,3)=0.
   c(2,1)=c(1,2)
   c(2,2)=c(1,1)
   c(2,3)=0.
   c(3,1)=0.
   c(3,2)=0.
   c(3,3)=g0

!  {e}=[B]{u}
   a=((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.   ! area of element
   a2=a*2.0

   b(1,1)=(y2-y3)/a2
   b(1,2)=0.
   b(1,3)=(y3-y1)/a2
   b(1,4)=0.
   b(1,5)=(y1-y2)/a2
   b(1,6)=0.

   b(2,1)=0.
   b(2,2)=(x3-x2)/a2
   b(2,3)=0.
   b(2,4)=(x1-x3)/a2
   b(2,5)=0.
   b(2,6)=(x2-x1)/a2

   b(3,1)=b(2,2)
   b(3,2)=b(1,1)
   b(3,3)=b(2,4)
   b(3,4)=b(1,3)
   b(3,5)=b(2,6)
   b(3,6)=b(1,5)

!  [k] = At[B][C][B]
   att0=a*tt0

   do 30 i=1,3
      do 32 j=1,6
         do 34 k=1,3
            etemp(i,j)=etemp(i,j)+c(i,k)*b(k,j)
34       continue
32    continue
30 continue

   do 40 i=1,6
      do 42 j=1,6
         do 44 k=1,3
!           b(k,i) is used due to transpose
            ekb(i,j)=ekb(i,j)+b(k,i)*etemp(k,j)*att0
44       continue
42    continue
40 continue

!  assign [ekb] to [18x18]
   inew(1)=1
   inew(2)=2
   inew(3)=7
   inew(4)=8
   inew(5)=13
   inew(6)=14
   do i=1,6
      ii=inew(i)
      do j=1,6
         jj=inew(j)
         ek(ii,jj) = ekb(i,j)
      enddo
   enddo

   END SUBROUTINE elmstfCST



   SUBROUTINE ROTATE  (XLL,XMM,XNN,EKB,EK)
!  ROTATE stiffness matrix
!  MAkes a 2-D coordinate transformations.

   USE global_parameters
   IMPLICIT NONE
   
   REAL(KIND=CGREAL) XLL,XMM,XNN
   REAL(KIND=CGREAL) EK(18,18),EKB(18,18)
   REAL(KIND=CGREAL) rt(3,3),r(3,3),ektemp(18,18)
   REAL(KIND=CGREAL) ddd
   INTEGER :: i,j,k, ii,jj,kk,in,jn,j1,j2
   ddd=sqrt(1-XNN**2)

   if (abs(XNN) .gt. 0.9999) then
      r(1,1)  = +xnn
      r(1,2)  = +0.0
      r(1,3)  = -0.0
      r(2,1)  = -0.0
      r(2,2)  =  1.0
      r(2,3)  =  0.0
      r(3,1)  =  0.0
      r(3,2)  =  0.0
      r(3,3)  =  xnn
   elseif (abs(XLL) .gt. 0.9999) then
      r(1,1)  = +0.0
      r(1,2)  = +0.0
      r(1,3)  = -1.0
      r(2,1)  = -0.0
      r(2,2)  = +xll
      r(2,3)  =  0.0
      r(3,1)  = +xll
      r(3,2)  =  0.0
      r(3,3)  =  0.0
   else
      r(1,1)  = +xll*xnn/ddd    ! Refer to Doyle "Nonlinear book", pp 159
      r(1,2)  = +xmm*xnn/ddd
      r(1,3)  = -ddd
      r(2,1)  = -xmm/ddd
      r(2,2)  = +xll/ddd
      r(2,3)  = 0.0
      r(3,1)  = xll
      r(3,2)  = xmm
      r(3,3)  = xnn
   endif

!  take [Rtrans][k][R] using the nature of [R] to speed computation.
!  [k] is sectioned off into 6 3x3s then multiplied [rtrans][k][r]
!  Refer to Doyle "Nonlinear book", pp 156-157

   IF (ABS(XNN) .GT. 0.9999) THEN
      do i=0,5
         do j=0,5
            ii=i*3
            jj=j*3
            ek(ii+1,jj+1) = ekb(ii+1,jj+1)
            ek(ii+1,jj+2) = ekb(ii+1,jj+2)*r(1,1)
            ek(ii+1,jj+3) = ekb(ii+1,jj+3)
            ek(ii+2,jj+1) = ekb(ii+2,jj+1)*r(1,1)
            ek(ii+2,jj+2) = ekb(ii+2,jj+2)
            ek(ii+2,jj+3) = ekb(ii+2,jj+3)*r(1,1)
            ek(ii+3,jj+1) = ekb(ii+3,jj+1)
            ek(ii+3,jj+2) = ekb(ii+3,jj+2)*r(1,1)
            ek(ii+3,jj+3) = ekb(ii+3,jj+3)
         enddo
      enddo
      return
   ELSEIF (ABS(XLL) .GT. 0.9999) THEN
      do i=0,5
         do j=0,5
            ii=i*3
            jj=j*3
            ek(ii+1,jj+1) = ekb(ii+3,jj+3)
            ek(ii+1,jj+2) = ekb(ii+3,jj+2)
            ek(ii+1,jj+3) =-ekb(ii+3,jj+1)*r(2,2)
            ek(ii+2,jj+1) = ekb(ii+2,jj+3)
            ek(ii+2,jj+2) = ekb(ii+2,jj+2)
            ek(ii+2,jj+3) =-ekb(ii+2,jj+1)*r(2,2)
            ek(ii+3,jj+1) =-ekb(ii+1,jj+3)*r(2,2)
            ek(ii+3,jj+2) =-ekb(ii+1,jj+2)*r(2,2)
            ek(ii+3,jj+3) = ekb(ii+1,jj+1)
         enddo
      enddo
      return
   ENDIF

!  get transpose

   do in=1,3
      do jn=1,3
         rt(jn,in)=r(in,jn)
      enddo
   enddo

   do 22 i=0,5
      do 22 j=0,5
         j1=i*3
         j2=j*3
!        [k][R]
         do 23 k=1,3
         do 23 ii=1,3
            ektemp(j1+k,j2+ii)=0.0
            do 23 jj=1,3
               ektemp(j1+k,j2+ii)=ektemp(j1+k,j2+ii)+ekb(j1+k,j2+jj)*r(jj,ii)
 23       continue

!        [Rtrans][k]
         do 24 k=1,3
         do 24 ii=1,3
            ek(j1+k,j2+ii)=0.0
            do 24 jj=1,3
               ek(j1+k,j2+ii)=ek(j1+k,j2+ii)+rt(k,jj)*ektemp(j1+jj,j2+ii)
 24      continue
 22     continue

   END SUBROUTINE ROTATE




   SUBROUTINE TRANS3D (L,M,N,EK,BETA)
!  This subroutine makes 3-D coordinate transformations.

   USE global_parameters
   IMPLICIT NONE

   REAL(KIND=CGREAL) :: ek(12,12),rt(3,3),r(3,3),ktemp(12,12)
   REAL(KIND=CGREAL) :: M,N,L,BETA,pi
   REAL(KIND=CGREAL) :: sb,cb,d
   INTEGER :: in,jn,i,j,k,j1,j2,ii,jj

   pi=4.0*atan(1.0)

   sb=sin(beta*pi/180)
   cb=cos(beta*pi/180)
   d=sqrt(1-n**2)

!  if (abs(l).ge. 0.995 .and. abs(beta).le. 0.01) return
   IF (ABS(N).GT.0.995) THEN
      r(1,1)  =  0.0
      r(1,2)  =  0.0
      r(1,3)  =  n
      r(2,1)  = -n*sb
      r(2,2)  =  cb
      r(2,3)  =  0.0
      r(3,1)  = -n*cb
      r(3,2)  = -sb
      r(3,3)  =  0.0
   ELSE
      r(1,1)  =  l
      r(1,2)  =  m
      r(1,3)  =  n
      IF (ABS(BETA) .LE. .01) THEN
              r(2,1)  =  -m/d
              r(2,2)  =  l/d
              r(2,3)  =  0.0
              r(3,1)  =  -l*n/d
              r(3,2)  =  -m*n/d
              r(3,3)  =  d
      ELSE
              r(2,1)  =  -(m*cb+l*n*sb)/d
              r(2,2)  =  (l*cb-m*n*sb)/d
              r(2,3)  =  d*sb
              r(3,1)  =  (m*sb-l*n*cb)/d
              r(3,2)  =  -(l*sb+m*n*cb)/d
              r(3,3)  =  d*cb
      ENDIF
   ENDIF

   DO IN=1,3
      DO JN=1,3
         rt(jn,in)=r(in,jn)
      ENDDO
   ENDDO

!  take [Rtrans][K][R] using the nature of [R] to speed computation.
!  k is sectioned off into 3x3s then multiplied [rtrans][k][r]

   DO 22 I=0,3
      DO 22 J=0,3
         do 23 k=1,3
         do 23 ii=1,3
            j1=i*3
            j2=j*3
            ktemp(j1+k,j2+ii)=0.0
            do 23 jj=1,3
               ktemp(j1+k,j2+ii)=ktemp(j1+k,j2+ii)+ek(j1+k,j2+jj)*r(jj,ii)
 23      continue

         do 24 k=1,3
            do 24 ii=1,3
            ek(j1+k,j2+ii)=0.0
            do 24 jj=1,3
            ek(j1+k,j2+ii)=ek(j1+k,j2+ii)+rt(k,jj)*ktemp(j1+jj,j2+ii)
 24      continue

 22 CONTINUE

   END SUBROUTINE TRANS3D


   SUBROUTINE ASSEMBCOL(AA,A,I1,J1,K1,iBody)
!  ASSEMBle element matrices in COLumn form

   USE global_parameters
   USE fea_unstructure_surface
   USE flow_parameters, ONLY : nPtsMax

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: AA((nPtsMax*nodeDoF)**2), A(18,18)
   INTEGER :: I1,J1,K1
   INTEGER :: iBody

   INTEGER :: IDOF(18)
   INTEGER :: I,J,IEQN1,IEQN2,JBAND,ILOC

!  Set idof to posible DoF number of each nodes
   DO I= 1,6
      IDOF(I)    = (I1-1)*6 + I
      IDOF(I+6)  = (J1-1)*6 + I
      IDOF(I+12) = (K1-1)*6 + I
   ENDDO

!  Store the values for individual array in global array
   DO 20 I = 1,18

      IEQN1 = JBC(IDOF(I))

      IF (IEQN1 .GT. 0) THEN
         DO J= I,18
            IEQN2 = JBC(IDOF(J))
            IF (IEQN2 .GT. 0) THEN
               IF (IEQN1 .GT. IEQN2) THEN
                  jband= (ieqn1-ieqn2)+1
                  iloc = NLOC(ieqn1)
                  aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
!                if (iloc +jband-1 .eq. 102) print *, 'aa(102)=', aa(iloc +jband-1)

!                 aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
               ELSE
                  jband= (ieqn2-ieqn1)+1
                  iloc = nloc(ieqn2)
                  aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
!                if (iloc +jband-1 .eq. 102) print *, 'aa(102)=', aa(iloc +jband-1)
!                 aa(ieqn1,jband) = aa(ieqn1,jband) + a(i,j)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

 20 CONTINUE

   END SUBROUTINE ASSEMBCOL



   SUBROUTINE STOREcol(AA,iFile,iBody)
!  STORE COLumn oriented matrix

   USE global_parameters
   USE fea_unstructure_surface
   USE flow_parameters, ONLY : nPtsBodyMarker
   
   IMPLICIT NONE
   
   INTEGER :: iFile, iBody
   INTEGER :: i,iloc,imax,j

   REAL(KIND=CGREAL) :: AA(nPtsBodyMarker(iBody))

   do i=1,Neq
      iloc = nloc(i)
      imax = iprof(i)
      write(iFile) ( AA(iloc+j-1), j=1,imax )
   enddo

   END SUBROUTINE STOREcol 



!  ELeMent STiFFness for Membrane with Rotation Triangle
   SUBROUTINE elmstfMRT(e0,g0,tt0,pL0,alpha0,beta0,a,x1,x2,x3,y1,y2,y3,ek,ekb)
   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE
   
   REAL(KIND=CGREAL) :: e0,g0,tt0,pL0,alpha0,beta0,a
   REAL(KIND=CGREAL) :: x1,x2,x3,y1,y2,y3
   REAL(KIND=CGREAL) :: ek(18,18),ekb(9,9)
   REAL(KIND=CGREAL) :: xt(3),yt(3),dmt(3,3),sm(9,9),ektemp(9,9)
   REAL(KIND=CGREAL) :: alpha,beta,f,fbeta
   REAL(KIND=CGREAL) :: znu0, xkap
   INTEGER :: inew(9)
   INTEGER ::  Lst(12)
   INTEGER :: i,j,m,ii,jj
   INTEGER :: status
   
   data Lst /1,2,4,5,7,8,3,6,9,10,11,12/

!ccccccccc!!!!!!!!!!!!!! change plane stress or strain !!!!!!!!
!        znu=(e0/g0/2.0)-1.0
!        coef=e0*tt0/(1-znu*znu)
!        dmt(1,1) = 1.0*coef
!        dmt(1,2) = znu*coef
!        dmt(1,3) = 0.0
!        dmt(2,1) = znu*coef
!        dmt(2,2) = 1.0*coef
!        dmt(2,3) = 0.0
!        dmt(3,1) = 0.0
!        dmt(3,2) = 0.0
!        dmt(3,3) = (1-znu)*0.5*coef

         znu0 = (e0/g0/2.0)-1.0
!        plane stress
         if (pL0 .gt. 0) xkap=(3.-znu0)/(1.+znu0)
!        plane strain
         if (pL0 .lt. 0) xkap=3.-4.*znu0

!        {s}=[C]{e}
!        dmt is general D matrix, xkap is used, equivalent to /nv version
         dmt(1,1)=g0*tt0/(xkap-1.)*(xkap+1.)
         dmt(1,2)=g0*tt0/(xkap-1.)*(3.-xkap)
         dmt(1,3)=0.
         dmt(2,1)=dmt(1,2)
         dmt(2,2)=dmt(1,1)
         dmt(3,2)=0.
         dmt(3,1)=0.
         dmt(3,2)=0.
         dmt(3,3)=g0*tt0

         xt(1)=x1
         xt(2)=x2
         xt(3)=x3
         yt(1)=y1
         yt(2)=y2
         yt(3)=y3

         f=1.0
         alpha=alpha0   ! Ref: Bergan & Felippa 1985
         beta =beta0

         do i=1,9
            do j=1,9
               ekb(i,j)=0.0
            enddo
         enddo
         
         m=9
         call sm3mb(xt,yt,dmt,alpha,f,Lst,ekb,m,status)

         if (beta .gt. 0.0) then
             fbeta=f*beta
             call sm3mh(xt,yt,dmt,fbeta,Lst,ekb,m,status)
         endif

!        assign [ekb] to [18x18]
         inew(1)=1
         inew(2)=2
         inew(3)=6
         inew(4)=7
         inew(5)=8
         inew(6)=12
         inew(7)=13
         inew(8)=14
         inew(9)=18
         do i=1,9
            ii=inew(i)
            do j=1,9
               jj=inew(j)
               ek(ii,jj) = ekb(i,j)
            enddo
         enddo

      END SUBROUTINE elmstfMRT


      SUBROUTINE sm3mb(x,y,dm,alpha,f,Ls,sm,m,status)
!     From Bergan and Felippa, CMAME, 1985
         USE global_parameters, ONLY : CGREAL
         IMPLICIT NONE
         
!     Basic stiffness
         REAL(KIND=CGREAL) :: x(3),y(3),dm(3,3),sm(9,9)
         REAL(KIND=CGREAL) :: alpha,f
         INTEGER :: Ls(12),m,n
         INTEGER :: i,j,k,L
         INTEGER :: status
         
         REAL(KIND=CGREAL) :: p(9,3)
         REAL(KIND=CGREAL) :: area2,c
         REAL(KIND=CGREAL) :: d11,d12,d13,d22,d23,d33
         REAL(KIND=CGREAL) :: x21,x32,x13,y21,y32,y13
         REAL(KIND=CGREAL) :: x12,x23,x31,y12,y23,y31
         REAL(KIND=CGREAL) :: s1,s2,s3
        
         x21 = x(2) - x(1)
         x32 = x(3) - x(2)
         x13 = x(1) - x(3)
         x12 = -x21
         x23 = -x32
         x31 = -x13

         y21 = y(2) - y(1)
         y32 = y(3) - y(2)
         y13 = y(1) - y(3)
         y12 = -y21
         y23 = -y32
         y31 = -y13

         area2 = y21*x13 - x21*y13

!        The following p matrix is L matrix in Doyle's book p139, or in Bergan &
!        Felippa 1985, Equ. 4.63, but the sequence has been changed. 
!        alpha (Drilling parameter) related entries are moved to row 7~9.

         p(1,1) = y23
         p(2,1) = 0.0
         p(3,1) = y31
         p(4,1) = 0.0
         p(5,1) = y12
         p(6,1) = 0.0
         p(1,2) = 0.0
         p(2,2) = x32
         p(3,2) = 0.0
         p(4,2) = x13
         p(5,2) = 0.0
         p(6,2) = x21
         p(1,3) = x32
         p(2,3) = y23
         p(3,3) = x13
         p(4,3) = y31
         p(5,3) = x21
         p(6,3) = y12
         n=6

         if (alpha .ne. 0.0) then
             p(7,1) = y23 * (y13-y21)*alpha/6.0
             p(7,2) = x32 * (x31-x12)*alpha/6.0
             p(7,3) =       (x31*y13-x12*y21)*alpha/3.0  !See Doyle "Nonliear book", pp 139, or Felippa 
             p(8,1) = y31 * (y21-y32)*alpha/6.0
             p(8,2) = x13 * (x12-x23)*alpha/6.0
             p(8,3) =       (x12*y21-x23*y32)*alpha/3.0
             p(9,1) = y12 * (y32-y13)*alpha/6.0
             p(9,2) = x21 * (x23-x31)*alpha/6.0
             p(9,3) =       (x23*y32-x31*y13)*alpha/3.0
             n=9
         endif

         c= 0.5*f/area2
         d11 = c*dm(1,1)
         d22 = c*dm(2,2)
         d33 = c*dm(3,3)
         d12 = c*dm(1,2)
         d13 = c*dm(1,3)
         d23 = c*dm(2,3)

         do j=1,n
            L  = ls(j)
            s1 = d11*p(j,1) + d12*p(j,2) + d13*p(j,3)
            s2 = d12*p(j,1) + d22*p(j,2) + d23*p(j,3)
            s3 = d13*p(j,1) + d23*p(j,2) + d33*p(j,3)
            do i=1,j
               k = ls(i)
               sm(k,L) = sm(k,L) + (s1*p(i,1) + s2*p(i,2) + s3*p(i,3))
               sm(L,k) = sm(k,L)
            enddo
         enddo
         
      END SUBROUTINE sm3mb


      SUBROUTINE sm3mh(x,y,dm,f,Ls,sm,m,status)
!        higher stiffness
         USE global_parameters, ONLY : CGREAL
         IMPLICIT NONE

         REAL(KIND=CGREAL) :: x(3),y(3),dm(3,3),f,sm(9,9)
         INTEGER :: Ls(12)
         INTEGER :: status
         
         REAL(KIND=CGREAL) :: xc(3), yc(3), xm(3),ym(3)
         REAL(KIND=CGREAL) :: bh(3,3),gt(9,9),hh(3,9)
         REAL(KIND=CGREAL) :: sqh(3,3), t(9), qx(3,3),qy(3,3)
         REAL(KIND=CGREAL) :: p(9,3)
         REAL(KIND=CGREAL) :: alpha,area,area2,c
         REAL(KIND=CGREAL) :: a1j, a2j, a3j, b1j, b2j,b3j
         REAL(KIND=CGREAL) :: d11,d12,d13,d22,d23,d33,jxx,jxy,jyy
         REAL(KIND=CGREAL) :: x0, y0, xi,yi
         REAL(KIND=CGREAL) :: cj, sj,dl,dx,dy
         REAL(KIND=CGREAL) :: s1,s2,s3,s4,s5,s6
         REAL(KIND=CGREAL) :: gti(9,9),wk(9,9)

         INTEGER :: iperm(9)
         INTEGER :: m
         INTEGER :: i,j,k,L
        
         area2 = (y(2) -y(1))*(x(1)-x(3)) - (x(2)-x(1))*(y(1)-y(3))
!*       write(*,*)'@@ AREA x 2: ',area2

         x0 = (x(1)+x(2)+x(3))/3.0
         y0 = (y(1)+y(2)+y(3))/3.0
         area=0.5*area2
         c = 1.0/sqrt(area)

         xc(1) = c*(x(1)-x0)
         xc(2) = c*(x(2)-x0)
         xc(3) = c*(x(3)-x0)
         yc(1) = c*(y(1)-y0)
         yc(2) = c*(y(2)-y0)
         yc(3) = c*(y(3)-y0)
         xm(1) = 0.5*(xc(2)+xc(3))
         xm(2) = 0.5*(xc(3)+xc(1))
         xm(3) = 0.5*(xc(1)+xc(2))
         ym(1) = 0.5*(yc(2)+yc(3))
         ym(2) = 0.5*(yc(3)+yc(1))
         ym(3) = 0.5*(yc(1)+yc(2))

!        form G^T in GT and initialize HH

          do i=1,9
             do j=1,6
                gt(j,i) = 0
             enddo
             hh(1,i)=0
             hh(2,i)=0
             hh(3,i)=0
         enddo

         d11 = f*dm(1,1)
         d22 = f*dm(2,2)
         d33 = f*dm(3,3)
         d12 = f*dm(1,2)
         d13 = f*dm(1,3)
         d23 = f*dm(2,3)
         jxx = -2.0*(xc(1)*xc(2) + xc(2)*xc(3) + xc(3)*xc(1))/3.0
         jxy =      (xc(1)*yc(1) + xc(2)*yc(2) + xc(3)*yc(3))/3.0
         jyy = -2.0*(yc(1)*yc(2) + yc(2)*yc(3) + yc(3)*yc(1))/3.0
         
         do 2500 j=1,3
            dx = xm(j) - xc(j)
            dy = ym(j) - yc(j)
            dl = sqrt(dx*dx + dy*dy)
            cj = dx/dl
            sj = dy/dl

!           !!!a2j b2j different than paper
            a1j = -0.5*sj*cj**2
            a2j =  0.5*cj**3
            b2j = -0.5*sj**3
            b3j =  0.5*sj**2*cj
            a3j = -(b2j + a1j + a1j)
            b1j = -(b3j + b3j + a2j)

            gt(1,2*j-1) =       1.
            gt(2,2*j  ) =   1.
            gt(3,2*j-1) =  -yc(j)
            gt(3,2*j  )=  xc(j)
            gt(3,  j+6)=   c
            gt(4,2*j-1)=  xc(j)
            gt(6,2*j-1)=  yc(j)
            gt(5,2*j  )=  yc(j)
            gt(6,2*j  )=  xc(j)
            hh(j,j+6) =     1.
            qx(j,1)     =      a1j
            qx(j,2)     =      b2j
            qx(j,3)     =     -2.0*b3j
            qy(j,1)     =      a2j
            qy(j,2)     =      b3j
            qy(j,3)     =     -2.0*a1j
            s1 =         d11*qx(j,1)    + d12*qx(j,2)   + d13*qx(j,3)
            s2 =         d12*qx(j,1)    + d22*qx(j,2)   + d23*qx(j,3)
            s3 =         d13*qx(j,1)    + d23*qx(j,2)   + d33*qx(j,3)
            s4 =         d11*qy(j,1)    + d12*qy(j,2)   + d13*qy(j,3)
            s5 =         d12*qy(j,1)    + d22*qy(j,2)   + d23*qy(j,3)
            s6 =         d13*qy(j,1)    + d23*qy(j,2)   + d33*qy(j,3)
            do i = 1,3
              xi =        xc(i)
              yi =        yc(i)
              gt(j+6,2*i-1) =    a1j*xi*xi + 2.*a2j*xi*yi + a3j*yi*yi
              gt(j+6,2*i)   =    b1j*xi*xi + 2.*b2j*xi*yi + b3j*yi*yi
              gt(j+6,i+6)   =   -c*(cj*xi+sj*yi)
            enddo

            do i=1,j
               sqh(i,j) = jxx*( qx(i,1)*s1+qx(i,2)*s2+qx(i,3)*s3)  &
                  + jxy*( qx(i,1)*s4+qx(i,2)*s5+qx(i,3)*s6   &
                         +qy(i,1)*s1+qy(i,2)*s2+qy(i,3)*s3)  &
                  + jyy*( qy(i,1)*s4+qy(i,2)*s5+qy(i,3)*s6)
            enddo
 2500   continue

!       Factor G' and backsolve to obtain H
!       Form physical stiffness and add to incoming SM
        do i=1,9
           do j=1,9
              gti(i,j)=gt(i,j)
           enddo
        enddo
        
        call ainver(gt,9,iperm,wk)
        
        do i=1,3
           do j=1,9
              hh(i,j)=gt(j,i+6)
           enddo
        enddo
        
        do 3000 j = 1,9
           L = ls(j)
           s1 = sqh(1,1)*hh(1,j) + sqh(1,2)*hh(2,j) + sqh(1,3)*hh(3,j)
           s2 = sqh(1,2)*hh(1,j) + sqh(2,2)*hh(2,j) + sqh(2,3)*hh(3,j)
           s3 = sqh(1,3)*hh(1,j) + sqh(2,3)*hh(2,j) + sqh(3,3)*hh(3,j)
           do i = 1,j
              k = Ls(i)
              sm(k,L) = sm(k,L) + (s1*hh(1,i) + s2*hh(2,i) + s3*hh(3,i))
              sm(L,k) = sm(k,L)
           enddo
 3000   continue

 81     format(1x,20(g11.5,1x))
!
      END SUBROUTINE sm3mh



