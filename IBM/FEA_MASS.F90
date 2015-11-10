   SUBROUTINE FEA_MASS(iBody)
!  FORM MASS matrix: lumped or consistent

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE
!   INTEGER :: nmprop(totNumTriElem(iBody))

!   REAL(KIND=CGREAL) :: qms(nPtsBodyMarker(iBody))
   REAL(KIND=CGREAL) :: em(18,18), emb(18,18), em12(12,12)
   REAL(KIND=CGREAL) :: rh0,area,tt0
   REAL(KIND=CGREAL) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
   REAL(KIND=CGREAL) :: axy,ayz,azx
   REAL(KIND=CGREAL) :: dxb,dyb,dzb
   REAL(KIND=CGREAL) :: xll,xmm,xnn
   REAL(KIND=CGREAL) :: xb1,yb1,zb1,xb2,yb2,zb2,xb3,yb3,zb3
   REAL(KIND=CGREAL) :: temp8

   INTEGER :: iBody
   INTEGER :: iecho
   INTEGER :: ibandm,imss,ilump
   INTEGER :: i,i1,j1,k1,nnp

   CHARACTER*40 STR40

   iecho = 2

!  Zero array before assembling
   mss(:)=0.0

!  Form the element form matrix, and assemble
   DO 50 i=1,totNumTriElem(iBody)

      i1 = triElemNeig(iBody,1,i)
      j1 = triElemNeig(iBody,2,i)
      k1 = triElemNeig(iBody,3,i)

!      tt0 = PROPERTY(ELMatType(i),3)
!      rh0 = PROPERTY(ELMatType(i),4)

!     Eelment type is set to 1 temporarily \\Wanh
      tt0 = PROPERTY(1,3)
      rh0 = PROPERTY(1,4)

      x1 = xBodyMarker(iBody,i1)
      y1 = yBodyMarker(iBody,i1)
      z1 = zBodyMarker(iBody,i1)

      x2 = xBodyMarker(iBody,j1)
      y2 = yBodyMarker(iBody,j1)
      z2 = zBodyMarker(iBody,j1)

      x3 = xBodyMarker(iBody,k1)
      y3 = yBodyMarker(iBody,k1)
      z3 = zBodyMarker(iBody,k1)

!     Determine vector area
      axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
      ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
      azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
      area=sqrt( axy*axy + ayz*ayz + azx*azx)
      xll=ayz/area
      xmm=azx/area
      xnn=axy/area

!     Transform element global coords to local X-Y
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

!     Calculate the mass matrix and assemble 
      emb(:,:)=0.0

      ilump = 2

      call elmmasCST(rh0,area,tt0,emb,ilump)
      call elmmasMRT(rh0,area,tt0,emb,xb1,xb2,xb3,yb1,yb2,yb3)
      call elmmasPLT(rh0,area,tt0,emb,xb1,xb2,xb3,yb1,yb2,yb3,ilump)

!      print *, 'emb(1,1)   after PLT=', emb(1,1)
!      print *, 'emb(18,18) after PLT=', emb(18,18)


      call rotate(xll,xmm,xnn,emb,em)

!      if (ilump .eq. 1) then
!         call assembLUM(mss,em,i1,j1,k1)
!      elseif (ilump .eq. 2) then

       call assembCOL(mss,em,i1,j1,k1,iBody)

!      endif
 50 CONTINUE 

!! STORE  mass matrices on disk in case
!         if (iecho .eq. 2) then
!            if (ilump.eq.1) then    
!                write(iout,*)'MASS: diagonal'
!                write(iout,22) (mss(i), i= 1, neq )
!            else
!               write(iout,'(a)')'MASS: COLumn form '
!               do 11 i=1,neq  
!                  iloc=nloc(i)
!                  jmax=iprof(i)
!                  write(iout,22) (mss(iloc+j-1),  j=1,jmax)
! 11            continue
!            endif
! 22         format(1x,6(g13.6))
!         endif
!         if (iecho .eq. 3) then
!            if (ilump.eq.1) then    
!                write(iout,'(a)')'MASS: lumped'
!                write(iout,*) neq,ilump
!                write(iout,'(1x,20(g12.6,1x))') (mss(i), i= 1, neq )
!            else
!               write(iout,'(a)')'MASS: band storage'
!               call storeBND_f(mss,maxstiff,iprof,nloc,iout,neq,iband,qms)
!            endif
!         endif
!         if (iecho .eq. 4) then
!            read(ikbd,'(a)') str40 
!             write(ilog,*) str40,'                 :: filename'
!            open(unit=itmp,file=str40 ,form='unformatted')
!            rewind itmp
!            if (ilump.eq.1) then    
!                write(itmp) neq,ilump
!                write(itmp) (mss(i), i= 1, neq )
!            else
!               call storeBND(mss,maxstiff,iprof,nloc,itmp,neq,iband,qms)
!            endif
!            close (itmp)
!         endif
!

!         if (ilump.eq.1) then    
!             ibandm = 1 
!         else
             ibandm = IBAND
!         endif
!         rewind(imss)
!         if (ilump.eq.1) then    
!             write(imss) neq, ibandm  
!             do i=1,neq
!                write(imss) mss(i)
!             enddo
!         else
!             call STOREcol(mss,imss,iBody)
!         endif
!
         write(*,*) '@@ FORMMASS:  Formed  [M]   OK'

!        print *, 'mss(1)=',mss(1)
   END SUBROUTINE FEA_MASS


   SUBROUTINE elmmasCST(rho,area,th,em,ilump)
!  ELeMent MASs matrix for Constant Strain Triangle
!
   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: rho,area,th
   REAL(KIND=CGREAL) :: em(18,18), emb(9,9)
   REAL(KIND=CGREAL) :: roat
   INTEGER :: I,J,II,JJ,inew(9),ilump

   emb(:,:) = 0.0

   if (ilump.eq.1) then
!     contibutions to lumped mass matrix
      roat = rho*area*th/3.0
      DO I=1,9
         emb(i,i) = roat
      ENDDO

   elseif (ilump.eq.2) then
!     consistent mass matrix
      roat = rho*area*th/12.0
      emb(1,1) = roat*2
      emb(1,4) = roat
      emb(1,7) = roat
      emb(4,1) = roat
      emb(4,4) = roat*2
      emb(4,7) = roat
      emb(7,1) = roat
      emb(7,4) = roat
      emb(7,7) = roat*2
      emb(2,2) = roat*2
      emb(2,5) = roat
      emb(2,8) = roat
      emb(5,2) = roat
      emb(5,5) = roat*2
      emb(5,8) = roat
      emb(8,2) = roat
      emb(8,5) = roat
      emb(8,8) = roat*2
      emb(3,3) = roat*2
      emb(3,6) = roat
      emb(3,9) = roat
      emb(6,3) = roat
      emb(6,6) = roat*2
      emb(6,9) = roat
      emb(9,3) = roat
      emb(9,6) = roat
      emb(9,9) = roat*2
   endif

!  assign to full element matrix
   inew(1)=1                    !u1
   inew(2)=2                    !v1
   inew(3)=3                    !w1
   inew(4)=7                    !u2
   inew(5)=8
   inew(6)=9
   inew(7)=13
   inew(8)=14
   inew(9)=15
   do i=1,9
      ii=inew(i)
      do j=1,9
         jj=inew(j)
         em(ii,jj) = emb(i,j)
      enddo
   enddo

   END SUBROUTINE elmmasCST


!  ELeMent MASs matrix for Moment Rotation Triangle
   SUBROUTINE elmmasMRT(rho,AREA,TH,EM,X1,X2,X3,Y1,Y2,Y3)

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   REAL(KIND=CGREAL) :: rho,AREA,TH
   REAL(KIND=CGREAL) :: EM(18,18), emb(9,9)
   REAL(KIND=CGREAL) :: X1,Y1,X2,Y2,X3,Y3

   INTEGER :: inew(9)

   emb(:,:) = 0.0

!     local: u v . . . phiz ......
!      if (ilump.eq.1) then
!           contibutions to lumped mass matrix
!            roat = rho*area*th/3.0
!            rad  = sqrt(area/3.0)
!            alpha=1.0e-1
!            emb(1,1) = roat
!            emb(2,2) = roat
!            emb(3,3) = roat*rad*rad*alpha
!            emb(4,4) = roat
!            emb(5,5) = roat
!            emb(6,6) = roat*rad*rad*alpha
!            emb(7,7) = roat
!            emb(8,8) = roat
!            emb(9,9) = roat*rad*rad*alpha
!
!      elseif (ilump.eq.2) then
!           consistent mass matrix
            roat = rho*area*th/12.0
            rad  = sqrt(area/3.0)
            alpha=1.0e-1
            emb(1,1) = roat*2
            emb(1,4) = roat
            emb(1,7) = roat
            emb(4,1) = roat
            emb(4,4) = roat*2
            emb(4,7) = roat
            emb(7,1) = roat
            emb(7,4) = roat
            emb(7,7) = roat*2
            emb(2,2) = roat*2
            emb(2,5) = roat
            emb(2,8) = roat
            emb(5,2) = roat
            emb(5,5) = roat*2
            emb(5,8) = roat
            emb(8,2) = roat
            emb(8,5) = roat
            emb(8,8) = roat*2
            emb(3,3) = roat*rad*rad*alpha
            emb(3,6) = 0.0  
            emb(3,9) = 0.0  
            emb(6,3) = 0.0  
            emb(6,6) = roat*rad*rad*alpha
            emb(6,9) = 0.0  
            emb(9,3) = 0.0  
            emb(9,6) = 0.0  
            emb(9,9) = roat*rad*rad*alpha
!      endif

!     assign to full element matrix
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
            em(ii,jj) = emb(i,j)
         enddo
      enddo

   END SUBROUTINE elmmasMRT


   SUBROUTINE elmmasPLT(rho,area,th,em,x1,x2,x3,y1,y2,y3,ilump)
!  ELeMent MASs matrix for PLaTe: triangle

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   REAL(KIND=CGREAL) :: rho,area,th,x1,x2,x3,y1,y2,y3
   REAL(KIND=CGREAL) :: em(18,18), emb(9,9)
   INTEGER :: ilump

   REAL(KIND=CGREAL) :: A(9,9),CC(9,9),ca(9,9),yn(9,9)
   REAL(KIND=CGREAL) :: roat,alpha
   INTEGER :: indx(9), inew(9)
   INTEGER :: i,j

   COMMON /EQb/ DELTA 
!
   emb(:,:) = 0.0

   IF (ilump == 1) THEN
!     contibutions to lumped mass matrix
      roat = rho*area*th/3.0
      alpha=1.0e-6
      emb(1,1) = roat
      emb(2,2) = roat*alpha
      emb(3,3) = roat*alpha
      emb(4,4) = roat
      emb(5,5) = roat*alpha
      emb(6,6) = roat*alpha
      emb(7,7) = roat
      emb(8,8) = roat*alpha
      emb(9,9) = roat*alpha

   ELSEIF (ilump.eq.2) THEN

      DELTA1=(X2-X1)*(Y3-Y1)
      DELTA2=(Y2-Y1)*(X3-X1)
      DELTA=.5*(DELTA1-DELTA2)
      A1=X3-X2
      A2=X1-X3
      A3=X2-X1
      B1=Y2-Y3
      B2=Y3-Y1
      B3=Y1-Y2

!     CONSTRUCT MATRIX [A] 
      A(:,:) = 0.0

      A(1,1)=2.*DELTA 
      A(2,1)=A1 
      A(2,2)=A2 
      A(2,3)=A3 
      A(2,6)=A3 
      A(2,7)=A2 
      A(3,1)=-B1
      A(3,2)=-B2
      A(3,3)=-B3
      A(3,6)=-B3
      A(3,7)=-B2
      A(4,2)=2.*DELTA 
      A(5,1)=A1 
      A(5,2)=A2 
      A(5,3)=A3 
      A(5,4)=A1 
      A(5,8)=A3 
      A(6,1)=-B1
      A(6,2)=-B2
      A(6,3)=-B3
      A(6,4)=-B1
      A(6,8)=-B3
      A(7,3)=2.*DELTA 
      A(8,1)=A1 
      A(8,2)=A2 
      A(8,3)=A3 
      A(8,5)=A2 
      A(8,9)=A1 
      A(9,1)=-B1
      A(9,2)=-B2
      A(9,3)=-B3
      A(9,5)=-B2
      A(9,9)=-B1
      A(:,:)=A(:,:)/(2.*DELTA)

!****************************************************************** 
!    INVERT MATRIX [A]  after  [A] CONTAINS THE INVERSION 
     ia = 9
     call ainver(A,ia,indx,yn)
!****************************************************************** 
      K0 = 0
      K1 = 1
      K2 = 2
      k3 = 3
      k4 = 4
!     Lumped mass matrix
      roat = rho*th

!     Duplication
      hrs002  = hrs(k0,k0,k2,delta)
      hrs011  = hrs(k0,k1,k1,delta)
      hrs013  = hrs(k0,k1,k3,delta) 
      hrs020  = hrs(k0,k2,k0,delta)
      hrs022  = hrs(k0,k2,k2,delta) 
      hrs024  = hrs(k0,k2,k4,delta)
      hrs031  = hrs(k0,k3,k1,delta)
      hrs033  = hrs(k0,k3,k3,delta) 
      hrs042  = hrs(k0,k4,k2,delta) 
      hrs101  = hrs(k1,k0,k1,delta)
      hrs103  = hrs(k1,k0,k3,delta) 
      hrs110  = hrs(k1,k0,k1,delta)

      hrs110  = hrs(k1,k1,k0,delta)
      hrs112  = hrs(k1,k1,k2,delta)
      hrs114  = hrs(k1,k1,k4,delta) 
      hrs121  = hrs(k1,k2,k1,delta)
      hrs123  = hrs(k1,k2,k3,delta)
      hrs130  = hrs(k1,k3,k0,delta) 
      hrs132  = hrs(k1,k3,k2,delta) 
      hrs141  = hrs(k1,k4,k1,delta) 
      hrs200  = hrs(k2,k0,k0,delta)
      hrs202  = hrs(k2,k0,k2,delta) 
      hrs204  = hrs(k2,k0,k4,delta) 
      hrs211  = hrs(k2,k1,k1,delta)
      hrs213  = hrs(k2,k1,k3,delta)
      hrs220  = hrs(k2,k2,k0,delta)
      hrs222  = hrs(k2,k2,k2,delta)
      hrs231  = hrs(k2,k3,k1,delta)
      hrs240  = hrs(k2,k4,k0,delta)
      hrs301  = hrs(k3,k0,k1,delta)
      hrs303  = hrs(k3,k0,k3,delta) 
      hrs310  = hrs(k3,k1,k0,delta) 
      hrs312  = hrs(k3,k1,k2,delta)
      hrs321  = hrs(k3,k2,k1,delta) 
      hrs330  = hrs(k3,k3,k0,delta) 
      hrs402  = hrs(k4,k0,k2,delta) 
      hrs411  = hrs(k4,k1,k1,delta) 
      hrs420  = hrs(k4,k2,k0,delta) 

      cc(1,1) = hrs200       
      cc(1,2) = hrs110       
      cc(1,3) = hrs101       
      cc(1,4) = hrs220        + hrs211 /2.0
      cc(1,5) = hrs112        + hrs211 /2.0
      cc(1,6) = hrs301        + hrs211 /2.0
      cc(1,7) = hrs310        + hrs211 /2.0
      cc(1,8) = hrs121        + hrs211 /2.0
      cc(1,9) = hrs202        + hrs211 /2.0
      cc(2,2) = hrs020       
      cc(2,3) = hrs011       
      cc(2,4) = hrs130        + hrs121 /2.0
      cc(2,5) = hrs022        + hrs121 /2.0
      cc(2,6) = hrs211        + hrs121 /2.0
      cc(2,7) = hrs220        + hrs121 /2.0
      cc(2,8) = hrs031        + hrs121 /2.0
      cc(2,9) = hrs112        + hrs121 /2.0
      cc(3,3) = hrs002       
      cc(3,4) = hrs121        + hrs112 /2.0
      cc(3,5) = hrs013        + hrs112 /2.0
      cc(3,6) = hrs202        + hrs112 /2.0
      cc(3,7) = hrs211        + hrs112 /2.0
      cc(3,8) = hrs022        + hrs112 /2.0
      cc(3,9) = hrs103        + hrs112 /2.0
      cc(4,4) = hrs240 + hrs222  /4.0+   hrs231       
      cc(4,5) = hrs132 + hrs222  /4.0+  (hrs123  +hrs231)/2
      cc(4,6) = hrs321 + hrs222  /4.0+  (hrs312  +hrs231)/2
      cc(4,7) = hrs330 + hrs222  /4.0+  (hrs321  +hrs231)/2
      cc(4,8) = hrs141 + hrs222  /4.0+  (hrs132  +hrs231)/2
      cc(4,9) = hrs222 + hrs222  /4.0+  (hrs213  +hrs231)/2
      cc(5,5) = hrs024 + hrs222 /4.0+ (hrs123 +hrs123)/2
      cc(5,6) = hrs213 + hrs222 /4.0+ (hrs312 +hrs123)/2
      cc(5,7) = hrs222 + hrs222 /4.0+ (hrs321 +hrs123)/2
      cc(5,8) = hrs033 + hrs222 /4.0+ (hrs132 +hrs123)/2
      cc(5,9) = hrs114 + hrs222 /4.0+ (hrs213 +hrs123)/2
      cc(6,6) = hrs402 + hrs222 /4.0+    (hrs312 +hrs312 )/2
      cc(6,7) = hrs411 + hrs222 /4.0+    (hrs321 +hrs312 )/2
      cc(6,8) = hrs222 + hrs222 /4.0+    (hrs132 +hrs312 )/2
      cc(6,9) = hrs303 + hrs222 /4.0+    (hrs213 +hrs312 )/2
      cc(7,7) = hrs420 + hrs222 /4.0+    (hrs321 +hrs321 )/2
      cc(7,8) = hrs231 + hrs222 /4.0+    (hrs132 +hrs321 )/2
      cc(7,9) = hrs312 + hrs222 /4.0+    (hrs213 +hrs321 )/2
      cc(8,8) = hrs042 + hrs222 /4.0+    (hrs132 +hrs132 )/2
      cc(8,9) = hrs123 + hrs222 /4.0+    (hrs213 +hrs132 )/2
      cc(9,9) = hrs204 + hrs222 /4.0+    (hrs213 +hrs213 )/2

!     Impose symmetry
      do i=1,9
         do j=i+1,9
            cc(j,i) = cc(i,j)
         enddo
      enddo

      cc(:,:) = cc(:,:)*roat

!     Form product    [A]trans x ([CC] x [A])
      call DAxB (cc,9,9,A,9,9,ca)
      call DtAxB (A,9,9,ca,9,9,emb)

   ENDIF

!  Assign to full element matrix
   inew(1)=3                    !w1
   inew(2)=4                    !phix1
   inew(3)=5                    !phiy1
   inew(4)=9                    !w2
   inew(5)=10
   inew(6)=11
   inew(7)=15
   inew(8)=16
   inew(9)=17
   do i=1,9
      ii=inew(i)
      do j=1,9
         jj=inew(j)
         em(ii,jj) = emb(i,j)
      enddo
   enddo
!
   END SUBROUTINE elmmasPLT     



   FUNCTION HRS(IR,IS,IT,delta)
!****************************************************************** 
!  COMPUTE THE INTEGRAL OF (L1**IR)*(L2**IS)*(L3**IT)
!  OVER THE REGION OF EACH ELEMENT
!****************************************************************** 

   FACR=1. 
   FACS=1.0
   FACT=1.0
   FACRST=1.0

   IRST=IR+IS+IT+2.

   DO 5 I=1,IR 
 5 FACR=FACR*I 

   DO 10 J=1,IS
 10  FACS=FACS*J 

   DO 15 K=1,IT
 15       FACT=FACT*K 

   DO 20 L=1,IRST
 20       FACRST=FACRST*L 

   HRS=2.*DELTA*FACR*FACS*FACT/FACRST

   RETURN
   END 


!  Double precision [A] times [B]
   subroutine DAxB (A,nra,nca,B,nrb,ncb,C)
   real*8  A(nra,nca),B(nrb,ncb),C(nra,ncb)
   real*8  sum
!
   do i=1,nra
      do j=1,ncb
         sum=0.0
         do k=1,nca
            sum=sum+A(i,k)*B(k,j)
         enddo
         C(i,j)=sum
      enddo
   enddo

   return
   END SUBROUTINE DAxB 


!  Double precision [A]trans times [B]
   SUBROUTINE DtAxB (A,nra,nca,B,nrb,ncb,C) 
   real*8    A(nra,nca),B(nrb,ncb),C(nca,ncb)
   real*8  sum

   do i=1,nca
      do j=1,ncb
         sum=0.0
         do k=1,nra
            sum=sum+a(k,i)*B(k,j)
         enddo
         C(i,j)=sum
      enddo
   enddo

   END SUBROUTINE DtAxB 



! SUBROUTINE ELMMAS
!     This subroutine determines the mass matrix for the system.
!
      subroutine elmmasFRM(  rho, area, length, zix, em,iglobal,ilump)

         real rho, area, length,zix
         real*8   em(12,12)

         do 10 i = 1,12
            do 12 j = 1,12
               em(i,j) = 0.0
 12         continue
 10      continue

         if (ilump.ne.1) goto 30

!        contibutions to lumped mass matrix
         roal = rho*area*length/2.0
!
         em(1,1) = roal
         em(2,2) = roal
         em(3,3) = roal
         em(4,4) = roal*zix/area
         em(5,5) = roal*length*length/48
         em(6,6) = roal*length*length/48
            em(7,7)   =   em(1,1)
            em(8,8)   =   em(2,2)
            em(9,9)   =   em(3,3)
            em(10,10) =   em(4,4)
            em(11,11) =   em(5,5)
            em(12,12) =   em(6,6)
         return
!
 30   continue
!        contributions to consistent  mass matrix
         roala = rho*area*length/6.0
         roalb = rho*area*length/420.0
         roalc = rho*area*length/12.0
!
         em(1,1) =  roala*2.0
         em(1,7) =  roala

         em(2,2) =  roalb*156.0
         em(2,6) =                      + roalb*22.0*length
         em(2,8) =   roalb*54.0
         em(2,12) =                      - roalb*13.0*length

         em(3,3) =  roalb*156.0
         em(3,5) =                      - roalb*22.0*length
         em(3,9) =   roalb*54.0
         em(3,11) =                      + roalb*13.0*length

         em(4,4) =                        roala*2.0*zix/area
         em(4,10) =                       roala*zix/area

         em(5,5) =                        roalb*4.0*length*length
         em(5,9) =                       -roalb*13.0*length
         em(5,11) =                       -roalb*3.0*length*length

         em(6,6) =                        roalb*4.0*length*length
         em(6,8) =                        roalb*13.0*length
         em(6,12) =                       -roalb*3.0*length*length

       em(7,7)   =   em(1,1)
       em(8,8)   =   em(2,2)
       em(9,9)   =   em(3,3)
       em(10,10) =   em(4,4)
       em(11,11) =   em(5,5)
       em(12,12) =   em(6,6)
       em(8,12)   =  -em(2,6)
       em(9,11)   =  -em(3,5)

!           impose  symmetry
            do 20 i= 2, 12
               do 22 j= 1, i-1
                  em(i,j) = em(j,i)
 22            continue
 20         continue
!
      return
      end
! 


!     integrals of shape functions
      subroutine kgmat(dum1,dum2,dum3,dum4,dum5,dum6,dmat,delta)
         real*8 dmat(9,9)

            d1d1=dum1
            d2d2=dum2
            d3d3=dum3
            d1d2=dum4
            d1d3=dum5
            d2d3=dum6

         k0 = 0
         k1 = 1
         k2 = 2
         k3 = 3
         k4 = 4

              h14=hrs(k4,k0,k0,delta)
              h24=hrs(k0,k4,k0,delta)
              h34=hrs(k0,k0,k4,delta)

              h13h2=hrs(k3,k1,k0,delta)
              h1h23=hrs(k1,k3,k0,delta)
              h23h3=hrs(k0,k3,k1,delta)
              h2h33=hrs(k0,k1,k3,delta)
              h1h33=hrs(k1,k0,k3,delta)
              h13h3=hrs(k3,k0,k1,delta)

              h12h22=hrs(k2,k2,k0,delta)
              h22h32=hrs(k0,k2,k2,delta)
              h12h32=hrs(k2,k0,k2,delta)

              h12=hrs(k2,k0,k0,delta)
              h22=hrs(k0,k2,k0,delta)
              h32=hrs(k0,k0,k2,delta)

              h1h2=hrs(k1,k1,k0,delta)
              h2h3=hrs(k0,k1,k1,delta)
              h1h3=hrs(k1,k0,k1,delta)

              h12h2h3=hrs(k2,k1,k1,delta)
              h1h22h3=hrs(k1,k2,k1,delta)
              h1h2h32=hrs(k1,k1,k2,delta)

!..........Elements of Geometric matrix

            dmat(1,1)=d1d1*hrs(0,0,0,delta)   
            dmat(1,2)=d1d2*hrs(0,0,0,delta)
            dmat(1,3)=d1d3*hrs(0,0,0,delta)
            dmat(1,4)=(h22+0.5*h2h3)*d1d1+(2.*h1h2+0.5*h1h3)*d1d2+ (0.5*h1h2)*d1d3
            dmat(1,5)=(0.5*h2h3)*d1d1+(h32+0.5*h1h3)*d1d2+(2.0*h2h3+0.5*h1h2)*d1d3
            dmat(1,6)=(2.0*h1h3+0.5*h2h3)*d1d1+(0.5*h1h3)*d1d2+ &
                  (h12+0.5*h1h2)*d1d3
            dmat(1,7)=(2.0*h1h2+0.5*h2h3)*d1d1+(h12+0.5*h1h3)*d1d2+ &
                   (0.5*h1h2)*d1d3
            dmat(1,8)=(0.5*h2h3)*d1d1+(2.0*h2h3+0.5*h1h3)*d1d2+  &
                  (h22+0.5*h1h2)*d1d3
            dmat(1,9)=(h32+0.5*h2h3)*d1d1+(0.5*h1h3)*d1d2+ &
                   (2.0*h1h3+0.5*h1h2)*d1d3
!
            dmat(2,2)=d2d2*hrs(0,0,0,delta)
            dmat(2,3)=d2d3*hrs(0,0,0,delta)
            dmat(2,4)=(h22+0.5*h2h3)*d1d2+(2.*h1h2+0.5*h1h3)*d2d2+ &
                   (0.5*h1h2)*d2d3
            dmat(2,5)=(0.5*h2h3)*d1d2+(h32+0.5*h1h3)*d2d2+   &
                   (2.0*h2h3+0.5*h1h2)*d2d3
            dmat(2,6)=(2.0*h1h3+0.5*h2h3)*d1d2+(0.5*h1h3)*d2d2+  &
                   (h12+0.5*h1h2)*d2d3
            dmat(2,7)=(2.0*h1h2+0.5*h2h3)*d1d2+(h12+0.5*h1h3)*d2d2+ &
                   (0.5*h1h2)*d2d3
            dmat(2,8)=(0.5*h2h3)*d1d2+(2.0*h2h3+0.5*h1h3)*d2d2+  &
                   (h22+0.5*h1h2)*d2d3
            dmat(2,9)=(h32+0.5*h2h3)*d1d2+(0.5*h1h3)*d2d2+  &
                   (2.0*h1h3+0.5*h1h2)*d2d3

            dmat(3,3)=d3d3*hrs(0,0,0,delta)
            dmat(3,4)=(h22+0.5*h2h3)*d1d3+(2.*h1h2+0.5*h1h3)*d2d3+  &
                   (0.5*h1h2)*d3d3
            dmat(3,5)=(0.5*h2h3)*d1d3+(h32+0.5*h1h3)*d2d3+  &
                   (2.0*h2h3+0.5*h1h2)*d3d3
            dmat(3,6)=(2.0*h1h3+0.5*h2h3)*d1d3+(0.5*h1h3)*d2d3+  &
                   (h12+0.5*h1h2)*d3d3
            dmat(3,7)=(2.0*h1h2+0.5*h2h3)*d1d3+(h12+0.5*h1h3)*d2d3+  &
                   (0.5*h1h2)*d3d3
            dmat(3,8)=(0.5*h2h3)*d1d3+(2.0*h2h3+0.5*h1h3)*d2d3+  &
                   (h22+0.5*h1h2)*d3d3
            dmat(3,9)=(h32+0.5*h2h3)*d1d3+(0.5*h1h3)*d2d3+  &
                   (2.0*h1h3+0.5*h1h2)*d3d3

          dmat(4,4)=d1d1*h14+d1d1*h23h3+(4.*d1d2+d1d3)*h1h23+  &
                   (4.0*d2d2+0.25*d3d3+2.0*d2d3)*h12h22+  &
                   (0.25*d1d1)*h22h32+(0.25*d2d2)*h12h32+  &
                   (2.*d2d2+0.5*d2d3)*h12h2h3+  &
                   (3.*d1d2+0.5*d1d3)*h1h22h3+(0.5*d1d2)*h1h2h32
            dmat(4,5)=(0.5*d1d3)*h1h23+(0.5*d1d1+2.*d1d3)*h23h3+  &
                  (0.5*d2d2)*h1h33+(0.5*d1d2)*h2h33+  &
                 (0.25*d3d3+d2d3)*h12h22+(0.25*d1d1+d1d2+d1d3)*h22h32+  &
                  (0.25*d2d2)*h12h32+(d2d2+0.5*d2d3)*h12h2h3+ &
                  (d3d3+1.5*d1d2+4.*d2d3+0.5*d1d3)*h1h22h3+ &
                  (2.*d2d2+0.5*d1d2+1.5*d2d3)*h1h2h32
            dmat(4,6)=(0.5*d2d3)*h13h3+(0.5*d1d1)*h23h3+  &
                   (0.5*d3d3+d2d3)*h13h2+(0.5*d1d3)*h1h23+  &
                   (0.25*d3d3+d2d3+d1d3)*h12h22+(0.25*d1d1)*h22h32+  &
                   (0.25*d2d2+d1d2)*h12h32+  &
                   (d2d2+4.*d1d2+0.5*d2d3+1.5*d1d3)*h12h2h3+  &
                   (2.*d1d1+1.5*d1d2+0.5*d1d3)*h1h22h3+  &
                   (d1d1+0.5*d1d2)*h1h2h32  
          dmat(4,7)=(2.*d1d1+0.5*d1d3)*h1h23+(0.5*d1d1)*h23h3+  &
                   (2.*d2d2+0.5*d2d3)*h13h2+(0.5*d2d2)*h13h3+  &
                   (0.25*d3d3+5.*d1d2+d2d3+d1d3)*h12h22+  &
                  (0.25*d1d1)*h22h32+(0.25*d2d2)*h12h32+  &
                   (d2d2+1.5*d1d2+0.5*d2d3)*h12h2h3+   &
                   (d1d1+1.5*d1d2+0.5*d1d3)*h1h22h3+(0.5*d1d2)*h1h2h32
          dmat(4,8)=d1d3*h24+(0.5*d1d1+2.*d1d2+0.5*d1d3)*h23h3+   &
                  (0.5*d3d3+2.*d2d3+0.5*d1d3)*h1h23+          &
                  (0.25*d3d3+d2d3)*h12h22+(0.25*d1d1+d1d2)*h22h32+  &
                   (0.25*d2d2)*h12h32+(d2d2+0.5*d2d3)*h12h2h3+  &
                   (4.*d2d2+1.5*d1d2+1.5*d2d3+0.5*d1d3)*h1h22h3+  &
                   (d2d2+0.5*d1d2)*h1h2h32
           dmat(4,9)=(0.5*d1d1)*h2h33+(0.5*d1d1)*h23h3+(0.5*d1d2)*h1h33+ &
                  (0.5*d1d3)*h1h23+(0.25*d3d3+d2d3)*h12h22+  &
                   (1.25*d1d1)*h22h32+(d2d3+0.25*d2d2)*h12h32+   &
                   (d3d3+4.5*d2d3+d2d2)*h12h2h3+  &
                   (1.5*d1d2+2.5*d1d3)*h1h22h3+   &
                   (2.5*d1d2+1.5*d1d3)*h1h2h32    

            dmat(5,5)=d2d2*h34+d2d2*h1h33+(d1d2+4.*d2d3)*h2h33+  &
                   (0.25*d3d3)*h12h22+    &
                   (0.25*d1d1+4.*d3d3+2.*d1d3)*h22h32+  &
                   (0.25*d2d2)*h12h32+(0.5*d2d3)*h12h2h3+ &
                   (2.*d3d3+0.5*d1d3)*h1h22h3+   &
                   (0.5*d1d2+3.*d2d3)*h1h2h32
            dmat(5,6)=(0.5*d2d2+2.*d1d2)*h1h33+(0.5*d3d3)*h13h2+  &
                   (0.5*d2d3)*h13h3+(0.5*d1d2)*h2h33+   &
                  (0.25*d3d3)*h12h22+(0.25*d1d1+d1d3)*h22h32+   &
                   (0.25*d2d2+d1d2+d2d3)*h12h32+   &
                   (2.*d3d3+0.5*d2d3+1.5*d1d3)*h12h2h3+  &
                   (d3d3+0.5*d1d3)*h1h22h3+  &
                   (d1d1+0.5*d1d2+1.5*d2d3+4.0*d1d3)*h1h2h32

           dmat(5,7)=(0.5*d2d2)*h13h3+(0.5*d2d2)*h1h33+(0.5*d1d2)*h2h33+ &
                  (0.5*d2d3)*h13h2+(0.25*d3d3+d1d3)*h12h22+  &
                   (0.25*d1d1+d1d3)*h22h32+(1.25*d2d2)*h12h32+  &
                   (1.5*d1d2+2.5*d2d3)*h12h2h3+  &
                   (d1d1+d3d3+4.5*d1d3)*h1h22h3+  &
                   (2.5*d1d2+1.5*d2d3)*h1h2h32
            dmat(5,8)=(2.*d2d2+0.5*d1d2)*h2h33+(0.5*d2d2)*h1h33+  &
                  (2.*d3d3+0.5*d1d3)*h23h3+(0.5*d3d3)*h1h23+  &
                   (0.25*d3d3)*h12h22+  &
                   (0.25*d1d1+d1d2+5.*d2d3+d1d3)*h22h32+  &
                   (0.25*d2d2)*h12h32+(0.5*d2d3)*h12h2h3+  &
                   (d3d3+1.5*d2d3+0.5*d1d3)*h1h22h3+  &
                   (d2d2+0.5*d1d2+1.5*d2d3)*h1h2h32
            dmat(5,9)=d1d2*h34+(0.5*d1d1+0.5*d1d2+2.*d1d3)*h2h33+  &
                   (0.5*d2d2+0.5*d1d2+2.*d2d3)*h1h33+  &
                   (0.25*d3d3)*h12h22+(0.25*d1d1+d1d3)*h22h32+   &
                   (0.25*d2d2+d2d3)*h12h32+(d3d3+0.5*d2d3)*h12h2h3+  &
                   (d3d3+0.5*d1d3)*h1h22h3+  &
                   (4.*d3d3+0.5*d1d2+1.5*d2d3+1.5*d1d3)*h1h2h32

            dmat(6,6)=d3d3*h14+d3d3*h13h2+(d2d3+4.*d1d3)*h13h3+  &
                   (0.25*d3d3)*h12h22+(0.25*d1d1)*h22h32+  &
                   (0.25*d2d2+4.*d1d1+2.*d1d2)*h12h32+  &
                   (0.5*d2d3+3.*d1d3)*h12h2h3+(0.5*d1d3)*h1h22h3+   &
                   (2.*d1d1+0.5*d1d2)*h1h2h32
            dmat(6,7)=d2d3*h14+(0.5*d2d2+0.5*d2d3+2.*d1d2)*h13h3+  &
                   (0.5*d3d3+0.5*d2d3+2.*d1d3)*h13h2+  &
                   (0.25*d3d3+d1d3)*h12h22+(0.25*d1d1)*h22h32+  &
                   (0.25*d2d2+d1d2)*h12h32+  &
                  (4.*d1d1+1.5*d1d2+0.5*d2d3+1.5*d1d3)*h12h2h3+  &
                   (d1d1+0.5*d1d3)*h1h22h3+(d1d1+0.5*d1d2)*h1h2h32
           dmat(6,8)=(0.5*d3d3)*h1h23+(0.5*d3d3)*h13h2+(0.5*d2d3)*h13h3+  &
                    (0.5*d1d3)*h23h3+(1.25*d3d3)*h12h22+  &
                   (0.25*d1d1+d1d2)*h22h32+(0.25*d2d2+d1d2)*h12h32+  &
                    (2.5*d2d3+1.5*d1d3)*h12h2h3+  &
                    (1.5*d2d3+2.5*d1d3)*h1h22h3+   &
                    (d1d1+d2d2+4.5*d1d2)*h1h2h32
          dmat(6,9)=(2.*d1d1+0.5*d1d2)*h1h33+(0.5*d1d1)*h2h33+  &
                   (2.*d3d3+0.5*d2d3)*h13h3+(0.5*d3d3)*h13h2+  &
                   (0.25*d3d3)*h12h22+(0.25*d1d1)*h22h32+  &
                   (0.25*d2d2+d1d2+d2d3+5.*d1d3)*h12h32+  &
                   (d3d3+0.5*d2d3+1.5*d1d3)*h12h2h3+ &
                   (0.5*d1d3)*h1h22h3+(d1d1+0.5*d1d2+1.5*d1d3)*h1h2h32

           dmat(7,7)=d2d2*h14+d2d2*h13h3+(4.*d1d2+d2d3)*h13h2+  &
                   (4.*d1d1+0.25*d3d3+2.*d1d3)*h12h22+  &
                   (0.25*d1d1)*h22h32+(0.25*d2d2)*h12h32+  &
                   (3.*d1d2+0.5*d2d3)*h12h2h3+  &
                   (2.*d1d1+0.5*d1d3)*h1h22h3+(0.5*d1d2)*h1h2h32
          dmat(7,8)=(0.5*d2d2)*h13h3+(0.5*d3d3+2.*d1d3)*h1h23+  &
                  (0.5*d2d3)*h13h2+(0.5*d1d3)*h23h3+  &
                 (0.25*d3d3+d2d3+d1d3)*h12h22+(0.25*d1d1+d1d2)*h22h32+  &
                  (0.25*d2d2)*h12h32+  &
                  (2.*d2d2+1.5*d1d2+0.5*d2d3)*h12h2h3+  &
                  (d1d1+4.*d1d2+1.5*d2d3+0.5*d1d3)*h1h22h3+   &
                  (d2d2+0.5*d1d2)*h1h2h32
          dmat(7,9)=(0.5*d1d1)*h2h33+(0.5*d2d2+2.*d2d3)*h13h3+  &
                   (0.5*d2d3)*h13h2+(0.25*d3d3+d1d3)*h12h22+ &
                   (0.25*d1d1)*h22h32+ &
                   (0.25*d2d2+d1d2+d2d3)*h12h32+  &
                   (d3d3+1.5*d1d2+0.5*d2d3+4.*d1d3)*h12h2h3+  &
                   (d1d1+0.5*d1d3)*h1h22h3+  &
                  (2.*d1d1+0.5*d1d2+1.5*d1d3)*h1h2h32

           dmat(8,8)=d3d3*h24+d3d3*h1h23+(4.*d2d3+d1d3)*h23h3+  &
                 (0.25*d3d3)*h12h22+  &
                  (0.25*d1d1+4.*d2d2+2.*d1d2)*h22h32+ &
                 (0.25*d2d2)*h12h32+(0.5*d2d3)*h12h2h3+ &
                  (3.*d2d3+0.5*d1d3)*h1h22h3+ &
                  (2.*d2d2+0.5*d1d2)*h1h2h32
         dmat(8,9)=(0.5*d1d1+2.*d1d2)*h2h33+(0.5*d3d3)*h1h23+  &
                  (0.5*d1d2)*h1h33+(0.5*d1d3)*h23h3+  &
                  (0.25*d3d3)*h12h22+  &
                  (0.25*d1d1+d1d2+d1d3)*h22h32+  &
                  (0.25*d2d2+d2d3)*h12h32+(d3d3+0.5*d2d3)*h12h2h3+  &
                  (2.*d3d3+1.5*d2d3+0.5*d1d3)*h1h22h3+  &
                  (d2d2+0.5*d1d2+4.*d2d3+1.5*d1d3)*h1h2h32

           dmat(9,9)=d1d1*h34+d1d1*h2h33+(d1d2+4.*d1d3)*h1h33+  &
                  (0.25*d3d3)*h12h22+(0.25*d1d1)*h22h32+ &
                  (0.25*d2d2+2.*d2d3+4.*d3d3)*h12h32+ &
                (2.*d3d3+0.5*d2d3)*h12h2h3+(0.5*d1d3)*h1h22h3+  &
                (0.5*d1d2+3.*d1d3)*h1h2h32

!           impose symmetry
            do 10 i=1,9
               do 10 j=1,9
                  dmat(j,i)=dmat(i,j)
10          continue
      return
      end


   SUBROUTINE ELM_ASSIGN(DISPFUL,IGLOBAL,nmprop,ippp,iBody)
!  ELeMement ASSIGNment  to get stresses etc

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface
   USE flow_parameters, ONLY : nPtsBodyMarker
   REAL (KIND=CGREAL):: dispful(nPtsBodyMarker(iBody),6)
   REAL (KIND=CGREAL):: strain(18),stress(18), force(18)
   REAL (KIND=CGREAL):: disploc(3,6), temp(6)
   REAL (KIND=CGREAL):: delt(12),nload(12)
   REAL (KIND=CGREAL):: znu
   REAL(KIND=CGREAL) :: xll,xmm,xnn  !Input
   REAL(KIND=CGREAL) :: dxb,dyb,dzb
   REAL(KIND=CGREAL) :: x1,x2,y1,y2,z1,z2
   
   INTEGER :: i1,j1,k1

!  integer nmprop(nel)

!  For each element, calculate the strain, stress at centroid
!   write(*,*)'@@ loop over elements ',nel
!   write(*,*)' '

   do 50 i=1,nel
      i1 = triElemNeig(iBody,1,i)
      j1 = triElemNeig(iBody,2,i)
      k1 = triElemNeig(iBody,3,i)

      IF (ELETYPE .EQ. 2) THEN
!        Frame
!         e0 = PROPERTY(ELMatType(i),1)
!         g0 = PROPERTY(ELMatType(i),2)
!         a0 = PROPERTY(ELMatType(i),3)
!         r0 = PROPERTY(ELMatType(i),4)
!         b0 = PROPERTY(ELMatType(i),5)
!        Element type is set to 1 temporarily \\Wanh
         e0 = PROPERTY(1,1)
         g0 = PROPERTY(1,2)
         a0 = PROPERTY(1,3)
         r0 = PROPERTY(1,4)
         b0 = PROPERTY(1,5)

!         zix0 = PROPERTY(ELMatType(i),6)
!         ziy0 = PROPERTY(ELMatType(i),7)
!         ziz0 = PROPERTY(ELMatType(i),8)
!
         dx = xBodyMarker(iBody,j1)-xBodyMarker(iBody,i1)
         dy = yBodyMarker(iBody,j1)-yBodyMarker(iBody,i1)
         dz = zBodyMarker(iBody,j1)-zBodyMarker(iBody,i1)
         xl=sqrt(dx*dx+dy*dy+dz*dz) 
         xll=dx/xl 
         xmm=dy/xl 
         xnn=dz/xl 

!        obtain local disps
         do kk=1,6
            delt(kk  )=dispful(i1,kk)
            delt(kk+6)=dispful(j1,kk)
         enddo
!        The following 9 lines are commented out for later development //Wanh
!         call trns3dv(delt,xll,xmm,xnn,nload,b0) 
!         do kk=1,6
!            disploc(1,kk) = nload(kk)
!            disploc(2,kk) = nload(kk+6)
!         enddo
!
!         call convertFRM(disploc,xl,  &
!                                  zix0,ziy0,ziz0,a0,e0,g0,  &
!                                  strain,stress,force,iglobal )

      ELSEIF (ELETYPE .EQ. 1) THEN
!        Plate
!         e0 = PROPERTY(ELMatType(i),1)
!         g0 = PROPERTY(ELMatType(i),2)
!         t0 = PROPERTY(ELMatType(i),3)
!         r0 = PROPERTY(ELMatType(i),4)
!         pL0 = PROPERTY(ELMatType(i),5)
!        Element type is set to 1 temporarily \\Wanh
         e0 = PROPERTY(1,1)
         g0 = PROPERTY(1,2)
         t0 = PROPERTY(1,3)
         r0 = PROPERTY(1,4)
         pL0 = PROPERTY(1,5)
         
!         zip0 = PROPERTY(ELMatType(i),6)
         znu = e0/(2.0*g0) - 1.0
         zip0 = e0*t0*t0*t0/(1-znu*znu)/12.0  !uniform plate
!         zia0 = PROPERTY(ELMatType(i),7)
!         zib0 = PROPERTY(ELMatType(i),8)
         zia0 = 1.5
         zib0 = 0.5

         x1 = xBodyMarker(iBody,i1)
         x2 = xBodyMarker(iBody,j1)
         x3 = xBodyMarker(iBody,k1)

         y1 = yBodyMarker(iBody,i1)
         y2 = yBodyMarker(iBody,j1)
         y3 = yBodyMarker(iBody,k1)

         z1 = zBodyMarker(iBody,i1)
         z2 = zBodyMarker(iBody,j1)
         z3 = zBodyMarker(iBody,k1)

!        determine vector area
         axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
         ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
         azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
         area=sqrt( axy*axy + ayz*ayz + azx*azx)
         xll=ayz/area
         xmm=azx/area
         xnn=axy/area

!        transform element coords to local X-Y
         xb1=x1
         yb1=y1
         zb1=z1
         call rotvec(xll,xmm,xnn,(x2-x1),(y2-y1),(z2-z1),dxb,dyb,dzb)
         xb2=x1+dxb
         yb2=y1+dyb
         zb2=z1+dzb
         call rotvec(xll,xmm,xnn,(x3-x1),(y3-y1),(z3-z1),dxb,dyb,dzb) 
         xb3=x1+dxb
         yb3=y1+dyb
         zb3=z1+dzb

!        rotate displacements to local coords
         nnp = nPtsBodyMarker(iBody)
         call ROTDISP(xll,xmm,xnn,dispful,disploc,nnp,i1,j1,k1,iBody)

!        we only need the in-plane stresses
         alpha=zia0
         beta =zib0

!        The following 4 lines are commented out for later development //Wanh         
!         call convertMRT(disploc, area, &
!                         e0,t0,g0,pl0,alpha,beta, &
!                         xb1,xb2,xb3,yb1,yb2,yb3, &
!                         strain, stress,force )
      endif

!     Bottom distinction between elements
!     SAVE centroidal values for stability analysis
      IF (ELETYPE .EQ. 2) THEN
!        for frame
         temp(1) = stress(1)*a0
         temp(2) = stress(2)*a0
         temp(3) = stress(3)*a0
      ELSEIF (ELETYPE .EQ. 1) THEN
         temp(1) = (stress(1)+stress(7)+stress(13))/3.0
         temp(2) = (stress(2)+stress(8)+stress(14))/3.0
         temp(3) = (stress(3)+stress(9)+stress(15))/3.0
      ENDIF
      write(igeo) (temp(j), j=1,3)

 50 continue !  End of loop over elements


   END SUBROUTINE ELM_ASSIGN



! SUBROUTINE ELMgeom
!     This subroutine calculates the element stiffness matrices.
!
       subroutine elmgeomFRM( length, eg,s,iglobal   )
          real     length
          real*8   eg(12,12)

!         initialize all eg elements to zero
          do 90 i=1,12
             do 90 j=1,12
90              eg(i,j)=0.0

!         Stiffness matrix in local coordinates

          emlenz = s/(30.0*length)
          alpha=s*1.0e-3
          beta = 1.0


          if (iglobal .eq. 11 .OR. iglobal .eq. 21  &
                             .OR. iglobal .eq. 31) then
              eg(1,1)   =   alpha
              eg(2,2)   =   36*emlenz 
              eg(3,3)   =   36*emlenz   
              eg(4,4)   =   alpha
              eg(5,5)   =   0.0
              eg(6,6)   =   0.0

          else
              eg(1,1)   =   alpha
              eg(2,2)   =   36*emlenz 
              eg(3,3)   =   36*emlenz   
              eg(4,4)   =   alpha
              eg(5,5)   =   4.0*emlenz*length*length
              eg(6,6)   =   4.0*emlenz*length*length

              eg(2,6)   =   3.0*emlenz*length
              eg(3,5)   =  -3.0*emlenz*length

          endif

          eg(7,7)   =   eg(1,1)
          eg(8,8)   =   eg(2,2)
          eg(9,9)   =   eg(3,3)
          eg(10,10) =   eg(4,4)
          eg(11,11) =   eg(5,5)
          eg(12,12) =   eg(6,6)

          eg(1,7)   =   -eg(1,1)
          eg(2,8)   =   -eg(2,2)
          eg(2,12)  =    eg(2,6)
          eg(3,9)   =   -eg(3,3)
          eg(3,11)  =    eg(3,5)
          eg(4,10)  =   -eg(4,4)
          eg(5,9)   =   -eg(3,5)
          eg(5,11)  =   -eg(5,5)/4.0
          eg(6,8)   =   -eg(2,6)
          eg(6,12)  =   -eg(6,6)/4.0

          eg(8,12)  =   -eg(2,6)
          eg(9,11)  =   -eg(3,5)

!         impose the symmetry

          do 10 i= 1, 12
             do 10 j= i, 12
 10             eg(j,i) = eg(i,j)

!         check diagonal terms
          do i=1,12
             if ( abs(eg(i,i)) .lt. 1.0e-12) eg(i,i)=1.0e-12
          enddo

      return 
      end



   SUBROUTINE ROTDISP(xll,xmm,xnn,dispful,disploc,nnp,i1,j1,k1,iBody)
!  ROTate DISPlacements to local coords
   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface
   USE flow_parameters, ONLY : nPtsBodyMarker
   
   REAL(KIND=CGREAL) :: xll,xmm,xnn
   REAL(KIND=CGREAL) :: dispful(nPtsBodyMarker(iBody),6)
   REAL (KIND=CGREAL):: disploc(3,6)

   INTEGER :: i,i1,j1,k1,nnp

   REAL (KIND=CGREAL):: temp(6) 


   do kk=1,6
      temp(kk)=dispful(i1,kk)
   enddo

   call rotvec(xll,xmm,xnn,temp(1),temp(2),temp(3),   &
               disploc(1,1),disploc(1,2),disploc(1,3))
   call rotvec(xll,xmm,xnn,temp(4),temp(5),temp(6),  &
               disploc(1,4),disploc(1,5),disploc(1,6))

   do kk=1,6
      temp(kk)=dispful(j1,kk)
   enddo

   call rotvec(xll,xmm,xnn,temp(1),temp(2),temp(3),    &
               disploc( 2,1),disploc( 2,2),disploc( 2,3))
   call rotvec(xll,xmm,xnn,temp(4), temp(5), temp(6),  &
               disploc( 2,4),disploc( 2,5),disploc( 2,6))

   do kk=1,6
      temp(kk)=dispful(k1,kk)
   enddo

   call rotvec(xll,xmm,xnn,temp(1),temp(2),temp(3),    &
               disploc( 3,1),disploc( 3,2),disploc( 3,3) )
   call rotvec(xll,xmm,xnn,temp(4), temp(5), temp(6),  &
               disploc( 3,4),disploc( 3,5),disploc( 3,6) )

   END SUBROUTINE ROTDISP

