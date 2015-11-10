   SUBROUTINE uduCOL(a,imult)
!  UDU decomposition of COLumn profiled system

   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: a(MAXSTIFF)
   INTEGER :: imult

   REAL(KIND=CGREAL) :: temp,sum
   INTEGER :: i,j,k,im1,jm1,j2,j3,is,io,jcol
   INTEGER :: iloc,iloc1,iloci,ilocj

   if (a(1) .eq. 0.0d0) then
      imult = 0
      return
   endif

   do 10 j=2,NEQ
      jm1=j-1
      j2=j-iprof(j)+1
      if (j2.lt.1) then
         j2=1
      endif

!     off-diagonal terms
      if (jm1.eq.1) then
         is=j
         io=1
         iloc = NLOC(is) + io - 1
         sum=a(iloc)
!        sum=a(j,1)
      else
         do 20 i=j2+1,jm1
               im1=i-1
               is=j
               io=j-i+1
               iloc = NLOC(is) + io - 1
               sum=a(iloc   )
!              sum=a(i,j-i+1)
!
               j3=i-IPROF(i)+1
               jcol=j3
               if (j3 .lt. j2) then
                  jcol=j2
               endif
!              do 21 k=j2,im1
               do 21 k=jcol,im1
                  is=i
                  io=i-k+1
                  iloci = NLOC(is) + io - 1
                  is=j
                  io=j-k+1
                  ilocj = NLOC(is) + io - 1
                  sum=sum-a(iloci  )*a(ilocj  )
!                 sum=sum-a(k,i-k+1)*a(k,j-k+1)
                  imult=imult+1
 21            continue
               a(iloc   )=sum
!              a(i,j-i+1)=sum
 20      continue
         is=j
         io=1
         iloc = NLOC(is) + io - 1
         sum=a(iloc   )
!        sum=a(j,1)
      endif

!     diagonal terms
      do 30 k=j2,jm1
         is=j
         io=j-k+1
         ilocj = NLOC(is) + io - 1
         is=k
         io=1
         iloc1 = NLOC(is) + io - 1
         temp=a(ilocj  )/a(iloc1)
         sum=sum-temp*a(ilocj  )
         a(ilocj  )=temp
!        temp=a(k,j-k+1)/a(k,1)
!        sum=sum-temp*a(k,j-k+1)
!        a(k,j-k+1)=temp
         imult=imult+2
 30   continue
      if (sum.eq.0.0d0) then
         imult=0
         return
      endif

      is = j
      io = 1
      iloc = NLOC(is) + io - 1
      a(iloc) = sum
!     a(j,1) = sum
 10  continue

   END SUBROUTINE uduCOL


   SUBROUTINE bakCOL(a,maxstore,b,imult)
!  BAcK solver of COLumn profiled system, return wk
   USE global_parameters
   USE boundary_arrays
   USE unstructured_surface_arrays
   USE fea_unstructure_surface


!   REAL(KIND=CGREAL) :: a(maxstore),b(nPtsMax*nodeDoF)
   REAL(KIND=CGREAL) :: a(maxstore),b(NEQ)
   
   INTEGER :: maxstore,imult
   REAL(KIND=CGREAL) :: sum

!  forward substitutions
   DO 10 i=1,NEQ
!     j=i-iband+1
      j=i-iprof(i)+1
      if (i.le.iprof(i) ) then
         j=1
      endif
      jb=i-iband+1
      jbb=jb
      if (i.le.iband    ) then
         jbb=1
      endif
      sum=b(i)
      km1=i-1
      if (j.gt.km1) then
         wk(i)=sum
      else
         do 11 k=j,km1
            is=i
            io=i-k+1
            iloc = nloc(is) + io - 1
            sum=sum-a(iloc   )*wk(k)
!           sum=sum-a(k,i-k+1)*wk(k)
            imult=imult+1
 11      continue
         wk(i)=sum
      endif
 10 ENDDO

!  middle terms
   DO I=1,NEQ
      is=i
      io=1
      iloc = nloc(is) + io - 1
      wk(i)=wk(i)/a(iloc )
!     wk(i)=wk(i)/a(i,1)
      imult=imult+1
   ENDDO

!  backward substitution
   do 50 i1=1,neq
      i=neq-i1+1
      j=i+iprof2(i) -1
      if (j.gt.neq) then
         j=neq
      endif
      jb=i+iband-1
      jbb=jb
      if (jb.gt.neq) then
         jbb=neq
      endif
      sum=wk(i)
      k2=i+1
      if (k2.gt.j) then
         wk(i)=sum
      else
         do 40 k=k2,j
            is=k
            io=k-i+1
            if (io .gt. iprof(is)) goto 40
            iloc = nloc(is) + io - 1
            sum=sum-a(iloc   )*wk(k)
!           sum=sum-a(i,k-i+1)*wk(k)
            imult=imult+1
 40      continue
         wk(i)=sum
      endif
 50 continue

   END SUBROUTINE bakCOL


!   SUBROUTINE rdCOL(aa,iprof,nloc)
!!  ReaD COLumn oriented storage
!   include 'commons.std'
!
!   integer  iprof(neq), nloc(neq)
!   real*8   aa(neq*iband)

!   write(*,*)'@@ reading <<stadyn.stf<<'
!   rewind(istf)
!   read(istf) j1,j2  
!   write(ilog,*)'@@ neq iband ',j1,j2  
!   write(iout,*)'@@ reading COL form'  
!   write(iout,*)'@@ neq iband ',j1,j2  
!   do 8 i=1,neq  
!               iloc = nloc(i)
!               imax = iprof(i)
!               read (istf   ) ( aa(iloc + j-1 ), j=1,imax )
!               write(iout,22) ( aa(iloc + j-1 ), j=1,imax )
! 22            format(1x,6(g13.6))
! 8 continue

!   END SUBROUTINE rdCOL


!     UDU solution of banded system using UDU decomposition
      subroutine udu(a,neq,iband,imult)
          real*8 a(neq,iband), temp,sum

          if (a(1,1) .eq. 0.0  ) then
              imult=0
              return
          endif
          if (neq .eq. 1) then
              imult=1
              return
          endif

          do 10 j=2,neq
             jm1=j-1
             j2=j-iband+1
             if (j2.lt.1) then
                 j2=1
             endif

!            off-diagonal terms
             if (jm1.eq.1) then
                 sum=a(j,1)
             else
                do 20 i=j2+1,jm1
                   im1=i-1
                   sum=a(i,j-i+1)
                   do 21 k=j2,im1
                      sum=sum-a(k,i-k+1)*a(k,j-k+1)
                      imult=imult+1
 21                continue
                   a(i,j-i+1)=sum
 20             continue
                sum=a(j,1)
             endif

!            diagonal terms
             do 30 k=j2,jm1
                temp=a(k,j-k+1)/a(k,1)
                sum=sum-temp*a(k,j-k+1)
                a(k,j-k+1)=temp
                   imult=imult+2
 30          continue
             if (sum.eq.0.0  ) then
                 imult=0
                 return
             endif
             a(j,1)=sum
 10        continue

      return
      end

      subroutine bak(a,b,neq,iband,wk,imult)
          real*8 a(neq,iband),b(neq),wk(neq), sum

!        forward substitutions
          do 10 i=1,neq
             j=i-iband+1
             if (i.le.iband) then
                 j=1
             endif
             sum=b(i)
             km1=i-1
             if (j.gt.km1) then
                wk(i)=sum
             else
                do 11 k=j,km1
                   sum=sum-a(k,i-k+1)*wk(k)
                   imult=imult+1
 11             continue
                wk(i)=sum
             endif
 10       continue

!         middle terms
          do 30 i=1,neq
             wk(i)=wk(i)/a(i,1)
                   imult=imult+1
 30       continue

!         backward substitution
          do 50 i1=1,neq
             i=neq-i1+1
             j=i+iband-1
             if (j.gt.neq) then
                 j=neq
             endif
             sum=wk(i)
             k2=i+1
             if (k2.gt.j) then
                wk(i)=sum
             else
                do 40 k=k2,j
                   sum=sum-a(i,k-i+1)*wk(k)
                   imult=imult+1
 40             continue
                wk(i)=sum
             endif
 50        continue

      end

! ------------------------------------------------------------------- !

       SUBROUTINE AINVER(a,n,indx,yn)
       USE global_parameters

       REAL(KIND=CGREAL) :: a(n,n), yn(n,n)
       INTEGER :: n,indx(n)

       do 20  i = 1,n
          do j = 1,n
             yn(i,j) = 0.0
          enddo
          yn(i,i) = 1.0
 20    continue

       np=n
       call ludcmp(a,n,np,indx,d)

       do 30  j = 1,n
          call lubksb(a,n,np,indx,yn(1,j))
 30    continue

       do 40  j = 1,n
          do 40  i = 1,n
             a(i,j) = yn(i,j)
 40    continue

       END SUBROUTINE AINVER




   SUBROUTINE ludcmp(a,n,np,indx,d)

   USE global_parameters

   parameter (nmax=100,tiny=1.0e-20)
   REAL(KIND=CGREAL) :: a(np,np), vv(nmax)
   INTEGER :: indx(n),np

      d = 1.0
      do 12  i = 1,n
         aamax = 0.0
         do j = 1,n
            if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
         enddo

         if (aamax .eq. 0.0) then !pause 'Singular Matrix'
            write (*,*) 'Singular Matrix'
         endif

         vv(i) = 1.0/aamax
 12   continue

      do 19  j = 1,n
         do 14  i = 1,j-1
            sum = a(i,j)
            do 13  k = 1,i-1
               sum = sum - a(i,k)*a(k,j)
 13         continue
            a(i,j) = sum
 14      continue
         aamax = 0.0
         do 16  i = j,n
            sum = a(i,j)
            do 15  k = 1,j-1
               sum = sum - a(i,k)*a(k,j)
 15         continue
            a(i,j) = sum
            dum = vv(i)*abs(sum)
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            endif
 16      continue
         if (j .ne. imax) then
            do 17  k = 1,n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
 17         continue
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j) .eq. 0.0) a(j,j) = tiny
         if (j .ne. n) then 
            dum = 1.0/a(j,j)
            do 18  i = j+1,n
               a(i,j) = a(i,j)*dum
 18         continue
         endif
 19   continue

      return
   END SUBROUTINE ludcmp

   SUBROUTINE lubksb(a,n,np,indx,b)

   USE global_parameters

   REAL(KIND=CGREAL) :: a(np,np), b(n)
   INTEGER :: indx(n)

      ii = 0
      do 12  i = 1,n
         ll = indx(i)
         sum = b(ll)
         b(ll) = b(i)
         if (ii .ne. 0) then
            do 11  j = ii , i-1
               sum = sum - a(i,j)*b(j)
 11         continue
         else if (sum .ne. 0.0) then
            ii = i
         endif
         b(i) = sum
 12   continue
      do 14 i = n,1,-1
         sum = b(i)
         if (i .lt. n) then
            do 13  j = i+1,n
               sum = sum - a(i,j)*b(j)
 13         continue
         endif
         b(i) = sum/a(i,i)
 14   continue

   END SUBROUTINE lubksb




