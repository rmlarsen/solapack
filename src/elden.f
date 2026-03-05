c     Copyright Rasmus Munk Larsen, Stanford University, 2003

c
c****************************************************************************
c

      subroutine elden(uplo, n, nrhs, D, E, lambda, RHS, ldrhs, resnorm,
     *     work)

      implicit none
      integer n, nrhs, ldrhs
      character*1 uplo
      double precision D(*),E(*),RHS(ldrhs,*),resnorm(*),work(*)
      double precision lambda
      
c  
c ELDEN: Compute QR (QL) factorization of a regularized n x n upper (lower) 
c        bidiagonal matrix B with diagonal elements D(1),...,D(n), and
c        off-diagonal elements E(1),...,E(n-1).
c
c        This is done by applying a series of 2d-1 Givens rotations 
c        G_1,G_2,...,G_2n-1, such that
c
c        G_2n-1 * G_2n-2...G_1 * (    B     ) = transpose(Q) * (    B     )
c                                (lambda*I_n)                  (lambda*I_n)
c                                             = ( R )
c                                               ( 0 )
c
c        where R is an k x k upper (lower) bidiagonal matrix. 
c        This is illustrated below for upper bidiagonal B:
c
c                      ( D(1) E(1)             | RHS(1,1) ... RHS(1,nrhs) )
c                      (      D(2) E(2)        | RHS(2,1) ... RHS(2,nrhs) )
c                      (           ...         |  :               :       )
c                      (                E(n-1) |  :               :       )
c  G_1 * G_2..G_2n-1 * (                D(n)   | RHS(n,1) ... RHS(n,nrhs) )
c                      ( lambda                |  0       ...      0      )
c                      (      lamdba           |  0       ...      0      )
c                      (           ...         |  :                :      )
c                      (                lambda |  0       ...      0      )
c
c              ( D(1)' E(1)'               | RHS(1,1)' ... RHS(1,nrhs)' )
c              (       D(2)' E(2)'         | RHS(2,1)' ... RHS(2,nrhs)' )
c              (             ...           |  :               :         )
c              (                   E(n-1)' |  :               :         )
c    ---->     (                   D(n)'   | RHS(n,1)' ... RHS(n,nrhs)' )
c              (  0    ...   ...   0       |  r(1,1)   ...  r(1,nrhs)   )
c              (  :     0          :       |  r(2,1)   ...  r(2,nrhs)   )
c              (  :           .    :       |  :         :      :        )
c              (  0     ...        0       |  r(n,1)   ...  r(1,nrhs)   )
c
c        D is overwritten with R(1,1),R(2,2),...,R(n,n), and E is overwritten
c        with R(1,2),R(2,3),...,R(n-1,n).
c        If nrhs>0 the matrix RHS is overwritten by
c
c          RHS' <- Gn-1 * Gn-2 * ...  * G1 * RHS = transpose(Q) * RHS,
c
c        i.e. after the QR-factorization, the solutions to the regularized
c        linear least squares problems
c
c            min { || B X(:,i) - RHS(:,i) ||_2^2 + 
c                  lambda^2 || X(:,i) ||_2^2 } , i=1,...,nrhs
c
c        are simply found by solving the following bidiagonal system with 
c        multiple right-hand sides
c 
c                      R X = RHS'
c
c        On exit resnorm(i) = || r(:,i) ||_2^2, i.e. the square of the residual
c        norm corresponding to the ith right-hand size.
c
c        See also DBDQR, DBDSOLVE and DBDMM.
c
c     (C) Rasmus Munk Larsen, JILA, 1997
c

      integer i,j,cc1,ss1,cc2,ss2
      double precision a,b,s,c,t,t1,lam,gamma
      external drotg,drot,dscal,lsame
      logical lsame


      if (lambda.eq.0.0) return
c
c     Noting to do. System is already on bidiagonal form.
c

      lam = lambda
      do i=1,nrhs
         resnorm(i) = 0.0
      enddo
      do i=1,nrhs+4*(n+1)
        work(i) = 0.0
      enddo
c
c Set offsets in work array for 'arrays'  cc1(*),ss1(*),
c cc2(*) and ss2(*) which are used for collecting the Givens rotations.
c
      cc1 = nrhs
      ss1 = cc1 + n+1
      cc2 = ss1 + n+1
      ss2 = cc2 + n+1
c
c Compute QR factorization. Collect Givens rotations in cc1(*),ss1(*),
c cc2(*) and ss2(*) and apply them in parallel to RHS afterwards.
c
      if (lsame(uplo,'L')) then
c
c     Matrix is LOWER bidiagonal.
c
         do i=n,1,-1
c
C     inlined call givens(D(i),lam,c,s)
c
            if (lam.eq.0.0) then
               c = 1.0
               s = 0.0
            else
               if (abs(lam).gt.abs(D(i))) then
                  t = -D(i)/lam 
                  t1 = dsqrt(1.0+t**2)
                  D(i) = -lam*t1
                  s = 1.0/t1
                  c = s*t
               else
                  t = -lam/D(i) 
                  t1 = dsqrt(1.0+t**2)
                  D(i) = D(i)*t1
                  c = 1.0/t1
                  s = c*t
               endif
            endif            
c
c end inlined call givens(D(i),lam,c,s)
c
            work(cc1+i)=c
            work(ss1+i)=s
c            D(i) = c*D(i) - s*lam
            if (i.gt.1) then
               gamma = s*E(i-1)
               E(i-1) = c*E(i-1)
c
C     inlined call givens(lambda,gamma,c,s)
c
               if (gamma.eq.0.0) then
                  c = 1.0
                  s = 0.0
               else
                  if (abs(gamma).gt.abs(lambda)) then
                     t = -lambda/gamma 
                     t1 = dsqrt(1.0+t**2)
                     lam = -gamma*t1
                     s = 1.0/t1
                     c = s*t
                  else
                     t = -gamma/lambda
                     t1 = dsqrt(1.0+t**2)
                     lam = lambda*t1
                     c = 1.0/t1
                     s = c*t
                  endif
               endif
c
c     end inlined call givens(lambda,gamma,c,s)
c
               work(cc2+i) = c
               work(ss2+i) = s
c               lam = c*lambda - s*gamma            
            endif
         enddo
         
c
c     Apply the orthogonal transformations in parallel 
c     to the right-hand sides.
c
c$doacross  local(i,j,a,b,t1,t),
c$& share(RHS,work,resnorm,cc1,ss1,cc2,ss2,nrhs,n)
cc$PAR DOALL private(i,a,b,t1,t), readonly(cc1,cc2,ss1,ss2,n,work)
cc$PAR DOALL shared(RHS,resnorm)
         do j=1,nrhs
            t1 = 0d0
            t = 0d0
            do i=n,2,-1
               a = work(cc1+i)*RHS(i,j) - work(ss1+i)*t1
               b = work(ss1+i)*RHS(i,j) + work(cc1+i)*t1
               RHS(i,j) = a
               t1 = b
               t = t + (work(cc2+i)*t1)**2
               t1 = -work(ss2+i)*t1
            enddo
            a = work(cc1+1)*RHS(1,j) - work(ss1+1)*t1
            b = work(ss1+1)*RHS(1,j) + work(cc1+1)*t1
            RHS(1,j) = a
            t1 = b
            resnorm(j) = t + t1**2
         enddo 

      else
c
c     Matrix is UPPER bidiagonal.
c
         do i=1,n
c
C     inlined call givens(D(i),lam,c,s)
c
            if (lam.eq.0.0) then
               c = 1.0
               s = 0.0
            else
               if (abs(lam).gt.abs(D(i))) then
                  t = -D(i)/lam 
                  t1 = dsqrt(1.0+t**2)
                  D(i) = -lam*t1
                  s = 1.0/t1
                  c = s*t
               else
                  t = -lam/D(i) 
                  t1 = dsqrt(1.0+t**2)
                  D(i) = D(i)*t1
                  c = 1.0/t1
                  s = c*t
               endif
            endif            
c
c end inlined call givens(D(i),lam,c,s)
c
            work(cc1+i)=c
            work(ss1+i)=s
c            D(i) = c*D(i) - s*lam
            if (i.lt.n) then
               gamma = s*E(i)
               E(i) = c*E(i)
c
C     inlined call givens(lambda,gamma,c,s)
c
               if (gamma.eq.0.0) then
                  c = 1.0
                  s = 0.0
               else
                  if (abs(gamma).gt.abs(lambda)) then
                     t = -lambda/gamma 
                     t1 = dsqrt(1.0+t**2)
                     lam = -gamma*t1
                     s = 1.0/t1
                     c = s*t
                  else
                     t = -gamma/lambda 
                     t1 = dsqrt(1.0+t**2)
                     lam = lambda*t1
                     c = 1.0/t1
                     s = c*t
                  endif
               endif

c     end inlined call givens(lambda,gamma,c,s)
c
               work(cc2+i) = c
               work(ss2+i) = s
c               lam = c*lambda - s*gamma            
            endif
         enddo

c
c     Apply the orthogonal transformations in parallel 
c     to the right-hand sides.
c
c$doacross  local(i,j,a,b,t,t1), 
c$&     share(RHS,work,resnorm,cc1,ss1,cc2,ss2,nrhs,n)
cc$PAR DOALL private(i,a,b,t1,t), readonly(cc1,cc2,ss1,ss2,n,work)
cc$PAR DOALL shared(RHS,resnorm)
         do j=1,nrhs
            t1 = 0d0
            t = 0d0
            do i=1,n-1
               a = work(cc1+i)*RHS(i,j) - work(ss1+i)*t1
               b = work(ss1+i)*RHS(i,j) + work(cc1+i)*t1
               RHS(i,j) = a
               t1 = b
               t = t + (work(cc2+i)*t1)**2
               t1 = -work(ss2+i)*t1
            enddo
            a = work(cc1+n)*RHS(n,j) - work(ss1+n)*t1
            b = work(ss1+n)*RHS(n,j) + work(cc1+n)*t1
            RHS(n,j) = a
            t1 = b
            resnorm(j) = t + t1**2
         enddo                   
      endif    
      return
      end
c
c****************************************************************************
c
      subroutine dbdqr(side, trans, n, k, D, E, C , ldc, work)
      implicit none
      character*1 side,trans
      integer n,k,ldc
      double precision D(*),E(*),C(ldc,*),work(*)
c
c DBDQR: Compute the R factor in the QR-factorization of a (n+1) x n lower
c        bidiagonal matrix B with diagonal elements D(1),...,D(n), and
c        sub-diagonals E(1),...,E(n). This is done by applying a series
c        of Givens rotations G_1,G_2,...,G_n, such that
c
c        G_n * G_n-1..G_1 * B = transpose(Q) * B = R
c
c        where R is an k x k upper bidiagonal matrix.
c
c        D is overwritten with R(1,1),R(2,2),...,R(n,n), and E is overwritten
c        with R(1,2),R(2,3),...,R(n-1,n). E(n) is set to zero.
c
c        If k>0 the matrix C is overwritten with:
c
c                        side = 'L'       side = 'R'
c        trans = 'N':    Q * C              C * Q
c        trans = 'T':    Q^T * C            C * Q^T
c
c        If side='L', then C is assumed to be (n+1) x k.
c        If side='R', then C is assumed to be k x (n+1).
c
c        Overwriting C can be useful when solving the linear least-squares 
c        problem with k right-hand sides
c
c                min || B X - C ||_2 .
c
c         If C is overwritten with Q^T * C then the solution X is found
c         after the call to DBDQR by solving the n x n upper bidiagonal system 
c
c                      R X = C(1:n,1:k)
c 
c         C(n+1,i) contains the residual corresponding to the i'th 
c         right-hand side.
c              
c
c     (C) Rasmus Munk Larsen, JILA, 1997
c
      integer i,j
      integer ic,is
      double precision ss,cc,t,t1
      double precision a,b
      logical lsame
      external lsame

      ic = 0
      is = n

      do i=1,n-1
c
c     inlined call givens(D(i),E(i),c,s)
c
         if (E(i).eq.0.0) then
            cc = 1.0
            ss = 0.0
         else
            if (abs(E(i)).gt.abs(D(i))) then
               t = -D(i)/E(i) 
               t1 = dsqrt(1.0+t**2)
               D(i) = -E(i)*t1
               ss = 1.0/t1
               cc = ss*t
            else
               t = -E(i)/D(i) 
               t1 = dsqrt(1.0+t**2)
               D(i) = D(i)*t1
               cc = 1.0/t1
               ss = cc*t
            endif
         endif
c     
c     end inlined call givens(D(i),E(i),c,s)
c
         work(ic+i) = cc
         work(is+i) = ss
c         D(i) = -ss*E(i)+cc*D(i)
         E(i) = -ss*D(i+1)
         D(i+1) = cc*D(i+1)
      enddo 
c      call givens(D(n),E(n),cc,ss)
c
c     inlined call givens(D(n),E(n),cc,ss)
c
         if (E(n).eq.0.0) then
            cc = 1.0
            ss = 0.0
         else
            if (abs(E(n)).gt.abs(D(n))) then
               t = -D(n)/E(n) 
               t1 = dsqrt(1.0+t**2)
               D(n) = -E(n)*t1
               ss = 1.0/t1
               cc = ss*t
            else
               t = -E(n)/D(n) 
               t1 = dsqrt(1.0+t**2)
               D(n) = D(n)*t1
               cc = 1.0/t1
               ss = cc*t
            endif
         endif
c     
c     end inlined call givens(D(n),E(n),cc,ss)
c
      work(ic+n) = cc
      work(is+n) = ss
c      D(n) = -ss*E(n)+cc*D(n)
      E(n) = 0.0

      if (k.gt.0) then
c
c  Apply Givens rotations to C: 
c
         if (lsame(side,'L') .and. lsame(trans,'T')) then
c         C <- transpose(Q) * C
c$doacross share(C,ic,is,k,n), local(i,j,a,b)
cc$PAR DOALL private(i,a,b), readonly(ic,is,n,work)
cc$PAR DOALL shared(C)
            do j=1,k
               do i=1,n
                  a = work(ic+i)*C(i,j) - work(is+i)*C(i+1,j)
                  b = work(is+i)*C(i,j) + work(ic+i)*C(i+1,j)
                  C(i,j) = a
                  C(i+1,j) = b
               enddo
            enddo 
         else if (lsame(side,'L') .and. lsame(trans,'N')) then
c        C <- Q * C
            write (*,*) 'side=L,  trans=N'
c$doacross share(C,ic,is,k,n), local(i,j,a,b)
cc$PAR DOALL private(i,a,b), readonly(ic,is,n,work)
cc$PAR DOALL shared(C)
            do j=1,k
               do i=n,1,-1
                  a = work(ic+i)*C(i,j) + work(is+i)*C(i+1,j)
                  b = -work(is+i)*C(i,j) + work(ic+i)*C(i+1,j)
                  C(i,j) = a
                  C(i+1,j) = b
               enddo
            enddo             
         else if (lsame(side,'R') .and. lsame(trans,'T')) then
c        C <- C * transpose(Q)
            do i=n,1,-1
               ss = work(is+i)
               cc = work(ic+i)
c$doacross share(C,ic,is,k,n,ss,cc), local(j,a,b)
cc$PAR DOALL private(j,a,b), readonly(ss,cc,i,ic,is,n,work)
cc$PAR DOALL shared(C)
               do j=1,k
                  a = cc*C(j,i) + ss*C(j,i+1)
                  b = -ss*C(j,i) + cc*C(j,i+1)
                  C(j,i) = a
                  C(j,i+1) = b
               enddo
            enddo             
         else if (lsame(side,'R') .and. lsame(trans,'N')) then
c        C <- C * Q
            do i=1,n
               ss = work(is+i)
               cc = work(ic+i)
c$doacross share(C,ic,is,k,n,ss,cc), local(j,a,b)
cc$PAR DOALL private(j,a,b), readonly(ss,cc,i,ic,is,n,work)
cc$PAR DOALL shared(C)
               do j=1,k
                  a = cc*C(j,i) - ss*C(j,i+1)
                  b = ss*C(j,i) + cc*C(j,i+1)
                  C(j,i) = a
                  C(j,i+1) = b
               enddo
            enddo             
         endif
      endif
      return
      end
c
c****************************************************************************
c
      subroutine dbdsolve(uplo, n, nrhs, D, E, RHS, ldrhs, work)
      implicit none
      character*1 uplo
      integer n,nrhs,ldrhs
      double precision D(*),E(*),RHS(ldrhs,*),work(*)
c
c     DBDSOLVE: solve system with bidiagonal matrix B with diagonal elements 
c               D(1),...,D(n), and off-diagonal elements E(1),...,E(n-1)
c               and multiple righthand sides.
c     
c
c     (C) Rasmus Munk Larsen, JILA, 1997
c
      integer i,j
      external lsame
      logical lsame

C compute reciprocals of diagonal elements to avoid division in inner loop.
      do i=1,n
         work(i) = 1D0/D(i)
      enddo
      if (lsame(uplo,'l')) then
c
c     Matrix is LOWER bidiagonal.
c
c$doacross local(i,j), share(D,E,RHS,n,nrhs)
cc$PAR DOALL private(i), readonly(E,work,n)
cc$PAR DOALL shared(RHS)
         do j=1,nrhs
            RHS(1,j) = RHS(1,j)*work(1)
            do i=2,n
               RHS(i,j) = (RHS(i,j) - E(i-1) * RHS(i-1,j)) * work(i)
            enddo
         enddo
      else
c
c     Matrix is UPPER bidiagonal.
c
c$doacross local(i,j), share(D,E,RHS,n,nrhs)
cc$PAR DOALL private(i), readonly(E,work,n)
cc$PAR DOALL shared(RHS)
         do j=1,nrhs
            RHS(n,j) = RHS(n,j)*work(n)
            do i=n-1,1,-1
               RHS(i,j) = (RHS(i,j) - E(i) * RHS(i+1,j)) * work(i)
            enddo
         enddo
      endif
      return
      end
c
c****************************************************************************
c
      subroutine dbdmm(uplo, m, n, D, E, A, lda, B,ldb)
      implicit none
      character*1 uplo
      integer m,n,lda,ldb
      double precision D(*),E(*),A(lda,*),B(ldb,*)
c
c     DBDMM: Multiply an m x m bidiagonal with an m x n general
c            matrix A store the result in th m x n general matrix B.
c
c     (C) Rasmus Munk Larsen, JILA, 1997
c
      integer i,j
      external lsame
      logical lsame

      if (lsame(uplo,'l')) then
c$doacross share(A,B,D,E,m,n),local(i,j)
cc$PAR DOALL private(i), readonly(A,E,D)
cc$PAR DOALL shared(B)
         do j=1,n
            B(1,j) = A(1,j)*D(1)
            do i=2,m
               B(i,j) = A(i-1,j)*E(i-1) + A(i,j)*D(i)
            enddo
         enddo
      else
c$doacross local(i,j), share(A,B,D,E,m,n)
cc$PAR DOALL private(i), readonly(A,E,D)
cc$PAR DOALL shared(B)
         do j=1,n
            do i=1,m-1
               B(i,j) = D(i)*A(i,j) + E(i)*A(i+1,j)
            enddo
            B(m,j) = D(m)*A(m,j) 
         enddo
      endif
      return
      end

c
c****************************************************************************
c
      subroutine dbdmm_ovwr(uplo, m, n, D, E, A, lda)
      implicit none
      character*1 uplo
      integer m,n,lda
      double precision D(*),E(*),A(lda,*)
c
c     DBDMM: Multiply an m x m bidiagonal B with an m x n general
c            matrix A and overwrite A with the result B*A.
c
c     (C) Rasmus Munk Larsen, JILA, 1997
c
      integer i,j
      external lsame
      logical lsame
      double precision Ai,oldAi

      if (lsame(uplo,'l')) then
c$doacross share(A,D,E,m,n),local(i,j,Ai,oldAi)
cc$PAR DOALL private(i,Ai,oldAi), readonly(E,D,m)
cc$PAR DOALL shared(A)
         do j=1,n
            oldAi = A(1,j)
            A(1,j) = A(1,j)*D(1)
            do i=2,m
               Ai = A(i,j)
               A(i,j) = oldAi*E(i-1) + A(i,j)*D(i)
               oldAi = Ai
            enddo
         enddo
      else
c$doacross local(i,j), share(A,D,E,m,n)
cc$PAR DOALL private(i), readonly(E,D,m)
cc$PAR DOALL shared(A)
         do j=1,n
            do i=1,m-1
               A(i,j) = D(i)*A(i,j) + E(i)*A(i+1,j)
            enddo
            A(m,j) = D(m)*A(m,j) 
         enddo
      endif
      return
      end

