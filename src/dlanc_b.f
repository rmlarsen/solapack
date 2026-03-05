c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine dlanc_b(m, n, k, ITER, APROD, alpha, beta, 
     c     U,  ldu, V, ldv, C, ldc, work)

      implicit none
      integer m, n, k, ITER, ldu, ldv, ldc
      double precision U(ldu,*), V(ldv,*), C(ldc,*), 
     c     alpha(*), beta(*), work(*)
      external APROD
c
c----------------------------------------------------------------------
c     DLANC_B: Lanczos bidiagonalization.
c----------------------------------------------------------------------
c   Purpose:
c       Computes ITER steps of the Lanczos bidiagonalization of a 
c       rectangular matrix A. If A is m x n DLANC_B computes matrices 
c       U (m x ITER+1),  V (n x ITER) and B (ITER+1 x ITER), such that
c
c                A ~= U * B * V^T 
c
c       Matrix B is lower bidiagonal with diagonal elements 
c       alpha(1),alpha(2),...alpha(ITER) and sub-diagonal elements 
c       beta(1),beta(2),...,beta(ITER) and the largest singular values
c       of B are good approximations of the singular values of A. Z
c       The columns of U and V are orthogonal and contain the left 
c       and right Lanczos vectors respectively.
c       Matrix A is only referenced via the user supplied subroutine 
c       APROD which computes the product of A or A^T times a vector. 
c       This is useful when A is either large, sparse or structured,
c       or if, for some other reason, the effect of A working on a 
c       vector can be computed efficiently without forming A 
c       explicitly.
c----------------------------------------------------------------------
c   Arguments:
c
c       REORTH   (input) INTEGER
c                Determines the type of reorthogonalization used for 
c                enforcing orthogonality among the two sets of Lanczos
c                vectors {u_i} and {v_i}, i=,1...,ITER, i.e. the 
c                columns of U and V.
c                Possible values:
c                  0 : No reorthogonalization.
c                  1 : Classical Gram-Schmidt reorthogonalization.
c                  2 : Modified Gram-Schmidt reorthogonalization.
c
c       M        (input) INTEGER
c                The number of rows of the matrix A.
c
c       N        (input) INTEGER
c                The number of columns of the matrix A.
c
c       ITER     (input) INTEGER    
c                Number of steps of the Lanczos bidiagonalization to compute.
c                
c       APROD    (input) EXTERNAL (subroutine)
c                User supplied subroutine that computes the product of a vector
c                by the matrix A or its transpose. APROD must be on the form
c
c                   SUBROUTINE APROD(TRANS, M, N, ALPHA, X, BETA, Y)
c                   CHARACTER*1 TRANS
c                   INTEGER M, N
c                   DOUBLE PRECISION A, B, X(*), Y(*)
c
c                If trans='N' APROD should compute
c                      y <- alpha * A * x + beta * y
c
c                If trans='T' APROD should compute
c                      y <- alpha * A^T * x + beta * y
c       
c       ALPHA    (output) DOUBLE PRECISION, dimension ITER.
c                The ITER diagonal elements of the B.
c
c       BETA     (output) DOUBLE PRECISION, dimension ITER.
c                The ITER sub-diagonal elements of the B.
c 
c       U        (input/output) DOUBLE PRECISION, dimension (LDU,ITER+1).
c                On input U(1:M,1) contains the starting vector for the 
c                Lanczos bidiagonalization. If the computed bidiagonalization 
c                is to be used in solving a least-squares problem, it is 
c                advisable to set the starting vector equal to the right-hand
c                side. If || U(1:M,1) ||_2 = 0 a random starting vector is
c                generated.
c                On exit the columns of U contain the left Lanczos vectors.
c
c       LDU      (input) INTEGER
c                Leading dimension of the array U.
c 
c       V        (output) DOUBLE PRECISION, dimension (LDV,ITER).
c                On exit the columns of V contain the right Lanczos vectors.
c
c       LDV      (input) INTEGER
c                Leading dimension of the array V.
c
c       C        (input/output) DOUBLE PRECISION array, dimension (LDC,K)
c                On output the array C is overwritten with V^T * C.
c
c       LDC      (input) INTEGER
c                Leading dimension of the array C.
c
c       WORK     (workspace) DOUBLE PRECISION array, dimension 
c                ITER*(N+M+K+2) + M + MAX(M,N) if REORTH = 'H'
c                ITER*K + MAX(M,N) otherwise.
c
c----------------------------------------------------------------------
c
c     (C) Rasmus Munk Larsen, JILA, 1997
c
      logical lsame
      double precision dlamch,rand,pdnrm2
      external pdnrm2,pdscal,rand
      external dlamch,lsame      
c     
      integer itmpC,iwork
      integer it,i,j,timenum1,timenum2
      
      real stoptimer, reorth_time, dtime, tarr(2),t,tsum,treo
      external stoptimer,starttimer,getnewtimer
      

      call startnewtimer(timenum2)
      call getnewtimer(timenum1)      
      t=dtime(tarr)
      if ((m.le.0).or.(n.le.0).or.(ITER.le.0)) return

      itmpC = 1 
      iwork = itmpC+k*ITER
c     
c     If u is zero then generate a random starting vector
      beta(1) = pdnrm2(m,u,1)
      if (beta(1).eq.0d0) then
         write (*,*) 'DLANC_B: WARNING! Generating random starting' //
     c        'vector.'
         do i=1,m
            u(i,1) = rand()
         enddo
         beta(1) = pdnrm2(m,u,1)         
      endif

c     Make sure starting vector has unit norm.
      if (beta(1).ne.1D0) then
         call pdscal(m,1D0 / beta(1), u, 1)
      endif
      
c     
c     Start Lanczos bidiagonalization iteration.
c     
      reorth_time = 0e0
      tsum = 0e0
      treo = 0e0

      do it=1,ITER
         if (mod(it,10).eq.0) write (*,*) 'IT = ',it
         if (it.eq.1) then
            call pdzero(n,v,1)
            call APROD('t',m,n,1D0,u,0D0,v)
         else
            call pdcopy(n,v(1,it-1),1,v(1,it),1)
            call APROD('t',m,n,1D0,u(1,it),-beta(it-1),v(1,it))
         endif
c     
c     *********** REORTHOGONALIZATION OF V_IT ***************
c     
         tsum=tsum+dtime(tarr)
         call starttimer(timenum1)

         call mod_gram_schmidt(n,it-1,v,ldv,v(1,it))
         alpha(it) = pdnrm2(n,v(1,it),1)

         treo = treo + dtime(tarr)
         reorth_time = reorth_time + stoptimer(timenum1)         
         call pdscal(n,1D0/alpha(it),v(1,it),1)

c     ******************************************************

c     Compute next row in V^T*C
         if (k.gt.0) then
            call dgemv('T',n,k,1D0,C,ldc,v(1,it),1,0D0,
     c           work(itmpC+it-1),ITER)
         endif

         call pdcopy(m,u(1,it),1,u(1,it+1),1)
         call APROD('n',m,n,1D0,v(1,it),-alpha(it),u(1,it+1))
            
            
c     *********** REORTHOGONALIZATION OF U_IT ***************
         call starttimer(timenum1)
         tsum=tsum+dtime(tarr)

         call mod_gram_schmidt(m,it,u,ldu,u(1,it+1))
         beta(it) = pdnrm2(m,u(1,it+1),1)

         treo=treo+dtime(tarr)
         reorth_time = reorth_time + stoptimer(timenum1)
         call pdscal(m,1D0/beta(it),u(1,it+1),1)
      enddo         

c     C <- tmpC

      do i=1,k
         do j=1,ITER
            C(j,i) = work(itmpC+ITER*(i-1)+j-1)
         enddo
      enddo
      write (*,*) 'Time spent reorthogonalizing = ',reorth_time,treo
      write (*,*) 'Total time in dlanc_b        = ',
     c     stoptimer(timenum2),tsum+treo
      return
      end

