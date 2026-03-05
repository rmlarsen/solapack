c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine dgemm_ovwr(transa,m,n,k,alpha,A,lda,beta,B,ldb,
     c     dwork,ldwork)
c
c     compute B <- alpha*op(A)*B + beta*B
c
c     Uses dwork as temporary storage, processing columns of B in
c     blocks that fit in the workspace.
c
      implicit none
      character*1 transa
      integer m,n,k,lda,ldb,ldwork
      double precision alpha,beta,A(lda,*),B(ldb,*),dwork(ldwork)
      integer i,j,blocksize,ncols
      external dgemm,dcopy,dscal,daxpy

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.m) stop 'Too little workspace in dgemm_ovwr'
      if (m.gt.ldb) stop 'm>ldb in DGEMM_OVWR'
      blocksize = ldwork/m

      do i=1,n,blocksize
         ncols = min(blocksize, n-i+1)
         call dgemm(transa,'N',m,ncols,k,alpha,A,lda,
     c              B(1,i),ldb,0D0,dwork,m)
         if (beta.eq.0D0) then
            do j=0,ncols-1
               call dcopy(m, dwork(j*m+1), 1, B(1,i+j), 1)
            enddo
         else
            call dscal(m*ncols, beta, B(1,i), 1)
            do j=0,ncols-1
               call daxpy(m, 1D0, dwork(j*m+1), 1, B(1,i+j), 1)
            enddo
         endif
      enddo
      return
      end


      subroutine dgemm_ovwr_left(transb,m,n,k,alpha,A,lda,beta,B,ldb,
     c     dwork,ldwork)
c
c     compute  A <- alpha*A*op(B)
c
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldwork
      double precision alpha,beta,A(lda,*),B(ldb,*),dwork(ldwork)
      integer i,j,blocksize,nrows
      external dgemm,dcopy

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.n) stop 'Too little workspace in dgemm_ovwr_left'
      blocksize = ldwork/n

      do i=1,m,blocksize
         nrows = min(blocksize, m-i+1)
         call dgemm('n',transb,nrows,n,k,alpha,A(i,1),lda,
     c              B,ldb,0d0,dwork,nrows)
         do j=0,n-1
            call dcopy(nrows, dwork(j*nrows+1), 1, A(i,j+1), 1)
         enddo
      enddo
      return
      end

