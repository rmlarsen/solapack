
c
c****************************************************************************
c     Modified Gram-Schmidt orthogonalization using BLAS.
c
c     Orthogonalizes vnew against the k columns of V:
c       for i = 1, ..., k:
c         s = dot(V(:,i), vnew)
c         vnew = vnew - s * V(:,i)
c

      subroutine  mod_gram_schmidt(n,k,V,ldv,vnew)
      implicit none
      integer n,k,ldv
      double precision V(ldv,*),vnew(*)
      integer i
      double precision s,ddot
      external ddot,daxpy

c     Check for quick return
      if ((k.le.0).or.(n.le.0)) return
      do i=1,k
         s = ddot(n, V(1,i), 1, vnew, 1)
         call daxpy(n, -s, V(1,i), 1, vnew, 1)
      enddo
      end
