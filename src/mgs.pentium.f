
c
c****************************************************************************
C     This simple version of MGS is faster on Pentium machines with the GNU
c     G77 compiler.

      subroutine  mod_gram_schmidt(n,k,V,ldv,vnew)
      implicit none
      integer n,k,ldv
      double precision V(ldv,*),vnew(*)
      integer i,j
      double precision s

c     Check for quick return
      if ((k.le.0).or.(n.le.0)) return
      do i=1,k
         s = 0d0
         do j=1,n
            s = s + V(j,i)*vnew(j)
         enddo
         do j=1,n
            vnew(j) = vnew(j) - s*V(j,i)
         enddo
      enddo
      end
