c
c****************************************************************************
c
C This version of MGS is faster on MIPS and UltraSPARC.

      subroutine mod_gram_schmidt(n,k,V,ldv,vnew)
c     
c     Modified Gram-Schmidt orthogonalization:
c     Orthogalizes vnew against the k vectors in V by the
c     iterative process     
c     
c     FOR i= [s_1:e_1 s_2:e_2 ... s_l:e_l] DO          
c       vnew = vnew - DOT( V(:,i), vnew ) * V(:,i) 
c
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n,k,ldv
      double precision V(ldv,*),vnew(*)
c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero
      parameter(one = 1.0, zero = 0.0)
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j
      double precision vn0,newcoef,coef

c     Check for quick return
      if ((k.le.0).or.(n.le.0)) return
c     Select the next block of columns from V
      coef = zero
c$doacross local(j) shared(V,vnew,n) reduction(coef)
c$& mp_schedtype=simple
c!$OMP PARALLEL DO local(j) shared(V,vnew,n) reduction(+:coef)
cc$PAR DOALL readonly(vnew,V), reduction(coef)
      do j=1,n
         coef = coef + V(j,1)*vnew(j)
      enddo         
c   interleaved (software pipelined) loops improve performance
c   of inner loop on machines with fused multiply-add.
      do i=2,k
         newcoef = zero
c$doacross local(j,vn0) shared(vnew,V,i,n,coef) reduction(newcoef)
c$& mp_schedtype=simple
c!$OMP PARALLEL DO local(j,vn0) shared(vnew,V,i,n,coef)
c!$&  reduction(+:newcoef)
cc$PAR DOALL private(vn0), readonly(i,coef,V), shared(vnew), reduction(newcoef)
         do j=1,n
            vn0 = vnew(j) - coef*V(j,i-1)
            newcoef = newcoef + vn0*V(j,i) 
            vnew(j) = vn0
         enddo
         coef = newcoef
      enddo
c$doacross local(j) shared(vnew,V,k,coef,n)
c$& mp_schedtype=simple
c!$OMP PARALLEL DO local(j) shared(vnew,V,k,coef,n)
cc$PAR DOALL readonly(k,coef,V), shared(vnew)
      do j=1,n
         vnew(j) = vnew(j) - coef*V(j,k)
      enddo
      end
