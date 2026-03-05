c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine dkron_matvec(trans,errscale,alpha,x,beta,y)
c     
c     Compute the matrix vector product 
c     
c     y <- alpha (F1 @ G1 + F2 @ G2) * x + beta y
c     
c     where '@' denotes a 'pseudo-Kronecker' product, and '*' is normal
c     matrix multiplication. Fi and Gi could describe the radial and
c     the lattitudinal behavior of the kernels for the helioseismic 2-D
c     rotational inversion problem respectively.

      implicit none
      include '2dsola.n.include.f'

      character*1 trans, errscale
      double precision alpha,beta
      double precision x(*),y(*)

      real time1,timearray(2)
      real dtime
      external dtime

      integer k,iblk,ylocal1,ylocal2,ylocal3
      integer inl_first,n_nl,ilm_first,n_lm,xlocal,k1,k2
      integer lwrk
      parameter(lwrk=2*maxtheta*maxmodes+
     c     1+maxnlm+maxpts)
      double precision wrk(lwrk), xtmp(maxnlm)
c     
c     External subroutines and functions called
c     
      logical lsame
      external lsame,dgemm,dscal

c*******************dkron_matvec begins here **********************

      if (lsame(trans,'t')) then         
c     
c*******************TRANSPOSED **********************
c     
c     y <- beta*y
c     
         if (beta.eq.0d0) then
            call pdzero(N_points,y,1)
         else
            call pdscal(N_points,beta,y,1)
         endif

         call pdcopy(M_kers,x,1,xtmp,1)
         if (lsame(errscale,'y')) then
            call covar_solve('t',xtmp)
         endif
C     
         if ((N_rad*(M_kers+M_lm*N_theta)).gt.
     c        (N_theta*(M_kers+M_nl*N_rad))) then
c     
c     Case: y <- VEC( (G * X) * F')
c     
            ylocal1 = 1
            ylocal2 = ylocal1+N_theta*M_nl
            xlocal  = ylocal2+N_theta*M_nl
            call pdzero(2*N_theta*M_nl,wrk,1)
c     
c     ylocal1 <- G1*X,  ylocal2 <- G2*X
c     
            do iblk=1,nlblocks
               do k=1,nblocks(iblk)
                  k1=iblock(k,iblk)
                  k2=iblock(k+1,iblk)-1
                  inl_first = inl(k1)
                  ilm_first = ilm(k1)
                  n_nl = inl(k2) - inl_first+1
                  n_lm = ilm(k2) - ilm_first+1
                  call  dgemm('n','n',N_theta,n_nl,n_lm,1D0,
     c                 g1(1,ilm_first), 
     c                 maxtheta,xtmp(k1),n_lm, 1D0,
     c                 wrk(ylocal1+(inl_first-1)*N_theta),N_theta)
                  call  dgemm('n','n',N_theta,n_nl,n_lm,1D0,
     c                 g2(1,ilm_first), 
     c                 maxtheta,xtmp(k1),n_lm, 1D0,
     c                 wrk(ylocal2+(inl_first-1)*N_theta),N_theta)
               enddo
            enddo
c     
c     y <- alpha*(ylocal1*transpose(F1) + ylocal2*transpose(F2)) + y
c     
            call dgemm('n','t',N_theta,N_rad,M_nl,alpha,
     c           wrk(ylocal1),N_theta,f1,maxrad,1D0,y,N_theta)
            call dgemm('n','t',N_theta,N_rad,M_nl,alpha,
     c           wrk(ylocal2),N_theta,f2,maxrad,1D0,y,N_theta)
         else
c     
c     Case: y <- VEC( G * (X * F'))
c     
            ylocal1 = 1
            ylocal2 = ylocal1+N_rad*M_lm
            xlocal = ylocal2+N_rad*M_lm
            call pdzero(2*N_rad*M_lm,wrk,1)
c     
c     ylocal1 <-( X*transpose(F1),   ylocal2 <- X*transpose(F2))
c     
            do iblk=1,nlblocks
               do k=1,nblocks(iblk)
                  k1=iblock(k,iblk)
                  k2=iblock(k+1,iblk)-1
                  inl_first = inl(k1)
                  ilm_first = ilm(k1)
                  n_nl = inl(k2) - inl_first+1
                  n_lm = ilm(k2) - ilm_first+1
                  call  dgemm('n','t',n_lm,N_rad,n_nl,1D0,xtmp(k1),
     c                 n_lm, f1(1,inl_first),maxrad,1D0,
     c                 wrk(ylocal1+ilm_first-1),M_lm)
                  call  dgemm('n','t',n_lm,N_rad,n_nl,1D0,xtmp(k1),
     c                 n_lm, f2(1,inl_first),maxrad,1D0,
     c                 wrk(ylocal2+ilm_first-1),M_lm)
               enddo
            enddo
c     
c     y  <- alpha * (G1*ylocal1 + G2*ylocal2) + y
c     
            call dgemm('n','n',N_theta,N_rad,M_lm,alpha,g1,maxtheta,
     c           wrk(ylocal1),M_lm, 1D0,y,N_theta)
            call dgemm('n','n',N_theta,N_rad,M_lm,alpha,g2,maxtheta,
     c           wrk(ylocal2),M_lm, 1D0,y,N_theta)
         endif
c     
      else
c     
c**************************NOT TRANSPOSED ***************************
c     
c     
         if  ((N_rad*(M_kers+M_lm*N_theta)).lt.
     c        (N_theta*(M_kers+M_nl*N_rad))) then
c     
c     Case: y <- VEC( (G' * X) * F )
c     
            ylocal1 = 1
            ylocal2 = ylocal1+N_rad*M_lm
            ylocal3 = ylocal2+N_rad*M_lm
c
c     ylocal1 <- alpha*transpose(G1)*X, ylocal2 <- alpha*transpose(G2)*X
c     
            call dgemm('t','n',M_lm,N_rad,N_theta,alpha,g1,maxtheta,x,
     c           N_theta, 0D0,wrk(ylocal1),M_lm)
            call dgemm('t','n',M_lm,N_rad,N_theta,alpha,g2,maxtheta,x,
     c           N_theta, 0D0,wrk(ylocal2),M_lm)
c     
c     y = ylocal1*F1 + ylocal2*F2
c     
            do iblk=1,nlblocks
               do k=1,nblocks(iblk)
                  k1=iblock(k,iblk)
                  k2=iblock(k+1,iblk)-1
                  inl_first = inl(k1)
                  ilm_first = ilm(k1)
                  n_nl = inl(k2) - inl_first+1
                  n_lm = ilm(k2) - ilm_first+1
                  call  dgemm('n','n', n_lm, n_nl, N_rad, 1D0,
     c                 wrk(ylocal1+ilm_first-1), M_lm, f1(1,inl_first),
     c                 maxrad, 0D0, xtmp(k1), n_lm)
                  call  dgemm('n','n', n_lm, n_nl, N_rad, 1D0,
     c                 wrk(ylocal2+ilm_first-1), M_lm, f2(1,inl_first),
     c                 maxrad, 1D0, xtmp(k1), n_lm)
               enddo
            enddo
         else
c     
c     Case: y <- VEC( G' * (X * F) )
c     
            ylocal1 = 1
            ylocal2 = ylocal1+N_theta*M_nl
            ylocal3 = ylocal2+N_theta*M_nl
c     
c     ylocal1 <- alpha * X * F1, ylocal2 <- alpha * X * F2
c     
            time1=dtime(timearray)
            call dgemm('n','n',N_theta,M_nl,N_rad,alpha,x,N_theta,f1,
     c           maxrad,0D0,wrk(ylocal1),N_theta)
            call dgemm('n','n',N_theta,M_nl,N_rad,alpha,x,N_theta,f2,
     c           maxrad,0D0,wrk(ylocal2),N_theta)
c     
c     y <- transpose(G1)*ylocal1 + transpose(G2)*ylocal2
c     
            do iblk=1,nlblocks
               do k=1,nblocks(iblk)
                  k1=iblock(k,iblk)
                  k2=iblock(k+1,iblk)-1
                  inl_first = inl(k1)
                  ilm_first = ilm(k1)
                  n_nl = inl(k2) - inl_first+1
                  n_lm = ilm(k2) - ilm_first+1
                  call  dgemm('t','n', n_lm, n_nl, N_theta, 1D0,
     c                 g1(1,ilm_first), maxtheta, wrk(ylocal1+
     c                 (inl_first-1)*N_theta), N_theta, 0D0,
     c                 xtmp(k1), n_lm)
                  
                  call  dgemm('t','n', n_lm, n_nl, N_theta, 1D0,
     c                 g2(1,ilm_first), maxtheta, wrk(ylocal2+
     c                 (inl_first-1)*N_theta), N_theta, 1D0, 
     c                 xtmp(k1), n_lm)
               enddo
            enddo
         endif
         if (lsame(errscale,'y')) then
            call covar_solve('n',xtmp)
         endif
         call pdaxpby(M_kers,1d0,xtmp,1,beta,y,1)
      endif
      return 
      end

