c     Copyright Rasmus Munk Larsen, Stanford University, 2003
c     
c************************************************************************
c     
      subroutine set_chol_sigma(bandwidth)
      implicit none
      include '2dsola.n.include.f'
      integer i,info,idx,nexp,j,k,bandwidth
      real stoptimer
c     
c     Compute cholesky factorization of covariance matrices M_i for each 
c     multiplet i=1,...,M_nl such that M_i = C_i * C_i^T.
c     Set off-diagonal elements outside 'bandwidth' to zero, e.i.
c     m_{ij} where |i-j| > bandwidth is set to zero.
c     
      call starttimer(6)
      write (*,*) 'Factorizing covariance matrix. Bandwidth =',bandwidth
      write (*,*) 'M_nl, nnblocks =',M_nl,nnblocks
      write (*,*) 'inblock(1),inblock(nnblocks) =',inblock(1),
     c     inblock(nnblocks)

c$doacross  local(i,j,k,idx,nexp,info)
c$& mp_schedtype=simple
cc$PAR DOALL private(j,k,idx,nexp,info), shared(chol_sigma)
cc$PAR DOALL readonly(inblock,modeset)
      do i=1,nnblocks
         idx = inblock(i)
         nexp = modeset(idx,6)
         
c     write (*,*) 'l,n = ',modeset(idx,1),modeset(idx,2)
         do j=1,nexp
            do k=1,nexp
               if (abs(k-j).ge.bandwidth) chol_sigma(j,k,i) = 0.0
            enddo
         enddo
c     
c     Compute Cholesky factorization.
c     
         call dpotrf('l',nexp,chol_sigma(1,1,i),maxexp,info)

         if (info.ne.0) then
            write (*,*) 'Error in cholesky factorization, i,(l,n) = ',
     c           modeset(idx,1), modeset(idx,2)
            write(*,*) 'i,Info = ',i,info
            stop 'Covariance matrix not positive definite.'
         endif
      enddo
      write (*,*) 'Time spent in set_chol_sigma = ',stoptimer(6)
      return
      end
      

      subroutine covar_solve(trans,x)
      implicit none
      character*1 trans
      include '2dsola.n.include.f'
      double precision x(*)
      integer i,nexp,idx
c     
c     Multiply vector of splittings by the inverse Cholesky factors 
c     of covariance matrices,
c     i.e. solve 
c     
c     op(C_i) y(nblock(i):(nblock(i+1)-1)) = x(nblock(i):(nblock(i+1)-1)) 
c     
c     and place the result in x. 
c     If trans='n' the op(C_i) = c_i,
c     If trans = 't' then op(C_i) = transpose(C_i)
c     
      if (icovar.eq.1) then

c$doacross  local(i,idx,nexp)
c$& mp_schedtype=simple    
cc$PAR DOALL private(idx,nexp), shared(chol_sigma,x)
cc$PAR DOALL readonly(inblock,modeset)
         do i=1,nnblocks
            idx = inblock(i)
            nexp = modeset(idx,6)
            call dtrsv('l',trans,'n',nexp,chol_sigma(1,1,i),maxexp,
     c           x(idx),1)
         enddo
      else
c$doacross  local(i)
c$& mp_schedtype=simple
cc$PAR DOALL shared(x), readonly(invsigma_split)
         do i=1,M_kers
            x(i) = x(i) * invsigma_split(i)
         enddo
      endif
      return
      end
c     
c     
      subroutine read_covariance(ilun)
      implicit none
      include '2dsola.n.include.f'
      integer ilun
      
      integer i,j,k,l,n,nexp,idx
      
      do i=1,M_nl
         read(ilun,*) l,n
         idx = inblock(i)
         nexp = modeset(idx,6)
c     write(*,*) 'l,n,i,inblock(i),nexp = ',l,n,inblock(i),nexp
         if (l.ne.modeset(idx,1) .or. n.ne.modeset(idx,2)) then
            write (*,*) 'Covariance file at i,(l,n) = ',i,l,n
            write (*,*) 'is incompatible with '
            write (*,*) 'modeset. Expected (l,n) = ',
     c           modeset(idx,1),modeset(idx,2) 
            stop 'Incompatible Modeset and Covariance file'
         endif
         do j=1,nexp
            read(ilun,*) (chol_sigma(k,j,i),k=1,4)
            read(ilun,*) (chol_sigma(k,j,i),k=5,8)
            read(ilun,*) (chol_sigma(k,j,i),k=9,12)
            read(ilun,*) (chol_sigma(k,j,i),k=13,16)
            read(ilun,*) (chol_sigma(k,j,i),k=17,18)
c     write (*,*) (chol_sigma(k,j,i),k=1,18)
         enddo
      enddo
      return
      end

