c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine set_derivker
c     
      implicit none
      include '2dsola.n.include.f'
      integer i,j,k1,k2
      double precision surfterm,cs,Omega_S(maxtheta),f1old,f2old,tmp
c
c     Set up integrated kernels using the trpezoidal rule
c
      do i=1,M_nl
         f1old = 2*f1(1,i)
         f2old = 2*f2(1,i)
         f1(1,i) = 0d0
         f2(1,i) = 0d0
         do j=2,N_rad-1
            tmp = f1(j-1,i)+0.5*(f1(j,i)+f1old)
            f1old = f1(j,i)
            f1(j,i) = tmp
            tmp = f2(j-1,i)+0.5*(f2(j,i)+f2old)
            f2old = f2(j,i)
            f2(j,i) = tmp
         enddo         
         f1(j,i) = f1(j-1,i)+0.5*f1old + f1(j,i)
         f2(j,i) = f2(j-1,i)+0.5*f2old + f2(j,i)
      enddo
c
c     Set up surface rotation using simple 3-term expansion.
c
      do j=1,N_theta
         cs = cos(rthetw(j,2))
         Omega_S(j) = 455.0 - 51.2*cs**2 - 84*cs**4
      enddo
c
c     Calculate surface term
c
      do i=1,M_kers
         surfterm = 0d0
         k1 = inl(i)
         k2 = ilm(i)
         do j=1,N_theta
            surfterm = surfterm + Omega_S(j)* 
     c           (f1(N_rad,k1)*g1(j,k2) + f2(N_rad,k1)*g2(j,k2))
         enddo
c         write(*,*) 'surface term, split = ',surfterm,modeset(i,7)
         modeset(i,7) = -modeset(i,7) + surfterm
      enddo      
      return
      end
c     
c***************************************************************
c     
      subroutine solve_constraint(alpha, akh1, sqrtw)
c     
      implicit none
      include '2dsola.n.include.f'
      double precision alpha, akh1(*), sqrtw(*)
c
c     Local variables
c
      double precision dtmp
c
c     External procedures
c      
      double precision pdnrm2
      external pdscal,pdcopy,dkron_matvec,pdnrm2,pdset
c     
c     Compute kernel integrals
c     *** STEP 2.2 ***
c     
      write(6,*) 'Computing kernel integrals'
      
      call pdset(N_points,1d0,work,1)
      call dkron_matvec('n','y',1D0,work,0D0,hhvec(3))

      hhvec(1)=0D0
      hhvec(2)=0D0
c     
c     Construct householder reflector
c
      write(6,*) 'Constructing householder reflector'
      
c      call gen_hh(M_kers, hhvec(3), hhvec(2), hhvec(1))

      alpha = pdnrm2(M_kers,hhvec(3),1)
      if (alpha.eq.0d0) stop 'Integrals of kernels are all zero!?'
      hhvec(2) = 1d0/(alpha*(alpha+dabs(hhvec(3))))
      hhvec(1) = -dsign(alpha,hhvec(3))
      hhvec(3) = hhvec(3) - hhvec(1)
      alpha = 1D0/hhvec(1)
      write (*,*) 'alpha = ',alpha
c
c     Compute alpha times the first column of the householder reflector 
c
c            alpha*h_1 = alpha * ( e_1 - tau * v(1) * v)
c
      dtmp = -hhvec(2)*hhvec(3)/hhvec(1)
      call pdaxpby(M_kers,dtmp,hhvec(3),1,0d0,work,1)
      work(1) = work(1)+alpha
c
c     Compute akh1 = alpha * W^(1/2) * K * h1
c
      call dkron_matvec('t','y',1D0,work,0D0,akh1)
      call pdaxty(N_points,1d0,invsqrtw,1,akh1,1)
      end
c     
c***************************************************************
c     
      subroutine set_u1(icase,u1,akh1,sqrtw)
      implicit none
      include '2dsola.n.include.f'
      
c     Set starting vector for Lanczos iteration to be the right-hand side 
c     corresponding to at target function at (\r_0,\theta_0) = (0.9,60 deg)
c     with widths (\Delta r, \Delta \theta) = (0.02, 0.1)

      integer icase
      double precision u1(*),akh1(*),sqrtw(*),rand
      integer i
      external rand

      call set_target(0.70D0,0.7D0,0.04D0,0.1D0,u1,icase,0)
      do  i=1,N_points
         u1(i) = sqrtw(i) * u1(i) - akh1(i)
      enddo
      end
c     
c***************************************************************
c     
      subroutine transform_rhs(splits,V,ldv) 
      implicit none
      include '2dsola.n.include.f'
      integer ldv
      double precision V(ldv,*)
      double precision splits(*)
c
c     Apply Householder transform to rhs (splittings) 
c     
      write (6,*) 'transform_rhs: Applying householder reflector to rhs'
      call dgemv('T',M_kers,1,1d0,splits,maxnlm,hhvec(3),1,0d0,work,1)
      call dger(M_kers,1,-hhvec(2),hhvec(3),1,work,1,splits,maxnlm)
c
c     Project splittings on orthogonal (Lanczos) basis:
c       splits(2:m) <- V^T * splits(2:m)      
c
      call dgemv('T',M_kers-1,N_iter,1d0,V,ldv,splits(2),1,0d0,work,1)
      call dcopy(M_kers-1, work, 1, splits(2),1)      
      end
c     
c***************************************************************
c     
      subroutine set_bidiag(U, ldu, bidiag, V,ldv) 
      implicit none
      include '2dsola.n.include.f'
      integer ldu, ldv
      double precision U(ldu,*), V(ldv,*)
      double precision bidiag(maxpts,*)
c
c     Compute Lanczos bidiagonalization
c
c         A * V_k = U_{k+1} * B_k  (k=N_iter)
c
c     by performing N_iter steps of the Golub-Kahan (Lanczos) 
c     bidiagonalization.
c
c     On entry the first  of U contains the starting 
c     vector for the Lanczos iteration.
c
      external daprod
      double precision dummy(1)
c     
      if (lwork .lt. (N_iter+M_kers+N_points)) then
         stop 'lwork too small in set_bidiag'
      endif
c     
c     Compute Lanczos bidiagonalization
c
      write (*,*) 'Calling dlanc_b'
      call dlanc_b(N_points,M_kers-1,0,N_iter,daprod,
     c     bidiag(1,1), bidiag(1,2),U,ldu,V,ldv,dummy,
     c     1, work)
c     
c     Compute QR-factorization B_k = Q_k R_k
c
      write (*,*) 'Calling dbdqr'
      call dbdqr('R','N',N_iter,N_points,bidiag(1,1),
     c     bidiag(1,2),U,ldu,work)
      return
      end
c     
c********************************************************************
c     
      subroutine set_transtarget(icase,sqrtw,target,ldt,akh1,U,ldu,
     c     dwork,ldwork)
c     
      implicit none
      include '2dsola.n.include.f'
      integer icase,ldu,ldt,ldwork
      double precision target(ldt,*),U(ldu,*)
      double precision sqrtw(*),akh1(*),dwork(ldwork)
      integer i,j
c     
c     set right-hand side, transpose(U) * (sqrt(W) * target - alpha*K*H_1)
c     
      if(icase.lt.100) then
c$doacross local(i,j),  shared(sqrtw,target,akh1,N_targets,N_points)
cc$PAR DOALL private(j), shared(target), readonly(sqrtw,akh1)
         do  j=1,N_targets
            do  i=1,N_points
               target(i,j) = sqrtw(i) * target(i,j) - akh1(i)
            enddo
         enddo          
      else
c$doacross  local(i,j),  shared(sqrtw,target,N_targets,N_points)
cc$PAR DOALL private(j), shared(target), readonly(sqrtw)
         do  j=1,N_targets
            do  i=1,N_points
               target(i,j) = target(i,j) * sqrtw(i)
            enddo
         enddo         
      end if
c     
c     Multiply target by transpose(U) 
c     
      write (6,*) 'set_transtarget: target  <- transpose(U) * target' 
      call dgemm_ovwr('T',N_iter,N_targets,N_points,1D0,U,ldu,0D0,
     c     target,ldt,dwork,ldwork)

      return
      end
c     
c**********************************************************************
c     
      subroutine solve_bidiag(target,ldt,bidiag,lambda,resnorm,solnorm,
     c     alpha,icase)
c     
      implicit none
      include '2dsola.n.include.f'
      integer icase,ldt
      double precision target(ldt,*),bidiag(maxpts,*)
      double precision lambda(*),resnorm(*),solnorm(*)
      double precision alpha
c     
c     Local Variables
c     
      character*1 uplo
      double precision ddot
      external dcopy, elden, dbdsolve,ddot,daxpy

      integer first,last,i,n,j
      integer numblocks,start(maxtargets), end(maxtargets)
      double precision thistrdoff
      double precision mywork(6*maxpts+maxtargets+1)
 
      if (N_targets.gt.0) then
c     
c     Solve for \xi_1, \xi_2,...,\xi_q using Eldens algorithm 
c     
         uplo = 'u'
         n = N_iter
         first = 1
         last = 1
         numblocks = 0
         thistrdoff = lambda(1)         
c     
c     WHILE (last <= N_targets)
 50      if (last.le.N_targets) then           
c     
c     For every group with the same value of lambda do...
c     
c     DO
 100        last = last + 1
            if (last.gt.N_targets) goto 150
            if (.not.(thistrdoff.eq.lambda(last))) goto 150
            goto 100
c     WHILE (last <= N_targets) AND (thistrdoff = trdoff(last)
 150        numblocks = numblocks+1 
            start(numblocks) = first
            end(numblocks) = last
            first = last
            if (first.le.N_targets) thistrdoff = lambda(first)
            goto 50
         endif

cc$doacross if (n*N_targets.gt.10000),
cc$& local(j,i,first,last,mywork),
cc$& shared(bidiag,target,resnorm,solnorm,lambda,alpha),
cc$& shared(n,uplo,numblocks,N_targets,start,end)
         do j=1,numblocks
            first = start(j)
            last = end(j)
            call dcopy(n,bidiag(1,1),1,mywork(1),1)
            call dcopy(n,bidiag(1,2),1,mywork(n+1),1)         
c     
c     Compute a QR-factorization of the augmented bidiagonal system
c     using Elden's algorithm.
c     
            write (*,*) 'Calling Elden'
            call elden(uplo,n,last-first,mywork(1),mywork(n+1), 
     c           sqrt(lambda(first)), target(1,first),ldt, 
     c           resnorm(first),mywork(2*n+1))
c     
c     Solve for \xi_first,...,\xi_last and overwrite
c     target(:,first:last-1) with solution
c     
            call dbdsolve(uplo,n,last-first,mywork(1),
     c           mywork(n+1),target(1,first),ldt,mywork(2*n+1))
c     
c     Compute solution norm ( = err. magnification) and 
c     residual norms (deviation from target).
c     
c$doacross local(i) share(target,first,last,n,solnorm)
cc$PAR DOALL shared(solnorm), readonly(n,target)
            do i=first,last-1
               solnorm(i) = ddot(n,target(1,i),1,target(1,i),1)
            enddo
            call daxpy(last-first,-lambda(first),solnorm(first),
     c           1, resnorm(first),1)
         enddo
         if (icase.lt.100) then
            do i=1,N_targets
               solnorm(i) = solnorm(i)+alpha*alpha
            enddo
         endif
         do i=1,N_targets
            resnorm(i) = sqrt(max(0d0,resnorm(i)))
            solnorm(i) = sqrt(max(0d0,solnorm(i)))
         enddo
      endif
      return
      end
c     
c*******************************************************************
c     
      subroutine set_rot(icase,rot,XI,ldxi,splits)
c     
c     set rotation rate, based on \Xi and transformed data beta
c     Note that inversion is assumed to be done for error-weighted case
c     
      implicit none
      include '2dsola.n.include.f'
      integer icase,ldxi
      double precision rot(*), XI(ldxi,*), splits(*)
      double precision alpha
      integer i,n
c     
c     Compute the rotation rate
c     
      write (6,*) 'set_rot: M_kers,N_points,N_targets,alpha = ',
     c     M_kers,N_points,N_targets,1D0/hhvec(1)

      write (6,*) 'set_rot: f_est = transpose(XI)*beta'
      n = N_iter
      call dgemv('T',n,N_targets,1D0,XI,ldxi,splits(2),1,
     c     0D0,rot,1)
      
      if (icase.lt.100) then
         alpha = splits(1)/hhvec(1)
         do i=1,N_targets
            rot(i) = rot(i)+alpha
         enddo
      endif
      return
      end
c     
c*****************************************************************
c     
      subroutine set_avker(icase, U, ldu, bidiag, akh1, XI, ldxi,
     c     sqrtw)
c     
c     Overwrite solution array XI with averaging kernels.
c     XI <- W^(1/2) * (alpha*K*h_1 + U*B*XI)
c     
      implicit none
      include '2dsola.n.include.f'
      integer icase,ldu,ldxi
      double precision U(ldu,*), bidiag(maxpts,*), akh1(*), sqrtw(*)
      double precision XI(ldxi,*)

      integer i,j
      external dbdmm_ovwr,daxpy

      write (6,*) 'set_avker:  XI <-- B * XI'     
      call dbdmm_ovwr('u',N_iter,N_targets,bidiag(1,1),bidiag(1,2),
     c     XI,ldxi)

      if (N_iter.lt.N_points) then
c$doacross local(i,j), share(XI,bidiag,N_iter,N_targets)
cc$PAR DOALL private(i), shared(XI)
         do j=1,N_targets
            do i=N_iter+1,N_points
               XI(i,j) = 0D0
            enddo
         enddo
      endif
c     
c     Multiply by target by U 
c     
      write (6,*) 'set_avker: XI <-- U * XI'
      call dgemm_ovwr('N',N_points,N_targets,N_iter,1D0,U,ldu,0D0,XI,
     c     ldxi,work,lwork)
c     
c    Add alpha*K*h_1 to each column of avker and divide out quadrature weights.
c     
      if (icase.lt.100) then
c$doacross local(i,j), shared(akh1,XI,sqrtw,N_targets,N_points)
cc$PAR DOALL private(i), shared(XI), readonly(invsqrtw,akh1)
         do j=1,N_targets
            do i=1,N_points
               XI(i,j) = 2*XI(i,j)*invsqrtw(i)+(akh1(i)-
     c              XI(i,j))*invsqrtw(i)
            enddo
         enddo
      else
c$doacross  local(i,j), shared(XI,sqrtw,N_targets,N_points)
cc$PAR DOALL private(i), shared(XI), readonly(sqrtw)
         do j=1,N_targets
            do i=1,N_points
               XI(i,j) = XI(i,j)*invsqrtw(i)
            enddo
         enddo         
      endif
      return
      end
c     
c*******************************************************************
c     
      subroutine set_correlation(fcorr, icase, trdoff, targetparms,rot,
     c     sigma_rot, depart, XI, ldxi)
c     
      implicit none
      include '2dsola.n.include.f'
      character*256 fcorr
      integer icase,ldxi, num
      double precision XI(ldxi,*),trdoff(*),targetparms(maxtargets,*),
     c     rot(*),sigma_rot(*),depart(*)
      double precision alpha
      integer first,i,j,blocksize,idx
      double precision ddot,invnorm(maxtargets), alpha2

      if (icase.lt.100) then
         alpha = 1D0/hhvec(1)
      else
         alpha = 0D0
      endif
      alpha2 = alpha*alpha      

c$doacross local(i)
cc$PAR DOALL shared(invnorm), readonly(alpha2,XI,N_iter)
      do i=1,N_targets
         invnorm(i) = 1d0/sqrt(alpha2+ddot(N_iter,XI(1,i),1,XI(1,i),1))
      enddo

      blocksize = lwork/N_targets
      first = 1
      do while (first.le.N_targets)  
         num = min(blocksize,N_targets-first+1)
         call dgemm('t','n',N_targets,num,N_iter,1d0,XI,ldxi,
     c        XI(1,first),ldxi,0d0,work,N_targets)
         do i=first,first+num-1
            do j=1,N_targets
               idx = j+(i-first)*N_targets
               work(idx) = (alpha2 + work(idx))*invnorm(i)*invnorm(j)
            enddo
         enddo
         if (first.eq.1) then
            call c_writeavker(ibyteswap,num,N_targets,icase,
     c           trdoff(first), targetparms(first,1), maxtargets, 
     c           rot(first), sigma_rot(first),depart(first), work, 
     c           N_targets, fcorr)
         else
            call c_appendavker(ibyteswap,num,N_targets,icase,
     c           trdoff(first), targetparms(first,1), maxtargets, 
     c           rot(first), sigma_rot(first),depart(first), work, 
     c           N_targets, fcorr)
         endif
         first = first+blocksize
      enddo
      return
      end
c     
c*******************************************************************
c     
      subroutine set_coeff_disk(fcoeff,icase, trdoff, targetparms,rot,
     c     sigma_rot, depart, V,ldv,XI,ldxi)
c     
c     Set inversion coefficients explicitly.
c     Overwrite V with coefficients: 
c     V <- H*[alpha*ones(1,ntargets); V*XI] = alpha*h1 + H2*V*XI
c
      implicit none
      include '2dsola.n.include.f'
      character*256 fcoeff
      integer icase,ldv,ldxi
      integer blksz, lwrk, num
      double precision V(ldv,*), XI(ldxi,*)
      double precision alpha, trdoff(*),targetparms(maxtargets,*),
     c     rot(*),sigma_rot(*),depart(*)
      integer i, first,j
      external dgemm, dgemv, dger
      
      blksz = lwork/M_kers-1
      lwrk = blksz*M_kers
      if (icase.lt.100) then
         alpha = 1D0/hhvec(1)
      else
         alpha = 0D0
      endif

      first = 1
      do while (first.le.N_targets)
         num = min(blksz,N_targets-first+1)
c
c     work <- [alpha*ones(1,ntargets); V*XI]
c
         do i=1,num
            work(1+(i-1)*M_kers) = alpha
         enddo
         call dgemm('N','N',M_kers-1,num,N_iter,1D0,V(2,1),
     c        ldv,XI(1,first),ldxi,0d0,work(2),M_kers)
c
c     Apply Householder transform
c     
c         write (6,*) 'set_coeff: q <- H * (alpha; q^hat)'
         call dgemv('T',M_kers,num,1D0,work,M_kers,hhvec(3),1,0D0,
     c        work(lwrk+1),1)
         call dger(M_kers,num,-hhvec(2),hhvec(3),1,work(lwrk+1),1,
     c        work,M_kers)   
c
c     Take out error weighting
c
cc$doacross  local(i,j)
cc$PAR DOALL shared(work), readonly(invsigma_split)
         do i=1,num
            do j=1,M_kers
               work(j+(i-1)*M_kers) =  work(j+(i-1)*M_kers) * 
     c              invsigma_split(j)
            enddo
         enddo
c
c     Write result to disk:
c         
         if (first .eq. 1) then
            call c_writeavker(ibyteswap,num,M_kers,icase,
     c           trdoff(first), targetparms(first,1), maxtargets, 
     c           rot(first), sigma_rot(first), depart, work, M_kers, 
     c           fcoeff)
         else
            call c_appendavker(ibyteswap,num,M_kers,icase,
     c           trdoff(first), targetparms(first,1), maxtargets, 
     c           rot(first), sigma_rot(first), depart, work, M_kers, 
     c           fcoeff)
         endif
         first = first+blksz
      enddo
      return
      end
c     
c************************************************************************
c     
      subroutine set_departure(avker, target, depart)
c     
c     evaluates departure integral for 2D SOLA inversion
c     
      implicit none
      include '2dsola.n.include.f'
      integer n
      double precision avker(*), target(*),depart
c     
      depart=0
      do n=1,N_points
         depart=depart+rthetw(n,3)*(avker(n)-target(n))**2
      enddo
c     
      return
      end
c
c*********************************************************************
c
      subroutine daprod(trans,m,n,alpha,x,beta,y)
c
c     Compute y <- alpha*op(A)*x + beta*y
c      
      implicit none
      include '2dsola.n.include.f'
      character*1 trans
      integer m,n
      double precision alpha,beta,a
      double precision x(*),y(*)
c
c     Local variables
c
      integer i,tmp,tmp1
c
c     External subroutines called
c      
      double precision pddot
      logical lsame
      external lsame,dkron_matvec,pddot
c
c     Check that dimensions are correct
c
      if ((.not.m.eq.N_points).or.(.not.n.eq.(M_kers-1))) then
         write (*,*) 'Error in daprod: Wrong dimensions'
         write (*,*) 'm,n = ',m,n
         stop 'wrong dimensions in daprod'
      endif
      
      if (lsame(trans,'n')) then
c
c     Compute y <- alpha *  W^(1/2) * Kt * H2 * x + beta y
c         
         tmp = 1
         tmp1 = tmp+M_kers
c
c     Apply Householder reflector to vector directly without calling app_hh.
c
         a = hhvec(2)*pddot(M_kers-1,x,1,hhvec(4),1)
         work(tmp)=-a*hhvec(3)
c$doacross share(work,x,a,hhvec)
cc$PAR DOALL shared(work), readonly(x,a,hhvec)
         do i=1,M_kers-1
            work(tmp+i) = x(i)-a*hhvec(3+i)
         enddo
         call dkron_matvec('t','y',alpha,work(tmp),0D0,work(tmp1))
         if (beta.eq.0D0) then
c$doacross share(y,work,invsqrtw)
cc$PAR DOALL shared(y), readonly(work,invsqrtw)
            do i=1,N_points
               y(i) = work(tmp1+i-1)*invsqrtw(i)
            enddo
         else
c$doacross share(y,work,invsqrtw,beta)
cc$PAR DOALL shared(y), readonly(work,invsqrtw,beta)
            do i=1,N_points
               y(i) = work(tmp1+i-1)*invsqrtw(i)+beta*y(i)
            enddo
         endif
      else
c
c     Compute y <- alpha * H2^T * Kt^T * W^(1/2) * x + beta y
c         
         tmp = 1
         tmp1 = tmp+N_points
c$doacross share(x,work,invsqrtw)
cc$PAR DOALL shared(work), readonly(x,invsqrtw)
         do i=1,N_points
            work(tmp+i-1) = x(i)*invsqrtw(i)
         enddo
c     tmp1 = alpha * Kt^T * W^(1/2) * x
         call dkron_matvec('n','y',alpha,work(tmp),0D0,work(tmp1))
c
c     Apply Householder reflector to vector directly without calling app_hh.
c
         a = hhvec(2)*pddot(M_kers,work(tmp1),1,hhvec(3),1)
         if (beta.eq.0D0) then
c$doacross share(y,work,a,hhvec)
cc$PAR DOALL shared(y), readonly(work,a,hhvec)
            do i=1,M_kers-1
               y(i) = work(tmp1+i)-a*hhvec(3+i)
            enddo
         else
c$doacross share(y,work,a,hhvec,beta)
cc$PAR DOALL shared(y), readonly(work,a,hhvec,beta)
            do i=1,M_kers-1
               y(i) = work(tmp1+i)-a*hhvec(3+i)+beta*y(i)
            enddo
         endif
      endif
      return 
      end
