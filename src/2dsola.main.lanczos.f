c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine sola2d_lanczos(U,V,target)
c     
c     2D SOLA inversion programme, for rotational inversion
c     Takes kernels in the form of A matrix as produced by RLS
c     code, together with mesh and integration weight, and
c     carries out SOLA inversion.
c     
c     Original version 2/4/96
c     
c     Modified 14/8/96, to allow setting r_0 dependent target width
c     (using scaling with sound speed, obtained from reading amdl file).
c     
c     Modified 18/8/96, including option of inverting for radial gradient
c     Modified 18/8/96, introducing target with zero second moment
c     Modified 20/8/96, implementing new algorithm based on bidiagonalization
c     Modified 12/2/97, including L-curve analysis.
c     Modified june 1997, replacing central algorithms with Lanczos bidiagonalization
c     Modified  9/4/98, including covariance matrices for a-coef's in inversion.
c     Modified August 2003, major cleanup to remove unnecessary features and
c                           enable byteswapping I/O.
c     Rasmus Munk Larsen, Computer Sci. Dep., University of Aarhus.

c     PROGRAM BEGINS HERE: 

c     Try avoid shooting ourselves in the foot...
      implicit none
c     
c     U, V and target arrays are allocated and parsed to sola2d_lanczos by
c     the C wrapper.
c     NOTICE: in set_avker, averaging kernels overwrite transformed 
c     coefficients stored in array 'target'. 
c
      double precision U(*),V(*),target(*)
c     
c     Include global constants and common blocks.
c     
      include '2dsola.n.include.f'

c     Names of data files.
      character*256 fker, fmesh, fdata, fcovar, fsol, favker, fcoeff,
     *     fbidiag, ftranstarget, famdl, fcorr, ftarget

c     Arrays for storing bidiagonal and transformed kernel data.
      double precision bidiag(maxpts,2),akh1(maxnlm)
c     Square root of quadrature weights.
      double precision sqrtw(maxpts)
c     Frequency splittings data.
      double precision splits(maxnlm)
c     Target parameters
      double precision  targetparms(maxtargets,4),
     *     trdoff(maxtargets), depart(maxtargets),
     *     fitparms(maxtargets,4),cogparms(maxtargets,6)
c     Solution (rotation rate) and error estimates
      double precision rot(maxtargets),sigma_rot(maxtargets)
      common/maincom/bidiag,akh1,sqrtw,splits,targetparms,
     c     trdoff,depart,fitparms,cogparms

c     Control parameters (what to do, read-write control etc.)
      integer iprecision,irscal,irdbidiag, iwrbidiag, irdtranstarget, 
     *     icase, iwrtranstarget, itargettype,
     *     icalcsol, icalcavker, icalccoeff, icalccorr, itargt, 
     *     idocalc, bandwidth,iwrtarget

c     Misc. local variables (counters and the like).
      double precision alpha,xrscal
      integer i,j,idummy,idummy1,tnum,tnum1,ndata

c     External and intrinsic functions:
      real stoptimer, dtime 
      external dtime,stoptimer,starttimer
c
c**********************Program begin here ************************************
c
      call startnewtimer(tnum)
c     
c     If ibyteswap=0 then unformatted data files are read in as they are.
c     Otherwise they are byte-swapped to convert from little to big
c     endian format (or vice versa).
c     
      ibyteswap=0
      write(6,*) 'Enter ibyteswap'
      read(5,*) ibyteswap
      write(6,*) ibyteswap
c     
c     Set icase for SOLA. 
c     icase = icase0 + 100*icase1. 
c     Here icase0 determines target function. Possible values:
c     icase = 0: simple Gaussian target
c     icase = 1: flags for setting kernel for radial derivative.
c     itargettype=0: simple Gaussian target
c     itargettype=1:  Gaussian*sin(theta)**2
c     itargettype=2:  sinc function.
c     
      icase = 0
      itargettype = 0
      write(6,*) 'Enter icase, itargettype'
      read(5,*) icase,itargettype
      write(6,*) icase,itargettype
c     
c     Set number of iterations and reorthogonalization flag for
c     Lanczos bidiagonalization.
c     
      N_iter = 100
      write(6,*) 'Enter N_iter'
      read(5,*) N_iter
      write(6,*) N_iter
c     
c     Set flags for input or output of bidiagonal factorization
c     of kernel matrix and transformed target functions.
c     
      irdbidiag=0
      iwrbidiag=1
      irdtranstarget=0
      iwrtranstarget=1
c     
      write(6,*) 'Enter irdbidiag, iwrbidiag, irdtranstarget,'  //
     c     ' iwrtranstarget'
      read(5,*) irdbidiag, iwrbidiag, irdtranstarget, 
     *     iwrtranstarget
      write(6,*) irdbidiag, iwrbidiag, irdtranstarget, 
     *     iwrtranstarget
c     
c     Set iprecision to 0/1 for reading mesh and kernels in 
c     single/double precision.
c     
      iprecision = 0
      write (6,*) 'Enter iprecision'
      read(5,*) iprecision
      write (6,*) iprecision
c     
c     Set icovar=1 to read covarience matrix for splittings from 
c     file fcovar.
c     Set icovar=0 to assume uncorreated noise and use variance estimates given
c     in modeset file.
      icovar = 0
      fcovar = 'XXX'
      bandwidth = 20
      write (6,*) 'Enter icovar, bandwidth'
      read(5,*) icovar,bandwidth
      write (6,*) icovar,bandwidth
      write (6,*) 'Enter fcovar'
      read(5,'(a)') fcovar
      write (6,*) fcovar
c     
c     Set flags determining what to compute. 
c     
      icalcsol=1
      icalcavker=1
      icalccoeff=0
      icalccorr=1
c     
      write(6,*) 'Enter icalcsol, icalcavker, icalccorr, ' //
     c     'icalccoeff, iwrtarget'
      read(5,*)  icalcsol, icalcavker, icalccorr, icalccoeff,
     c     iwrtarget
      write(6,*) icalcsol, icalcavker, icalccorr, icalccoeff,
     c     iwrtarget      
c     
c     Set possible parameters for scaling of radial target width
c     If irscal = 1, scale with sound speed, relative to value at
c     xrscal
c     
      irscal=0
      xrscal=0.7
c     
      write(6,*) 'Enter irscal, xrscal'
      read(5,*) irscal, xrscal
      write(6,*) irscal, xrscal

      write(6,*) 'Enter amdl file to set sound speed for scaling'
      read(5,'(a)') famdl
      write(6,'(a)') famdl
c     
c     Check if there is anything to do...
c     
      if ((icalcsol.eq.0).and.(icalcavker.eq.0).and.
     c     (icalccoeff.eq.0).and.(icalccorr.eq.0)) then
         idocalc = 0
      else
         idocalc = 1
      endif
      if (idocalc.eq.0) then
         stop 'Nothing to do!'
      endif
c     
c     open files for input and output
c     
      write(6,*) 'Enter file name for kernel input'
      read(5,'(a)') fker
      write(6,'(a)') fker
c     
      write(6,*) 'Enter file name for mesh input'
      read(5,'(a)') fmesh
      write(6,'(a)') fmesh
c     
      write(6,*) 'Enter file name for data input'
      read(5,'(a)') fdata
      write(6,'(a)') fdata
c     
      write(6,*) 'Enter file name for solution output'
      read(5,'(a)') fsol
      write(6,'(a)') fsol
      if (icalcsol.eq.1) then
         open(ifsol,file=fsol,status='unknown')
      endif
c     
      write(6,*) 'Enter file name for averaging kernel output'
      read(5,'(a)') favker
      write(6,'(a)') favker
c     
      write(6,*) 'Enter file name for correlation function output'
      read(5,'(a)') fcorr
      write(6,'(a)') fcorr
c     
      write(6,*) 'Enter file name for inversion coefficients output'
      read(5,'(a)') fcoeff
      write(6,'(a)') fcoeff
c     
      write(6,*) 'Enter file name for factorization input or output'
      read(5,'(a)') fbidiag
      write(6,'(a)') fbidiag
c
      write(6,*) 
     *     'Enter file name for transformed target input or output'
      read(5,'(a)') ftranstarget
      write(6,'(a)') ftranstarget
c     
      write(6,*) 'Enter file name for target output'
      read(5,'(a)') ftarget
      write(6,'(a)') ftarget
c
c**************************Start reading data *************************
c     
      call initkronecker(iprecision,fdata,fker,fmesh,fcovar,bandwidth)
      write(*,*) 'M_kers, N_targets, N_iter = ',M_kers, N_targets,
     c     N_iter 
      ndata = M_kers

      if (icase.lt.0) then
         call set_derivker
         icase = 0
      endif
c     
c     Compute square root of quadrature weights - a quantity used in 
c     many places.
c     
      do i=1,N_points
         sqrtw(i) = sqrt(rthetw(i,3))
      enddo
c     
c     Scale the splittings by the errors.  
c
      if (modeset(1,8).eq.0.0) then
c     errorfree data?
         write(*,*) 'ASSUMING error free data!'
         do i=1,M_kers
            splits(i) = modeset(i,7)
            invsigma_split(i) = 1d0
         enddo
      elseif (icovar.eq.0) then
         do i=1,M_kers
            splits(i) = modeset(i,7)/modeset(i,8)
            invsigma_split(i) = 1d0/modeset(i,8)
         enddo
      else
         call dcopy(M_kers,modeset(1,7),1,splits,1)
         do i=1,M_kers
            invsigma_split(i) = 1d0/modeset(i,8)
         enddo             
         call covar_solve('n',splits)
      endif
c     
c     Set up targets and tradeoff parameters
c
      call def_targets(5,itargt,targetparms,trdoff)
      write (6,*) N_targets, ' sets of target parameters have been' //
     c     ' read'

********************Start computational part *************************
c     
c     Compute kernel integrals and set up Householder reflector to
c     remove linear constraint
c     
      call solve_constraint(alpha, akh1,sqrtw)
c     
c     test for reading or setting  bidiagonal factorization
c     
      if(irdbidiag.eq.0) then
         call set_u1(icase,U,akh1,sqrtw)

         call startnewtimer(tnum1)
         call set_bidiag(U,N_points,bidiag, V(2), M_kers)
         write(*,*) 'Time spent in bidiag = ',stoptimer(tnum1)

         if(iwrbidiag.eq.1) then
            write(6,112) fbidiag
            call c_writebidiag(ibyteswap,N_points,M_kers,N_iter,bidiag,
     c           maxpts, hhvec, maxnlm+2, U,N_points, V(2), M_kers,
     c           fbidiag)
         endif
c     
      else
c     
c     read factorization
c     
         write (6,*) 'Reading factorization from ',fbidiag

         call c_readbidiag(ibyteswap,N_points,M_kers,N_iter,bidiag,
     c        maxpts, hhvec, maxnlm+2, U,N_points, V(2), M_kers,
     c        fbidiag)
         write(6,*) 'N_points,M_kers, N_iter = ',N_points,M_kers,N_iter
c     test for same number of kernels and data
         if(M_kers.ne.ndata) then
            write(6,110) ndata,M_kers,fdata,fker,fbidiag
            stop 'incompatible bidiagonalization'
         endif
      endif
c     
c     Set up target functions:
c     Possibly scale target width in radius with sound speed
c     
      if (irscal.eq.1) then
         write (6,*) 'Scaling target width in radius with sound speed.'
         do i=1,N_targets
            call scale_rwidth(ibyteswap,famdl,targetparms(i,1),xrscal,
     c           targetparms(i,3))
         enddo
      endif
      write (6,*) 'Allocating target'
c     
c     Begin reading or computing transformed target functions
c     
      if(irdtranstarget.eq.0)  then
         write (*,*) 'Setting up targets'
         call startnewtimer(tnum1)
c
c     Set up targets
c
         call set_targets( targetparms,target,N_points,icase,
     c        itargettype)
         write(*,*) 'Time spent in set_targets = ',stoptimer(tnum1)
         write (6,*) N_targets,' target functions have been set up'
         call set_COG(target,N_points,cogparms)
         if (iwrtarget.eq.1) then
            call c_writeavker(ibyteswap, N_targets, N_points, icase, 
     c           trdoff, targetparms, maxtargets, cogparms(1,1), 
     c           cogparms(1,2), depart, target,N_points, favker)
         endif
c     
c     Compute transformed targets
c     
         call startnewtimer(tnum1)
         call set_transtarget(icase,sqrtw,target,N_points,akh1,U,
     c        N_points, work,lwork)
         write(*,*) 'Time spent in set_transtargets = ',
     c        stoptimer(tnum1)
         
         if(iwrtranstarget.eq.1) then
            write(6,113) ftranstarget
            call c_wrmatr(ibyteswap,N_iter,N_targets,target,N_points, 
     c           ftranstarget)
         endif
c     
      else
c     
c     Read transformed targets
c     
         write (6,*) 'Reading transformed targets from ',ftranstarget
         call c_rdmatr(ibyteswap,idummy,idummy1,target,N_points, 
     c        ftranstarget)
      endif
c
c     Compute the RLS solution to the damped bidiagonal system
c     
      write (*,*) 'Solving bidiagonal system.'
      call startnewtimer(tnum1)
      call solve_bidiag(target,N_points,bidiag,trdoff,depart,
     c     sigma_rot, alpha,icase)
      write(*,*) 'Time spent in solve_bidiag = ',stoptimer(tnum1) 
c     
c     Transform splittings according to the bidiagonalization
c     and compute rotation rate by transforming the solution
c     to the reduced system back to the original space.
c     
      call startnewtimer(tnum1)
      call transform_rhs(splits,V(2),M_kers) 
      call set_rot(icase,rot,target,N_points,splits)
      write(*,*) 'rot=',(rot(i),i=1,min(10,N_targets))
      write(*,*) 'Time spent in set_rot = ',stoptimer(tnum1)
c
c     Compute solution error correlation
c      
      if (icalccorr.eq.1) then
         write (*,*) 'Calculating solution error correlation.'
         call startnewtimer(tnum1)
         call set_correlation(fcorr, icase, trdoff, 
     c        targetparms, rot, sigma_rot, depart, target, N_points)
         write(*,*) 'Time spent in set_correlation = ',
     c        stoptimer(tnum1)
      endif
c
c     Compute inversion coefficients
c
      if (icalccoeff.eq.1) then
         write (*,*) 'Calculating inversion coefficients.'
         call startnewtimer(tnum1)
         call set_coeff_disk(fcoeff, icase, trdoff, targetparms, rot, 
     c        sigma_rot, depart, V,M_kers,target,N_points)
         write(*,*) 'Time spent in set_coeff = ',stoptimer(tnum1)         
      endif

      if (icalcavker.eq.1) then
c
c     Compute averaging kernels
c
         call startnewtimer(tnum1)
         call set_avker(icase, U, N_points, bidiag, akh1,
     c        target, N_points, sqrtw)
         write(*,*) 'Time spent in set_avker = ',stoptimer(tnum1)

         call startnewtimer(tnum1)
         write (*,*) 'Writing averaging kernels to ',favker
         call c_writeavker(ibyteswap, N_targets, N_points, icase,
     c        trdoff, targetparms, maxtargets, rot, sigma_rot,
     c        depart, target,N_points, favker)
c     
c     If averaging kernels were calculated then determine actual 
c     target positions and widths.
c     
         write (*,*) 'Calculating center of gravity ' // 
     c        'of averaging kernels.'
         call startnewtimer(tnum1)
         call set_COG(target,N_points,cogparms) 
         write(*,*) 'Time spent in set_COG = ',stoptimer(tnum1)
      endif
c     
c     Output inferred rotation rate with formal error bars
c     and measured widths of averaging kernels
c
      if (icalcsol.eq.1) then
         write (*,*) 'Writing solution to ',fsol
         write(ifsol,115) fdata, fker
         do i=1,N_targets
            write(ifsol,130) trdoff(i),(targetparms(i,j),j=1,4), 
     c           rot(i), sigma_rot(i)
            if (icalcavker.eq.1) then
               write(ifsol,131) cogparms(i,1),cogparms(i,2),
     c              (fitparms(i,j),j=1,4)
            else
               write(ifsol,131) -1d0,-1d0,-1d0,-1d0,-1d0,-1d0
            endif
         enddo
         call flush(6)
      endif
      write(*,*) 'Normal termination'
      write(*,*) 'Total execution time = ',stoptimer(tnum),' seconds.'

 110  format(//' ***** Error in 2D SOLA: ndata = ',i7,
     *     ' ne M_kers =',i7/
     *     '       Data   file: ',a60/
     *     '       kernel file: ',a60/
     *     '       bidiag file: ',a60/
     *     '       Execution terminated')
 112  format('/ Write bidiag factorization matrix to file'/
     *     2x,a77)
 113  format('/ Write transformed targets matrix to file'/
     *     2x,a77)
 115  format('# 2D SOLA results'/'#'/'# Data   file:'/'# ',a75/
     *     '# Kernel file:'/'# ',a75/'#'/
     *     '# trade-off, r0, theta0, deltar, deltat, Omega/2 pi,',
     *     'error'/
     *     '# r_cg, theta_cg, r_fit, theta_fit, deltar_fit, ',
     *     'deltat_fit:'/'#')
 130  format(1pe11.3,0p4f10.5,1p2e13.5)
 131  format(6e13.5)
      end
