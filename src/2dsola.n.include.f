c     Copyright Rasmus Munk Larsen, Stanford University, 2003
c	
c     Parameters defining the maximal problem size:
c     
c     maxrad:         max number of discretization points in radius.
c     maxtheta:       max number of discretization points in latitude.
c     maxpts:         max total number of discretization points in (r,theta).
c     maxl:           maximal degree of mode.
c     maxlm:          max number of (l,m) combinations in modeset.
c     maxnl/maxmodes: max number of (n,l) combinations (multiplets) in modeset.
c     maxnlm:         max total number of modes (n,l,m).
c
      integer maxrad,maxtheta,maxpts,maxnlm,maxtargets,maxsplits
      integer maxiter,maxmodes,maxnl,maxlm,maxl,maxexp
      parameter(maxrad = 103,maxtheta = 53, maxpts = maxrad*maxtheta)
      parameter(maxmodes=5001,maxl=301,maxnl=maxmodes,maxexp=36)
      parameter(maxlm=(maxl+1)*(maxl+2)/2,maxnlm=45001)
      parameter(maxtargets=maxpts+2, maxsplits=1, maxiter=2000)
c
c     Parameters for L-curve analysis.
c     maxtrdoff: max number of points on the L-curve.
c
      integer maxtrdoff
      parameter(maxtrdoff=100)
c
c     Problem dimensions.
c
      integer M_kers,M_nl,M_lm
      integer N_points,N_rad,N_theta
      integer N_targets,N_iter
      common/cdim/M_kers,M_nl,M_lm,N_points,N_rad,
     c     N_theta,N_targets,N_iter

c
c     Data defining the kernel matrix, mesh and modeset.
c
c     modeset       :  
c     rthetw        :  Mesh points and quadrature (integration) weights.
c                      The i'th meshpoint is
c                         (r,theta) = (rthetw(i,1),rthetw(i,2)) 
c                      and rthetw(i,3) gives the quadrature eights for this
c                      point.
c     f1, f2        :  Matrices describing the radial behavior 
c                      of the kernels.
c     inl           :  Index array for f1, f2. The radial behavior of the 
c                      i'th kernel is described by f1(*,inl(i)) and
c                      f2(*,inl(i)).
c     g1, g2        :  Matrices describing the latitudinal behavior 
c                      of the kernels.
c     ilm           :  Index array for g1, g2. The latitudinal behavior of the 
c                      i'th kernel is described by g1(*,ilm(i)) and
c                      g2(*,ilm(i)).
c     chol_sigma    :  chol_sigma(*,*,i), i=1,...,_nl is the Cholesky factor of the 
c                      error-covariance matrix for the i'th multiplet.
c     invsigma_split:  The reciprocal of the error estimates of the splittings. 
c                      (= The reciprocal of the sqrt of the diagonal elements of chol_sigma.)
c     invsqrtw      :  The reciprocal of the square root of the quadrature
c                      (integration) weights.
c     ilblock       :  [ilblock(i):ilblock(i+1)-1] is the i'th block
c                      of modes with the same value of l
c     nlblocks      :  Gives the number of such blocks, i.e. the length 
c                      of ilblock.
c     inblock       :  [inblock(i):inblock(i+1)-1} is the block consisting of 
c                      modes from the i'th multiplet.
c     nnblocks      :  The number of blocks in inblock, must be equal to M_nl.
c     iblock        :  [iblock(i,j):iblock(i,j+1)-1] is the j'th block
c                      of modes with the same number of splittings/a-coefs in
c                      the i'th block of modes with the same value of l.
c                      In other words, iblock(i,*) describes at subdivision 
c                      the block [ilblock(i):ilblock(i+1)-1]
c     nblocks       :  nblocks(i) gives the number of subdivisions of 
c                      [ilblock(i):ilblock(i+1)-1],  i.e. the length 
c                      of iblock(i,*).
c
      double precision modeset(maxnlm,8),rthetw(maxpts,3)
      double precision f1(maxrad,maxmodes),f2(maxrad,maxmodes)
      double precision g1(maxtheta,maxlm),g2(maxtheta,maxlm)
      double precision invsigma_split(maxnlm),invsqrtw(maxpts)
      double precision hhvec(maxnlm+2),Dw(maxpts)
      double precision rmesh(maxrad),tmesh(maxtheta)
      double precision invsigma_nl(maxmodes),invsigma_lm(maxlm)
      double precision chol_sigma(maxexp,maxexp,maxmodes)
      integer icovar, inl(maxnlm),ilm(maxnlm), ilblock(10*maxl)
      integer nlblocks,inblock(maxnlm), nnblocks,ibyteswap
      integer iblock(2*maxl,10*maxl),nblocks(10*maxl)

      common/cfg/invsigma_nl,invsigma_lm,Dw,rmesh,tmesh,modeset,rthetw,
     c     hhvec,invsqrtw,chol_sigma,invsigma_split,f1,f2,g1,g2,ilm,
     c     inl,iblock,nblocks,inblock,nnblocks,ilblock,nlblocks
      common/comctrl/icovar,ibyteswap
c     
c     Common workspace.
c      
      integer lwork
      parameter(lwork = 64*maxnlm)
      double precision work(lwork)
c
      common/cwrk/work

c
c     I/O unit numbers for data files.
c
      integer iflog,ifker, ifmesh, ifdata, ifcovar, ifsol, ifavker, 
     *     ifcoeff, ifbidiag, iftranstarget, iftransrhs, ifamdl,
     *     iflcurve, ifres,ifrestrict,ifconv,ifcorr
      parameter(iflog=27, ifker=2, ifmesh=3, ifdata=4, ifsol=10)
      parameter(ifavker=11, ifcoeff=12, ifres=13, ifconv=14)
      parameter(ifrestrict=15,ifcorr=16,ifbidiag=21, iftranstarget=22)
      parameter(iftransrhs=23, iflcurve=24, ifamdl=25, ifcovar=26)
c
c     Directory and logical unit number for temporary files
      integer iftmp
      parameter (iftmp=27)
      character*(*) tmpdir
      parameter (tmpdir='/tmp')
