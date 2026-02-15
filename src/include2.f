c     Copyright Rasmus Munk Larsen, Stanford University, 2003

c     extra parameters for 2d inversions

      parameter (maxr1=402,maxt=1251,maxt1=202)
c     parameter (maxr1=102,maxt=1251,maxt1=42)
      parameter (maxlm=(maxl+1)*(maxl+2)/2,maxnlm=50001)
c..      parameter (maxlm=(maxl+1)*(maxl+2)/2,maxnlm=3100)
      parameter (maxa=maxexp+36,maxla=5*maxl)
      parameter (maxrin=210,maxtin=60)
      parameter (maxrav=220,maxtav=220)
c     parameter (maxrav=102,maxtav=101)

c     maxr1: maximum number of radial gridpoints for inversions
c     maxt: maximum number of gridpoints in theta for setting up plm's
c     maxt1: maximum number of gridpoints in theta for inversion
c     maxlm: maximum number of l,m values
c     maxnlm: maximum total number of modes (m's and a's)
c     maxa: maximum number of a coefficients to use (even exclusive)
c     maxla: maximum number of combinations of l's and maximum a's,
c            most often maxla=maxl. Must be set higher if different
c            types of a-coeffecients are mixed.
c     maxrin,maxtin: maximum number of points where omega is given
c     maxrav,maxtav: maximum number of points to set averaging kernels

      integer llm(maxlm),mlm(maxlm),itlm(maxlm),it1lm(maxlm)
      integer lnlm(maxnlm),mnlm(maxnlm),itnlm(maxnlm),it1nlm(maxnlm)
      integer ilm(maxnlm),iln(maxnlm)

      common/clm/lmax,nlm,nnlm,llm,mlm,itlm,it1lm,
     c           lnlm,mnlm,itnlm,it1nlm,ilm,iln

c     lmax: maximum l set
c     nlm: number of (l,m)'s set
c     nnlm: total number of splittings
c     llm(i) and mlm(i): l and (m or k) of the i'th G
c     itlm(i): type of the i'th G, it=0 for m, it=1 for ak(m/l),
c              it=2 for ak(m/L) and it=3 for orthogonal
c     it1lm(i): max number of k's used for the i'th G
c     the nlm arrays are like the lm arrays but for the i'th splitting
c     ilm(i): number in G list of the i'th splitting
c     iln(i): number in mode list of the i'th splitting
c     xnlm(i)=xlm(ilm(i)), must hold for i=1,nnlm
