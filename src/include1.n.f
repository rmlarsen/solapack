c     Copyright Rasmus Munk Larsen, Stanford University, 2003

c     parameter (maxrad=602,maxr0=maxrad,maxmodes=1500,maxl=260)
c     parameter (maxrad=602,maxr0=maxrad,maxmodes=1000,maxl=65)
      parameter (maxrad=602,maxr0=maxrad,maxmodes=4000,maxl=310)
      parameter (maxrrd=5000)
      parameter (maxexp=36)
      character*(*) tmpdir
      parameter (tmpdir='/tmp')
      character*24 timestamp

c     maxrad: maximum number of gridpoints for kernels
c     maxr0: maximum number of target radii for averaging kernels
c            and coeffecients
c     maxmodes: maximum number of modes (n,l) to be used
c     maxl: highest l to be used
c     maxrrd: maximum no. of mesh points allowed in eigenfunctions
c             on file
c     maxexp: maximum number of expansion coefficients (theta,m)

      integer nmode(maxmodes),lmode(maxmodes)
      real freq(maxmodes),radmesh(maxrad)
      double precision dradmesh(maxrad)
      character*80 suff

      common/ccase/iexpcase,nexp
      common/cid/nmodes,nmode,lmode,freq
      common/commesh/dradmesh,radmesh,nrad
      common/csuff/suff
      common/ctime/timestamp

c     nmode,lmode: identification of modes
c     freq: frequency of mode (muHz)
c     iexpcase: case of expansions
c            0 for latitude independent rotation
c            1 for a sub i (m/l) and orthogonal polynomials
c            2 for a sub i (m/L) and orthogonal polynomials
c            3 for a sub i (truly orthogonal) and orthogonal polynomials
c     nexp: number of expansion coefficients
c     suff: suffix to add to filenames
