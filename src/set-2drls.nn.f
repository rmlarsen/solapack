c     Copyright Rasmus Munk Larsen, Stanford University, 2003


      program test2d
c
c  set up and output A matrix for 2-D inversion
c  Note A matrix is output in single precision
c
      implicit double precision (a-h,o-z)
      character*256 cfile,fsplit, afile, frthet, fsum
      include 'include1.n.f'
      include 'include2.f'

      parameter (maxpoints=maxr1*maxt1)

      integer lker(maxnlm),nker(maxnlm),mker(maxnlm),itker(maxnlm),
     *        it1ker(maxnlm)
      real frqker(maxnlm), splker(maxnlm), errker(maxnlm)
      double precision rpoints(maxpoints),tpoints(maxpoints),
     c     weights(maxpoints)

      common/incom/lker,nker,mker,itker,it1ker,frqker, splker, errker,
     c     rpoints,tpoints,weights
c
      write(6,*) 'maxpoints, maxnlm =',maxpoints, maxnlm
c
c  set up file names: 
c  control file (note this contains file name for eigenfunction
c                control file)
c  file of splittings
c  file for A matrix output
c  file for r and theta point and weight output
c  file for summary of mode set and indices
c
      read(5,'(a)') cfile
      read(5,'(a)') fsplit
      read(5,'(a)') afile
      read(5,'(a)') frthet
      read(5,'(a)') fsum
c
c  enter 1 to scale with errors, 0 not to scale
c
      read(5,*) ierrfl
c
      open(31,file=fsum,status='unknown')
      write(31,*) 'mode file: ',fsplit(1:60)
      write(31,*) ierrfl, '   ierrfl'
c
      open (10,file=fsplit,status='old',form='formatted')
c
      i=1
   10 read (10,*,end=15) lker(i),nker(i),frqker(i),mker(i),
     *                   itker(i),it1ker(i),splker(i), errker(i)
      i=i+1
      go to 10
c
   15 nkers=i-1
      write(6,110) nkers
      write(31,*) nkers,'   nkers'
      close (10)
      call set2dker(afile,nkers,lker,nker,mker,itker,it1ker,cfile,
     c              npoints,rpoints,tpoints,weights)

c     
      open (10,file=frthet,status='unknown',form='unformatted')
      write (10) npoints
      write (10) (rpoints(j),j=1,npoints),
     *     (tpoints(j),j=1,npoints),
     *     (weights(j),j=1,npoints)
      close (10)
      write (6,*) 'nkers, npoints = ', nkers, npoints
      close(31)
      stop
  110 format(//' No of modes read in:',i7)

      end
