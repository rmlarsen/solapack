c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine set2dker(afile,nkers,lker,nker,mker,itker,it1ker
     c     ,cfile,npoints,rpoints,tpoints,weights)

c     procedure to set up the A matrices

      implicit double precision (a-h,o-z)
      include 'include1.n.f'
      include 'include2.f'

      integer lker(maxnlm),nker(maxnlm),mker(maxnlm)
      integer itker(maxnlm),it1ker(maxnlm),npoints
      character*80 cfile, afile
      real rpoints(*),tpoints(*)
      real weights(*)

      double precision tmesh1(maxt1),tweights(maxt1)
      double precision rmesh1(maxr1),rweights(maxr1)
      character*80 cefunc
      real dummy(2),dtime,etime,time1
      real f1(maxmodes,maxr1),f2(maxmodes,maxr1)
      real g1(maxlm,maxt1),g2(maxlm,maxt1)
      integer nnlm1(maxnlm),llist(maxnlm)
      logical found
      integer ijk(maxr1,maxt1)
      common/fgcom/f1,f2,g1,g2

c     nmode,lmode: identification of modes
c     freq: frequency of mode (muHz)
c     radmesh: mesh in radius for kernels and eigenfunctions
c     cefunc: control file for reading of eigenfunctions

c     read parameters etc.
      write (6,*) 'reading parameters from file'
      write (6,'(A80)') cfile
      write (6,*) '------------------'
      call catfile(cfile)
      write (6,*) '------------------'
      open (11,file=cfile,form='formatted')
      call clnfile(11,10)
      close(11)
      write (6,*) 'control file for reading of eigenfunctions'
      read (10,'(A80)') cefunc
      write (6,*) 'number of points to skip in r for setting up f'
      read (10,*) nskipr
      write (6,*) 'number of points in x for setting up the plm''s'
      read (10,*) nsetp
      write (6,*) 'number of points to skip in x for setting up g'
      read (10,*) nskipt
      write (6,*) 'irdf, iwrf, irdg, iwrg'
      read (10,*) irdf,iwrf,irdg,iwrg
c     write (6,*) '0 for piecewise linear omega,'
c     write (6,*) '1 for expansion in cos^2n theta'
c     read (10,*) icaseb
      icaseb=0
      close(10)

      nnlm=nkers
      do 5,i=1,nnlm
        lnlm(i)=lker(i)
        nnlm1(i)=nker(i)
        mnlm(i)=mker(i)
        itnlm(i)=itker(i)
        it1nlm(i)=it1ker(i)
 5    continue
      call setilm
      lmax1=-1
      do 10,i=1,nnlm
        if (lnlm(i).gt.lmax1) lmax1=lnlm(i)
 10   continue
      i1=0
      do 80,l=0,lmax1
        nl=0
        do 20,i=1,nnlm
          if (lnlm(i).eq.l) then
            nl=nl+1
            llist(nl)=i
          end if
 20     continue
        i2=i1+1
        do 60,i=1,nl
          il=llist(i)
          found=.false.
          do 30,ii=i2,i1
            if (nnlm1(il).eq.nmode(ii)) then
              found=.true.
              ifound=ii
            end if
  30      continue
          if (found) then
            iln(il)=ifound
          else
            i1=i1+1
            lmode(i1)=l
            nmode(i1)=nnlm1(il)
c           freq(i1)=fnlm(il)
            iln(il)=i1
          end if
 60     continue
 80   continue
      nmodes=i1
      write (6,*) nmodes,nnlm
      pi=4.0D0*atan(1.0D0)
      pi2=pi/2.0D0
      write (6,*) 'start reading eigenfunctions'
      time1=dtime(dummy)
      call readefun(cefunc)
      write (6,*) 'eigenfunctions have been read, time=',
     c            dtime(dummy),dummy(1),dummy(2)
      nrad1=(nrad-1)/nskipr+1
      write (*,*) 'nrad,nskipr,nrad1 = ',nrad,nskipr,nrad1
c
c     Set up F matrices.
c
      write (6,*) 'start setting up f matrix'
      write (6,*) '#nl, #rpoints = ',nmodes,nrad1 
      time1=dtime(dummy)
      call setf(nskipr,f1,f2)
      write (6,*) 'f matrix has been set up, time=',
     c     dtime(dummy),dummy(1),dummy(2)

      nsetp1=(nsetp-1)/nskipt+1
c
c     Set up G matrices.
c
      write (6,*) 'start setting up g matrix'
      write (6,*) '#lm, #tpoints = ',nlm,nsetp1 
      time1=dtime(dummy)
      call setg(icaseb,nsetp,nsetp1,g1,g2,tmesh1)
      write (6,*) 'g matrix has been set up, time=',
     c     dtime(dummy),dummy(1),dummy(2)

c
c     Set up grid and integration weights.
c
      time1=dtime(dummy)
      do 100,i=1,nrad
        dradmesh(i)=radmesh(i)
 100  continue
      do 110,i=1,nrad1
        rmesh1(i)=dradmesh(nskipr*(i-1)+1)
 110  continue
      call dsetint(1,nrad1,rmesh1,rweights)
      call dsetint(1,nsetp1,tmesh1,tweights)
      npoints=nrad1*nsetp1
      i=0
      do j=1,nrad1
        do k=1,nsetp1
          i=i+1
          ijk(j,k)=i
          rpoints(i)=rmesh1(j)
          tpoints(i)=tmesh1(k)
          weights(i)=rweights(j)*tweights(k)
        end do
      end do

      
c
c     Write modeinfo to 'afile'
c
      open (21,file=afile,form='unformatted',
     c     status='unknown')

      write (21) nnlm,npoints
      write (21) nnlm,(iln(i),i=1,nnlm),(ilm(i),i=1,nnlm)
c
c     Write F matrices to 'afile'
c
      write (6,*) 'start writing f matrix'
      time1=dtime(dummy)
      write (21) nrad1,nmodes
      write (21) ((f1(i,j),j=1,nrad1),i=1,nmodes)
      write (21) ((f2(i,j),j=1,nrad1),i=1,nmodes)
      write (6,*) 'f matrix has been written, time=',
     c     dtime(dummy),dummy(1),dummy(2)
c
c     Write G matrices to 'afile'
c
      write (6,*) 'start writing g matrix'
      time1=dtime(dummy)
      write (21) nsetp1,nlm
      write (21) ((g1(i,j),j=1,nsetp1),i=1,nlm)
      write (21) ((g2(i,j),j=1,nsetp1),i=1,nlm)
      close(21)
        write (6,*) 'g matrix has been written, time=',
     c              dtime(dummy),dummy(1),dummy(2)

      write (6,*) 'total time used:',etime(dummy)
      write (6,*) 'user time:      ',dummy(1)
      write (6,*) 'system time:    ',dummy(2)
      write(31,*) nsetp1, nrad1, 
     *  '    output no. of points in theta and r'
      write(31,*) nsetp, nrad, '    orig. no. of points in theta and r'
      end

