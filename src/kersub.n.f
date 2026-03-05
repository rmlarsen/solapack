c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine readefun(cfile)

c     reads eigenfunctions from one or more .z files

      include 'include1.n.f'

      real efuncr(maxrad,maxmodes),efunct(maxrad,maxmodes)
      real rm1(maxrrd),rho(maxrrd)
      character*80 cfile,filename(10),lnfile,modelfile
      real cs(50),oldmesh(maxrrd),func(2,maxrrd),strmesh(maxrrd),
     *  orgmesh(maxrrd),strfunc(2,maxrrd)
      double precision dcs(50),doldmesh(maxrrd),dfunc(2,maxrrd)
      logical found(maxmodes),lstop

      common/cfunc/efuncr,efunct

c     nmode,lmode: identification of modes
c     freq: frequency of mode (muHz)
c     radmesh: Mesh in radius for kernels and eigenfunctions
c     efuncr and efunct: Eigenfunctions as read from .z file
c     cfile: control file
c     iselect: 0 for reading the eigenfunctions corresponding
c                to lmode(i),nmode(i),i=1 to nmodes
c              1 for reading all kernels with lmin<=l<=lmax,...
c              2 for reading list of l,n from file
c     isetfrq: 0 for not setting frequencies, 1 for setting them
c                (only sets for iselect=0 or 2)
c     lmin,lmax,nmin,nmax,fmin,fmax: window in l,n,f for iselect=1
c     lnfile: name of file to read l,n list from
c     nskip: only every nskip radial meshpoints from the original
c            eigenfunctions are used
c     nrsmsh:  for nrsmsh = 1, when reading from two or more files,
c               reset meshfor subsequent files based on mesh 
c               for initial file
c               (Note: this has not been done consistently
c               for itype = 2, where model file must be read)
c     itype: 0 for .z files with real numbers
c            1 for .z files with double precision numbers
c            2 for GONG type files
c     iadd: 0 for not adding extra gridpoint
c           1 for adding extra gridpoint

      write (6,*) 'reading parameters from file'
      write (6,'(A80)') cfile
      write (6,*) '--------------------'
      call catfile(cfile)
      write (6,*) '--------------------'
      open (10,file=cfile,form='formatted')
      call clnfile(10,11)
      close(10)
      read (11,*) itype,iadd
      if ((iadd.ne.0).and.(iadd.ne.1)) stop
      read (11,*) iselect, isetfrq
      read (11,*) lmin,lmax,nmin,nmax,fmin,fmax
      read (11,*) lnfile
      read (11,*) nskip, nrsmsh
      read (11,*) modelfile
      read (11,*) nfiles
      do 10,i=1,nfiles
        read (11,'(A80)') filename(i)
 10   continue
      close(11)
      if (iselect.eq.2) then
c       for consistency
        isetfrq=1
        nmodes=0
        open (10,file=lnfile,form='formatted')
 11     read (10,*,end=13) lhelp,nhelp
        nmodes=nmodes+1
        if (nmodes.gt.maxmodes) then
          write (6,*) 'too many modes, only dimensioned for' //
     c          'maxmodes=', maxmodes,' modes'
          close(10)
          stop
        end if
        lmode(nmodes)=lhelp
        nmode(nmodes)=nhelp
        goto 11
 13     continue
        close(10)
c       from now on pretend iselect=0
        iselect=0
      end if
      if (iselect.eq.0) then
        do 15,i=1,nmodes
          found(i)=.false.
 15     continue
      else
        nmodes=0
      end if
      nrst1=0
      write (6,*) 'reading eigenfunctions from file(s):'
      do 50,nfile=1,nfiles
        write (6,'(A80)') filename(nfile)
        if (iselect.eq.1) then
          write (6,*) 'modes selected:'
          write (6,*)
          write (6,*) '  i,  l,  n,  f'
        end if
        if (itype.eq.0) then
          open (10,file=filename(nfile),form='unformatted',
     c          status='old')
          read (10) nr1,(oldmesh(i),i=1,nr1)
        else if (itype.eq.1) then
          open (10,file=filename(nfile),form='unformatted',
     c          status='old')
          read (10) nr1,(doldmesh(i),i=1,nr1)
          write(6,*) 'nr1,nmodes = ',nr1,nmodes
          if (nr1.gt.maxrrd) then
             write (6,*) 'too many mesh points in eigenfunctions,'
             write (6,*) 'asking for nr1 = ',nr1
             write (6,*) 'only dimensioned for maxrrd = ',maxrrd
             close(10)
             stop
          end if
          do i=1,nr1
            oldmesh(i)=doldmesh(i)
          end do
        else
          open (10,file=modelfile,form='formatted',status='old')
          read (10,*) nr1
          write(6,*) 'nr1,nmodes = ',nr1,nmodes
          if (nr1.gt.maxrrd) then
             write (6,*) 'too many mesh points in eigenfunctions,'
             write (6,*) 'asking for nr1 = ',nr1
             write (6,*) 'only dimensioned for maxrrd = ',maxrrd
             close(10)
             stop
          end if
          read (10,*) (dum,oldmesh(i),rho(i),dum,i=1,nr1)
          close(10)
          open (10,file=filename(nfile),form='formatted',status='old')
        end if
c
c  test for storing mesh for later interpolation
c
        nrorg1=nr1
        if(nrsmsh.eq.1) then
          if(nfile.eq.1) then
            do i=1,nr1
              strmesh(i)=oldmesh(i)
            end do
            nrst1=nr1
          else
            do i=1,nr1
              orgmesh(i)=oldmesh(i)
            end do
            do i=1,nrst1
              oldmesh(i)=strmesh(i)
            end do
            nr1=nrst1
          end if
        end if
        nr2=nr1
        if (iadd.eq.1) then
          do i=nr1,2,-1
            oldmesh(i+1)=oldmesh(i)
            rho(i+1)=rho(i)
          end do
          oldmesh(2)=0.5*oldmesh(3)
          rho(2)=rho(1)
          nr2=nr1+1
        end if
        nrad=(nr2-1)/nskip+1
        if (nrad.gt.maxrad) then
          write (6,*) 'too many meshpoints'
          write (6,*) 'asked for',nrad,' meshpoints but only'
          write (6,*) 'dimensioned for maxrad=',maxrad,' meshpoints'
          stop
        end if
        do 17,i=0,nrad-1
          radmesh(i+1)=oldmesh(nskip*i+1)
          dradmesh(i+1)=oldmesh(nskip*i+1)
 17     continue
        if ((itype.eq.0).or.(itype.eq.1)) then
          do 18,i=1,nrad
            if (radmesh(i).eq.0.0) then
              rm1(i)=0.0
            else
              rm1(i)=1.0/sqrt(radmesh(i))
            end if
 18       continue
        else
          do 19,i=1,nrad
            rm1(i)=radmesh(i)*sqrt(rho(i))
 19       continue
        end if
 20     continue
        if (itype.eq.0) then
          read (10,end=40) cs,((func(i,n),i=1,2),n=1,nrorg1)
          l=nint(cs(18))
          n=nint(cs(19))
          f=1000*cs(27)
        else if (itype.eq.1) then
          read (10,end=40) dcs,((dfunc(i,n),i=1,2),n=1,nrorg1)
          do i=1,nrorg1
            func(1,i)=dfunc(1,i)
            func(2,i)=dfunc(2,i)
          end do
          l=nint(dcs(18))
          n=nint(dcs(19))
          f=1000*dcs(27)
        else
          read (10,*,end=40) l,n,f
          f=1000*f
          read (10,*) ((func(i,j),i=1,2),j=1,nrorg1)
        end if
c
c  test for resetting mesh
c
        if(nrsmsh.eq.1.and.nfile.gt.1) then
          do i=1,2
            do j=1,nrorg1
              strfunc(i,j)=func(i,j)
            end do
          end do
          do i=1,nr1
            call lir(strmesh(i),orgmesh,func(1,i),strfunc,2,2,
     *            nrorg1,i,inter)
          end do
        end if
        p=max(l-1,1)
        if (iadd.eq.1) then
          do i=nr1,2,-1
            func(1,i+1)=func(1,i)
            func(2,i+1)=func(2,i)
          end do
          func(1,2)=(0.5**p)*func(1,3)
          func(2,2)=(0.5**p)*func(2,3)
        end if
        if (iselect.eq.0) then
          do 30,i=1,nmodes
            if ((lmode(i).eq.l).and.(nmode(i).eq.n)) then
              found(i)=.true.
              do 31,j=0,nrad-1
                efuncr(j+1,i)=func(1,nskip*j+1)*rm1(j+1)
                efunct(j+1,i)=func(2,nskip*j+1)*rm1(j+1)
 31           continue
              if (isetfrq.eq.1) freq(i)=f
              goto 30
            end if
 30       continue
        else
          if ((l.ge.lmin).and.(l.le.lmax).and.
     c        (n.ge.nmin).and.(n.le.nmax).and.
     c        (f.ge.fmin).and.(f.le.fmax)) then
            nmodes=nmodes+1
            write (6,'(3i4,f7.1)') nmodes,l,n,f
            lmode(nmodes)=l
            nmode(nmodes)=n
            freq(nmodes)=f
            do 35,j=0,nrad-1
              efuncr(j+1,nmodes)=func(1,nskip*j+1)*rm1(j+1)
              efunct(j+1,nmodes)=func(2,nskip*j+1)*rm1(j+1)
 35         continue
          end if
        end if
        goto 20
 40     close(10)
 50   continue
      write (6,*) nr1,nr2,nrad,nskip
      if (iselect.eq.0) then
        lstop=.false.
        do 60,i=1,nmodes
          if (.not.found(i)) then
            write (6,*) 'missing mode with (l,n)=',lmode(i),nmode(i)
            lstop=.true.
          end if
 60     continue
        if (lstop) stop
      end if
      write(6,*) 'exiting from readefunc'
c     do j=1,nrad
c       write (6,*) radmesh(j),
c    c              (efuncr(j,i),i=1,nmodes),(efunct(j,i),i=1,nmodes)
c     end do
      return
      end

      subroutine setf(nskipr,f1,f2)

c     sets up integrals of eigenfunctions and basefunctions for
c     2d calculations

      implicit double precision (a-h,o-z)
      include 'include1.n.f'
      include 'include2.f'

      real efuncr(maxrad,maxmodes),efunct(maxrad,maxmodes)
      double precision drmesh(maxrad),rmesh1(maxr1),weights(maxrad)
      double precision k1(0:maxrad),k2(0:maxrad)
      double precision w1(-maxrad:maxrad,maxr1),whelp(0:maxrad)
      real f1(maxmodes,maxr1),f2(maxmodes,maxr1)

      common/cfunc/efuncr,efunct
      common/setfcom/w1,whelp

      do 10,i=1,nrad
        drmesh(i)=radmesh(i)
 10   continue
      nrad1=(nrad-1)/nskipr+1
      do 15,i=1,nrad1
        rmesh1(i)=drmesh(nskipr*(i-1)+1)
 15   continue
      call dsetint(1,nrad,drmesh,weights)
      drmesh(1)=0.0
c     do 20,i=2,nrad
c       weights(i)=weights(i)/drmesh(i)
c20   continue
      sum=0.0D0
      do 80,i=1,nrad1
        do 30,j=-maxrad,maxrad
          w1(j,i)=0.0D0
 30     continue
        i0=nskipr*(i-1)+1
        r0=drmesh(i0)
        if (i.gt.1) then
          call dsetint(1,nskipr+1,drmesh(i0-nskipr),whelp)
          d=1.0D0/(r0-drmesh(i0-nskipr))
          do 40,j=1-nskipr,0
            w1(j,i)=w1(j,i)+(1.0D0-(r0-drmesh(i0+j))*d)*whelp(j+nskipr)
 40       continue
        end if
        if (i.lt.nrad1) then
          call dsetint(1,nskipr+1,drmesh(i0),whelp)
          d=1.0D0/(drmesh(i0+nskipr)-drmesh(i0))
          do 60,j=0,nskipr-1
            w1(j,i)=w1(j,i)+(1.0D0-(drmesh(i0+j)-r0)*d)*whelp(j)
 60       continue
        end if
        do 70,j=-nskipr,nskipr
          sum=sum+w1(j,i)
 70     continue
 80   continue
      write (6,*) 'integration error in r',sum-(drmesh(nrad)-drmesh(1))
      do 200,i=1,nmodes
        factor=0.0D0
        do 110,j=2,nrad
          factor=factor+weights(j)*(efuncr(j,i)**2+efunct(j,i)**2)
 110    continue
        factor=1.0D0/factor
        l=lmode(i)
        capl2=l*(l+1.0)
        c1=2.0/sqrt(capl2)
        c2=factor/capl2
        k1(1)=0.0D0
        k2(1)=0.0D0
        do 120,j=2,nrad
          k1(j)=(efuncr(j,i)*(efuncr(j,i)-c1*efunct(j,i)))*
c    c          factor/drmesh(j)
     c          factor
c         k2(j)=c2*efunct(j,i)**2/drmesh(j)
          k2(j)=c2*efunct(j,i)**2
 120    continue
        sum1=0.0D0
        sum2=0.0D0
        do 150,j=1,nrad1
          j1=nskipr*(j-1)+1
          f1j=0.0D0
          f2j=0.0D0
          do 140,k=max(1-j1,-nskipr),min(nrad-j1,nskipr)
            f1j=f1j+w1(k,j)*k1(j1+k)
            f2j=f2j+w1(k,j)*k2(j1+k)
 140      continue
          f1(i,j)=f1j
          f2(i,j)=f2j
          sum1=sum1+f1j
          sum2=sum2+f2j
 150    continue
 200  continue

      end

      subroutine setilm

c     setls ilm, llm, mlm, itlm and it1lm from lnlm, mnlm, itnlm and it1nlm
c     such that xnlm(i)=xlm(ilm(i)) for i=1,nnlm
c     also sets lmax and nlm

      implicit double precision (a-h,o-z)
      include 'include1.n.f'
      include 'include2.f'

      integer llist(maxnlm)
      logical found
      
      lmax=-1
      do 10,i=1,nnlm
        if (lnlm(i).gt.lmax) lmax=lnlm(i)
 10   continue
      i1=0
      do 100,l=0,lmax
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
            if ((mnlm(il).eq.mlm(ii)).and.
     c          (itnlm(il).eq.itlm(ii)).and.
     c          (it1nlm(il).eq.it1lm(ii))) then
              found=.true.
              ifound=ii
            end if
  30      continue
          if (found) then
            ilm(il)=ifound
          else
            i1=i1+1
            llm(i1)=l
            mlm(i1)=mnlm(il)
            itlm(i1)=itnlm(il)
            it1lm(i1)=it1nlm(il)
            ilm(il)=i1
          end if
 60     continue
 100  continue
      nlm=i1

      end

      subroutine setg(icaseb,nsetp,nset1,g1,g2,tmesh1)

c     sets up integrals of plms and basefunctions for
c     2d calculations
c     represents the plms on nsetp points
c     for icaseb=0 use basefunctions for a piecewice linear function
c     with nset1 points in latitude (endpoints included)
c     for icaseb=1 use cos(theta)^(2(n-1)) as basefunctions for n=1,nset1
c     for icaseb=2 use orthogonal polynomials as basefunctions for n=1,nset1

      implicit double precision (a-h,o-z)
      include 'include1.n.f'
      include 'include2.f'

      double precision weights(maxt),tmesh(maxt),tmesh1(maxt1)
      double precision xmesh(maxt),xmesh1(maxt1),smesh(maxt)
      double precision smesh1(maxt),smesh3(maxt),sc2mesh(maxt)
      double precision acoeff(maxa,maxl,maxla)
      integer ila(maxl,maxa,3),llist(maxnlm)
      real g1(maxlm,maxt1),g2(maxlm,maxt1)
      double precision w1(-maxt:maxt,maxt1),whelp(0:maxt)
      double precision plm(maxt,0:maxl),dplm(maxt,0:maxl)
      double precision g1l(maxt),g2l(maxt)
      double precision c(0:maxexp-1,0:maxexp-1)
      integer indx(0:maxl),itmp
      logical setl(0:maxl)

      common/setgcom/w1,whelp,plm,dplm,g1l,g2l,acoeff

      do i=1,maxt1
         do j=1,maxlm
            g1(j,i) = 0e0
            g2(j,i) = 0e0
         enddo
      enddo
      pi=4.0D0*atan(1.0D0)
      pi2=pi/2.0D0
      do 10,i=1,nsetp
        tmesh(i)=pi2*(i-1)/(nsetp-1.0D0)
        xmesh(i)=cos(tmesh(i))
        smesh(i)=sin(tmesh(i))
        if (smesh(i).eq.0.0D0) then
          smesh1(i)=0.0D0
        else
          smesh1(i)=1.0D0/smesh(i)
        end if
        smesh3(i)=smesh(i)**3
        sc2mesh(i)=2.0D0*smesh(i)*xmesh(i)
 10   continue
c     nset1=(nsetp-1)/nskipt+1
      if (nset1.le.1) then
c got to set nskipt to something
        nskipt=1
      else
        nskipt=(nsetp-1)/(nset1-1)
        write (6,*) 'nset1 = ', nset1, ' nsetp = ',nsetp
        if ((mod(nset1-1,2).ne.0).or.(mod(nskipt,2).ne.0)) then
          write (6,*) 'both resulting number of latitude meshpoints-1'
          write (6,*) 'and skip should be divisible by 2'
          stop
        end if
      end if
      do 15,i=1,nset1
        tmesh1(i)=tmesh(nskipt*(i-1)+1)
        xmesh1(i)=xmesh(nskipt*(i-1)+1)
 15   continue
      call dsetint(2,nsetp,tmesh,weights)
      if (icaseb.eq.0) then
        sum=0.0D0
        do 80,i=1,nset1
          do 30,j=-maxt,maxt
            w1(j,i)=0.0D0
 30       continue
          i0=nskipt*(i-1)+1
          x0=tmesh(i0)
          if (i.gt.1) then
            call dsetint(2,nskipt+1,tmesh(i0-nskipt),whelp)
            d=1.0D0/(x0-tmesh(i0-nskipt))
            do 40,j=1-nskipt,0
              w1(j,i)=w1(j,i)+(1.0D0-(x0-tmesh(i0+j))*d)*whelp(j+nskipt)
 40         continue
          end if
          if (i.lt.nset1) then
            call dsetint(2,nskipt+1,tmesh(i0),whelp)
            d=1.0D0/(tmesh(i0+nskipt)-tmesh(i0))
            do 60,j=0,nskipt-1
              w1(j,i)=w1(j,i)+(1.0D0-(tmesh(i0+j)-x0)*d)*whelp(j)
 60         continue
          end if
          do 70,j=-nskipt,nskipt
            sum=sum+w1(j,i)
 70       continue
 80     continue
        write (6,*) 'integration error in theta',
     c              sum-(tmesh(nsetp)-tmesh(1))
      else if (icaseb.eq.1) then
        do 110,i=1,nset1
          i1=2*(i-1)
          sum=0.0
          do 100,j=1,nsetp
            w1(j,i)=weights(j)*xmesh(j)**i1
            sum=sum+w1(j,i)
 100      continue
 110    continue
      else if (icaseb.eq.2) then
        call ortho(c)
        do 130,i=1,nset1
          do 120,j=1,nsetp
            w1(j,i)=0.0
 120      continue
 130    continue
        sum=0.0
        do 210,j=1,nsetp
          xj2=xmesh(j)**2
          do 160,i=0,nexp-1
            wx=0.0
            do 150,k=0,nexp-1
              wx=wx+c(k,i)*xj2**k
 150        continue
            w1(j,i+1)=weights(j)*wx
 160      continue
c         w1(j,1)=weights(j)
c         w1(j,2)=weights(j)*(xj2-0.2)
c         w1(j,3)=weights(j)*(xj2**2-2./3.*xj2+1./21.)
 210    continue
      else
      end if
      do 220,i=1,nlm
        setl(llm(i))=.true.
 220  continue
c     do 310,k=1,lmax
c       do 300,l=1,maxa
      do 310,k=1,maxa
        do 300,l=1,lmax
          ila(l,k,1)=0
          ila(l,k,2)=0
          ila(l,k,3)=0
 300    continue
 310  continue
      ila1=0
      do 350,i=1,nlm
        if (itlm(i).ne.0) then
          l=llm(i)
          it=itlm(i)
          it1=it1lm(i)
          if (ila(l,it1,it).eq.0) then
            ila1=ila1+1
            ila(l,it1,it)=ila1
            call mtoac(it,it1,l,acoeff(1,1,ila1))
          end if
        end if
 350  continue
      do 1000,l=0,lmax
        indx(l)=l
 1000 continue 
      e1=-1
      e2=-1
      do 1300,m=1,lmax
        rm2=m*m
        call setplm(lmax,m,nsetp,indx,xmesh,maxt,plm,dplm)
        do 1250,l=m,lmax
          if (setl(l)) then
            nl=0
            do 1050,i=1,nlm
              if (llm(i).eq.l) then
                nl=nl+1
                llist(nl)=i
              end if
 1050       continue
            do 1100,j=1,nsetp
              p=plm(j,l)
              p2=p*p
              dp=dplm(j,l)
              g1l(j)=p2*smesh(j)
              g2l(j)=dp**2*smesh3(j)+p*dp*sc2mesh(j)+rm2*p2*smesh1(j)
 1100       continue

            if (icaseb.eq.0) then
              sum1=0.0D0
              sum2=0.0D0
              do 1130,j=1,nset1
                j1=nskipt*(j-1)+1
                g1j=0.0D0
                g2j=0.0D0
                if (m.eq.0) then
                  do 1110,k=max(1-j1,-nskipt),min(nsetp-j1,nskipt)
                    j1k=j1+k
                    p=plm(j1k,l)
                    dp=dplm(j1k,l)
                    w=w1(k,j)
                    g1j=g1j+w*g1l(j1k)
                    g2j=g2j+w*(dp**2*smesh3(j1k)+p*dp*sc2mesh(j1k))
 1110             continue
                else
                  do 1120,k=max(1-j1,-nskipt),min(nsetp-j1,nskipt)
                    j1k=j1+k
                    w=w1(k,j)
                    g1j=g1j+w*g1l(j1k)
                    g2j=g2j+w*g2l(j1k)
 1120             continue
                end if
                sum1=sum1+2.0D0*g1j
                sum2=sum2+2.0D0*g2j
                g1j=2.0D0*m*g1j
                g2j=2.0D0*m*g2j
                do 1125,i=1,nl
                  i1=llist(i)
                  if (itlm(i1).eq.0) then
                    if (mlm(i1).eq.m) then
                      g1(i1,j)=g1j
                      g2(i1,j)=g2j
                    end if
                  else
                     k=mlm(i1)
                     it=itlm(i1)
                     it1=it1lm(i1)
                     if (m.eq.0) then
                        ac = 0D0
                     else
                        ac = acoeff(k,m,ila(l,it1,it))
                     endif
                     g1(i1,j)=g1(i1,j)+g1j*ac
                     g2(i1,j)=g2(i1,j)+g2j*ac
                  end if
 1125           continue
 1130         continue
            else if ((icaseb.eq.1).or.(icaseb.eq.2)) then
              sum1=0.0D0
              sum2=0.0D0
              do 1160,j=1,nset1
                g1j=0.0
                g2j=0.0
                do 1150,k=1,nsetp
                  w=w1(k,j)
                  g1j=g1j+w*g1l(k)
                  g2j=g2j+w*g2l(k)
 1150           continue
                if (j.eq.1) then
                  sum1=sum1+2.0D0*g1j
                  sum2=sum2+2.0D0*g2j
                end if
                g1j=2.0D0*m*g1j
                g2j=2.0D0*m*g2j
                do 1155,i=1,nl
                  i1=llist(i)
                  if (itlm(i1).eq.0) then
                    if (mlm(i1).eq.m) then
                      g1(i1,j)=g1j
                      g2(i1,j)=g2j
                    end if
                  else
                    k=mlm(i1)
                    it=itlm(i1)
                    it1=it1lm(i1)
                    itmp = ila(l,it1,it)
                    ac=acoeff(k,m,itmp)
                    g1(i1,j)=g1(i1,j)+g1j*ac
                    g2(i1,j)=g2(i1,j)+g2j*ac
                  end if
 1155           continue
 1160         continue
            else
            end if
            if (m.ne.0) then
              if (abs(sum1-1.0D0).gt.e1) then
                e1=abs(sum1-1.0D0)
                le1=l
                me1=m
              end if
              x=l*l+l-1
              if (abs((sum2-x)/x).gt.e2) then
                e2=abs((sum2-x)/x)
                le2=l
                me2=m
              end if
            end if
          end if
 1250   continue
 1300 continue
      write (6,*) 'maximum error in 1. integral:',e1,' at l,m=',le1,me1
      write (6,*) 'maximum error in 2. integral:',e2,' at l,m=',le2,me2

      end


      subroutine msort(nnm,mindex,nna,aindex)

      include 'include1.n.f'
      include 'include2.f'

      integer mindex(maxnlm),nm(maxl),im(maxl),aindex(maxnlm)

      do 50,m=1,maxl
        nm(m)=0
 50   continue
      nnm=0
      nna=0
      do 100,i=1,nnlm
        if (itnlm(i).eq.0) then
          nnm=nnm+1
          nm(mnlm(i))=nm(mnlm(i))+1
        else
          nna=nna+1
        end if
 100  continue
      im(1)=0
      do 200,m=2,maxl
        im(m)=im(m-1)+nm(m-1)
 200  continue
      ima=0
      do 300,i=1,nnlm
        if (itnlm(i).eq.0) then
          m=mnlm(i)
          im(m)=im(m)+1
          mindex(im(m))=i
        else
          ima=ima+1
          aindex(ima)=i
        end if
 300  continue

      end

      subroutine ortho(c)

      implicit double precision (a-h,o-z)

      include 'include1.n.f'

      parameter (m1=maxexp-1)

      double precision c(0:maxexp-1,0:maxexp-1)
      double precision pc(0:2*maxexp-1,0:2*maxexp-1)
      double precision p(0:m1,0:m1),p1(0:m1,0:m1),d(0:m1,0:m1)
      double precision b(0:m1),pb1p1(0:m1,0:m1),work(0:m1)
      integer ipvt(0:m1)

      n1=nexp-1
      do 20,i=0,2*n1
        do 10,j=0,2*n1
          pc(j,i)=0.0
 10     continue
 20   continue
      pc(0,0)=1.
      pc(1,1)=1.
      do 40,i=2,2*nexp-1
        c1=(2.0D0*i-1.0D0)/i
        c2=(i-1.0D0)/i
        do 30,j=0,i
          if (i.gt.0.and.j.gt.0) pc(i,j)=pc(i,j)+c1*pc(i-1,j-1)
          if (i.gt.1) pc(i,j)=pc(i,j)-c2*pc(i-2,j)
 30     continue
 40   continue
      do 60,i=0,n1
        do 50,j=0,n1
          p(i,j)=pc(2*j,2*i)
          p1(i,j)=pc(2*j,2*i)
          d(i,j)=pc(2*j+1,2*i+1)
 50     continue
 60   continue
      call dgeco(p1,maxexp,nexp,ipvt,rcond,work)
      call dgedi(p1,maxexp,nexp,ipvt,det,work,1)
      b(0)=1.
      do 70,i=1,n1
        b(i)=-b(i-1)*(2.0D0*i-1.0D0)*2.0D0*i/(4.0D0*i*i)
 70   continue
      do 130,i=0,n1
        do 120,j=0,n1
          cij=0.0
          do 110,k=0,n1
            cij=cij+p(i,k)*p1(k,j)/b(k)
 110      continue
          pb1p1(i,j)=cij
 120    continue
 130  continue
      do 230,i=0,n1
        do 220,j=0,n1
          cij=0.0
          do 210,k=0,n1
            cij=cij+pb1p1(i,k)*d(k,j)
 210      continue
          c(i,j)=cij
 220    continue
 230  continue
      do 260,j=0,n1
        cjj=c(j,j)
        do 240,i=0,j
          c(i,j)=c(i,j)/cjj
 240    continue
        do 250,i=j+1,n1
          c(i,j)=0.0
 250    continue
 260  continue
      end
