c     Copyright Rasmus Munk Larsen, Stanford University, 2003


      subroutine readspl(cfile,split,sigma)

c     read splittings from file

      implicit double precision (a-h,o-z)
      include 'include1.f'

      real split(maxmodes,maxexp),sigma(maxmodes,maxexp)
      real split1(101),sigma1(101)
      character*80 cfile,splitfile

c     cfile: name of file to read parameters from
c     nexp: number of expansion coefficients to read
c     nmode,lmode: identification of modes:
c     freq: frequency of mode (muHz)
c     split and sigma: observed splittings and errors on those (nHz)
c     nmodes: number of modes read in
c     nexp1: number of expansion coefficients in file
c     ieven=1 if even splittings are provided in file too
c     if ierr=0 all errors are set to err
c     otherwise errors are read from the splitting file
c     splitfile: name of splitting file

      write (6,*) 'start reading splittings, using control file:'
      write (6,*) cfile
      write (6,*) '--------------------'
      call catfile(cfile)
      write (6,*) '--------------------'
      open (10,file=cfile,form='formatted')
      call clnfile(10,11)
      close(10)
      read(11,*) nexp1,ieven,ierr,err
      read (11,'(A80)') splitfile
      close(11)
      write (6,*) 'splittings read from file: '
      write (6,'(A80)') splitfile
      write (6,*) 
      write (6,*) 
     c  '  i   l   n   freq  (split(j),j=1,nexp) (sigma(j),j=1,nexp)'
      write (6,*)
      nmodes=0
      open (10,file=splitfile,form='formatted')
 10   continue
        if (ierr.eq.0) then
          read (10,*,end=40) rl,rn,f,(split1(i),i=1,nexp1)
          do 20,i=1,nexp1
            sigma1(i)=err
 20       continue
        else
          read (10,*,end=40) rl,rn,f,(split1(i),i=1,nexp1),
     c                     (sigma1(i),i=1,nexp1)
        end if
        nmodes=nmodes+1
        if (nmodes.gt.maxmodes) then
          write (6,*) 'too many modes, only dimensioned for maxmodes=',
     c                maxmodes,' modes'
          close(10)
          stop
        end if
        lmode(nmodes)=nint(rl)
        nmode(nmodes)=nint(rn)
        freq(nmodes)=f
        do 30,i=1,nexp
          split(nmodes,i)=split1(1+(ieven+1)*(i-1))
          sigma(nmodes,i)=sigma1(1+(ieven+1)*(i-1))
 30     continue
        write (6,'(3i4,f8.1,100f7.2)') nmodes,nint(rl),nint(rn),f,
     c    (split(nmodes,i),i=1,nexp),(sigma(nmodes,i),i=1,nexp)
      goto 10
 40   continue
      
      end

      subroutine setspl(nker,ker,omega,split)

c     sets integrals of omega*ker, doesn't do expansions
c     both the kernels and omega must be given on radmesh(i),i=1,nrad

      implicit double precision (a-h,o-z)
      include 'include1.f'

      real omega(maxrad)
      real ker(maxrad,maxmodes)
      double precision weights(maxrad)
      real split(maxmodes)

      do 5,i=1,nrad
        dradmesh(i)=radmesh(i)
 5    continue
      call dsetint(1,nrad,dradmesh,weights)
      do 10,j=1,nrad
        weights(j)=weights(j)*omega(j)
 10   continue
      do 20,i=1,nker
        sum=0.0D0
        do 30,j=1,nrad
          sum=sum+ker(j,i)*weights(j)
 30     continue
        split(i)=sum
 20   continue

      end

      subroutine dsetspl(nker,ker,omega,split)

c     sets integrals of omega*ker, doesn't do expansions
c     both the kernels and omega must be given on radmesh(i),i=1,nrad

      implicit double precision (a-h,o-z)
      include 'include1.f'

      double precision omega(maxrad)
      real ker(maxrad,maxmodes)
      double precision split(maxmodes),weights(maxrad)

      call dsetint(1,nrad,dradmesh,weights)
      do 10,j=1,nrad
        weights(j)=weights(j)*omega(j)
 10   continue
      do 20,i=1,nker
        sum=0.0D0
        do 30,j=1,nrad
          sum=sum+ker(j,i)*weights(j)
 30     continue
        split(i)=sum
 20   continue

      end

      subroutine mtoa(icase,na,l,split,a)

c     calculates a coeffecients a(i),i=1,na given splittings
c     split(m),m=1,l 
c     expands in m/l for icase=1, in m/L for icase=2
c     and in truly orthogonal functions for icase=3

      implicit double precision (a-h,o-z)
      include 'include1.f'
      include 'include2.f'
      parameter (maxb=2*maxa-1)

      real split(maxl),a(maxa)
      double precision mlmesh(maxl),p(maxl,0:maxb),dp(maxl,0:maxb)
      integer indx(0:maxl)
      double precision b(maxa,maxa),right(maxa),help(maxa)

      na1=min(l,na)
      na2=2*na1-1
      l2=l*(l+1)
      capl=sqrt(l2+0.0D0)
      xl=l
      if (icase.eq.2) xl=capl
      x=1.0D0/xl
      do 100,m=1,l
        mlmesh(m)=m*x
 100  continue
      do 120,i=0,na2
        indx(i)=i
 120  continue
      call setplm(na2,0,l,indx,mlmesh,maxl,p,dp)
c     the normal normalization of the p's is restored
      do 150,i=1,na2,2
        do 140,m=1,l
          p(m,i)=p(m,i)/sqrt(i+0.5D0)
 140    continue
 150  continue
      if (icase.eq.3) then
c     orthogonalize the polynomials
        do 190,i=1,na2,2
          do 170,j=1,i-2,2
            c=ddot1(l,p(1,i),p(1,j))/ddot1(l,p(1,j),p(1,j))
            do 160,m=1,l
              p(m,i)=p(m,i)-c*p(m,j)
 160        continue
 170      continue
          c=p(l,i)
          do 180,m=1,l
            p(m,i)=p(m,i)/c
 180      continue
 190    continue
      end if
      do 250,i=1,na1
        i1=2*i-1
        ri=0.0D0
        do 210,k=1,l
          ri=ri+p(k,i1)*split(k)
 210    continue
        right(i)=ri
        do 240,j=1,na1
          j1=2*j-1
          bji=0.0D0
          do 230,k=1,l
            bji=bji+p(k,i1)*p(k,j1)
 230      continue
          b(j,i)=bji
 240    continue
 250  continue
      call dpoco(b,maxa,na1,rcond,help,info)
      call dposl(b,maxa,na1,right)
      do 300,i=1,na1
        a(i)=right(i)*x
 300  continue

      end

      subroutine mtoac(icase,na,l,coeff)

c     calculates a coeffecients coeff(i,m),i=1,na and m=1,l such
c     that a(i)=sum_{m=1}^l coeff(i,m) split(m)
c     expands in m/l for icase=1, in m/L for icase=2
c     and in truly orthogonal functions for icase=3

      implicit double precision (a-h,o-z)
      include 'include1.f'
      include 'include2.f'
      parameter (maxb=2*maxa-1)

      double precision coeff(maxa,maxl)
      double precision mlmesh(maxl),p(maxl,0:maxb),dp(maxl,0:maxb)
      integer indx(0:maxl)
      double precision b(maxa,maxa),right(maxa),help(maxa)

      na1=min(l,na)
      na2=2*na1-1
      l2=l*(l+1)
      capl=sqrt(l2+0.0D0)
      xl=l
      if (icase.eq.2) xl=capl
      x=1.0D0/xl
      do 100,m=1,l
        mlmesh(m)=m*x
 100  continue
      do 120,i=0,na2
        indx(i)=i
 120  continue
      call setplm(na2,0,l,indx,mlmesh,maxl,p,dp)
c     the normal normalization of the p's is restored
      do 150,i=1,na2,2
        do 140,m=1,l
          p(m,i)=p(m,i)/sqrt(i+0.5D0)
 140    continue
 150  continue
      if (icase.eq.3) then
c     orthogonalize the polynomials
        do 190,i=1,na2,2
          do 170,j=1,i-2,2
            c=ddot1(l,p(1,i),p(1,j))/ddot1(l,p(1,j),p(1,j))
            do 160,m=1,l
              p(m,i)=p(m,i)-c*p(m,j)
 160        continue
 170      continue
          c=p(l,i)
          do 180,m=1,l
            p(m,i)=p(m,i)/c
 180      continue
 190    continue
      end if
      do 250,i=1,na1
        i1=2*i-1
        do 240,j=1,na1
          j1=2*j-1
          bji=0.0D0
          do 230,k=1,l
            bji=bji+p(k,i1)*p(k,j1)
 230      continue
          b(j,i)=bji
 240    continue
 250  continue
      call dpoco(b,maxa,na1,rcond,help,info)
c     write (6,*) icase,l,rcond
      do 310,m=1,l
        do 270,j=1,na1
          right(j)=p(m,2*j-1)
 270    continue
        call dposl(b,maxa,na1,right)
        do 300,i=1,na1
          coeff(i,m)=right(i)*x
 300    continue
 310  continue

      end

      subroutine atom(icase,na,l,a,split)

c     calculates splittings split(m),m=1,l given a coeffecients
c     a(2i-1),i=1,na
c     expands in m/l for icase=1, in m/L for icase=2
c     and in truly orthogonal functions for icase=3

      implicit double precision (a-h,o-z)
      include 'include1.f'
      include 'include2.f'
      parameter (maxb=2*maxa-1)

      real split(maxl),a(maxa)
      double precision mlmesh(maxl),p(maxl,0:maxb),dp(maxl,0:maxb)
      integer indx(0:maxl)

      na1=min(l,na)
      na2=2*na1-1
      l2=l*(l+1)
      capl=sqrt(l2+0.0D0)
      xl=l
      if (icase.eq.2) xl=capl
      x=1.0D0/xl
      do 100,m=1,l
        mlmesh(m)=m*x
 100  continue
      do 120,i=0,na2
        indx(i)=i
 120  continue
      call setplm(na2,0,l,indx,mlmesh,maxl,p,dp)
c     the normal normalization of the p's is restored
      do 150,i=1,na2,2
        do 140,m=1,l
          p(m,i)=p(m,i)/sqrt(i+0.5D0)
 140    continue
 150  continue
      if (icase.eq.3) then
c     orthogonalize the polynomials
        do 190,i=1,na2,2
          do 170,j=1,i-2,2
            c=ddot1(l,p(1,i),p(1,j))/ddot1(l,p(1,j),p(1,j))
            do 160,m=1,l
              p(m,i)=p(m,i)-c*p(m,j)
 160        continue
 170      continue
          c=p(l,i)
          do 180,m=1,l
            p(m,i)=p(m,i)/c
 180      continue
 190    continue
      end if
      do 250,m=1,l
        sm=0.0D0
        do 230,i=1,na1
          i1=2*i-1
          sm=sm+p(m,i1)*a(i)
 230    continue
        split(m)=sm*xl
 250  continue

      end

      subroutine atomc(icase,na,l,coeff)

c     calculates a coeffecients coeff(i,m),i=1,na and m=1,l such
c     split(m)=sum_{i=1}^na coeff(i,m) a(i)
c     expands in m/l for icase=1, in m/L for icase=2
c     and in truly orthogonal functions for icase=3

      implicit double precision (a-h,o-z)
      include 'include1.f'
      include 'include2.f'
      parameter (maxb=2*maxa-1)

      double precision coeff(maxa,maxl)
      double precision mlmesh(maxl),p(maxl,0:maxb),dp(maxl,0:maxb)
      integer indx(0:maxl)

      na1=min(l,na)
      na2=2*na1-1
      l2=l*(l+1)
      capl=sqrt(l2+0.0D0)
      xl=l
      if (icase.eq.2) xl=capl
      x=1.0D0/xl
      do 100,m=1,l
        mlmesh(m)=m*x
 100  continue
      do 120,i=0,na2
        indx(i)=i
 120  continue
      call setplm(na2,0,l,indx,mlmesh,maxl,p,dp)
c     the normal normalization of the p's is restored
      do 150,i=1,na2,2
        do 140,m=1,l
          p(m,i)=p(m,i)/sqrt(i+0.5D0)
 140    continue
 150  continue
      if (icase.eq.3) then
c     orthogonalize the polynomials
        do 190,i=1,na2,2
          do 170,j=1,i-2,2
            c=ddot1(l,p(1,i),p(1,j))/ddot1(l,p(1,j),p(1,j))
            do 160,m=1,l
              p(m,i)=p(m,i)-c*p(m,j)
 160        continue
 170      continue
          c=p(l,i)
          do 180,m=1,l
            p(m,i)=p(m,i)/c
 180      continue
 190    continue
      end if
      do 250,m=1,l
        do 230,i=1,na2,2
          coeff(i,m)=p(m,i)*xl
 230    continue
 250  continue

      end
