c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine catfile(filename)

      character*80 filename
      character*80 line

      open(14,file=filename,form='formatted')
 10   read(14,'(A80)',end=20) line
      write(*,'(A80)') line
      goto 10
 20   continue
      close(14)
      end

      subroutine nextline(unit,line)

c     reads lines from unit until it finds a line not starting
c     with a # and returns the line

      integer unit
      character*(*) line

 10   read (unit,'(A80)') line
      if (line(1:1).eq.'#') goto 10

      end

      subroutine nxtline1(unit,line)

c     reads lines from unit until it finds a line not starting
c     with a # and returns the line

      integer unit
      character*(*) line

 10   read (unit,'(A1000)') line
      if (line(1:1).eq.'#') goto 10

      end

      subroutine clnfile(iunit1,iunit2)

c     takes non-comment lines from iunit1 and places them in
c     a scratch file iunit2

      character*80 line

      open(iunit2,status='scratch',form='formatted')
 10   read(iunit1,'(a80)',end=20) line
      if (line(1:1).eq.'#') goto 10
      write(iunit2,'(a80)') line
      goto 10
 20   continue
      rewind(iunit2)
      end

      subroutine addsuff(s1,suf,s2)

      character*80 s1,suf,s2

      i=index(s1,' ')
      i1=index(suf,' ')
      s2=s1
      s2(i:i+i1-2)=suf(1:i1-1)

      end

      subroutine setint(method,nmesh,mesh,intw)

c     sets integration weights for mesh into intw
c     uses trapezoidal rule for method=1
c     2'nd order, exact up to and including 1'st degree polynomials
c     for method=2 does richardson extrapolation from (nmesh+1)/2 and
c     nmesh meshpoints for nmesh odd and an approximation to this for
c     nmesh even, equivalent to simpsons rule for nmesh odd and evenly
c     spaced mesh
c     4'th order for odd # of meshpoints, 3'rd for even # of meshpoints
c     generally exact for up to and including 1'st degree polynomials
c     but exact for up to and including 3'rd degree polynomials for
c     odd # of meshpoints and evenly spaced mesh (simpson case)
c     for method=3 does a fit to a 3'rd degree polynomial for each
c     subinterval from the nearest 4 points and integrates this
c     polynomial over the subinterval
c     always 4'th order and exact up to 3'rd degree polynomials
c     less accurate than method 2 for high degree polynomials

      real mesh(nmesh),intw(nmesh)


      if (method.eq.1) then
        do 10,i=1,nmesh
          i1=max(1,i-1)
          i2=min(nmesh,i+1)
          intw(i)=0.5E0*(mesh(i2)-mesh(i1))
 10     continue
      else if (method.eq.2) then
        c=2.0E0/3.0E0
        do 20,i=1,nmesh
          i1=max(1,i-1)
          i2=min(nmesh,i+1)
          intw(i)=c*(mesh(i2)-mesh(i1))
 20     continue
        if (mod(nmesh,2).eq.1) then
          c=-1.0E0/6.0E0
          do 25,i=1,nmesh,2
            i1=max(1,i-2)
            i2=min(nmesh,i+2)
            intw(i)=intw(i)+c*(mesh(i2)-mesh(i1))
 25       continue
        else
          c=-1.0E0/12.0E0
          do 30,i=1,nmesh+1,2
            i0=min(nmesh,max(1,i))
            i1=max(1,i-2)
            i2=min(nmesh,i+2)
            intw(i0)=intw(i0)+c*(mesh(i2)-mesh(i1))
 30       continue
          c=-1.0E0/12.0E0
          do 40,i=0,nmesh,2
            i0=min(nmesh,max(1,i))
            i1=max(1,i-2)
            i2=min(nmesh,i+2)
            intw(i0)=intw(i0)+c*(mesh(i2)-mesh(i1))
 40       continue
        end if
      else if (method.eq.3) then
        do 100,i=1,nmesh
          intw(i)=0.0E0
 100    continue
        do 120,i=1,nmesh-1
          i1=i-1
          i2=i
          i3=i+1
          i4=i+2
          if (i.eq.1) i1=i+3
          if (i.eq.nmesh-1) i4=i-2
          h1=mesh(i2)-mesh(i1)
          h2=mesh(i3)-mesh(i2)
          h3=mesh(i4)-mesh(i3)
          h21=mesh(i3)-mesh(i1)
          h32=mesh(i4)-mesh(i2)
          h321=mesh(i4)-mesh(i1)
          intw(i1)=intw(i1)-h2**3*(h2+2*h3)/(12*h1*h21*h321)
          intw(i2)=intw(i2)+h2*(4*h1*h2+h2*h2+6*h1*h3+2*h2*h3)/
     c                      (12*h1*h32)
          intw(i3)=intw(i3)+h2*(4*h3*h2+h2*h2+6*h1*h3+2*h2*h1)/
     c                      (12*h3*h21)
          intw(i4)=intw(i4)-h2**3*(h2+2*h1)/(12*h3*h32*h321)
 120    continue
      end if
      end

      subroutine interp(method,nin,xin,yin,nout,xout,yout)

c     interpolates function given by (xin(i),yin(i),i=1,nin)
c     to a new mesh given by (xout(i),i=1,nout) into (yout(i),i=1,nout)
c     both xin and xout must be increasing, and nin.ge.2
c     uses linear interpolation for method=1
c     extrapolates linearly outside interval given by xin

      real xin(nin),yin(nin),xout(nout),yout(nout)

      n=method
      n=max(n,1)
      n=min(n,nin-1)
      np1=n+1
      nd2=n/2
      if (n.eq.1) then
        n2=2
        do 10,i=1,nout
 20       if ((xout(i).gt.xin(n2)).and.(n2.lt.nin)) then
            n2=n2+1
            goto 20
          end if
          yout(i)=yin(n2-1)+(yin(n2)-yin(n2-1))*
     c            (xout(i)-xin(n2-1))/(xin(n2)-xin(n2-1))
 10     continue
      else
        if (2*(n/2).eq.n) then
c n is even so find nearest point in input to the target output point
          n0=1
          xi0=xin(n0)
          xi1=xin(n0+1)
          do 110,i=1,nout
            xo=xout(i)
 105        if ((n0.lt.nin).and.(abs(xo-xi0).gt.abs(xo-xi1))) then
              n0=n0+1
              xi0=xi1
              xi1=xin(n0+1)
              goto 105
            end if
            n1=min(max(1,n0-nd2),nin-n)
c           write (6,*) n0,n1,xi0,xi1,xo,(xin(j),j=n1,n1+np1-1)
            call rpolint(xin(n1),yin(n1),np1,xo,yout(i),dummy)
 110      continue
        else
c n is odd
          n2=2
          do 210,i=1,nout
            xo=xout(i)
 220        if ((xo.gt.xin(n2)).and.(n2.lt.nin)) then
              n2=n2+1
              goto 220
            end if
            n1=min(max(1,n2-1-nd2),nin-n)
            call rpolint(xin(n1),yin(n1),np1,xo,yout(i),dummy)
 210      continue
        end if
      end if

      end

      subroutine dsetint(method,nmesh,mesh,intw)

c     double precision version
c     sets integration weights for mesh into intw
c     uses trapezoidal rule for method=1
c     2'nd order, exact up to and including 1'st degree polynomials
c     for method=2 does richardson extrapolation from (nmesh+1)/2 and
c     nmesh meshpoints for nmesh odd and an approximation to this for
c     nmesh even, equivalent to simpsons rule for nmesh odd and evenly
c     spaced mesh
c     4'th order for odd # of meshpoints, 3'rd for even # of meshpoints
c     generally exact for up to and including 1'st degree polynomials
c     but exact for up to and including 3'rd degree polynomials for
c     odd # of meshpoints and evenly spaced mesh (simpson case)
c     for method=3 does a fit to a 3'rd degree polynomial for each
c     subinterval from the nearest 4 points and integrates this
c     polynomial over the subinterval
c     always 4'th order and exact up to 3'rd degree polynomials
c     less accurate than method 2 for high degree polynomials

      implicit double precision (a-h,o-z)

      double precision mesh(nmesh),intw(nmesh)

      if (method.eq.1) then
        do 10,i=1,nmesh
          i1=max(1,i-1)
          i2=min(nmesh,i+1)
          intw(i)=0.5D0*(mesh(i2)-mesh(i1))
 10     continue
      else if (method.eq.2) then
        c=2.0D0/3.0D0
        do 20,i=1,nmesh
          i1=max(1,i-1)
          i2=min(nmesh,i+1)
          intw(i)=c*(mesh(i2)-mesh(i1))
 20     continue
        if (mod(nmesh,2).eq.1) then
          c=-1.0D0/6.0D0
          do 25,i=1,nmesh,2
            i1=max(1,i-2)
            i2=min(nmesh,i+2)
            intw(i)=intw(i)+c*(mesh(i2)-mesh(i1))
 25       continue
        else
          c=-1.0D0/12.0D0
          do 30,i=1,nmesh+1,2
            i0=min(nmesh,max(1,i))
            i1=max(1,i-2)
            i2=min(nmesh,i+2)
            intw(i0)=intw(i0)+c*(mesh(i2)-mesh(i1))
 30       continue
          c=-1.0D0/12.0D0
          do 40,i=0,nmesh,2
            i0=min(nmesh,max(1,i))
            i1=max(1,i-2)
            i2=min(nmesh,i+2)
            intw(i0)=intw(i0)+c*(mesh(i2)-mesh(i1))
 40       continue
        end if
      else if (method.eq.3) then
        do 100,i=1,nmesh
          intw(i)=0.0D0
 100    continue
        do 120,i=1,nmesh-1
          i1=i-1
          i2=i
          i3=i+1
          i4=i+2
          if (i.eq.1) i1=i+3
          if (i.eq.nmesh-1) i4=i-2
          h1=mesh(i2)-mesh(i1)
          h2=mesh(i3)-mesh(i2)
          h3=mesh(i4)-mesh(i3)
          h21=mesh(i3)-mesh(i1)
          h32=mesh(i4)-mesh(i2)
          h321=mesh(i4)-mesh(i1)
          intw(i1)=intw(i1)-h2**3*(h2+2*h3)/(12*h1*h21*h321)
          intw(i2)=intw(i2)+h2*(4*h1*h2+h2*h2+6*h1*h3+2*h2*h3)/
     c                      (12*h1*h32)
          intw(i3)=intw(i3)+h2*(4*h3*h2+h2*h2+6*h1*h3+2*h2*h1)/
     c                      (12*h3*h21)
          intw(i4)=intw(i4)-h2**3*(h2+2*h1)/(12*h3*h32*h321)
 120    continue
      end if

      end

      subroutine dinterp(method,nin,xin,yin,nout,xout,yout)

c     interpolates function given by (xin(i),yin(i),i=1,nin)
c     to a new mesh given by (xout(i),i=1,nout) into (yout(i),i=1,nout)
c     both xin and xout must be increasing, and nin.ge.2
c     uses linear interpolation for method=1
c     extrapolates linearly outside interval given by xin

      implicit double precision (a-h,o-z)

      double precision xin(nin),yin(nin),xout(nout),yout(nout)

c     write (6,*) 1
      n=method
      n=max(n,1)
      n=min(n,nin-1)
      np1=n+1
      nd2=n/2
      if (n.eq.1) then
        n2=2
        do 10,i=1,nout
c       write (6,*) n,i,nout,n2,xout(i),xin(n2),nin,xin(n2-1)
 20       if ((xout(i).gt.xin(n2)).and.(n2.lt.nin)) then
            n2=n2+1
            goto 20
          end if
          yout(i)=yin(n2-1)+(yin(n2)-yin(n2-1))*
     c            (xout(i)-xin(n2-1))/(xin(n2)-xin(n2-1))
 10     continue
      else
        if (2*(n/2).eq.n) then
c n is even so find nearest point in input to the target output point
          n0=1
          xi0=xin(n0)
          xi1=xin(n0+1)
          do 110,i=1,nout
            xo=xout(i)
 105        if ((n0.lt.nin).and.(abs(xo-xi0).gt.abs(xo-xi1))) then
              n0=n0+1
              xi0=xi1
              xi1=xin(n0+1)
              goto 105
            end if
            n1=min(max(1,n0-nd2),nin-n)
c           write (6,*) n0,n1,xi0,xi1,xo,(xin(j),j=n1,n1+np1-1)
            call dpolint(xin(n1),yin(n1),np1,xo,yout(i),dummy)
 110      continue
        else
c n is odd
          n2=2
          do 210,i=1,nout
            xo=xout(i)
 220        if ((xo.gt.xin(n2)).and.(n2.lt.nin)) then
              n2=n2+1
              goto 220
            end if
            n1=min(max(1,n2-1-nd2),nin-n)
            call dpolint(xin(n1),yin(n1),np1,xo,yout(i),dummy)
 210      continue
        end if
      end if

      end

      double precision function ddot1(n,a,b)

      double precision a(n),b(n),dot

      dot=0.0D0
      do 100,i=1,n
        dot=dot+a(i)*b(i)
 100  continue
      ddot1=dot
      end

      subroutine polint(xa,ya,n,x,y,dy)
      parameter (nmax=10)
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end

      subroutine dpolint(xa,ya,n,x,y,dy)
      parameter (nmax=10)
      implicit double precision (a-h,o-z)
      dimension xa(n),ya(n),c(nmax),d(nmax)
c     do i=1,n
c       write (6,*) i,xa(i),ya(i)
c     end do
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end


      subroutine rpolint(xa,ya,n,x,y,dy)
      parameter (nmax=10)
      implicit real (a-h,o-z)
      dimension xa(n),ya(n),c(nmax),d(nmax)
c     do i=1,n
c       write (6,*) i,xa(i),ya(i)
c     end do
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end


      subroutine dsetip(method,nin,xin,nout,xout,n,ifirst,maxc,coeff)

c     sets coeffecients for interpolation from at net given by
c     (xin(i),i=1,nin) to a new mesh given by (xout(i),i=1,nout)
c     both xin and xout must be increasing, and nin.ge.2
c     uses n'th order interpolation for method=n and extrapolates
c     outside interval given by xin
c     n and coeff is set such that
c     yout(i)=sum_{j=0}^{n} y(ifirst(i)+j)*coeff(j,i)

      implicit double precision (a-h,o-z)

      parameter (maxmethod=10)

      double precision xin(nin),xout(nout),coeff(0:maxc,nout)
      double precision yin(0:maxmethod)
      integer ifirst(nout)

      n=method
      n=max(n,1)
      n=min(n,nin-1)
      np1=n+1
      nd2=n/2
      if (2*(n/2).eq.n) then
c n is even so find nearest point in input to the target output point
        n0=1
        xi0=xin(n0)
        xi1=xin(n0+1)
        do 150,i=1,nout
          xo=xout(i)
 110      if ((n0.lt.nin).and.(abs(xo-xi0).gt.abs(xo-xi1))) then
            n0=n0+1
            xi0=xi1
            xi1=xin(n0+1)
            goto 110
          end if
          n1=min(max(1,n0-nd2),nin-n)
          ifirst(i)=n1
          do 130,j=0,n
            do 120,k=0,n
              yin(k)=0.0D0
 120        continue
            yin(j)=1.0D0
            call dpolint(xin(n1),yin,np1,xo,yo,dummy)
            coeff(j,i)=yo
 130      continue
 150    continue
      else
c n is odd
        n2=2
        do 250,i=1,nout
          xo=xout(i)
 210      if ((xo.gt.xin(n2)).and.(n2.lt.nin)) then
            n2=n2+1
            goto 210
          end if
          n1=min(max(1,n2-1-nd2),nin-n)
          ifirst(i)=n1
          do 230,j=0,n
            do 220,k=0,n
              yin(k)=0.0D0
 220        continue
            yin(j)=1.0D0
            call dpolint(xin(n1),yin,np1,xo,yo,dummy)
            coeff(j,i)=yo
 230      continue
 250    continue
      end if

      end

      subroutine setip(method,nin,xin,nout,xout,n,ifirst,maxc,coeff)

c     sets coeffecients for interpolation from at net given by
c     (xin(i),i=1,nin) to a new mesh given by (xout(i),i=1,nout)
c     both xin and xout must be increasing, and nin.ge.2
c     uses n'th order interpolation for method=n and extrapolates
c     outside interval given by xin
c     n and coeff is set such that
c     yout(i)=sum_{j=0}^{n} y(ifirst(i)+j)*coeff(j,i)

      parameter (maxmethod=10)
     
      real xin(nin),xout(nout),coeff(0:maxc,nout)
      real yin(0:maxmethod)
      integer ifirst(nout)

      n=method
      n=max(n,1)
      n=min(n,nin-1)
      np1=n+1
      nd2=n/2
      if (2*(n/2).eq.n) then
c n is even so find nearest point in input to the target output point
        n0=1
        xi0=xin(n0)
        xi1=xin(n0+1)
        do 150,i=1,nout
          xo=xout(i)
 110      if ((n0.lt.nin).and.(abs(xo-xi0).gt.abs(xo-xi1))) then
            n0=n0+1
            xi0=xi1
            xi1=xin(n0+1)
            goto 110
          end if
          n1=min(max(1,n0-nd2),nin-n)
          ifirst(i)=n1
          do 130,j=0,n
            do 120,k=0,n
              yin(k)=0.0
 120        continue
            yin(j)=1.0
            call polint(xin(n1),yin,np1,xo,yo,dummy)
            coeff(j,i)=yo
 130      continue
 150    continue
      else
c n is odd
        n2=2
        do 250,i=1,nout
          xo=xout(i)
 210      if ((xo.gt.xin(n2)).and.(n2.lt.nin)) then
            n2=n2+1
            goto 210
          end if
          n1=min(max(1,n2-1-nd2),nin-n)
          ifirst(i)=n1
          do 230,j=0,n
            do 220,k=0,n
              yin(k)=0.0
 220        continue
            yin(j)=1.0
            call polint(xin(n1),yin,np1,xo,yo,dummy)
            coeff(j,i)=yo
 230      continue
 250    continue
      end if

      end


      subroutine   lir(z,zi,y,yi,ii,id,nt,l,inter)
c     subroutine  lir1(z,zi,y,yi,ii,id,nt,l,inter)
c
c
c
c                interpolation/extrapolation routine
c
c
c     for a such that z=zi(a),  sets y(i)=yi(i,a), i=1,ii
c
c     zi(n),yi(i,n) must be supplied for n=1,nt and i=1,ii
c     id is first dimension of yi
c
c     inter is set to 1 for interpolation and 0 for extrapolation
c
c     if l.le.1, scan to find the zi(n) which immediately bound z
c                starts at n=1
c     if l.gt.1, scan starts from value of n from previous call of lir
c
c
c     lir use cubic interpolation/extrapolation unless nt.lt.4
c     lir1 use linear interpolation/extrapolation
c
c     note
c     ****
c     most of the computation is performed in single precision
c
c
c
      dimension zi(1),y(1),yi(1),a(4)
      data n/-1/
c
      il=0
      go to 1
      entry lir1(z,zi,y,yi,ii,id,nt,l,inter)
      il=1
    1 continue
      ir=1
c
c     check nt and reset il if necessary
      if(nt.lt.2) go to 101
      if(nt.lt.4) il=1
c
c     addressing constants
      inter=1
      ir1=ir-1
      ird=ir*id
      iir=(ii-1)*ir+1
      j=(nt-1)*ir+1
      diff=zi(j)-zi(1)
c
c     set index for start of search
      n=(n-2)*ir+1
      if(l.le.1.or.n.lt.1) n=1
c
c     determine position of z within zi
    2 if(n.gt.j) go to 8
      if(diff) 4,102,3
    3 if(zi(n)-z) 5,6,9
    4 if(zi(n)-z) 9,6,5
    5 n=n+ir
      go to 2
c
c     set y when z lies on a mesh point
    6 j=(n-1)*id
      do 7 i=1,iir
      y(i)=yi(i+j)
    7 if(y(i).eq.0.) y(i+ir1)=0.
      go to 30
c
c     control when z does not lie on a mesh point
    8 inter=0
    9 if(n.le.1) inter=0
      if(il.eq.1) go to 20
c
c     cubic interpolation/extrapolation
c
c     pivotal point (m) and point (k) closest to z
   10 m=n
      k=3
      if(n.gt.1+ir) go to 11
      m=1+ir+ir
      k=n
   11 if(n.lt.j) go to 12
      m=j-ir
      k=4
c
c
c     weighting factors
   12 y1=zi(m-ir*2)
      y2=zi(m-ir)
      y3=zi(m)
      y4=zi(m+ir)
c
      z1=z-y1
      z2=z-y2
      z3=z-y3
      z4=z-y4
c
   13 z12=z1*z2
      z34=z3*z4
c
   14 a(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
      a(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
      a(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
      a(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))
c
c     correct a(k)
   15 diff=a(1)+a(2)+a(3)+a(4)
      a(k)=(1.d0+a(k))-diff
c
c     compute y
   16 m=(m-1)/ir-3
      m=m*ird
      do 18 i=1,iir
      k=i+m
      yy=0.
      do 17 j=1,4
      k=k+ird
      diff=yi(k)
   17 yy=yy+a(j)*diff
      y(i)=yy
   18 if(y(i).eq.0.) y(i+ir1)=0.
      go to 30
c
c     linear interpolation/extrapolation
   20 if(n.eq.1) n=1+ir
      if(n.gt.j) n=j
      z1=zi(n)
      y1=(z1-z)/(z1-zi(n-ir))
      y2=1.0-y1
      j=(n-1)*id
      m=j-ird
      do 21 i=1,iir,ir
      y(i)=y1*yi(i+m)+y2*yi(i+j)
   21 if(y(i).eq.0.) y(i+ir1)=0.
c
c     reset n
   30 n=(n+ir-1)/ir
      return
c
c
c     diagnostics
  101 write(6,1001) nt
      return
  102 write(6,1002) zi(1),nt,zi(j)
      return
c
 1001 format(/1x,10('*'),5x,'there are fewer than two data points in',
     *      ' lir     nt =',i4,5x,10('*')/)
 1002 format(/1x,10('*'),5x,'extreme values of independent variable',
     *      ' equal in lir',5x,10('*')/16x,'zi(1) =',1pe13.5,',   ',
     *       'zi(',i4,') =',1pe13.5/)
c
      end
