c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      subroutine setplm(l,m,nx,indx,x,nplm,plm,dplm)

c     sets plm(i,indx(l))=plm(x(i),l)
c     sets dplm(i,indx(l))=dplm(x(i),l)/dx for l=m,l and i=1,nx
c     the plm's are normalized such that int from -1 to 1 plm sup 2 dx=1

      implicit double precision (a-h,o-z)

      parameter (nmax=10001,eps=1.0D-12,lbig=1000)

c     lbig is used not to get indexing errors when debugging

      double precision x(nx),plm(nplm,0:lbig),dplm(nplm,0:lbig)
      double precision x1(nmax),s(nmax),x2(nmax),x3(nmax)
      integer indx(0:lbig)

      m2=m**2
      im=indx(m)
      do 10,i=1,nx
        x1(i)=min(1.0D0-eps,max(-1.0D0+eps,x(i)))
        s(i)=sqrt(1.0D0-x1(i)*x1(i))
        plm(i,im)=sqrt((2*m+1)/2.0D0)
        x2(i)=1.0D0/(x1(i)*x1(i)-1.0D0)
        x3(i)=x1(i)*x2(i)
 10   continue
      fact=-1.0D0
      do 30,i=1,m
        fact=fact+2.0D0
        c=fact/sqrt(2.0D0*i*(2*i-1))
        do 20,j=1,nx
          plm(j,im)=-plm(j,im)*c*s(j)
 20     continue
 30   continue
      if (l.gt.m) then
        im1=indx(m+1)
c       c=2m+1
        c=sqrt(2*m+3.0D0)
        do 50,i=1,nx
          plm(i,im1)=x1(i)*c*plm(i,im)
 50     continue
      end if
      do 110,l0=m+2,l
        l02=l0*l0
        il0=indx(l0)
        il1=indx(l0-1)
        il2=indx(l0-2)
        c1=sqrt((4*l02-1.0D0)/(l02-m2))
        c2=sqrt(((2*l0+1.0D0)*(l0+m-1)*(l0-m-1))/((2*l0-3)*(l02-m2)))
        do 100,i=1,nx
          plm(i,il0)=c1*x1(i)*plm(i,il1)-c2*plm(i,il2)
 100    continue
 110  continue
      do 200,i=1,nx
        dplm(i,im)=m*x3(i)*plm(i,im)
 200  continue
      do 220,l0=m+1,l
        rl0=l0
        il0=indx(l0)
        il1=indx(l0-1)
        c=sqrt((2*l0+1.0D0)*(l0*l0-m2)/(2*l0-1))
        do 210,i=1,nx
          dplm(i,il0)=rl0*x3(i)*plm(i,il0)-c*x2(i)*plm(i,il1)
 210    continue
 220  continue
      end

