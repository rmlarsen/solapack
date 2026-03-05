c     Copyright Rasmus Munk Larsen, Stanford University, 2003

c     
c**********************************************************************
c     
      subroutine set_COG(avker, ldav, cogparms)
      implicit none
      include '2dsola.n.include.f'
      integer ldav
      double precision avker(ldav,*), cogparms(maxtargets,*)
      double precision rsum,tsum,rsum_old,tsum_old,pi
      integer i,j,k,lowidx,upidx,idx
      double precision wint(maxpts)     
      double precision lininterp,integral
      integer idebug
      parameter(idebug=0)
      external integral,lininterp

      pi    = 2.0*acos(0d0)
c     
c     Set quadrature weights
      call set_qweights(1,wint)

      do i=1,N_targets                  
c     
c     Find center of gravity of averaging kernels.
c     
         rsum = 0d0
         tsum = 0d0
         do j=1,N_points
            rsum = rsum + avker(j,i) * rthetw(j,1) * wint(j)
            tsum = tsum + avker(j,i) * rthetw(j,2) * wint(j)
         enddo
         cogparms(i,1) = rsum
         cogparms(i,2) = tsum
c     
c     Determine lower and upper quartile points
c
         rsum = 0d0
         lowidx = 0
         upidx = 0
         do j=1,N_rad
            rsum_old = rsum
            do k=1,N_theta
               idx = (j-1)*N_theta+k
               rsum = rsum + avker(idx,i) * wint(idx)
            enddo
            if (rsum.gt.0.25 .and. lowidx.eq.0)  then
               lowidx = (j-1)*N_theta+1
c     Interpolate between gridpoints.
               if (j.eq.1) then
                  cogparms(i,3) = lininterp(0d0,0d0,rsum,
     c                 rthetw(lowidx,1),0.25d0)
               else
                  cogparms(i,3) = lininterp(rsum_old,
     c                 rthetw(lowidx-N_theta,1),
     c                 rsum,rthetw(lowidx,1),0.25d0)
               endif
            else
               if (rsum.gt.0.75 .and. upidx.eq.0) then
                  upidx = (j-1)*N_theta+1
c     Interpolate between gridpoints.
                  if (j.eq.1) then
                     cogparms(i,4) = lininterp(0d0,0d0,rsum,
     c                    rthetw(upidx,1),0.75d0)
                  else
                     cogparms(i,4) = lininterp(rsum_old,
     c                    rthetw(upidx-N_theta,1),
     c                    rsum,rthetw(upidx,1),0.75d0)
                  endif
c                 cogparms(i,4) = rthetw(upidx,1)
                  goto 10
               endif
            endif
         enddo

 10      tsum = 0d0
         lowidx = 0
         upidx = 0
         do k=1,N_theta
            tsum_old = tsum
            do j=1,N_rad
               idx = (j-1)*N_theta+k
               tsum = tsum + avker(idx,i) * wint(idx)
            enddo
            if (tsum.gt.0.25 .and. lowidx.eq.0)  then
               lowidx = k
c     Interpolate between gridpoints.
               if (k.eq.1) then
                  cogparms(i,5) = lininterp(0d0,0d0,tsum,
     c                 rthetw(lowidx,2),0.25d0)
               else
                  cogparms(i,5) = lininterp(tsum_old,
     c                 rthetw(lowidx-1,2),
     c                 tsum,rthetw(lowidx,2),0.25d0)
               endif
            endif
            if (tsum.gt.0.75 .and. upidx.eq.0) then 
               upidx = k
c     Interpolate between gridpoints.
               if (k.eq.1) then
                  cogparms(i,6) = lininterp(0d0,0d0,tsum,
     c                 rthetw(upidx,2), 0.75d0)
               else
                  cogparms(i,6) = lininterp(tsum_old,
     c                 rthetw(upidx-1,2),tsum,
     c                 rthetw(upidx,2),0.75d0)
               endif
               goto 20
            endif
         enddo
 20      if (idebug.eq.1) then
         write(*,*) '----------------------------------------------' //
     c           '---------------------------------'
         write(*,*) 'COG parameters are'
         write(*,*) 'r_cg       =',cogparms(i,1),',  theta_cg   =',
     c        cogparms(i,2)*180/pi
         write(*,*) 'r-width    =',cogparms(i,4)-cogparms(i,3),
     c        ',  theta-width=',
     c        (cogparms(i,6)-cogparms(i,5))*180/pi
         write(*,*) 'r_lower    =',cogparms(i,3),',  r_upper    =',
     c        cogparms(i,4)
         write(*,*) 'theta_lower=',cogparms(i,5)*180/pi,
     c         ',  theta_upper=', cogparms(i,6)*180/pi
         write (*,*) 'Integral outside main lobe = ',
     c        1.0-integral(avker(1,i),
     c        cogparms(i,3),cogparms(i,4),cogparms(i,5),cogparms(i,6))
         write(*,*) '------------------------------------------' //
     c        '-------------------------------------'
      endif
      enddo
      return
      end

c     
c**********************************************************************
c     

      subroutine set_qweights(thetarule,weights)     

c     Set integration weights.
c
      implicit none
      include '2dsola.n.include.f'
      integer thetarule
      double precision weights(*)
      double precision rweights(maxrad),tweights(maxtheta)
      double precision r(maxrad),theta(maxtheta),dt,m
      integer i,j

      write (*,*) 'N_theta,N_rad = ',N_theta,N_rad 
      do i=1,N_rad
         r(i) = rthetw((i-1)*N_theta+1,1)
      enddo
      rweights(1) = 0.5*(r(2)-r(1))
      do i=2,N_rad-1
         rweights(i) = 0.5*(r(i+1)-r(i-1))
      enddo
      rweights(N_rad) = 0.5*(r(N_rad)-r(N_rad-1))
      do i=1,N_theta
         theta(i) =  rthetw(i,2)
      enddo
      dt = theta(2)-theta(1)

      if (thetarule.eq.1) then
c     
c     trapezoidal rule.
c
         tweights(1) = 0.5*dt
         do i=2,N_theta-1
            tweights(i) = dt
         enddo
         tweights(N_theta) = 0.5*dt
      else
c
c     Simpsons rule.
         tweights(1) = dt/3.0
         tweights(N_theta) =  dt/3.0
         do i=2,N_theta-1,2
            tweights(i) = 4.0*dt/3.0
         enddo
         do i=3,N_theta-2,2
            tweights(i) = 2.0*dt/3.0
         enddo
      endif
c
      do i=1,N_rad
         do j=1,N_theta
            weights(j+(i-1)*N_theta)=r(i)*tweights(j)*rweights(i)
         enddo
      enddo

      m = abs(weights(1)-rthetw(1,3))
      j = 1
      do i=1,N_rad*N_theta
         if (abs(weights(i)-rthetw(i,3)).gt.m) then
            m = abs(weights(i)-rthetw(i,3))
            j = i
         endif
      enddo
      write (*,*) 'j, max error = ',j,m
      return
      end
      
c
c**********************************************************************
c
      double precision function integral(y,rmin,rmax,tmin,tmax)
      implicit none
      include '2dsola.n.include.f'
      double precision y(*),rmin,rmax,tmin,tmax,sum
      integer j

      sum = 0
      do j=1,N_points
         if (rthetw(j,1).ge.rmin .and. rthetw(j,1).le.rmax .and.
     c        rthetw(j,2).ge.tmin .and. rthetw(j,2).le.tmax) then
            sum = sum + y(j) *  rthetw(j,3)
         endif
      enddo
      integral = sum
      return
      end
c     
c**********************************************************************
c     

      double precision function  lininterp(a,ya,b,yb, x) 
      implicit none
      double precision a,ya,b,yb, x
      lininterp = (ya*(b-x)+yb*(x-a))/(b-a)
      end

