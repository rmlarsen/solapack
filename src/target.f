c     Copyright Rasmus Munk Larsen, Stanford University, 2003

c     
c***************************************************************
c     
      subroutine def_targets(ids,itargt,targetparms,trdoff)
c     
c     Set up target parameters for SOLA inversion, defined by 
c     parameters read from unit ids.
c     
c     The definition of the targets is controlled by the parameter
c     itargt read first:
c     
c     itargt = 0: read individual target positions and widths, and
c     trade-off parameters.
c     The subsequent lines have the format
c     x0, theta0, deltax, deltatheta, trdoff
c     itargt > 0: set up targets based on kernel mesh or otherwise
c     on square grid.
c     The subsequent lines have the format
c     deltax, deltatheta, tradeoff   (assumed to be fixed)
c     x1, x2, theta1, theta2, nxtar, nttar
c     
c     In this case has been implemented:
c     
c     itargt = 1: choose points in kernel mesh with x 
c     between x1 and x2, theta between theta1 and theta2
c     skipping every nxtar point in radius and every nttar point
c     in latitude.
c     
c     NOTE: If itargt=1 then N_rad and N_theta must be defined upon entry.
c     
      implicit double precision (a-h, o-z)
      include '2dsola.n.include.f'
      double precision  targetparms(maxtargets,*),trdoff(*)
      integer i,j
c     
c     read control parameter
c     
c     read(ids,end=90,err=90) itargt
      write(*,*) 'Input itargt'
      read(ids,*) itargt
      write(6,*) itargt
      if ((itargt.lt.0).or.(itargt.gt.1)) goto 90
c     
      if(itargt.eq.0) then
c     
c     start loop over target and tradeoff parameters
c     
         N_targets = 0
 10      N_targets = N_targets + 1
         read(ids,*,end=20) (targetparms(N_targets,j),j=1,4), 
     c        trdoff(N_targets)
c     write(6,*)  (targetparms(N_targets,j),j=1,4), trdoff(N_targets)
c     
c     loop back for new target parameters
c     
         if (N_targets.gt.maxtargets) then
            write (6,*) 'Too many targets.'
            write (6,*) 'Increase maxtargets and re-compile.'
            stop 'Out of workspace'
         endif
         go to 10
c     
 20      continue
         N_targets = N_targets-1       
         return
c     
      else
c     
         read(ids,*,end=80,err=90) deltxr,delttr,trdofr
         write(6,*) deltxr,delttr,trdofr
         read(ids,*,end=80,err=90) x1, x2, theta1, theta2, nxtar, nttar
         write (6,*) x1, x2, theta1, theta2, nxtar, nttar
c     
      end if
c     
c     test for version of setting up
c     
      if(itargt.eq.1) then
c     
c     set targets from kernel mesh
c     
         if (nxtar.lt.1) nxtar = 1
         if (nttar.lt.1) nttar = 1
         N_targets=0
         do i=1,N_rad,nxtar
            do j=1,N_theta,nttar
               n = N_theta*(i-1)+j
               if(x1.le.rthetw(n,1) .and. rthetw(n,1).le.x2 .and.
     *              theta1.le.rthetw(n,2) .and. rthetw(n,2).le.theta2)
     *              then
                  N_targets=N_targets+1
                  targetparms(N_targets,1)=rthetw(n,1)
                  targetparms(N_targets,2)=rthetw(n,2)
                  targetparms(N_targets,3)=deltxr
                  targetparms(N_targets,4)=delttr
                  trdoff(N_targets)=trdofr
                  if (N_targets.gt.maxtargets) then
                     write (6,*) 'Too many targets.'
                     write (6,*) 'Increase maxtargets and re-compile.'
                     stop 'Out of workspace'
                  endif
               end if
            enddo
         enddo
c     
      else
         write(6,110) itargt
         stop 'target not defined'
      end if
      return
 80   write(6,130) itargt
      return
c     
c     diagnostics
c     
 90   write(6,120) itargt
      stop 'input error'
c     
 110  format(///' ***** Error in s/r def_targets. itargt =',i5,
     *     ' not implemented')
 120  format(///' ***** Input error in s/r def_targets, for itargt =',
     *     i5)
 130  format(///' ***** end of file reached in def_targets, for ' //
     *     ' itargt =',i5)
      end
c     
c     
c***************************************************************
c     
      subroutine set_target(x0, theta0, deltax, 
     *     deltat, target, icase, itype)
c     
c     Set target function for 2D rotational SOLA inversion,
c     at location (x0, theta0) and widths (deltax, deltat).
c     
c     icase = icase0 + 100*icase1. 
c     Here icase0 determines target function.
c     icase1 = 1 flags for setting kernel for radial derivative.
c     
c     The current version only allows the following target functions:
c     icase0 = 0: Gaussian target, symmetrized around the
c     equator and normalized to have integral 1 over half-disk,
c     i.e., integral 0.5 over quarter disk.
c     icase1 = 0: modified Gaussian target, with zero second moment in
c     both radial and latitude directions
c     
      implicit none
      include '2dsola.n.include.f'
      integer icase
      double precision x0,theta0, deltax, deltat
      double precision target(*)
c     
c     Local variables
c     
      integer i,icase0,itype
      double precision dx,dx2,dt,dt2,dts,dts2,PI,C,x,t,xi
      double precision integral,starget(maxpts),
     c     weight(maxpts),delta,rmax,rmin,piby2

      parameter(PI = 3.14159265358979323846D0, C = 2.0D0/3.0D0)

      external daxpy,dscal,ddot
      double precision ddot
      
      piby2 = PI/2.0
      rmax = 0.998
      rmin = 0.1

      if (itype.eq.0) then
         icase0=mod(icase,100)        
         do i=1,N_points
            dx=(rthetw(i,1)-x0)/deltax
            dx2=dx*dx
            dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
            dt2=dt*dt
            dts=(rthetw(i,2)-(PI-theta0))* rthetw(i,1)/deltat
            dts2=dts*dts    
            target(i) = -dx2-dt2
            starget(i) = -dx2-dts2
         enddo

         do i=1,N_points
            starget(i) = exp(starget(i))
            target(i) = exp(target(i))
         enddo

c     
c     calculate integral of symmetrized target.
c     
         integral = 0d0
         if(icase0.eq.0) then
            do i=1,N_points
               integral=integral + rthetw(i,3)*(starget(i)+target(i))
            enddo
         else            
            do i=1,N_points
               dx=(rthetw(i,1)-x0)/deltax
               dx2=dx*dx
               dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
               dt2=dt*dt
               dts=(rthetw(i,2)-(PI-theta0))* rthetw(i,1)/deltat
               dts2=dts*dts    
               target(i)=(1-C*dx2) * ( (1-C*dt2) * target(i) + 
     c              (1-C*dts2) * starget(i) )
               integral=integral + rthetw(i,3)*target(i)
            enddo
         end if
c     
c     normalize target
c     
         do i=1,N_points
            target(i) = target(i)/integral
            starget(i) = starget(i)/integral
         enddo
c     
c     
c     test for setting kernel for rotation rate or its derivative
c     
         
         if(icase.eq.100) then
c     
c     reset target, but including normalization, for derivative
            do  i=1,N_points                     
               x = rthetw(i,1)
               t = rthetw(i,2)
               dx=(x - x0)/deltax
               dx2=dx*dx
               dt=(t - theta0)/deltat
               dt2=dt*dt
               dts=(t-(PI-theta0))/deltat
               dts2=dts*dts 
               if (x.ne.0d0) then
                  xi = 1d0/x
               else
                  xi = 0d0
               endif
               target(i) = (2*dx/deltax+2*dt2*x-xi)*target(i)+
     c              (2*dx/deltax + 2*dts2*x - xi)*starget(i)
            enddo
         else if (icase.eq.200) then
c     reset target, but including normalization, for 2. derivative
            do  i=1,N_points
               x = rthetw(i,1)
               t = rthetw(i,2)
               dx=(x-x0)/deltax
               dx2=dx*dx
               dt=(t-theta0)/deltat
               dt2=dt*dt
               dts=(t-(PI-theta0))/ deltat
               dts2=dts*dts    
               if (x.ne.0d0) then
                  xi = 1d0/x
               else
                  xi = 0d0
               endif
               target(i) = ( 2*dt2*(2*dt2*x*x-3) +
     c              2/(deltax*deltax)*(2*dx2-2*deltax*dx*xi-1) +
     c              8*dx/deltax*dt2*x)*target(i) +
     c              ( 2*dts2*(2*dts2*x*x-3) +
     c              2/(deltax*deltax)*(2*dx2-2*deltax*dx*xi-1) +
     c              8*dx/deltax*dts2*x)*starget(i)
            enddo
         else
            do i=1,N_points
               target(i) = target(i) + starget(i)
            enddo
         end if
      elseif (itype.eq.2) then
         integral = 0d0
         do i=1,N_points
            dx=2.0*dsign(max(abs((rthetw(i,1)-x0)/deltax),1d-8),
     c           (rthetw(i,1)-x0))
            dt=2.0*dsign(max(abs((rthetw(i,2)-theta0)*rthetw(i,1)/
     c           deltat), 1d-8), rthetw(i,2)-theta0)   
            dts=2.0*dsign(max(abs((rthetw(i,2)-(PI-theta0))*
     c           rthetw(i,1)/
     c           deltat), 1d-8), rthetw(i,2)-(PI-theta0))   
            target(i) = sin(dx)/dx*(sin(dt)/dt+sin(dts)/dts)
c     target(i) = sin(sqrt(dx**2+dt**2))/sqrt(dx**2+dt**2) +
c     c            sin(sqrt(dx**2+dts**2))/sqrt(dx**2+dts**2)
            integral = integral + target(i)*rthetw(i,3)
         enddo
         do i=1,N_points
            target(i) = target(i)/integral
         enddo
      elseif (itype.eq.3) then
         delta = max(deltat,theta0)
         do i =1,N_points
            dx=(rthetw(i,1)-x0)/deltax
            dx2=dx*dx
            dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
            dt2=dt*dt
            dts=(rthetw(i,2)-(PI-theta0))* rthetw(i,1)/deltat
            dts2=dts*dts    
            target(i) = -dx2-dt2
            starget(i) = -dx2-dts2
            if (rthetw(i,2).lt.delta*pi/2.0) then
               weight(i) = sin(rthetw(i,2)/delta)**2
            else
               weight(i) = 1d0
            endif
         enddo
         
         do i=1,N_points
            starget(i) = exp(starget(i))
            target(i) = exp(target(i))
         enddo

         integral = 0d0
         do i=1,N_points
            target(i) = (target(i) + starget(i))*weight(i)
            integral = integral + rthetw(i,3)*target(i)
         enddo
         call dscal(N_points,1D0/integral,target,1)         
      elseif (itype.eq.4) then
         delta = 0.2
         do i=1,N_points
            dx=(rthetw(i,1)-x0)/deltax
            dx2=dx*dx
            dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
            dt2=dt*dt
            dts=(rthetw(i,2)-(PI-theta0))* rthetw(i,1)/deltat
            dts2=dts*dts    
            target(i) = -dx2-dt2
            starget(i) = -dx2-dts2
            if (rthetw(i,2).lt.delta) then
               weight(i) = sin(piby2*rthetw(i,2)/delta)**4
            else
               weight(i) = 1d0
            endif
            if (rthetw(i,1).gt.rmax) then
               if (rthetw(i,1).gt.1d0) then
                  weight(i) = 0d0
               else
                  weight(i) = weight(i)*sin(piby2*(1.0-rthetw(i,1)) / 
     c                 (1.0-rmax))**2
               endif
            elseif  (rthetw(i,1).lt.rmin)  then
               if (rthetw(i,1).eq.0d0) then
                  weight(i) = 0d0
               else
                  weight(i) = weight(i)*sin(piby2*rthetw(i,1) /
     c                 rmin)**4
               endif
            endif
         enddo
         
         do i=1,N_points
            starget(i) = exp(starget(i))
            target(i) = exp(target(i))
         enddo

         integral=0d0
         do i=1,N_points
            target(i) = (target(i) + starget(i))*weight(i)
            integral=integral + rthetw(i,3)*target(i)
         enddo
         call dscal(N_points,1D0/integral,target,1)         
      endif
      return
      end
c     
c***************************************************************
c     
c     
      subroutine set_targets(targetparms, target, ldt, icase,itype)
c     
c     set  'N_targets' targets for 2D rotational SOLA inversion,
c     at location (x0(:), theta0(:)) and widths (deltax(:), deltat(:))
c     
c     icase = icase0 + 100*icase1. 
c     Here icase0 determines target function.
c     icase1 = 1 flags for setting kernel for radial derivative.
c     
c     The current version only allows the following target functions:
c     icase0 = 0: Gaussian target, symmetrized around the
c     equator and normalized to have integral 1 over half-disk,
c     i.e., integral 0.5 over quarter disk.
c     icase1 = 0: modified Gaussian target, with zero second moment in
c     both radial and latitude directions
c     
      implicit none
      include '2dsola.n.include.f'
      integer icase,ldt
      double precision targetparms(maxtargets,*)
      double precision target(ldt,*)
c     Local variables
      integer i,j,icase0,itype
      double precision x0,theta0,deltax,deltat,xi,x,t
      double precision dx,dx2,dt,dt2,dts,dts2,PI,C
      double precision integral,starget(maxpts),power,
     c     weight(maxpts),delta,rmax,rmin,piby2

      parameter(PI = 3.14159265358979323846D0, C = 2.0D0/3.0D0)

      external daxpy,dscal,ddot
      double precision ddot
      
      if (N_targets.gt.0) then
         piby2 = PI/2.0
         icase0=mod(icase,100)
         rmax = 0.998
         rmin = 0.1
         
c$doacross local(i,j,dx,dx2,dt,dt2,dts,dts2,integral,starget,weight),
c$& local(x0,theta0,deltax,deltat,delta,power,x,t,xi),
c$& shared(target,rthetw,targetparms,N_targets,N_points,icase) 
cc$PAR DOALL private(j,dx,dx2,dt,dt2,dts,dts2,integral,starget,weight)
cc$PAR DOALL private(x0,theta0,deltax,deltat,delta,power,x,t,xi)
cc$PAR DOALL shared(target)
cc$PAR DOALL readonly(rthetw,targetparms,N_points,icase)
         do j=1,N_targets
            x0 = targetparms(j,1)
            theta0 = targetparms(j,2)
            deltax = targetparms(j,3)
            deltat = targetparms(j,4)
c            write (*,*) 'x0 ,theta0, deltax ,deltat ',x0 ,theta0, deltax ,deltat
            if (itype.eq.0) then
               do i=1,N_points
                  dx=(rthetw(i,1)-x0)/deltax
                  dx2=dx*dx
                  dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
                  dt2=dt*dt
                  dts=(rthetw(i,2)-(PI-theta0))* rthetw(i,1)/deltat
                  dts2=dts*dts    
                  target(i,j) = -dx2-dt2
                  starget(i) = -dx2-dts2
               enddo
               
               do i=1,N_points
                  starget(i) = exp(starget(i))
                  target(i,j) = exp(target(i,j))
               enddo
c     
c     calculate integral of symmetrized target.
c     
               integral=0d0
               if(icase0.eq.0) then
                  do i=1,N_points
                     integral=integral + rthetw(i,3)*(starget(i) + 
     c                    target(i,j))
                  enddo
               else            
                  do i=1,N_points
                     dx=(rthetw(i,1)-x0)/deltax
                     dx2=dx*dx
                     dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
                     dt2=dt*dt
                     dts=(rthetw(i,2)-(PI-theta0))* rthetw(i,1)/deltat
                     dts2=dts*dts    
                     target(i,j)=(1-C*dx2)*( (1-C*dt2)*target(i,j) + 
     c                    (1-C*dts2)*starget(i) )
                     integral=integral + rthetw(i,3)*target(i,j)
                  enddo
               end if
c     
c     normalize target
c     
c     write (*,*) 'integral = ',integral
               do i=1,N_points
                  target(i,j) = target(i,j)/integral
                  starget(i) = starget(i)/integral
               enddo
c     
c     test for setting kernel for rotation rate or its derivative
c     
               if (icase.eq.100) then
c     
c     reset target, but including normalization, for 1. derivative
c     
c                  write(*,*) 'blah100'
                  do  i=1,N_points                     
                     x = rthetw(i,1)
                     t = rthetw(i,2)
                     dx=(x - x0)/deltax
                     dx2=dx*dx
                     dt=(t - theta0)/deltat
                     dt2=dt*dt
                     dts=(t-(PI-theta0))/deltat
                     dts2=dts*dts 
                     if (x.ne.0d0) then
                        xi = 1d0/x
                     else
                        xi = 0d0
                     endif
                     target(i,j) = (2*dx/deltax+2*dt2*x-xi) * 
     c                    target(i,j)+
     c                    (2*dx/deltax + 2*dts2*x - xi)*starget(i)
                  enddo

               else if (icase.ge.200) then
c     
c     reset target, but including normalization, for 2. derivative
c     
c                  write(*,*) 'blah200'
                  do  i=1,N_points
                     x = rthetw(i,1)
                     t = rthetw(i,2)
                     dx=(x-x0)/deltax
                     dx2=dx*dx
                     dt=(t-theta0)/deltat
                     dt2=dt*dt
                     dts=(t-(PI-theta0))/ deltat
                     dts2=dts*dts    
                     if (x.ne.0d0) then
                        xi = 1d0/x
                     else
                        xi = 0d0
                     endif
                     target(i,j) = ( 2*dt2*(2*dt2*x*x-3) +
     c                    2/(deltax*deltax)*(2*dx2-2*deltax*dx*xi-1) +
     c                    8*dx/deltax*dt2*x)*target(i,j) +
     c                    ( 2*dts2*(2*dts2*x*x-3) +
     c                    2/(deltax*deltax)*(2*dx2-2*deltax*dx*xi-1) +
     c                    8*dx/deltax*dts2*x)*starget(i)
                  enddo
               else
                  do i=1,N_points
                     target(i,j) = target(i,j) + starget(i)
                  enddo
               end if
            elseif (itype.eq.2) then
               integral = 0d0
               do i=1,N_points
                  dx=2.0*dsign(max(abs((rthetw(i,1)-x0)/deltax),1d-8),
     c                 (rthetw(i,1)-x0))
                  dt=2.0*dsign(max(abs((rthetw(i,2)-theta0)*rthetw(i,1)
     c                 /deltat), 1d-8), rthetw(i,2)-theta0)   
                  dts=2.0*dsign(max(abs((rthetw(i,2)-(PI-theta0)) * 
     c                 rthetw(i,1)/deltat),1d-8), 
     c                 rthetw(i,2)-(PI-theta0))   
c     dx=(rthetw(i,1)-x0)/deltax
c     dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
c     dts=(rthetw(i,2)+theta0 - PI)*rthetw(i,1)/deltat
                  target(i,j) = sin(dx)/dx*(sin(dt)/dt+sin(dts)/dts)
c     target(i,j) = sin(sqrt(dx**2+dt**2))/sqrt(dx**2+dt**2) +
c     c                 sin(sqrt(dx**2+dts**2))/sqrt(dx**2+dts**2)
                  integral = integral + target(i,j)*rthetw(i,3)
               enddo
               do i=1,N_points
                  target(i,j) = target(i,j)/integral
               enddo
            elseif (itype.eq.3) then
               delta = max(deltat,theta0)
               do i=1,N_points
                  dx=(rthetw(i,1)-x0)/deltax
                  dx2=dx*dx
                  dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
                  dt2=dt*dt
                  dts=(rthetw(i,2)-(PI-theta0))* rthetw(i,1)/deltat
                  dts2=dts*dts    
                  target(i,j) = -dx2-dt2
                  starget(i) = -dx2-dts2
               enddo
               
               do i=1,N_points
                  weight(i) = sin(rthetw(i,2))
                  starget(i) = exp(starget(i))
                  target(i,j) = exp(target(i,j))
               enddo

               integral=0d0
               do i=1,N_points
                  target(i,j) = (target(i,j) + starget(i)) * 
     c                 weight(i)**(1.5)
                  integral=integral + rthetw(i,3)*target(i,j)
               enddo
               call dscal(N_points,1D0/integral,target(1,j),1)

            elseif (itype.eq.4) then
               delta = 0.25
               power = 9
               do i=1,N_points
                  dx=(rthetw(i,1)-x0)/deltax
                  dx2=dx*dx
                  dt=(rthetw(i,2)-theta0)*rthetw(i,1)/deltat
                  dt2=dt*dt
                  dts=(rthetw(i,2)-(PI-theta0))* rthetw(i,1)/deltat
                  dts2=dts*dts    
                  target(i,j) = -dx2-dt2
                  starget(i) = -dx2-dts2
                  if (rthetw(i,2).lt.delta) then
                     weight(i) = sin(piby2*rthetw(i,2)/delta)**power
                  else
                     weight(i) = 1d0
                  endif
                  if (rthetw(i,1).gt.rmax) then
                     if (rthetw(i,1).gt.1d0) then
                        weight(i) = 0d0
                     else
                        weight(i) = weight(i)*sin(piby2* (1.0 - 
     c                       rthetw(i,1)) / (1.0-rmax))**2
                     endif
                  elseif  (rthetw(i,1).lt.rmin)  then
                     if (rthetw(i,1).eq.0d0) then
                        weight(i) = 0d0
                     else
                        weight(i) = weight(i)*sin(piby2*rthetw(i,1) / 
     c                       rmin)**4
                     endif
                  endif
               enddo
               
               do i=1,N_points
                  starget(i) = exp(starget(i))
                  target(i,j) = exp(target(i,j))
               enddo

               integral=0d0
               do i=1,N_points
                  target(i,j) = (target(i,j) + starget(i))*weight(i)
                  integral=integral + rthetw(i,3)*target(i,j)
               enddo
               call dscal(N_points,1D0/integral,target(1,j),1)
            endif
         enddo
      endif
      return
      end
c     
c     
c***************************************************************
c     
      subroutine scale_rwidth(ibyteswap,amdlfile,x0,xrscal,deltax)
c     
c     Scale target width provided in deltax with sound speed at
c     target location x0, relative to sound speed at r/R = xrscal.
c     ids is unit number of file containing amdl model, which must
c     have been opened before call. Model is read in first call.
c     
c     Original version: 14/8/96
c     
      implicit double precision (a-h, o-z)
      character*256 amdlfile
      parameter(nnmax=4000,ia=8)
      double precision cgrav,crscal
      integer init
      dimension x(nnmax),aa(ia,nnmax),data(8),cint(1)
c     
      save
      data init / 0 /
      data cgrav / 6.67232d-8 /

c     
c     test for initializing model, setting dimensionless sound speed into
c     aa(6,.)
c     
      if(init.eq.0) then
         call c_rdamdl(ibyteswap,x,aa,data,nn,nmod,ivar,
     c        icry,ia,amdlfile)
         call flush(6)
         if(icry.ne.0) stop 'scale_radwidth'
c     
         do n=2,nn
            aa(6,n)=x(n)*sqrt(aa(1,n)/aa(2,n))
         enddo
         aa(6,1)=sqrt(data(2)*data(3)*aa(3,1)/(cgrav*data(1)*data(4)))
c     write(81,'(0pf10.6,1p3e13.5)') 
c     *    (x(n),aa(1,n),aa(2,n),aa(6,n),n=1,nn)
c     
c     find value at reference radius
c     
         lint=0
         call dlir(xrscal,x,cint,aa(6,1),1,ia,nn,lint,inter)
         crscal=cint(1)
c     write(6,*) 'crscal set to',crscal
         
         init=1
c     
      end if
c     
c     interpolate to actual target location
c     
      
      call dlir(x0,x,cint,aa(6,1),1,ia,nn,lint,inter)
      fact=cint(1)/crscal
      deltax=fact*deltax
c     write(6,*) 'fact, deltax set to',fact,deltax
      
      return
      end
c     
c***************************************************************



      subroutine   dlir(z,zi,y,yi,ii,id,nt,l,inter)
c     subroutine  dlir1(z,zi,y,yi,ii,id,nt,l,inter)
c     
c     
c     
c     interpolation/extrapolation routine
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
c     starts at n=1
c     if l.gt.1, scan starts from value of n from previous call of dlir
c     
c     
c     dlir use cubic interpolation/extrapolation unless nt.lt.4
c     dlir1 use linear interpolation/extrapolation
c     
c     
c     Double precision version.
c     ++++++++++++++++++++++++
c     
c     Dated 10/3/90.
c     
c     Note: this double precision version of the routine has same name
c     as single precision version, but is distinguished by its file name
c     
c     
      implicit double precision(a-h,o-z)
      dimension zi(*),y(*),yi(*),a(4)
      data n/-1/
c     
      il=0
      go to 1
      entry dlir1(z,zi,y,yi,ii,id,nt,l,inter)
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
 7       if(y(i).eq.0.) y(i+ir1)=0.
         go to 30
c     
c     control when z does not lie on a mesh point
 8       inter=0
 9       if(n.le.1) inter=0
         if(il.eq.1) go to 20
c     
c     cubic interpolation/extrapolation
c     
c     pivotal point (m) and point (k) closest to z
         m=n
         k=3
         if(n.gt.1+ir) go to 11
         m=1+ir+ir
         k=n
 11      if(n.lt.j) go to 12
         m=j-ir
         k=4
c     
c     
c     weighting factors
 12      y1=zi(m-ir*2)
         y2=zi(m-ir)
         y3=zi(m)
         y4=zi(m+ir)
c     
         z1=z-y1
         z2=z-y2
         z3=z-y3
         z4=z-y4
c     
         z12=z1*z2
         z34=z3*z4
c     
         a(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
         a(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
         a(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
         a(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))
c     
c     correct a(k)
         diff=a(1)+a(2)+a(3)+a(4)
         a(k)=(1.d0+a(k))-diff
c     
c     compute y
         m=(m-1)/ir-3
         m=m*ird
         do 18 i=1,iir
            k=i+m
            yy=0.
            do 17 j=1,4
               k=k+ird
               diff=yi(k)
 17            yy=yy+a(j)*diff
               y(i)=yy
 18            if(y(i).eq.0.) y(i+ir1)=0.
               go to 30
c     
c     linear interpolation/extrapolation
 20            if(n.eq.1) n=1+ir
               if(n.gt.j) n=j
               z1=zi(n)
               y1=(z1-z)/(z1-zi(n-ir))
               y2=1.0-y1
               j=(n-1)*id
               m=j-ird
               do 21 i=1,iir,ir
                  y(i)=y1*yi(i+m)+y2*yi(i+j)
 21               if(y(i).eq.0.) y(i+ir1)=0.
c     
c     reset n
 30               n=(n+ir-1)/ir
                  return
c     
c     
c     diagnostics
 101              write(6,1001) nt
                  return
 102              write(6,1002) zi(1),nt,zi(j)
                  return
c     
 1001   format(/1x,10('*'),5x,'there are fewer than two data points in',
     *                 ' dlir     nt =',i4,5x,10('*')/)
 1002   format(/1x,10('*'),5x,'extreme values of independent variable',
     *    ' equal in dlir',5x,10('*')/16x,'zi(1) =',1pe13.5,',   ',
     *    'zi(',i4,') =',1pe13.5/)
c     
        end
