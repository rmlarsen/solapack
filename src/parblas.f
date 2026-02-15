c     Copyright Rasmus Munk Larsen, Stanford University, 2003

      double precision function pdnrm2(n, x, incx)
      implicit none
      integer n, incx
      double precision x(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0)) then    
         sum = 0d0
c$doacross local(i) shared(x,n,incx) reduction(sum)
cc$PAR DOALL readonly(x,incx), reduction(sum)
         do i=1,n
            sum = sum + x(1+(i-1)*incx)**2
         enddo
         pdnrm2 = sqrt(sum)
      else
         pdnrm2 = 0d0
      endif   
      return
      end
c
c****************************************************************************
c 
      

      subroutine pdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$doacross local(i) shared(x,alpha) 
cc$PAR DOALL shared(x), readonly(alpha)
            do i=1,n
               x(i) = alpha*x(i)
            enddo
         else
c$doacross local(i) shared(x,n,incx,alpha)
cc$PAR DOALL  shared(x), readonly(alpha,incx)
            do i=1,n
               x(1+(i-1)*incx) = alpha*x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end
      
c
c****************************************************************************
c 
         

      subroutine pdcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$doacross local(i) shared(x,y,n) 
cc$PAR DOALL  shared(y), readonly(x)
            do i=1,n
               y(i) = x(i)
            enddo
         else
c$doacross local(i) shared(x,y,n,incx,incy)
cc$PAR DOALL  shared(y), readonly(x,incx,incy)
            do i=1,n
               y(1+(i-1)*incy) = x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 
      subroutine pdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha,x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$doacross local(i) shared(x,y,n) 
cc$PAR DOALL  shared(y), readonly(x,alpha)
            do i=1,n
               y(i) = alpha*x(i) + y(i)
            enddo
         else
c$doacross local(i) shared(x,y,n,incx,incy)
cc$PAR DOALL  shared(y), readonly(alpha,x,incx,incy)
            do i=1,n
               y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c              y(1+(i-1)*incy)
            enddo
         endif
      endif
      return
      end
      
c
c****************************************************************************
c 
         
      subroutine pdzero(n, x, incx)
      implicit none
      integer n, incx
      double precision x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$doacross local(i) shared(x,n) 
cc$PAR DOALL  shared(x)
            do i=1,n
               x(i) = 0d0
            enddo
         else
c$doacross local(i) shared(x,n,incx)
cc$PAR DOALL  shared(x), readonly(incx)
            do i=1,n
               x(1+(i-1)*incx) = 0d0
            enddo
         endif
      endif
      return
      end


      subroutine pdset(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$doacross local(i) shared(x,n,alpha) 
cc$PAR DOALL  shared(x), readonly(alpha)
            do i=1,n
               x(i) = alpha
            enddo
         else
c$doacross local(i) shared(x,n,incx,alpha)
cc$PAR DOALL  shared(x), readonly(alpha,incx)
            do i=1,n
               x(1+(i-1)*incx) = alpha
            enddo
         endif
      endif
      return
      end


c
c****************************************************************************
c 
     
         
      double precision function pddot(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = 0d0
c$doacross local(i) shared(x,y,n) reduction(sum)
cc$PAR DOALL  readonly(x,y), reduction(sum)
            do i=1,n
               sum = sum + x(i) * y(i)
            enddo
         else
            sum = 0d0
c$doacross local(i) shared(x,y,n,incx,incy) reduction(sum)
cc$PAR DOALL  readonly(x,y,incx,incy), reduction(sum)
            do i=1,n
               sum = sum + x(1+(i-1)*incx) * y(1+(i-1)*incy)
            enddo
         endif
         pddot = sum
      else
         pddot = 0d0
      endif   
      return
      end


  
      subroutine pdaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,beta,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero .and. beta.eq.zero) then
         if (incy.eq.1) then
c$doacross shared(y) 
cc$PAR DOALL shared(y)
            do i=1,n
               y(i) = zero
            enddo
         else
c$doacross shared(y,incy)
cc$PAR DOALL shared(y), readonly(incy)
            do i=1,n
               y(1+(i-1)*incy) = zero
            enddo
         endif
         
      else if (alpha.eq.zero .and. beta.ne.zero) then
         
         call pdscal(n,beta,y,incy)

      else if (alpha.ne.zero .and. beta.eq.zero) then

         if (alpha.eq.one) then
            call pdcopy(n,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$doacross shared(alpha,x,y) 
cc$PAR DOALL shared(y), readonly(alpha,x) 
               do i=1,n
                  y(i) = alpha*x(i)
               enddo
            else
c$doacross shared(alpha,x,y,incx,incy)
cc$PAR DOALL shared(y), readonly(alpha,x,incx,incy)
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)
               enddo
            endif
         endif

      else

         if (beta.eq.one) then
c DAXPY
            call pdaxpy(n,alpha,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$doacross shared(alpha,x,y,n) 
cc$PAR DOALL shared(y), readonly(alpha,x,n) 
               do i=1,n
                  y(i) = alpha*x(i) + beta*y(i)
               enddo
            else
c$doacross shared(alpha,beta,x,y,n,incx,incy)
cc$PAR DOALL shared(y), readonly(alpha,beta,x,n,incx,incy)
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c                 beta*y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end



      subroutine pdaxty(n,alpha,x,incx,y,incy)
c
c     Y = alpha*X*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero) then
         if (incy.eq.1) then
c$doacross shared(y) 
cc$PAR DOALL shared(y)
            do i=1,n
               y(i) = zero
            enddo
         else
c$doacross shared(y,incy)
cc$PAR DOALL shared(y), readonly(incy)
            do i=1,n
               y(1+(i-1)*incy) = zero
            enddo
         endif
         
      else if (alpha.ne.zero) then

         if (alpha.eq.one) then
            if (incx.eq.1 .and. incy.eq.1) then
c$doacross shared(alpha,x,y) 
cc$PAR DOALL shared(y), readonly(x) 
               do i=1,n
                  y(i) = x(i)*y(i)
               enddo
            else
c$doacross shared(alpha,x,y,incx,incy)
cc$PAR DOALL shared(y), readonly(alpha,x,incx,incy)
               do i=1,n
                  y(1+(i-1)*incy) = x(1+(i-1)*incx)*y(1+(i-1)*incy)
               enddo
            endif

         else
            if (incx.eq.1 .and. incy.eq.1) then
c$doacross shared(alpha,x,y) 
cc$PAR DOALL shared(y), readonly(alpha,x) 
               do i=1,n
                  y(i) = alpha*x(i)*y(i)
               enddo
            else
c$doacross shared(alpha,x,y,incx,incy)
cc$PAR DOALL shared(y), readonly(alpha,x,incx,incy)
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)*
     c                 y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end
