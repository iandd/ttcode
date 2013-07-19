cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine odeint(ny,nlambda,nmu,x1,x2,y1,y2,lambda,mu,disp,h,
     +  hmax,stepmax,tol)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Integrates the set of real ordinary differential equations of a
c real variable, y'(x)=f(x,y), using a 4th/5th-order embedded Runge-
c Kutta method with adaptive stepsize control (Numerical Recipes,
c Section 16.2).  The equations may contain real eigenvalues and
c real parameters.  The function f(x,y) is defined in the subroutine
c "derivs".  The real function "error" supplies a dimensionless local
c error estimate based on the local error in y.  The subroutine
c "display" can be used to display the solution during the integration.

c ny (integer): order of the system (number of functions y)

c nlambda (integer): number of real eigenvalues lambda

c nmu (integer): number of real parameters nmu

c x1 (real): initial value of x

c x2 (real): final value of x

c y1 (real array): initial values of y

c y2 (real array): final values of y (to be returned)

c lambda (real array): real eigenvalues

c mu (real array): real parameters

c disp (logical): set if the subroutine "display" is to be called
c during the integration

c h (real): initial stepsize (negative if x2<x1)

c hmax (real): maximum permitted stepsize (negative if x2<x1)

 5    Ac stepmax (integer): maximum number of (successful) steps permitted

c tol (real): tolerance to be achieved

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer*4 nymax
      parameter (nymax=100)

      logical disp
      integer*4 ny,nlambda,nmu,stepmax
      real*8 x1,x2,mu(nmu),h,hmax,tol
      real*8 y1(ny),y2(ny),lambda(nlambda)

      logical backwards,last,good
      integer*4 i,j
      real*8 x
      real*8 y(nymax)

      if (x1.gt.x2) then
        backwards=.true.
      else if (x1.lt.x2) then
        backwards=.false.
      else
        write (12,'(a13)') 'odeint: x1=x2'
        stop
      endif

      x=x1
      do i=1,ny
        y(i)=y1(i)
      enddo

      do i=1,stepmax
        last=.false.
        if (disp) call display(ny,nlambda,nmu,x,y,lambda,mu)
        if (((.not.backwards).and.((x+h).ge.x2)).or.
     +    (backwards.and.((x+h).le.x2))) then
          h=x2-x
          last=.true.
        endif
        call goodstep(ny,nlambda,nmu,x,y,lambda,mu,h,hmax,tol,good)
        if (last.and.good) then
          do j=1,ny
            y2(j)=y(j)
          enddo
          if (disp) call display(ny,nlambda,nmu,x,y,lambda,mu)
          return
        endif
      enddo

      write (12,'(a22)') 'odeint: too many steps'
      stop

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine goodstep(ny,nlambda,nmu,x,y,lambda,mu,h,hmax,tol,
     +  good)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer*4 nymax
      parameter (nymax=100)

      logical good
      integer*4 ny,nlambda,nmu
      real*8 x,mu(nmu),h,hmax,tol
      real*8 y(ny),lambda(nlambda)

      integer*4 i
      real*8 e,error,xnext
      real*8 ynext(nymax),yerror(nymax)

      good=.true.
 1    call rkstep(ny,nlambda,nmu,x,y,lambda,mu,h,ynext,yerror)
      e=(error(ny,nlambda,nmu,x,y,lambda,mu,h,ynext,yerror)/tol)
      if (e.lt.1.0d0) then
        x=x+h
        do i=1,ny
          y(i)=ynext(i)
        enddo
        h=h*dmin1(5.0d0,dmax1(0.1d0,0.9d0*(e**(-0.2d0))))
        if (dabs(h).gt.dabs(hmax)) h=hmax
        return
      else
        good=.false.
        h=h*dmin1(5.0d0,dmax1(0.1d0,0.9d0*(e**(-0.2d0))))
        xnext=x+h
        if (xnext.eq.x) then
          write (12,'(a25)') 'qstep: stepsize underflow'
          stop
        endif
        goto 1
      endif

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rkstep(ny,nlambda,nmu,x,y,lambda,mu,h,ynext,yerror)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer*4 nymax
      parameter (nymax=100)

      integer*4 ny,nlambda,nmu
      real*8 x,mu(nmu),h
      real*8 y(ny),lambda(nlambda),ynext(ny),yerror(ny)

      integer*4 i
      real*8 a1,a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53,b54,
     +  b61,b62,b63,b64,b65,c1,c2,c3,c4,c5,c6,d1,d2,d3,d4,d5,d6
      real*8 ytemp(nymax),f1(nymax),f2(nymax),f3(nymax),f4(nymax),
     +  f5(nymax),f6(nymax)

      parameter (a1=0.0d0,
     +           a2=(1.0d0/5.0d0),
     +           a3=(3.0d0/10.0d0),
     +           a4=(3.0d0/5.0d0),
     +           a5=1.0d0,
     +           a6=(7.0d0/8.0d0),
     +           b21=(1.0d0/5.0d0),
     +           b31=(3.0d0/40.0d0),
     +           b32=(9.0d0/40.0d0),
     +           b41=(3.0d0/10.0d0),
     +           b42=-(9.0d0/10.0d0),
     +           b43=(6.0d0/5.0d0),
     +           b51=-(11.0d0/54.0d0),
     +           b52=(5.0d0/2.0d0),
     +           b53=-(70.0d0/27.0d0),
     +           b54=(35.0d0/27.0d0),
     +           b61=(1631.0d0/55296.0d0),
     +           b62=(175.0d0/512.0d0),
     +           b63=(575.0d0/13824.0d0),
     +           b64=(44275.0d0/110592.0d0),
     +           b65=(253.0d0/4096.0d0),
     +           c1=(37.0d0/378.0d0),
     +           c2=0.0d0,
     +           c3=(250.0d0/621.0d0),
     +           c4=(125.0d0/594.0d0),
     +           c5=0.0d0,
     +           c6=(512.0d0/1771.0d0),
     +           d1=c1-(2825.0d0/27648.0d0),
     +           d2=c2-0.0d0,
     +           d3=c3-(18575.0d0/48384.0d0),
     +           d4=c4-(13525.0d0/55296.0d0),
     +           d5=c5-(277.0d0/14336.0d0),
     +           d6=c6-(1.0d0/4.0d0))

      call derivs(ny,nlambda,nmu,x+a1*h,y,lambda,mu,f1)
      do i=1,ny
        ytemp(i)=y(i)+h*b21*f1(i)
      enddo
      call derivs(ny,nlambda,nmu,x+a2*h,ytemp,lambda,mu,f2)
      do i=1,ny
        ytemp(i)=y(i)+h*(b31*f1(i)+b32*f2(i))
      enddo
      call derivs(ny,nlambda,nmu,x+a3*h,ytemp,lambda,mu,f3)
      do i=1,ny
        ytemp(i)=y(i)+h*(b41*f1(i)+b42*f2(i)+b43*f3(i))
      enddo
      call derivs(ny,nlambda,nmu,x+a4*h,ytemp,lambda,mu,f4)
      do i=1,ny
        ytemp(i)=y(i)+h*(b51*f1(i)+b52*f2(i)+b53*f3(i)+b54*f4(i))
      enddo
      call derivs(ny,nlambda,nmu,x+a5*h,ytemp,lambda,mu,f5)
      do i=1,ny
        ytemp(i)=y(i)+h*(b61*f1(i)+b62*f2(i)+b63*f3(i)+b64*f4(i)+
     +    b65*f5(i))
      enddo
      call derivs(ny,nlambda,nmu,x+a6*h,ytemp,lambda,mu,f6)
      do i=1,ny
        ynext(i)=y(i)+h*(c1*f1(i)+c2*f2(i)+c3*f3(i)+c4*f4(i)+c5*f5(i)+
     +    c6*f6(i))
        yerror(i)=h*(d1*f1(i)+d2*f2(i)+d3*f3(i)+d4*f4(i)+d5*f5(i)+
     +    d6*f6(i))
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE derivs(ny_in,m,y,f)
      use global_constants
      use set_parameters
      use central_entropy
      IMPLICIT NONE
      integer :: ny_in
      double precision :: m,y(ny_in),f(ny_in)
      double precision :: rho,epsilonpp
      double precision :: xkapP,xkapR,xkapA
      double precision :: dxkapP_dT,dxkapA_dT
      double precision :: gradconv,gradad,gradrad,grad,beta
      double precision :: saumon_rho
      double precision :: mixinglength
c Define f here...
c     m = mass
c     y(1) = r  f(1) = dr/dm
c     y(2) = l  f(2) = dl/dm
c     y(3) = P  f(3) = dP/dm
c     y(4) = T  f(4) = dT/dm
!-- density
      rho = saumon_rho(y(4),y(3),0.7d0,0.3d0)
!-- opacity
      call BurrowsOPC(y(4),rho,xkapP,
     $     xkapR,xkapA,xkap(-1,-1,-1,*),
     $     dxkapP_dT,dxkapA_dT,-1,-1,-1)
!-- nuclear energy
c      call ppchain(epsilonpp,rho,y(4),XX)
      epsilonpp = 0.d0
!-- dr/dm
      f(1) = 1.d0/(4.d0*PI*(y(1)**2.d0)*rho)
!-- dl/dm
      f(2) = epsilonpp-y(4)*DS_Dt
!-- dP/dm
      f(3) = -GRAV*m/(4.d0*PI*(y(1)**4.d0))
!-- set beta =1 (ie no radiation pressure)
      beta = 1.d0
!-- dT/dm (check which gradient to use...)
      gradad = (1.d0 + ((1.d0-beta)*(4.d0+beta)/(beta**2.d0)))/
     %     ((5.d0/2.d0) + 4.d0*((1.d0-beta)*(4.d0+beta)/(beta**2.d0)))
      gradrad = 3.d0*xkapR*y(2)*y(3)/(16.d0*PI*GRAV*ARAD*FD*
     +     (y(4)**4.d0)*m)
      gradconv =  mixinglength(m,ny,y,gradrad,gradad,XX,YY)
      if(gradrad.lt.gradconv) then
         grad=gradrad
c         print *,'radgrad',gradrad,gradconv
c         stop
      endif
      if(gradrad.ge.gradconv) then
         grad=gradconv
      endif
      f(4) = -(GRAV*m*y(4)/(4.d0*PI*(y(1)**4.d0)*y(3)))*grad
c      write(6,'(a,6(e14.6))') 'der:l',y(2),f(2),f(4)
c      write(6,'(a,6(e14.6))') 'der',f(1),f(2),f(3),f(4)
      RETURN
      END SUBROUTINE derivs


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function error(ny,nlambda,nmu,x,y,lambda,mu,h,ynext,yerror)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer*4 ny,nlambda,nmu
      real*8 error,x,mu(nmu),h
      real*8 y(ny),lambda(nlambda),ynext(ny),yerror(ny)

      integer*4 i
      real*8 ymax,yemax

c Define error estimate appropriate to application...

      ymax=0.0d0
      yemax=0.0d0
      do i=1,ny
        ymax=dmax1(ymax,dabs(y(i)))
        yemax=dmax1(yemax,dabs(yerror(i)))
      enddo
      error=yemax/ymax

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine display(ny,nlambda,nmu,x,y,lambda,mu)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer*4 ny,nlambda,nmu
      real*8 x,mu(nmu)
      real*8 y(ny),lambda(nlambda)

c Write required variables...

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
