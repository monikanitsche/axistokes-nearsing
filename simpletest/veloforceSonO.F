c Contains
c      SUBROUTINE veloforceSonO : compute single layer at target pt x,y
c                                 due to sheet at xs,ys

#include "flags.h"
      SUBROUTINE veloforceSonO(x,y,n,xs,ys,ns,f1,f2,
     +       xd0,xddd0,ydd0,f1d0,f2dd0,
     +       xdn,xdddn,yddn,f1dn,f2ddn,
     +       u,v)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc     Date: 1/28/18
c
c     Computes velocity u,v at sheet x,y,n induced by sheet xs,ys,ns,f1,f2
c     ASSUME: 1) xs(0),xs(n) on axis. 
c             2) x(0),x(n) on axis
c             3) x,y and xs ys are two different sheets! so no local corrections
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,k,n,ns
c INPUT/OUTPUT
      MODE x(0:*),y(0:*),xs(0:*),ys(0:*),f1(0:*),f2(0:*),
     +      xd0,xdn,xddd0,xdddn,ydd0,yddn,f1d0,f1dn,f2dd0,f2ddn,
     +      u(0:*),v(0:*)
      MODE F(0:nmax),E(0:nmax),ksq(0:nmax),rho1sq,rho2 
      MODE xi,xsq,fact,sumx,a2,a4,b,c,c2,c3,c4,eth,i10,i11,i30,i31,
     +     i32,hold,m11,m21,m12,m22,left,right,corrus1,corrvs1,mm,delalf
      MODE corrus2,corrvs2,right3,left3,temp

      delalf = pi/ns
c     Note: THROUGHOUT, in original code, had f1=kappa*yp,  f2=-kappa*xd

c     Initialize single layer velocity components us, vs
      do j = 0,n
         u(j) = 0
         v(j) = 0
      enddo

c     COMPUTE VELOENDPTS v(0), v(n) . Note: us(0)=us(n)=0  copied from IntSingleSelfForce
      do k = 1,ns-1      !add kth component, see bubbles/math/Gs.axis
         xi   = ys(k) - y(0)
         xsq  = xs(k)**2
         fact = xs(k)/dsqrt(xsq + xi**2)**3
         v(0)= v(0) + fact*( f1(k)*xs(k)*xi + f2(k)*(2*xi**2+xsq) )
      enddo
      v(0) = v(0)*pi2*delalf 
c     Add  O(h^2) contrib from 1st deriv, skip O(h^4) contrib from 3rd deriv
c     ADD O(h^4) contribution on 12/30/20, see Msendderiv.nb
      xi   = ys(0) - y(0)
      left= 4*pi*f2(0)*xd0/abs(xi)
      temp=xd0/abs(xi)
      left3=pi2*( f2dd0*6*temp + f1d0*6*temp**2 +2*f2(0)*
     +   ( -6*temp**3+xddd0/abs(xi)-3*ydd0*temp/abs(xi) ))
      xi   = ys(ns) - y(0)
      right= 4*pi*f2(ns)*xdn/abs(xi)
      temp=xdn/abs(xi)
      right3=pi2*( f2ddn*6*temp + f1dn*6*temp**2 +2*f2(ns)*
     +   ( -6*temp**3+xdddn/abs(xi)-3*yddn*temp/abs(xi) ))
      v(0)= v(0)-delalf**2/12*(right-left)+delalf**4/720*(right3-left3)
c      print*,y(0),v(0),v(0)+delalf**4/720*(right3-left3)
c Check: good!! using y(0)=pm0.8d0-4 (!!!), print*,ns,v(0), stop
 
      do k = 1,ns-1
         xi   = ys(k) - y(n)
         xsq  = xs(k)**2
         fact = xs(k)/dsqrt(xsq + xi**2)**3
         v(n) = v(n) + fact*(f1(k)*xs(k)*xi+f2(k)*(2*xi**2+xsq))
      enddo
      v(n) = v(n)*pi2*delalf
c     Add  O(h^2) contrib from 1st deriv, skip O(h^4) contrib from 3rd deriv
      xi   = y(n) - ys(0)
      left= 4*pi*f2(0) *xd0   /abs(xi)
      temp=xd0/abs(xi)
      left3=pi2*( f2dd0*6*temp + f1d0*6*temp**2 +2*f2(0)*
     +   ( -6*temp**3+xddd0/abs(xi)-3*ydd0*temp/abs(xi) ))
      xi   = y(n)-ys(ns) 
      right= 4*pi*f2(ns)*xdn/abs(xi)
      temp=xdn/abs(xi)
      right3=pi2*( f2ddn*6*temp + f1dn*6*temp**2 +2*f2(ns)*
     +   ( -6*temp**3+xdddn/abs(xi)-3*yddn*temp/abs(xi) ))
      v(n)= v(n)-delalf**2/12*(right-left)
     + +delalf**4/720*(right3-left3)
c Check: good!! used y(n)=-0.8d0, print*,ns,v(n), stop

c     COMPUTE VELOINT u(j),v(j),j=1,n-1
c     Adding contribution of kth term in sum to jth velocity
c     points on axis (k=0,ns) do not contribute
      do  j = 1,n-1
        do  k = 1,ns-1
           sumx = xs(k) + x(j)
           xi   = ys(k) - y(j)
           ksq(k) = 4*xs(k)*x(j)/( xi**2 + sumx**2)
        enddo
 
        call ellip(ksq(1),F(1),E(1),ns-2)
 
        do k = 1,ns-1
          sumx = xs(k) + x(j)
          xi   = ys(k) - y(j)

          mm =1-ksq(k)
          eth = E(k)/mm

          a2=x(j)**2 + xs(k)**2 + xi**2 
          a4=a2*a2
          b=2*x(j)*xs(k)
          c2=sumx**2 + xi**2
          c=sqrt(c2)
          c3=c*c2
          c4=c2*c2

          i10 = 4*F(k)/c
          i11 = 4/(b*c)*( a2*F(k) - c2*E(k) )
          i30 = 4/c3*eth
          i31 = 4/(b*c3)*(a2*eth-c2*F(k))
          i32 = 4/(b**2*c3)*( a4*eth-2*a2*c2*F(k)+c4*E(k) )

          hold = i11+ (xs(k)**2+x(j)**2)*i31- xs(k)*x(j)*(i30+i32)
          m11 =xs(k)*hold 
          m12 =xs(k)*xi*( xs(k)*i31-x(j)*i30 )
          u(j) = u(j)+ m11*f1(k)+m12*f2(k)
c for f=kappa.n f1=yd*kap, f2=-xd*kap

          m21 =xs(k)*xi*( xs(k)*i30-x(j)*i31 )
          hold = i10+xi**2*i30
          m22 =xs(k)*hold
          v(j) = v(j)+m21*f1(k)+m22*f2(k)
        enddo
      enddo

      do j = 1,n-1
         u(j) = u(j)*delalf
         v(j) = v(j)*delalf
      enddo

      do j=1,n-1
c       Computing and adding endcorrection from 1st derivatives
c        call compendsf(x(j),y(j),ys,ns,xd,xddd,ydd,f2,f2dd,f1d,
c     +     corrus1,corrvs1,corrus2,corrvs2,delalf)
        call compendsf(x(j),y(j),
     +     xd0,xdn,xddd0,xdddn,ys(0),ys(ns),ydd0,yddn,
     +     f1d0,f1dn,f2(0),f2(ns),f2dd0,f2ddn,
     +     corrus1,corrvs1,corrus2,corrvs2,delalf)
        u(j) = u(j) + corrus1 + corrus2
        v(j) = v(j) + corrvs1 + corrvs2
      enddo

      do j=0,n
        u(j)= -u(j)/pi8
        v(j)= -v(j)/pi8
      enddo

      return
      end
