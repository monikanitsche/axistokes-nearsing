#include "flags.h"
      SUBROUTINE velodSonO(xs,ys,ns,u,v,x,y,n,udb,vdb,
     +   xd,yd,xddd0,ydd0,ud0,vdd0,xdddn,yddn,udn,vddn)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes the double layer component at x,y induced by xs,ys,u,v
c           ud=int Q11*u+Q12*v dalf 
c
c
c     common all derivatives of x,y, precomputed earlier
c     input: xs,ys,u,v,ns  values on sheet, size 0:ns
c            xd,yd,xddd0,ydd0,ud0,vdd0,xdddn,yddn,udn,vddn :
c                  derivatives of xs,ys,u,v on sheet, at endpoints j0,ns
c            x,y,n target point, size 0:n, 0 on axis
c     output: udb,vdb double layer velo at target points 0:n
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,k,ns,n
      MODE xs(0:*),ys(0:*),u(0:*),v(0:*),x(0:*),y(0:*),udb(0:*),
     +  vdb(0:*),xd(0:*),yd(0:*),xddd0,ydd0,ud0,vdd0,xdddn,yddn,udn,vddn
      MODE ksq(0:nmax),F(0:nmax),E(0:nmax),delalf
      MODE xi,xsq,fact,sumx,a2,a4,a6,b,b2,b3,c,c2,c3,c4,c5,c6,
     +       eth,efh,ac,cb,k2,sing,i50,i51,i52,i53,
c     +       q11,q12,q21,q22,
     +       q11,q12,q21,q22,m,left,right,efhp,right3,left3
      MODE corrud1,corrvd1,corrud2,corrvd2,oneminksq(0:nmax)
      MODE eminusone

      delalf = pi/ns

      do j = 0,n
         udb(j) = 0
         vdb(j) = 0
      enddo

c     SPECIAL CASE 1 j = 0 (were x,y on axis)
      do k = 1,ns-1
         xi   = ys(k) - y(0)
         xsq  = xs(k)**2
         fact = xs(k)/sqrt(xsq + xi**2)**3
         fact = 6*xi*fact/(xsq + xi**2)
         vdb(0)= vdb(0)-fact*(u(k)*xs(k)+v(k)*xi)*(xs(k)*yd(k)-xi*xd(k))
      enddo
      vdb(0) = vdb(0)*pi2*delalf

c     Add O(h^2) contrib from 1st deriv at endpoints of S1
      xi   = ys(0) - y(0)
      left = 12*pi*v(0)*xd(0)**2/(xi*abs(xi))
      left3 = pi2*( 6*v(0)*xd(0)*(-15*xd(0)**3+4*xddd0*xi**2
     +                            -12*xi*xd(0)*ydd0)/(xi**3*abs(xi))
     +       +36*ud0*xd(0)**3/(xi**2*abs(xi))
     +       +18*vdd0*xd(0)**2/(xi*abs(xi)) )
      xi   = ys(ns) - y(0)
      right= 12*pi*v(ns)*xd(ns)**2/(xi*abs(xi))
      right3 = pi2*( 6*v(ns)*xd(ns)*(-15*xd(ns)**3+4*xdddn*xi**2
     +                         -12*xi*xd(ns)*yddn)/(xi**3*abs(xi))
     +       +36*udn*xd(ns)**3/(xi**2*abs(xi))
     +       +18*vddn*xd(ns)**2/(xi*abs(xi)) )
      vdb(0)= vdb(0) -delalf**2/12*(right-left)
     + +delalf**4/720*(right3-left3)

c     SPECIAL CASE 2 j = n (were x,y on axis)
      do k = 1,ns-1
         xi   = ys(k) - y(n)
         xsq  = xs(k)**2
         fact = xs(k)/sqrt(xsq + xi**2)**3
         fact = 6*xi*fact/(xsq + xi**2)
         vdb(n)= vdb(n)-fact*(u(k)*xs(k)+v(k)*xi)*(xs(k)*yd(k)-xi*xd(k))
      enddo
      vdb(n) = vdb(n)*pi2*delalf

c     Add O(h^2) contrib from 1st deriv at endpoints of S1
      xi   = ys(0)-y(n) 
      left=12*pi*v(0)*xd(0)**2/(xi*abs(xi))
      left3 = pi2*( 6*v(0)*xd(0)*(-15*xd(0)**3+4*xddd0*xi**2
     +                            -12*xi*xd(0)*ydd0)/(xi**3*abs(xi))
     +       +36*ud0*xd(0)**3/(xi**2*abs(xi))
     +       +18*vdd0*xd(0)**2/(xi*abs(xi)) )
      xi   = ys(ns) - y(n)
      right=12*pi*v(ns)*xd(ns)**2/(xi*abs(xi))
      right3 = pi2*( 6*v(ns)*xd(ns)*(-15*xd(ns)**3+4*xdddn*xi**2
     +                         -12*xi*xd(ns)*yddn)/(xi**3*abs(xi))
     +       +36*udn*xd(ns)**3/(xi**2*abs(xi))
     +       +18*vddn*xd(ns)**2/(xi*abs(xi)) )
      vdb(n)= vdb(n) -delalf**2/12*(right-left)
     + +delalf**4/720*(right3-left3)

c     REMAINING CASES  j = 1,n-1
 
      do j = 1,n-1
c         if (j.eq.8) print*,x(j),y(j),'in SonO 0'
 
c        Calculate velo at x(j) due to xs(k) ring (k=0,ns dont contribute)
         do  k   = 1,ns-1
            sumx = xs(k) + x(j)
            xi   = ys(k) - y(j)
 
            ksq(k) = 4*xs(k)*x(j)/( xi**2 + sumx**2)
            oneminksq(k) = (xi**2+(xs(k)-x(j))**2)/( xi**2 + sumx**2)
         enddo
 
c USED FOR TESTING ROUTINES, COMPARISON WITH MATLAB
c         ksq(1)=1-0.0001q0
c         ksq(2)=1-0.00001q0
c         ksq(3)=1-0.000001q0
c         ksq(4)=1-0.0000001q0
c         ksq(5)=1-0.00000001q0
c         ksq(6)=1-0.0000000000001q0
c         call ellip(ksq(1),F(1),E(1),ns-2)
c         oneminksq(1)=0.0001q0
c         oneminksq(2)=0.00001q0
c         oneminksq(3)=0.000001q0
c         oneminksq(4)=0.0000001q0
c         oneminksq(5)=0.00000001q0
c         oneminksq(6)=0.0000000000001q0

c         call ellip(ksq(1),F(1),E(1),ns-2)
         call ellip1mink(ksq(1),oneminksq(1),F(1),E(1),ns-2)
 
         do k = 1,ns-1
            sumx = xs(k) + x(j)
            xi   = ys(k) - y(j)

c            m =1-ksq(k)
            m =oneminksq(k)
            eth = E(k)/m
            efh = (2*(1+m)/m*E(k)-F(k))/(3*m)
            efhp = (2*(1+m)/m*E(k)-F(k))/(3*m) -2/(3*m**2)
            eminusone = E(k)-1

c coefficients ec1,ec2,ec3,ed1,ed2,ed3 set in params.h
            if (m.lt.1.e-4) eminusone=m*(ec1+ec2*m+ec3*m**2)
     +                  -m*log(m)*(ed1+ed2*m+ed3*m**2)
            efhp=(2*(1+m)/m*eminusone-F(k)+2)/(3*m)  !efhp

            a2=x(j)**2 + xs(k)**2 + xi**2 
            a4=a2*a2
            a6=a4*a2
            b=2*x(j)*xs(k)
            b2=b**2
            b3=b**3
            c2=sumx**2 + xi**2
            c=sqrt(c2)
            c3=c*c2
            c4=c2*c2
            c5=c4*c
            c6=c3*c3
 
            ac=(1+m)/2
            cb=2/(1-m)
            k2=ksq(k)

            sing = 8/(3*c5*m**2)
            i50= 4/c5*(efhp) 
            i51= 4/c5*(cb*(ac*efhp-eth)+4/(3*m*(1-m)) ) 
            i52=4/c5*(cb**2*(ac**2*efhp-2*ac*eth+F(k))+8/(3*m*(1-m)**2))
            i53= 4/c5*( cb**3*(ac**3*efhp-3*ac**2*eth+3*ac*F(k)-E(k)) 
     +            +4*(3+m**2)/(3*m*(1-m)**3) )

            q11=-6*xs(k)*( xs(k)**3*i51 - xs(k)**2*x(j)*(i50+2*i52) +
     +                    xs(k)*x(j)**2*(i53+2*i51) - x(j)**3*i52
     +                    + (xs(k)-x(j))**3*sing)
                
            q12=-6*xs(k)*xi*((xs(k)**2+x(j)**2)*i51-xs(k)*x(j)*(i50+i52)
     +                    + (xs(k)-x(j))**2*sing)
            q21= q12
            q22= -6*xs(k)*xi**2*( xs(k)*i51 -x(j)*i50+(xs(k)-x(j))*sing)

c            if (j.eq.8)
c     + print*,k,q11,q12,q21,q22,yd(k),xd(k),u(k),v(k)
            udb(j) = udb(j) + (q11*yd(k)-q12*xd(k)) * u(k)
     +                      + (q21*yd(k)-q22*xd(k)) * v(k)

c TEST: if (j.eq.960) then, if (k.eq.384) then, print i50,i51,e,f,eth,efh,efhp, 
c in double and quad before and after replacing E-1. reduced 1e-8 to 1e-12

            q11=-6*xs(k)*xi*( x(j)**2*i52 +xs(k)**2*i50-2*xs(k)*x(j)*i51
     +                    + (x(j)-xs(k))**2*sing)
            q12=-6*xs(k)*xi**2*( xs(k)*i50 -x(j)*i51 +(xs(k)-x(j))*sing)
            q21= q12
            q22= -6*xs(k)*xi**3*(i50+sing)

            vdb(j) = vdb(j) + (q11*yd(k)-q12*xd(k)) * u(k)
     +                      + (q21*yd(k)-q22*xd(k)) * v(k)
         enddo
      enddo
 
      do j = 1,n-1
         udb(j) = udb(j)*delalf
         vdb(j) = vdb(j)*delalf
c         if (j.eq.8) print*,udb(j),vdb(j),'in SonO 1'
      enddo
 
      do j=1,n-1
         call compenddv(x(j),y(j),delalf,
     +      xd(0),xddd0,ys(0),ydd0,ud0,v(0),vdd0,
     +      xd(ns),xdddn,ys(ns),yddn,udn,v(ns),vddn,
     +      corrud1,corrvd1,corrud2,corrvd2)
        udb(j) = udb(j) + corrud1  + corrud2
        vdb(j) = vdb(j) + corrvd1  + corrvd2
c         if (j.eq.8) print*,udb(j),vdb(j),'in SonO 2'
      enddo
 
      do j = 0,n
         udb(j) = udb(j)/pi8
         vdb(j) = vdb(j)/pi8
c         if (j.eq.8) print*,udb(j),vdb(j),'in SonO 3'
      enddo

      return
      end


      SUBROUTINE velodSonPnok(xs,ys,ns,u,v,k0,xx,yy,upt,vpt,
     +   xd,yd,xddd0,ydd0,ud0,vdd0,xdddn,yddn,udn,vddn,gkself)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes the double layer integral induced by xs,ys,u,v,ns at
c     one point P(xx,yy), WITHOUT contribution of the kth component
c     ASSUME k0.ne.0 or ns
c
c     output: upt,vpt double layer velo at target points 0:n
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,k,ns,k0
      MODE xs(0:*),ys(0:*),u(0:*),v(0:*),xx,yy,upt,vpt,
     +  xd(0:*),yd(0:*),xddd0,ydd0,ud0,vdd0,xdddn,yddn,udn,vddn
      MODE ksq(0:nmax),F(0:nmax),E(0:nmax),delalf
      MODE xi,xsq,fact,sumx,a2,a4,a6,b,b2,b3,c,c2,c3,c4,c5,c6,
     +       eth,efh,ac,cb,k2,sing,i50,i51,i52,i53,
c     +       q11,q12,q21,q22,m,left,right,efhp,right3,left3
     +       m,left,right,efhp,right3,left3,
     +       q11u(0:nmax),q12u(0:nmax),q21u(0:nmax),q22u(0:nmax),
     +       q11v(0:nmax),q12v(0:nmax),q21v(0:nmax),q22v(0:nmax)
      MODE corrud1,corrvd1,corrud2,corrvd2,gkself

      delalf = pi/ns

      upt = 0
      vpt = 0

c     ASSUME x(j),y(j) NOT ON AXIS
 
c     xx=x(j)
c     yy=y(j)
c      do j = 1,n-1
 
c        Calculate velo at x(j) due to xs(k) ring (k=0,ns dont contribute)
         do  k   = 1,ns-1
            sumx = xs(k) + xx
            xi   = ys(k) - yy
 
            ksq(k) = 4*xs(k)*xx/( xi**2 + sumx**2)
         enddo
 
         call ellip(ksq(1),F(1),E(1),ns-2)
 
         do k = 1,ns-1
            sumx = xs(k) + xx
            xi   = ys(k) - yy

            m =1-ksq(k)
            eth = E(k)/m
            efh = (2*(1+m)/m*E(k)-F(k))/(3*m)
            efhp = (2*(1+m)/m*E(k)-F(k))/(3*m) -2/(3*m**2)

            a2=xx**2 + xs(k)**2 + xi**2 
            a4=a2*a2
            a6=a4*a2
            b=2*xx*xs(k)
            b2=b**2
            b3=b**3
            c2=sumx**2 + xi**2
            c=sqrt(c2)
            c3=c*c2
            c4=c2*c2
            c5=c4*c
            c6=c3*c3
 
            ac=(1+m)/2
            cb=2/(1-m)
            k2=ksq(k)

            sing = 8/(3*c5*m**2)
            i50= 4/c5*(efhp) 
            i51= 4/c5*(cb*(ac*efhp-eth)+4/(3*m*(1-m)) ) 
            i52=4/c5*(cb**2*(ac**2*efhp-2*ac*eth+F(k))+8/(3*m*(1-m)**2))
            i53= 4/c5*( cb**3*(ac**3*efhp-3*ac**2*eth+3*ac*F(k)-E(k)) 
     +            +4*(3+m**2)/(3*m*(1-m)**3) )

            q11u(k)=-6*xs(k)*( xs(k)**3*i51 - xs(k)**2*xx*(i50+2*i52) +
     +                    xs(k)*xx**2*(i53+2*i51) - xx**3*i52
     +                    + (xs(k)-xx)**3*sing)
                
            q12u(k)=-6*xs(k)*xi*((xs(k)**2+xx**2)*i51-xs(k)*xx*(i50+i52)
     +                    + (xs(k)-xx)**2*sing)
            q21u(k)= q12u(k)
            q22u(k)= -6*xs(k)*xi**2*( xs(k)*i51 -xx*i50+(xs(k)-xx)*sing)

c            if (k.ne.k0) then
c            upt = upt + (q11*yd(k)-q12*xd(k)) * u(k)
c     +                      + (q21*yd(k)-q22*xd(k)) * v(k)
c            else
c              gkself=q11*yd(k)*u(k)
c              gkself=-q12*xd(k)*u(k)
c              gkself=q21*yd(k)*v(k)
c              gkself=-q22*xd(k)*v(k)
c            endif
c
            q11v(k)=-6*xs(k)*xi*( xx**2*i52 +xs(k)**2*i50-2*xs(k)*xx*i51
     +                    + (xx-xs(k))**2*sing)
            q12v(k)=-6*xs(k)*xi**2*( xs(k)*i50 -xx*i51 +(xs(k)-xx)*sing)
            q21v(k)= q12v(k)
            q22v(k)= -6*xs(k)*xi**3*(i50+sing)
c
c            if (k.ne.k0) then
c            vpt = vpt + (q11*yd(k)-q12*xd(k)) * u(k)
c     +                      + (q21*yd(k)-q22*xd(k)) * v(k)
c            else
c              gkself=q11*yd(k)*u(k)
c              gkself=-q12*xd(k)*u(k)
c              gkself=q21*yd(k)*v(k)
c              gkself=-q22*xd(k)*v(k)
c            endif
      enddo
      do k = 1,k0-1
            upt = upt + (q11u(k)*yd(k)-q12u(k)*xd(k)) * u(k)
     +                + (q21u(k)*yd(k)-q22u(k)*xd(k)) * v(k)
            vpt = vpt + (q11v(k)*yd(k)-q12v(k)*xd(k)) * u(k)
     +                + (q21v(k)*yd(k)-q22v(k)*xd(k)) * v(k)
      enddo
      do k = ns-1,k0+1,-1
            upt = upt + (q11u(k)*yd(k)-q12u(k)*xd(k)) * u(k)
     +                      + (q21u(k)*yd(k)-q22u(k)*xd(k)) * v(k)
            vpt = vpt + (q11v(k)*yd(k)-q12v(k)*xd(k)) * u(k)
     +                      + (q21v(k)*yd(k)-q22v(k)*xd(k)) * v(k)
      enddo
      gkself= q11u(k0)*yd(k)*u(k)
c      gkself= q11v(k0)*yd(k)*u(k)
c     gkself=-q12v(k0)*xd(k)*u(k)
c     gkself= q21v(k0)*yd(k)*v(k)
c      gkself=-q22v(k0)*xd(k)*v(k)
 
c      do j = 1,n-1
         upt = upt*delalf
         vpt = vpt*delalf
c      enddo
 
c      do j=1,n-1
         call compenddv(xx,yy,delalf,
     +      xd(0),xddd0,ys(0),ydd0,ud0,v(0),vdd0,
     +      xd(ns),xdddn,ys(ns),yddn,udn,v(ns),vddn,
     +      corrud1,corrvd1,corrud2,corrvd2)
        upt = upt + corrud1  + corrud2
        vpt = vpt + corrvd1  + corrvd2
c      enddo
 
         upt = upt/pi8
         vpt = vpt/pi8

      return
      end
