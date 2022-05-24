#include "flags.h"
      SUBROUTINE velodSonP(xx,yy,alf,xs,ys,ns,u,v,xd,yd,
     +      xddd0,ydd0,ud0,vdd0,xdddn,yddn,udn,vddn,
     +      uu,vv)
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
c INPUT
      integer ns
      MODE xx,yy
      MODE alf(0:*),xs(0:*),ys(0:*),u(0:*),v(0:*),xd(0:*),yd(0:*)
      MODE xddd0,ydd0,ud0,vdd0,xdddn,yddn,udn,vddn
c OUTPUT
      MODE uu,vv
c LOCAL
      integer k,k0
      MODE ksq(0:nmax),F(0:nmax),E(0:nmax),delalf
      MODE xi,sumx,a2,a4,a6,b,b2,b3,c,c2,c3,c4,c5,c6,
     +       eth,efh,ac,cb,sing,i50,i51,i52,i53,
     +       q11,q12,q21,q22,m,efhp,eminusone
      MODE corrud1,corrvd1,corrud2,corrvd2,oneminksq(0:nmax)
      MODE alfb,xb,yb,dist,corrqu,corrqv,corrbu,corrbv,xjunk
      MODE ud(0:nmax),ud2(0:nmax),vd(0:nmax),vd2(0:nmax),xdb,ydb,
     +       restu,restv,hk,rest1
      logical yes,nearaxis,largecorr

      delalf = pi/ns

      uu=0
      vv=0
 
c      print*,xx,yy,'in velod 0'
c     Calculate velo at xx due to xs(k) ring (k=0,ns dont contribute)
      do  k   = 1,ns-1
            sumx = xs(k) + xx
            xi   = ys(k) - yy
 
            ksq(k) = 4*xs(k)*xx/( xi**2 + sumx**2)
            oneminksq(k) = (xi**2+(xs(k)-xx)**2)/( xi**2 + sumx**2)
      enddo
 
c     call ellip(ksq(1),F(1),E(1),ns-2)
      call ellip1mink(ksq(1),oneminksq(1),F(1),E(1),ns-2)
 
      do k = 1,ns-1
            sumx = xs(k) + xx
            xi   = ys(k) - yy

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

            sing = 8/(3*c5*m**2)
            i50= 4/c5*(efhp) 
            i51= 4/c5*(cb*(ac*efhp-eth)+4/(3*m*(1-m)) ) 
            i52=4/c5*(cb**2*(ac**2*efhp-2*ac*eth+F(k))+8/(3*m*(1-m)**2))
            i53= 4/c5*( cb**3*(ac**3*efhp-3*ac**2*eth+3*ac*F(k)-E(k)) 
     +            +4*(3+m**2)/(3*m*(1-m)**3) )

            q11=-6*xs(k)*( xs(k)**3*i51 - xs(k)**2*xx*(i50+2*i52) +
     +                    xs(k)*xx**2*(i53+2*i51) - xx**3*i52
     +                    + (xs(k)-xx)**3*sing)
                
            q12=-6*xs(k)*xi*((xs(k)**2+xx**2)*i51-xs(k)*xx*(i50+i52)
     +                    + (xs(k)-xx)**2*sing)
            q21= q12
            q22= -6*xs(k)*xi**2*( xs(k)*i51 -xx*i50+(xs(k)-xx)*sing)

            uu = uu + (q11*yd(k)-q12*xd(k)) * u(k)
     +                      + (q21*yd(k)-q22*xd(k)) * v(k)
c            print*,k,q11,q12,q21,q22,yd(k),xd(k),u(k),v(k)

c TEST: if (j.eq.960) then, if (k.eq.384) then, print i50,i51,e,f,eth,efh,efhp, 
c in double and quad before and after replacing E-1. reduced 1e-8 to 1e-12

            q11=-6*xs(k)*xi*( xx**2*i52 +xs(k)**2*i50-2*xs(k)*xx*i51
     +                    + (xx-xs(k))**2*sing)
            q12=-6*xs(k)*xi**2*( xs(k)*i50 -xx*i51 +(xs(k)-xx)*sing)
            q21= q12
            q22= -6*xs(k)*xi**3*(i50+sing)

            vv = vv + (q11*yd(k)-q12*xd(k)) * u(k)
     +                      + (q21*yd(k)-q22*xd(k)) * v(k)
      enddo
 
      uu = uu*delalf
      vv = vv*delalf
c      print*,uu,vv,'in velod 1'
 
      call compenddv(xx,yy,delalf,
     +      xd(0),xddd0,ys(0),ydd0,ud0,v(0),vdd0,
     +      xd(ns),xdddn,ys(ns),yddn,udn,v(ns),vddn,
     +      corrud1,corrvd1,corrud2,corrvd2)
      uu = uu + corrud1  + corrud2
      vv = vv + corrvd1  + corrvd2
c      print*,uu,vv,'in velod 2'

      uu = uu/pi8
      vv = vv/pi8
c      print*,uu,vv,'in velod 3'

      yes=.false.
      nearaxis=.false.
      corrbu=0
      corrbv=0
      corrqu=0
      corrqv=0
     
      if ( (xx**2+(yy-ys(0))**2).lt.0.2375**2) nearaxis=.true. !.2375=.95/4
      call x2cl2circ(xx,yy,ns,yes,k0,alfb,xb,yb,dist)
      if (yes) then   
         xjunk=0
         largecorr=.false.
         call correctqnearsheetnew(alf,ns,xx,yy,alfb,xb,yb,dist,
     +             corrqu,corrqv,xdb,ydb,largecorr,k0,restu,restv,
     +             hk,rest1,xjunk,xjunk,xjunk,xjunk,xjunk,xjunk)
         if (nearaxis) 
     +        call compncorrBq(xx,yy,alf,ns,xd(0),ys(0),ud0,v(0),
     +                vdd0,alfb,xb,yb,xdb,ydb,dist,k0,
     +                corrbu,corrbv)
      else
         if (nearaxis) 
     +      call compBq(xx,yy,alf,ns,
     +                  xd(0),ys(0),ud0,v(0),vdd0,
     +                  corrbu,corrbv)
      endif
c      print*,corrbu,corrbv,' in velod 4'

      uu=uu+corrqu+corrbu
      vv=vv+corrqv+corrbv

      return
      end
