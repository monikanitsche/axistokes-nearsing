#include "flags.h"
      SUBROUTINE veloforceSonP(xpt,ypt,alf,xs,ys,ns,f1,f2,
     +       xd0,xddd0,ydd0,f1d0,f2dd0,
     +       xdn,xdddn,yddn,f1dn,f2ddn,
     +       upt,vpt)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc     Date: 1/28/18
c
c     Computes velocity u,v at sheet x,y,n induced by sheet xs,ys,ns,f1,f2
c     ASSUME: 1) xs(0),xs(n) on axis. 
c             2) x(0),x(n) on axis
c             3) x,y and xs ys are two different sheets! so no local corrections
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns
      MODE xpt,ypt,alf(0:*),xs(0:*),ys(0:*),f1(0:*),f2(0:*),
     +       xd0,xddd0,ydd0,f1d0,f2dd0,xdn,xdddn,yddn,f1dn,f2ddn
c OUTPUT
      MODE upt,vpt
c LOCAL
      integer k,k0
      MODE F(0:nmax),E(0:nmax),ksq(0:nmax)
      MODE xi,sumx,a2,a4,b,c,c2,c3,c4,eth,i10,i11,i30,i31,i32,
     +     hold,m11,m21,m12,m22,mm,delalf
      MODE alfb,corrmu,corrmv,corrbu,corrbv,dist,xb,yb,xdb,ydb
      MODE corrus1,corrvs1,corrus2,corrvs2
      logical yes,nearaxis

      delalf = pi/ns

c     Initialize single layer velocity components us, vs
      upt = 0
      vpt = 0

c     Adding contribution of kth term in sum to jth velocity
c     points on axis (k=0,ns) do not contribute
      do  k = 1,ns-1
           sumx = xs(k) + xpt
           xi   = ys(k) - ypt
           ksq(k) = 4*xs(k)*xpt/( xi**2 + sumx**2)
      enddo
 
      call ellip(ksq(1),F(1),E(1),ns-2)
 
      do k = 1,ns-1
          sumx = xs(k) + xpt
          xi   = ys(k) - ypt

          mm =1-ksq(k)
          eth = E(k)/mm

          a2=xpt**2 + xs(k)**2 + xi**2 
          a4=a2*a2
          b=2*xpt*xs(k)
          c2=sumx**2 + xi**2
          c=sqrt(c2)
          c3=c*c2
          c4=c2*c2

          i10 = 4*F(k)/c
          i11 = 4/(b*c)*( a2*F(k) - c2*E(k) )
          i30 = 4/c3*eth
          i31 = 4/(b*c3)*(a2*eth-c2*F(k))
          i32 = 4/(b**2*c3)*( a4*eth-2*a2*c2*F(k)+c4*E(k) )

          hold = i11+ (xs(k)**2+xpt**2)*i31- xs(k)*xpt*(i30+i32)
          m11 =xs(k)*hold 
          m12 =xs(k)*xi*( xs(k)*i31-xpt*i30 )
          upt = upt+ m11*f1(k)+m12*f2(k)
c for f=kappa.n f1=yd*kap, f2=-xd*kap

          m21 =xs(k)*xi*( xs(k)*i30-xpt*i31 )
          hold = i10+xi**2*i30
          m22 =xs(k)*hold
          vpt = vpt+m21*f1(k)+m22*f2(k)
      enddo

      upt = upt*delalf
      vpt = vpt*delalf

c     Computing and adding endcorrection from 1st and 3rd derivatives
      call compendsf(xpt,ypt,
     +     xd0,xdn,xddd0,xdddn,ys(0),ys(ns),ydd0,yddn,
     +     f1d0,f1dn,f2(0),f2(ns),f2dd0,f2ddn,
     +     corrus1,corrvs1,corrus2,corrvs2,delalf)
      upt = upt + corrus1 + corrus2
      vpt = vpt + corrvs1 + corrvs2

      upt= -upt/pi8
      vpt= -vpt/pi8

c add corrections if near interface or near axis 
      yes=.false.
      nearaxis=.false.
      corrmu=0
      corrmv=0
      corrbu=0
      corrbv=0

      if ( (xpt**2+(ypt-ys(0))**2).lt.0.2375**2) nearaxis=.true. 
      call x2cl2circ(xpt,ypt,ns,yes,k0,alfb,xb,yb,dist)
      if (yes) then   
         call correctmnearsheet(xpt,ypt,alf,ns,k0,alfb,xb,yb,dist,
     +             corrmu,corrmv,xdb,ydb)
         if (nearaxis) 
     +      call compncorrBm(xpt,ypt,alf,ns,xd0,ys(0),f1d0,f2(0),
     +             f2dd0,alfb,xb,yb,xdb,ydb,dist,k0,
     +             corrbu,corrbv)
      else
         if (nearaxis) 
     +      call compBm(xpt,ypt,alf,ns,xd0,ys(0),f1d0,f2(0),f2dd0,
     +             corrbu,corrbv)
      endif

      upt=upt+corrmu+corrbu
      vpt=vpt+corrmv+corrbv

      return
      end
