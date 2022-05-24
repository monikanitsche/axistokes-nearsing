#include "flags.h"

      SUBROUTINE compncorrBm(xx,yy,alf,ns,xd0,y0,f1d0,f200,f2d20,
     +   alfb,xb,yb,xdb,ydb,dist,k0,corrbu,corrbv)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Called when point near axis is too close to sheet
c     Output: corrections E[B_Q]_beta to get the exact integral of 
c     THE APPROXIMATION B on curve xa
c     using onesided chebichev points (trapezoid rule in beta)
c     INPUT: xp=xproj,yp=yproj,xd,yd  these define the approximate curve xa
c
c     f1d0,f200,f2d20 : value of f1d,f2,f2d2 at alf=0
c     f10,f11,f12,f13 : coefficients of expansion of f1 about alfb
c     f20,f21,f22,f23 : coefficients of expansion of f2 about alfb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c     INPUT
      integer ns,k0,k0bet
      MODE xx,yy,alf(0:*),xd0,y0,f1d0,f200,f2d20,
     +     alfb,xb,yb,xdb,ydb,dist
c     OUTPUT
      MODE corrbu,corrbv
c     LOCAL
      integer jmax
      MODE bintu,bintv,btrapu,btrapv
      MODE alfend,cc1,cc2,dd1,dd2,
     +     corrbintu,corrbintv,corrbtrapu,corrbtrapv
      MODE bet(0:nmax)
      MODE betb,delalf,xjunk,
     +     dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +     f10,f11,f12,f13,f20,f21,f22,f23,
     +     c1bet,c2bet,c3bet,c4bet,d1bet,d2bet,d3bet,d4bet,
     +     f10bet,f11bet,f12bet,f13bet,f20bet,f21bet,f22bet,f23bet

      call intbm(xx,yy,ns,xd0,y0,f1d0,f200,f2d20,bintu,bintv,
     +           alfend,alfb,xb,yb,xdb,ydb,cc1,cc2,dd1,dd2,jmax)

      call trapbm(xx,yy,alf,jmax,xd0,y0,f1d0,f200,f2d20,
     +           cc1,cc2,dd1,dd2,btrapu,btrapv)


c     sets derivatives at alfb of B(beta) AND B(alf)
      call derivsmatpBbet(
c  input
     +     xx,yy,ns,alfb,xb,yb,xdb,ydb,
     +     xd0,y0,f1d0,f200,f2d20,alfend,cc1,cc2,dd1,dd2,
c output coefficients of xs,ys,f1,f2 (alf)
     +     dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +     f10,f11,f12,f13,f20,f21,f22,f23,
c output coefficients of xs,ys,f1,f2 (bet)
     +     bet,betb,c1bet,c2bet,c3bet,c4bet,d1bet,d2bet,d3bet,d4bet,
     +     f10bet,f11bet,f12bet,f13bet,f20bet,f21bet,f22bet,f23bet)

      k0bet=ns/2  !not needed, only for Gj-Hj
      call correctionmatxnew(
c input
     +        bet,betb,ns,dist,k0bet,
     +        xb,xx,c1bet,c2bet,c3bet,c4bet,dy,d1bet,d2bet,d3bet,d4bet,
     +        f10bet,f11bet,f12bet,f13bet,f20bet,f21bet,f22bet,f23bet,
c output
     +        corrbintu,corrbintv)

c     compute corrections  for B_M(alf) 
c      call correctionmatx(alf,jmax,xx,yy,xp,yp,dist,   !USE WITH V3n
c     +       alfb,corrbtrapu,corrbtrapv,
c     +       xdp,xddp,xdddp,ydp,yddp,ydddp,
c     +       f1p,f1dp,f1ddp,f2p,f2dp,f2d2p)
      call correctionmatxnew(
c input
     +        alf,alfb,jmax,dist,k0,
     +        xb,xx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +        f10,f11,f12,f13,f20,f21,f22,f23,
c output
     +        corrbtrapu,corrbtrapv)

      bintu=bintu+corrbintu
      bintv=bintv+corrbintv
      btrapu=btrapu+corrbtrapu
      btrapv=btrapv+corrbtrapv
      corrbu=bintu-btrapu
      corrbv=bintv-btrapv

      return
      end


      SUBROUTINE derivsmatpBbet(
c  input
     +     xx,yy,ns,alfb,xb,yb,xdb,ydb,
     +     xd0,y0,f1d0,f200,f2d20,alfend,cc1,cc2,dd1,dd2,
c output coefficients of xs,ys,f1,f2 (alf)
     +     dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +     f10,f11,f12,f13,f20,f21,f22,f23,
c output coefficients of xs,ys,f1,f2 (bet)
     +     bet,betb,c1bet,c2bet,c3bet,c4bet,d1bet,d2bet,d3bet,d4bet,
     +     f10bet,f11bet,f12bet,f13bet,f20bet,f21bet,f22bet,f23bet)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,ns
c INPUT      
      MODE xx,yy,alfb,xb,yb,xdb,ydb,
     +     xd0,y0,f1d0,f200,f2d20,alfend,cc1,cc2,dd1,dd2
c OUTPUT      
      MODE dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    f10,f11,f12,f13,f20,f21,f22,f23,
     +    betb,c1bet,c2bet,c3bet,c4bet,d1bet,d2bet,d3bet,d4bet,
     +    f10bet,f11bet,f12bet,f13bet,f20bet,f21bet,f22bet,f23bet
      MODE bet(0:*)
c LOCAL      
      MODE f1,f1d,f1d2,f1d3,f2,f2d,f2d2,f2d3,delbet,
     +     alfd,alfd2,alfd3,alfd4,
     +     xd2,xd3,xd4,yd2,yd3,yd4,
     +     f1dbet,f1d2bet,f1d3bet,f2dbet,f2d2bet,f2d3bet


      delbet=pi/ns
      do k=0,ns
         bet(k)=k*delbet
      enddo

c     THESE ARE DERIVATIVES WRT ALFA
      f1    = f1d0*alfb
      f1d   = f1d0
      f1d2  = 0
      f1d3  = 0

      f2     = f200+f2d20*alfb**2/2
      f2d    = f2d20*alfb
      f2d2   = f2d20
      f2d3   = 0

      xd2 = 6*cc1*alfb   +20*cc2*alfb**3
      xd3 = 6*cc1        +60*cc2*alfb**2
      yd2 = 2*dd1        +12*dd2*alfb**2
      yd3 =               24*dd2*alfb

      dx =xb-xx
      c1 = xdb
      c2 = xd2/2
      c3 = xd3/6
      dy=yb-yy
      d1 = ydb
      d2 = yd2/2
      d3 = yd3/6

      f10=f1
      f11=f1d
      f12=f1d2/2
      f13=0
      f20=f2
      f21=f2d
      f22=f2d2/2
      f23=0

#ifdef THIRD
      xd4=0
      yd4=0
      c4=0
      d4=0
#endif
#ifdef FOURTH
      xd4 = 120*cc2*alfb
      yd4 =  24*dd2
      c4=xd4/24
      d4=yd4/24
#endif


c     NOW FIND DERIVATIVES WRT BETA
      betb  = 2*acos(1-alfb/alfend)      !projection xp=x(betb)
      alfd  =  sin(betb/2)*alfend/2   !dalfdbet at beta=betb
      alfd2 =  cos(betb/2)*alfend/4
      alfd3 = -sin(betb/2)*alfend/8
      alfd4 = -cos(betb/2)*alfend/16


c derivatives of x(alf(beta)) wrt beta
c derivatives of f1(alf(beta)) wrt beta
      c1bet  = xdb*alfd
      c2bet  = (xd2*alfd**2+xdb*alfd2)/2
      c3bet  = (xd3*alfd**3 + 3*xd2*alfd*alfd2 + xdb*alfd3)/6

      d1bet  = ydb*alfd
      d2bet  = (yd2*alfd**2+ydb*alfd2)/2
      d3bet  = (yd3*alfd**3 + 3*yd2*alfd*alfd2 + ydb*alfd3)/6

      f1dbet  = f1d*alfd
      f1d2bet = f1d2*alfd**2+f1d*alfd2
      
      f2dbet  = f2d*alfd
      f2d2bet = f2d2*alfd**2+f2d*alfd2

c product fbet=f*alfdbet: use alfdbet=alfdp+alfddp*bet+alfdddp/2*bet^2
      f10bet  = f1*alfd
      f11bet  = f1*alfd2 +   f1dbet*alfd
      f12bet  = (f1*alfd3 + 2*f1dbet*alfd2 + f1d2bet*alfd)/2

      f20bet  = f2*alfd
      f21bet  = f2*alfd2  +   f2dbet*alfd
      f22bet  = (f2*alfd3 + 2*f2dbet*alfd2 + f2d2bet*alfd)/2

#ifdef THIRD
      c4bet=0
      d4bet=0
      f13bet=0
      f23bet=0
#endif
#ifdef FOURTH
      c4bet = (xd4*alfd**4 + 6*alfd**2*alfd2*xd3+
     +         4*xd2*alfd*alfd3 + 3*xd2*alfd2**2 +xdb*alfd4)/24
      d4bet = (yd4*alfd**4 + 6*alfd**2*alfd2*yd3+
     +         4*yd2*alfd*alfd3 + 3*yd2*alfd2**2 +ydb*alfd4)/24

      f1d3bet = f1d3*alfd**3 + 3*f1d2*alfd*alfd2 + f1d*alfd3
      f2d3bet = f2d3*alfd**3 + 3*f2d2*alfd*alfd2 + f2d*alfd3

      f13bet = (f1d3*alfd**3 + 3*f1d2*alfd*alfd2 + f1d*alfd3)/6
      f23bet = (f2d3*alfd**3 + 3*f2d2*alfd*alfd2 + f2d*alfd3)/6
#endif

      return
      end



      SUBROUTINE compBm(xx,yy,alf,ns,xd0,y0,f1d0,f20,f2d20,
     +        corrbu,corrbv)
#include "params.h"
c INPUT
      integer ns
      MODE xx,yy,alf(0:*),xd0,y0,f1d0,f20,f2d20
c OUTPUT
      MODE corrbu,corrbv
c LOCAL
      integer jmax,k0
      MODE bintu,bintv,btrapu,btrapv
      MODE alfb,xb,yb,xdb,ydb
      MODE alfend,cc1,cc2,dd1,dd2
      MODE xk,xik,uxk,vxk,uyk,vyk

      call findalfp(alfb,xx,yy,xb,yb,xdb,ydb)
      call intbm(xx,yy,ns,xd0,y0,f1d0,f20,f2d20,bintu,bintv,
     +           alfend,alfb,xb,yb,xdb,ydb,cc1,cc2,dd1,dd2,jmax)

      call trapbm(xx,yy,alf,jmax,xd0,y0,f1d0,f20,f2d20,
     +     cc1,cc2,dd1,dd2,btrapu,btrapv)

      corrbu=bintu-btrapu
      corrbv=bintv-btrapv
      return
      end

c        write(8,*)xx,yy,u(j),v(j),corrmu,corrmv,
c     +   (bintmu+corrbintmu)-(btrapmu+corrbtrapmu),
c     +   (bintmv+corrbintmv)-(btrapmv+corrbtrapmv),zero,zero
