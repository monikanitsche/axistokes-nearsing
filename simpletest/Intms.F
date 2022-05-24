c Contains
c SUBROUTINE IntQs: comp ud,vd at x,y due to xs,ys and then corrects
c SUBROUTINE findclosest: (not used) generic routine to find closest pt
c                        on interface
c SUBROUTINE initvars: init vars to .false., =0, and vals of xmax,width 
c SUBROUTINE fixroundoff: routine to replace hj-gj by rest if largecorr
c                        (not used, not needed here)

#include "flags.h"
      SUBROUTINE IntMs(xs,ys,ns,f1,f2,x,y,n,u,v)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Computes single layer integral usg,vsg 
c    induced by sheet xs,ys,ns with density u,v (SHEET 1)
c    at sheet x,y,n (SHEET 2). Note:
c         x: radial 
c         y: axial
c    Assumes both sheets have j=0,n on axis
c    Includes corrections 
c         near interface (xs,ys), 
c         near axis, 
c         for small r 
c    Prints results in files set in drivergr
c
c    VARIABLES
c     global derivs1.h: derivatives of SHEET 1
c     input:  xs,ys,ns (position of SHEET 1)
c             alf    (parametrization of SHEET 1)
c             u,v    (density of SHEET 1)
c             x,y,n  target points (SHEET 2) must be different from SHEET 1 
c     output: u,v  (sgle layer integral at x,y on SHEET 2 induced by SHEET 1)
c
c     Variables:
c       subscript s: sheet
c       subscript a: approximating sheet (for B)
c       subscript p or b: projection 
c       subscript d: derivatives
c
c       xd,xdd,xddd, yd,ydd,yddd: derivatives of xs,ys at alfj j=0,ns
c       xp,yp,alfb,xdp,xddp,xdddp,etc : projection: values at alfb
c       x0,y0,xd0,yd0,xdd0,ydd0,etc : values at axis alf=0
c       xn,yn,xdn,ydn,xddn,yddn,etc : values at axis alfn=pi
c
c       xsa,ysa,alf : approximating sheet to get B
c
c OUTLINE
c compute T[G_j] (velod udb(j)) j=0,...,n
c correct endpoints j=0 and j=n
c determine if nearsheet (yes), nearaxis, smallr
c if yes
c     compute correction E[H]
c     determine if large corr (not needed here)
c        if largecorr, add rest 
c     if nearaxis, add E[B]
c else
c     if nearaxis, add E[B]
c endif
c if smallr, use extrapolation
c
c print: T[G], correction E[G], correction E[B]. Fullcorr T[G]+E[G]+E[B]
c if smallr, then print extrapolated entry in last column
c plot data using ploterrgriddb with ploterrgridsuball3.m
c (properly reads in data in 1600q4 which is in a different format)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
#include "derivs1.h"    
C     INPUT
      integer ns,n
      MODE xs(0:*),ys(0:*),f1(0:*),f2(0:*),x(0:*),y(0:*)
C     OUTPUT
      MODE u(0:*),v(0:*)
C     LOCAL
      integer j,k0
      MODE xx,yy,alfb,xb,yb,dist,corrmu,corrmv,corrbu,corrbv,
     +       uextrap,vextrap
      MODE ud3b,zd4b,xjunk,zero,width,corrsize,xmax
      MODE f1d(0:nmax),f1d2(0:nmax),f2d(0:nmax),f2d2(0:nmax),xdb,ydb
      logical yes,nearaxis,smallr

c compute derivatives of "densities", needed below
      call compderiv4(f1,f1d,f1d2,f2,f2d,f2d2,ns)

c compute velocity at (x,y) induced by (xs,ys) using T4[G] or T6[G]     
      call veloforceSonO(x,y,n,xs,ys,ns,f1,f2,
     +   xd(0), xddd(0), ydd(0), f1d(0), f2d2(0),
     +   xd(ns),xddd(ns),ydd(ns),f1d(ns),f2d2(ns),
     +   u,v)

c compute correction for target points on axis (j=0)
      alfb=0
      yy=y(0)
      dist=abs(yy-ys(0))
      corrmv=0
      if (dist.lt.0.4) then
         zd4b=2*(ys(2)-4*ys(1)+3*ys(0))/delalf**4  !z even
         call correctionmendsh(alf,ns,yy,ys(0),dist,alfb,corrmv,
     +     xd(0),xddd(0),ydd(0),zd4b,f1d(0),f2(0),f2d2(0))
      endif
      zero=0
      write(8,1000)x(j),y(j), u(j),v(j),
     +        zero,corrmv, zero,zero, zero,zero

c compute corrections for target points j=1..n-1 not on axis
      do j = 1,n-1  
         xx=x(j)
         yy=y(j)

c        Initialize variables
         call initvarsm(ns,yes,nearaxis,smallr,corrmu,corrmv,
     +                  corrbu,corrbv,width,xmax,uextrap,vextrap)

c        Is (xx,yy) too close to where interface meets axis?
         if ( (xx**2+(yy-ys(0))**2).lt.0.2375**2) nearaxis=.true. !.2375=.95/4

c        Is (xx,yy) too close to sheet?
c        call x2close2S(j,xx,yy,xs,ys,ns,yes,k0,alfb,xproj,yproj,dist,dcut,dels)
         call x2cl2circ(xx,yy,ns,yes,k0,alfb,xb,yb,dist)

c        Is radius xx too small, close to sheet?
         if ((xx.lt.xmax).and.(abs(yy-ys(0)).lt.width)) smallr=.true. 

         if (yes) then   
            call correctmnearsheet(xx,yy,alf,ns,k0,alfb,xb,yb,dist,
     +             corrmu,corrmv,xdb,ydb)

           if (nearaxis) 
     +        call compncorrBm(xx,yy,alf,ns,xd(0),ys(0),f1d(0),f2(0),
     +             f2d2(0),alfb,xb,yb,xdb,ydb,dist,k0,
     +             corrbu,corrbv)
         else
            if (nearaxis) 
     +        call compBm(xx,yy,alf,ns,xd(0),ys(0),f1d(0),f2(0),f2d2(0),
     +             corrbu,corrbv)
         endif
         if (smallr) then
            call compsmallrm(xx,yy,xmax,alf,xs,ys,ns,f1,f2,
     +                   xd(0),xddd(0), ydd(0), f1d(0), f2d2(0),
     +                   xd(ns),xddd(ns),ydd(ns),f1d(ns),f2d2(ns),
     +                   uextrap,vextrap)          
c            print*,uextrap,vextrap
         endif

c print data for plotting
c         if (largecorr) then 
c            write(8,1000)xx,yy, udb(j),vdb(j),
c     +           corrmu,corrmv, corrbu,corrbv,
c     +           ucorrected+corrbu,vcorrected+corrbv
         write(8,1000)xx,yy, u(j),v(j),
     +        corrmu,corrmv, corrbu,corrbv, uextrap,vextrap
      enddo

#ifdef QUAD
1000  format(2f14.10,8e29.20)
#endif
#ifdef DBLE
1000  format(2f14.10,8e23.14)
#endif

      return
      end


      SUBROUTINE initvarsm(ns,yes,nearaxis,smallr,corrmu,corrmv,
     +                  corrbu,corrbv,width,xmax,uextrap,vextrap)
#include "params.h"
      integer ns
      MODE corrmu,corrmv,width,corrbu,corrbv,xmax,uextrap,vextrap
      logical yes,nearaxis,smallr

      yes=.false.
      nearaxis=.false.
      smallr=.false.
      corrmu=0
      corrmv=0
      corrbu=0
      corrbv=0
      uextrap=0
      vextrap=0
      if (ns.eq.100) then
c           width=0.24
c           xmax=0.03
           width=0.1
           xmax=0.02
      elseif (ns.eq.200) then
c           width=0.12
c           xmax=0.03
           width=0.04
           xmax=0.02
      elseif (ns.eq.400) then
c           width=0.12
c           xmax=3*pi/ns
           width=0.0
           xmax=-1
      elseif (ns.eq.800) then
c           xmax=4*pi/ns
           width=0.0
           xmax=-1
      endif
      return
      end

