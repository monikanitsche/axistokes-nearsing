c Contains
c SUBROUTINE IntQs: comp ud,vd at x,y due to xs,ys and then corrects
c SUBROUTINE findclosest: (not used) generic routine to find closest pt
c                        on interface
c SUBROUTINE initvars: init vars to .false., =0, and vals of xmax,width 
c SUBROUTINE fixroundoff: routine to replace hj-gj by rest if largecorr
c                        (not used, not needed here)

#include "flags.h"
      SUBROUTINE IntQs(xs,ys,ns,u,v,x,y,n,udb,vdb)
c      SUBROUTINE IntDb2on1(xs,ys,ns,u,v,x,y,n,ud,vd)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Computes dble layer integral udb,vdb 
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
c     output: udb,vdb  (dble layer integral at x,y on SHEET 2 induced by SHEET 1)
c
c     Variables:
c       subscript s: sheet
c       subscript a: approximating sheet (for B)
c       subscript p or b: projection 
c       subscript d: derivatives
c       subscript db: double layer as in udb,vdb
c       subscript sg: single layer as in usg,vsg
c       subscript dbs: self-induced double layer velo on sheet, as in udbs,vdbs
c
c       xd,xdd,xddd, yd,ydd,yddd: derivatives of xs,ys at alfj j=0,ns
c       xp,yp,alfb,xdp,xddp,xdddp,etc : projection: values at alfb
c       udbsp,vdbsp: self-induced dble layer velo at alfb
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
      MODE xs(0:*),ys(0:*),u(0:*),v(0:*),x(0:*),y(0:*)
C     OUTPUT
      MODE udb(0:*),vdb(0:*)
C     LOCAL
      integer j,k0
      MODE xx,yy,alfb,xb,yb,dist,corrqu,corrqv,corrbu,corrbv,
     +       ucorrected,vcorrected,uextrap,vextrap
      MODE ud3b,zd4b,xjunk,zero,width,corrsize
      MODE ud(0:nmax),ud2(0:nmax),vd(0:nmax),vd2(0:nmax),xdb,ydb,
     +       restu,restv,hk,rest1,xmax
      logical yes,nearaxis,smallr,largecorr

c compute derivatives of "densities", needed below
      call compderiv4(u,ud,ud2,v,vd,vd2,ns)

c compute velocity at (x,y) induced by (xs,ys) using T4[G] or T6[G]     
      call velodSonO(xs,ys,ns,u,v,x,y,n,udb,vdb,xd,yd,
     +   xddd(0),ydd(0),ud(0),vd2(0),xddd(ns),ydd(ns),ud(ns),vd2(ns))

c compute correction for target points on axis (j=0)
      j=0
      alfb=0
      yy=y(j)
      dist=abs(yy-ys(0))
      corrqv=0
      if (dist.lt.0.4) then
         ud3b=(u(2)-2*u(1))/delalf**3           !u odd
         zd4b=2*(ys(2)-4*ys(1)+3*ys(0))/delalf**4  !z even
         call correctionqendsh(alf,ns,yy,ys(0),dist,alfb,corrqv,
     +     xd(0),xddd(0),ydd(0),zd4b,ud(0),ud3b,v(0),vd2(0))
      endif
      zero=0
      write(9,1000)x(j),y(j), udb(j),vdb(j),
     +        zero,corrqv, zero,zero, zero,zero

c compute corrections for target points j=1..n-1 not on axis
      do j = 1,n-1  
         xx=x(j)
         yy=y(j)

c        Initialize variables
         call initvars(ns,yes,nearaxis,smallr,largecorr,
     +                  corrbu,corrbv,corrqu,corrqv,width,xmax)

c        Is (xx,yy) too close to where interface meets axis?
         if ( (xx**2+(yy-ys(0))**2).lt.0.2375**2) nearaxis=.true. !.2375=.95/4

c        Is (xx,yy) too close to sheet?
c        call x2close2S(j,xx,yy,xs,ys,ns,yes,k0,alfb,xproj,yproj,dist,dcut,dels)
         call x2cl2circ(xx,yy,ns,yes,k0,alfb,xb,yb,dist)

c        Is radius xx too small, close to sheet?
         xmax=min(0.03,4*delalf)  !this was used by mistake for 
c                                 !those .r3 results! worked like a charm
c                                 !should comment out in future
         if ((xx.lt.xmax).and.(abs(yy-ys(0)).lt.width)) smallr=.true. 

         if (yes) then   
            xjunk=0
            call correctqnearsheetnew(alf,ns,xx,yy,alfb,xb,yb,dist,
     +             corrqu,corrqv,xdb,ydb,largecorr,k0,restu,restv,
     +             hk,rest1,xjunk,xjunk,xjunk,xjunk,xjunk,xjunk)

#ifdef ROUNDOFFFIX
c not needed for current runs since corrsize is never too large (d not tiny) 
           corrsize= max(abs(corrqu),abs(corrqv))
           if ((corrsize.ge.10).and.(.not.nearaxis)) largecorr=.true. 
           if (largecorr) 
     +        call fixroundoff(alf,xs,ys,ns,u,v,xd,yd,
     +               xddd(0), ydd(0), ud(0), vd2(0),
     +               xddd(ns),ydd(ns),ud(ns),vd2(ns),
     +               xx,yy,alfb,xb,yb,dist,k0,xdb,ydb,
     +               ucorrected,vcorrected)
#endif
c          if nearaxis and tooclose 
c          (need corrections of B, alfp known)
           if (nearaxis) 
     +        call compncorrBq(xx,yy,alf,ns,xd(0),ys(0),ud(0),v(0),
     +             vd2(0),alfb,xb,yb,xdb,ydb,dist,k0,
c     +             bintu,bintv,btrapu,btrapv)
     +             corrbu,corrbv)
         else
c           if nearaxis but not tooclose 
c           (no corrections of B, need to find alfp)
            if (nearaxis) 
     +         call compBq(xx,yy,alf,ns,
     +                     xd(0),ys(0),ud(0),v(0),vd2(0),
     +             corrbu,corrbv)
         endif
         if (smallr) then
            call compsmallr(xx,yy,xmax,alf,xs,ys,ns,u,v,xd,yd,
     +                   xddd(0), ydd(0), ud(0), vd2(0),
     +                   xddd(ns),ydd(ns),ud(ns),vd2(ns),
     +                   uextrap,vextrap)          
        endif

c print data for plotting
         if (largecorr) then 
            write(9,1000)xx,yy, udb(j),vdb(j),
     +           corrqu,corrqv, corrbu,corrbv,
     +           ucorrected+corrbu,vcorrected+corrbv
         elseif (smallr) then 
             write(9,1000)xx,yy, udb(j),vdb(j),
     +        corrqu,corrqv, corrbu,corrbv, uextrap,vextrap
         else
c           write(9,1000)xx,yy,udb(j),vdb(j),udb(j)+corrqu,vdb(j)+corrqv,
c     +      bintu,bintv,btrapu,btrapv
           write(9,1000)xx,yy, udb(j),vdb(j),
     +        corrqu,corrqv, corrbu,corrbv, zero,zero
         endif
      enddo

#ifdef QUAD
1000  format(2f14.10,8e29.20)
#endif
#ifdef DBLE
1000  format(2f14.10,8e23.14)
#endif

      return
      end


      SUBROUTINE findclosest(j0,xx,yy,xs,ys,ns)
c find closest on coarsest mesh, for sake of estimating errors
      implicit none
      integer j,j0,ns,jstep
      MODE xx,yy,xs(0:*),ys(0:*)
      MODE dist, distnew
      dist=1
      jstep=ns/100
c      print*,jstep
c      do j=1,ns/3
      do j=jstep,ns/3,jstep
        distnew=(xs(j)-xx)**2+(ys(j)-yy)**2
c        print*,j,dist,distnew,xs(j),xx
        if ((distnew.ge.dist).or.(xs(j).ge.xx)) then
          j0=j
          goto 100
        endif
        dist=distnew
      enddo
      print*,'WARNING','j0 too large'
      j0=ns/3
100   continue
      return
      end


      SUBROUTINE initvars(ns,yes,nearaxis,smallr,largecorr,
c     +                  bintu,bintv,btrapu,btrapv,corrqu,corrqv,width)
     +                  corrbu,corrbv,corrqu,corrqv,width,xmax)
#include "params.h"
      integer ns
      MODE corrqu,corrqv,width,corrbu,corrbv,xmax
      logical yes,nearaxis,smallr,largecorr

      yes=.false.      
      nearaxis=.false.
      smallr=.false.
      largecorr=.false.
c      bintu=0
c      bintv=0
c      btrapu=0
c      btrapv=0
      corrbu=0
      corrbv=0
      corrqu=0
      corrqv=0
      if (ns.eq.100) then
           width=0.24
           xmax=0.03
      elseif (ns.eq.200) then
           width=0.12
           xmax=0.03
      elseif (ns.eq.400) then
           width=0.12
           xmax=3*pi/ns
      elseif (ns.eq.800) then
           width=0.07
           xmax=4*pi/ns
      endif
      return
      end



      SUBROUTINE fixroundoff(alf,xs,ys,ns,u,v,xd,yd,
     +               xddd0, ydd0, ud0, vd20,
     +               xdddn, yddn, udn, vd2n,
     +               xx,yy,alfb,xb,yb,dist,k0,xdb,ydb,
     +               ucorrected,vcorrected)
#include "params.h"
c INPUT
      integer ns,k0
      MODE alf(0:*),xs(0:*),ys(0:*),u(0:*),v(0:*),xd(0:*),yd(0:*),
     +     xddd0,ydd0,ud0,vd20,xdddn,yddn,udn,vd2n,
     +     xx,yy,alfb,xb,yb,dist,xdb,ydb
c OUTPUT
      MODE ucorrected,vcorrected
c LOCAL
      MODE xk,xik,uxk,vxk,uyk,vyk,restu,restv,hk,rest1,gkself,
     +  corrqunk,corrqvnk,udnk,vdnk,delalf
      logical largecorr

      delalf=pi/ns
      xk=xs(k0)
      xik=ys(k0)-yy
      uxk=u(k0)*xd(k0)
      vxk=v(k0)*xd(k0)
      uyk=u(k0)*yd(k0)
      vyk=v(k0)*yd(k0)
      largecorr=.true.

      call correctqnearsheetnew(alf,ns,xx,yy,alfb,xb,yb,dist,
     +            corrqunk,corrqvnk,xdb,ydb,largecorr,k0,restu,restv,
     +            hk,rest1,xk,xik,uxk,vxk,uyk,vyk)

      call velodSonPnok(xs,ys,ns,u,v,k0,xx,yy,udnk,vdnk,
     +               xd,yd,xddd0,ydd0,ud0,vd20,xdddn,
     +               yddn,udn,vd2n,gkself)
c     To test: include print statements in testrest.F HERE
      ucorrected=udnk+corrqunk+restu*delalf
      vcorrected=vdnk+corrqvnk+restv*delalf
      return
      end



c*******************
c USED FOR PLOTS OF U,VD,UDD ON AXIS (in compaxisplot)
c        call compaxisq(alf,yy,ys(0),dist,alfb,
c     +       xs,ys,ns,u,v,xd,yd,
c     +       xddd(0),xddd(ns),ydd(0),ydd(ns),
c     +       ud(0),ud(ns),vd2(0),vd2(ns),vax,udax,vddax)
c        write(10,*)0,yy,vax,udax,vddax
c*******************
