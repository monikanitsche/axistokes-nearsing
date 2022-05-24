#include "flags.h"
      SUBROUTINE correctmnearsheet(xx,yy,alf,ns,k0,alfb,xb,yb,dist,
     +         corrmu,corrmv,c1,d1) 
c     +         xd,xdd,xddd,yd,ydd,yddd,
c     +         f1,f1d,f1dd,f2,f2d,f2dd,
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Called when too close to sheet
c     Output: corrections E[H_Q] and E[H_M]
c     2/1/21: rewritten to minimize common blocks,
c             replaced hermite...by derivsatalfb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     INPUT
      integer k0,ns
      MODE xx,yy,alf(0:*),alfb,xb,yb,dist
c     OUTPUT      
      MODE corrmu,corrmv,c1,d1
c     LOCAL
      MODE dx,dy,c2,c3,c4,d2,d3,d4,
     +     f10,f11,f12,f13,f20,f21,f22,f23

c      sets values of all variables at alfb (=old hermite)
       call derivsmatalfp(
c input
     +    xx,yy,alfb,xb,yb,
c output
     +    dx,dy,c1,c2,c3,c4,d1,d2,d3,d4,
     +    f10,f11,f12,f13,f20,f21,f22,f23)
c       print*,k0,xb,xx,c1,c2,c3
c      print*,d1,d2,d3
c      print*,f10,f11,f12,f20,f21,f22

c     compute corrections for Mu,Mv
      call correctionmatxnew(
c input
     +        alf,alfb,ns,dist,k0,
     +        xb,xx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +        f10,f11,f12,f13,f20,f21,f22,f23,
c output
     +        corrmu,corrmv)
c      print*,corrmu,corrmv

      return
      end


      SUBROUTINE derivsmatalfp(xx,yy,alfb,xb,yb,
     +    dx,dy,c1,c2,c3,c4,d1,d2,d3,d4,
     +    f10,f11,f12,f13,f20,f21,f22,f23)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     interpolates variables at alf(kk-2:kk+3) to alfb
c     using HERMITE cubic interpolant,needs 6 points
c
c     for now: use known exact values
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     INPUT      
      MODE xx,yy,alfb,xb,yb
c     OUTPUT      
      MODE dx,dy,c1,c2,c3,c4,d1,d2,d3,d4,
     +    f10,f11,f12,f13,f20,f21,f22,f23
C     LOCAL
      MODE rad,f1,f1d,f1d2,f1d3,f2,f2d,f2d2,f2d3,
     +     xd,xd2,xd3,xd4,yd,yd2,yd3,yd4

c      call fithermite(alf(kk-2),alfp,vdd(kk-2),c)

       rad = 0.95q0

c values at basepoint alfb
       f1   =  sin(alfb)
       f1d  =  cos(alfb)
       f1d2 = -sin(alfb)
       f1d3 = -cos(alfb)

       f2   =  cos(alfb)
       f2d  = -sin(alfb)
       f2d2 = -cos(alfb)
       f2d3 =  sin(alfb)

       xd   =  rad*cos(alfb)
       xd2  = -rad*sin(alfb)
       xd3  = -rad*cos(alfb)
       xd4  =  rad*sin(alfb)

       yd   =  rad*sin(alfb)
       yd2  =  rad*cos(alfb)
       yd3  = -rad*sin(alfb)
       yd4  = -rad*cos(alfb)

       dx=xb-xx
       c1=xd
       c2=xd2/2
       c3=xd3/6

       dy=yb-yy
       d1=yd
       d2=yd2/2
       d3=yd3/6

       f10=f1
       f11=f1d
       f12=f1d2/2

       f20=f2
       f21=f2d
       f22=f2d2/2

#ifdef THIRD
       c4=0
       d4=0
       f13=0
       f23=0
#endif
#ifdef FOURTH
       c4=xd4/24
       d4=yd4/24
       f13=f1d3/6
       f23=f2d3/6
#endif

      return
      end
