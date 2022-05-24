c Contains:
c SUBROUTINE correctqnearsheetnew
c           outputs corrections qu,qv
c SUBROUTINE derivsqatalfpnew
c           computes coefficients needed for correqtionqatxnew
c           currently simply by using exact solutions
#include "flags.h"

      SUBROUTINE correctqnearsheetnew(alf,ns,xx,yy,alfb,xb,yb,dist,
     +      corrqu,corrqv,c1,d1,tooclose,k0,restu,restv,
     +      hk,rest1,xk,xik,uxk,vxk,uyk,vyk)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Called when too close to sheet
c     Output: corrections E[H_Q] and E[H_M]
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     INPUT
      integer ns,k0
      MODE alf(0:*),xx,yy,alfb,xb,yb,dist,xk,xik,
     +      uxk,vxk,uyk,vyk
      logical tooclose
c     OUTPUT
      MODE corrqu,corrqv,c1,d1,restu,restv
c     LOCAL
      MODE dx,dy,c2,c3,c4,d2,d3,d4,
     +  ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +  uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3,hk,rest1

c      sets coefficients of expansions of all functions about alfb 
c      for now: short and quick, not general
       call derivsqatalfpnew(
c input
     +  xx,yy,alfb,xb,yb,  
c output
     +  dx,dy,c1,c2,c3,c4,d1,d2,d3,d4,
     +  ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +  uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3)

c     compute corrections  for Qu,Qv   
      call correctionqatxnew(
c input
     +      alf,alfb,ns,dist,k0,xk,xik,
     +      xb,xx,c1,c2,c3,c4,dy,d1,d2,d3,d4, 
     +      ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3, 
     +      uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3, 
     +      uxk,vxk,uyk,vyk,
c output
     +      corrqu,corrqv,tooclose,restu,restv,hk,rest1) 

      return
      end


      SUBROUTINE derivsqatalfpnew(xx,yy,alfb,xb,yb,  !input
c output
     +  dx,dy,c1,c2,c3,c4,d1,d2,d3,d4,
     +  ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +  uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3)
c(alf,alfp,kk,xd,xdd,xddd,yd,ydd,yddd,
c    +     u,ud,udd,v,vd,vdd,
c    +     xdp,xddp,xdddp,ydp,yddp,ydddp,
c    +     up,udp,uddp,vp,vdp,vddp)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     interpolates variables at alf(kk-2:kk+3) to alfp
c     using HERMITE cubic interpolant,needs 6 points
c
c     input: alf,alfp,kk,xd,xdd,xddd,yd,ydd,yddd 
c     output: xdp,xddp,xdddp,ydp,yddp,ydddp values at alfp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c INPUT
      MODE xx,yy,alfb,xb,yb  !input
c OUTPUT
      MODE dx,dy,c1,c2,c3,c4,d1,d2,d3,d4,
     +  ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +  uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3
c LOCAL
      MODE rad,xd,xd2,xd3,xd4,yd,yd2,yd3,yd4,u,ud,ud2,ud3,v,vd,vd2,
     +  vd3

       rad=0.95q0
c      xb  = rad*sin(alfb)
       xd  = rad*cos(alfb)
       xd2 =-rad*sin(alfb)
       xd3 =-rad*cos(alfb)

c      yb  =-rad*cos(alfb)-4
       yd  = rad*sin(alfb)
       yd2 = rad*cos(alfb)
       yd3 =-rad*sin(alfb)

       u   = sin(alfb)
       ud  = cos(alfb)
       ud2 =-sin(alfb)
       v   = cos(alfb)
       vd  =-sin(alfb)
       vd2 =-cos(alfb)

c four third order method set these derivatives to zero
#ifdef THIRD
       xd4=0
       yd4=0
       ud3=0
       vd3=0
#endif

#ifdef FOURTH
       xd4 = rad*sin(alfb)
       yd4 =-rad*cos(alfb)
       ud3 =-cos(alfb)
       vd3 = sin(alfb)
#endif

       call compallcoeffuys(
c input
     +    xb,xx,xd,xd2,xd3,xd4,
     +    yb,yy,yd,yd2,yd3,yd4,
     +    u,ud,ud2,ud3,v,vd,vd2,vd3,
c output
     +    dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +    uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3)

      return
      end
