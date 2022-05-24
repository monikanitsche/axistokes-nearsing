cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Contains
c SUBROUTINE compallcoeffuys : turns derivatives into coefficients
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include "flags.h"
      SUBROUTINE compallcoeffuys(
c input: all derivatives of x,y,u,v at basepoint (either alfp or betb)
     +    x,xx,xd,xd2,xd3,xd4,
     +    y,yy,yd,yd2,yd3,yd4,
     +    u,ud,ud2,ud3,v,vd,vd2,vd3,
c ouput: all coefficients of x,y,ux,uy,vx,vy of expansions
     +    dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +    uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c turns derivatives into coefficients
c input: all derivatives of x,y,u,v at basepoint (either alfp or betb)
c ouput: all coefficients of x,y,ux,uy,vx,vy of expansions
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     INPUT
      MODE
     +    x,xx,xd,xd2,xd3,xd4,
     +    y,yy,yd,yd2,yd3,yd4,
     +    u,ud,ud2,ud3,v,vd,vd2,vd3
c     OUTPUT
      MODE
     +    dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +    uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3

      dx=x-xx
      c1=xd
      c2=xd2/2
      c3=xd3/6
c      print*,real(xd),real(xd2),real(xd3),real(xd4),'herex'
c      print*,real(c1),real(c2),real(c3),real(c4),'herec'
      dy=y-yy
      d1=yd
      d2=yd2/2
      d3=yd3/6

      ux0=u*xd
      ux1=ud*xd+u*xd2
      ux2=(ud2*xd+2*ud*xd2+u*xd3)/2

      vx0=v*xd
      vx1=vd*xd+v*xd2
      vx2=(vd2*xd+2*vd*xd2+v*xd3)/2

      uy0=u*yd
      uy1=ud*yd+u*yd2
      uy2=(ud2*yd+2*ud*yd2+u*yd3)/2

      vy0=v*yd
      vy1=vd*yd+v*yd2
      vy2=(vd2*yd+2*vd*yd2+v*yd3)/2

#ifdef THIRD
      c4=0
      d4=0
      ux3=0
      vx3=0
      uy3=0
      vy3=0
#endif

#ifdef FOURTH
      c4=xd4/24
      d4=yd4/24
      ux3=(ud3*xd+3*ud2*xd2+3*ud*xd3+u*xd4)/6
      vx3=(vd3*xd+3*vd2*xd2+3*vd*xd3+v*xd4)/6
      uy3=(ud3*yd+3*ud2*yd2+3*ud*yd3+u*yd4)/6
      vy3=(vd3*yd+3*vd2*yd2+3*vd*yd3+v*yd4)/6
#endif

      return
      end


