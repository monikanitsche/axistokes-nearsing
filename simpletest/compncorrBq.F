ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Contains
c SUBROUTINE compncorrB : outputs bintu+corr, btrapu+corr 
c SUBROUTINE derivsqatpBbetnew: outputs all coefficients of x,y,uys 
c                               in expansions about alfp and about betb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "flags.h"
      SUBROUTINE compncorrBq(xx,yy,alf,ns,xd0,y0,ud0,v0,vd20,
     +     alfb,xb,yb,xdb,ydb,dist,k0,
     +     corrbu,corrbv)
c     +     bintu,bintv,btrapu,btrapv)
c     +     k0,xk,xik,uxk,vxk,uyk,vyk)
c   later: corrbu=(bintu+corrbintu) - (btrapu+corrtbrapu)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computes E[Bbet] = (intB+corrIntB)-(trapB+corrtrapB)
c Input:
c    xx,yy  target point
c    alf,ns parametrization of xs,ys
c    xd0,y0,ud0,v0,vd20: values on axis used to obtain quadratic
c                        or cubic approximations of xs,ys,u,v
c    alfb,xb,yb,xdb,ydb: values at projection=base point
c Output:
c   for now, corrected intbu =intbu+corrintbu
c            corrected trapbu=trapbu+corrtrapbu
c   later: corrbu=(intbu+corrintbu) - (trapbu+corrtrapbu)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns,k0
      MODE xx,yy,alf(0:*),xd0,y0,ud0,v0,vd20,
     +       alfb,xb,yb,xdb,ydb,dist
c OUTPUT
      MODE corrbu,corrbv
c LOCAL
      integer jmax,k0bet
      MODE bintu,bintv,btrapu,btrapv
      MODE alfend,cc1,cc2,dd1,dd2,
     +      uxk,vxk,uyk,vyk,xk,xik,
     +      corrbintu,corrbintv,corrbtrapu,corrbtrapv
      MODE bet(0:nmax),betb,
     +    dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +    uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3,
     +    c1bet,c2bet,c3bet,c4bet,d1bet,d2bet,d3bet,d4bet,
     +    ux0bet,ux1bet,ux2bet,ux3bet,vx0bet,vx1bet,vx2bet,vx3bet,
     +    uy0bet,uy1bet,uy2bet,uy3bet,vy0bet,vy1bet,vy2bet,vy3bet,
     +      uxkbet,vxkbet,uykbet,vykbet
      MODE restu,restv,hk,term1,xikbet,xkbet,delalf,xjunk
      MODE disttopt,btrapunk,btrapvnk,corrbtrapunk,corrbtrapvnk
      logical tooclose

      call intbq(xx,yy,ns,xd0,y0,ud0,v0,vd20,bintu,bintv,
     +           alfend,alfb,xb,yb,xdb,ydb,cc1,cc2,dd1,dd2,jmax,
     +           k0bet,xkbet,xikbet,uxkbet,vxkbet,uykbet,vykbet)

      call trapbq(xx,yy,alf,jmax,xd0,y0,ud0,v0,vd20,
     +     cc1,cc2,dd1,dd2,btrapu,btrapv)

      call derivsqatpBbetnew(xx,yy,ns,alfb,xb,yb,xdb,ydb,  !input
     +    xd0,y0,ud0,v0,vd20,alfend,cc1,cc2,dd1,dd2,
c output coefficients of xs,ys,ux,uy,vx,vy (alf)
     +    dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +    uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3,
c output coefficients of xs,ys,uy*dalf,ux*dalf,vx*dalf,vy*dalf (bet)
     +    bet,betb,c1bet,c2bet,c3bet,c4bet,d1bet,d2bet,d3bet,d4bet,
     +    ux0bet,ux1bet,ux2bet,ux3bet,vx0bet,vx1bet,vx2bet,vx3bet,
     +    uy0bet,uy1bet,uy2bet,uy3bet,vy0bet,vy1bet,vy2bet,vy3bet)

      tooclose=.false.
      corrbintu=0
      corrbintv=0
      corrbtrapu=0
      corrbtrapv=0
c     compute corrections  for B_Q(beta) 
c      if (dist.le.5*delsbet) call correctionqatxnew(

c    NOTE xx=x0
c      xk=x(alf(k))
c      xik=y(alf(k))-yy
      call correctionqatxnew(
     +      bet,betb,ns,dist,k0bet,xkbet,xikbet,
     +      xb,xx,c1bet,c2bet,c3bet,c4bet,dy,d1bet,d2bet,d3bet,d4bet,
     +      ux0bet,ux1bet,ux2bet,ux3bet,vx0bet,vx1bet,vx2bet,vx3bet,
     +      uy0bet,uy1bet,uy2bet,uy3bet,vy0bet,vy1bet,vy2bet,vy3bet,
     +      uxkbet,vxkbet,uykbet,vykbet,
     +      corrbintu,corrbintv,tooclose,restu,restv,  !output
     + hk,term1) !output

      xjunk=0
      call correctionqatxnew(
c     +      alf,alfb,jmax,dist,k0,xk,xik,
     +      alf,alfb,jmax,dist,k0,xjunk,xjunk,
     +      xb,xx,c1,c2,c3,c4,dy,d1,d2,d3,d4, !input
     +      ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3, !input
     +      uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3, !input
c     +      uxk,vxk,uyk,vyk,
     +      xjunk,xjunk,xjunk,xjunk,
     +      corrbtrapu,corrbtrapv,tooclose,restu,restv,
     + hk,term1) !output

      bintu=bintu+corrbintu
      bintv=bintv+corrbintv
      btrapu=btrapu+corrbtrapu
      btrapv=btrapv+corrbtrapv
      corrbu=bintu-btrapu
      corrbv=bintv-btrapv

      return
      end


      SUBROUTINE derivsqatpBbetnew(xx,yy,ns,alfp,x,y,xd,yd,  !input
     +    xd0,y0,ud0,v0,vd20,alfend,cc1,cc2,dd1,dd2,
c output coefficients of xs,ys,ux,uy,vx,vy (alf)
     +    dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +    uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3,
c output coefficients of xs,ys,uy*dalf,ux*dalf,vx*dalf,vy*dalf (bet)
     +    bet,betp,c1bet,c2bet,c3bet,c4bet,d1bet,d2bet,d3bet,d4bet,
     +    ux0bet,ux1bet,ux2bet,ux3bet,vx0bet,vx1bet,vx2bet,vx3bet,
     +    uy0bet,uy1bet,uy2bet,uy3bet,vy0bet,vy1bet,vy2bet,vy3bet)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Sets derivatives at alp of of all functions of alf (MUST BE 
c                  consistent with definition of B(bet) in intbq)
c     Then finds derivatives at betp where alfp=alf(betp)
c     of of all as functions of bet. obtained by chain rule
c
c     output: 
c     coefficients of expansions in alf about alfp
c        x(alf)=dx+c1*alf+c2*alf^2+...+c4*alf^4
c        y(alf)=dy+d1*alf+d2*alf^2+...+d4*alf^4
c        ux(alf)=u(alf)*dx(alf)=ux0+ux1*alf+...
c        uy(alf)=u(alf)*dy(alf)=uy0+uy1*alf+...
c        vx(alf)=v(alf)*dx(alf)=vx0+vx1*alf+...
c        vy(alf)=v(alf)*dy(alf)=vy0+vy1*alf+...
c     coefficients of expansions in bet about betp
c        x(bet)=dx+c1*bet+c2*bet^2+...+c4*bet^4
c        y(bet)=dy+d1*bet+d2*bet^2+...+d4*bet^4
c        ux(bet)=u(bet)*dx(bet)=ux0bet+ux1bet*bet+...
c        uy(bet)=u(bet)*dy(bet)=uy0bet+uy1bet*bet+...
c        vx(bet)=v(bet)*dx(bet)=vx0bet+vx1bet*bet+...
c        vy(bet)=v(bet)*dy(bet)=vy0bet+vy1bet*bet+...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,ns
c INPUT      
      MODE xx,yy,alfp,x,y,xd,yd,
     +    xd0,y0,ud0,v0,vd20,alfend,cc1,cc2,dd1,dd2
c OUTPUT 
      MODE dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +    uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3,
     +    bet(0:*),betp,c1bet,c2bet,c3bet,c4bet,d1bet,d2bet,d3bet,d4bet,
     +    ux0bet,ux1bet,ux2bet,ux3bet,vx0bet,vx1bet,vx2bet,vx3bet,
     +    uy0bet,uy1bet,uy2bet,uy3bet,vy0bet,vy1bet,vy2bet,vy3bet
c     LOCAL      
      MODE u,ud,ud2,ud3,v,vd,vd2,vd3,xd2,xd3,xd4,yd2,yd3,yd4,
     +  delbet,alfd,alfd2,alfd3,alfd4,
     +  udbet,ud2bet,ud3bet,vdbet,vd2bet,vd3bet,
     +  xdbet,xd2bet,xd3bet,xd4bet,ydbet,yd2bet,yd3bet,yd4bet

      delbet=pi/ns
      do k=0,ns
         bet(k)=k*delbet
      enddo

c     THESE ARE DERIVATIVES WRT ALFA at alfp
      u   = ud0*alfp
      ud  = ud0
      ud2 = 0
      ud3 = 0

      v   = v0+vd20*alfp**2/2
      vd  = vd20*alfp
      vd2 = vd20
      vd3 = 0

      xd2 = 6*cc1*alfp   +20*cc2*alfp**3
      xd3 = 6*cc1        +60*cc2*alfp**2
c      yddp =      ydd0       +   12*dd1*alfp**2+30*dd2*alfp**4
c      ydddp=                     24*dd1*alfp  +120*dd2*alfp**3
      yd2 = 2*dd1        +12*dd2*alfp**2
      yd3 =               24*dd2*alfp

#ifdef THIRD
      xd4=0
      yd4=0
#endif

#ifdef FOURTH
      xd4 =              120*cc2*alfp
      yd4 =               24*dd2
#endif

      call compallcoeffuys(
c input
     +    x,xx,xd,xd2,xd3,xd4,
     +    y,yy,yd,yd2,yd3,yd4,
     +    u,ud,ud2,ud3,v,vd,vd2,vd3,
c output
     +    dx,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +    ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +    uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3)

c     NOW FIND DERIVATIVES WRT BETA
      betp  = 2*acos(1-alfp/alfend)      !projection xp=x(betp)
      alfd  =  sin(betp/2)*alfend/2   !dalfdbet at bet=betp
      alfd2 =  cos(betp/2)*alfend/4
      alfd3 = -sin(betp/2)*alfend/8
      alfd4 = -cos(betp/2)*alfend/16

c      delsbet=dels*alfd
c derivatives of x(alf(bet)) wrt bet
c derivatives of u(alf(bet)) wrt bet
c v=v0, y=y0 are constant, xp0=xp0
      xdbet  = xd*alfd
      xd2bet = xd2*alfd**2+xd*alfd2
      xd3bet = xd3*alfd**3 + 3*xd2*alfd*alfd2 + xd*alfd3
      ydbet  = yd*alfd
      yd2bet = yd2*alfd**2+yd*alfd2
      yd3bet = yd3*alfd**3 + 3*yd2*alfd*alfd2 + yd*alfd3
      udbet  = ud*alfd
      ud2bet = ud2*alfd**2+ud*alfd2
      vdbet  = vd*alfd
      vd2bet = vd2*alfd**2+vd*alfd2

#ifdef THIRD
      ud3bet=0
      vd3bet=0
      xd4bet=0
      yd4bet=0
#endif

#ifdef FOURTH
      ud3bet = ud3*alfd**3+ 3*ud2*alfd*alfd2 + ud*alfd3
      vd3bet = vd3*alfd**3+ 3*vd2*alfd*alfd2 + vd*alfd3
      xd4bet = alfd**4*xd4 +6*alfd**2*alfd2*xd3
     +  +4*alfd*alfd3*xd2 +3*alfd2**2*xd2 +xd*alfd4
      yd4bet = alfd**4*yd4 +6*alfd**2*alfd2*yd3
     +  +4*alfd*alfd3*yd2 +3*alfd2**2*yd2 +yd*alfd4
#endif

      call compallcoeffuys(
c input
     +    x,xx,xdbet,xd2bet,xd3bet,xd4bet,
     +    y,yy,ydbet,yd2bet,yd3bet,yd4bet,
     +    u,udbet,ud2bet,ud3bet,v,vdbet,vd2bet,vd3bet,
c output
     +    dx,c1bet,c2bet,c3bet,c4bet,dy,d1bet,d2bet,d3bet,d4bet,
     +    ux0bet,ux1bet,ux2bet,ux3bet,vx0bet,vx1bet,vx2bet,vx3bet,
     +    uy0bet,uy1bet,uy2bet,uy3bet,vy0bet,vy1bet,vy2bet,vy3bet)

      return
      end


      SUBROUTINE intbq(xx,yy,ns,xd0,y0,ud0,v0,vdd0,bintu,bintv,
     +                 alfend,alfp,xp,yp,xdp,ydp,cc1,cc2,dd1,dd2,jmax,
     +                 k0,xk,xik,uxk,vxk,uyk,vyk)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute "exact" integral of b for qu,qv, to test
c     needed to get the error E[b]
c     uses ns tchebb pts on [0,alfend] plus 2nd and 4th order corrections
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns                !length of xs, and of xsa(bet)
      MODE xx,yy,               !target point
     +     xd0,y0,ud0,v0,vdd0,  !values of xs,ys at alf=0
     +     alfp,xp,yp,xdp,ydp   !values at basepoint
c OUTPUT
      integer jmax         ! index of upper bound alfend=alf(jmax)
      integer k0           ! index of point x(k0),xi(k0) closest to xb,dy
      MODE bintu,bintv,           !exact integrals of Bu,Bv
     +     alfend,cc1,cc2,dd1,dd2 !coeff in approximation xsa,ysa
      MODE xk,xik,uxk,vxk,uyk,vyk !vals at k0,needed for "rest", not used
c LOCAL
      integer k
      MODE alf(0:nmax),alfd(0:nmax),delbet,xsa(0:nmax),xia(0:nmax),
     +    sumx,m,b,c,c2,c3,c5,eth,efh,efhp,ac,cb,sing,i50,i51,i52,i53,
     +    q11,q12,q21,q22,beta,ua(0:nmax),va(0:nmax),
     +    xda(0:nmax),yda(0:nmax),
     +    r,alfd2l,dalfql,dbetql,d3betql,dbetqr,d3betqr,
     +    xi0,aa,bb,distmin,distk
      MODE ksq(0:nmax),F(0:nmax),E(0:nmax),qu(0:nmax),qv(0:nmax)
      MODE eminusone,oneminksq(0:nmax)

      delbet = pi/ns  !=delalf
      jmax=nint(4*xx/(pi/100))*ns/100
c      jmax=19*ns/100
c      alfend=max(10*delbet,jmax*delbet)
c      print*,4*xx,alfend,jmax
      jmax=max(1*ns/100,jmax)
      alfend=jmax*delbet
c      print*,alfend
c****************
c      read from testxapprox.m

      cc1=(5*xp/alfp-4*xd0-xdp)/(2*alfp**2)
      cc2=(xdp+2*xd0-3*xp/alfp)/(2*alfp**4)

      dd1=(2*(yp-y0)/alfp-ydp/2)/alfp
      dd2=(ydp/2-(yp-y0)/alfp)/alfp**3
c      aa=(ydp-ydd0*alfp)/alfp**3
c      bb=((yp-y0)-ydd0*alfp**2/2)/alfp**4
c      dd2=(aa-4*bb)/(2*alfp**2)
c      dd1=-(aa-6*bb)/2
      

**************
c      xda=xd0   !first try straight line
c      yda=0        !0th order approximation
c      va=v0        !0th order approximation
      xi0=y0-yy    !0th order approximation, basepoint

c     compute b(xx,yy) where yy.ne.0 (point not on axis)

      distmin=1
      do k=1,ns+2
        beta=k*delbet
        alf(k)=(1-cos(beta/2))*alfend
        alfd(k)=sin(beta/2)*alfend/2
        alfd2l=alfend/4
c        alf(k)=(1-cos(beta))*alfend/2
c        alfd(k)=sin(beta)*alfend/2
c        alfd2l=alfend/2

c USE EXACT DERIVATIVE OF XSA,YSA-YY=XIA
        xsa(k)=xd0*alf(k) +   cc1*alf(k)**3  +   cc2*alf(k)**5
        xda(k)=xd0        + 3*cc1*alf(k)**2  + 5*cc2*alf(k)**4
        xia(k)=(y0-yy) +  dd1*alf(k)**2 +  dd2*alf(k)**4    !through alfp
        yda(k)=         2*dd1*alf(k)    +4*dd2*alf(k)**3    !through alfp
c        xia(k)=(y0-yy)+ydd0*alf(k)**2/2+ dd1*alf(k)**4+   dd2*alf(k)**6
c        yda(k)=        ydd0*alf(k)    +4*dd1*alf(k)**3+ 6*dd2*alf(k)**5
        ua(k) =ud0*alf(k)  !1st order approximation
        va(k) =v0 +vdd0*alf(k)**2/2  !1st order approximation

        sumx = xsa(k) + xx
        ksq(k) = 4*xsa(k)*xx/(xia(k)**2 + sumx**2)
        oneminksq(k) = (xia(k)**2+(xsa(k)-xx)**2)/(xia(k)**2 + sumx**2)

        distk= xia(k)**2+ (xsa(k)-xx)**2
        if (distk.lt.distmin) then
                k0=k
                distmin=distk
        endif
      enddo
      xk=xsa(k0)
      xik=xia(k0)
      uxk=ua(k0)*xda(k0)
      vxk=va(k0)*xda(k0)
      uyk=ua(k0)*yda(k0)
      vyk=va(k0)*yda(k0)
c      print*,'k0bet',k0,xk,xik
 
c      call ellip(ksq(1),F(1),E(1),ns+1)
      call ellip1mink(ksq(1),oneminksq(1),F(1),E(1),ns+1)
 
      do k = 1,ns+2
c       if (k.eq.40) print*,xsa(k),F(k),E(k),ksq(k)
         sumx = xsa(k) + xx

c         m =1-ksq(k)
         m=oneminksq(k)
         eth = E(k)/m
         efh = (2*(1+m)/m*E(k)-F(k))/(3*m)
c         efhp = (2*(1+m)/m*E(k)-F(k))/(3*m) -2/(3*m**2)

         eminusone=E(k)-1
        if (m.lt.1.e-4) eminusone=m*(ec1+ec2*m+ec3*m**2)
     +                  -m*log(m)*(ed1+ed2*m+ed3*m**2)

         efhp=(2*(1+m)/m*eminusone-F(k)+2)/(3*m)  !efhp

         b=2*xx*xsa(k)
         c2=sumx**2 + xia(k)**2
         c=sqrt(c2)
         c3=c*c2
         c5=c2*c3
 
         ac=(1+m)/2
         cb=2/(1-m)
c         sing=0
c            i50= 4/c5*efh
c            i51= 4/(c3*b)*(ac*efh-eth) 
c            i52= 4/(c*b**2)*(ac**2*efh-2*ac*eth+F(k))
c            i53= 4*c/b**3*(ac**3*efh-3*ac**2*eth+3*ac*F(k)-E(k))

cc            i50= 4/c5*efh
cc            i51= 4/c5*cb*(ac*efh-eth) 
cc            i52= 4/c5*cb**2*(ac**2*efh-2*ac*eth+F(k))
cc            i53= 4/c5*cb**3*(ac**3*efh-3*ac**2*eth+3*ac*F(k)-E(k))

         sing = 8/(3*c5*m**2)
         i50= 4/c5*(efhp) 
         i51= 4/c5*(cb*(ac*efhp-eth)+4/(3*m*(1-m)) ) 
         i52=4/c5*(cb**2*(ac**2*efhp-2*ac*eth+F(k))+8/(3*m*(1-m)**2))
         i53= 4/c5*( cb**3*(ac**3*efhp-3*ac**2*eth+3*ac*F(k)-E(k)) 
     +        +4*(3+m**2)/(3*m*(1-m)**3) )

         q11=-6*xsa(k)*( xsa(k)**3*i51 - xsa(k)**2*xx*(i50+2*i52) +
     +       xsa(k)*xx**2*(i53+2*i51) - xx**3*i52 + (xsa(k)-xx)**3*sing)
                
         q12=-6*xsa(k)*xia(k)*((xsa(k)**2+xx**2)*i51-xsa(k)*xx*(i50+i52)
     +                    + (xsa(k)-xx)**2*sing)
         q21= q12
         q22= -6*xsa(k)*xia(k)**2*( xsa(k)*i51 -xx*i50+(xsa(k)-xx)*sing)

         qu(k)=(q11*yda(k)-q12*xda(k))*ua(k)
     +       + (q21*yda(k)-q22*xda(k))*va(k)


         q11=-6*xsa(k)*xia(k)*( xx**2*i52 +xsa(k)**2*i50-2*xsa(k)*xx*i51
     +                    + (xx-xsa(k))**2*sing)
         q12=-6*xsa(k)*xia(k)**2*( xsa(k)*i50 -xx*i51 +(xsa(k)-xx)*sing)
         q21= q12
         q22= -6*xsa(k)*xia(k)**3*(i50+sing)

         qv(k)=(q11*yda(k)-q12*xda(k))*ua(k)
     +        +(q21*yda(k)-q22*xda(k))*va(k)

         qu(k)=qu(k)*alfd(k)
         qv(k)=qv(k)*alfd(k)

      enddo
c         q(ns)=0
c         leftv(j)=(-3*q(0)+4*q(1)-q(2))/(2*delalf)
c         rightv(j)=(3*q(ns)-4*q(ns-1)+q(ns-2))/(2*delalf)
 
      bintu=0
      bintv=0
      do k = 1,ns-1
         bintu = bintu+ qu(k)
         bintv = bintv+ qv(k)
      enddo
      bintu=bintu*delbet
      bintv=bintv*delbet

c add corrections using 1st and 3rd derivatives at left and right
c B(bet)=B(alf(bet))*alf'(bet)
c At left:
c first derivative dB/dbet=0 at left, since alf=0
c third derivatives dB^3/dbet^3 only requires dB/dalf at left, see notes
c At right:
c use finite difference approximations: 4th order for q', 2nd order for q'''
c see "Chapter12.pdf", "fd.pdf",  for reference

      r =sqrt(xi0**2+xx**2)

      dbetql=0   !since alfd=0 at left
      dalfql=-12*(pi*v0*xd0**2*xi0**2*xx)/r**5  !including c1,c2,d2,d2
      d3betql=3*dalfql*alfd2l**2 
      dbetqr=(-qu(ns+2)+8*qu(ns+1)-8*qu(ns-1)+qu(ns-2))/(12*delbet)
      d3betqr=(qu(ns+2)-2*qu(ns+1)+2*qu(ns-1)-qu(ns-2))/(2*delbet**3)
      bintu= bintu+qu(ns)*delbet/2
     +    - delbet**2*(dbetqr -dbetql )/12
     +    + delbet**4*(d3betqr-d3betql)/720

      dbetql=0   !since alfd=0 at left
      dalfql=12*(pi*v0*xd0**2*xi0**3)/r**5  !unchanged 
      d3betql=3*dalfql*alfd2l**2 
      dbetqr=(-qv(ns+2)+8*qv(ns+1)-8*qv(ns-1)+qv(ns-2))/(12*delbet)
      d3betqr=(qv(ns+2)-2*qv(ns+1)+2*qv(ns-1)-qv(ns-2))/(2*delbet**3)
      bintv= bintv+qv(ns)*delbet/2
     +    - delbet**2*(dbetqr -dbetql )/12
     +    + delbet**4*(d3betqr-d3betql)/720
 
      bintu = bintu/pi8
      bintv = bintv/pi8

      return
      end


      SUBROUTINE trapbq(xx,yy,alf,jmax,xd0,y0,ud0,v0,vdd0,
     +  cc1,cc2,dd1,dd2,trapbu,trapbv)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute trapezoid approximation T_6 of B_q for qu,qv, to test
c     needed to get the error E[b]
c     uses ns tchebb pts on [0,alfend] plus 2nd and 4th order corrections
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,jmax,ns
      MODE xx,yy,xd0,y0,ud0,v0,trapbu,trapbv,cc1,cc2,dd1,dd2,alf(0:*)
      MODE delalf,xi0,va(0:nmax),xsa(0:nmax),xia(0:nmax),xda(0:nmax),
     +    yda(0:nmax),ua(0:nmax),sumx,m,b,c,c2,c3,c5,eth,efh,efhp,ac,cb,
     +    i50,i51,i52,i53,q11,q12,q21,q22,r,dql,d3ql,dqr,d3qr,f1,f2,f3
      MODE ksq(0:nmax),F(0:nmax),E(0:nmax),qu(0:nmax),qv(0:nmax)
      MODE xddd0,ydd0,vdd0,alfa,dfine,sing,rsq,d3qln
      MODE eminusone, oneminksq(0:nmax)

      ns=jmax
      delalf = alf(jmax)/jmax  !=delalf
      dfine=delalf/100
c      print*,delalf
c      stop
c      va=v0        !0th order approximation
      xi0=y0-yy    !0th order approximation, basepoint
c      vdd0=0
      xddd0=6*cc1
      ydd0=2*dd1

c     compute b(xx,yy) where yy.ne.0 (point not on axis)

      do k=1,ns
        xsa(k)=xd0*alf(k) +   cc1*alf(k)**3  +   cc2*alf(k)**5
        xda(k)=xd0        + 3*cc1*alf(k)**2  + 5*cc2*alf(k)**4
        xia(k)=(y0-yy) +  dd1*alf(k)**2 +  dd2*alf(k)**4  
        yda(k)=         2*dd1*alf(k)    +4*dd2*alf(k)**3   
c        xia(k)=(y0-yy)+ydd0*alf(k)**2/2+ dd1*alf(k)**4+   dd2*alf(k)**6
c        yda(k)=        ydd0*alf(k)    +4*dd1*alf(k)**3+ 6*dd2*alf(k)**5
        ua(k) =ud0*alf(k)  
        va(k) =v0 +vdd0*alf(k)**2/2  

        sumx = xsa(k) + xx
        ksq(k) = 4*xsa(k)*xx/( xia(k)**2 + sumx**2)
        oneminksq(k) = (xia(k)**2+(xsa(k)-xx)**2)/( xia(k)**2 + sumx**2)
      enddo
      do k=ns+1,ns+5
        alfa=alf(ns)+(k-ns-3)*dfine
        xsa(k)=xd0*alfa +   cc1*alfa**3  +   cc2*alfa**5
        xda(k)=xd0        + 3*cc1*alfa**2  + 5*cc2*alfa**4
        xia(k)=(y0-yy) +  dd1*alfa**2 +  dd2*alfa**4  
        yda(k)=         2*dd1*alfa    +4*dd2*alfa**3   
c        xia(k)=(y0-yy)+ydd0*alfa**2/2+ dd1*alfa**4+   dd2*alfa**6
c        yda(k)=        ydd0*alfa    +4*dd1*alfa**3+ 6*dd2*alfa**5
        ua(k) =ud0*alfa  
        va(k) =v0 +vdd0*alfa**2/2  

        sumx = xsa(k) + xx
        ksq(k) = 4*xsa(k)*xx/( xia(k)**2 + sumx**2)
        oneminksq(k) = (xia(k)**2+(xsa(k)-xx)**2)/(xia(k)**2 + sumx**2)
      enddo
 
c      call ellip(ksq(1),F(1),E(1),ns+4)
      call ellip1mink(ksq(1),oneminksq(1),F(1),E(1),ns+4)
 
      do k = 1,ns+5
         sumx = xsa(k) + xx

c         m =1-ksq(k)
         m =oneminksq(k)

         eth = E(k)/m
         efh = (2*(1+m)/m*E(k)-F(k))/(3*m)
c         efhp = (2*(1+m)/m*E(k)-F(k))/(3*m) -2/(3*m**2)

         eminusone = E(k)-1
         if (m.lt.1.e-4) eminusone=m*(ec1+ec2*m+ec3*m**2)
     +                  -m*log(m)*(ed1+ed2*m+ed3*m**2)

         efhp=(2*(1+m)/m*eminusone-F(k)+2)/(3*m)  !efhp


         b=2*xx*xsa(k)
         c2=sumx**2 + xia(k)**2
         c=sqrt(c2)
         c3=c*c2
         c5=c2*c3
 
         ac=(1+m)/2
         cb=2/(1-m)

         sing = 8/(3*c5*m**2)
         i50= 4/c5*(efhp) 
         i51= 4/c5*(cb*(ac*efhp-eth)+4/(3*m*(1-m)) ) 
         i52=4/c5*(cb**2*(ac**2*efhp-2*ac*eth+F(k))+8/(3*m*(1-m)**2))
         i53= 4/c5*( cb**3*(ac**3*efhp-3*ac**2*eth+3*ac*F(k)-E(k)) 
     +        +4*(3+m**2)/(3*m*(1-m)**3) )

c         sing=0
c         i50= 4/c5*efh
c         i51= 4/(c3*b)*(ac*efh-eth) 
c         i52= 4/(c*b**2)*(ac**2*efh-2*ac*eth+F(k))
c         i53= 4*c/b**3*(ac**3*efh-3*ac**2*eth+3*ac*F(k)-E(k))

         q11=-6*xsa(k)*( xsa(k)**3*i51 - xsa(k)**2*xx*(i50+2*i52) +
     +       xsa(k)*xx**2*(i53+2*i51) - xx**3*i52 + (xsa(k)-xx)**3*sing)
                
         q12=-6*xsa(k)*xia(k)*((xsa(k)**2+xx**2)*i51-xsa(k)*xx*(i50+i52)
     +                    + (xsa(k)-xx)**2*sing)
         q21= q12
         q22= -6*xsa(k)*xia(k)**2*( xsa(k)*i51 -xx*i50+(xsa(k)-xx)*sing)

         qu(k)= (q11*yda(k)-q12*xda(k))*ua(k)
     +         +(q21*yda(k)-q22*xda(k))*va(k)


         q11=-6*xsa(k)*xia(k)*( xx**2*i52 +xsa(k)**2*i50-2*xsa(k)*xx*i51
     +                    + (xx-xsa(k))**2*sing)
         q12=-6*xsa(k)*xia(k)**2*( xsa(k)*i50 -xx*i51 +(xsa(k)-xx)*sing)
         q21= q12
         q22= -6*xsa(k)*xia(k)**3*(i50+sing)

         qv(k)=(q11*yda(k)-q12*xda(k))*ua(k)
     +        +(q21*yda(k)-q22*xda(k))*va(k)
      enddo
 
      trapbu=0
      trapbv=0
      do k = 1,ns-1
         trapbu = trapbu+ qu(k)
         trapbv = trapbv+ qv(k)
      enddo
      trapbu=trapbu*delalf
      trapbv=trapbv*delalf

c      xi0 = y0-yy
      rsq = xx**2+xi0**2
      r = sqrt(rsq)

      dql=-12*(pi*v0*xd0**2*xi0**2*xx)/r**5  !including c1,c2,d2,d2
      f1 = 2*( 3*xx**2*ud0*ydd0 -            
     +     xi0**2*(3*xd0*vdd0+4*xddd0*v0) )
      f2 = 6*xd0*xi0*(2*xx**2-3*xi0**2)/rsq  
     +               *(xd0*ud0-2*v0*ydd0)
      f3 = 15*xd0**3*xi0**2*v0*(-3*xx**2+4*xi0**2)/rsq**2 
      d3ql=6*pi*xx*xd0/r**5*(f1+f2+f3)
c      d3qln=(qu(3)-3*qu(2)+3*qu(1))/delalf**3
c      print*,'d3q=',d3ql,d3qln,(d3ql-d3qln)
c     this is correct finite difference but for n=100 values are 
c     completely off. large 3rd derivatives.

      dqr=(-qu(ns+5)+8*qu(ns+4)-8*qu(ns+2)+qu(ns+1))/(12*dfine)
      d3qr=(qu(ns+5)-2*qu(ns+4)+2*qu(ns+2)-qu(ns+1))/(2*dfine**3)
      trapbu= trapbu+qu(ns)*delalf/2
     +    - delalf**2*(dqr -dql )/12
     +    + delalf**4*(d3qr-d3ql)/720

      dql=12*(pi*v0*xd0**2*xi0**3)/r**5  !unchanged 
      f1 = 2*( -3*xx**2*ud0*ydd0 +
     +     xi0**2*(3*xd0*vdd0+4*xddd0*v0) )
      f2 = 6*xd0*xi0*(3*xx**2-2*xi0**2)/rsq
     +               *(-xd0*ud0+2*v0*ydd0)
      f3 = 15*xd0**3*xi0**2*v0*(5*xx**2-2*xi0**2)/rsq**2
      d3ql = 6*pi*xi0*xd0/r**5*(f1+f2+f3)
c      d3ql=(qv(3)-3*qv(2)+3*qv(1))/delalf**3
      dqr=(-qv(ns+5)+8*qv(ns+4)-8*qv(ns+2)+qv(ns+1))/(12*dfine)
      d3qr=(qv(ns+5)-2*qv(ns+4)+2*qv(ns+2)-qv(ns+1))/(2*dfine**3)
c      dqr=(-qv(ns+2)+8*qv(ns+1)-8*qv(ns-1)+qv(ns-2))/(12*delalf)
c      d3qr=(qv(ns+2)-2*qv(ns+1)+2*qv(ns-1)-qv(ns-2))/(2*delalf**3)
      trapbv= trapbv+qv(ns)*delalf/2
     +    - delalf**2*(dqr -dql )/12
     +    + delalf**4*(d3qr-d3ql)/720
 
      trapbu = trapbu/pi8
      trapbv = trapbv/pi8

      return
      end


      SUBROUTINE trapbqnok(k0,xx,yy,alf,jmax,xd0,y0,ud0,v0,vdd0,
     +  cc1,cc2,dd1,dd2,trapbu,trapbv)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute trapezoid approximation T_6 of B_q for qu,qv, to test
c     needed to get the error E[b]
c     uses ns tchebb pts on [0,alfend] plus 2nd and 4th order corrections
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,jmax,ns,k0
      MODE xx,yy,xd0,y0,ud0,v0,trapbu,trapbv,cc1,cc2,dd1,dd2,alf(0:*)
      MODE delalf,xi0,va(0:nmax),xsa(0:nmax),xia(0:nmax),xda(0:nmax),
     +    yda(0:nmax),ua(0:nmax),sumx,m,b,c,c2,c3,c5,eth,efh,efhp,ac,cb,
     +    i50,i51,i52,i53,q11,q12,q21,q22,r,dql,d3ql,dqr,d3qr,f1,f2,f3
      MODE ksq(0:nmax),F(0:nmax),E(0:nmax),qu(0:nmax),qv(0:nmax)
      MODE xddd0,ydd0,vdd0,alfa,dfine,sing,rsq,d3qln

      ns=jmax
      delalf = alf(jmax)/jmax  !=delalf
      dfine=delalf/100
c      print*,delalf
c      stop
c      va=v0        !0th order approximation
      xi0=y0-yy    !0th order approximation, basepoint
c      vdd0=0
      xddd0=6*cc1
      ydd0=2*dd1

c     compute b(xx,yy) where yy.ne.0 (point not on axis)

      do k=1,ns
        xsa(k)=xd0*alf(k) +   cc1*alf(k)**3  +   cc2*alf(k)**5
        xda(k)=xd0        + 3*cc1*alf(k)**2  + 5*cc2*alf(k)**4
        xia(k)=(y0-yy) +  dd1*alf(k)**2 +  dd2*alf(k)**4  
        yda(k)=         2*dd1*alf(k)    +4*dd2*alf(k)**3   
c        xia(k)=(y0-yy)+ydd0*alf(k)**2/2+ dd1*alf(k)**4+   dd2*alf(k)**6
c        yda(k)=        ydd0*alf(k)    +4*dd1*alf(k)**3+ 6*dd2*alf(k)**5
        ua(k) =ud0*alf(k)  
        va(k) =v0 +vdd0*alf(k)**2/2  

        sumx = xsa(k) + xx
        ksq(k) = 4*xsa(k)*xx/( xia(k)**2 + sumx**2)
      enddo
      do k=ns+1,ns+5
        alfa=alf(ns)+(k-ns-3)*dfine
        xsa(k)=xd0*alfa +   cc1*alfa**3  +   cc2*alfa**5
        xda(k)=xd0        + 3*cc1*alfa**2  + 5*cc2*alfa**4
        xia(k)=(y0-yy) +  dd1*alfa**2 +  dd2*alfa**4  
        yda(k)=         2*dd1*alfa    +4*dd2*alfa**3   
c        xia(k)=(y0-yy)+ydd0*alfa**2/2+ dd1*alfa**4+   dd2*alfa**6
c        yda(k)=        ydd0*alfa    +4*dd1*alfa**3+ 6*dd2*alfa**5
        ua(k) =ud0*alfa  
        va(k) =v0 +vdd0*alfa**2/2  

        sumx = xsa(k) + xx
        ksq(k) = 4*xsa(k)*xx/( xia(k)**2 + sumx**2)
      enddo
 
      call ellip(ksq(1),F(1),E(1),ns+4)
 
      do k = 1,ns+5
         sumx = xsa(k) + xx

         m =1-ksq(k)
         eth = E(k)/m
         efh = (2*(1+m)/m*E(k)-F(k))/(3*m)
         efhp = (2*(1+m)/m*E(k)-F(k))/(3*m) -2/(3*m**2)

         b=2*xx*xsa(k)
         c2=sumx**2 + xia(k)**2
         c=sqrt(c2)
         c3=c*c2
         c5=c2*c3
 
         ac=(1+m)/2
         cb=2/(1-m)

         sing = 8/(3*c5*m**2)
         i50= 4/c5*(efhp) 
         i51= 4/c5*(cb*(ac*efhp-eth)+4/(3*m*(1-m)) ) 
         i52=4/c5*(cb**2*(ac**2*efhp-2*ac*eth+F(k))+8/(3*m*(1-m)**2))
         i53= 4/c5*( cb**3*(ac**3*efhp-3*ac**2*eth+3*ac*F(k)-E(k)) 
     +        +4*(3+m**2)/(3*m*(1-m)**3) )

c         sing=0
c         i50= 4/c5*efh
c         i51= 4/(c3*b)*(ac*efh-eth) 
c         i52= 4/(c*b**2)*(ac**2*efh-2*ac*eth+F(k))
c         i53= 4*c/b**3*(ac**3*efh-3*ac**2*eth+3*ac*F(k)-E(k))

         q11=-6*xsa(k)*( xsa(k)**3*i51 - xsa(k)**2*xx*(i50+2*i52) +
     +       xsa(k)*xx**2*(i53+2*i51) - xx**3*i52 + (xsa(k)-xx)**3*sing)
                
         q12=-6*xsa(k)*xia(k)*((xsa(k)**2+xx**2)*i51-xsa(k)*xx*(i50+i52)
     +                    + (xsa(k)-xx)**2*sing)
         q21= q12
         q22= -6*xsa(k)*xia(k)**2*( xsa(k)*i51 -xx*i50+(xsa(k)-xx)*sing)

         qu(k)= (q11*yda(k)-q12*xda(k))*ua(k)
     +         +(q21*yda(k)-q22*xda(k))*va(k)


         q11=-6*xsa(k)*xia(k)*( xx**2*i52 +xsa(k)**2*i50-2*xsa(k)*xx*i51
     +                    + (xx-xsa(k))**2*sing)
         q12=-6*xsa(k)*xia(k)**2*( xsa(k)*i50 -xx*i51 +(xsa(k)-xx)*sing)
         q21= q12
         q22= -6*xsa(k)*xia(k)**3*(i50+sing)

         qv(k)=(q11*yda(k)-q12*xda(k))*ua(k)
     +        +(q21*yda(k)-q22*xda(k))*va(k)
      enddo
 
      trapbu=0
      trapbv=0
      do k = 1,k0-1
         trapbu = trapbu+ qu(k)
         trapbv = trapbv+ qv(k)
      enddo
      do k = k0+1,ns-1
         trapbu = trapbu+ qu(k)
         trapbv = trapbv+ qv(k)
      enddo
      trapbu=trapbu*delalf
      trapbv=trapbv*delalf

c      xi0 = y0-yy
      rsq = xx**2+xi0**2
      r = sqrt(rsq)

      dql=-12*(pi*v0*xd0**2*xi0**2*xx)/r**5  !including c1,c2,d2,d2
      f1 = 2*( 3*xx**2*ud0*ydd0 -            
     +     xi0**2*(3*xd0*vdd0+4*xddd0*v0) )
      f2 = 6*xd0*xi0*(2*xx**2-3*xi0**2)/rsq  
     +               *(xd0*ud0-2*v0*ydd0)
      f3 = 15*xd0**3*xi0**2*v0*(-3*xx**2+4*xi0**2)/rsq**2 
      d3ql=6*pi*xx*xd0/r**5*(f1+f2+f3)
c      d3qln=(qu(3)-3*qu(2)+3*qu(1))/delalf**3
c      print*,'d3q=',d3ql,d3qln,(d3ql-d3qln)
c     this is correct finite difference but for n=100 values are 
c     completely off. large 3rd derivatives.

      dqr=(-qu(ns+5)+8*qu(ns+4)-8*qu(ns+2)+qu(ns+1))/(12*dfine)
      d3qr=(qu(ns+5)-2*qu(ns+4)+2*qu(ns+2)-qu(ns+1))/(2*dfine**3)
      trapbu= trapbu+qu(ns)*delalf/2
     +    - delalf**2*(dqr -dql )/12
     +    + delalf**4*(d3qr-d3ql)/720

      dql=12*(pi*v0*xd0**2*xi0**3)/r**5  !unchanged 
      f1 = 2*( -3*xx**2*ud0*ydd0 +
     +     xi0**2*(3*xd0*vdd0+4*xddd0*v0) )
      f2 = 6*xd0*xi0*(3*xx**2-2*xi0**2)/rsq
     +               *(-xd0*ud0+2*v0*ydd0)
      f3 = 15*xd0**3*xi0**2*v0*(5*xx**2-2*xi0**2)/rsq**2
      d3ql = 6*pi*xi0*xd0/r**5*(f1+f2+f3)
c      d3ql=(qv(3)-3*qv(2)+3*qv(1))/delalf**3
      dqr=(-qv(ns+5)+8*qv(ns+4)-8*qv(ns+2)+qv(ns+1))/(12*dfine)
      d3qr=(qv(ns+5)-2*qv(ns+4)+2*qv(ns+2)-qv(ns+1))/(2*dfine**3)
c      dqr=(-qv(ns+2)+8*qv(ns+1)-8*qv(ns-1)+qv(ns-2))/(12*delalf)
c      d3qr=(qv(ns+2)-2*qv(ns+1)+2*qv(ns-1)-qv(ns-2))/(2*delalf**3)
      trapbv= trapbv+qv(ns)*delalf/2
     +    - delalf**2*(dqr -dql )/12
     +    + delalf**4*(d3qr-d3ql)/720
 
      trapbu = trapbu/pi8
      trapbv = trapbv/pi8

      return
      end



c      xi = yn-yy
c      rsq = xx**2+xi**2
c      r = dsqrt(rsq)
c      fact = 12*pi*vn*xdn**2/r**5
c      gpur = -fact*xi**2*xx
c      gpvr =  fact*xi**3

C     FOR BOTH U AND V COMPONENTS
c      f1 = 2*( 3*xx**2*udn*yddn -
c     +     xi**2*(3*xdn*vddn+4*xdddn*vn) )
c      f2 = 6*xdn*xi*(2*xx**2-3*xi**2)/rsq
c     +               *(xdn*udn-2*vn*yddn)
c      f3 = 15*xdn**3*xi**2*vn*(-3*xx**2+4*xi**2)/rsq**2
c      gpppur = 6*pi*xx*xdn/r**5*(f1+f2+f3)
c      f1 = 2*( -3*xx**2*udn*yddn +
c     +     xi**2*(3*xdn*vddn+4*xdddn*vn) )
c      f2 = 6*xdn*xi*(3*xx**2-2*xi**2)/rsq
c     +               *(-xdn*udn+2*vn*yddn)
c      f3 = 15*xdn**3*xi**2*vn*(5*xx**2-2*xi**2)/rsq**2
c      gpppvr = 6*pi*xi*xdn/r**5*(f1+f2+f3)



      SUBROUTINE compBq(xx,yy,alf,ns,
     +                  xd0,y0,ud0,v0,vd20,
c     +                  bintu,bintv,btrapu,btrapv)
     +     corrbu,corrbv)
c   later: corrbu=(bintu) - (btrapu)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computes E[Bbet] = intB-trapB  (no correction)
c Input:
c    xx,yy  target point
c    alf,ns parametrization of xs,ys
c    xd0,y0,ud0,v0,vd20: values on axis used to obtain quadratic
c                        or cubic approximations of xs,ys,u,v
c    alfb,xb,yb,xdb,ydb: values at projection=base point
c Output:
c   for now, corrected intbu =intbu+corrintbu
c            corrected trapbu=trapbu+corrtrapbu
c   later: corrbu=(intbu+corrintbu) - (trapbu+corrtrapbu)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns
      MODE xx,yy,alf(0:*),
     +     xd0,y0,ud0,v0,vd20
c OUTPUT
      MODE corrbu,corrbv
c LOCAL
      integer jmax,k0
      MODE bintu,bintv,btrapu,btrapv
      MODE alfb,xb,yb,xdb,ydb
      MODE alfend,cc1,cc2,dd1,dd2
      MODE xk,xik,uxk,vxk,uyk,vyk

      call findalfp(alfb,xx,yy,xb,yb,xdb,ydb)
      call intbq(xx,yy,ns,xd0,y0,ud0,v0,vd20,bintu,bintv,
     +           alfend,alfb,xb,yb,xdb,ydb,cc1,cc2,dd1,dd2,jmax,
     +           k0,xk,xik,uxk,vxk,uyk,vyk)

      call trapbq(xx,yy,alf,jmax,xd0,y0,ud0,v0,vd20,
     +     cc1,cc2,dd1,dd2,btrapu,btrapv)
      corrbu=bintu-btrapu
      corrbv=bintv-btrapv
      return
      end


      SUBROUTINE findalfp(alfp,xx,yy,xproj,yproj,xdp,ydp)
c find projection on circle, use for all n to test convergence
      implicit none
c INPUT
      MODE xx,yy
c OUTPUT
      MODE alfp,xproj,yproj,xdp,ydp
c LOCAL
      MODE rad

      rad=0.95q0
      xproj = xx*rad/sqrt(xx**2+(yy+4)**2)
      yproj = (yy+4)*rad/sqrt(xx**2+(yy+4)**2)-4
      alfp=asin(xproj/rad)
      xdp=rad*cos(alfp)
      ydp=rad*sin(alfp)
      return
      end

