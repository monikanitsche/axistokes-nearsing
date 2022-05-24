c whenever there is x/rho2, expand together, to get better accuracy near
c axis. so expand numerators/x. to third order only

c Contains
c SUBROUTINES 
c     correctionqatxnew: given all coefficients at basepoint,
c                        outputs corrections for intQu intQv
c     corrqu11: all these output corrections and local term 
c     corrqu12: for the corresponding components
c     corrqu22:
c     corrqv11:
c     corrqv12:
c     corrqv22:
c FUNCTIONS :These output the correction for the 3 terms in each component
c      corrqrho4
c      corrqrho2
c      corrqlog
c SUBROUTINES 
c     compcoeffnew: computes coefficients of some factors used 
c                   by all expansions, such as
c                   e3,e4,e5,cm,1/r2 1/r2^3 1/r2^5, 1/r2x..,, all "small")
#include "flags.h"
      SUBROUTINE correctionqatxnew(alf,alfb,ns,d,k0,xk,xik,
     +      xb,x0,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +      ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +      uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3,
     +      uxk,vxk,uyk,vyk,
     +      corru,corrv,largecorr,restu,restv,hk,termout)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    alfb : basepoint
c    xb,xd,xdd,xddd,.. values AT BASEPOINT alfb
c    x0,y0: target point (y0 not needed, in dy=yb-y0)
c
c Use expansions in bubbles/math/math2022:
c        expqu11redo2v2.pdf
c        expqu12redo2v2.pdf
c        expqu22redo2v2.pdf
c        expqv11redo2v2.pdf
c        expqv12redo2v2.pdf
c        expqv22redo2v2.pdf
c to find corrections for Qu11*u, Qu12*v, Qu12*u, Qu22*v
c
c Note: x=xb+c1*alf+c2*alf^2+c3*alf^3+c4*alf^4
c       xi=dy+d1*alf+d2*alf^2+d3*alf^3+d4*alf^4
c       f=f0+f1*alf+f2*alf^2+f3*alf^3
c dx=xb-x0

c Note2: Q11*u =mnum*f/r4+ mnum*f/r2 * 1/r1^2  + lognum*f/r2^3 * log(r1^2)
c
c        r1^2=|x-x0|^2 = (x-x0)^2+xi^2, xi=y-y0  dy=yb-y0
c        r2^2=         = (x+x0)^2+xi^2           dx=xb-x0
c
c        rho^2=d^2+calf^2
c
c Throughout use: x0>0
c
c corrnok intH - sumHk, k.ne.k0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
C     INPUT
      integer ns,k0
      MODE alf(0:*),alfb,d,xk,xik,
     +      xb,x0,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +      ux0,ux1,ux2,ux3,vx0,vx1,vx2,vx3,
     +      uy0,uy1,uy2,uy3,vy0,vy1,vy2,vy3,
     +      uxk,vxk,uyk,vyk
      logical largecorr
C     OUTPUT
      MODE corru,corrv,restu,restv,hk,termout
C3     LOCAL
      MODE corrs(nb),trapsing(nb),hk2,cm
      MODE corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8,
     +       rest1,rest2,rest3,rest4,rest5,rest6,rest7,rest8

      call compcoeffnew(xb,x0,d,c1,c2,c3,c4,dy,d1,d2,d3,d4,xk,xik,cm)
      call compallcorrnew(alf,alfb,ns,d,cm,corrs,k0,largecorr,trapsing)

      call corrqu11(x0,uy0,uy1,uy2,uy3,corrs,
     +       largecorr,xk,xik,uyk,trapsing,corr1,rest1,hk2) 
      hk=hk2
      termout=rest1
c      print*,rest1*pi/800,'rest1'
      call corrqu12(x0,ux0,ux1,ux2,ux3,corrs,
     +       largecorr,xk,xik,uxk,trapsing,corr2,rest2,hk2) 
      call corrqu12(x0,vy0,vy1,vy2,vy3,corrs,
     +       largecorr,xk,xik,vyk,trapsing,corr3,rest3,hk2) 
      call corrqu22(x0,vx0,vx1,vx2,vx3,corrs,
     +       largecorr,xk,xik,vxk,trapsing,corr4,rest4,hk2) 

      call corrqv11(x0,uy0,uy1,uy2,uy3,corrs,
     +       largecorr,xk,xik,uyk,trapsing,corr5,rest5,hk2) 
c      hk=hk2
c      termout=rest5
      call corrqv12(x0,ux0,ux1,ux2,ux3,corrs,
     +       largecorr,xk,xik,uxk,trapsing,corr6,rest6,hk2) 
c      hk=-hk2
c      termout=-rest6
      call corrqv12(x0,vy0,vy1,vy2,vy3,corrs,
     +       largecorr,xk,xik,vyk,trapsing,corr7,rest7,hk2) 
c      hk=hk2
c      termout=rest7
      call corrqv22(x0,vx0,vx1,vx2,vx3,corrs,
     +       largecorr,xk,xik,vxk,trapsing,corr8,rest8,hk2) 
c      hk=-hk2
c      termout=-rest8

      corru=(corr1 -corr2 +corr3 -corr4)/pi8
      corrv=(corr5 -corr6 +corr7 -corr8)/pi8

      restu=(rest1 -rest2 +rest3 -rest4)/pi8
      restv=(rest5 -rest6 +rest7 -rest8)/pi8

      return
      end


      SUBROUTINE corrqu11(x0,f0,f1,f2,f3,corrs,
     +       largecorr,xk,xik,f0k,trapsing,corr,rest,hk) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  h^3 method: numer6 is not used by corrqrho
c  have x to 3rd order, w to 2nd, all small products to 5th
c
c  1/rho4: should call prod2(fact,irx,q), prod255(q,a3,piece1)
c  1/rho2: should call prod1(fact,ir3x,q), prod133(q,a3,piece1)
c  logrho: should call prod1(fact,ir5x,q)
c
c  right now, since we have more derivatives, use
c  1/rho4: prod3(fact,irx,q), prod366(q,a3,piece1)
c  1/rho2: prod2(fact,ir3x,q), prod244(q,a3,piece1)
c  logrho: prod2(fact,ir5x,q)
c
c  which are the correct calls for h^4 method.
c  with d4=0, prod446 same as prod335, except that 6th comp .ne. 0
c
c  using prod3(fact,irx,q) instead of prod2(fact,irx,q)
c  seems to cause larger errors near axis. 
c  but lower errors away from axis CHECK THIS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*),xk,xik,f0k,trapsing(*)
      logical largecorr
c OUTPUT 
      MODE corr,rest,hk
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),q(0:6),small(0:6),
     +       piece1(0:6),piece2(0:6),piece3(0:6),piece4(0:6),
     +       corrqrho4,corrqrho2,corrqlog,corr1,corr2,corr3
      MODE hk1,hk2,hk3,piece5(0:6)
      MODE denom,Fk,reglog,reglog1,rest0num,rest1num,rest2num
      
c 1/RHO4 --------------------------------------------------
      lead = -2/x0

c Piece 1: small=delr^3=a3,   fact=(x+x0)^3
      do j=0,4
        fact(j)=x3(j) + 3*x2(j)*x0 + 3*x1(j)*x02+ones(j)*x03
      enddo
      call prodrho4(fact,irx,a3,piece1)

c Piece 2: small=-delr2*xi2=a2b2,  fact=(x+x0)^2
      do j=0,4
        fact(j)=-(x2(j) + 2*x1(j)*x0 + ones(j)*x02)
      enddo
      call prodrho4(fact,irx,a2b2,piece2)

c Piece 3: small=delr*xi^4=ab4,  fact=-(x+x0)
      do j=0,4
        fact(j) = -(x1(j) + ones(j)*x0 )
      enddo
      call prodrho4(fact,irx,ab4,piece3)

c Piece 4: small=xi^6=b6,  fact=1
      call prodrho4(ones,irx,b6,piece4)

c Piece 1 + Piece2
      do j=0,6
        numer(j)=piece1(j)+piece2(j)+piece3(j)+piece4(j)
      enddo
      corr1=corrqrho4(lead,numer,f0,f1,f2,f3,e3,e4,e5,corrs) 
      hk1=corrqrho4(lead,numer,f0,f1,f2,f3,e3,e4,e5,trapsing) 

c 1/RHO2 --------------------------------------------------
      lead = -3/(2*x0)

c Piece 1: small=delr^2=a2, fact=3*x^4+10*x^3*x0+12*x^2*x0^2+6*x*x0^3+x0^4
      do j=0,3
        fact(j)=3*x4(j)+10*x3(j)*x0+12*x2(j)*x02+6*x1(j)*x03+ones(j)*x04
      enddo
      call prodrho2(fact,ir3x,a2,piece1)

c Piece 2: small=delr*xi^2=ab2,  fact=
      do j=0,3
        fact(j)=5*x3(j)+13*x2(j)*x0+11*x1(j)*x02+3*ones(j)*x03
      enddo
      call prodrho2(fact,ir3x,ab2,piece2)

c Piece 3: small=xi^4=b4,  fact=
      do j=0,3
        fact(j) = -3*(x2(j)+4*x1(j)*x0+3*ones(j)*x02)
      enddo
      call prodrho2(fact,ir3x,b4,piece3)

c Piece 4: small=-5xi^6=-5b6,  fact=1
      do j=0,4
         small(j) = -5*b6(j)
      enddo
      call prodrho2(ones,ir3x,small,piece4)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)+piece2(j)+piece3(j)+piece4(j)
      enddo

      corr2=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      hk2=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,trapsing)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = -3/(16*x0)

c Piece 1: small=delr=a1, fact=..
      do j=0,2
        fact(j)=x5(j)+25*x4(j)*x0+78*x3(j)*x02+94*x2(j)*x03
     +          +49*x1(j)*x04+9*ones(j)*x05
      enddo
      call prodlog(fact,ir5x,a1,piece1)

c Piece 2: small=xi^2=b2,  fact=
      do j=0,2
        fact(j)=-(17*x4(j)+80*x3(j)*x0+142*x2(j)*x02
     +            +112*x1(j)*x03+33*ones(j)*x04)
      enddo
      call prodlog(fact,ir5x,b2,piece2)

c Piece 3: small=xi^4=b4,  fact=-3*(...)
      do j=0,2
        fact(j) = -3*(11*x2(j)+24*x1(j)*x0+13*ones(j)*x02)
      enddo
      call prodlog(fact,ir5x,b4,piece3)

c Piece 4: small=-15xi^6=-15b6,  fact=1
      do j=0,2
         small(j) = -15*b6(j)
      enddo
      call prodlog(ones,ir5x,small,piece4)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 1 (skip xi^2)
      fact(0)=-3*(x1(0)+x0)**3*(13*x3(0)+77*x2(0)*x0+31*x1(0)*x02+7*x03)
      denom=64*x1(0)*x0*r2**7
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces 
      do j=0,2
        numer(j)=piece1(j)+piece2(j)+piece3(j)+piece4(j)+piece5(j)/lead
      enddo
      corr3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      hk3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,trapsing)
      
      corr=corr1+corr2+corr3
      hk=hk1+hk2+hk3

      rest=0
#ifdef ROUNDOFFFIX
c REST (value of Gk-Hk at alf(k0)) ---------------------------------------
c G-H is found in mathematica (H approx F2/rho1^4 + F1/rho1^2 +F0*log(rho1sq))
c it has 2 components
c     logterm =-F0*log(rho2sq)  
c     rest = rest0 + rest1*m*logm)
c     
c     drop high order xik terms. Careful: xik=O(h), not O(d). Also m=O(h^2)
c     drop m^2logm = O(h^4logh) terms, which are mult by h. Thus O(h^5logh)
c     BELOW: logr2k=log(r2ksq)

      if (largecorr) then

        lead =-3
        denom= 16*xk*x0*r2k5

        Fk = xk6+24*xk5*x0+xk4*(53*x02-17*xik2)-
     -       3*(x02+xik2)**2*(3*x02+5*xik2)+16*xk3*(x03-5*x0*xik2)-
     -       xk2*(45*x04+142*x02*xik2+33*xik4)-
     -       8*xk*(5*x05+14*x03*xik2+9*x0*xik4)
        reglog=-lead*Fk/denom*logr2k

        Fk=-3*(xk+x0)**3*(13*xk3+77*xk2*x0+31*xk*x02+7*x03)
        reglog1=-Fk/(4*denom)*m*logr2k

        rest0num=3*(xk+x0)**3*((47+8*log2)*xk3+(67+168*log2)*xk2*x0+
     -                     (109-104*log2)*xk*x02+3*(11-24*log2)*x03)
     -           -3*xik2*(xk+x0)**2*((-113+136*log2)*xk2
     +                   +2*(-95+184*log2)*xk*x0+ 3*(-27+88*log2)*x02)
     +           -9*xik4*(xk+x0)*((-27+88*log2)*xk+(-21+104*log2)*x0)
c dropped xik6 term O(h^7)
        rest1num=(xk+x0)**3*((-47+312*log2)*xk3+3*(-221+616*log2)*xk2*x0
     -                   +3*(-23+248*log2)*xk*x02+(11+168*log2)*x03)
     +           +xik2*(xk+x0)**2*((-25+648*log2-162*logm)*xk2
     -                   +2*(-11+600*log2-150*logm)*xk*x0
     +                   +(47+456*log2-114*logm)*x02)
c dropped xik4*m*logm term O(h^7*logh)

c        rest2num=3*(xk+x0)**3*((-833+1776*log2-444*logm)*xk3+
c     -           3*(-3519+6160*log2-1540*logm)*xk2*x0+
c     -           3*(-577+1264*log2-316*logm)*xk*x02+
c     -           (-191+528*log2-132*logm)*x03)
c dropped xik2*m2*logm term O(h^7*logh)

        rest=reglog + reglog1
     +              + rest0num/(2*denom) 
     +              + m*rest1num/(8*denom)
c     +              + m**2*rest2num/(512*denom)  

c Note: rest is necessarily O(h^4) since F2/r4, F1/r2, F0logm,mlogm
c are only approximated to fourth order. so the formula in mathematica
c used to find rest is only 4th order.

c multiply by density at xk
        rest=f0k*rest
      endif
#endif
      return
      end
      

      SUBROUTINE corrqu12(x0,f0,f1,f2,f3,corrs,
     +       largecorr,xk,xik,f0k,trapsing,corr,rest,hk) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*),xk,xik,f0k,trapsing(*)
      logical largecorr
c OUTPUT 
      MODE corr,rest,hk
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),small(0:6),
     +       piece1(0:6),piece2(0:6),piece3(0:6),
     +       corrqrho4,corrqrho2,corrqlog,corr1,corr2,corr3
      MODE hk1,hk2,hk3,piece5(0:6)
      MODE denom,Fk,reglog,reglog1,rest0num,rest1num
      
c 1/RHO4 --------------------------------------------------
      lead = -4/x0

c Piece 1: small=delr^2*xi=a2b,   fact=(x+x0)^2
      do j=0,4
        fact(j)=x2(j) + 2*x1(j)*x0 + x02*ones(j)  
      enddo
      call prodrho4(fact,ir,a2b,piece1)

c Piece 2: small=-xi^5,  fact=1
      do j=0,6
        small(j)=-b5(j)
      enddo
      call prodrho4(ones,ir,small,piece2)

c Sum Pieces
      do j=0,6
        numer(j)=piece1(j)+piece2(j)
      enddo
      corr1=corrqrho4(lead,numer,f0,f1,f2,f3,e3,e4,e5,corrs) 
c      hk1=corrqrho4(lead,num0,num1,num2,num3,num4,num5,num6,
c     +                f0,f1,f2,f3,e3,e4,e5,trapsing) 

c 1/RHO2 --------------------------------------------------
      lead = -3/x0

c Piece 1: fact=(x+x0)^2,  small=delr^2*xi+4*xi^3,  
c Use that fact is the same as previous
      do j=0,4
         small(j) = a2b(j) + 4*b3(j)
      enddo
      call prodrho2(fact,ir3,small,piece1)

c Piece 2: small=3xi^5=b5,  fact=1
      do j=0,4
         small(j) = 3*b5(j)
      enddo
      call prodrho2(ones,ir3,small,piece2)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)+piece2(j)
      enddo

      corr2=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
c      hk2=corrqrho2(lead,num0,num1,num2,num3,num4,
c     +       ir30,ir31,ir32,ir33,f0,f1,f2,f3,e3,e4,trapsing)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = -3/(8*x0)
c Piece 1: small=xi,  fact=5*x^4+32*x^3*x0+54*x^2*x0^2+32*x*x0^3+5*x0^4
      do j=0,2
         fact(j) =  5*x4(j) + 32*x3(j)*x0 + 54*x2(j)*x02 
     +           + 32*x1(j)*x03 +5*x04*ones(j)
      enddo
      call prodlog(fact,ir5,b1,piece1)

c Piece 2: small=8xi^3  fact=(x+x0)^2=x^2+2*x*x0+x0^2
      do j=0,2
         fact(j) =  8*(x2(j) + 2*x1(j)*x0 + ones(j)*x02)
      enddo
      call prodlog(fact,ir5,b3,piece2)

c Piece 3: 3xi^5
      do j=0,2
         small(j) = 3*b5(j)
      enddo
      call prodlog(ones,ir5,small,piece3)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order xi (skip xi^3)
      fact(0)=-3*b1(0)*(x1(0)+x0)**2*(3*x2(0)+26*x1(0)*x0+3*x02)
      denom=32*x0*r2**7
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces
      do j=0,2
        numer(j)=piece1(j)+piece2(j)+piece3(j)+piece5(j)/lead
      enddo

      corr3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
c      hk3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,trapsing)
      
      corr=corr1+corr2+corr3

      rest=0
#ifdef ROUNDOFFFIX
c REST (value of Gk-Hk at alf(k0)) ---------------------------------------
c G-H is found in mathematica (H approx F2/rho1^4 + F1/rho1^2 +F0*log(rho1sq))
c it has 2 components
c     logterm =-F0*log(rho2sq)  
c     rest = rest0 + m(rest11 + rest12*logm)
c     
c     drop high order xik terms. Careful: xik=O(h), not O(d). Also m=O(h^2)
c     drop m^2logm = O(h^4logh) terms, which are mult by h. Thus O(h^5logh)
c     BELOW: logr2k=log(r2ksq)

      if (largecorr) then

        lead =-3
        denom= 8*x0*r2k5

        Fk = xik*(5*xk4+32*xk3*x0+5*x04+8*x02*xik2+3*xik4+
     -            16*xk*x0*(2*x02+xik2)+xk2*(54*x02+8*xik2))
        reglog=-lead*Fk/denom*logr2k

        Fk= -3*xik*(xk+x0)**2*(3*xk2+26*xk*x0+3*x02)
        reglog1=-Fk/(4*denom)*m*logr2k

        rest0num=3*xik*(xk+x0)**2*((7+40*log2)*xk2
     +                  +2*(-7+88*log2)*xk*x0+(7+40*log2)*x02)
     +            +xik3*48*(1+4*log2)*(xk+x0)**2
c dropped xik5 O(h^6)

        rest1num= xik*(xk+x0)**2*((-29+72*log2-18*logm)*xk2+
     -       2*(-163+312*log2-78*logm)*xk*x0+(-29+72*log2-18*logm)*x02)
c dropped xik3*mlogm O(h^6logh)

        rest=reglog+reglog1 + rest0num/(2*denom) + m*rest1num/(8*denom)

c multiply by density at xk
        rest=f0k*rest
      endif
#endif
      return
      end
      

      SUBROUTINE corrqu22(x0,f0,f1,f2,f3,corrs,
     +       largecorr,xk,xik,f0k,trapsing,corr,rest,hk) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*),xk,xik,f0k,trapsing(*)
      logical largecorr
c OUTPUT 
      MODE corr,rest,hk
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),q(0:6),p(0:6),
     +       piece1(0:6),piece2(0:6),piece3(0:6),
     +       corrqrho4,corrqrho2,corrqlog,corr1,corr2,corr3
      MODE hk1,hk2,hk3,piece5(0:6)
      MODE denom,Fk,reglog,reglog1,rest0num,rest1num
      
c 1/RHO4 --------------------------------------------------
      lead = -8/x0

c Piece 1: small=delr*xi^2=ab2,  fact=x(x+x0)=x^2+x*x0
      do j=0,4
        fact(j)=x2(j) + x1(j)*x0 
      enddo
      call prodrho4(fact,ir,ab2,piece1)

c Piece 2: small=xi^4,  fact=x
      call prodrho4(x1,ir,b4,piece2)

c Sum Pieces
      do j=0,6
        numer(j)=piece1(j)+piece2(j)
      enddo
      corr1=corrqrho4(lead,numer,f0,f1,f2,f3,e3,e4,e5,corrs) 
c      hk1=corrqrho4(lead,num0,num1,num2,num3,num4,num5,num6,
c     +                f0,f1,f2,f3,e3,e4,e5,trapsing) 

c 1/RHO2 --------------------------------------------------
      lead = 6/x0

c Piece 1: small=xi^2=b2,  fact=x^3+4*x^2*x0+3*x*x0^2
      do j=0,3
         fact(j) = x3(j) + 4*x2(j)*x0 + 3*x1(j)*x02
      enddo
      call prodrho2(fact,ir3,b2,piece1)

c Piece 2: small=xi^4=b4,  fact=x
      call prodrho2(x1,ir3,b4,piece2)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j) +piece2(j)
      enddo

      corr2=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
c      hk2=corrqrho2(lead,num0,num1,num2,num3,num4,
c     +       ir30,ir31,ir32,ir33,f0,f1,f2,f3,e3,e4,trapsing)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = -3/(4*x0)
c piece 1: small=xi^2=b2,  fact=x^3+8x^2*x0+7*x*x0^2
      do j=0,2
         fact(j) =  x3(j) + 8*x2(j)*x0 + 7*x1(j)*x02 
      enddo
      call prodlog(fact,ir5,b2,piece1)

c piece 2: small=xi^4=b4,  fact=x=x1
      call prodlog(x1,ir5,b4,piece2)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 2 (skip xi^4)
      fact(0)= -3*b2(0)*x1(0)*(x2(0)+12*x1(0)*x0+11*x02)
      denom=16*x0*r2**7
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom
c Sum Pieces
      do j=0,2
        numer(j)=piece1(j) +piece2(j)  +piece5(j)/lead
      enddo

      corr3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
c      hk3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,trapsing)
      
      corr=corr1+corr2+corr3

      rest=0
#ifdef ROUNDOFFFIX
c REST (value of Gk-Hk at alf(k0)) ---------------------------------------
c G-H is found in mathematica (H approx F2/rho1^4 + F1/rho1^2 +F0*log(rho1sq))
c it has 2 components
c     logterm =-F0*log(rho2sq)  
c     rest = rest0 + m(rest11 + rest12*logm)
c     
c     drop high order xik terms. Careful: xik=O(h), not O(d). Also m=O(h^2)
c     drop m^2logm = O(h^4logh) terms, which are mult by h. Thus O(h^5logh)
c     BELOW: logr2term=log(r2ksq)/r2k**5  

      if (largecorr) then
        lead =-3

        denom= 4*x0*r2k5

        Fk = xk*xik2*(xk2+8*xk*x0+7*x02+xik2)
        reglog=-lead*Fk/denom*logr2k

        Fk= -3*xik2*xk*(xk2+12*xk*x0+11*x02)
        reglog1=-Fk/(4*denom)*m*logr2k

        rest0num=xik2*3*xk*(xk+x0)*((-1+8*log2)*xk+(-15+56*log2)*x0)
     +          +xik4*3*(-1+8*log2)*xk
        rest1num=xik2*xk*(xk+x0)*( (-11+24*log2-6*logm)*xk
     +                            +(-145+264*log2-66*logm)*x0)
c dropped xik4*mlogm O(h^7logh)

        rest=reglog+reglog1 + rest0num/(2*denom) + m*rest1num/(8*denom)

c multiply by density at xk
        rest=f0k*rest
      endif
#endif
      return
      end

      
      SUBROUTINE corrqv11(x0,f0,f1,f2,f3,corrs,
     +       largecorr,xk,xik,f0k,trapsing,corr,rest,hk) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*),xk,xik,f0k,trapsing(*)
      logical largecorr
c OUTPUT 
      MODE corr,rest,hk
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),q(0:6),small(0:6),
     +       piece1(0:6),piece2(0:6),piece3(0:6),piece4(0:6),
     +       corrqrho4,corrqrho2,corrqlog,corr1,corr2,corr3
      MODE hk1,hk2,hk3,piece5(0:6)
      MODE denom,Fk,reglog,reglog1,rest0num,rest1num
      
c 1/RHO4 --------------------------------------------------
      lead = -4

c Piece 1: small=delr^2*xi=a2b,   fact=(x+x0)^2
      do j=0,4
        fact(j)= x2(j) + 2*x1(j)*x0 + ones(j)*x02
      enddo
      call prodrho4(fact,irx,a2b,piece1)
c      call prod3(fact,irx,q)              !q=fact*ir
c      call prod366(q,a2b,piece1)         !piece1=fact*ir*a2b
c      call prod4(fact,irx,q)              !q=fact*ir
c      call prod466(q,a2b,piece1)         !piece1=fact*ir*a2b

c Piece 2: small=delr*xi^3=ab3,  fact=-2(x+x0)
      do j=0,4
        fact(j) = -2*(x1(j) +  ones(j)*x0)
      enddo
      call prodrho4(fact,irx,ab3,piece2)

c Piece 3: small=xi^5=b5,  fact=1
      call prodrho4(ones,irx,b5,piece3)

c Sum Pieces 
      do j=0,6
        numer(j)=piece1(j)+piece2(j)+piece3(j)
      enddo
      corr1=corrqrho4(lead,numer,f0,f1,f2,f3,e3,e4,e5,corrs) 
c      hk1=corrqrho4(lead,num0,num1,num2,num3,num4,num5,num6,
c     +                f0,f1,f2,f3,e3,e4,e5,trapsing) 

c 1/RHO2 --------------------------------------------------
      lead = 3

c Piece 1: small=xi=b1,  fact=-(...)
      do j=0,3
        fact(j) = -( 5*x4(j)+8*x3(j)*x0-2*x2(j)*x02
     +              -8*x1(j)*x03-3*ones(j)*x04)
      enddo
      call prodrho2(fact,ir3x,b1,piece1)

c Piece 2: small=xi^3=b3,  fact=2(...)
      do j=0,3
        fact(j) = 2*(x2(j)+4*x1(j)*x0+3*ones(j)*x02)
      enddo
      call prodrho2(fact,ir3x,b3,piece2)

c Piece 3: small=3xi^5=3b5,  fact=1
      do j=0,4
         small(j) = 3*b5(j)
      enddo
      call prodrho2(ones,ir3x,small,piece3)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)+piece2(j)+piece3(j)
      enddo

      corr2=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
c      hk2=corrqrho2(lead,num0,num1,num2,num3,num4,
c     +       ir30,ir31,ir32,ir33,f0,f1,f2,f3,e3,e4,trapsing)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = 3.q0/8

c Piece 1: small=xi=b1,  fact=
      do j=0,2
        fact(j)=  19*x4(j)+48*x3(j)*x0+42*x2(j)*x02
     +           +16*x1(j)*x03+3*ones(j)*x04
      enddo
      call prodlog(fact,ir5x,b1,piece1)
c      call prod2(fact,ir5x,q)              !q=fact*ir5x
c      call prod2(q,b1,piece1)      !piece2=fact*ir5x*small

c Piece 2: small=xi^3=b3,  fact=2*(...)
      do j=0,2
        fact(j) = 2*(5*x2(j)+8*x1(j)*x0+3*ones(j)*x02)
      enddo
      call prodlog(fact,ir5x,b3,piece2)

c Piece 3: small=3xi^5=3b5,  fact=1
      do j=0,2
         small(j) = 3*b5(j)
      enddo
      call prodlog(ones,ir5x,small,piece3)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 1 (skip xi^2)
      fact(0)=b1(0)*3*(x1(0)+x0)**2*(25*x2(0)+6*x1(0)*x0+x02)
      denom=32*x1(0)*r2**7
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces 
      do j=0,2
        numer(j)=piece1(j)+piece2(j)+piece3(j)  +piece5(j)/lead
      enddo
      corr3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      hk3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,trapsing)
      
      corr=corr1+corr2+corr3

      rest=0
#ifdef ROUNDOFFFIX
c REST (value of Gk-Hk at alf(k0)) ---------------------------------------
c G-H is found in mathematica (H approx F2/rho1^4 + F1/rho1^2 +F0*log(rho1sq))
c it has 2 components
c     logterm =-F0*log(rho2sq)/r2^5  
c     rest = rest0 + m(rest11 + rest12*logm)
c     
c     drop high order xik terms. Careful: xik=O(h), not O(d). Also m=O(h^2)
c     drop m^2logm = O(h^4logh) terms, which are mult by h. Thus O(h^5logh)
c     BELOW: logr2k=log(r2ksq)

      if (largecorr) then
        lead =3
        denom= 8*xk*r2k5

        Fk = xik*(19*xk4+48*xk3*x0+16*xk*x0*(x02+xik2)+3*(x02+xik2)**2+
     -        2*xk2*(21*x02+5*xik2))
        reglog=-lead*Fk/denom*logr2k

        Fk= 3*xik*(xk+x0)**2*(25*xk2+6*xk*x0+x02)
        reglog1=-Fk/(4*denom)*m*logr2k

        rest0num= -3*xik*(xk+x0)**2*((-23+152*log2)*xk2
     +                           +2*(7+40*log2)*xk*x0+3*(3+8*log2)*x02)
     +            -6*xik3*(xk+x0)*((7+40*log2)*xk+3*(3+8*log2)*x0)
c dropped xik5 O(h^6)
        rest1num=-xik*((xk+x0)**2*((-319+600*log2-150*logm)*xk2
     -        +2*(-29+72*log2-18*logm)*xk*x0+(-7+24*log2-6*logm)*x02))
c dropped xik3*mlogm O(h^6logh)

        rest=reglog+reglog1 + rest0num/(2*denom) + m*rest1num/(8*denom)

c multiply by density at xk
        rest=f0k*rest
      endif
#endif
      return
      end


      SUBROUTINE corrqv12(x0,f0,f1,f2,f3,corrs,
     +       largecorr,xk,xik,f0k,trapsing,corr,rest,hk) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*),xk,xik,f0k,trapsing(*)
      logical largecorr
c OUTPUT 
      MODE corr,rest,hk
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),q(0:6),
     +       piece1(0:6),piece2(0:6),piece3(0:6),piece4(0:6),
     +       corrqrho4,corrqrho2,corrqlog,corr1,corr2,corr3
      MODE hk1,hk2,hk3,piece5(0:6)
      MODE denom,Fk,reglog,reglog1,rest0num,rest1num
      
c 1/RHO4 --------------------------------------------------
      lead = 8

c Piece 1: small=delr*xi^2=ab2,   fact=-(x+x0)
      do j=0,4
        fact(j)= -( x1(j) + ones(j)*x0 )
      enddo
      call prodrho4(fact,ir,ab2,piece1)

c Piece 2: small=xi^4=b4,  fact=1
      call prodrho4(ones,ir,b4,piece2)

c Sum Pieces 
      do j=0,6
        numer(j)=piece1(j)+piece2(j)
      enddo
      corr1=corrqrho4(lead,numer,f0,f1,f2,f3,e3,e4,e5,corrs) 
c      hk1=corrqrho4(lead,num0,num1,num2,num3,num4,num5,num6,
c     +                f0,f1,f2,f3,e3,e4,e5,trapsing) 

c 1/RHO2 --------------------------------------------------
      lead = -6

c Piece 1: small=xi^2=b2,  fact=2(...)
      do j=0,3
        fact(j) = 3*x2(j) +4*x1(j)*x0 +ones(j)*x02
      enddo
      call prodrho2(fact,ir3,b2,piece1)

c Piece 2: small=xi^4=b4,  fact=1
      call prodrho2(ones,ir3,b4,piece2)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)+piece2(j)
      enddo

      corr2=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
c      hk2=corrqrho2(lead,num0,num1,num2,num3,num4,
c     +       ir30,ir31,ir32,ir33,f0,f1,f2,f3,e3,e4,trapsing)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = 3.q0/4

c Piece 1: small=xi^2=b2,  fact=
      do j=0,2
        fact(j) = 7*x2(j) +8*x1(j)*x0 +ones(j)*x02
      enddo
      call prodlog(fact,ir5,b2,piece1)

c Piece 2: small=xi^4=b4,  fact=1
      call prodlog(ones,ir5,b4,piece2)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 2 (skip xi^4)
      fact(0)=b2(0)*3*(11*x2(0)+12*x1(0)*x0+x02)
      denom=16*r2**7
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces 
      do j=0,2
        numer(j)=piece1(j)+piece2(j)  +piece5(j)/lead 
      enddo
      corr3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      hk3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,trapsing)
      
      corr=corr1+corr2+corr3

      rest=0
#ifdef ROUNDOFFFIX
c REST (value of Gk-Hk at alf(k0)) ---------------------------------------
c G-H is found in mathematica (H approx F2/rho1^4 + F1/rho1^2 +F0*log(rho1sq))
c it has 2 components
c     logterm =-F0*log(rho2sq)/r2^5  
c     rest = rest0 + m(rest11 + rest12*logm)
c     
c     drop high order xik terms. Careful: xik=O(h), not O(d). Also m=O(h^2)
c     drop m^2logm = O(h^4logh) terms, which are mult by h. Thus O(h^5logh)
c     BELOW: logr2k=log(r2ksq)

      if (largecorr) then
        lead =3
        denom= 4*xk*r2k5

        Fk = xik2*(7*xk2+8*xk*x0+x02+xik2)
        reglog=-lead*Fk/denom*logr2k

        Fk= 3*xik2*(11*xk2+12*xk*x0+x02)
        reglog1=-Fk/(4*denom)*m*logr2k

        rest0num=-3*xik2*(xk+x0)*((-15+56*log2)*xk+(-1+8*log2)*x0)
     +             +xik4*(3-24*log2)
c dropped xik6 term O(h^7)

        rest1num=-xik2*((xk+x0)*((-145+264*log2-66*logm)*xk
     +                             +(-11+24*log2-6*logm)*x0))
c dropped xik4*m*logm term O(h^7*logh)

        rest=reglog+reglog1 + rest0num/(2*denom) + m*rest1num/(8*denom)

c multiply by density at xk
        rest=f0k*rest
      endif
#endif
      return
      end



      SUBROUTINE corrqv22(x0,f0,f1,f2,f3,corrs,
     +       largecorr,xk,xik,f0k,trapsing,corr,rest,hk) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*),xk,xik,f0k,trapsing(*)
      logical largecorr
c OUTPUT 
      MODE corr,rest,hk
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),q(0:6),
     +       piece1(0:6),piece2(0:6),piece3(0:6),piece4(0:6),
     +       corrqrho4,corrqrho2,corrqlog,corr1,corr2,corr3
      MODE hk1,hk2,hk3,piece5(0:6)
      MODE denom,Fk,reglog,reglog1,rest0num,rest1num
      
c 1/RHO4 --------------------------------------------------
      lead = -16

c Piece 1: small=xi^3=b3,   fact=x
      call prodrho4(x1,ir,b3,piece1)

c Sum Pieces 
      do j=0,6
        numer(j)=piece1(j)
      enddo
      corr1=corrqrho4(lead,numer,f0,f1,f2,f3,e3,e4,e5,corrs) 
c      hk1=corrqrho4(lead,num0,num1,num2,num3,num4,num5,num6,
c     +                f0,f1,f2,f3,e3,e4,e5,trapsing) 

c 1/RHO2 --------------------------------------------------
      lead = -12

c Piece 1: small=xi^3=b3,  fact=x
      call prodrho2(x1,ir3,b3,piece1)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)
      enddo

      corr2=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
c      hk2=corrqrho2(lead,num0,num1,num2,num3,num4,
c     +       ir30,ir31,ir32,ir33,f0,f1,f2,f3,e3,e4,trapsing)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = 9.q0/2

c Piece 1: small=xi^3=b3,  fact=x
      call prodlog(x1,ir5,b3,piece1)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 3 (skip xi^2)
      fact(0)=15*x1(0)*b3(0)
      denom=8*r2**7
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces 
      do j=0,2
        numer(j)=piece1(j)  +piece5(j)/lead 
      enddo
      corr3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      hk3=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,trapsing)
      
      corr=corr1+corr2+corr3

      rest=0
#ifdef ROUNDOFFFIX
c REST (value of Gk-Hk at alf(k0)) ---------------------------------------
c G-H is found in mathematica (H approx F2/rho1^4 + F1/rho1^2 +F0*log(rho1sq))
c it has 2 components
c     logterm =-F0*log(rho2sq)/r2^5  
c     rest = rest0 + m(rest11 + rest12*logm)
c     
c     drop high order xik terms. Careful: xik=O(h), not O(d). Also m=O(h^2)
c     drop m^2logm = O(h^4logh) terms, which are mult by h. Thus O(h^5logh)
c     BELOW: logr2k=log(r2ksq)

      if (largecorr) then
        lead =9
        denom= 2*r2k5

        Fk = xk*xik3
        reglog=-lead*Fk/denom*logr2k

        Fk=15*xk*xik3
        reglog1=-Fk/(4*denom)*m*logr2k

        rest0num= -3*(-7+24*log2)*xk*xik3
        rest1num= xik3*(67-120*log2+30*logm)*xk
        rest=reglog+reglog1 + rest0num/(2*denom) + m*rest1num/(8*denom)

c multiply by density at xk
        rest=f0k*rest
      endif
#endif
      return
      end
      

      MODE FUNCTION corrqrho4(lead,num,f0,f1,f2,f3,e3,e4,e5,corr)
      implicit none
c INPUT      
      MODE lead,num(0:*),f0,f1,f2,f3,e3,e4,e5,corr(*)
c LOCAL      
      MODE g0,g1,g2,g3,g4,g5,g6,
     +   c20,c21,c22,c23,c24,c25,c26,c33,c34,c35,c36,c37,c38,c46,c47,
     +   c48,c49,c410,c59,c510,c511,c512

c Use f to second order and multiply mnum*f=t to O(d^6+...+alf^6)
c (only need f to third order since mnum is O(d^3) ) 
      g0= num(0)*f0
      g1= num(1)*f0 + num(0)*f1
      g2= num(2)*f0 + num(1)*f1 + num(0)*f2
      g3= num(3)*f0 + num(2)*f1 + num(1)*f2 + num(0)*f3
      g4= num(4)*f0 + num(3)*f1 + num(2)*f2 + num(1)*f3
      g5= num(5)*f0 + num(4)*f1 + num(3)*f2 + num(2)*f3
      g6= num(6)*f0 + num(5)*f1 + num(4)*f2 + num(3)*f3

c      print*,lead,' lead'
c      print*,num4,' num4'
c      print*,lead*g4,' lead*num/r2*uy 4'
c      print*,'lead*um x/r2*vx'
c      print*,lead*g0
c      print*,lead*g1
c      print*,lead*g2
c      print*,lead*g3
c      print*,lead*g4
c      print*,lead*g5

c Compute coefficients of g/(r1sq+e3*alf^3+e4*alf^4...) in expansion
c Here: cjk = coefficient of Hjk = alf^k/r1sq**j
c
c Remember: g0 = O(d^2), g1=O(d), I am not including terms 
c                                 smaller than needed !!!

c Coefficients of rhosq**2 = 1/r4 terms, need O(alf^6 ... d^6) on top
      c20 = g0    !coeff of     1/r4
      c21 = g1    !coeff of alf^1/r4
      c22 = g2    !coeff of alf^2/r4
      c23 = g3    !coeff of alf^3/r4
      c24 = g4    !coeff of alf^4/r4
      c25 = g5    !coeff of alf^5/r4
      c26 = g6    !coeff of alf^6/r4

c Coefficients of rhosq**3 = 1/r6 terms, need O(alf^8... d^8) on top
      c33 = -2*                 e3*g0  !coeff of alf^3/r6
      c34 = -2*(        e4*g0 + e3*g1) !coeff of alf^4/r6
      c35 = -2*(e5*g0 + e4*g1 + e3*g2) !coeff of alf^5/r6
      c36 = -2*(e5*g1 + e4*g2 + e3*g3) !coeff of alf^6/r6
      c37 = -2*(e5*g2 + e4*g3 + e3*g4) !coeff of alf^7/r6 (skip g0 term)
      c38 = -2*(e5*g3 + e4*g4 + e3*g5) !coeff of alf^8/r6 (skip g0, g1 terms)

c Coefficients of rhosq**4 = 1/r8 terms, need O(alf^10... d^10) on top
      c46 = 3*e3**2*g0              !coeff of alf^6/r8
      c47 = 3*e3*(2*e4*g0 + e3*g1)  !coeff of alf^7/r8
      c48 = 3*e3*(2*e4*g1 + e3*g2)  !coeff of alf^8/r8
      c49 = 3*e3*(2*e4*g2 + e3*g3)  !coeff of alf^9/r8  
      c410 =3*e3*(2*e4*g3 + e3*g4)  !coeff of alf^10/r8  

c Coefficients of rhosq**5 = 1/r10 terms, need O(alf^12... d^12) on top
      c59  = -4*e3**3*g0    !coeff of alf^9/r10
      c510 = -4*e3**3*g1    !coeff of alf^10/r10
      c511 = -4*e3**3*g2    !coeff of alf^11/r10
      c512 = -4*e3**3*g3    !coeff of alf^12/r10

c Compute correction for M11 m11num*f/r2 /r1sq
      corrqrho4 = lead*(c20*corr(9) ! 1/arg2  d3 h/d4 =h/d
     +                +c21*corr(10)   ! u/arg2 d2 h/d3=h/d
     +                +c22*corr(11)   ! u2/arg2 d h/d2=h/d
     +                +c23*corr(12)  ! u3/arg2  h/d
     +                +c24*corr(13)  ! u4/arg2  h
     +                +c25*corr(14)  ! u5/arg2  hd
     +                +c33*corr(16)  ! u3/arg3  d3h/d3 =h
     +                +c34*corr(17)  ! u4/arg3  d2h/d2 =h
     +                +c35*corr(18)  ! u5/arg3  dh/d  =h
     +                +c36*corr(19)  ! u6/arg3  h
     +                +c37*corr(20)  ! u7/arg3  hd
     +                +c46*corr(22)  ! u6/arg4  d3h/d2 =hd
     +                +c47*corr(23)  ! u7/arg4  d2h/d =hd
     +                +c48*corr(24)  ! u8/arg4  d h =hd
     +                +c49*corr(25)  ! u9/arg4  hd
#ifdef FOURTH
     +                +c26*corr(15)  ! u6/arg2  hd2 SKIP FOR O(h3)
     +                +c38*corr(21)  ! u8/arg3  hd2  SKIP FOR O(h3)
     +                +c410*corr(26) ! u10/arg4 hd2  SKIP FOR O(h3)
     +                +c59*corr(27)  ! u9/arg5 d3*h/d=hd2 SKIPfrO(h3)
     +                +c510*corr(28) ! u10/arg5 d2*h=hd2 SKIPforO(h3)
     +                +c511*corr(29) ! u11/arg5 d*hd=hd2 SKIPforO(h3)
     +                +c512*corr(30) ! u12/arg5 hd2  SKIP FOR O(h3)
#endif
     +  )

c     TEST : 
c     for qu22, corrqrho4 agrees to 1.e-17 with original code
c      print*,'first 5'
c      print*,'all rho4'
c      print*,lead*c20*corr(9)
c      print*,lead*c21*corr(10)   
c      print*,lead*c22*corr(11)
c      print*,lead*c23*corr(12)  
c      print*,lead*c24*corr(13)  
c       print*,lead*c25*corr(14)  ! u5/arg2  hd
c       print*,lead*c33*corr(16)  ! u3/arg3  d3h/d3 =h
c       print*,lead*c34*corr(17)  ! u4/arg3  d2h/d2 =h
c       print*,lead*c35*corr(18)  ! u5/arg3  dh/d  =h
c       print*,lead*c36*corr(19)  ! u6/arg3  h
c       print*,lead*c37*corr(20)  ! u7/arg3  hd
c       print*,lead*c46*corr(22)  ! u6/arg4  d3h/d2 =hd
c       print*,lead*c47*corr(23)  ! u7/arg4  d2h/d =hd
c       print*,lead*c48*corr(24)  ! u8/arg4  d h =hd
c       print*,lead*c49*corr(25)  ! u9/arg4  hd

      return
      end


      MODE FUNCTION corrqrho2(lead,num,f0,f1,f2,f3,e3,e4,corr)
c for all rho2 terms, mnum is O(d^2) : integrable, not PV!
c Note: For Qs 1/rh0**2 terms are multiplied by 1/r2^3
      implicit none
c INPUT
      MODE lead,num(0:*),f0,f1,f2,f3,e3,e4,corr(*)
c LOCAL
      MODE g0,g1,g2,g3,g4,
     +  c10,c11,c12,c13,c14,c23,c24,c25,c26,c36,c37,c38,pi8

c Use f to second order and multiply mnum*f=t to O(d^4+...+alf^4)
c (only need f to second order since mnum is O(d^2) ) Actually: 
c  only need f to first order since mnum is O(d^3)! but will proceed as before )
      g0= num(0)*f0
      g1= num(1)*f0 + num(0)*f1
      g2= num(2)*f0 + num(1)*f1 + num(0)*f2 
      g3= num(3)*f0 + num(2)*f1 + num(1)*f2 + num(0)*f3
      g4= num(4)*f0 + num(3)*f1 + num(2)*f2 + num(1)*f3

c      print*,lead,' lead'
c      print*,'lead*num x/r2*vx'
c      print*,lead*g0
c      print*,lead*g1
c      print*,lead*g2
c      print*,lead*g3
c      print*,lead*g4
c Compute coefficients of g/(r1sq+e3*alf^3+e4*alf^4...) in expansion
c Here: cjk = coefficient of Hjk = alf^k/r1sq**j
c
c Remember: g0 = O(d^2), g1=O(d), I am not including terms 
c                                 smaller than needed !!!

c Coefficients of 1/rhosq = 1/r2 terms, need O(alf^4 ... d^4) on top
      c10 = g0    !coeff of     1/r2
      c11 = g1    !coeff of alf^1/r2
      c12 = g2    !coeff of alf^2/r2
      c13 = g3    !coeff of alf^3/r2
      c14 = g4    !coeff of alf^4/r2

c Coefficients of 1/rhosq**2 = 1/r4 terms, need O(alf^6... d^6) on top
      c23 =       -e3*g0    !coeff of alf^3/r4
      c24 = -e4*g0-e3*g1    !coeff of alf^4/r4
      c25 = -e4*g1-e3*g2    !coeff of alf^5/r4 (skip g0 term)
      c26 = -e4*g2-e3*g3    !coeff of alf^6/r4 (skip g0, g1 terms)

c Coefficients of 1/rhosq**3 = 1/r6 terms, need O(alf^8... d^8) on top
      c36 = e3**2*g0  !coeff of alf^6/r6
      c37 = e3**2*g1  !coeff of alf^7/r6  (skip g0 term)
      c38 = e3**2*g2  !coeff of alf^8/r6  (skip g0,g1 term)

c Compute correction for M11 m11num*f/r2 /r1sq
      corrqrho2 = lead*(c10*corr(4)   ! 1/arg d2 h/d2 =h
     +                +c11*corr(5)   ! u/arg d h/d =h
     +                +c12*corr(6)   ! u2/arg h
     +                +c13*corr(7)  ! u3/arg hd
     +                +c23*corr(12)  ! u3/arg2 d2 h/d =hd
     +                +c24*corr(13)  ! u4/arg2 d h =hd
     +                +c25*corr(14)  ! u5/arg2 hd
#ifdef FOURTH
     +                +c14*corr(8)  ! u4/arg  hd2 SKIP FOR Oh3
     +                +c26*corr(15) ! u6/arg2 hd2 SKIP FOR Oh3
     +                +c36*corr(19) ! u6/arg3 d2*h=hd2 SKIPfrOh3
     +                +c37*corr(20) ! u7/arg3 d*hd=hd2 SKIPfrOh3
     +                +c38*corr(21) ! u8/arg3 hd2 SKIP FOR Oh3
#endif
     +  )

c      pi8=32*atan(1.q0)
c      print*,'all'
c      print*,lead*c10*corr(4)/pi8 !' %1/arg'
c      print*,lead*c11*corr(5)/pi8 !' %u/arg'
c      print*,lead*c12*corr(6)/pi8 !' %u2/arg'
c      print*,lead*c13*corr(7)/pi8 !'%u3/arg'
c      print*,lead*c23*corr(12)/pi8 !'%u3/arg2'
c      print*,lead*c24*corr(13)/pi8 !'%u4/arg2'
c      print*,lead*c25*corr(14)/pi8 !'%u5/arg2'
c      print*,lead*c14*corr(8)/pi8  ! u4/arg  hd2 SKIP FOR Oh3
c      print*,lead*c26*corr(15)/pi8 ! u6/arg2 hd2 SKIP FOR Oh3
c      print*,lead*c36*corr(19)/pi8 ! u6/arg3 d2*h=hd2 SKIPfrOh3
c      print*,lead*c37*corr(20)/pi8 ! u7/arg3 d*hd=hd2 SKIPfrOh3
c      print*,lead*c38*corr(21)/pi8 ! u8/arg3 hd2 SKIP FOR Oh3

      return
      end


      MODE FUNCTION corrqlog(lead,num,f0,f1,f2,f3,e3,e4,corr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Note: For Qs, log(rho**2) is multiplied by 1/r2^5
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c INPUT
      MODE lead,num(0:*),f0,f1,f2,f3,e3,e4,corr(*)
c LOCAL
      MODE g0,g1,g2,g3,c00,c01,c02,c03,c13,c14,c26,pi8

c Use f to second order and multiply lognum*f=t to O(d^2+...+alf^2)
      g0= num(0)*f0
      g1= num(0)*f1 + num(1)*f0
      g2= num(0)*f2 + num(1)*f1 + num(2)*f0
c      g3= num(0)*f3 + num(1)*f2 + num(2)*f1 + num(3)*f0
      g3=0

c coefficients of log(rhosq)  (rhosq=rho2)
      c00 = g0    !coeff of       log(rhosq)
      c01 = g1    !coeff of alf^1*log(rhosq)
      c02 = g2    !coeff of alf^2*log(rhosq)
      c03 = g3    !coeff of alf^2*log(rhosq)
c coefficients of 1/rhosq = 1/rhosq terms, need O(alf^4... d^4) on top
      c13 =        e3*g0    !coeff of alf^3/r2
      c14 =  e4*g0+e3*g1    !coeff of alf^4/r2
c coefficients of 1/rhosq**2 = 1/rhosq^2 terms, need O(alf^6... d^6) on top
      c26 = -e3**2*g0/2  !coeff of alf^6/r4

      corrqlog = lead*( c00*corr(1)  ! logarg h
     +                +c01*corr(2)   ! u*logarg h2
     +                +c13*corr(7)   ! u3/arg hd 
#ifdef FOURTH
     +                +c02*corr(3)   ! u2*logarg h3 SKIP FOR Oh3
c     +                +c03*corr(31)  ! u3*logarg h4 ADD for Oh4
     +                +c14*corr(8)   ! u4/arg hd2 SKIP FOR Oh3
     +                +c26*corr(15)  ! u6/arg2 hd2 SKIP FOR Oh3
#endif
     +   )

c      pi8=32*atan(1.q0)
c      print*,'all'
c      print*,lead*c00*corr(1)/pi8,'% logarg'
c      print*,lead*c01*corr(2)/pi8,'% ulogarg'
c      print*,lead*c13*corr(7)/pi8,'% u3/arg'
c      print*,lead*c02*corr(3)/pi8,'% u2logarg'
c      print*,lead*c14*corr(8)/pi8,'% u4/arg'
c      print*,corrqlog,c00,corr(1),lead

      return
      end


      SUBROUTINE compcoeffnew(xb,x0,d,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +       xk,xik,cm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Revised April 23, 2022
c This routine computes all coefficients that are common to many
c subroutines, including the coefficients
c
c      e3,e4,e5,cm
c the coefficients of the factors in the denominator (r2=rho2)
c      1/r2   = ir     1/(x*r2)   =irx
c      1/r2^3 = ir3    1/(x*r2^3) =ir3x
c      1/r2^5 = ir5    1/(x*r2^5) =ir5x
c and the coefficients of the "small" terms
c
c      a1=r-r0           ab=(r-r0)*xi
c      a2=(r-r0)^2       ab2=(r-r0)*xi^2
c      a3=(r-r0)^3       ab3=(r-r0)*xi^3
c      b1=xi             ab4=(r-r0)*xi^4
c      b2=xi^2           a2b=(r-r0)^2*xi
c      b3=xi^3           a2b2=(r-r0)^2*xi^2
c      b4=xi^4
c      b5=xi^5
c      b6=xi^6
c
c Each one of the small terms gets 7, 5 or 3 coefficients 
c (to order 6,4 or 2, see below)

c Each of the other factors gets 4 coefficients (to order 3)

c These are all the coefficients needed in the 4th order method. 
c Note that for O(h^4)
c
c    Q2=F2/rho^4*sigma    requires triple roots in F2 
c                         and those to c6*alf^6, sigma to c3*alf^3
c    Q1=F1/rho^2*sigma    has only double roots in F1 (in this axi
c                         Stokes case there happen to be no simple roots) 
c                         so we require all the double roots to
c                         c4*alf^4, sigma to c2*alf^2
c    Q0=F0*log(rho)*sigma requires F0,sigma to c3*alf^2 (slightly more
c                         than necessary
c 
c For convenience, there are all put in COMMON block
c
c All the coefficients are either computed in Mathematica, see expmiscredo.pdf
c or obtained by multiplying out polynomials with known coefficients.This
c latter approach is equivalent to the coefficients obtained w Mathematica
c by expanding for example (r-r0)^3 with known expansions for r, but
c more transparent, and requires less operations.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT       
      MODE xb,x0,d,c1,c2,c3,c4,dy,d1,d2,d3,d4,xk,xik
c OUTPUT : all variables in coeff.h and      
      MODE cm
c LOCAL
      MODE dx,xd,xd2,xd3,xd4,yd,yd2,yd3,yd4,
     +  rs0,rs1,rs2,rs3,rs4,ix(0:6),r1sqk,r2sqk

      r1sqk=(xk-x0)**2+xik**2
      r2sqk=(xk+x0)**2+xik**2

      r2k5=sqrt(r2sqk)**5
      logr2k=log(r2sqk)

      m=r1sqk/r2sqk
      logm=log(m)

      log2=log(2.q0)
      x02  = x0**2
      x03  = x0**3
      x04  = x0**4
      x05  = x0**5
      xk2  = xk**2
      xk3  = xk**3
      xk4  = xk**4
      xk5  = xk**5
      xk6  = xk**6
      xik2 = xik**2
      xik3 = xik**3
      xik4 = xik**4

      dx=xb-x0
      xd=c1
      xd2=2*c2
      xd3=6*c3
      xd4=24*c4
      yd=d1
      yd2=2*d2
      yd3=6*d3
      yd4=24*d4

      dsq=d**2
      cmsq=(xd**2+yd**2)+dx*xd2+dy*yd2
      cm=sqrt(cmsq)
      e3=xd*xd2+yd*yd2 + (dx*xd3+dy*yd3)/3
      e4=(xd2**2+yd2**2)/4+(dx*xd4+dy*yd4)/12+(xd3*xd+yd3*yd)/3
      e5=(2*(xd2*xd3+yd2*yd3) + (xd4*xd+yd4*yd))/12

c coefficients of 1
      ones(0)=1
      ones(1)=0
      ones(2)=0
      ones(3)=0
      ones(4)=0
      ones(5)=0
      ones(6)=0

c coefficients of rho2**2 (see expmiscredo.pdf
c      rs0= dx**2 + dy**2 + 4*dx*x0 + 4*x0**2
c      rs1= 2*(d1*dy+c1*(dx+2*x0))
c      rs2= c1**2+d1**2+2*c2*dx+2*d2*dy+4*c2*x0
c      rs3= 2*(c1*c2+d1*d2+c3*dx+d3*dy+2*c3*x0)
      rs0= dy**2 + (x0+xb)**2
      rs1= 2*(d1*dy+c1*(x0+xb))
      rs2= c1**2+d1**2+2*(d2*dy+c2*(x0+xb))
      rs3= 2*(c1*c2+d1*d2+d3*dy+c3*(x0+xb))
      rs4= c2**2+2*c1*c3+d2**2+2*d1*d3+2*d4*dy+2*c4*x0+2*c4*xb

      r2=sqrt(rs0)

c coefficients of 1/r2
      ir(0)=1/r2
      ir(1)=-rs1/(2*r2**3)
      ir(2)= (3*rs1**2 - 4*rs0*rs2)/(8*r2**5)
      ir(3)=(-5*rs1**3 + 12*rs0*rs1*rs2 - 8*rs0**2*rs3)/(16*r2**7)
      ir(4)=(35*rs1**4-120*rs0*rs1**2*rs2+96*rs0**2*rs1*rs3+
     -       16*rs0**2*(3*rs2**2-4*rs0*rs4))/(128*r2**9)

c coefficients of 1/r2^3
      ir3(0)= 1/r2**3
      ir3(1)= (-3*rs1)/(2*r2**5)
      ir3(2)= 3*(5*rs1**2 - 4*rs0*rs2)/(8*r2**7)
      ir3(3)= (-35*rs1**3 + 60*rs0*rs1*rs2 - 24*rs0**2*rs3)/(16*r2**9)
      ir3(4)=(315*rs1**4-840*rs0*rs1**2*rs2+480*rs0**2*rs1*rs3+
     -        48*rs0**2*(5*rs2**2-4*rs0*rs4))/(128*r2**11)

c coefficients of 1/r2^5
      ir5(0)=  1/r2**5
      ir5(1)= -5*rs1/(2*r2**7)
      ir5(2)=  5*(7*rs1**2 - 4*rs0*rs2)/(8*r2**9)
      ir5(3)= -5*(21*rs1**3 - 28*rs0*rs1*rs2 + 8*rs0**2*rs3)/(16*r2**11)
      ir5(4)= (5*(231*rs1**4-504*rs0*rs1**2*rs2+224*rs0**2*rs1*rs3+
     -         16*rs0**2*(7*rs2**2-4*rs0*rs4)))/(128*r2**13)

c coefficients of 1/x
      ix(0) =  1/xb
      ix(1) = -c1/xb**2
      ix(2) =  (c1**2 - c2*xb)/xb**3
      ix(3) = -(c1**3 - 2*c1*c2*xb + c3*xb**2)/xb**4
      ix(4) = (c1**4-3*c1**2*c2*xb+2*c1*c3*xb**2+xb**2*(c2**2-c4*xb))
     +         /xb**5
c      call prod3(ir,ix,irx)   ! 1/(x*r2)
c      call prod2(ir3,ix,ir3x) ! 1/(x*r2^3)
c      call prod2(ir5,ix,ir5x) ! 1/(x*r2^5)
      call prod4(ir,ix,irx)   ! 1/(x*r2)
      call prod4(ir3,ix,ir3x) ! 1/(x*r2^3)
      call prod4(ir5,ix,ir5x) ! 1/(x*r2^5)

c coefficients of x
      x1(0)=xb
      x1(1)=c1
      x1(2)=c2
      x1(3)=c3
      x1(4)=c4
c coefficients of x2=x^2
c      call prod3(x1,x1,x2)  !x2=x^2
      call prod446(x1,x1,x2)  !x2=x^2
c      call prod3(x1,x2,x3)  !x3=x^2
c      call prod3(x2,x2,x4)  !x4=x^2
c      call prod3(x2,x3,x5)  !x5=x^2
      call prod4(x1,x2,x3)  !x3=x^2
      call prod4(x2,x2,x4)  !x4=x^2
      call prod4(x2,x3,x5)  !x5=x^2
c coefficients of a1=r-r0=x-x0
      a1(0)=dx
      a1(1)=c1
      a1(2)=c2
      a1(3)=c3
      a1(4)=c4
      call prod446(a1,a1,a2)    ! a2=delr^2
      call prod466(a1,a2,a3)  ! a3=delr^3
c coefficients of b1=xi
      b1(0)=dy
      b1(1)=d1
      b1(2)=d2
      b1(3)=d3
      b1(4)=d4
c      call prod4(b1,b1,b2)    ! b2=xi^2
      call prod446(b1,b1,b2)    ! b2=xi^2
      call prod466(b1,b2,b3)  ! b3=xi^3
      call prod6(b2,b2,b4)  ! b4=xi^4
      call prod6(b3,b2,b5)  ! b5=xi^5
      call prod6(b3,b3,b6)  ! b6=xi^6
      call prod446(a1,b1,ab1) ! ab1=delr*xi
      call prod466(a1,b2,ab2) ! ab2=delr*xi^2
      call prod466(a1,b3,ab3) ! ab3=delr*xi^3
      call prod466(a1,b4,ab4) ! ab4=delr*xi^4
      call prod466(b1,a2,a2b) ! a2b=delr^2*xi
      call prod6(a2,b2,a2b2)! a2b2=delr^2*xi^2

c      a2b20=a20*b20
c      a2b21=a20*b21 + a21*b20
c      a2b22=a20*b22 + a21*b21 + a22*b20
c      a2b23=a20*b23 + a21*b22 + a22*b21 + a23*b20
c      a2b24=a20*b24 + a21*b23 + a22*b22 + a23*b21 + a24*b20
c      a2b25=          a21*b24 + a22*b23 + a23*b22 + a24*b21
c      a2b26=                    a22*b24 + a23*b23 + a24*b22

      return
      end


      subroutine prod1(x,y,z)
c     computes z=x*y where x,y,z are poly or degree 2
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2)= 0
      z(3)=0
      z(4)=0
      z(5)=0
      z(6)=0
      return 
      end


      subroutine prod2(x,y,z)
c     computes z=x*y where x,y,z are poly or degree 2
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3)=0
      z(4)=0
      z(5)=0
      z(6)=0
      return 
      end

      subroutine prod3(x,y,z)
c     computes z=x*y where x,y,z are poly or degree 3
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3) = x(0)*y(3) + x(1)*y(2) + x(2)*y(1) + x(3)*y(0)
      z(4)=0
      z(5)=0
      z(6)=0
      return 
      end

      subroutine prod4(x,y,z)
c     computes z=x*y where x,y,z are degree 4
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3) = x(0)*y(3) + x(1)*y(2) + x(2)*y(1) + x(3)*y(0)
      z(4) = x(0)*y(4) + x(1)*y(3) + x(2)*y(2) + x(3)*y(1) + x(4)*y(0)
      z(5)=0
      z(6)=0
      return 
      end

      subroutine prod5(x,y,z)
c     computes z=x*y where x,y,z are degree 4
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0)= x(0)*y(0)
      z(1)= x(0)*y(1)+x(1)*y(0)
      z(2)= x(0)*y(2)+x(1)*y(1)+x(2)*y(0)
      z(3)= x(0)*y(3)+x(1)*y(2)+x(2)*y(1)+x(3)*y(0)
      z(4)= x(0)*y(4)+x(1)*y(3)+x(2)*y(2)+x(3)*y(1)+x(4)*y(0)
      z(5)= x(0)*y(5)+x(1)*y(4)+x(2)*y(3)+x(3)*y(2)+x(4)*y(1)+x(5)*y(0)
      z(6)=0
      return 
      end

      subroutine prod466(x,y,z)
c     computes z=x*y where x,y,z are degree 4
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0)= x(0)*y(0)
      z(1)= x(0)*y(1)+x(1)*y(0)
      z(2)= x(0)*y(2)+x(1)*y(1)+x(2)*y(0)
      z(3)= x(0)*y(3)+x(1)*y(2)+x(2)*y(1)+x(3)*y(0)
      z(4)= x(0)*y(4)+x(1)*y(3)+x(2)*y(2)+x(3)*y(1)+x(4)*y(0)
      z(5)= x(0)*y(5)+x(1)*y(4)+x(2)*y(3)+x(3)*y(2)+x(4)*y(1)
      z(6)= x(0)*y(6)+x(1)*y(5)+x(2)*y(4)+x(3)*y(3)+x(4)*y(2)
      return 
      end

      subroutine prod6(x,y,z)
c     computes z=x*y where x,y,z are degree 4
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0)= x(0)*y(0)
      z(1)= x(0)*y(1)+x(1)*y(0)
      z(2)= x(0)*y(2)+x(1)*y(1)+x(2)*y(0)
      z(3)= x(0)*y(3)+x(1)*y(2)+x(2)*y(1)+x(3)*y(0)
      z(4)= x(0)*y(4)+x(1)*y(3)+x(2)*y(2)+x(3)*y(1)+x(4)*y(0)
      z(5)= x(0)*y(5)+x(1)*y(4)+x(2)*y(3)+x(3)*y(2)+x(4)*y(1)+x(5)*y(0)
      z(6)= x(0)*y(6)+x(1)*y(5)+x(2)*y(4)+x(3)*y(3)+x(4)*y(2)+x(5)*y(1)
     +     +x(6)*y(0)
      return 
      end

      subroutine prod122(x,y,z)
c     computes z=x*y where x is degree 3 and y,z are degree 4
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) 
      z(3)= 0
      z(4)=0
      z(5)=0
      z(6)=0
      return 
      end

      subroutine prod133(x,y,z)
c     computes z=x*y where x is degree 3 and y,z are degree 4
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) 
      z(3) = x(0)*y(3) + x(1)*y(2) 
      z(4)=0
      z(5)=0
      z(6)=0
      return 
      end

      subroutine prod233(x,y,z)
c     computes z=x*y where x is degree 3 and y,z are degree 4
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3) = x(0)*y(3) + x(1)*y(2) + x(2)*y(1) 
      z(4)=0
      z(5)=0
      z(6)=0
      return 
      end

      subroutine prod344(x,y,z)
c     computes z=x*y where x is degree 3 and y,z are degree 4
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3) = x(0)*y(3) + x(1)*y(2) + x(2)*y(1) + x(3)*y(0)
      z(4) = x(0)*y(4) + x(1)*y(3) + x(2)*y(2) + x(3)*y(1)
      z(5)=0
      z(6)=0
      return 
      end

      subroutine prod446(x,y,z)
c     computes z=x*y where x,y are degree 4, z is degree 6
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3) = x(0)*y(3) + x(1)*y(2) + x(2)*y(1) + x(3)*y(0)
      z(4) = x(0)*y(4) + x(1)*y(3) + x(2)*y(2) + x(3)*y(1) + x(4)*y(0)
      z(5) =             x(1)*y(4) + x(2)*y(3) + x(3)*y(2) + x(4)*y(1)
      z(6) =                         x(2)*y(4) + x(3)*y(3) + x(4)*y(2)
      return 
      end

      subroutine prod366(x,y,z)
c     computes z=x*y where x,y are degree 4, z is degree 6
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3) = x(0)*y(3) + x(1)*y(2) + x(2)*y(1) + x(3)*y(0)
      z(4) = x(0)*y(4) + x(1)*y(3) + x(2)*y(2) + x(3)*y(1) 
      z(5) = x(0)*y(5) + x(1)*y(4) + x(2)*y(3) + x(3)*y(2)
      z(6) = x(0)*y(6) + x(1)*y(5) + x(2)*y(4) + x(3)*y(3)
      return 
      end

      subroutine prod255(x,y,z)
c     computes z=x*y where x,y are degree 4, z is degree 6
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3) = x(0)*y(3) + x(1)*y(2) + x(2)*y(1)
      z(4) = x(0)*y(4) + x(1)*y(3) + x(2)*y(2)
      z(5) = x(0)*y(5) + x(1)*y(4) + x(2)*y(3) 
      z(6)=0
      return 
      end

      subroutine prod244(x,y,z)
c     computes z=x*y where x,y are degree 4, z is degree 6
      implicit none
      MODE x(0:*),y(0:*),z(0:*)
      z(0) = x(0)*y(0)
      z(1) = x(0)*y(1) + x(1)*y(0)
      z(2) = x(0)*y(2) + x(1)*y(1) + x(2)*y(0)
      z(3) = x(0)*y(3) + x(1)*y(2) + x(2)*y(1)
      z(4) = x(0)*y(4) + x(1)*y(3) + x(2)*y(2)
      z(5)= 0
      z(6)=0
      return 
      end


      SUBROUTINE prodrho4(fact,denom,small,piece)
      implicit none
c INPUT/OUTPUT
      MODE fact(0:*),denom(0:*),small(0:*),piece(0:*)
c LOCAL
      MODE q(0:6)

#ifdef FOURTH
      call prod3(fact,denom,q)              !q=fact*ir
      call prod366(q,small,piece)         !piece1=fact*ir*a2b
c      call prod4(fact,irx,q)              !q=fact*ir
c      call prod466(q,a3,piece1)         !piece1=fact*ir*a2b
#endif
#ifdef THIRD
      call prod2(fact,denom,q)              !q=fact*ir
      call prod255(q,small,piece)         !piece1=fact*ir*a2b
#endif
      return
      end


      SUBROUTINE prodrho2(fact,denom,small,piece)
      implicit none
c INPUT/OUTPUT
      MODE fact(0:*),denom(0:*),small(0:*),piece(0:*)
c LOCAL
      MODE q(0:6)

#ifdef FOURTH
c      call prod2(fact,denom,q)              !q=fact*ir
c      call prod244(q,small,piece)         !piece1=fact*ir*a2b
      call prod3(fact,denom,q)        !q=fact*ir3x
      call prod344(q,small,piece)      !piece2=fact*ir3x*small
#endif
#ifdef THIRD
      call prod2(fact,denom,q)              !q=fact*ir
      call prod233(q,small,piece)         !piece1=fact*ir*a2b
#endif
      return
      end


      SUBROUTINE prodlog(fact,denom,small,piece)
      implicit none
c INPUT/OUTPUT
      MODE fact(0:*),denom(0:*),small(0:*),piece(0:*)
c LOCAL
      MODE q(0:6)

#ifdef FOURTH
      call prod2(fact,denom,q)              !q=fact*ir
      call prod2(q,small,piece)         !piece1=fact*ir*a2b
#endif
#ifdef THIRD
      call prod1(fact,denom,q)              !q=fact*ir
      call prod1(q,small,piece)         !piece1=fact*ir*a2b
#endif
      return
      end
