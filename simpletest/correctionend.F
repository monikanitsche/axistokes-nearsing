c 06/11/20: change to GLOBAL CORRECTIONS, skip middle term
c 02/17/22: Switch x->r, y->z,        DONE
c           Note f1 = f2inpaper, f2=f1inpaper

#include "flags.h"
      SUBROUTINE correctionmend(alf,ns,z0,zb,d,alfp,corrv,corrud,
     +    corrvdd,rdb,rd3b,zd2b,zd4b,f1db,f2b,f2d2b)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Input: alfp and all derivatives at alfp 
c            the point (r0=0,z0) off sheet on axis at which velo is to be evaluated
c            the normal projection (0,zproj) of that point onto sheet 
c            the distance d from (0,z0) to (0,zb)
c
c     Output: corrv,corrud,corrvdd corrections for v,ud,vdd at (0,z0)
c             coefficients given in GsderAxisExp.nb/pdf
c
c Note: all integrands for r0=0 are odd, so only c1,c3,c5,c7,etc
c       coefficients, that is odd powers of u^k in numerator
c Note
c        corr(1)   u/rt3      ! u/(1+u2)^3/2
c        corr(2)   u3/rt3      !u3/(1+u2)^3/2
c        corr(3)   u/rt5      ! u/(1+u2)^5/2
c        corr(4)   u3/rt5      !u3/(1+u2)^5/2
c        corr(5)   u5/rt5      !u5/(1+u2)^5/2
c        corr(6)   u/rt7       !
c        corr(7)   u3/rt7      !
c        corr(8)   u5/rt7      !u5/(1+u2)^7/2
c        corr(9)   u7/rt7      !u7/(1+u2)^7/2
c        corr(10)  u/rt9      !
c        corr(11)  u3/rt9     !
c        corr(12)  u5/rt9     !
c        corr(13)  u7/rt9     !
c        corr(14)  u5/rt11    !
c        corr(15)  u7/rt11    !
c
c  Note: call zb=zbase, rdb=rdbase, xib=xibase etc
c  SUBSCRIPT B REFERS TO BASE=PROJECTION
c
c  use O(h^2) approximations of all derivatives
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,ns
C     INPUT/OUTPUT
      MODE z0,zb,d,alfp,corrv,corrud,corrvdd,alf(0:*)
      MODE rdb,rd3b,zd2b,zd4b,f1db,f2b,f2d2b
C     LOCAL
      MODE cmsq,cm,corr(20),corrl(20),xib,c1,c3,c5,c7,e4

      xib=zb-z0
      cmsq=rdb**2+xib*zd2b
      cm=sqrt(cmsq)
      e4 = rdb*rd3b/3 + zd2b**2/4  +xib*zd4b/12

      call compallcorrend(alf,alfp,d,cm,ns,corr,corrl)

c v at r=0
      c1= 2*xib**2*f2b*rdb
      c3= f2b*rdb**3 + rdb*(f1db*rdb+2*f2b*zd2b)*xib
     +         + (3*f2d2b*rdb+f2b*rd3b)*xib**2/3
      c5=(rdb*(f2d2b*rdb**2+f1db*rdb*zd2b+f2b*(rdb*rd3b+zd2b**2)))/2
      corrv= c1*corrl(1)+c3*corrl(2) +c5*corrl(19) !3/2 terms
     +      -3*e4/2*(c1*corr(5) +c3*corr(16))   !5/2 terms
c      corrv=2*pi*corrv; corrv=-corrv/pi8
      corrv=-corrv/4

c ud at r=0
      c1= -2*f2b*rdb*xib**3
      c3=  f2b*rdb**3*xib - rdb*(2*f1db*rdb+3*f2b*zd2b)*xib**2
     +     -  (3*f2d2b*rdb + f2b*rd3b)*xib**3/3
      c5=  rdb**3*(2*f1db*rdb+f2b*zd2b)/2
     +    +rdb*(f2d2b*rdb**2 + f2b*rdb*rd3b
     +         -4*f1db*rdb*zd2b - 3*f2b*zd2b**2)*xib/2
      corrud= c1*corr(3)+c3*corr(4) +c5*corr(5) !5/2 terms
     +       -5*e4/2*(c1*corr(8) +c3*corr(9))      !7/2 terms
c      corrud=pi*corrud; corrud=-corrud/pi8
      corrud=-corrud/8
c      print*,corrud,corr(1),corr(2),corr(3),c1 !corr(4),corr(5) !c1*corr(3),c3*corr(4),c5*corr(5)

c vdd at r=0
      c1= -8*f2b*rdb*xib**4
      c3=  8*f2b*rdb**3*xib**2 - 4*rdb*(3*f1db*rdb+4*f2b*zd2b)*xib**3
     +    -4*(f2d2b*rdb+f2b*rd3b/3)*xib**4
      c5= f2b*rdb**5+rdb**3*(3*f1db*rdb+8*f2b*zd2b)*xib
     +    +2*rdb*(   2*f2d2b*rdb**2+2*f2b*rd3b
     +              -9*f1db*rdb*zd2b-6*f2b*zd2b**2   )*xib**2
      c7=rdb**3/6*( 3*f2d2b*rdb**2  +5*f2b*rdb*rd3b
     +             +9*f1db*rdb*zd2b +12*f2b*zd2b**2)
c     +    + rdb/6*(3*f1d3b*rdb**3 + 6*f1db*rdb*(2rdb*rd3b-9*zd2b**2) 
c     +    + 4*(6*f2d2b*rdb**2*zd2b + 6*f2b*rdb*rd3b*zd2b 
c     +    - 6*f2b*zd2b**3+f2b*rdb**2*zd4b))*xib
    
      corrvdd=  c1*corr(6) +c3*corr(7) +c5*corr(8) +c7*corr(9) !7/2 terms
     +    -7*e4/2*(c1*corr(12) +c3*corr(13)+c5*corr(17)) !9/2 terms
c      corrvdd=pi*corrvdd; corrvdd=-corrvdd/pi8
      corrvdd=-corrvdd/8

      return
      end


      SUBROUTINE correctionmendsh(alf,ns,z0,zb,d,alfp,corrv,
     +    rdb,rd3b,zd2b,zd4b,f1db,f2b,f2d2b)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     short version of above to return corrv only
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,ns
C     INPUT/OUTPUT
      MODE z0,zb,d,alfp,corrv,alf(0:*)
      MODE rdb,rd3b,zd2b,zd4b,f1db,f2b,f2d2b
C     LOCAL
      MODE cmsq,cm,e4,corr(20),corrl(20),xib,c1,c3,c5

      xib=zb-z0
      cmsq=rdb**2+xib*zd2b
      cm=sqrt(cmsq)
      e4 = rdb*rd3b/3 + zd2b**2/4  +xib*zd4b/12

      call compallcorrend(alf,alfp,d,cm,ns,corr,corrl)

c v at r=0
      c1= 2*xib**2*f2b*rdb
      c3= f2b*rdb**3 + rdb*(f1db*rdb+2*f2b*zd2b)*xib
     +         + (3*f2d2b*rdb+f2b*rd3b)*xib**2/3
      c5=(rdb*(f2d2b*rdb**2+f1db*rdb*zd2b+f2b*(rdb*rd3b+zd2b**2)))/2
      corrv= c1*corrl(1)+c3*corrl(2) +c5*corrl(19) !3/2 terms
     +      -3*e4/2*(c1*corr(5) +c3*corr(16))   !5/2 terms
c      corrv=2*pi*corrv; corrv=-corrv/pi8
      corrv=-corrv/4

      return
      end


      SUBROUTINE correctionqend(alf,ns,z0,zb,d,alfp,corrv,
     +    corrud,corrvdd,rdb,rd3b,zd2b,zd4b,udb,ud3b,vb,vd2b)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Input: alfp and all derivatives at alfp 
c            the point (r=0,z0) off sheet at which velo is to be evaluated
c            the normal projection (0,zb) of that point onto sheet 
c            the distance d from (0,z0) to (0,zb)
c
c     SUBSCRIPT B FOR BASE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,ns
C     INPUT/OUTPUT
      MODE z0,zb,d,alfp,corrv,corrud,corrvdd,alf(0:*)
      MODE rdb,rd3b,zd2b,zd4b,udb,ud3b,vb,vd2b
C     LOCAL
      MODE cmsq,cm,corr(20),corrl(20),xib,c1,c3,c5,c7,c9,e4
C     COMMON
c      common/alfpendvals/rdb,rd3b,zd2b,udb,vb,vd2b

      xib=zb-z0
      cmsq=rdb**2+xib*zd2b
      cm=sqrt(cmsq)
      e4 = rdb*rd3b/3 + zd2b**2/4  +xib*zd4b/12

      call compallcorrend(alf,alfp,d,cm,ns,corr,corrl)

c v at r=0 
      c1= 12*xib**3*vb*rdb**2
      c3= 6*xib**2*rdb**2*(2*udb*rdb+vb*zd2b)
     +   +2*xib**3*rdb*(4*vb*rd3b+3*vd2b*rdb)
      c5= -3*vb*rdb**2*zd2b**2*xib
     +  +rdb/2*(4*rdb**2*ud3b+16*rd3b*vb*zd2b
     +          +rdb*(20*rd3b*udb+6*vd2b*zd2b-vb*zd4b))*xib**2
      c7=-3*rdb**2*zd2b**2/2*(2*rdb*udb+vb*zd2b)

      corrv= c1*corrl(3)+c3*corrl(4)+c5*corrl(5) +c7*corrl(16) !5/2 terms
     +    -2.5q0*e4*(c1*corrl(8)+c3*corrl(9))                  !7/2 terms
c      corrv=pi*corrv/pi8
      corrv=corrv/8

c ud at r=0 
      c1= -12*rdb**2*vb*xib**4
      c3=  18*rdb**4*vb*xib**2
     +    -6*rdb**2*(3*rdb*udb+vb*zd2b)*xib**3
     +    -2*rdb*(4*rd3b*vb+3*rdb*vd2b)*xib**4
      c5=  6*rdb**4*(2*rdb*udb+vb*zd2b)*xib
     +    +3*rdb**2*(6*rdb*rd3b*vb+3*rdb**2*vd2b
     +                -rdb*udb*zd2b+3*vb*zd2b**2)*xib**2
     +    -rdb*(3*rdb**2*ud3b+10*rd3b*vb*zd2b
     +          +rdb*(15*rd3b*udb+3*vd2b*zd2b-vb*zd4b))*xib**3
      c7= -3*rdb**4*vb*zd2b**2/2 
     +    +rdb**2/2*(4*rdb**3*ud3b+15*vb*zd2b**3
     +    +rdb*zd2b*(20*rd3b*vb+21*udb*zd2b)
     +    +rdb**2*(28*rd3b*udb+6*vd2b*zd2b-vb*zd4b))*xib

      corrud= c1*corr(6)+c3*corr(7)+c5*corr(8)+c7*corr(9)  !7/2 terms
     +    -3.5q0*e4*(c1*corr(12)+c3*corr(13)+c5*corr(17))              !9/2 terms
c      corrud=pi*corrud/pi8
      corrud=corrud/8

c vdd at r=0
      c1= -60*rdb**2*vb*xib**5
      c3= 150*rdb**4*vb*xib**3-30*rdb**2*(4*rdb*udb+vb*zd2b)*xib**4
     +     -2*rdb*(20*rd3b*vb+15*rdb*vd2b+6*udb*zd2b)*xib**5
      c5=  45*rdb**4*(2*rdb*udb+3*vb*zd2b)*xib**2 
     +     +3*rdb**2*( 50*rdb*rd3b*vb+25*rdb**2*vd2b
     +                -28*rdb*udb*zd2b+30*vb*zd2b**2)*xib**3
     +     -5*rdb/2*(8*rdb**2*ud3b +12*zd2b*(2*rd3b*vb +udb*zd2b)
     +              +rdb*(40*rd3b*udb +6*vd2b*zd2b -3*vb*zd4b))*xib**4
      c7=   3*rdb**4*zd2b/2*(32*rdb*udb+15*vb*zd2b)*xib
     +     +3*rdb**2/4*(20*rdb**3*ud3b +140*vb*zd2b**3
     +                  +4*rdb*zd2b*(55*rd3b*vb+18*udb*zd2b)
     +   +5*rdb**2*(28*rd3b*udb+18*vd2b*zd2b+vb*zd4b))*xib**2

c     +   +rdb**2/4*(165*vb*zd2b**4+6rdb*zd2b**2*(35*rd3b*vb+38*udb*zd2b)
c     +              +2*rdb**3*(16*ud3b*zd2b+udb*zd4b)
c     +              +5*rdb**2*zd2b*(56*rd3b*udb+9*vd2b*zd2b-3*vb*zd4b)
c     +             )*xib
      c9=3*rdb**4*zd2b**2/4*(2*rdb*udb-5*vb*zd2b)

      corrvdd= c1*corr(10)+c3*corr(11)+c5*corr(12)  !9/2 terms
     +        +c7*corr(13)+c9*corr(17)              !9/2 terms
     +    -4.5q0*e4*(c1*corr(14)+c3*corr(15)+c5*corr(18)) !11/2 terms
c      corrvdd=pi*corrvdd/pi8
      corrvdd=corrvdd/8
      return
      end


      SUBROUTINE correctionqendsh(alf,ns,z0,zb,d,alfp,corrv,
     +    rdb,rd3b,zd2b,zd4b,udb,ud3b,vb,vd2b)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Input: alfp and all derivatives at alfp 
c            the point (r=0,z0) off sheet at which velo is to be evaluated
c            the normal projection (0,zb) of that point onto sheet 
c            the distance d from (0,z0) to (0,zb)
c
c     SUBSCRIPT B FOR BASE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,ns
C     INPUT/OUTPUT
      MODE z0,zb,d,alfp,corrv,alf(0:*)
      MODE rdb,rd3b,zd2b,zd4b,udb,ud3b,vb,vd2b
C     LOCAL
      MODE cmsq,cm,corr(20),corrl(20),xib,c1,c3,c5,c7,e4
C     COMMON
c      common/alfpendvals/rdb,rd3b,zd2b,udb,vb,vd2b

      xib=zb-z0
      cmsq=rdb**2+xib*zd2b
      cm=sqrt(cmsq)
      e4 = rdb*rd3b/3 + zd2b**2/4  +xib*zd4b/12

      call compallcorrend(alf,alfp,d,cm,ns,corr,corrl)

c v at r=0 
      c1= 12*xib**3*vb*rdb**2
      c3= 6*xib**2*rdb**2*(2*udb*rdb+vb*zd2b)
     +   +2*xib**3*rdb*(4*vb*rd3b+3*vd2b*rdb)
      c5= -3*vb*rdb**2*zd2b**2*xib
     +  +rdb/2*(4*rdb**2*ud3b+16*rd3b*vb*zd2b
     +          +rdb*(20*rd3b*udb+6*vd2b*zd2b-vb*zd4b))*xib**2
      c7=-3*rdb**2*zd2b**2/2*(2*rdb*udb+vb*zd2b)

      corrv= c1*corrl(3)+c3*corrl(4)+c5*corrl(5) +c7*corrl(16) !5/2 terms
     +    -2.5q0*e4*(c1*corrl(8)+c3*corrl(9))                  !7/2 terms
c      corrv=pi*corrv/pi8
      corrv=corrv/8

      return
      end
