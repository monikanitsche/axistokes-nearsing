Contains:
c SUBROUTINE compendsf: computes contributions from 1st and 3rd
c                       derivatives to T6[G], G=Gus,Gvs
c SUBROUTINE compenddv: computes contributions from 1st and 3rd
c                       derivatives to T6[G], G=Gud,Gvd

#include "flags.h"
      SUBROUTINE compendsf(xx,yy,
     +    xd0,xdn,xddd0,xdddn,ys0,ysn,ydd0,yddn,
     +    f1d0,f1dn,f20,f2n,f2dd0,f2ddn,
     +    corrus1,corrvs1,corrus2,corrvs2,delalf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Computes the correction for the single layer velocities introduced
c     by the first third derivatives
c             corrus1 = - h^2/12[Bp-Gp]
c             corrus2 = + h^4/720[Bppp-Gppp]
c     Note 12/22/20: this should be instead
c             corrus1 = - h^2/12[Gp]
c             corrus2 = + h^4/720[Gppp]
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c  #include "constants.h"
c INPUT/OUTPUT
      MODE xx,yy,
     +    xd0,xdn,xddd0,xdddn,ys0,ysn,ydd0,yddn,
     +    f1d0,f1dn,f20,f2n,f2dd0,f2ddn,
     +    corrus1,corrvs1,corrus2,corrvs2,delalf
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      MODE x2,x4,xi,xi2,xi4,root,root2,gpul,gpur,gpvl,gpvr,left,right,
     +   gpppul,gpppur,gpppvl,gpppvr,m11pp,m12p,m12ppp,m21pp,m22p,m22ppp

      x2 = xx**2
      x4 = x2**2

c LEFT ENDPOINT
      xi = ys0-yy
      xi2 = xi**2
      xi4 = xi**4
      root2 = x2+xi2
      root = sqrt(root2)
c      gpul = -pi2*f20 *xd0 *xx*xi/root**3 
c      gpvl =  pi2*f20 *xd0 *(x2+2*xi2)/root**3 
c  SEE bubbles/math/Msderiv.nb
c     f1 odd => f1=f1pp=0
c     f2 even => f2p=f2ppp=0
c     m11=0, m11p=0, m11ppp=0 => (m11f1)'  =m11f1'=0
c                                (m11f1)'''=2m11''f1'+m11f1'''=2m11''f1'
c     m12=0, m12pp=0  => (m12f2)'  = m12'f2
c                        (m12f2)'''= m12'''f2+3m12'f2''
      m11pp = pi2*xd0**2*xx*(x2-2*xi2)/root**5
      m12p  =-pi2*xd0*xx*xi/root**3
      m12ppp=-pi*xx*( 9*xd0**3*(x2-4*xi2)*xi 
     +  + 2*xddd0*xi*root2**2
     +  + 6*xd0*(x4-x2*xi2-2*xi4)*ydd0 )/root**7
      gpul = m12p*f20 
      gpppul=3*m11pp*f1d0 + m12ppp*f20 + 3*m12p*f2dd0
c     m21=0, m21p=0, m21ppp=0 => (m21f1)'  =m21f1'=0
c                                (m21f1)'''=2m21''f1'+m21f1'''=2m21''f1'
c     m22=0, m22pp=0  => (m22f2)'  = m22'f2
c                        (m22f2)'''= m22'''f2+3m22'f2''
      m21pp =-pi2*xd0**2*xi*(x2-2*xi2)/root**5
      m22p  = pi2*xd0*(x2+2*xi2)/root**3
      m22ppp= pi*(
     +    2*xddd0*root2**2*(x2+2*xi2)
     +  + 3*xd0**3*(x4+8*x2*xi2-8*xi4)
     +  + 6*xd0*xi*(x4-x2*xi2-2*xi4)*ydd0 )/root**7
      gpvl = m22p*f20 
      gpppvl=3*m21pp*f1d0 + m22ppp*f20 + 3*m22p*f2dd0


c RIGHT ENDPOINT

      xi = ysn-yy
      xi2 = xi**2
      xi4 = xi**4
      root2 = x2+xi2
      root = sqrt(root2)
c      gpur = -pi2*f2n *xd(ns) *xx*xi/root**3 
c      gpvr =  pi2*f2n *xd(ns) *(x2+2*xi2)/root**3  
c  SEE bubbles/math/Msderiv.nb
c     f1 odd => f1=f1pp=0
c     f2 even => f2p=f2ppp=0
c     m11=0, m11p=0, m11ppp=0 => (m11f1)'  =m11f1'=0
c                                (m11f1)'''=2m11''f1'+m11f1'''=2m11''f1'
c     m12=0, m12pp=0  => (m12f2)'  = m12'f2
c                        (m12f2)'''= m12'''f2+3m12'f2''
      m11pp = pi2*xdn**2*xx*(x2-2*xi2)/root**5
      m12p  =-pi2*xdn*xx*xi/root**3
      m12ppp=-pi*xx*( 9*xdn**3*(x2-4*xi2)*xi 
     +  + 2*xdddn*xi*root2**2
     +  + 6*xdn*(x4-x2*xi2-2*xi4)*yddn )/root**7
      gpur = m12p*f2n 
      gpppur=3*m11pp*f1dn + m12ppp*f2n + 3*m12p*f2ddn
c     m21=0, m21p=0, m21ppp=0 => (m21f1)'  =m21f1'=0
c                                (m21f1)'''=2m21''f1'+m21f1'''=2m21''f1'
c     m22=0, m22pp=0  => (m22f2)'  = m22'f2
c                        (m22f2)'''= m22'''f2+3m22'f2''
      m21pp =-pi2*xdn**2*xi*(x2-2*xi2)/root**5
      m22p  = pi2*xdn*(x2+2*xi2)/root**3
      m22ppp= pi*(
     +    2*xdddn*root2**2*(x2+2*xi2)
     +  + 3*xdn**3*(x4+8*x2*xi2-8*xi4)
     +  + 6*xdn*xi*(x4-x2*xi2-2*xi4)*yddn )/root**7
      gpvr = m22p*f2n 
      gpppvr=3*m21pp*f1dn + m22ppp*f2n + 3*m22p*f2ddn

      left  = gpul
      right = gpur
      corrus1 = -delalf**2/12*(right-left)
      left  = gpvl
      right = gpvr
      corrvs1 = -delalf**2/12*(right-left)


      left  = gpppul
      right = gpppur
      corrus2 = delalf**4/720*(right-left)
      left  = gpppvl
      right = gpppvr
      corrvs2 = delalf**4/720*(right-left)

      return
      end


      SUBROUTINE compenddv(xx,yy,delalf,
     +   xd0,xddd0,y0,ydd0,ud0,v0,vdd0,
     +   xdn,xdddn,yn,yddn,udn,vn,vddn,
     +   corrud1,corrvd1,corrud2,corrvd2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Computes the correction for the double layer velocities introduced
c     by the first third derivatives
c             corrud1 = - h^2/12[Bp-Gp]
c             corrud2 = + h^4/720[Bppp-Gppp]
c    Note 12/22/20: This should be, instead, check today
c             corrud1 = - h^2/12[Gp]
c             corrud2 = + h^4/720[Gppp]
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j
      MODE xx,yy,delalf,
     +   xd0,xddd0,y0,ydd0,ud0,v0,vdd0,
     +   xdn,xdddn,yn,yddn,udn,vn,vddn,
     +   corrud1,corrvd1,corrud2,corrvd2
      MODE xi,r,rsq,fact,gpul,gpur,gpvl,gpvr,gpppul,gpppur,
     +    gpppvl,gpppvr,left,right,f1,f2,f3

      xi = y0-yy
      rsq = xx**2+xi**2
      r = sqrt(rsq)
      fact = 12*pi*v0*xd0**2/r**5   !Q22
      gpul = -fact*xi**2*xx
      gpvl =  fact*xi**3

C     FOR BOTH U AND V COMPONENTS
      f1 = 2*( 3*xx**2*ud0*ydd0 -            !Q11+Q22
     +     xi**2*(3*xd0*vdd0+4*xddd0*v0) )
c      f4 = 6*xd0*xi*2*xx**2/rsq   
c     +               *(xd0*ud0-2*v0*ydd0)
c      f2 = 6*xd0*xi*(-3*xi**2)/rsq    
c     +               *(xd0*ud0-2*v0*ydd0)
      f2 = 6*xd0*xi*(2*xx**2-3*xi**2)/rsq    !Q12=Q21
     +               *(xd0*ud0-2*v0*ydd0)
      f3 = 15*xd0**3*xi**2*v0*(-3*xx**2+4*xi**2)/rsq**2 !Q22

      gpppul = 6*pi*xx*xd0/r**5*(f1+f2+f3)
      f1 = 2*( -3*xx**2*ud0*ydd0 +
     +     xi**2*(3*xd0*vdd0+4*xddd0*v0) )
c      f4 = 6*xd0*xi*(3*xx**2)/rsq
c     +               *(-xd0*ud0+2*v0*ydd0)
c      f2 = 6*xd0*xi*(-2*xi**2)/rsq
c     +               *(-xd0*ud0+2*v0*ydd0)
      f2 = 6*xd0*xi*(3*xx**2-2*xi**2)/rsq
     +               *(-xd0*ud0+2*v0*ydd0)
      f3 = 15*xd0**3*xi**2*v0*(5*xx**2-2*xi**2)/rsq**2
      gpppvl = 6*pi*xi*xd0/r**5*(f1+f2+f3)

      xi = yn-yy
      rsq = xx**2+xi**2
      r = sqrt(rsq)
      fact = 12*pi*vn*xdn**2/r**5
      gpur = -fact*xi**2*xx
      gpvr =  fact*xi**3

C     FOR BOTH U AND V COMPONENTS
      f1 = 2*( 3*xx**2*udn*yddn -
     +     xi**2*(3*xdn*vddn+4*xdddn*vn) )
c      f4 = 6*xdn*xi*(2*xx**2)/rsq
c     +               *(xdn*udn-2*vn*yddn)
c      f2 = 6*xdn*xi*(-3*xi**2)/rsq
c     +               *(xdn*udn-2*vn*yddn)
      f2 = 6*xdn*xi*(2*xx**2-3*xi**2)/rsq
     +               *(xdn*udn-2*vn*yddn)
      f3 = 15*xdn**3*xi**2*vn*(-3*xx**2+4*xi**2)/rsq**2

      gpppur = 6*pi*xx*xdn/r**5*(f1+f2+f3)

      f1 = 2*( -3*xx**2*udn*yddn +
     +     xi**2*(3*xdn*vddn+4*xdddn*vn) )
c      f4 = 6*xdn*xi*(3*xx**2)/rsq
c     +               *(-xdn*udn+2*vn*yddn)
c      f2 = 6*xdn*xi*(-2*xi**2)/rsq
c     +               *(-xdn*udn+2*vn*yddn)
      f2 = 6*xdn*xi*(3*xx**2-2*xi**2)/rsq
     +               *(-xdn*udn+2*vn*yddn)
      f3 = 15*xdn**3*xi**2*vn*(5*xx**2-2*xi**2)/rsq**2
      gpppvr = 6*pi*xi*xdn/r**5*(f1+f2+f3)

      left  = gpul
      right = gpur
      corrud1 = -delalf**2/12*(right-left)
      left  = gpvl
      right = gpvr
      corrvd1 = -delalf**2/12*(right-left)
      left  = gpppul
      right = gpppur
      corrud2 = delalf**4/720*(right-left)
      left  = gpppvl
      right = gpppvr
      corrvd2 = delalf**4/720*(right-left)

      return
      end
