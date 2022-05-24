#include "flags.h"
      SUBROUTINE compaxism(alf,z0,rs,zs,ns,f1,f2,
     +  rd0,rdn,rddd0,rdddn,zdd0,zddn,f1d0,f1dn,f2dd0,f2ddn,v,ud,vdd)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Input: 
c         z0 - value of z of target point on axis (r=0) 
c               at which to evaluate all velocities
c         alf - parameter for the sheet
c         rs,zs,ns,f1,f2 - description of sheet
c         rd0,rddd0,zdd0,f1d0,f2dd0 - values of derivs at left endpt
c         rdn,rdddn,zddn,f1dn,f2ddn - values of derivs at right endpt
c      Output:
c         v,ud,vdd - values of axial and radial velo at (0,z0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns
      real*8 alf(0:*),z0,rs(0:*),zs(0:*),f1(0:*),f2(0:*),
     +      rd0,rdn,rddd0,rdddn,zdd0,zddn,f1d0,f1dn,f2dd0,f2ddn
c OUTPUT
      real*8 v,ud,vdd
c LOCAL
      integer j,k,n
      real*8 xi,xisq,r,rsq,fact,alfb,zd4b,
     +     left,right,delalf,corrv,corrud,corrvdd
      real*8 right3,left3,temp,distl,distr

      delalf = pi/ns

c     COMPUTE v 
      v=0
      do k = 1,ns-1      
         xi   = zs(k) - z0
         rsq  = rs(k)**2
         fact = rs(k)/dsqrt(rsq + xi**2)**3
         v= v + fact*( f1(k)*rs(k)*xi + f2(k)*(2*xi**2+rsq) )
      enddo
      v = v*pi2*delalf 
c     Add O(h^2) contrib from 1st deriv
c     Add O(h^4) contrib from 3rd deriv on 12/30/20, see Msendderiv.nb
      xi   = zs(0) - z0
      left= 4*pi*f2(0)*rd0/abs(xi)
      temp=rd0/abs(xi)
      left3=pi2*( f2dd0*6*temp + f1d0*6*temp**2 +2*f2(0)*
     +   ( -6*temp**3+rddd0/abs(xi)-3*zdd0*temp/abs(xi) ))
      xi   = zs(ns) - z0
      right= 4*pi*f2(ns)*rdn/abs(xi)
      temp=rdn/abs(xi)
      right3=pi2*( f2ddn*6*temp + f1dn*6*temp**2 +2*f2(ns)*
     +   ( -6*temp**3+rdddn/abs(xi)-3*zddn*temp/abs(xi) ))
      v= v-delalf**2/12*(right-left)+delalf**4/720*(right3-left3)
c      v= v-delalf**2/12*(right-left)
      v=-v/pi8

c     COMPUTE ud
      ud=0
      do k = 1,ns-1      
         xi   = zs(k) - z0
         rsq  = rs(k)**2
         fact = rs(k)/dsqrt(rsq + xi**2)**5
         ud= ud + fact*( f1(k)*rs(k)+f2(k)*xi )*( rsq-2*xi**2 ) 
      enddo
      ud = ud*pi*delalf 
c     Add O(h^2) contrib from 1st deriv
c     SKIP O(h^4). For 3rd deriv,  see Msendderiv.nb ??
      xi   = zs(0) - z0
      left= -2*pi*f2(0)*rd0/(abs(xi)*xi)
c      left3=...
      xi   = zs(ns) - z0
      right= -2*pi*f2(ns)*rdn/(abs(xi)*xi)
c      right3=...
      ud= ud-delalf**2/12*(right-left)
      ud=-ud/pi8  !MINUS SIGN NOT CORRECT NEED TO FIX EVERYWHERE

c     COMPUTE vdd
      vdd=0
      do k = 1,ns-1      
         xi   = zs(k) - z0
         xisq = xi**2
         r    = rs(k)
         rsq  = rs(k)**2
         fact = r/dsqrt(rsq + xisq)**7
         vdd= vdd + fact*(3*f1(k)*r*xi*(rsq-4*xisq)
     +                     +f2(k)*(rsq**2+8*rsq*xisq-8*xisq**2 ) )
      enddo
      vdd = vdd*pi*delalf 
c     Add  O(h^2) contrib from 1st deriv, skip O(h^4) contrib from 3rd deriv
      xi   = zs(0) - z0
      xisq = xi**2
      left= -8*pi*f2(0)*rd0/(abs(xi)*xisq)
      xi   = zs(ns) - z0
      xisq = xi**2
      right= -8*pi*f2(ns)*rdn/(abs(xi)*xisq)
      vdd= vdd-delalf**2/12*(right-left)
      vdd=-vdd/pi8

      distl=abs(z0-zs(0))
      distr=abs(z0-zs(ns))
      if (distl.le.20*delalf) then
              alfb=0
        zd4b=2*(zs(2)-4*zs(1)+3*zs(0))/delalf**4  !z even
        call correctionmend(alf,ns,z0,zs(0),distl,alfb,corrv,corrud,
     +        corrvdd,rd0,rddd0,zdd0,zd4b,f1d0,f2(0),f2dd0)
        v   =   v + corrv
        ud  =  ud + corrud
        vdd = vdd + corrvdd
      endif
      if (distr.le.5*delalf) then
              alfb=pi
        zd4b=2*(zs(ns-1)-4*zs(ns-1)+3*zs(ns))/delalf**4  !z even
        call correctionmend(alf,ns,z0,zs(ns),distr,alfb,
     +       corrv,corrud,corrvdd,rdn,rdddn,zddn,zd4b,f1dn,f2(ns),f2ddn)
        v   =   v + corrv
        ud  =  ud + corrud
        vdd = vdd + corrvdd
      endif

      return
      end


      SUBROUTINE compaxismsh(alf,z0,rs,zs,ns,f1,f2,
     +  rd0,rdn,rddd0,rdddn,zdd0,zddn,f1d0,f1dn,f2dd0,f2ddn,v)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c returns only v on axis
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns
      real*8 alf(0:*),z0,rs(0:*),zs(0:*),f1(0:*),f2(0:*),
     +      rd0,rdn,rddd0,rdddn,zdd0,zddn,f1d0,f1dn,f2dd0,f2ddn
c OUTPUT
      real*8 v
c LOCAL
      integer j,k,n
      real*8 xi,xisq,r,rsq,fact,alfb,zd4b,
     +     left,right,delalf,corrv
      real*8 right3,left3,temp,distl,distr

      delalf = pi/ns

c     COMPUTE v 
      v=0
      do k = 1,ns-1      
         xi   = zs(k) - z0
         rsq  = rs(k)**2
         fact = rs(k)/dsqrt(rsq + xi**2)**3
         v= v + fact*( f1(k)*rs(k)*xi + f2(k)*(2*xi**2+rsq) )
      enddo
      v = v*pi2*delalf 
c     Add O(h^2) contrib from 1st deriv
c     Add O(h^4) contrib from 3rd deriv on 12/30/20, see Msendderiv.nb
      xi   = zs(0) - z0
      left= 4*pi*f2(0)*rd0/abs(xi)
      temp=rd0/abs(xi)
      left3=pi2*( f2dd0*6*temp + f1d0*6*temp**2 +2*f2(0)*
     +   ( -6*temp**3+rddd0/abs(xi)-3*zdd0*temp/abs(xi) ))
      xi   = zs(ns) - z0
      right= 4*pi*f2(ns)*rdn/abs(xi)
      temp=rdn/abs(xi)
      right3=pi2*( f2ddn*6*temp + f1dn*6*temp**2 +2*f2(ns)*
     +   ( -6*temp**3+rdddn/abs(xi)-3*zddn*temp/abs(xi) ))
      v= v-delalf**2/12*(right-left)+delalf**4/720*(right3-left3)
c      v= v-delalf**2/12*(right-left)
      v=-v/pi8


      distl=abs(z0-zs(0))
      distr=abs(z0-zs(ns))
      if (distl.le.20*delalf) then
              alfb=0
        zd4b=2*(zs(2)-4*zs(1)+3*zs(0))/delalf**4  !z even
        call correctionmendsh(alf,ns,z0,zs(0),distl,alfb,
     +      corrv,rd0,rddd0,zdd0,zd4b,f1d0,f2(0),f2dd0)
        v   =   v + corrv
      endif
      if (distr.le.5*delalf) then
              alfb=pi
        call correctionmendsh(alf,ns,z0,zs(ns),distr,alfb,
     +       corrv,rdn,rdddn,zddn,f1dn,f2(ns),f2ddn)
        v   =   v + corrv
      endif

      return
      end



      SUBROUTINE compaxisq(alf,zz,rs,zs,ns,u,v,rd,zd,
     +  rddd0,rdddn,zdd0,zddn,ud0,udn,vdd0,vddn,
     +  vax,udax,vddax)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns
      real*8 alf(0:*),zz,rs(0:*),zs(0:*),u(0:*),v(0:*),
     +      rd(0:*),zd(0:*),
     +      rddd0,rdddn,zdd0,zddn,ud0,udn,vdd0,vddn
c OUTPUT
      real*8 vax,udax,vddax
c LOCAL
      integer j,k,n
      MODE xi,xisq,r,rsq,fact,alfb,
     +     left,right,delalf,corrv,corrud,corrvdd
      MODE right3,left3,temp,distl,distr,ud3b,zd4b

      delalf = pi/ns

c     COMPUTE v 
      vax=0

      do k = 1,ns-1
         xi   = zs(k) - zz
         xisq = xi**2
         r=rs(k)
         rsq  = r**2
         fact = r*xi/dsqrt(rsq + xisq)**5
         vax= vax+fact*(u(k)*r+v(k)*xi)*(r*zd(k)-xi*rd(k))
      enddo
      vax = -12*pi*vax*delalf

c     Add O(h^2,h^4) contrib from 1st deriv at endpoints of S1
      xi   = zs(0) - zz
      left = 12*pi*v(0)*rd(0)**2/(xi*abs(xi))
      left3 = pi2*( 6*v(0)*rd(0)*(-15*rd(0)**3+4*rddd0*xi**2
     +                            -12*xi*rd(0)*zdd0)/(xi**3*abs(xi))
     +       +36*ud0*rd(0)**3/(xi**2*abs(xi))
     +       +18*vdd0*rd(0)**2/(xi*abs(xi)) )
      xi   = zs(ns) - zz
      right= 12*pi*v(ns)*rd(ns)**2/(xi*abs(xi))
      right3 = pi2*( 6*v(ns)*rd(ns)*(-15*rd(ns)**3+4*rdddn*xi**2
     +                         -12*xi*rd(ns)*zddn)/(xi**3*abs(xi))
     +       +36*udn*rd(ns)**3/(xi**2*abs(xi))
     +       +18*vddn*rd(ns)**2/(xi*abs(xi)) )
      vax= vax -delalf**2/12*(right-left)
     + +delalf**4/720*(right3-left3)
      vax=vax/pi8

c     COMPUTE ud
      udax=0
      do k = 1,ns-1
         xi   = zs(k) - zz
         xisq = xi**2
         r=rs(k)
         rsq  = r**2
         fact = r/dsqrt(rsq + xisq)**7
         udax = udax+fact*(
     +         zd(k)*r* (u(k)*r*(  rsq-4*xisq)+v(k)*xi*(2*rsq-3*xisq))
     +        -rd(k)*xi*(u(k)*r*(2*rsq-3*xisq)+v(k)*xi*(3*rsq-2*xisq)) )
      enddo
      udax = -6*pi*udax*delalf

c     Add O(h^2,h^4) contrib from 1st deriv at endpoints of S1
      xi    = zs(0) - zz
      left  = -12*pi*v(0)*rd(0)**2/(xi**2*abs(xi))
      xi    = zs(ns) - zz
      right = -12*pi*v(ns)*rd(ns)**2/(xi**2*abs(xi))
      udax= udax -delalf**2/12*(right-left)
      udax= udax/pi8

c     COMPUTE vdd
      vddax=0
      do k = 1,ns-1
         xi   = zs(k) - zz
         xisq = xi**2
         r=rs(k)
         rsq  = r**2
         fact = r*xi/dsqrt(rsq + xisq)**9
         vddax= vddax+fact*(
     +        zd(k)*(7*r**4*u(k)+15*r**3*v(k)*xi-26*rsq*u(k)*xisq
     +                                   -20*r*v(k)*xi**3+2*u(k)*xi**4)
     +     -5*rd(k)*xi*(3*r**3*u(k)+5*rsq*v(k)*xi
     +                                   -4*r*u(k)*xisq-2*v(k)*xi**3) 
     +  )
      enddo
      vddax = -6*pi*vddax*delalf

c     Add O(h^2,h^4) contrib from 1st deriv at endpoints of S1
      xi   = zs(0) - zz
      left = -60*pi*v(0)*rd(0)**2/(xi**3*abs(xi))
      xi   = zs(ns) - zz
      right= -60*pi*v(ns)*rd(ns)**2/(xi**3*abs(xi))
      vddax= vddax -delalf**2/12*(right-left)
      vddax=vddax/pi8

      distl=abs(zz-zs(0))
      distr=abs(zz-zs(ns))
      if (distl.le.60*delalf) then
        alfb=0
        ud3b=(u(2)-2*u(1))/delalf**3           !u odd
        zd4b=2*(zs(2)-4*zs(1)+3*zs(0))/delalf**4  !z even

        call correctionqend(alf,ns,zz,zs(0),distl,alfb,corrv,corrud,
     +             corrvdd,rd(0),rddd0,zdd0,zd4b,ud0,ud3b,v(0),vdd0)
        vax   =   vax + corrv
        udax  =  udax + corrud
        vddax = vddax + corrvdd
      endif
      if (distr.le.5*delalf) then
              alfb=pi
        call correctionqend(alf,ns,zz,zs(ns),distr,alfb,
     +       corrv,corrud,corrvdd,rd(ns),rdddn,zddn,udn,v(ns),vddn)
        vax   =   vax + corrv
        udax  =  udax + corrud
        vddax = vddax + corrvdd
      endif

      return
      end


      SUBROUTINE compaxisqsh(alf,zz,rs,zs,ns,u,v,rd,zd,
     +  rddd0,rdddn,zdd0,zddn,ud0,udn,vdd0,vddn,vax)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Computes only axial velocity at r=0,z=zz
c     Throughout code:  u=dr/dt  radial velocity odd
c                       v=dz/dt  axial velocity  even
c     (unlike paper)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns
      real*8 alf(0:*),zz,rs(0:*),zs(0:*),u(0:*),v(0:*),
     +      rd(0:*),zd(0:*),
     +      rddd0,rdddn,zdd0,zddn,ud0,udn,vdd0,vddn
c OUTPUT
      real*8 vax
c LOCAL
      integer j,k,n
      real*8 xi,xisq,r,rsq,fact,alfb,
     +     left,right,delalf,corrv
      real*8 right3,left3,temp,distl,distr,ud3b,zd4b

      delalf = pi/ns

c     COMPUTE v 
      vax=0

      do k = 1,ns-1
         xi   = zs(k) - zz
         xisq = xi**2
         r=rs(k)
         rsq  = r**2
         fact = r*xi/dsqrt(rsq + xisq)**5
         vax= vax+fact*(u(k)*r+v(k)*xi)*(r*zd(k)-xi*rd(k))
      enddo
      vax = -12*pi*vax*delalf

c     Add O(h^2,h^4) contrib from 1st,3rd deriv at endpoints of S1
      xi   = zs(0) - zz
      left = 12*pi*v(0)*rd(0)**2/(xi*abs(xi))
      left3 = pi2*( 6*v(0)*rd(0)*(-15*rd(0)**3+4*rddd0*xi**2
     +                            -12*xi*rd(0)*zdd0)/(xi**3*abs(xi))
     +       +36*ud0*rd(0)**3/(xi**2*abs(xi))
     +       +18*vdd0*rd(0)**2/(xi*abs(xi)) )
      xi   = zs(ns) - zz
      right= 12*pi*v(ns)*rd(ns)**2/(xi*abs(xi))
      right3 = pi2*( 6*v(ns)*rd(ns)*(-15*rd(ns)**3+4*rdddn*xi**2
     +                         -12*xi*rd(ns)*zddn)/(xi**3*abs(xi))
     +       +36*udn*rd(ns)**3/(xi**2*abs(xi))
     +       +18*vddn*rd(ns)**2/(xi*abs(xi)) )
      vax= vax -delalf**2/12*(right-left)
     + +delalf**4/720*(right3-left3)
      vax=vax/pi8

      distl=abs(zz-zs(0))
      distr=abs(zz-zs(ns))
      if (distl.le.60*delalf) then
        alfb=0
        ud3b=(u(2)-2*u(1))/delalf**3           !u odd
        zd4b=2*(zs(2)-4*zs(1)+3*zs(0))/delalf**4  !z even
        call correctionqendsh(alf,ns,zz,zs(0),distl,alfb,corrv,
     +             rd(0),rddd0,zdd0,zd4b,ud0,ud3b,v(0),vdd0)

        vax   =   vax + corrv
      endif
      if (distr.le.5*delalf) then
              stop 'fix it'
              alfb=pi
c        call correctionqendsh(alf,ns,zz,zs(ns),distr,alfb,corrv,
c     +           rd(ns),rdddn,zddn,zd4n,udn,ud3b,v(ns),vddn)
        vax   =   vax + corrv
      endif

      return
      end
