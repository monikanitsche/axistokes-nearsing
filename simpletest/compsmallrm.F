#include "flags.h"
      SUBROUTINE compsmallrm(rr,zz,rmax,alf,rs,zs,ns,f1,f2,
     +           rd0,rd30,zd20,f1d0,f2d20,
     +           rdn,rd3n,zd2n,f1dn,f2d2n,uextrap,vextrap)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Computes u,v for small r=O(h) near axis-intersect by extrapolation
c     Extrapolates along arc parallel to interface.
c     Third order extrapolation obtained using
c     u (odd) :  udax=du/dr
c     v (even):  vax=v vddax=d^2v/dr^2 at r=0 uxax=du/dx at r=0 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer ns
      MODE rr,zz,rd0,rd30,zd20,f1d0,f2d20,rdn,rd3n,zd2n,f1dn,f2d2n
      MODE alf(0:*),rs(0:*),zs(0:*),f1(0:*),f2(0:*)
c OUTPUT
      MODE uextrap,vextrap,vextrap0
c LOCAL
      MODE rad,delalf,theta,zcirc,d,z0,delz,delr,hloc,
     +   zplus,zmin,vaxplus,vaxmin,udaxplus,udaxmin
      MODE vax,udax,vddax,vzax,udzax
      MODE c0,c1,c2,c3,c4,c5,c6,aa,bb,cc,dd,thetmax,thetmaxplus,
     +   umax,umaxplus,udmax,vmax,vmaxplus,vdmax,
     +   rmax,zmax,rmaxplus,zmaxplus,fdd,radtarget,
     + thetmaxmin,vmaxmin,umaxmin,zmaxmin,rmaxmin,radmax

c find horizontal signed distance D to interface =>  z0=zs(0)-D,
c and thus Delta z= zz-z0. Also Delta r=rr

      rad= 0.95q0
      delalf=pi/ns

c     point on circular sheet of radius rad at r=rr has axial comp zcirc
c      theta=asin(rr/rad)
c      zcirc=-4-rad*cos(theta)
c      d=zcirc-zz
c      z0=zs(0)-d
c      radtarget=sqrt((zz+4)**2+rr**2)
c      theta = asin(rr/radtarget)
      radtarget=sqrt((zz+4)**2+rr**2)
      theta=asin(rr/radtarget)
      z0=-4-radtarget
      d=zs(0)-z0

      thetmax = asin(rmax/radtarget)
      rmax=radtarget*sin(thetmax)
      zmax=-4-radtarget*cos(thetmax)

c TESTING
c      call veloforceSonP(rr,zz,alf,rs,zs,ns,f1,f2,
c     +           rd0,rd30,zd20,f1d0,f2d20,
c     +           rdn,rd3n,zd2n,f1dn,f2d2n,umax,vmax)
c      stop
c ENDTESTING

      call veloforceSonP(rmax,zmax,alf,rs,zs,ns,f1,f2,
     +           rd0,rd30,zd20,f1d0,f2d20,
     +           rdn,rd3n,zd2n,f1dn,f2d2n,umax,vmax)

      thetmaxplus=thetmax+delalf/10
      rmaxplus=radtarget*sin(thetmaxplus)
      zmaxplus=-4-radtarget*cos(thetmaxplus)
      call veloforceSonP(rmaxplus,zmaxplus,alf,rs,zs,ns,f1,f2,
     +           rd0,rd30,zd20,f1d0,f2d20,
     +           rdn,rd3n,zd2n,f1dn,f2d2n,umaxplus,vmaxplus)
      thetmaxmin=thetmax-delalf/10
      rmaxmin=radtarget*sin(thetmaxmin)
      zmaxmin=-4-radtarget*cos(thetmaxmin)
      call veloforceSonP(rmaxmin,zmaxmin,alf,rs,zs,ns,f1,f2,
     +           rd0,rd30,zd20,f1d0,f2d20,
     +           rdn,rd3n,zd2n,f1dn,f2d2n,umaxmin,vmaxmin)
      udmax=(umaxplus-umaxmin)/(thetmaxplus-thetmaxmin)
      vdmax=(vmaxplus-vmaxmin)/(thetmaxplus-thetmaxmin)

c find values of derivatives on axis: vax,vddax,udax, and vzax
c here ud,vdd denote derrivatives wrt r on axis
c vz denote derrivatives wrt z on axis

      call compaxism(alf,z0,rs,zs,ns,f1,f2,
     +     rd0,rdn,rd30,rd3n,zd20,zd2n,f1d0,f1dn,f2d20,f2d2n,
     +     vax,udax,vddax)

      d=abs(z0-zs(0))
      hloc=min(delalf/100,d/10)
      zplus=z0+hloc
      zmin=z0-hloc
      call compaxismsh(alf,zplus,rs,zs,ns,f1,f2,
     +     rd0,rdn,rd30,rd3n,zd20,zd2n,f1d0,f1dn,f2d20,f2d2n,
     +     vaxplus)
      call compaxismsh(alf,zmin,rs,zs,ns,f1,f2,
     +     rd0,rdn,rd30,rd3n,zd20,zd2n,f1d0,f1dn,f2d20,f2d2n,
     +     vaxmin)
      vzax=(vaxplus-vaxmin)/(2*hloc)

c find polynomial u(theta)=c1*thet+c3*thet^3+c5*thet^5
      c1=udax*radtarget
      cc=umax-c1*thetmax
      dd=(udmax-c1)*thetmax
      aa=(5*cc-dd)/2
      bb=cc-aa
      c3=aa/thetmax**3
      c5=bb/thetmax**5
      uextrap= c1*theta+c3*theta**3+c5*theta**5

      c0=vax
      fdd=vddax*radtarget**2 + vzax*radtarget
      c2=fdd/2
      cc=vmax-c0-c2*thetmax**2;
      dd=(vdmax-2*c2*thetmax)*thetmax;
      aa=(6*cc-dd)/2
      bb=cc-aa
      c4=aa/thetmax**4
      c6=bb/thetmax**6
      vextrap= c0+ c2*theta**2+c4*theta**4+c6*theta**6

      return
      end

