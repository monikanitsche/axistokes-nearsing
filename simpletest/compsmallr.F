#include "flags.h"
      SUBROUTINE compsmallr(rr,zz,rmax,alf,rs,zs,ns,u,v,rd,zd,
     +           rd30,zd20,ud0,vd20,
     +           rd3n,zd2n,udn,vd2n,uextrap,vextrap)
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
      MODE rr,zz,rd30,zd20,ud0,vd20,rd3n,zd2n,udn,vd2n
      MODE alf(0:*),rs(0:*),zs(0:*),u(0:*),v(0:*),rd(0:*),zd(0:*)
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
c      zcirc=-4-rad*cos(theta)
c      d=zcirc-zz
c      z0=zs(0)-d
c      theta = asin(rr/radtarget)

      thetmax = asin(rmax/radtarget)
      rmax=radtarget*sin(thetmax)
      zmax=-4-radtarget*cos(thetmax)
c      thetmax=asin(rmax/rad)
c      zmax=-4-rad*cos(thetmax)-d
c      radmax=sqrt((zmax+4)**2+rmax**2)
c      thetmax=asin(rmax/radmax)

      call velodSonP(rmax,zmax,alf,rs,zs,ns,u,v,rd,zd,
     +           rd30,zd20,ud0,vd20,
     +           rd3n,zd2n,udn,vd2n,umax,vmax)
c      print*,rmax,zmax,umax


c      rr=0.012245549943937414d0  
c      zz=-4.9744186076484933d0
c      call velodSonP(rr,zz,alf,rs,zs,ns,u,v,rd,zd,
c     +           rd30,zd20,ud0,vd20,
c     +           rd3n,zd2n,udn,vd2n,umax,vmax)
c      print*,rr,zz
c      print*,thetmax,umax,vmax
c      print*,umax,vmax,'in compsmall fin'

      thetmaxplus=thetmax+delalf/10
      rmaxplus=radtarget*sin(thetmaxplus)
      zmaxplus=-4-radtarget*cos(thetmaxplus)
      call velodSonP(rmaxplus,zmaxplus,alf,rs,zs,ns,u,v,rd,zd,
     +           rd30,zd20,ud0,vd20,
     +           rd3n,zd2n,udn,vd2n,umaxplus,vmaxplus)
      thetmaxmin=thetmax-delalf/10
      rmaxmin=radtarget*sin(thetmaxmin)
      zmaxmin=-4-radtarget*cos(thetmaxmin)
      call velodSonP(rmaxmin,zmaxmin,alf,rs,zs,ns,u,v,rd,zd,
     +           rd30,zd20,ud0,vd20,
     +           rd3n,zd2n,udn,vd2n,umaxmin,vmaxmin)
      udmax=(umaxplus-umaxmin)/(thetmaxplus-thetmaxmin)
      vdmax=(vmaxplus-vmaxmin)/(thetmaxplus-thetmaxmin)
c      print*,thetmax,umax,udmax , ' thetmax,umax,udmax'
c      print*,thetmax,vmax,vdmax, ' thetmax,vmax,vdmax'

c find values of derivatives on axis: vax,vddax,udax, and vzax
c here ud,vdd denote derrivatives wrt r on axis
c vz denote derrivatives wrt z on axis

      call compaxisq(alf,z0,rs,zs,ns,u,v,rd,zd,
     +     rd30,rd3n,zd20,zd2n,ud0,udn,vd20,vd2n,vax,udax,vddax)

      d=abs(z0-zs(0))
      hloc=min(delalf/100,d/10)
      zplus=z0+hloc
      zmin=z0-hloc
c      call compaxisq(alf,zplus,rs,zs,ns,u,v,rd,zd,
c     +     rd30,rd3n,zd20,zd2n,ud0,udn,vd20,vd2n,vaxplus,udaxplus,vddax)
c      call compaxisq(alf,zmin,rs,zs,ns,u,v,rd,zd,
c     +     rd30,rd3n,zd20,zd2n,ud0,udn,vd20,vd2n,vaxmin,udaxmin,vddax)
      call compaxisqsh(alf,zplus,rs,zs,ns,u,v,rd,zd,
     +     rd30,rd3n,zd20,zd2n,ud0,udn,vd20,vd2n,vaxplus)
      call compaxisqsh(alf,zmin,rs,zs,ns,u,v,rd,zd,
     +     rd30,rd3n,zd20,zd2n,ud0,udn,vd20,vd2n,vaxmin)
      vzax=(vaxplus-vaxmin)/(2*hloc)
c      udzax=(udaxplus-udaxmin)/(2*hloc)

c set uetrap,vextrap
c      uextrap = udax*delr + udzax*delr*delz
c      c3=(umax-udax*alfmax)/alfmax**3
c      uextrap = udax*alf + c3*alf**3
c      vextrap0 = vax+vzax*delz + vddax*delr**2/2

c find polynomial u(theta)=c1*thet+c3*thet^3+c5*thet^5
      c1=udax*radtarget
      cc=umax-c1*thetmax
      dd=(udmax-c1)*thetmax
      aa=(5*cc-dd)/2
      bb=cc-aa
      c3=aa/thetmax**3
      c5=bb/thetmax**5
      uextrap= c1*theta+c3*theta**3+c5*theta**5
c      print*,c1,'udax'

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
c      print*,c0,c2,'c0,c2'

      return
      end

