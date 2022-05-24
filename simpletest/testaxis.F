#include "flags.h"
      PROGRAM testaxis
#include "params.h"
      integer ns,nr,j
      MODE rad,delalf,alf(0:nmax),xs(0:nmax),ys(0:nmax),
     +     xd(0:nmax),yd(0:nmax),u(0:nmax),v(0:nmax),
     +     xd30,xd3n,yd20,yd2n,ud0,udn,vd20,vd2n,r1,r2,delr,
     +  vax,udax,vddax,y(0:nmax),r
      MODE f1(0:nmax),f2(0:nmax),f1d0,f1dn,f2d20,f2d2n
      character*10 flend

      read(5,*)
      read(5,'(a25)')flend
      read(5,*)ns       !number of points on bubble
c      read(5,*)rad      !initial radius for outer bubble 1 rad1
c      read(5,*)y01       !initial center of outer bubble 1

      rad=0.95
      delalf = pi/ns
      do j=0,ns
        alf(j)=j*delalf
        xs(j)=  rad*sin(alf(j))
        ys(j)= -rad*cos(alf(j))
  
        xd(j)=  rad*cos(alf(j))
        yd(j)=  rad*sin(alf(j))

        u(j) = sin(alf(j))
        v(j) = cos(alf(j))

        f1(j) = sin(alf(j))
        f2(j) = cos(alf(j))
      enddo
      xd30 = -rad
      xd3n =  rad
      yd20 =  rad
      yd2n = -rad
      ud0  =   1
      udn  =  -1
      vd20 =  -1
      vd2n =   1
      f1d0  =   1
      f1dn  =  -1
      f2d20 =  -1
      f2d2n =   1

      r1=0.8d0
      r2=1.4d0
      nr=600
c equally spaced points about r=1
      delr = 2*(r2-r1)/(2*nr+3)
      print*,'a=['
      do j=0,nr+1
        r=r1+j*delr
        y(j)= -r
c        y=-.8075
        vax=0
        udax=0
        vddax=0
c        call compaxism(alf,y(j),xs,ys,ns,f1,f2,xd(0),xd(ns),xd30,xd3n,
c     +           yd20,yd2n,f1d0,f1dn,f2d20,f2d2n,vax,udax,vddax)
c        call compaxisq(alf,y(j),xs,ys,ns,u,v,xd,yd,
c     +     xd30,xd3n,yd20,yd2n,ud0,udn,vd20,vd2n,vax,udax,vddax)
c        print*, j,y(j),vax ,udax,vddax
        call compaxismsh(alf,y(j),xs,ys,ns,f1,f2,xd(0),xd(ns),xd30,xd3n,
     +           yd20,yd2n,f1d0,f1dn,f2d20,f2d2n,vax)
c        call compaxisqsh(alf,y(j),xs,ys,ns,u,v,xd,yd,
c     +     xd30,xd3n,yd20,yd2n,ud0,udn,vd20,vd2n,vax)
        print*, j,y(j),vax 
      enddo
      print*,'];'

      stop
      end
