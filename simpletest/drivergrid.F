#include "flags.h"
      program drivergrid1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     Computes and prints us,vs,udb,vd for test case1: 
c           circle with prescribed density.
c     computed on an annular grid
c
c     copied from drivertest1.F 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
#include "vars.h"      
      integer j,nr,k,ns(10),i
      MODE alfa,delalf,pih
      MODE rad,thet1d(0:nmax),thet1ddd(0:nmax),xs1(0:nmax),
     +  ys1(0:nmax),sinthet1(0:nmax),costhet1(0:nmax),
     +  thet1dd(0:nmax),kapmax1,vol01,vol1,
     +  xs2(0:nmax),ys2(0:nmax),rint(10),r,r1,r2,
     +  delthet,space,delr,thet
      MODE u(0:nmax),v(0:nmax)
      MODE udb(0:nmax),vdb(0:nmax),udbs(0:nmax),vdbs(0:nmax)
      MODE f11(0:nmax),f21(0:nmax)
      MODE uin(0:nmax),vin(0:nmax)
      character*25 flend
      character*4 lab(10)

      read(5,*)
      read(5,'(a25)')flend
      read(5,*)ns1       !number of points on bubble 
      read(5,*)rad      !initial radius for outer bubble 1 rad1
      read(5,*)y01       !initial center of outer bubble 1

c      open(9,file='datgriddb/info'//flend)
c      write(9,*)ns1, '   =n'
c#ifdef FOURTH
c      write(9,*)'FOURTH'
c#endif
c#ifdef THIRD
c      write(9,*)'THIRD'
c#endif
c#ifdef DBLE
c      write(9,*)'DBLE'
c#endif
c#ifdef QUAD
c      write(9,*)'QUAD'
c#endif
c#ifdef ROUNDOFFFIX
c      write(9,*)'ROUNDOFFFIX'
c#endif
c#ifdef NOROUNDOFFFIX
c      write(9,*)'NOROUNDOFFFIX'
c#endif
c      close(9)

c     SET xs1,ys1 exact to avoid errors, for now
      rad=0.95q0
      delalf=pi/ns1
      do j=0,ns1
        alfa=j*delalf
        thet1(j)=alfa
      enddo
      L1=pi*rad
      sp1=L1/pi
      vol01=4*pi*rad**3/3
      call compxy(ns1,L1,sp1,thet1,thet1d,thet1dd,thet1ddd,y01,
     +            sinthet1,costhet1,xs1,ys1,vol1,vol01)
      call compderivxy1(ns1,xs1,ys1,sp1,sinthet1,costhet1,thet1d,
     + thet1dd,thet1ddd,kapmax1)

      do j=0,ns1
        xs1(j)=     rad*sin(thet1(j))
        ys1(j)= -4 -rad*cos(thet1(j))
      enddo
c     Initialize density for both Ms and Qs
      do j=0,ns1
         f11(j)=sin(thet1(j))
         f21(j)=cos(thet1(j))
         uin(j)=sin(thet1(j))
         vin(j)=cos(thet1(j))
      enddo

c     SET grid xs2,ys2 : normal line or vertical line or grid
c     use even ng in theta direction this way symmetric not including top
c     j=1,ng/2 covers first quadrant, save only this, then copy to RHS
c     n=100, ng=1000: 10 points/ grid, 500 pts saved
c     to keep spacing in thet and r direction same use
c     mg= (rout-rin)/dels, dels=pi*rad/(ng+1), use even mg, only k=1,mg

      call init(rad,rint,lab,ns)

      pih=pi/2
      open(8,file='datgridsg/datsg'//flend)
      open(9,file='datgriddb/datdb'//flend)
c      open(10,file='dataxisq'//flend)
c      write(10,*)'a=['

      do i=1,6
c      do i=3,5
c      do i=4,4
        write(8,*)lab(i)
        write(9,*)lab(i)
        r1=rint(i)
        r2=rint(i+1)
        ns2=ns(i)
        delthet=pih/ns2
        space=rad*delthet
        nr=2*nint((r2-r1)/space/2)
        delr=(r2-r1)/(nr+1)
        do j=0,nr+1
          r=r1+j*delr
          do k=0,ns2
c          do k=0,10
            thet=k*delthet
            xs2(k)=r*sin(thet)
            ys2(k)=-(r*cos(thet)+4)
c            print*,k,xs2(k),ys2(k)
          enddo
c          stop
          xs2(ns2+1)=0
          ys2(ns2+1)=-4
          call IntMs(xs1,ys1,ns1,f11,f21,xs2,ys2,ns2+1,u,v)
          do k=0,ns2+1
             udbs(k)=0
             vdbs(k)=0
          enddo
         call IntQs(xs1,ys1,ns1,uin,vin,xs2,ys2,ns2+1,udb,vdb)
        enddo
        write(8,*)'];'
        write(9,*)'];'
      enddo
c        write(10,*)'];'

      stop
      end


      SUBROUTINE init(rad,rint,lab,ns)
      implicit none
      integer ns(*)
      MODE rint(*),rad
      character*4 lab(*) 

c      rint(1)=0.1*rad  !ns2=50 or 20
      rint(1)=0.01*rad  !ns2=50 or 20
      rint(2)=0.3*rad  !ns2=150 or 60
      rint(3)=0.6*rad  !ns2=500 or 200
      rint(4)=0.85*rad !ns2=1000 or 400
c      rint(4)=0.95*rad !ns2=1000 or 400
      rint(5)=1.2*rad  !ns2=500 or 200
c      rint(5)=1.05*rad  !ns2=500 or 200
      rint(6)=1.4*rad  !ns2=150 or 200
      rint(7)=1.7*rad  
      lab(1)='z1=['
      lab(2)='z2=['
      lab(3)='z3=['
      lab(4)='z4=['
      lab(5)='z5=['
      lab(6)='z6=['
c      ns(1)=  50 ! 20
      ns(1)=  70 ! 20
      ns(2)= 150 ! 60
      ns(3)= 500 !200
      ns(4)=1000 !400
c      ns(4)=200 !400
      ns(5)= 500 !200
      ns(6)= 150 ! 60

      return
      end
