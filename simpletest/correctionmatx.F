#include "flags.h"

      SUBROUTINE correctionmatx(alf,ns,xx,yy,xproj,yproj,delta,
     +        alfp,corrmu,corrmv,
     +        xdp,xddp,xdddp,ydp,yddp,ydddp,
     +        f1p,f1dp,f1ddp,f2p,f2dp,f2ddp)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Input: alfp and all derivatives at alfp 
c            the point xx,yy off sheet at which velo is to be evaluated
c            the normal projection xproj,yproj of that point onto sheet 
c            the distance dist=delta from xx,yy to xp,yp
c
c     Output: corruq[11,12,22]  vector containing 
c             corrvq[11,12,22]  the corrections for all Qs,
c             given by
c             (int[Qapprox]-T[Qapprox]) NOT over LOCAL interval centered
c             at alf(k0) which is near alfp): B=[alf(k0-10),alf(k0+10)]
c             Instead use GLOBAL interval [alf(0),alf(n)]
c
c            where Qapprox=H in notes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer ns
C     INPUT/OUTPUT
      MODE alf(0:*),xx,yy,xproj,yproj,delta,
     +        alfp,corrmu,corrmv,
     +        xdp,xddp,xdddp,ydp,yddp,ydddp,
     +        f1p,f1dp,f1ddp,f2p,f2dp,f2ddp
C     LOCAL
      MODE cmsq,cm,fact,coeffm1,coeffm2,coeffm3,coeffm4,corrs(22),
     + cfin0,cfin1,cfin2,cfin3,dfin0,dfin1,dfin2
      MODE corr(4)

      call compcoeffm(xx,yy,xproj,yproj,xdp,xddp,xdddp,
     +      ydp,yddp,ydddp)

      cmsq=(xdp**2+ydp**2)+(xproj-xx)*xddp+(yproj-yy)*yddp
      cm=sqrt(cmsq)
c      c=2*( (xproj-xx)*xdp + (yproj-yy)*ydp )

      fact = delta/cm
      coeffm1=1/(delta**2)
      coeffm2=1/(delta**4)
      coeffm3=1/(delta**6)
      coeffm4=1/(delta**8)

c      c = c/(delta*cm)
c      print*,'HERE, c=',c,delta,cm
c      print*,delta,c

      call compallcorr(alf,alfp,fact,ns,corrs)

c   dfinj : coefficient of H0j  logterms
c   cfinj : coefficient of H1j  u^j/(1+u^2) scaled by fact^j

      call compm11(xx,fact,cfin0,cfin1,cfin2,cfin3,dfin0,dfin1,dfin2,
     +     f1p,f1dp,f1ddp)
      corr(1)= (dfin0*corrs(1)+dfin1*corrs(2)+dfin2*corrs(3)) 
     +     +coeffm1*( cfin0*corrs(4)+cfin1*corrs(5)
     +              + cfin2*corrs(6)+cfin3*corrs(7) )

      call compm12(xx,fact,cfin0,cfin1,cfin2,cfin3,dfin0,dfin1,dfin2,
     +     f2p,f2dp,f2ddp)
      corr(2)= (dfin0*corrs(1)+dfin1*corrs(2)+dfin2*corrs(3)) 
     +     +coeffm1*( cfin0*corrs(4)+cfin1*corrs(5)
     +              + cfin2*corrs(6)+cfin3*corrs(7) )

      call compm21(xx,fact,cfin0,cfin1,cfin2,cfin3,dfin0,dfin1,dfin2,
     +     f1p,f1dp,f1ddp)
      corr(3)= (dfin0*corrs(1)+dfin1*corrs(2)+dfin2*corrs(3)) 
     +     +coeffm1*( cfin0*corrs(4)+cfin1*corrs(5)
     +              + cfin2*corrs(6)+cfin3*corrs(7) )

      call compm22(xx,fact,cfin0,cfin1,cfin2,cfin3,dfin0,dfin1,dfin2,
     +     f2p,f2dp,f2ddp)
      corr(4)= (dfin0*corrs(1)+dfin1*corrs(2)+dfin2*corrs(3)) 
     +     +coeffm1*( cfin0*corrs(4)+cfin1*corrs(5)
     +              + cfin2*corrs(6)+cfin3*corrs(7) )

      corrmu = -(corr(1)+corr(2))/pi8
      corrmv = -(corr(3)+corr(4))/pi8

      return
      end
