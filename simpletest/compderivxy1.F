#include "flags.h"
      subroutine compderivxy1(n,x,y,sp,sinthet,costhet,thetd,thetdd,
     +    thetddd,kapmax)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c   
c     Computes all necessary derivatives of x,y, including             c
c     kappa and its derivatives                                        c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
#include "derivs1.h"
      integer j,n
      MODE x(0:*),y(0:*),sinthet(0:*),costhet(0:*),thetd(0:*),
     +       thetdd(0:*),thetddd(0:*),sp,kapmax

      kapmax=0

      xddd0   = -sp*thetd(0)**2
c      xp(0)   =  sp*costhet(0)
c      yp(0)   =  0
c      xdd(0)  =  0
c      ydd(0)  =  sp*costhet(0)*thetp(0)
c      xddd(0)  =  xddd0
c      xddd(0)  =  0
      kap(0)  =  2*thetd(0)/sp
      ydddd0  =  sp*(-thetd(0)**3+thetddd(0))
      kapmax  =  max(abs(kap(0)),kapmax)
      alf(0)=0
      alf(n)=pi
      delalf=pi/n

      do j=0,n
        alf(j)=j*delalf
        xd(j)   =  sp*costhet(j)
        yd(j)   =  sp*sinthet(j)
        xdd(j)  = -sp*sinthet(j)*thetd(j)  ! -yp*thetp
        ydd(j)  =  sp*costhet(j)*thetd(j)  !  xp*thetp
        xddd(j) = - ydd(j)*thetd(j)-yd(j)*thetdd(j)
        yddd(j) =   xdd(j)*thetd(j)+xd(j)*thetdd(j)
      enddo
      do j=1,n-1
        kap(j)  =  thetd(j)/sp + sinthet(j)/x(j)
        kapmax  =  max(abs(kap(j)),kapmax)
      enddo

c      xp(n)  = sp*costhet(n)
c      yp(n)  = 0
c      xdd(n) = 0
c      ydd(n) = sp*costhet(n)*thetp(n)
      kap(n) = 2*thetd(n)/sp
      xdddn  =  sp*thetd(n)**2
      yddddn = sp*( thetd(n)**3-thetddd(n))
      kapmax  =  max(abs(kap(n)),kapmax)

      call compderiv2(kap,kapd,kapdd,n)

      return
      end

