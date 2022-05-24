#include "flags.h"
      subroutine compxy(n,L,sp,thet,thetp,thetpp,thetppp,y0,
     +                  sinthet,costhet,x,y,vol,vol0) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Solves xp=sp*cos(thet), where sp=const
c     Solves yp=sp*sin(thet), where sp=const
c     BUT: uses sp such that volume (computed to 6th order)=1
c          returns new sp
c     input: n,L,sp,thet,y0,
c     output: thetp,thetpp,thetppp,sinthet,costhet,x,y,vol,vol0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c#include "constants.h"
c#include "vars.h"
      integer j,n
      MODE x(0:*),y(0:*),L,thet(0:*),thetp(0:*),thetppp(0:*),y0,sp,
     +       ymid,sinthet(0:*),costhet(0:*),vol0,delalf
      MODE hold(0:nmax),thetpp(0:nmax),hsq,a1,a2,int1(0:nmax),
     +       int2(0:nmax),
     +       fpl,fpr,gpl,gpr,fpppl,fpppr,gpppl,gpppr,dfp,dgp,dfppp,dgppp
      MODE spt,vol

      delalf=pi/n

C     compute derivatives of thet
      do j=0,n
c        hold(j)=thet(j)-alf(j)
        hold(j)=thet(j)-j*delalf
      enddo
      call deriv3odd(hold,thetp,thetpp,thetppp,n)
      do j=0,n
         thetp(j)   = thetp(j)+1
         sinthet(j) = sin(thet(j))
         costhet(j) = cos(thet(j))
      enddo

C     integrate f=cos(thet), g=sin(thet) from 0 to alfj with trapezoid
      int1(0) = 0
      int2(0) = 0
      do j=1,n
        int1(j) = int1(j-1) + (costhet(j-1)+costhet(j))*delalf/2
        int2(j) = int2(j-1) + (sinthet(j-1)+sinthet(j))*delalf/2
      enddo

C     add corrections add endpoints using MacLaurin series
      hsq=delalf**2
      fpl = 0           !fp(0)
      gpl = thetp(0)    !gp(0)
      fpppl = 0           !fppp(0)
      gpppl = thetppp(0)-thetp(0)**3    !gpp(0)
      do j=1,n
        fpr = -sinthet(j)*thetp(j)
        gpr =  costhet(j)*thetp(j)
        a1 = thetppp(j)-thetp(j)**3
        a2 = 3*thetp(j)*thetpp(j)
        fpppr = - sinthet(j)*a1 - costhet(j)*a2
        gpppr = + costhet(j)*a1 - sinthet(j)*a2
        dfp = fpr-fpl
        dgp = gpr-gpl
        dfppp = fpppr-fpppl
        dgppp = gpppr-gpppl
c        int1(j) = int1(j) - hsq*( dfp/12)  !fourth order
c        int2(j) = int2(j) - hsq*( dgp/12)  !fourth order
        int1(j) = int1(j) - hsq*( dfp/12 - dfppp/720*hsq )
        int2(j) = int2(j) - hsq*( dgp/12 - dgppp/720*hsq )
      enddo
      do j=0,n
        x(j) = sp*int1(j)  !?? - alf(j)*salf/pi*int1(n)
        y(j) = sp*int2(j) + y0
      enddo

      call compvol(n,x,sp,sinthet,costhet,thetp,vol)

      spt=(vol0/vol)**(1/3.q0)
      sp=spt*sp
c      print*,vol0,vol,vol0/vol,spt
c      print*,spt-1,sp-0.95d0
      do j=0,n
        x(j) = spt*x(j)
        y(j) = spt*y(j) 
      enddo
      x(n)=0

      call compvol(n,x,sp,sinthet,costhet,thetp,vol)

      return
      end


      subroutine compvol(n,x,sp,sinthet,costhet,thetp,vol)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c compute volume using
c formula pi*[sum x_k^2 yp_k delalf 
c        + (delalf^4/720)*6*(xp^2(n)ypp(n)-xp^2(0)ypp(0)]
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,n
      MODE x(0:*),sinthet(0:*),costhet(0:*),thetp(0:*),sp,delalf
      MODE sum,vol,xp(0:nmax),yp(0:nmax),xp0,xpn,ypp0,yppn

      delalf=pi/n
      do j=0,n
        yp(j) = sp*sinthet(j)
      enddo
      xp0  = sp*costhet(0)
      xpn  = sp*costhet(n)
      ypp0 = xp0*thetp(0)
      yppn = xpn*thetp(n)

      sum=0
      do j=1,n-1
        sum=sum+x(j)**2*yp(j)
      enddo
      sum=sum*delalf+delalf**4*(xpn**2*yppn-xp0**2*ypp0)/120
      vol=pi*sum

      return
      end
