#include "flags.h"
      SUBROUTINE intbm(xx,yy,ns,xd0,y0,f1d0,f20,f2dd0,bintu,bintv,
     +  alfend,alfp,xp,yp,xdp,ydp,cc1,cc2,dd1,dd2,jmax)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computed "exact" integral of b for Mu,Mv, to test
c     needed to get the error E[b]
c     uses ns tchebb pts on [0,alfend] plus 2nd and 4th order corrections
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,ns,jmax
      MODE xx,yy,xd0,y0,f1d0,f20,f2dd0,bintu,bintv,
     +    alfend,alfp,xp,yp,xdp,ydp,cc1,cc2,dd1,dd2
      MODE alf(0:nmax),alfd(0:nmax),xsa(0:nmax),xia(0:nmax),
     +   f1a(0:nmax),f2a(0:nmax),delbet,beta,sumx,mm,eth,
     +   a2,a4,b,c2,c,c3,c4,i10,i11,i30,i31,i32,hold,m11,m12,m21,m22,
     +   r,alfd2l,dbetgl,d3betgl,dbetgr,d3betgr,xi0,
     +   m12p,dalfgul,m22p,dalfgvl,x2
      MODE ksq(0:nmax),F(0:nmax),E(0:nmax),gu(0:nmax),gv(0:nmax)

      delbet = pi/ns  !=delalf
      jmax=nint(4*xx/(pi/100))*ns/100
c      alfend=max(10*delbet,jmax*delbet)
      jmax=max(1*ns/100,jmax)
      alfend=jmax*delbet

c      read from testxapprox.m
      cc1=(5*xp/alfp-4*xd0-xdp)/(2*alfp**2)
      cc2=(xdp+2*xd0-3*xp/alfp)/(2*alfp**4)

      dd1=(2*(yp-y0)/alfp-ydp/2)/alfp
      dd2=(ydp/2-(yp-y0)/alfp)/alfp**3

      do k=1,ns+2
        beta=k*delbet
        alf(k)=(1-cos(beta/2))*alfend
        alfd(k)=sin(beta/2)*alfend/2
        alfd2l=alfend/4

c USE EXACT DERIVATIVE OF XSA,YSA-YY=XIA
        xsa(k)=xd0*alf(k) +   cc1*alf(k)**3  +   cc2*alf(k)**5
        xia(k)=(y0-yy) +  dd1*alf(k)**2 +  dd2*alf(k)**4    !through alfp
        f1a(k) =f1d0*alf(k)  
        f2a(k) =f20 +f2dd0*alf(k)**2/2  

        sumx = xsa(k) + xx
        ksq(k) = 4*xsa(k)*xx/( xia(k)**2 + sumx**2)
      enddo
 
      call ellip(ksq(1),F(1),E(1),ns+1)
 
      do k = 1,ns+2
          sumx = xsa(k) + xx

          mm =1-ksq(k)
          eth = E(k)/mm

          a2=xx**2 + xsa(k)**2 + xia(k)**2 
          a4=a2*a2
          b=2*xx*xsa(k)
          c2=sumx**2 + xia(k)**2
          c=sqrt(c2)
          c3=c*c2
          c4=c2*c2

          i10 = 4*F(k)/c
          i11 = 4/(b*c)*( a2*F(k) - c2*E(k) )
          i30 = 4/c3*eth
          i31 = 4/(b*c3)*(a2*eth-c2*F(k))
          i32 = 4/(b**2*c3)*( a4*eth-2*a2*c2*F(k)+c4*E(k) )

          hold = i11+ (xsa(k)**2+xx**2)*i31- xsa(k)*xx*(i30+i32)
          m11 =xsa(k)*hold 
          m12 =xsa(k)*xia(k)*( xsa(k)*i31-xx*i30 )
          gu(k) = m11*f1a(k)+m12*f2a(k)

          m21 =xsa(k)*xia(k)*( xsa(k)*i30-xx*i31 )
          hold = i10+xia(k)**2*i30
          m22 =xsa(k)*hold
          gv(k) = m21*f1a(k)+m22*f2a(k)
c
         gu(k)=gu(k)*alfd(k)
         gv(k)=gv(k)*alfd(k)
c         print*,k,xsa(k),xia(k),alfd(k),f1a(k),f2a(k)
      enddo
 
      bintu=0
      bintv=0
      do k = 1,ns-1
         bintu = bintu+ gu(k)
         bintv = bintv+ gv(k)
      enddo
      bintu=bintu*delbet
      bintv=bintv*delbet

      x2=xx**2
      xi0 = y0-yy
      r = sqrt(xx**2 + xi0**2)

      dbetgl=0   !since alfd=0 at left
      m12p  =-pi2*xd0*xx*xi0/r**3
      dalfgul = m12p*f20 
      d3betgl=3*dalfgul*alfd2l**2 
      dbetgr=(-gu(ns+2)+8*gu(ns+1)-8*gu(ns-1)+gu(ns-2))/(12*delbet)
      d3betgr=(gu(ns+2)-2*gu(ns+1)+2*gu(ns-1)-gu(ns-2))/(2*delbet**3)
      bintu= bintu+gu(ns)*delbet/2
     +    - delbet**2*(dbetgr -dbetgl )/12
     +    + delbet**4*(d3betgr-d3betgl)/720
c      print*,']; %d3bu',d3betgl
c      print*,dalfgul

      dbetgl=0   !since alfd=0 at left
      m22p  = pi2*xd0*(x2+2*xi0**2)/r**3
      dalfgvl = m22p*f20 
      d3betgl=3*dalfgvl*alfd2l**2 
      dbetgr=(-gv(ns+2)+8*gv(ns+1)-8*gv(ns-1)+gv(ns-2))/(12*delbet)
      d3betgr=(gv(ns+2)-2*gv(ns+1)+2*gv(ns-1)-gv(ns-2))/(2*delbet**3)
      bintv= bintv+gv(ns)*delbet/2
     +    - delbet**2*(dbetgr -dbetgl )/12
     +    + delbet**4*(d3betgr-d3betgl)/720
c      print*,'%d3bv',d3betgl
c      print*,dalfgvl
 
      bintu = -bintu/pi8
      bintv = -bintv/pi8

      return
      end


      SUBROUTINE trapbm(xx,yy,alf,jmax,xd0,y0,f1d0,f20,f2dd0,
     +  cc1,cc2,dd1,dd2,trapbu,trapbv)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute trapezoid approximation T_6 of B_q for gu,gv, to test
c     needed to get the error E[b]
c     uses ns tchebb pts on [0,alfend] plus 2nd and 4th order corrections
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,ns,jmax
c INPUT/OUTPUT
      MODE xx,yy,xd0,y0,f1d0,f20,f2dd0,trapbu,trapbv,
     +    cc1,cc2,dd1,dd2,alf(0:*)
      MODE delalf,alfd(0:nmax),xsa(0:nmax),xia(0:nmax),
     +   f1a(0:nmax),f2a(0:nmax),sumx,mm,eth,
     +   a2,a4,b,c2,c,c3,c4,i10,i11,i30,i31,i32,hold,m11,m12,m21,m22,
     +   r,alfd2l,dbetgl,d3betgl,dbetgr,d3betgr,
     +   m12p,m22p,m11pp,m12ppp,gpul,
     +   gpppul,dgur,d3gur,m21pp,m22ppp,gpvl,gpppvl,dgvr,d3gvr
      MODE ksq(0:nmax),F(0:nmax),E(0:nmax),gu(0:nmax),gv(0:nmax)
      MODE xddd0,ydd0,alfa,dfine,sing,rsq,x2,x4,xi0,xi2,xi4

      ns=jmax
      delalf = alf(jmax)/jmax  !=delalf
      dfine=delalf/100
      xddd0=6*cc1
      ydd0=2*dd1

c     compute b(xx,yy) where yy.ne.0 (point not on axis)

      do k=1,ns
        xsa(k)=xd0*alf(k) +   cc1*alf(k)**3  +   cc2*alf(k)**5
        xia(k)=(y0-yy) +  dd1*alf(k)**2 +  dd2*alf(k)**4  
        f1a(k) =f1d0*alf(k)  
        f2a(k) =f20 +f2dd0*alf(k)**2/2  

        sumx = xsa(k) + xx
        ksq(k) = 4*xsa(k)*xx/( xia(k)**2 + sumx**2)
      enddo
      do k=ns+1,ns+5
        alfa=alf(ns)+(k-ns-3)*dfine
        xsa(k)=xd0*alfa +   cc1*alfa**3  +   cc2*alfa**5
        xia(k)=(y0-yy) +  dd1*alfa**2 +  dd2*alfa**4  
        f1a(k) =f1d0*alfa  
        f2a(k) =f20 +f2dd0*alfa**2/2  

        sumx = xsa(k) + xx
        ksq(k) = 4*xsa(k)*xx/( xia(k)**2 + sumx**2)
      enddo
 
      call ellip(ksq(1),F(1),E(1),ns+4)
 
      do k = 1,ns+5
          sumx = xsa(k) + xx

          mm =1-ksq(k)
          eth = E(k)/mm

          a2=xx**2 + xsa(k)**2 + xia(k)**2 
          a4=a2*a2
          b=2*xx*xsa(k)
          c2=sumx**2 + xia(k)**2
          c=sqrt(c2)
          c3=c*c2
          c4=c2*c2

          i10 = 4*F(k)/c
          i11 = 4/(b*c)*( a2*F(k) - c2*E(k) )
          i30 = 4/c3*eth
          i31 = 4/(b*c3)*(a2*eth-c2*F(k))
          i32 = 4/(b**2*c3)*( a4*eth-2*a2*c2*F(k)+c4*E(k) )

          hold = i11+ (xsa(k)**2+xx**2)*i31- xsa(k)*xx*(i30+i32)
          m11 =xsa(k)*hold 
          m12 =xsa(k)*xia(k)*( xsa(k)*i31-xx*i30 )
          gu(k) = m11*f1a(k)+m12*f2a(k)

          m21 =xsa(k)*xia(k)*( xsa(k)*i30-xx*i31 )
          hold = i10+xia(k)**2*i30
          m22 =xsa(k)*hold
          gv(k) = m21*f1a(k)+m22*f2a(k)
      enddo
 
      trapbu=0
      trapbv=0
      do k = 1,ns-1
         trapbu = trapbu+ gu(k)
         trapbv = trapbv+ gv(k)
      enddo
      trapbu=trapbu*delalf
      trapbv=trapbv*delalf

      x2 = xx**2
      x4 = x2**2

      xi0 = y0-yy
      xi2 = xi0**2
      xi4 = xi0**4
      rsq = x2+xi2
      r = sqrt(rsq)
      m11pp = pi2*xd0**2*xx*(x2-2*xi2)/r**5
      m12p  =-pi2*xd0*xx*xi0/r**3
      m12ppp=-pi*xx*( 9*xd0**3*(x2-4*xi2)*xi0 
     +  + 2*xddd0*xi0*rsq**2
     +  + 6*xd0*(x4-x2*xi2-2*xi4)*ydd0 )/r**7
      gpul = m12p*f20 
      gpppul=3*m11pp*f1d0 + m12ppp*f20 + 3*m12p*f2dd0

      dgur=(-gu(ns+5)+8*gu(ns+4)-8*gu(ns+2)+gu(ns+1))/(12*dfine)
      d3gur=(gu(ns+5)-2*gu(ns+4)+2*gu(ns+2)-gu(ns+1))/(2*dfine**3)
      trapbu= trapbu+gu(ns)*delalf/2
     +    - delalf**2*(dgur -gpul )/12
     +    + delalf**4*(d3gur-gpppul)/720

      m21pp =-pi2*xd0**2*xi0*(x2-2*xi2)/r**5
      m22p  = pi2*xd0*(x2+2*xi2)/r**3
      m22ppp= pi*(
     +    2*xddd0*rsq**2*(x2+2*xi2)
     +  + 3*xd0**3*(x4+8*x2*xi2-8*xi4)
     +  + 6*xd0*xi0*(x4-x2*xi2-2*xi4)*ydd0 )/r**7
      gpvl = m22p*f20 
      gpppvl=3*m21pp*f1d0 + m22ppp*f20 + 3*m22p*f2dd0

      dgvr=(-gv(ns+5)+8*gv(ns+4)-8*gv(ns+2)+gv(ns+1))/(12*dfine)
      d3gvr=(gv(ns+5)-2*gv(ns+4)+2*gv(ns+2)-gv(ns+1))/(2*dfine**3)
      trapbv= trapbv+gv(ns)*delalf/2
     +    - delalf**2*(dgvr -gpvl )/12
     +    + delalf**4*(d3gvr-gpppvl)/720
 
      trapbu = -trapbu/pi8
      trapbv = -trapbv/pi8

      return
      end

