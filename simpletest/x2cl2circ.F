#include "flags.h"
      SUBROUTINE x2cl2circ(x,y,ns,yes,k,alfp,xproj,yproj,dist)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute exact values of d, xproj,yproj,alfp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer ns,k
      MODE x,y,xproj,yproj,dist,dcut,dcut2,alfp,dels,rad,delalf,y0
      logical yes

      rad=0.95q0  !for xs1, 0.6 for xs2
      dist=sqrt(x**2+(y+4)**2)  ! distance from (0,-4)
      dist=abs(dist-rad)          ! distance from circle
      delalf=pi/ns
      dels=rad*delalf

      y0=-4-.95q0
      dcut=5*dels
      dcut2=0.1*dels
c      if ((dist.le.dcut).and.(x.ge.0.6*abs(y-y0))) then 
      if (dist.le.dcut) then 
         yes =.true.
         xproj=x*rad/sqrt(x**2+(y+4)**2)
         yproj=(y+4)*rad/sqrt(x**2+(y+4)**2)-4
         alfp=asin(xproj/rad)
c         dist = sqrt( (x-xproj)**2 + (y-yproj)**2 )
         k=nint(alfp/delalf)
      endif

      return
      end


      SUBROUTINE x2close2S(j,x,y,xs,ys,ns,yes,kcl,alfp,xproj,yproj,
     +   dist,dcut,dels)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    xxxxxyyyyy  FILL IN
c     3) if dist<dcut, set yesint=true, return xproj,yproj,dist,dcut
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer k,kcl,ns,j
      MODE x,y,xproj,yproj,dist,dcut,alfp,xs(0:*),ys(0:*),delalf
      MODE ax,ay,bx,by,cx,cy,dels,la,lb,lc,adotb,adotc,fact,dalf
      MODE alf(0:nmax),c(6),rad,half
      logical yes

      half=1/2.q0
      rad=0.95q0  !for xs1, 0.6 for xs2

      yes=.false.
      delalf=pi/ns
      do k=0,ns
        alf(k)=k*delalf
      enddo

      dist=10
      do k=2,ns-2
        ax = xs(k+1)-xs(k-1)
        ay = ys(k+1)-ys(k-1)
        bx = x-xs(k+1)
        by = y-ys(k+1)
        cx = x-xs(k-1)
        cy = y-ys(k-1)
        la=sqrt(ax**2+ay**2) 
        dels=la/2
        lb=sqrt(bx**2+by**2)
        lc=sqrt(cx**2+cy**2)

c        if ( (lb+lc).le.5*dels ) then
        if ( (lb+lc).le.8*dels ) then
          adotb=-(ax*bx+ay*by)
          adotc=ax*cx+ay*cy
          if ((adotc.ge.0).and.(adotb.ge.0)) then
            dalf=adotc/la**2*2*delalf
            alfp=(k-1)*delalf+dalf
            fact=adotc/la**2
c            print*,fact
            if (fact.le.half) then
               kcl = k-1
               xproj=xs(k-1)+fact*2*(xs(k)-xs(k-1))
               yproj=ys(k-1)+fact*2*(ys(k)-ys(k-1))
c               print*,'x,y lin',xproj,yproj
c NEW 4/23/2019: use better approximation of xproj,yproj
c              EXACT VALUES, EXACT alfp
               xproj=x*rad/sqrt(x**2+(y+4)**2)
               yproj=(y+4)*rad/sqrt(x**2+(y+4)**2)-4
               alfp=asin(xproj/rad)
c               call fithermite(alf(kcl-2),alfp,xs(kcl-2),c)
c               xproj=c(1)
c               call fithermite(alf(kcl-2),alfp,ys(kcl-2),c)
c               yproj=c(1)
            else
              kcl=k
               xproj=xs(k)+2*(fact-half)*(xs(k+1)-xs(k))
               yproj=ys(k)+2*(fact-half)*(ys(k+1)-ys(k))
c               print*,'xsrad',sqrt(xs(k)**2+(ys(k)+4)**2) =0.6
c               print*,'x,y linr',xproj,yproj
c NEW 4/23/2019: use better approximation of xproj,yproj
c              EXACT VALUES, EXACT alfp
               xproj=x*rad/sqrt(x**2+(y+4)**2)
               yproj=(y+4)*rad/sqrt(x**2+(y+4)**2)-4
c               print*,'x,y exct',xproj,yproj
c               print*,alfp,sin(alfp)*rad,xproj
c               print*,asin(xproj/rad)
               alfp=asin(xproj/rad)
c               call fithermite(alf(kcl-2),alfp,xs(kcl-2),c)
c               xproj=c(1)
c               call fithermite(alf(kcl-2),alfp,ys(kcl-2),c)
c               yproj=c(1)
            endif
            dist = sqrt( (x-xproj)**2 + (y-yproj)**2 )
            dcut=4*dels
            if (dist.le.dcut) yes =.true.
            goto 20
          endif
        endif
      enddo

c Note: kcl is always the value on the left of interval
20    continue
      return
      end
