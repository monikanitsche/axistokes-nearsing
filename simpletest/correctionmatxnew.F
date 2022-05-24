c Contains
c SUBROUTINES 
c     correctionmatxnew: given all coefficients at basepoint,
c                        outputs corrections for intQu intQv
c     corrm11: all these output corrections and local term 
c     corrm12: for the corresponding components
c     corrm21:
c     corrm22:
c FUNCTIONS corrqrho4, corrqrho2, corrqlog in other file
c SUBROUTINES compcoeffnew also in other file

#include "flags.h"
      SUBROUTINE correctionmatxnew(alf,alfb,ns,d,k0,
     +      xb,x0,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +      f10,f11,f12,f13,f20,f21,f22,f23,
     +      corru,corrv)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    alfb : basepoint
c    xb,xd,xdd,xddd,.. values AT BASEPOINT alfb
c    x0,y0: target point (y0 not needed, in dy=yb-y0)
c
c Use expansions in bubbles/math/math2022:
c        expm11redov2.pdf
c        expm12redov2.pdf
c        expm21redov2.pdf
c        expm22redov2.pdf
c to find corrections for M11*f1, M12*f2, M21*f1, M22*f2
c
c Note: x=xb+c1*alf+c2*alf^2+c3*alf^3+c4*alf^4
c       xi=dy+d1*alf+d2*alf^2+d3*alf^3+d4*alf^4
c       f=f0+f1*alf+f2*alf^2+f3*alf^3
c dx=xb-x0

c Note2: M11*u = mnum*f/r2 * 1/r1^2  + lognum*f/r2^3 * log(r1^2)
c
c        r1^2=|x-x0|^2 = (x-x0)^2+xi^2, xi=y-y0  dy=yb-y0
c        r2^2=         = (x+x0)^2+xi^2           dx=xb-x0
c
c        rho^2=d^2+calf^2
c
c Throughout use: x0>0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
C     INPUT
      integer ns,k0
      MODE alf(0:*),alfb,d,xjunk,
     +      xb,x0,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +      f10,f11,f12,f13,f20,f21,f22,f23
      logical largecorr
C     OUTPUT
      MODE corru,corrv
C3     LOCAL
      MODE corrs(nb),trapsing(nb),cm
      MODE corr1,corr2,corr3,corr4

      largecorr=.false.
      xjunk=1
      call compcoeffnew(xb,x0,d,c1,c2,c3,c4,dy,d1,d2,d3,d4,
     +                  xjunk,xjunk,cm)
      call compallcorrnew(alf,alfb,ns,d,cm,corrs,k0,largecorr,trapsing)
c      print*,alfb,k0,corrs(1)

      call corrqm11(x0,f10,f11,f12,f13,corrs,corr1)
      call corrqm12(x0,f20,f21,f22,f23,corrs,corr2)
      call corrqm21(x0,f10,f11,f12,f13,corrs,corr3)
      call corrqm22(x0,f20,f21,f22,f23,corrs,corr4)

      corru= - (corr1 + corr2)/pi8
      corrv= - (corr3 + corr4)/pi8

c      print*,x0,f10,f11,f12
      return
      end


      SUBROUTINE corrqm11(x0,f0,f1,f2,f3,corrs,corr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*)
c OUTPUT 
      MODE corr
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),small(0:6),
     +       piece1(0:6),piece2(0:6),piece3(0:6),
     +       corrqrho2,corrqlog,corr1,corr2
      MODE piece5(0:6),denom
      
c 1/RHO2 --------------------------------------------------
      lead = 1/x0

c Piece 1: small=delr^2=a2, fact=(x+x0)^2
      do j=0,3
        fact(j)=x2(j)+2*x1(j)*x0+ones(j)*x02
      enddo
      call prodrho2(fact,ir,a2,piece1)

c Piece 2: small=-xi^4=-b4,  fact=1
      do j=0,4
         small(j) = -b4(j)
      enddo
      call prodrho2(ones,ir,small,piece2)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)+piece2(j)
      enddo

      corr1=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = -1/(4*x0)

c Piece 1: small=delr=a1, fact=5x^4 +8x^3*x0 +6x^2*x0^2 +8x*x0^3 +5*x0^4
      do j=0,2
        fact(j)=  5*x4(j)     + 8*x3(j)*x0    + 6*x2(j)*x02
     +          + 8*x1(j)*x03 + 5*ones(j)*x04
      enddo
      call prodlog(fact,ir3,ones,piece1)

c Piece 2: small=xi^2=b2,  fact=12x^2 + 16x*x0 + 12x0^2
      do j=0,2
        fact(j)=12*x2(j)+16*x1(j)*x0+12*ones(j)*x02
      enddo
      call prodlog(fact,ir3,b2,piece2)

c Piece 4: small=7xi^4=7b4,  fact=1
      do j=0,2
         small(j) = 7*b4(j)
      enddo
      call prodlog(ones,ir3,small,piece3)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 1 (skip xi^2)
c     13*x^4 + 80*x^3*x0 + 134*x^2*x0^2 + 80*x*x0^3 + 13*x0^4
      fact(0)=13*x4(0) +80*x3(0)*x0 +134*x2(0)*x02 +80*x1(0)*x03 +13*x04
      denom=-8*r2**5  !must equal actual denom after mult by lead
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces 
      do j=0,2
        numer(j)=piece1(j)+piece2(j)+piece3(j)+piece5(j)
      enddo
      corr2=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      
      corr=corr1+corr2

      return
      end
      


      SUBROUTINE corrqm12(x0,f0,f1,f2,f3,corrs,corr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*)
c OUTPUT 
      MODE corr
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),small(0:6),
     +       piece1(0:6),piece2(0:6),
     +       corrqrho2,corrqlog,corr1,corr2
      MODE piece5(0:6),denom
      
c 1/RHO2 --------------------------------------------------
      lead = 2/x0

c Piece 1: small=delr*xi=ab1, fact=x(x+x0)
      do j=0,3
        fact(j)=x2(j)+x1(j)*x0
      enddo
      call prodrho2(fact,ir,ab1,piece1)

c Piece 2: small=xi^3=b3,  fact=x
      call prodrho2(x1,ir,b3,piece2)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)+piece2(j)
      enddo

      corr1=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = 1/(2*x0)

c Piece 1: small=xi=b1, fact=x^3 + 4x^2*x0 + 3x*x0^2
      do j=0,2
        fact(j)= x3(j) +4*x2(j)*x0 +3*x1(j)*x02
      enddo
      call prodlog(fact,ir3,b1,piece1)

c Piece 2: small=xi^3=b3,  fact=x
      call prodlog(x1,ir3,b3,piece2)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 1 (skip xi^2)
      fact(0)=b1(0)*x1(0)*(x2(0)+8*x1(0)*x0+7*x02+b2(0))
      denom=8*r2**5  !must equal actual denom after mult by lead
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces 
      do j=0,2
        numer(j)=piece1(j)+piece2(j)+piece5(j)
      enddo
      corr2=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      
      corr=corr1+corr2

      return
      end
      
      

      SUBROUTINE corrqm21(x0,f0,f1,f2,f3,corrs,corr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*)
c OUTPUT 
      MODE corr
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),small(0:6),
     +       piece1(0:6),piece2(0:6),
     +       corrqrho2,corrqlog,corr1,corr2
      MODE piece5(0:6),denom
      
c 1/RHO2 --------------------------------------------------
      lead = 2

c Piece 1: small=delr*xi=ab1, fact=(x+x0)
      do j=0,3
        fact(j)=x1(j)+ones(j)*x0
      enddo
      call prodrho2(fact,ir,ab1,piece1)

c Piece 2: small=-xi^3=b3,  fact=1
      do j=0,4
        small(j)=-b3(j)
      enddo
      call prodrho2(ones,ir,small,piece2)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)+piece2(j)
      enddo

      corr1=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = -1.q0/2

c Piece 1: small=xi=b1, fact=..
      do j=0,2
        fact(j)= 3*x2(j) +4*x1(j)*x0 + ones(j)*x02
      enddo
      call prodlog(fact,ir3,b1,piece1)

c Piece 2: small=xi^3
      call prodlog(ones,ir3,b3,piece2)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 1 (skip xi^2)
      fact(0)=b1(0)*(7*x2(0)+8*x1(0)*x0+x02+b2(0))
      denom=8*r2**5  !must equal actual denom after mult by lead
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces 
      do j=0,2
        numer(j)=piece1(j)+piece2(j)+piece5(j)
      enddo
      corr2=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      
      corr=corr1+corr2

      return
      end
      
     

      SUBROUTINE corrqm22(x0,f0,f1,f2,f3,corrs,corr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#include "coeff.h"
c INPUT
      MODE x0,f0,f1,f2,f3,corrs(*)
c OUTPUT 
      MODE corr
c LOCAL
      integer j
      MODE lead,numer(0:6),fact(0:6),small(0:6),
     +       piece1(0:6),piece2(0:6),
     +       corrqrho2,corrqlog,corr1,corr2
      MODE piece5(0:6),denom
      
c 1/RHO2 --------------------------------------------------
      lead = 4

c Piece 1: small=xi^2=b2, fact=x
      call prodrho2(x1,ir,b2,piece1)

c Sum Pieces
      do j=0,4
        numer(j)=piece1(j)
      enddo

      corr1=corrqrho2(lead,numer,f0,f1,f2,f3,e3,e4,corrs)

c LOG TERM --------------------------------------------------
c compute numerator for g*log need g to alf^2, leaving alf^3->O(d^4)error
c First compute num=lognum to 2nd order =num0+num1*alf+num2*alf^2
      lead = -1

c Piece 1: small=xi=b1, fact=..
      do j=0,2
        fact(j)= 2*x3(j)+4*x2(j)*x0+2*x1(j)*x02
      enddo
      call prodlog(fact,ir3,ones,piece1)

c Piece 2: small=xi^2 fact = 3x
      do j=0,2
        fact(j)= 3*x1(j)
      enddo
      call prodlog(fact,ir3,b2,piece2)

c Piece 5: ADD MLOGM terms for fourth order: rho1^2/rho2^2*log(rho1^2/rho2^2)
c     (d^2+cm^2*alf^2)*logarg*fact/denom (demom=mathematicadenom*r2^2)
c     these terms have E[.]=O(h^3)*O(fact)
c     keep fact and denom only to order 1 (skip xi^2)
      fact(0)=x1(0)*(4*x2(0)+8*x1(0)*x0+4*x02+7*b2(0))
      denom=8*r2**5  !must equal actual denom after mult by lead
      piece5(0)=dsq*fact(0)/denom
      piece5(1)=0
      piece5(2)=cmsq*fact(0)/denom

c Sum Pieces 
      do j=0,2
        numer(j)=piece1(j)+piece2(j)+piece5(j)
      enddo
      corr2=corrqlog(lead,numer,f0,f1,f2,f3,e3,e4,corrs)
      
      corr=corr1+corr2

      return
      end
      


