C     Contains routines to compute derivatives of periodic functions
c     spectrally

#include "flags.h"
      subroutine deriv3odd(x,xp,xpp,xppp,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes derivatives of periodic functions x assuming that
c       x(0:2n) is odd
c     using spectral differentiation
c     By symmetry Fourier coefficients of y+ix all real
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,j1,j2,n,n2,jhold
      MODE x(0:*),xp(0:*),xpp(0:*),xppp(0:*),y(0:nmax)
      MODE coeff(0:4*nmax),sv(0:4*nmax),p1coeff(0:4*nmax),
     +       p2coeff(0:4*nmax),p3coeff(0:4*nmax),absval
 
      n2 = 2*n
      do j=0,n
        y(j)=0
      enddo
      call ftrans(y,x,n,coeff,sv,0)
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold =j
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) = 0 ! -jhold   *coeff(j2)  (by symmetry)
         p1coeff(j2) =  jhold   *coeff(j1)  
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = 0 ! -jhold**2*coeff(j2)
         p3coeff(j1) = 0 ! -jhold**3*coeff(j2)!should be +?
         p3coeff(j2) = -jhold**3*coeff(j1)
      enddo
      do j =n+1,n2-1
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold = j-n2
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) = 0 !-jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = 0 !-jhold**2*coeff(j2)
         p3coeff(j1) = 0 !-jhold**3*coeff(j2)
         p3coeff(j2) = -jhold**3*coeff(j1)
      enddo
 
      call dfft(n2,p1coeff,1,sv)
      call dfft(n2,p2coeff,1,sv)
      call dfft(n2,p3coeff,1,sv)
 
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         xp(j)   = p1coeff(j2)
         xpp(j)  = p2coeff(j2)
         xppp(j) = p3coeff(j2)
      enddo
 
      return
      end



      subroutine deriv3even(x,xp,xpp,xppp,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes derivatives of periodic functions x assuming that
c       x(0:2n) is even
c     using spectral differentiation
c     By symmetry Fourier coefficients of y+ix all real
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,j1,j2,n,n2,jhold
      MODE x(0:*),xp(0:*),xpp(0:*),xppp(0:*),y(0:nmax)
      MODE coeff(0:4*nmax),sv(0:4*nmax),p1coeff(0:4*nmax),
     +       p2coeff(0:4*nmax),p3coeff(0:4*nmax),absval
 
      n2 = 2*n
      do j=0,n
        y(j)=0
      enddo
      call ftrans(x,y,n,coeff,sv,0)
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold =j
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) = 0 ! -jhold   *coeff(j2)  (by symmetry)
         p1coeff(j2) =  jhold   *coeff(j1)  
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = 0 ! -jhold**2*coeff(j2)
         p3coeff(j1) = 0 ! -jhold**3*coeff(j2)!should be +?
         p3coeff(j2) = -jhold**3*coeff(j1)
      enddo
      do j =n+1,n2-1
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold = j-n2
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) = 0 !-jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = 0 !-jhold**2*coeff(j2)
         p3coeff(j1) = 0 !-jhold**3*coeff(j2)
         p3coeff(j2) = -jhold**3*coeff(j1)
      enddo
 
      call dfft(n2,p1coeff,1,sv)
      call dfft(n2,p2coeff,1,sv)
      call dfft(n2,p3coeff,1,sv)
 
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         xp(j)   = p1coeff(j1)
         xpp(j)  = p2coeff(j1)
         xppp(j) = p3coeff(j1)
      enddo
 
      return
      end


      subroutine deriv1even(x,xp,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes derivatives of periodic functions x assuming that
c       x(0:2n) is even
c     using spectral differentiation
c     By symmetry Fourier coefficients of y+ix all real
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,j1,j2,n,n2,jhold
      MODE x(0:*),xp(0:*),y(0:nmax)
      MODE coeff(0:4*nmax),sv(0:4*nmax),p1coeff(0:4*nmax),absval
 
      n2 = 2*n
      do j=0,n
        y(j)=0
      enddo
      call ftrans(x,y,n,coeff,sv,0)
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold =j
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) = 0 ! -jhold   *coeff(j2)  (by symmetry)
         p1coeff(j2) =  jhold   *coeff(j1)  
      enddo
      do j =n+1,n2-1
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold = j-n2
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) = 0 !-jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
      enddo
 
      call dfft(n2,p1coeff,1,sv)
 
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         xp(j)   = p1coeff(j1)
      enddo
 
      return
      end

 
 
      subroutine compderiv2(y,yp,ypp,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes derivatives of periodic y(alfa) assuming
c        y: even 
c     using spectral differentiation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,j1,j2,n,n2,jhold
      MODE y(0:*),yp(0:*),ypp(0:*)
      MODE coeff(0:4*nmax),sv(0:4*nmax),p1coeff(0:4*nmax),
     +       p2coeff(0:4*nmax),!p3coeff(0:4*nmax),p4coeff(0:4*nmax),
     +       absval,x(0:nmax)
 
      do j=0,n
         x(j) =0
      enddo
      n2 = 2*n
      call ftrans(y,x,n,coeff,sv,0)
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold =j
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) =  0 ! jhold   *coeff(j2) should be -?
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) =  0 ! -jhold**2*coeff(j2) should be +?
      enddo
      do j =n+1,n2-1
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold = j-n2
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) =  0 !-jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) =  0 !jhold**2*coeff(j2)
      enddo
 
      call dfft(n2,p1coeff,1,sv)
      call dfft(n2,p2coeff,1,sv)
 
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         yp(j)  = p1coeff(j1)
         ypp(j) = p2coeff(j1)
      enddo
 
      return
      end
 
 
      subroutine compderiv3(x,xp,y,yp,ypp,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes derivatives of periodic functions x,y assuming that
c       x:odd
c       y:even
c     using spectral differentiation. 
c     Because of symmetry y+ix=sum c_k e^{ikx} has purely real FCoeff c_k
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,j1,j2,n,n2,jhold
      MODE x(0:*),y(0:*),xp(0:*),yp(0:*),ypp(0:*)
      MODE coeff(0:4*nmax),sv(0:4*nmax),p1coeff(0:4*nmax),
     +       p2coeff(0:4*nmax),absval
 
      n2 = 2*n
      call ftrans(y,x,n,coeff,sv,0)
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold =j
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) = 0 !-jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = 0 !-jhold**2*coeff(j2)
      enddo
      do j =n+1,n2-1
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold = j-n2
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) = 0 !-jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = 0 !-jhold**2*coeff(j2)
      enddo
 
      call dfft(n2,p1coeff,1,sv)
      call dfft(n2,p2coeff,1,sv)
 
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         yp(j)  = p1coeff(j1)
         xp(j)  = p1coeff(j2)
         ypp(j) = p2coeff(j1)
      enddo
 
      return
      end
 

 
      subroutine compderiv1(x,xp,xpp,xppp,xpppp,y,yp,ypp,yppp,ypppp,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes derivatives of periodic functions x,y assuming that
c       x:odd
c       y:even
c     using spectral differentiation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,j1,j2,n,n2,jhold
      MODE x(0:*),y(0:*),xp(0:*),yp(0:*),xpp(0:*),ypp(0:*),xppp(0:*),
     + yppp(0:*),xpppp(0:*),ypppp(0:*)
      MODE coeff(0:4*nmax),sv(0:4*nmax),p1coeff(0:4*nmax),
     +       p2coeff(0:4*nmax),p3coeff(0:4*nmax),p4coeff(0:4*nmax),
     +       absval
 
      n2 = 2*n
      call ftrans(y,x,n,coeff,sv,0)
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold =j
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) =  jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = -jhold**2*coeff(j2)
         p3coeff(j1) = -jhold**3*coeff(j2)
         p3coeff(j2) = -jhold**3*coeff(j1)
         p4coeff(j1) =  jhold**4*coeff(j1)
         p4coeff(j2) =  jhold**4*coeff(j2)
      enddo
      do j =n+1,n2-1
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold = j-n2
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) =  jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = -jhold**2*coeff(j2)
         p3coeff(j1) = -jhold**3*coeff(j2)
         p3coeff(j2) = -jhold**3*coeff(j1)
         p4coeff(j1) =  jhold**4*coeff(j1)
         p4coeff(j2) =  jhold**4*coeff(j2)
      enddo
 
      call dfft(n2,p1coeff,1,sv)
      call dfft(n2,p2coeff,1,sv)
      call dfft(n2,p3coeff,1,sv)
      call dfft(n2,p4coeff,1,sv)
 
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         yp(j)  = p1coeff(j1)
         xp(j)  = p1coeff(j2)
         ypp(j) = p2coeff(j1)
         xpp(j) = p2coeff(j2)
         yppp(j) = p3coeff(j1)
         xppp(j) = p3coeff(j2)
         ypppp(j) = p4coeff(j1)
         xpppp(j) = p4coeff(j2)
      enddo
 
      return
      end

 
 
      subroutine compderiv4(x,xp,xpp,y,yp,ypp,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes derivatives of periodic functions x,y assuming that
c       x:odd
c       y:even
c     using spectral differentiation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
      integer j,j1,j2,n,n2,jhold
      MODE x(0:*),y(0:*),xp(0:*),yp(0:*),xpp(0:*),ypp(0:*)
      MODE coeff(0:4*nmax),sv(0:4*nmax),p1coeff(0:4*nmax),
     +       p2coeff(0:4*nmax),absval
 
      n2 = 2*n
      call ftrans(y,x,n,coeff,sv,0)
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold =j
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) =  jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = -jhold**2*coeff(j2)
c         p3coeff(j1) = -jhold**3*coeff(j2)
c         p3coeff(j2) = -jhold**3*coeff(j1)
c         p4coeff(j1) =  jhold**4*coeff(j1)
c         p4coeff(j2) =  jhold**4*coeff(j2)
      enddo
      do j =n+1,n2-1
         j1 = 2*j
         j2 = j1+1
         absval = sqrt(coeff(j1)**2+coeff(j2)**2)
         jhold = j-n2
         if (absval.lt.flevelx) jhold=0
         p1coeff(j1) =  jhold   *coeff(j2)
         p1coeff(j2) =  jhold   *coeff(j1)
         p2coeff(j1) = -jhold**2*coeff(j1)
         p2coeff(j2) = -jhold**2*coeff(j2)
c         p3coeff(j1) = -jhold**3*coeff(j2)
c         p3coeff(j2) = -jhold**3*coeff(j1)
c         p4coeff(j1) =  jhold**4*coeff(j1)
c         p4coeff(j2) =  jhold**4*coeff(j2)
      enddo
 
      call dfft(n2,p1coeff,1,sv)
      call dfft(n2,p2coeff,1,sv)
c      call dfft(n2,p3coeff,1,sv)
c      call dfft(n2,p4coeff,1,sv)
 
      do j =0,n
         j1 = 2*j
         j2 = j1+1
         yp(j)  = p1coeff(j1)
         xp(j)  = p1coeff(j2)
         ypp(j) = p2coeff(j1)
         xpp(j) = p2coeff(j2)
c         yppp(j) = p3coeff(j1)
c         xppp(j) = p3coeff(j2)
c         ypppp(j) = p4coeff(j1)
c         xpppp(j) = p4coeff(j2)
      enddo
 
      return
      end

 
 


