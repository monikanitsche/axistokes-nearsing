#include "flags.h"
      subroutine ellip(k2,F,E,n)
c
c     Computes the elliptic integrals F(k) and E(k)  using the
c     arithmetico - geometrical mean method (also see Burlisch
c     algorithm), where k2,F,E are vectors of length (0:n)
c
c     Input  :    k2  -  k**2
c
#include "params.h"
      integer j,l,n
      MODE a(0:nmax),b(0:nmax),c(0:nmax),E(0:*),
     +   fact(0:nmax),h(0:nmax),k2(0:*),sumc(0:nmax),F(0:*)
      logical flag
 
      do j = 0,n                      !vectorizes
         a(j) = 1
         b(j) = sqrt( 1 - k2(j) )
         sumc(j) = k2(j)
         fact(j) = 1
      enddo
 
      l = 0
 
800      continue
 
      flag = .false.
      l = l + 1
 
         do j = 0,n                  !vectorizes
            h(j) = a(j)
 
            c(j) = (a(j) - b(j))/2
            fact(j) = 2*fact(j)
            sumc(j) = sumc(j) + fact(j)*c(j)**2
 
            a(j) = (a(j) + b(j))/2
            if ( abs(h(j)-b(j)).gt.(ca*h(j))) flag = .true.
            b(j) = sqrt(h(j)*b(j))
        enddo
 
 
      if (flag) goto 800
 
c      print*,'l =',l
 
      do j = 0,n                     !vectorizes
         F(j) = pi2/(4*a(j))
         E(j) = (1 - sumc(j)/2) * F(j)
      enddo
 
      return
      end

      subroutine ellip1mink(k2,onemink2,F,E,n)
c
c     Computes the elliptic integrals F(k) and E(k)  using the
c     arithmetico - geometrical mean method (also see Burlisch
c     algorithm), where k2,F,E are vectors of length (0:n)
c
c     Input  :    k2  -  k**2
c
#include "params.h"
      integer j,l,n
      MODE a(0:nmax),b(0:nmax),c(0:nmax),E(0:*),
     +   fact(0:nmax),h(0:nmax),k2(0:*),sumc(0:nmax),F(0:*),
     +   onemink2(0:*)
      logical flag
 
      do j = 0,n                      !vectorizes
         a(j) = 1
c         b(j) = sqrt( 1 - k2(j) )
         b(j) = sqrt( onemink2(j) )
         sumc(j) = k2(j)
         fact(j) = 1
      enddo
 
      l = 0
 
800      continue
 
      flag = .false.
      l = l + 1
 
         do j = 0,n                  !vectorizes
            h(j) = a(j)
 
            c(j) = (a(j) - b(j))/2
            fact(j) = 2*fact(j)
            sumc(j) = sumc(j) + fact(j)*c(j)**2
 
            a(j) = (a(j) + b(j))/2
            if ( abs(h(j)-b(j)).gt.(ca*h(j))) flag = .true.
            b(j) = sqrt(h(j)*b(j))
        enddo
 
 
      if (flag) goto 800
 
c      print*,'l =',l
 
      do j = 0,n                     !vectorizes
         F(j) = pi2/(4*a(j))
         E(j) = (1 - sumc(j)/2) * F(j)
      enddo
 
      return
      end

