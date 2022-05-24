#include "flags.h"
      subroutine ftrans(x,y,n,coeff,sv,info)
c
c     Computes the Fourier Transform of the periodic 
c     vector (x+iy)(0,n) and stores it in coeff(0,2n)
c     where x is even y is odd
c 
      implicit none
      integer n,n2,j,j1,j2,info
      MODE x(0:*),y(0:*),sv(0:*),coeff(0:*)

      n2 = 2*n

      do j=0,n
         j1 = 2*j
         j2 = j1+1
         coeff(j1) = x(j)
         coeff(j2) = y(j)
      enddo

      do j=n2,n+1,-1
         j1 = 2*j
         j2 = j1+1
         coeff(j1) = coeff(2*(n2-j))
         coeff(j2) = -coeff(2*(n2-j)+1)
      enddo

      call dfft(n2,coeff,-1,sv)
      do j=0,2*n2-1
         coeff(j) = coeff(j)/n2
      enddo

c      if (info.eq.1) then
c      print*
c      print*,'b=['
c      do j=0,n2-1
c         j1 = 2*j
c         j2 = j1+1
c         write(*,'(i5,2e15.6)')j,coeff(j1),coeff(j2)
c      enddo
c      print*,'];'
c      endif
      return
      end
