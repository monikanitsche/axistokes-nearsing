#include "flags.h"
      SUBROUTINE compallcorrnew(alf,alfp,n,d,calf,corr,k0,tooclose,
     + trapsing)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Compute all corrections for all basisfunctions Hjk
c nb : number of basis functions, in params.h currently=35
c
c TESTED TO BE CORRECT! checked that convergence is O(h^4) and that
c values are correct, meaning that the scaling fscale is correct
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer n,k0
      MODE alf(0:*),alfp,d,calf
      logical tooclose
c OUTPUT
      MODE corr(*)
c LOCAL
      integer j,kl,kr,k
      MODE fact,h,delu,delu2,delu4,ul,ur,trap(nb),Fl(nb),Fr(nb),
     +   fpl(nb),fpr(nb),fp3l(nb),fp3r(nb),fscale(nb),inth(nb),corrtrap,
     +   trapsing(nb),trapnok(nb),one,half

      fact=d/calf

c      h=pi/n   !CAREFUL, CHANGED ON 2/10/21
      h=alf(n)/n
      delu=h/fact
      delu2=delu**2
      delu4=delu**4

c     initialize trapezoid integrals for all functions
      do j=1,31
        trap(j)=0
        trapsing(j)=0
      enddo

      one=1
      half=one/2
      kl=0
      kr=n
      if ((k0.eq.0).or.(k0.eq.n)) then
        call trapcomponent((alf(kl)-alfp)/fact,trap,half,d)
        call trapcomponent((alf(kr)-alfp)/fact,trap,half,d)
        do k=kl+1,kr-1  
           call trapcomponent((alf(k)-alfp)/fact,trap,one,d)
        enddo
        do j=1,31
          trapnok(j)=trap(j)*delu
          trap(j)=trapnok(j)
        enddo
      else
c       compute trapezoid approximation of all integrals 
        call trapcomponent((alf(kl)-alfp)/fact,trap,half,d)
        call trapcomponent((alf(kr)-alfp)/fact,trap,half,d)
c        do k=kl+1,kr-1  
c           call trapcomponent((alf(k)-alfp)/fact,trap,1.q0)
c        enddo
        do k=kl+1,k0-1
           call trapcomponent((alf(k)-alfp)/fact,trap,one,d)
        enddo
        do k=k0+1,kr-1
           call trapcomponent((alf(k)-alfp)/fact,trap,one,d)
        enddo
        k=k0
        call trapcomponent((alf(k)-alfp)/fact,trapsing,one,d)
        do j=1,31
          trapnok(j)=trap(j)*delu
          trap(j)=(trap(j)+trapsing(j))*delu
        enddo
      endif

c     compute corrections
      ul=(alf(kl)-alfp)/fact
      ur=(alf(kr)-alfp)/fact
      call allvals(ul,Fl,fpl,fp3l,d)
      call allvals(ur,Fr,fpr,fp3r,d)

      do j=1,31
        inth(j)    = Fr(j)-Fl(j)
        corrtrap= -delu2/12 *( fpr(j)- fpl(j))
     +            +delu4/720*(fp3r(j)-fp3l(j))
        corr(j) = inth(j)-(trap(j)+corrtrap)
c overwrite if too close!        
        if (tooclose) corr(j) = inth(j)-(trapnok(j)+corrtrap)
      enddo

c     scale corrections by fscale, given by change of variables: 
c         int f(alf)dalf  =  fscale*int f(u)du
c
c     for f=H0k, fscale = (d/c)^(k+1)
c     for f=Hjk, fscale = d^(k+1-2j)/c^(k+1)

      fscale(1) = d/calf              !H00      logu k=0    h
      fscale(2) = (d/calf)**2         !H01    u logu k=1    h2
      fscale(3) = (d/calf)**3         !H02   u^2logu k=2    h3
      fscale(31)= (d/calf)**4         !H03   u^3logu k=2    h3
      fscale(4) = 1   /(d*calf)       !H10    1/r2 j=1 k=0 M h 
      fscale(5) = 1   /calf**2        !H11    u/r2 j=1 k=1 M h 
      fscale(6) = d   /calf**3        !H12   u2/r2 j=1 k=2 M h
      fscale(7) = d**2/calf**4        !H13   u3/r2 j=1 k=3 M hd
      fscale(8) = d**3/calf**5        !H14   u4/r2 j=1 k=4 M hd2
      fscale(9) = 1   /(d**3*calf)    !H20    1/r4 j=2 k=0 Q h/d 
      fscale(10)= 1   /(d**2*calf**2) !H21    u/r4 j=2 k=1 Q h/d 
      fscale(11)= 1   /(d*calf**3)    !H22   u2/r4 j=2 k=2 Q h/d 
      fscale(12)= 1   /calf**4        !H23   u3/r4 j=2 k=3 Q h/d M hd
      fscale(13)= d   /calf**5        !H24   u4/r4 j=2 k=4 Q h   M hd
      fscale(14)= d**2/calf**6        !H25   u5/r4 j=2 k=5 Q hd  M hd
      fscale(15)= d**3/calf**7        !H26   u6/r4 j=2 k=6 Q hd2 M hd2
      fscale(16)= 1   /(d**2*calf**4) !H33   u3/r6 j=3 k=3 Q h 
      fscale(17)= 1   /(d*calf**5)    !H34   u4/r6 j=3 k=4 Q h 
      fscale(18)= 1   /calf**6        !H35   u5/r6 j=3 k=5 Q h 
      fscale(19)= d   /calf**7        !H36   u6/r6 j=3 k=6 Q h   M hd2
      fscale(20)= d**2/calf**8        !H37   u7/r6 j=3 k=7 Q hd  M hd2
      fscale(21)= d**3/calf**9        !H38   u8/r6 j=3 k=8 Q hd2 M hd2
      fscale(22)= 1   /(d*calf**7)    !H46   u6/r8 j=4 k=6 Q hd 
      fscale(23)= 1   /calf**8        !H47   u7/r8 j=4 k=7 Q hd 
      fscale(24)= d   /calf**9        !H48   u8/r8 j=4 k=8 Q hd 
      fscale(25)= d**2/calf**10       !H49   u9/r8 j=4 k=9 Q hd 
      fscale(26)= d**3/calf**11       !H410 u10/r8 j=4 k=10  Q hd2
      fscale(27)=    1/calf**10       !H59   u9/r10 j=5 k=9  Q hd2 
      fscale(28)=    d/calf**11       !H510 u10/r10 j=5 k=10 Q hd2 
      fscale(29)= d**2/calf**12       !H511 u11/r10 j=5 k=11 Q hd2 
      fscale(30)= d**3/calf**13       !H512 u12/r10 j=5 k=12 Q hd2

c h^4 method: 30 corrections
c h^3 method: 21 corrections
c h^2 method: 13 corrections

      do j=1,31
        corr(j)=fscale(j)*corr(j)
c        corrnok(j)=fscale(j)*corrnok(j)
c        print*,j,inth(j)*fscale(j)
        trapsing(j)=trapsing(j)*delu*fscale(j)
      enddo

      return
      end


      SUBROUTINE allvals(u,F,fp,fp3,d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computes F(u), fp(u), fp3(u) for all basis functions 
c input : u (scalar)
c output : F,fp,fp3 (vector)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c INPUT/OUTPUT
      MODE u,F(*),fp(*),fp3(*)  ,d
c LOCAL
      MODE u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,
     +       arg,arg2,arg3,arg4,arg5,arg6,arg7,arg8,
     +       atanu,logarg

      u2=u**2
      u3=u2*u
      u4=u2*u2
      u5=u2*u3
      u6=u2*u4
      u7=u3*u4
      u8=u4*u4
      u9=u8*u
      u10=u9*u
      u11=u10*u
      arg=1+u2
      arg2=arg**2
      arg3=arg**3
      arg4=arg**4
      arg5=arg**5
      arg6=arg**6
      arg7=arg**7
      arg8=arg**8
      atanu=atan(u)
      logarg=log(arg)

c H00c
      F(1)  = 2*atanu+u*(logarg-2)  + 2*log(d)*u
      fp(1) = 2*u/arg               
      fp3(1)= 4*u*(u2-3)/arg3
c H01c
      F(2)  = (-u2+arg*logarg)/2    + 2*log(d)*u2/2
      fp(2) =  2*u2/arg+logarg      + 2*log(d)
      fp3(2)= -2*(u4+6*u2-3)/arg3
c H02c
      F(3) = (6*u+u3*(-2+3*logarg)-6*atanu)/9  + 2*log(d)*u3/3
      fp(3)= 2*u*(u2/arg+logarg)               + 4*log(d)*u
      fp3(3) = 4*u*(6+3*u2+u4)/arg3
c H03
      F(31)  = (-(u2*(-2 + u2)) + 2*(-1 + u4)*logarg)/8
      fp(31) =  2*u4/arg + 3*u2*logarg 
      fp3(31)= (54*u2 + 60*u4 + 22*u6 + 6*arg3*logarg)/arg3
c H10c  
      F(4) = atanu
      fp(4)= -2*u/arg2
      fp3(4)= -24*u*(u2-1)/arg4
c H11c 
      F(5) = logarg/2
      fp(5)= (1-u2)/arg2
      fp3(5)= -6*(1-6*u2+u4)/arg4
c H12c 
      F(6) = u-atanu
      fp(6)= 2*u/arg2
      fp3(6)= 24*u*(u2-1)/arg4
c H13c
      F(7) = (u2-logarg)/2
      fp(7)= u2*(3+u2)/arg2
      fp3(7)= 6*(1-6*u2+u4)/arg4
c H14c
      F(8) = u/3*(u2-3)+atanu
      fp(8)= 2*u3*(2+u2)/arg2
      fp3(8)= -24*u*(u2-1)/arg4
c H20c
      F(9) =(u/arg+atanu)/2
      fp(9)= -4*u/arg3
      fp3(9)= 24*u*(3-5*u2)/arg5
c H21c
      F(10) = -1/arg/2
      fp(10)= (1-3*u2)/arg3
      fp3(10)= -12*(1-10*u2+5*u4)/arg5
c H22c
      F(11) = (-u/arg+atanu)/2
      fp(11)= -2*u*(u2-1)/arg3
      fp3(11)= -24*u*(2-5*u2+u4)/arg5
c H23c
      F(12) = (1/arg +logarg)/2
      fp(12)= -u2*(u2-3)/arg3
      fp3(12)= -6*(u6-15*u4+15*u2-1)/arg5
c H24c
      F(13) = u+u/arg/2-3*atanu/2
      fp(13)= 4*u3/arg3
      fp3(13)= 24*u*(1-5*u2+2*u4)/arg5
c H25c
      F(14) = (u2-1/arg-2*logarg)/2
      fp(14)= u4*(5+u2)/arg3
      fp3(14)= 12*u2*(5-10*u2+u4)/arg5
c H26c
      F(15) = -u/6*(12-2*u2+3/arg)+5*atanu/2
      fp(15)= 2*u5*(3+u2)/arg3
      fp3(15)= 24*u3*(5-3*u2)/arg5
c H33c
      F(16) = -(1+2*u2)/arg2/4
      fp(16)= -3*u2*(u2-1)/arg4
      fp3(16)= -6*(10*u6-45*u4+24*u2-1)/arg6
c H34c
      F(17) = (-u*(3+5*u2)/arg2 +3*atanu)/8
      fp(17)= -2*u3*(u2-2)/arg4
      fp3(17)= -24*u*(u6-9*u4+9*u2-1)/arg6
c H35c
      F(18) = ( (3+4*u2)/arg2 +2*logarg)/4
      fp(18)= -u4*(u2-5)/arg4
      fp3(18)= -6*u2*(u6-24*u4+45*u2-10)/arg6
c H36c
      F(19) = ( u*(15+25*u2+8*u4)/arg2 -15*atanu)/8
      fp(19)= 6*u5/arg4
      fp3(19)= 24*u3*(5-12*u2+3*u4)/arg6
c H37c
      F(20) = ( 2*u2-(5+6*u2)/arg2 -6*logarg )/4
      fp(20)= u6*(7+u2)/arg4
      fp3(20)= 6*u4*(35-42*u2+3*u4)/arg6
c H38c
      F(21)= (u*(-72+8*u2+6/arg2-39/arg) + 105*atanu)/24
      fp(21)=2*u7*(4+u2)/arg4
      fp3(21)=-48*u5*(-7+3*u2)/arg6
c H46c
      F(22) = -u*(15+40*u2+33*u4)/arg3/48+5*atanu/16
      fp(22)= -2*u5*(u2-3)/arg5
      fp3(22)= -24*u3*(-5+21*u2-13*u4+u6)/arg7
c H47c
      F(23) = (11+27*u2+18*u4)/12/arg3 +logarg/2
      fp(23)= -u6*(u2-7)/arg5
      fp3(23)= -6*u4*(-35+91*u2-33*u4+u6)/arg7
c H48c
      F(24) = (u*(105+280*u2+231*u4+48*u6)/arg3-105*atanu)/48
      fp(24)= 8*u7/arg5
      fp3(24)= 48*u5*(7-11*u2+2*u4)/arg7
c H49c
      F(25) = (-13-27*u2-9*u4+9*u6+3*u8)/arg3/6- 2*logarg 
      fp(25)= u8*(9+u2)/arg5
      fp3(25)= 24*u6*(21-18*u2+u4)/arg7
c H410c
      F(26) = ((u*(-315 - 840*u2 - 693*u4 - 144*u6 + 16*u8))/arg3 + 
     +    315*atanu)/48
      fp(26)= (2*u9*(5 + u2))/arg5
      fp3(26)= -240*u7*(-3 + u2)/arg7
c H59
      F(27) = (25 + 88*u2 + 108*u4 + 48*u6 + 12*arg4*logarg)/
     +  (24*arg4)
      fp(27)= -u8*(-9 + u2)/arg6
      fp3(27)= -6*u6*(-84 + 153*u2 - 42*u4 + u6)/arg8
c H510
      F(28) = (u*(315 + 1155*u2 + 1533*u4 + 837*u6 + 128*u8)/arg4 - 
     +    315*atanu)/128
      fp(28)= 10*u9/arg6
      fp3(28)= 120*u7*(6 - 7*u2 + u4)/arg8
c H511
      F(29) = -(77 + 248*u2 + 252*u4 + 48*u6 - 48*u8 - 
     +      12*u10)/arg4/24 - 5*logarg/2
      fp(29)= u10*(11 + u2)/arg6
      fp3(29)= 30*u8*(33 - 22*u2 + u4)/arg8
c H512
      F(30) = (u*(-3465 -12705*u2 -16863*u4 -9207*u6 -1408*u8 +128*u10)/
     -     arg4 + 3465*atanu)/384
      fp(30)= 2*u11*(6 + u2)/arg6
      fp3(30)= 120*u9*(11 - 3*u2)/arg8

      return
      end


      SUBROUTINE trapcomponent(u,trap,wt,d)
c computes one component f_(uj)*wt  where u=uj in the trapezoid rule
c for various functions f
      implicit none
c INPUT      
      MODE u,wt,d
c OUTPUT      
      MODE trap(*)
c LOCAL
      MODE u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,
     +   arg,arg1,arg2,arg3,arg4,arg5,logarg

      u2=u**2
      u3=u2*u
      u4=u3*u
      u5=u4*u
      u6=u5*u
      u7=u6*u
      u8=u7*u
      u9=u8*u
      u10=u9*u
      u11=u10*u
      u12=u11*u

      arg=1+u2
      arg1=arg/wt
      arg2=arg**2/wt
      arg3=arg**3/wt
      arg4=arg**4/wt
      arg5=arg**5/wt
c      logarg=log(arg)*wt
      logarg=log(d**2*arg)*wt

      trap(1) = trap(1) +   logarg  !H00
      trap(2) = trap(2) + u*logarg  !H01
      trap(3) = trap(3) + u2*logarg !H02
      trap(31)= trap(31)+ u3*logarg !H03
      trap(4) = trap(4) + 1 /arg1   !H10
      trap(5) = trap(5) + u /arg1   !H11
      trap(6) = trap(6) + u2/arg1   !H12
      trap(7) = trap(7) + u3/arg1   !H13
      trap(8) = trap(8) + u4/arg1   !H14
      trap(9) = trap(9) + 1 /arg2   !H20
      trap(10)= trap(10)+ u /arg2   !H21
      trap(11)= trap(11)+ u2/arg2   !H22
      trap(12)= trap(12)+ u3/arg2   !H23
      trap(13)= trap(13)+ u4/arg2   !H24
      trap(14)= trap(14)+ u5/arg2   !H25
      trap(15)= trap(15)+ u6/arg2   !H26
      trap(16)= trap(16)+ u3/arg3   !H33 
      trap(17)= trap(17)+ u4/arg3   !H34
      trap(18)= trap(18)+ u5/arg3   !H35
      trap(19)= trap(19)+ u6/arg3   !H36
      trap(20)= trap(20)+ u7/arg3   !H37
      trap(21)= trap(21)+ u8/arg3   !H38
      trap(22)= trap(22)+ u6/arg4   !H46
      trap(23)= trap(23)+ u7/arg4   !H47
      trap(24)= trap(24)+ u8/arg4   !H48
      trap(25)= trap(25)+ u9/arg4   !H49
      trap(26)= trap(26)+u10/arg4   !H410 
      trap(27)= trap(27)+ u9/arg5   !H59  
      trap(28)= trap(28)+u10/arg5   !H510 
      trap(29)= trap(29)+u11/arg5   !H511 
      trap(30)= trap(30)+u12/arg5   !H512 

      return
      end
