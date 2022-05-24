#include "flags.h"
      SUBROUTINE compallcorrend(alf,alfp,d,cm,n,corr,corrl)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c we use CHANGE OF VARIABLES u=alf/fact, fact=d/cm
c integrals of interest are of the form 
c
c         int alf^k/sqrt(d^2+cm^2 alf^2)^l dalf  
c
c        =  fact*fact^k/d^l int u^k/sqrt(1+u2)^l du
c
c here we compute errors T[Hkl]-I[Hkl]
c   corrl(8): uses 3rd derivative
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "params.h"
c INPUT
      integer n
      MODE alf(0:*),alfp,d,cm,corr(*),corrl(*)
c LOCAL
      integer j,jl,jr,k
      MODE fact,fact3,trapl
      MODE ul,ur,u,arg,fpl,fpr,Fl,Fr,h,h2,h4,trap(20),intf,fp3l,fp3r
      MODE u2,u3,u5,u7,rt3,rt5,rt7,rt9,rt11,u2l,u4l,u6l,u2r,u4r,u6r, !u8r,u8l
     +   rtl,rt3l,rt5l,rt7l,rt9l,rt11l,rt13l,u8l,u8r,u9,
     +   rtr,rt3r,rt5r,rt7r,rt9r,rt11r,rt13r,argl,argr
      MODE f13,f33,f53,f15,f35,f55,f75,f17,f37,f57,f77,
     +       f19,f39,f59,f79,
     +       f511,f711,f99,f911

      fact=d/cm

      h=pi/n
      do j=1,19
        trap(j)=0
      enddo

      jl=0
      jr=n

c add left endpoint
      k=jl
        u  =(alf(k)-alfp)/fact
        u2 =u**2
        u3 =u**3
        u5 =u**5
        u7 =u**7
        u9 =u**9
        arg =1+u2
        rt3 =2*sqrt(arg)**3
        rt5 =2*sqrt(arg)**5
        rt7 =2*sqrt(arg)**7
        rt9 =2*sqrt(arg)**9
        rt11=2*sqrt(arg)**11
        trap(1) = trap(1) +  u/rt3      ! u/rt3
        trap(2) = trap(2) + u3/rt3      !u3/rt3
        trap(19)= trap(19)+ u5/rt3      !u5/rt3
        trap(3) = trap(3) +  u/rt5      ! u/rt5
        trap(4) = trap(4) + u3/rt5      !u3/rt5
        trap(5) = trap(5) + u5/rt5      !u5/rt5
        trap(16)= trap(16)+ u7/rt5      !u7/rt5
        trap(6) = trap(6) +  u/rt7      ! u/rt7
        trap(7) = trap(7) + u3/rt7      !u3/rt7
        trap(8) = trap(8) + u5/rt7      !u5/rt7
        trap(9) = trap(9) + u7/rt7      !u7/rt7
        trap(10)= trap(10)+  u/rt9      ! u/rt9
        trap(11)= trap(11)+ u3/rt9      !u3/rt9
        trap(12)= trap(12)+ u5/rt9      !u5/rt9
        trap(13)= trap(13)+ u7/rt9      !u7/rt9
        trap(17)= trap(17)+ u9/rt9      !u7/rt9
        trap(14)= trap(14)+ u5/rt11     !u5/rt11
        trap(15)= trap(15)+ u7/rt11     !u7/rt11
        trap(18)= trap(18)+ u9/rt11     !u9/rt11
c add right endpoint
      k=jr
        u=(alf(k)-alfp)/fact
        u2=u**2
        u3=u**3
        u5=u**5
        u7=u**7
        u9=u**9
        arg=1+u2
        rt3 =2*sqrt(arg)**3
        rt5 =2*sqrt(arg)**5
        rt7 =2*sqrt(arg)**7
        rt9 =2*sqrt(arg)**9
        rt11=2*sqrt(arg)**11
        trap(1) = trap(1)+ u/rt3      ! u/(1+u2)^3/2
        trap(2) = trap(2)+u3/rt3      !u3/(1+u2)^3/2
        trap(19)= trap(19)+ u5/rt3      !u5/rt3
        trap(3) = trap(3)+ u/rt5      ! u/(1+u2)^5/2
        trap(4) = trap(4)+u3/rt5      !u3/(1+u2)^5/2
        trap(5) = trap(5)+u5/rt5      !u5/(1+u2)^5/2
        trap(16)= trap(16)+u7/rt5     !NEW
        trap(6) = trap(6)+u/rt7       !NEW
        trap(7) = trap(7)+u3/rt7      !NEW
        trap(8) = trap(8)+u5/rt7      !u5/(1+u2)^7/2
        trap(9) = trap(9)+u7/rt7      !u7/(1+u2)^7/2
        trap(10)= trap(10)+u/rt9      !NEW
        trap(11)= trap(11)+u3/rt9     !NEW
        trap(12)= trap(12)+u5/rt9     !NEW
        trap(13)= trap(13)+u7/rt9     !NEW
        trap(17)= trap(17)+ u9/rt9    !u9/rt9
        trap(14)= trap(14)+u5/rt11    !NEW
        trap(15)= trap(15)+u7/rt11    !NEW
        trap(18)= trap(18)+ u9/rt11     !u9/rt11
      
c add interior points
      do k=jl+1,jr-1
        u=(alf(k)-alfp)/fact
        u2=u**2
        u3=u**3
        u5=u**5
        u7=u**7
        u9=u**9
        arg=1+u2
        rt3 =sqrt(arg)**3
        rt5 =sqrt(arg)**5
        rt7 =sqrt(arg)**7
        rt9 =sqrt(arg)**9
        rt11=sqrt(arg)**11
        trap(1) = trap(1) + u/rt3      ! u/(1+u2)^3/2
        trap(2) = trap(2) +u3/rt3      !u3/(1+u2)^3/2
        trap(3) = trap(3) + u/rt5      ! u/(1+u2)^5/2
        trap(4) = trap(4) +u3/rt5      !u3/(1+u2)^5/2
        trap(19)= trap(19)+ u5/rt3      !u5/rt3
        trap(5) = trap(5) +u5/rt5      !u5/(1+u2)^5/2
        trap(16)= trap(16)+u7/rt5     !NEW
        trap(6) = trap(6) + u/rt7      !NEW
        trap(7) = trap(7) +u3/rt7      !NEW
        trap(8) = trap(8) +u5/rt7      !u5/(1+u2)^7/2
        trap(9) = trap(9) +u7/rt7      !u7/(1+u2)^7/2
        trap(10)= trap(10)+ u/rt9      !NEW
        trap(11)= trap(11)+u3/rt9      !NEW
        trap(12)= trap(12)+u5/rt9      !NEW
        trap(13)= trap(13)+u7/rt9      !NEW
        trap(17)= trap(17)+ u9/rt9    !u9/rt9
        trap(14)= trap(14)+u5/rt11     !NEW
        trap(15)= trap(15)+u7/rt11     !NEW
        trap(18)= trap(18)+ u9/rt11     !u9/rt11
      enddo
      do j=1,19
        trap(j)=trap(j)*h
      enddo

      h2=h**2
      h4=h2**2
      fact3=fact**3
      ul=(alf(jl)-alfp)/fact
        u2l=ul**2
        u4l=ul**4
        u6l=ul**6
        u8l=ul**8
        argl=1+u2l
        rtl=sqrt(argl)
        rt3l=rtl*argl
        rt5l=rt3l*argl
        rt7l=rt5l*argl
        rt9l=rt7l*argl
        rt11l=rt9l*argl
        rt13l=rt11l*argl
      ur=(alf(jr)-alfp)/fact
        u2r=ur**2
        u4r=ur**4
        u6r=ur**6
        u8r=ur**8
        argr=1+u2r
        rtr=sqrt(argr)
        rt3r=rtr*argr
        rt5r=rt3r*argr
        rt7r=rt5r*argr
        rt9r=rt7r*argr
        rt11r=rt9r*argr
        rt13r=rt11r*argr

c     u/rt3
      fpl= -3*u2l/rt5l+1/rt3l
      fp3l= -3*(3-24*u2l+8*u4l)/rt9l
      Fl = -1/rtl
      fpr= -3*u2r/rt5r+1/rt3r
      fp3r= -3*(3-24*u2r+8*u4r)/rt9r
      Fr = -1/rtr
      intf=Fr-Fl
      trap(1) = trap(1)-h2/12*(fpr-fpl)/fact  
     +  +h4/720*(fp3r-fp3l)/fact3
      corrl(1) = intf*fact-trap(1)
c     u3/rt3
      fpl= 3*u2l/rt5l
      fp3l= (6-63*u2l+36*u4l)/rt9l 
      Fl = (2+u2l)/rtl
      fpr= 3*u2r/rt5r
      fp3r= (6-63*u2r+36*u4r)/rt9r 
      Fr = (2+u2r)/rtr
      intf=Fr-Fl
      trap(2) = trap(2)-h2/12*(fpr-fpl)/fact  
     +  +h4/720*(fp3r-fp3l)/fact3
      corrl(2) = intf*fact-trap(2)

c     u5/rt3
      fpl= u4l*(5+2*u2l)/rt5l
      fp3l= (60*u2l-45*u4l)/rt9l 
      Fl = (-8-4*u2l+u4l)/(3*rtl)
      fpr= u4r*(5+2*u2r)/rt5r
      fp3r= (60*u2r-45*u4r)/rt9r 
      Fr = (-8-4*u2r+u4r)/(3*rtr)
      intf=Fr-Fl
      trap(19) = trap(19)-h2/12*(fpr-fpl)/fact  
     +  +h4/720*(fp3r-fp3l)/fact3
      corrl(19) = intf*fact-trap(19)

c     u/rt5
      fpl= (1-4*u2l)/rt7l 
      fp3l= -15*(1-12*u2l+8*u4l)/rt11l 
      Fl = -1/rt3l/3
      fpr= (1-4*u2r)/rt7r 
      fp3r= -15*(1-12*u2r+8*u4r)/rt11r 
      Fr = -1/rt3r/3
      intf=Fr-Fl
      trap(3) = trap(3)-h2/12*(fpr-fpl)/fact  
      trapl=trap(3) +h4/720*(fp3r-fp3l)/fact3
      corr(3) = intf*fact-trap(3)
      corrl(3) = intf*fact-trapl
c     u3/rt5
      fpl= u2l*(3-2*u2l)/rt7l
      fp3l= -3*(-2+39*u2l-56*u4l+8*u6l)/rt11l
      Fl = -(2+3*u2l)/rt3l/3
      fpr= u2r*(3-2*u2r)/rt7r
      fp3r= -3*(-2+39*u2r-56*u4r+8*u6r)/rt11r
      Fr = -(2+3*u2r)/rt3r/3
      intf=Fr-Fl
      trap(4) = trap(4)-h2/12*(fpr-fpl)/fact  
      trapl=trap(4) +h4/720*(fp3r-fp3l)/fact3
      corr(4) = intf*fact-trap(4)
      corrl(4) = intf*fact-trapl
c     u5/rt5
      fpl= 5*u4l/rt7l
      fp3l= 15*u2l*(4-13*u2l+4*u4l)/rt11l
      Fl = (8+12*u2l+3*u4l)/(3*rt3l)
      fpr= 5*u4r/rt7r
      fp3r= 15*u2r*(4-13*u2r+4*u4r)/rt11r
      Fr = (8+12*u2r+3*u4r)/(3*rt3r)
      intf=Fr-Fl
      trap(5) = trap(5)-h2/12*(fpr-fpl)/fact  
c     +  +h4/720*(fp3r-fp3l)/fact3
      trapl=trap(5) +h4/720*(fp3r-fp3l)/fact3
      corr(5) = intf*fact-trap(5)
      corrl(5) = intf*fact-trapl
c     u7/rt5
      fpl=u6l*(7+2*u2l)/rt7l 
      fpr=u6r*(7+2*u2r)/rt7r 
      fp3l= -105*u4l*(u2l-2)/rt11l
      fp3r= -105*u4r*(u2r-2)/rt11r
      Fl = (u6l-6*u4l-24*u2l-16)/rt3l/3
      Fr = (u6r-6*u4r-24*u2r-16)/rt3r/3
      intf=Fr-Fl
      trap(16) = trap(16)-h2/12*(fpr-fpl)/fact  
c     +  +h4/720*(fp3r-fp3l)/fact3
      trapl=trap(16) +h4/720*(fp3r-fp3l)/fact3
      corr(16) = intf*fact-trap(16)
      corrl(16) = intf*fact-trapl
c     u/rt7                            !NEW
      fpl= (1-6*u2l)/rt9l 
      fpr= (1-6*u2r)/rt9r 
      Fl = -1/(5*rt5l)
      Fr = -1/(5*rt5r)
      intf=Fr-Fl
      trap(6) = trap(6)-h2/12*(fpr-fpl)/fact  
      corr(6) = intf*fact-trap(6)
c     u3/rt7                           !NEW
      fpl= u2l*(3-4*u2l)/rt9l 
      fpr= u2r*(3-4*u2r)/rt9r 
      Fl = -(2+5*u2l)/(15*rt5l)
      Fr = -(2+5*u2r)/(15*rt5r)
      intf=Fr-Fl
      trap(7) = trap(7)-h2/12*(fpr-fpl)/fact  
      corr(7) = intf*fact-trap(7)
c     u5/rt7
      fpl= u4l*(5-2*u2l)/rt9l
      fp3l= 3*u2l*(20-115*u2l+88*u4l-8*u6l)/rt13l
      Fl = -(8+20*u2l+15*u4l)/(15*rt5l)
      fpr= u4r*(5-2*u2r)/rt9r
      fp3r= 3*u2r*(20-115*u2r+88*u4r-8*u6r)/rt13r
      Fr = -(8+20*u2r+15*u4r)/(15*rt5r)
      intf=Fr-Fl
      trap(8) = trap(8)-h2/12*(fpr-fpl)/fact  
c     +  +h4/720*(fp3r-fp3l)/fact3
      trapl=trap(8)+h4/720*(fp3r-fp3l)/fact3
      corr(8) = intf*fact-trap(8)
      corrl(8) = intf*fact-trapl
c     u7/rt7
      fpl= 7*u6l/rt9l 
      fp3l= 21*u4l*(10-19*u2l+4*u4l)/rt13l 
      Fl = (16+40*u2l+30*u4l+5*u6l)/(5*rt5l)
      fpr= 7*u6r/rt9r 
      fp3r= 21*u4r*(10-19*u2r+4*u4r)/rt13r 
      Fr = (16+40*u2r+30*u4r+5*u6r)/(5*rt5r)
      intf=Fr-Fl
      trap(9) = trap(9)-h2/12*(fpr-fpl)/fact  
      trapl=trap(9)+h4/720*(fp3r-fp3l)/fact3
      corr(9) = intf*fact-trap(9)
      corrl(9) = intf*fact-trapl
c     u/rt9                        !NEW
      fpl= (1-8*u2l)/rt11l 
      fpr= (1-8*u2r)/rt11r 
      Fl = -1/(7*rt7l)
      Fr = -1/(7*rt7r)
      intf=Fr-Fl
      trap(10) = trap(10)-h2/12*(fpr-fpl)/fact  
      corr(10) = intf*fact-trap(10)
c     u3/rt9                       !NEW
      fpl= u2l*(3-6*u2l)/rt11l 
      fpr= u2r*(3-6*u2r)/rt11r 
      Fl = -(2+7*u2l)/(35*rt7l)
      Fr = -(2+7*u2r)/(35*rt7r)
      intf=Fr-Fl
      trap(11) = trap(11)-h2/12*(fpr-fpl)/fact  
      corr(11) = intf*fact-trap(11)
c     u5/rt9                       !NEW
      fpl= u4l*(5-4*u2l)/rt11l 
      fpr= u4r*(5-4*u2r)/rt11r 
      Fl = -(8+28*u2l+35*u4l)/(105*rt7l)
      Fr = -(8+28*u2r+35*u4r)/(105*rt7r)
      intf=Fr-Fl
      trap(12) = trap(12)-h2/12*(fpr-fpl)/fact  
      corr(12) = intf*fact-trap(12)
c     u7/rt9                       !NEW
      fpl= u6l*(7-2*u2l)/rt11l 
      fpr= u6r*(7-2*u2r)/rt11r 
      Fl = -(16+56*u2l+70*u4l+35*u6l)/(35*rt7l)
      Fr = -(16+56*u2r+70*u4r+35*u6r)/(35*rt7r)
      intf=Fr-Fl
      trap(13) = trap(13)-h2/12*(fpr-fpl)/fact  
      corr(13) = intf*fact-trap(13)
c     u9/rt9                       !NEW
      fpl= 9*u8l/rt11l
      fpr= 9*u8r/rt11r
      Fl = (128+448*u2l+560*u4l+280*u6l+35*u8l)/(35*rt7l)
      Fr = (128+448*u2r+560*u4r+280*u6r+35*u8r)/(35*rt7r)
      intf=Fr-Fl
      trap(17) = trap(17)-h2/12*(fpr-fpl)/fact  
      corr(17) = intf*fact-trap(17)
c     u5/rt11                       !NEW
      fpl= u4l*(5-6*u2l)/rt13l 
      fpr= u4r*(5-6*u2r)/rt13r 
      Fl = -(8+36*u2l+63*u4l)/(315*rt9l)
      Fr = -(8+36*u2r+63*u4r)/(315*rt9r)
      intf=Fr-Fl
      trap(14) = trap(14)-h2/12*(fpr-fpl)/fact  
      corr(14) = intf*fact-trap(14)
c     u7/rt11                       !NEW
      fpl= u6l*(7-4*u2l)/rt13l 
      fpr= u6r*(7-4*u2r)/rt13r 
      Fl = -(16+72*u2l+126*u4l+105*u6l)/(315*rt9l)
      Fr = -(16+72*u2r+126*u4r+105*u6r)/(315*rt9r)
      intf=Fr-Fl
      trap(15) = trap(15)-h2/12*(fpr-fpl)/fact  
      corr(15) = intf*fact-trap(15)
c     u9/rt11                       !NEW
      fpl= u8l*(9-2*u2l)/rt13l
      fpr= u8r*(9-2*u2r)/rt13r
      Fl = -(128+576*u2l+1008*u4l+840*u6l+315*u8l)/(315*rt9l)
      Fr = -(128+576*u2r+1008*u4r+840*u6r+315*u8r)/(315*rt9r)
      intf=Fr-Fl
      trap(18) = trap(18)-h2/12*(fpr-fpl)/fact  
      corr(18) = intf*fact-trap(18)

c we use CHANGE OF VARIABLES u=alf/fact, fact=d/cm
c integrals of interest are of the form 
c            fact*fact^k/d^l int u^k/sqrt(1+u2)^l du
c    errors above equal errors in fact*int u^k/root(1+usq)^l du
c    complete the change of variables by multiplying by
c    fkl=fact^k/d^l=d^(k-l)/c^k=(1/d^(l-k)*c^k)
      f13  = 1/(cm*d**2)  
      f33  = 1/(cm**3) 
      f53  = d**2/(cm**5) 
      f15  = 1/(cm*d**4) 
      f35  = 1/(cm**3*d**2) 
      f55  = 1/(cm**5) 
      f75  = d**2/(cm**7) 
      f17  = 1/(cm*d**6) 
      f37  = 1/(cm**3*d**4) 
      f57  = 1/(cm**5*d**2) 
      f77  = 1/(cm**7) 
      f19  = 1/(cm*d**8) 
      f39  = 1/(cm**3*d**6) 
      f59  = 1/(cm**5*d**4) 
      f79  = 1/(cm**7*d**2) 
      f99  = 1/(cm**9) 
      f511 = 1/(cm**5*d**6) 
      f711 = 1/(cm**7*d**4) 
      f911 = 1/(cm**9*d**2) 
      corr(3) = corr(3)*f15     ! u/rt5      ! u/(1+u2)^5/2
      corr(4) = corr(4)*f35     ! u3/rt5
      corr(5) = corr(5)*f55      !u5/(1+u2)^5/2
      corr(16)= corr(16)*f75      !u5/(1+u2)^5/2
      corr(6) = corr(6)*f17      !u/rt7
      corr(7) = corr(7)*f37      !u3/rt7
      corr(8) = corr(8)*f57      !u5/(1+u2)^7/2
      corr(9) = corr(9)*f77      !u7/(1+u2)^7/2
      corr(10)= corr(10)*f19     !u/rt9
      corr(11)= corr(11)*f39     !u3/rt9
      corr(12)= corr(12)*f59     !u5/rt9
      corr(13)= corr(13)*f79     !u7/rt9
      corr(17)= corr(17)*f99     !u7/rt9
      corr(14)= corr(14)*f511    !u5/rt11
      corr(15)= corr(15)*f711    !u7/rt11
      corr(18)= corr(18)*f911    !u9/rt11

      corrl(1) = corrl(1)*f13     ! u/rt3      ! u/(1+u2)^3/2
      corrl(2) = corrl(2)*f33     ! u3/rt3      !u3/(1+u2)^3/2
      corr(19) = corr(19)*f53     ! u5/rt3
      corrl(3) = corrl(3)*f15     ! u/rt5      ! u/(1+u2)^5/2
      corrl(4) = corrl(4)*f35     ! u3/rt5
      corrl(5) = corrl(5)*f55      !u5/(1+u2)^5/2
      corrl(16)= corrl(16)*f75      !u5/(1+u2)^5/2
      corrl(8) = corrl(8)*f57      !u5/(1+u2)^7/2
      corrl(9) = corrl(9)*f77      !u7/(1+u2)^7/2

      return
      end

