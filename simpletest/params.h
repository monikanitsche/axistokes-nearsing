c     *** Parameters ***
      implicit none
      integer nmax,nmax2,nb
      integer ntop,nwall
      MODE pi,pi2,pi4,pi8,logterm,ca,nu,flevelx,log2,
     +  ed1,ed2,ed3,ec1,ec2,ec3
      MODE logtermw1,logtermw2,cosbetsq,csbet,radopening  !capnum
      parameter(pi= 3.14159265358979323846264338327950280q0,
     +   pi2=2*pi,pi4=4*pi,pi8=8*pi)
      parameter(log2=log(2.q0))
      parameter( ed1=1.q0/4, ed2=3.q0/32, ed3=15.q0/256,
     +  ec1=(4*log2-1)/4, ec2=(24*log2-13)/64, ec3=3*(5*log2-3)/64)
#ifdef QUAD
      parameter(ca=1e-20)
      parameter(flevelx=1e-20)
#endif
#ifdef DBLE
      parameter(ca=1e-14)
      parameter(flevelx=1e-13)
#endif
      parameter(nmax=60100,nmax2=2*nmax+1)
      parameter(nb=35) !number of basis functions

      common/constants/logtermw1,logtermw2,cosbetsq,csbet,logterm,
     +   radopening
      common/plateint/ntop,nwall
