UNAME = $(shell uname)

FFLAGS   = -g -Og -fbacktrace
# on papaya
#  F77      = pgf77
# on guava, and on shell
F77      = gfortran
CC       = gcc
CFLAGS   = -DLINUX
CPP      = cpp
CPPFLAGS = -DLINUX $(IPATH)

#             # Pattern Rules
%.o: %.F params.h derivs1.h vars.h flags.h coeff.h 
	$(F77) $(CPPFLAGS) $(FFLAGS) -c $*.F

#  Pattern Rules for special files requiring flags.h
#correctqnearsheetnewv2q.o: correctqnearsheetnewv2q.F flags.h 
#	$(F77) $(CPPFLAGS) $(FFLAGS) -c correctqnearsheetnewv2q.F

OBJTSTNEW = drivergrid.o \
	Intqs.o Intms.o velodSonO.o veloforceSonO.o compendcontrib.o \
	ellip.o trigderiv.o dfft.o ftrans.o compxyvol.o compderivxy1.o \
	x2cl2circ.o correctionend.o correctqnearsheet.o correctionqatx.o \
	compallcorr.o compallcorrend.o \
	compncorrBq.o compallcoeffuys.o compsmallr.o compaxis.o velodSonP.o \
	correctmnearsheet.o intbm.o correctionmatxnew.o compncorrBm.o \
	compsmallrm.o veloforceSonP.o


drivergrid.out : $(OBJTSTNEW) 
	$(F77) -o drivergrid.out $(OBJTSTNEW) -L/usr/lib -llapack -lblas

testaxis.out: testaxis.o compallcorrend.o correctionend.o compaxis.o
	$(F77) -o testaxis.out testaxis.o compallcorrend.o correctionend.o compaxis.o

clean : 
	rm -f *.o
