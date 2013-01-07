include ../../make.inc

#######################################################################
#  This is the makefile to create a library of the test matrix
#  generators used in LAPACK.  The files are organized as follows:
#
#     SCATGEN  -- Auxiliary routines called from both REAL and COMPLEX
#     DZATGEN  -- Auxiliary routines called from both DOUBLE PRECISION
#                 and COMPLEX*16
#     SMATGEN  -- Single precision real matrix generation routines
#     CMATGEN  -- Single precision complex matrix generation routines
#     DMATGEN  -- Double precision real matrix generation routines
#     ZMATGEN  -- Double precision complex matrix generation routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       libtmglib.a
#  and is created at the LAPACK directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################
ALLAUX  = lsamen.o

SCATGEN = slatm1.o slaran.o slarnd.o slaruv.o slabad.o slarnv.o
SLASRC  = slatb4.o slaset.o slartg.o
SMATGEN = slatms.o slatme.o slatmr.o \
	  slagge.o slagsy.o slarge.o slaror.o slarot.o slatm2.o slatm3.o
SINTRINSIC  = r_lg10.o r_sign.o pow_dd.o

DZATGEN = dlatm1.o dlaran.o dlarnd.o dlaruv.o dlabad.o dlarnv.o
DLASRC  = dlatb4.o dlaset.o dlartg.o
DMATGEN = dlatms.o dlatme.o dlatmr.o \
	  dlagge.o dlagsy.o dlarge.o dlaror.o dlarot.o dlatm2.o dlatm3.o
DINTRINSIC  = d_lg10.o d_sign.o pow_dd.o

CLASRC  = clatb4.o claset.o clartg.o clarnv.o clacgv.o csymv.o
CMATGEN = clatms.o clatme.o clatmr.o \
	  clagge.o clagsy.o clarge.o claror.o clarot.o clatm2.o clatm3.o \
	  claghe.o clarnd.o cdotc.o

ZLASRC  = zlatb4.o zlaset.o zlartg.o zlarnv.o zlacgv.o zsymv.o
ZMATGEN = zlatms.o zlatme.o zlatmr.o \
	  zlagge.o zlagsy.o zlarge.o zlaror.o zlarot.o zlatm2.o zlatm3.o \
	  zlaghe.o zlarnd.o zdotc.o

all: single double complex complex16

single: $(SMATGEN) $(SCATGEN) $(SLASRC) $(SINTRINSIC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(TMGLIB) $(SMATGEN) $(SCATGEN) \
	        $(SLASRC) $(SINTRINSIC) $(ALLAUX)
	$(RANLIB) ../$(TMGLIB)

double: $(DMATGEN) $(DZATGEN) $(DLASRC) $(DINTRINSIC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(TMGLIB) $(DMATGEN) $(DZATGEN) \
	        $(DLASRC) $(DINTRINSIC) $(ALLAUX)
	$(RANLIB) ../$(TMGLIB)

complex: $(CMATGEN) $(SCATGEN) $(CLASRC) $(SINTRINSIC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(TMGLIB) $(CMATGEN) $(SCATGEN) \
	        $(CLASRC) $(SINTRINSIC) $(ALLAUX)
	$(RANLIB) ../$(TMGLIB)

complex16: $(ZMATGEN) $(DZATGEN) $(ZLASRC) $(DINSTRINSIC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(TMGLIB) $(ZMATGEN) $(DZATGEN) \
	        $(ZLASRC) $(DINTRINSIC) $(ALLAUX)
	$(RANLIB) ../$(TMGLIB)

clean:
	rm -f *.o ../$(TMGLIB)

.c.o: ; $(CC) $(CFLAGS) $(CDEFS) -c $<
