include ../make.inc
HEADER = ../SRC

#######################################################################
#  This is the makefile to create a library for C-BLAS.
#  The files are organized as follows:
#
#       SBLAS1 -- Single precision real BLAS routines
#       CBLAS1 -- Single precision complex BLAS routines
#       DBLAS1 -- Double precision real BLAS routines
#       ZBLAS1 -- Double precision complex BLAS routines
#
#       CB1AUX -- Real BLAS routines called by complex routines
#       ZB1AUX -- D.P. real BLAS routines called by d.p. complex
#                 routines
#
#      ALLBLAS -- Auxiliary routines for Level 2 and 3 BLAS
#
#       SBLAS2 -- Single precision real BLAS2 routines
#       CBLAS2 -- Single precision complex BLAS2 routines
#       DBLAS2 -- Double precision real BLAS2 routines
#       ZBLAS2 -- Double precision complex BLAS2 routines
#
#       SBLAS3 -- Single precision real BLAS3 routines
#       CBLAS3 -- Single precision complex BLAS3 routines
#       DBLAS3 -- Double precision real BLAS3 routines
#       ZBLAS3 -- Double precision complex BLAS3 routines
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
#       blas.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################

SBLAS1 = isamax.o sasum.o saxpy.o scopy.o sdot.o snrm2.o \
	 srot.o sscal.o sswap.o
SBLAS2 = sgemv.o ssymv.o strsv.o sger.o ssyr2.o

DBLAS1 = idamax.o dasum.o daxpy.o dcopy.o ddot.o dnrm2.o \
	 drot.o dscal.o dswap.o
DBLAS2 = dgemv.o dsymv.o dtrsv.o dger.o dsyr2.o

CBLAS1 = icamax.o scasum.o caxpy.o ccopy.o scnrm2.o \
	 cscal.o cswap.o
CBLAS2 = cgemv.o chemv.o ctrsv.o cgerc.o cher2.o

ZBLAS1 = izamax.o dzasum.o zaxpy.o zcopy.o dznrm2.o \
	 zscal.o dcabs1.o zswap.o
ZBLAS2 = zgemv.o zhemv.o ztrsv.o zgerc.o zher2.o


all: single double complex complex16

single: $(SBLAS1) $(SBLAS2) $(SBLAS3)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(SBLAS1) $(ALLBLAS) $(SBLAS2) $(SBLAS3)
	$(RANLIB) $(BLASLIB)

double: $(DBLAS1) $(DBLAS2) $(DBLAS3)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(DBLAS1) $(ALLBLAS) $(DBLAS2) $(DBLAS3)
	$(RANLIB) $(BLASLIB)

complex: $(CBLAS1) $(CBLAS2) $(CBLAS3)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(CBLAS1) $(ALLBLAS) $(CBLAS2) $(CBLAS3)
	$(RANLIB) $(BLASLIB)

complex16: $(ZBLAS1) $(ZBLAS2) $(ZBLAS3)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(ZBLAS1) $(ALLBLAS) $(ZBLAS2) $(ZBLAS3)
	$(RANLIB) $(BLASLIB)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

clean:	
	rm -f *.o ../lib/libblas.a
