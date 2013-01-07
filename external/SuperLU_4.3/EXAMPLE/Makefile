include ../make.inc

#######################################################################
#  This makefile creates the example programs for the linear equation
#  routines in SuperLU.  The files are grouped as follows:
#
#       SLINEXM -- Single precision real example routines
#       DLINEXM -- Double precision real example routines
#       CLINEXM -- Double precision complex example routines
#       ZLINEXM -- Double precision complex example routines
#
#  Example programs can be generated for all or some of the four different
#  precisions.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single double
#  Alternatively, the command
#       make
#  without any arguments creates all four example programs.
#  The executable files are called
#       slinsol		slinsolx
#       dlinsol		dlinsolx
#       clinsol		clinsolx
#       zlinsol		zlinsolx
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################

HEADER   = ../SRC

SLINEXM		= slinsol.o
SLINEXM1	= slinsol1.o
SLINXEXM 	= slinsolx.o
SLINXEXM1 	= slinsolx1.o
SLINXEXM2	= slinsolx2.o
SITSOL	     	= sitersol.o sfgmr.o
SITSOL1	     	= sitersol1.o sfgmr.o

DLINEXM		= dlinsol.o
DLINEXM1	= dlinsol1.o
DLINXEXM 	= dlinsolx.o
DLINXEXM1 	= dlinsolx1.o
DLINXEXM2 	= dlinsolx2.o
SUPERLUEXM 	= superlu.o sp_ienv.o
DITSOL	     	= ditersol.o dfgmr.o
DITSOL1	     	= ditersol1.o dfgmr.o

CLINEXM 	= clinsol.o
CLINEXM1 	= clinsol1.o
CLINXEXM 	= clinsolx.o
CLINXEXM1 	= clinsolx1.o
CLINXEXM2 	= clinsolx2.o
CITSOL	     	= citersol.o cfgmr.o
CITSOL1	     	= citersol1.o cfgmr.o

ZLINEXM 	= zlinsol.o
ZLINEXM1 	= zlinsol1.o
ZLINXEXM 	= zlinsolx.o
ZLINXEXM1 	= zlinsolx1.o
ZLINXEXM2 	= zlinsolx2.o
ZITSOL	     	= zitersol.o zfgmr.o
ZITSOL1	     	= zitersol1.o zfgmr.o


all:    single double complex complex16

single:    slinsol slinsol1 slinsolx slinsolx1 slinsolx2 sitersol sitersol1
double:    dlinsol dlinsol1 dlinsolx dlinsolx1 dlinsolx2 superlu \
		ditersol ditersol1
complex:   clinsol clinsol1 clinsolx clinsolx1 clinsolx2 citersol citersol1
complex16: zlinsol zlinsol1 zlinsolx zlinsolx1 zlinsolx2 zitersol zitersol1

slinsol: $(SLINEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINEXM) $(LIBS) -lm -o $@

slinsol1: $(SLINEXM1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINEXM1) $(LIBS) -lm -o $@

slinsolx: $(SLINXEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINXEXM) $(LIBS) -lm -o $@

slinsolx1: $(SLINXEXM1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINXEXM1) $(LIBS) -lm -o $@

slinsolx2: $(SLINXEXM2) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINXEXM2) $(LIBS) -lm -o $@

sitersol: $(SITSOL) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SITSOL) $(LIBS) -lm -o $@

sitersol1: $(SITSOL1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SITSOL1) $(LIBS) -lm -o $@

dlinsol: $(DLINEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINEXM) $(LIBS) -lm -o $@

dlinsol1: $(DLINEXM1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINEXM1) $(LIBS) -lm -o $@

dlinsolx: $(DLINXEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINXEXM) $(LIBS) -lm -o $@

dlinsolx1: $(DLINXEXM1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINXEXM1) $(LIBS) -lm -o $@

dlinsolx2: $(DLINXEXM2) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINXEXM2) $(LIBS) -lm -o $@

superlu: $(SUPERLUEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SUPERLUEXM) $(LIBS) -lm -o $@

ditersol: $(DITSOL) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DITSOL) $(LIBS) -lm -o $@

ditersol1: $(DITSOL1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DITSOL1) $(LIBS) -lm -o $@

clinsol: $(CLINEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINEXM) $(LIBS) -lm -o $@

clinsol1: $(CLINEXM1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINEXM1) $(LIBS) -lm -o $@

clinsolx: $(CLINXEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINXEXM) $(LIBS) -lm -o $@

clinsolx1: $(CLINXEXM1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINXEXM1) $(LIBS) -lm -o $@

clinsolx2: $(CLINXEXM2) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINXEXM2) $(LIBS) -lm -o $@

citersol: $(CITSOL) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CITSOL) $(LIBS) -lm -o $@

citersol1: $(CITSOL1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CITSOL1) $(LIBS) -lm -o $@

zlinsol: $(ZLINEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINEXM) $(LIBS) -lm -o $@

zlinsol1: $(ZLINEXM1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINEXM1) $(LIBS) -lm -o $@

zlinsolx: $(ZLINXEXM) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINXEXM) $(LIBS) -lm -o $@

zlinsolx1: $(ZLINXEXM1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINXEXM1) $(LIBS) -lm -o $@

zlinsolx2: $(ZLINXEXM2) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINXEXM2) $(LIBS) -lm -o $@

zitersol: $(ZITSOL) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZITSOL) $(LIBS) -lm -o $@

zitersol1: $(ZITSOL1) $(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZITSOL1) $(LIBS) -lm -o $@

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o *linsol *linsol1 *linsolx *linsolx1 *linsolx2 \
	    superlu *itersol



