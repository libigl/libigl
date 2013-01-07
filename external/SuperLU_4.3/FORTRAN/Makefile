include ../make.inc

#######################################################################
#  This makefile creates the Fortran example interface to use the
#  C routines in SuperLU.
#######################################################################

HEADER   = ../SRC
LIBS	= $(SUPERLULIB) $(BLASLIB) -lm

# double real
F77EXM	= f77_main.o hbcode1.o c_fortran_dgssv.o

# double complex
ZF77EXM	= z_f77_main.o zhbcode1.o c_fortran_zgssv.o

all:	f77exm zf77exm

f77exm: $(F77EXM) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(F77EXM) $(LIBS) -o $@

zf77exm: $(ZF77EXM) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(ZF77EXM) $(LIBS) -o $@

c_fortran_zgssv.o: c_fortran_zgssv.c
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o f77exm zf77exm

