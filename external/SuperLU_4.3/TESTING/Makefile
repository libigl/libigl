include ../make.inc

#######################################################################
#  This makefile creates the test programs for the linear equation
#  routines in SuperLU.  The test files are grouped as follows:
#
#       ALINTST -- Auxiliary test routines
#       SLINTST -- Single precision real test routines
#       DLINTST -- Double precision real test routines
#       CLINTST -- Double precision complex test routines
#       ZLINTST -- Double precision complex test routines
#
#  Test programs can be generated for all or some of the four different
#  precisions.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single double
#  Alternatively, the command
#       make
#  without any arguments creates all four test programs.
#  The executable files are called
#       stest
#       dtest
#       ctest
#       ztest
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################

HEADER  = ../SRC

ALINTST = sp_ienv.o

SLINTST = sdrive.o sp_sconvert.o sgst01.o sgst02.o sgst04.o sgst07.o 

DLINTST = ddrive.o sp_dconvert.o dgst01.o dgst02.o dgst04.o dgst07.o 

CLINTST = cdrive.o sp_cconvert.o cgst01.o cgst02.o cgst04.o cgst07.o

ZLINTST = zdrive.o sp_zconvert.o zgst01.o zgst02.o zgst04.o zgst07.o

all: testmat single double complex complex16

testmat:
	(cd MATGEN; $(MAKE))

single: ./stest stest.out

./stest: $(SLINTST) $(ALINTST) $(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(SLINTST) $(ALINTST) \
        $(TMGLIB) $(SUPERLULIB) $(BLASLIB) -lm -o $@

stest.out: stest stest.csh
	@echo Testing SINGLE PRECISION linear equation routines 
	csh stest.csh

double: ./dtest dtest.out

./dtest: $(DLINTST) $(ALINTST) $(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(DLINTST) $(ALINTST) \
        $(TMGLIB) $(SUPERLULIB) $(BLASLIB) -lm -o $@

dtest.out: dtest dtest.csh
	@echo Testing DOUBLE PRECISION linear equation routines 
	csh dtest.csh

complex: ./ctest ctest.out

./ctest: $(CLINTST) $(ALINTST) $(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(CLINTST) $(ALINTST) \
	$(TMGLIB) $(SUPERLULIB) $(BLASLIB) -lm -o $@

ctest.out: ctest ctest.csh
	@echo Testing SINGLE COMPLEX linear equation routines 
	csh ctest.csh

complex16: ./ztest ztest.out

./ztest: $(ZLINTST) $(ALINTST) $(SUPERLULIB) $(TMGLIB)
	$(LOADER) $(LOADOPTS) $(ZLINTST) $(ALINTST) \
	$(TMGLIB) $(SUPERLULIB) $(BLASLIB) -lm -o $@

ztest.out: ztest ztest.csh
	@echo Testing DOUBLE COMPLEX linear equation routines 
	csh ztest.csh

##################################
# Do not optimize this routine   #
##################################
dlamch.o: dlamch.c ; $(CC) -c $<

timer.o:  timer.c ; $(CC) -O -c $<

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

clean:	
	rm -f *.o *test *.out

