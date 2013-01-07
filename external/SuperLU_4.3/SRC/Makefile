# makefile for sparse supernodal LU, implemented in ANSI C
include ../make.inc

#######################################################################
#  This is the makefile to create a library for SuperLU.
#  The files are organized as follows:
#
#       ALLAUX -- Auxiliary routines called from all precisions of SuperLU
#	LAAUX  -- LAPACK auxiliary routines called from all precisions
#	SLASRC -- LAPACK single precision real routines
#	DLASRC -- LAPACK double precision real routines
#	CLASRC -- LAPACK single precision complex routines
#	ZLASRC -- LAPACK double precision complex routines
#       SCLAUX -- LAPACK Auxiliary routines called from both real and complex
#       DZLAUX -- LAPACK Auxiliary routines called from both double precision
#                 and complex*16
#	SLUSRC -- Single precision real SuperLU routines
#       DLUSRC -- Double precision real SuperLU routines
#	CLUSRC -- Single precision complex SuperLU routines
#       ZLUSRC -- Double precision complex SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single double
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       superlu.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################

### LAPACK 
LAAUX 	= lsame.o xerbla.o
SLASRC 	= slacon.o 
DLASRC	= dlacon.o
CLASRC	= clacon.o scsum1.o icmax1.o
ZLASRC	= zlacon.o dzsum1.o izmax1.o
SCLAUX 	= slamch.o
DZLAUX 	= dlamch.o

### SuperLU 
ALLAUX 	= superlu_timer.o util.o memory.o get_perm_c.o mmd.o \
	  sp_coletree.o sp_preorder.o sp_ienv.o relax_snode.o \
	  heap_relax_snode.o colamd.o \
	  ilu_relax_snode.o ilu_heap_relax_snode.o mark_relax.o \
	  mc64ad.o qselect.o

SLUSRC = \
	sgssv.o sgssvx.o \
	ssp_blas2.o ssp_blas3.o sgscon.o  \
	slangs.o sgsequ.o slaqgs.o spivotgrowth.o \
	sgsrfs.o sgstrf.o sgstrs.o scopy_to_ucol.o \
	ssnode_dfs.o ssnode_bmod.o \
	spanel_dfs.o spanel_bmod.o \
	sreadhb.o sreadrb.o sreadtriple.o \
	scolumn_dfs.o scolumn_bmod.o spivotL.o spruneL.o \
	smemory.o sutil.o smyblas2.o \
	sgsisx.o sgsitrf.o sldperm.o \
	ilu_sdrop_row.o ilu_ssnode_dfs.o \
	ilu_scolumn_dfs.o ilu_spanel_dfs.o ilu_scopy_to_ucol.o \
	ilu_spivotL.o sdiagonal.o

DLUSRC = \
	dgssv.o dgssvx.o \
	dsp_blas2.o dsp_blas3.o dgscon.o \
	dlangs.o dgsequ.o dlaqgs.o dpivotgrowth.o  \
	dgsrfs.o dgstrf.o dgstrs.o dcopy_to_ucol.o \
	dsnode_dfs.o dsnode_bmod.o dpanel_dfs.o dpanel_bmod.o \
	dreadhb.o dreadrb.o dreadtriple.o \
	dcolumn_dfs.o dcolumn_bmod.o dpivotL.o dpruneL.o \
	dmemory.o dutil.o dmyblas2.o \
	dgsisx.o dgsitrf.o dldperm.o \
	ilu_ddrop_row.o ilu_dsnode_dfs.o \
	ilu_dcolumn_dfs.o ilu_dpanel_dfs.o ilu_dcopy_to_ucol.o \
	ilu_dpivotL.o ddiagonal.o 
        ## dgstrsL.o dgstrsU.o

CLUSRC = \
	scomplex.o cgssv.o cgssvx.o csp_blas2.o csp_blas3.o cgscon.o \
	clangs.o cgsequ.o claqgs.o cpivotgrowth.o  \
	cgsrfs.o cgstrf.o cgstrs.o ccopy_to_ucol.o \
	csnode_dfs.o csnode_bmod.o \
	cpanel_dfs.o cpanel_bmod.o \
	creadhb.o creadrb.o creadtriple.o \
	ccolumn_dfs.o ccolumn_bmod.o cpivotL.o cpruneL.o \
	cmemory.o cutil.o cmyblas2.o \
	cgsisx.o cgsitrf.o cldperm.o \
	ilu_cdrop_row.o ilu_csnode_dfs.o \
	ilu_ccolumn_dfs.o ilu_cpanel_dfs.o ilu_ccopy_to_ucol.o \
	ilu_cpivotL.o cdiagonal.o

ZLUSRC = \
	dcomplex.o zgssv.o zgssvx.o zsp_blas2.o zsp_blas3.o zgscon.o \
	zlangs.o zgsequ.o zlaqgs.o zpivotgrowth.o  \
	zgsrfs.o zgstrf.o zgstrs.o zcopy_to_ucol.o \
	zsnode_dfs.o zsnode_bmod.o \
	zpanel_dfs.o zpanel_bmod.o \
	zreadhb.o zreadrb.o zreadtriple.o \
	zcolumn_dfs.o zcolumn_bmod.o zpivotL.o zpruneL.o \
	zmemory.o zutil.o zmyblas2.o \
	zgsisx.o zgsitrf.o zldperm.o \
	ilu_zdrop_row.o ilu_zsnode_dfs.o \
	ilu_zcolumn_dfs.o ilu_zpanel_dfs.o ilu_zcopy_to_ucol.o \
	ilu_zpivotL.o zdiagonal.o

all:    single double complex complex16

single: $(SLUSRC) $(ALLAUX) $(LAAUX) $(SLASRC) $(SCLAUX)
	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) \
		$(SLUSRC) $(ALLAUX) $(LAAUX) $(SLASRC) $(SCLAUX)
	$(RANLIB) $(SUPERLULIB)

double: $(DLUSRC) $(ALLAUX) $(LAAUX) $(DLASRC) $(DZLAUX)
	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) \
                $(DLUSRC) $(ALLAUX) $(LAAUX) $(DLASRC) $(DZLAUX)
	$(RANLIB) $(SUPERLULIB)

complex: $(CLUSRC) $(ALLAUX) $(LAAUX) $(CLASRC) $(SCLAUX)
	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) \
		$(CLUSRC) $(ALLAUX) $(LAAUX) $(CLASRC) $(SCLAUX)
	$(RANLIB) $(SUPERLULIB)

complex16: $(ZLUSRC) $(ALLAUX) $(LAAUX) $(ZLASRC) $(DZLAUX)
	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) \
		$(ZLUSRC) $(ALLAUX) $(LAAUX) $(ZLASRC) $(DZLAUX)
	$(RANLIB) $(SUPERLULIB)


##################################
# Do not optimize these routines #
##################################
slamch.o: slamch.c ; $(CC) -c $(NOOPTS) $(CDEFS) $<
dlamch.o: dlamch.c ; $(CC) -c $(NOOPTS) $(CDEFS) $<
superlu_timer.o:  superlu_timer.c ; $(CC) -c $(NOOPTS) $<
##################################

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $<

clean:	
	rm -f *.o $(SUPERLULIB)
