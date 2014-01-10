include makeconf

SOURCEDIR    =  src
TESTDIR      =  test

#-------------------- RULES --------------------------------------

# Always remake everything if make is called from this top-level directory

all:	clean depend
	$(MAKE) -C $(SOURCEDIR) lib
	$(MAKE) -C $(TESTDIR)

distclean: clean
	$(RM) -f lib/libjmesh.a

clean:
	$(MAKE) -C $(SOURCEDIR) clean
	$(MAKE) -C $(TESTDIR) clean

depend:
	$(MAKE) -C $(SOURCEDIR) depend

backup :
	$(MAKE) -C $(SOURCEDIR) backup
