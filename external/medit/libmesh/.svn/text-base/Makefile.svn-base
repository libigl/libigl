# This is the Makefile for the libMesh library and helper
# applications.  This file is specific to the project.
# See the file Make.common for path configurations.
# Make.common is what should be included by applications
# using the library.

# Build the library and all utility applications by default:

default: all

# include the library options determined by configure
# (now included indirectly in Makefile.const - don't want to
# accidentally double-up on += variable modifications
# include Make.common

###############################################################################
include Makefile.const


###############################################################################
# Target:
#
target := $(mesh_library)

library: $(mesh_library)

bin/$(hosttype)_$(METHOD)/libmesh-config: contrib/bin/libmesh-config
	mkdir -p bin/$(hosttype)_$(METHOD)
	ln -sf ../../contrib/bin/libmesh-config bin/$(hosttype)_$(METHOD)/
	
all:: $(target) $(appbinfiles) bin/$(hosttype)_$(METHOD)/libmesh-config

###############################################################################

#
# 
#
.PHONY: all clean clobber default distclean library target \
	doc upload doc_upload log cvsweb TODO svnexpand examples

#
# static library
#
ifeq ($(findstring darwin,$(hostos)),darwin)
$(mesh_library_dir)/libmesh$(static_libext): $(objects)
	@$(shell mkdir -p $(mesh_library_dir))
	@echo "Linking "$@
	@libtool -static -o $(mesh_library) $(objects)
	@$(MAKE) -C contrib

else
$(mesh_library_dir)/libmesh$(static_libext): $(objects)
	@$(shell mkdir -p $(mesh_library_dir))
	@rm -f $(mesh_library)
	@echo "Linking "$@
	@$(AR) rv $(mesh_library) $(objects)
	@$(MAKE) -C contrib
endif
# shared library
#
$(mesh_library_dir)/libmesh$(shared_libext): $(objects)
	@$(shell mkdir -p $(mesh_library_dir))
	@echo "Linking "$@
	@$(libmesh_CXX) $(libmesh_CXXSHAREDFLAG) -o $(mesh_library) $(objects) $(libmesh_LDFLAGS)
	@$(MAKE) -C contrib

#
# Build just object files
#
obj: $(objects)

#
# Build the examples.  Can be run in parallel with "make -jN examples"
#
examples: $(mesh_library)
	@$(MAKE) -C examples

# Only link the examples.  Deprecated, use "make examples" instead.
#
# link_examples: $(mesh_library)
# 	@$(MAKE) -C examples link

#
# Run the examples
#
run_examples: $(mesh_library)
	@$(MAKE) -C examples run

#
# Clean the examples
#
clean_examples:
	@$(MAKE) -C examples clean

#
# Test the header files to make sure they all compile stand-alone
# Also a cheat here: removing an ex-configure-output so it isn't
# detected as an invalid header for people without fresh checkouts
#
test_headers:
	@rm -f include/base/libmesh_contrib_config.h
	@cd contrib/bin && ./test_headers.sh

#	
# Remove object files for the current mode
#
clean:
	@test ! -d contrib || $(MAKE) -C contrib clean
	@test ! -d examples || $(MAKE) -C examples clean
	@rm -f *~ include/*~ include/*/*~ src/*/*~ src/*/*.$(obj-suffix) doc/html/*~
	@rm -f .depend

#
# Make clean, remove all binaries and generated files.  Leaves libraries in-tact
#
clobber:
	@$(MAKE) clean
	@test ! -d contrib || $(MAKE) -C contrib clobber
	@test ! -d examples || $(MAKE) -C examples clobber
	@rm -rf config.status $(targ_dir) $(appbinfiles)

#
# Make clobber, remove documentation, removes all libraries & object files for all modes
# Should restore to a pristine state, except for files you
# have added
distclean:
	@rm -f src/*/*.$(hosttype).*.o
	@$(MAKE) clobber
	@test ! -d contrib || $(MAKE) -C contrib distclean
	@test ! -d examples || $(MAKE) -C examples distclean
	@rm -rf doc/man/man3
	@rm -rf doc/html/doxygen/*.html # split these up, otherwise command line gets too long
	@rm -rf doc/html/doxygen/*.php
	@rm -rf doc/html/doxygen/*.png
	@rm -rf doc/html/doxygen/*.gif
	@rm -rf doc/html/doxygen/*.map
	@rm -rf doc/html/doxygen/*.md5
	@rm -rf doc/html/doxygen/*.dot
	@rm -rf doc/html/doxygen/formula.repository doc/html/doxygen/graph_legend.dot
	@rm -rf doc/latex/doxygen
	@rm -rf doc/latex/*/*.aux doc/latex/*/*~ doc/latex/*/*.log doc/latex/*/*.out
	@rm -rf src/*/*.o
	@rm -rf lib/*_opt lib/*_dbg lib/*_pro lib/*_devel lib/*_oprof
	@rm -rf .depend



#
# all documentation
#
doc: doxygen_doc examples_doc

#
# doxygen documentation
#
doxygen_doc:
	$(doxygen) ./doc/Doxyfile
	@rm -rf doc/html/doxygen/*.map
	@rm -rf doc/html/doxygen/*.md5
	@rm -rf doc/html/doxygen/*.dot
	@rm -rf doc/html/doxygen/formula.repository

#
# examples documentation
#
examples_doc:
	contrib/bin/create_example_docs.sh

#
# Upload the web page to sourceforge.  We need a way to specify usernames
# other than $USER when connecting to sourceforge servers... Please set the
# environment variable: $LIBMESH_SVN_USER if your sourceforge username
# is different than whatever $USER is.
#
ifeq (x$(LIBMESH_SVN_USER),x) #unset
  upload_name=$(USER),libmesh@
else
  upload_name=$(LIBMESH_SVN_USER),libmesh@
endif

upload:
	chmod -R g+w ./doc/html/* ./doc/latex/*/*
	rsync -rltzve ssh --exclude '.svn' ./doc/html/ $(upload_name)web.sourceforge.net:/home/groups/l/li/libmesh/htdocs
	rsync -rltzve ssh --exclude '.svn' ./doc/latex/howto $(upload_name)web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -rltzve ssh --exclude '.svn' ./doc/latex/xda_format $(upload_name)web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	chmod -R g-w ./doc/html/* ./doc/latex/*/*


#
# Test Uploading (rsync -n) the web page to sourceforge.  Try this
# if you want to see what *would be* uploaded with a real rsync.
#
upload_test:
	chmod -R g+w ./doc/html/* ./doc/latex/*/*
	rsync -nrltzve ssh --exclude '.svn' ./doc/html/ $(upload_name)libmesh.sourceforge.net:/home/groups/l/li/libmesh/htdocs
	rsync -nrltzve ssh --exclude '.svn' ./doc/latex/howto $(upload_name)libmesh.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -nrltzve ssh --exclude '.svn' ./doc/latex/xda_format $(upload_name)libmesh.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	chmod -R g-w ./doc/html/* ./doc/latex/*/*


#
# Build and upload the documentation to sourceforge
#
doc_upload:
	@$(MAKE) doc
	@$(MAKE) upload


log: $(loggedfiles)
	cvs log $(loggedfiles) > cvs_log

#
# Web-based CVS documentation
# Please note: SourceForge provides a nice web interface to the CVS
# logs for the library.  You probably shouldn't need this target for
# anything.
cvsweb:
	./contrib/bin/cvs2html -f -p -o doc/cvshtml/index.html -v -a -b -n 2 -C crono.html

svnexpand:
        svn propset svn:keywords "Date Author Revision HeadURL Id" $(srcfiles) $(headerfiles) m4/*.m4

#
# Standalone applications.  Anything in the ./src/apps directory that ends in .C
# can be compiled with this rule.  For example, if ./src/apps/foo.C contains a main()
# and is a standalone program, then make bin/foo will work.
#
bin/%$(bin-suffix) : src/apps/%.$(obj-suffix) $(mesh_library)
	@echo "Linking $@"
	@$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(patsubst %.C,%.$(obj-suffix),$<) -o $@ $(libmesh_LIBS) $(libmesh_DLFLAGS) $(libmesh_LDFLAGS)

src/apps/%.$(obj-suffix) : src/apps/%.C
	@echo "Building $@"
	@$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c $< -o $(patsubst %.C,%.$(obj-suffix),$<) 

#
# In the contrib/bin directory, we run the test_headers.sh shell
# script.  This is a make rule for those tests.
#
contrib/bin/%.o : contrib/bin/%.C
	$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c $< -o $@

#
# Make a TODO list
#
TODO:
	@egrep -i '// *todo' $(srcfiles) $(headerfiles) \
	| perl -pi -e 's#(.*)(TODO:?)(\[.+\])#\3 \1\2\3#i;' \
	| perl -pi -e 's#\s*//\s*TODO:\s*(\[.+\])\s*#\n\1     #i;'


#
# Run static analysis tools
#
static_analysis:
	contrib/bin/cpplint.py --filter=-whitespace,-build/include,-runtime/int,-runtime/references,-readability/streams,-build/namespaces,-build/header_guard include/*/*.h src/*/*.C 2>&1 | tee cpplint.txt
	if type cppcheck >/dev/null 2>&1; then (for file in include/*/*.h src/*/*.C; do cppcheck `contrib/bin/libmesh-config --include` --enable=all -f $$file; done) 2>&1 | tee cppcheck.txt; fi


#
# Dependencies
#
.depend: $(srcfiles) $(appsrcfiles) $(csrcfiles) $(headerfiles)
	@$(perl) ./contrib/bin/make_dependencies.pl $(foreach i, $(wildcard include/ include/*/ contrib/*/ contrib/*/*/), -I./$(i)) "-S\$$(obj-suffix)" $(srcfiles) $(appsrcfiles) $(csrcfiles) > .depend



#
# Special rules for building some source files.  Must be *after* the 'all' target.
#

# This rule suppresses deprecated header warnings when compiling files with
# VTK headers.
src/mesh/vtk_io.$(obj-suffix) : src/mesh/vtk_io.C
	@echo "Compiling using special rule for vtk_io.C (in "$(mode)" mode) "$<"..."
	@$(libmesh_CXX) $(NODEPRECATEDFLAG) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c $< -o $@




###############################################################################
include .depend


# Local Variables:
# mode: makefile
# End:
