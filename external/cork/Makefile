# Warning: This makefile rampantly assumes it is being run in UNIX
#   Also, there may be some gnu tool chain assumptions: gcc and gmake?

# ***********
# * PREFACE *
# ***********------+
# | Subdirectories |
# +----------------+
# declare subdirectories of the source directory
SUBDIRECTORIES := util file_formats math isct mesh rawmesh accel

# make sure the subdirectories are mirrored in
# the obj/ debug/ and depend/ directories
# (HACK: use this dummy variable to get access to shell commands)
SHELL_HACK := $(shell mkdir -p bin lib include)
SHELL_HACK := $(shell mkdir -p $(addprefix obj/,$(SUBDIRECTORIES)))
SHELL_HACK := $(shell mkdir -p $(addprefix debug/,$(SUBDIRECTORIES)))
SHELL_HACK := $(shell mkdir -p $(addprefix depend/,$(SUBDIRECTORIES)))
# also make a directory to expose headers in
#SHELL_HACK := $(shell mkdir -p $(addprefix include/,$(SUBDIRECTORIES)))

# +----------+
# | Platform |
# +----------+-------------------
# | evaluates to 'Darwin' on mac
# | evaluates to ??? on ???
# +------------------------------
PLATFORM := $(shell uname)

-include makeConstants


# *********
# * TOOLS *
# *********--+
# | Programs |
# +----------+
CPP11_FLAGS := -std=c++11 
#include ../../build/Makefile.conf
#CC  := ${GG}
#CXX := ${GG}
ifeq ($(PLATFORM),Darwin) # on mac
    CPP11_FLAGS := $(CPP11_FLAGS)  -Wno-c++11-extensions
endif
RM  := rm
CP  := cp

# +--------------+
# | Option Flags |
# +--------------+
# Let the preprocessor know which platform we're on
ifeq ($(PLATFORM),Darwin)
  CONFIG  := $(CONFIG) -DMACOSX
else
  CONFIG  := $(CONFIG)
endif

# include paths (flattens the hierarchy for include directives)
INC       := -I src/ $(addprefix -I src/,$(SUBDIRECTORIES))
# get the location of GMP header files from makeConstants include file
GMP_INC_DIR := /opt/local/include/
GMPINC    := -I $(GMP_INC_DIR)
INC       := $(INC) $(GMPINC)

# use the second line to disable profiling instrumentation
# PROFILING := -pg
PROFILING :=
CCFLAGS   := -Wall $(INC) $(CONFIG) -O2 -DNDEBUG $(PROFILING)
CXXFLAGS  := $(CCFLAGS) $(CPP11_FLAGS)
CCDFLAGS  := -Wall $(INC) $(CONFIG) -ggdb
CXXDFLAGS := $(CCDFLAGS)

# Place the location of GMP libraries here
GMP_LIB_DIR := /opt/local/lib
GMPLD     := -L$(GMP_LIB_DIR) -lgmpxx -lgmp
# static version...
#GMPLD     := $(GMP_LIB_DIR)/libgmpxx.a $(GMP_LIB_DIR)/libgmp.a

LINK         := $(CXXFLAGS) $(GMPLD)
LINKD        := $(CXXDFLAGS) $(GMPLD)
ifeq ($(PLATFORM),Darwin)
  LINK  := $(LINK) -Wl,-no_pie
  LINKD := $(LINK) -Wl,-no_pie
endif


# ***********************
# * SOURCE DECLARATIONS *
# ***********************-----------------+
# | SRCS defines a generic bag of sources |
# +---------------------------------------+
MATH_SRCS    := 
UTIL_SRCS    := timer log
ISCT_SRCS    := empty3d quantization
MESH_SRCS    := 
RAWMESH_SRCS := 
ACCEL_SRCS   := 
FILE_SRCS    := files ifs off
SRCS         := \
    cork \
    $(addprefix math/,$(MATH_SRCS))\
    $(addprefix util/,$(UTIL_SRCS))\
    $(addprefix isct/,$(ISCT_SRCS))\
    $(addprefix mesh/,$(MESH_SRCS))\
    $(addprefix rawmesh/,$(RAWMESH_SRCS))\
    $(addprefix accel/,$(ACCEL_SRCS))\
    $(addprefix file_formats/,$(FILE_SRCS))


# +-----------------------------------+
# | HEADERS defines headers to export |
# +-----------------------------------+
MATH_HEADERS      := vec.h bbox.h ray.h
UTIL_HEADERS      := prelude.h memPool.h iterPool.h shortVec.h \
                     unionFind.h
ISCT_HEADERS      := unsafeRayTriIsct.h \
                     ext4.h fixext4.h gmpext4.h absext4.h \
                     quantization.h fixint.h \
                     empty3d.h \
                     triangle.h
RAWMESH_HEADERS   := rawMesh.h rawMesh.tpp
MESH_HEADERS      := mesh.h mesh.decl.h \
                     mesh.tpp mesh.topoCache.tpp \
                     mesh.remesh.tpp mesh.isct.tpp mesh.bool.tpp
ACCEL_HEADERS     := aabvh.h
FILE_HEADERS      := files.h
HEADERS           := \
    cork.h
#    $(addprefix math/,$(MATH_HEADERS))\
#    $(addprefix util/,$(UTIL_HEADERS))\
#    $(addprefix isct/,$(ISCT_HEADERS))\
#    $(addprefix mesh/,$(MESH_HEADERS))\
#    $(addprefix rawmesh/,$(RAWMESH_HEADERS))\
#    $(addprefix accel/,$(ACCEL_HEADERS))\
#    $(addprefix file_formats/,$(FILE_HEADERS))
HEADER_COPIES     := $(addprefix include/,$(HEADERS))

# +-----------------------------+
# | stand alone program sources |
# +-----------------------------+
MAIN_SRC := \
    $(SRCS) \
    main

# +---------------------------------------+
# | all sources for dependency generation |
# +---------------------------------------+
ALL_SRCS     := \
    $(SRCS)\
    main
DEPENDS := $(addprefix depend/,$(addsuffix .d,$(ALL_SRCS)))

# +--------------------------------+
# | Object Aggregation for Targets |
# +--------------------------------+

OBJ               := $(addprefix obj/,$(addsuffix .o,$(SRCS))) \
                     obj/isct/triangle.o
DEBUG             := $(addprefix debug/,$(addsuffix .o,$(SRCS))) \
                     obj/isct/triangle.o

MAIN_OBJ          := $(addprefix obj/,$(addsuffix .o,$(MAIN_SRC))) \
                     obj/isct/triangle.o
MAIN_DEBUG        := $(addprefix debug/,$(addsuffix .o,$(MAIN_SRC))) \
                     obj/isct/triangle.o

LIB_TARGET_NAME   := cork

# *********
# * RULES *
# *********------+
# | Target Rules |
# +--------------+
all: lib/lib$(LIB_TARGET_NAME).a includes \
     bin/off2obj bin/cork
debug: lib/lib$(LIB_TARGET_NAME)debug.a includes

lib/lib$(LIB_TARGET_NAME).a: $(OBJ)
	@echo "Bundling $@"
	@ar rcs $@ $(OBJ)

lib/lib$(LIB_TARGET_NAME)debug.a: $(DEBUG)
	@echo "Bundling $@"
	@ar rcs $@ $(DEBUG)

bin/cork: $(MAIN_OBJ)
	@echo "Linking cork command line tool"
	@$(CXX) -o bin/cork $(MAIN_OBJ) $(LINK)

bin/off2obj: obj/off2obj.o
	@echo "Linking off2obj"
	@$(CXX) -o bin/off2obj obj/off2obj.o $(LINK)

# +------------------------------+
# | Specialized File Build Rules |
# +------------------------------+

obj/isct/triangle.o: src/isct/triangle.c
	@echo "Compiling the Triangle library"
	@$(CC) -x c -O2 -DNO_TIMER \
               -DREDUCED \
               -DCDT_ONLY -DTRILIBRARY \
               -Wall -DANSI_DECLARATORS \
               -o obj/isct/triangle.o -c src/isct/triangle.c

# +------------------------------------+
# | Generic Source->Object Build Rules |
# +------------------------------------+
obj/%.o: src/%.cpp
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

debug/%.o: src/%.cpp
	@echo "Compiling $@"
	@$(CXX) $(CXXDFLAGS) -o $@ -c $<

# dependency file build rules
depend/%.d: src/%.cpp
	@$(CXX) $(CXXFLAGS) -MM $< | \
        sed -e 's@^\(.*\)\.o:@depend/$*.d debug/$*.o obj/$*.o:@' > $@

# +-------------------+
# | include copy rule |
# +-------------------+---------------------
# | This rule exists to safely propagate
# | header file dependencies to other
# | targets that depend on the common code
# +-----------------------------------------
includes: $(HEADER_COPIES)

include/%.h: src/%.h
	@echo "updating $@"
	@cp $< $@
#also support template implementation files
include/%.tpp: src/%.tpp
	@echo "updating $@"
	@cp $< $@

# +---------------+
# | cleaning rule |
# +---------------+
clean:
	-@$(RM) -r obj depend debug include bin lib
	-@$(RM) bin/off2obj
#	-@$(RM) gmon.out
	-@$(RM) lib/lib$(LIB_TARGET_NAME).a
	-@$(RM) lib/lib$(LIB_TARGET_NAME)debug.a

-include $(DEPENDS)
