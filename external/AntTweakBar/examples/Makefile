####### Compiler, tools and options


#---- Release
CXXCFG   	= -O3
LFLAGS		= -L../lib
OUT_DIR		= bin64
#---- Debug
#CXXCFG   	= -g -D_DEBUG
#LFLAGS		= -Wl -L../lib/debug64
#OUT_DIR	= debug64


CXX      	= gcc
CXXFLAGS 	= $(CXXCFG) -Wall -fno-strict-aliasing
INCPATH  	= -I../include -I/usr/local/include -I/usr/X11R6/include -I/usr/include
LIBS     	= -L/usr/X11R6/lib -lAntTweakBar -lGL -lGLU -lX11 -lXext -lXmu -lXrandr -lpthread -lm

DEL_FILE 	= rm -f
DEL_DIR  	= rmdir
NO_STDERR	= 2> /dev/null
EXP_SH      = '\#!/bin/sh'
EXP_PATH	= 'export LD_LIBRARY_PATH=`dirname $$0`/../../lib ; $$0.out'

####### Files


SRC_FILES	= TwSimpleGLFW.c TwSimpleGLUT.c TwSimpleSDL.c TwAdvanced1.cpp TwString.cpp TwDualGLUT.c


####### Build rules


#first:	depend all
first:  all

all: 	Makefile $(SRC_FILES)

	@echo "===== Build TwSimpleGLUT ===="
	$(CXX) $(CXXFLAGS) $(INCPATH) TwSimpleGLUT.c $(LFLAGS) -lglut $(LIBS) -o $(OUT_DIR)/TwSimpleGLUT.out
	@echo $(EXP_SH) > $(OUT_DIR)/TwSimpleGLUT
	@echo $(EXP_PATH) >> $(OUT_DIR)/TwSimpleGLUT
	@chmod +x $(OUT_DIR)/TwSimpleGLUT

	@echo "===== Build TwSimpleSDL ===="
	$(CXX) $(CXXFLAGS) $(INCPATH) TwSimpleSDL.c $(LFLAGS) -lSDL $(LIBS) -o $(OUT_DIR)/TwSimpleSDL.out
	@echo $(EXP_SH) > $(OUT_DIR)/TwSimpleSDL
	@echo $(EXP_PATH) >> $(OUT_DIR)/TwSimpleSDL
	@chmod +x $(OUT_DIR)/TwSimpleSDL

	@echo "===== Build TwSimpleGLFW ===="
	$(CXX) $(CXXFLAGS) $(INCPATH) TwSimpleGLFW.c $(LFLAGS) -lglfw $(LIBS) -o $(OUT_DIR)/TwSimpleGLFW.out
	@echo $(EXP_SH) > $(OUT_DIR)/TwSimpleGLFW
	@echo $(EXP_PATH) >> $(OUT_DIR)/TwSimpleGLFW
	@chmod +x $(OUT_DIR)/TwSimpleGLFW

	@echo "===== Build TwDualGLUT ===="
	$(CXX) $(CXXFLAGS) $(INCPATH) TwDualGLUT.c $(LFLAGS) -lglut $(LIBS) -o $(OUT_DIR)/TwDualGLUT.out
	@echo $(EXP_SH) > $(OUT_DIR)/TwDualGLUT
	@echo $(EXP_PATH) >> $(OUT_DIR)/TwDualGLUT
	@chmod +x $(OUT_DIR)/TwDualGLUT

	@echo "===== Build TwAdvanced1 ===="
	$(CXX) $(CXXFLAGS) $(INCPATH) TwAdvanced1.cpp $(LFLAGS) -lglfw $(LIBS) -o $(OUT_DIR)/TwAdvanced1.out
	@echo $(EXP_SH) > $(OUT_DIR)/TwAdvanced1
	@echo $(EXP_PATH) >> $(OUT_DIR)/TwAdvanced1
	@chmod +x $(OUT_DIR)/TwAdvanced1

	@echo "===== Build TwString ===="
	$(CXX) $(CXXFLAGS) $(INCPATH) TwString.cpp $(LFLAGS) -lglfw $(LIBS) -o $(OUT_DIR)/TwString.out
	@echo $(EXP_SH) > $(OUT_DIR)/TwString
	@echo $(EXP_PATH) >> $(OUT_DIR)/TwString
	@chmod +x $(OUT_DIR)/TwString

# append dependencies to this Makefile
#depend:
#	@echo "===== Make dependencies ====="
#	makedepend -Y
#	makedepend -a -Y -- $(CXXFLAGS) $(INCPATH) -- $(SRC_FILES) $(NO_STDERR)


# clean temporary files
clean:
	@echo "===== Clean ====="
	-$(DEL_FILE) *.o
	-$(DEL_FILE) *~ core *.core *.stackdump
	-$(DEL_FILE) debug64/*
	-$(DEL_DIR) debug64

	
####### DEPENDENCIES

TwSimpleGLFW.o: ../include/AntTweakBar.h
TwSimpleGLUT.o: ../include/AntTweakBar.h
TwDualGLUT.o:   ../include/AntTweakBar.h
TwSimpleSDL.o:  ../include/AntTweakBar.h
TwAdvanced1.o:  ../include/AntTweakBar.h
TwString.o:     ../include/AntTweakBar.h
