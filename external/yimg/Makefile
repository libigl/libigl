SOURCES = $(CXX_OBJECTS:%.o=%.cpp) $(C_OBJECTS:%.o=%.c) \
			$(MAIN_OBJECT:%.o=%.cpp)
LIB_OBJECTS = $(CXX_OBJECTS) $(C_OBJECTS)
CXX_OBJECTS = YImage.o
MAIN_OBJECT = showpng.o

ARCH = $(shell uname)

TARGET = showpng
LIBTARGET = libyimg.a

COMPILER = g++
#COMPILER = gfilt
WARNINGFLAGS = -Wall -Wshadow -Woverloaded-virtual -Wno-sign-compare
# WARNINGFLAGS += -Winit-self

# Does this help?  -falign-loops=16
COMPILERFLAGS = $(WARNINGFLAGS)

#COMPILERFLAGS += -g
COMPILERFLAGS += -O3
COMPILERFLAGS += -DNDEBUG

#COMPILERFLAGS += -fno-inline
#COMPILERFLAGS += -finline-functions
COMPILERFLAGS += -fkeep-inline-functions
#COMPILERFLAGS += -Winline

AR = /usr/bin/ar
ARFLAGS = -rus

MAKEDEPEND = g++ -MM
DEPENDFILE = Makefile.depend

INCLUDEFLAGS =
LINKTIMEFLAGS =
include Makefile.$(ARCH)

# clear out implicit rules
.SUFFIXES:

all: $(LIBTARGET) $(TARGET)

$(TARGET): $(MAIN_OBJECT) $(LIBTARGET)
	$(COMPILER) -o $@ \
	$(MAIN_OBJECT) \
	$(LIBTARGET) \
	$(FRAMEWORKS) \
	$(INCLUDEFLAGS) $(LINKTIMEFLAGS)
#	$(INCLUDEFLAGS) $(LINKTIMEFLAGS) $(LIBTARGET) $(LINKTIMEFLAGS) \


#$(TARGET): $(LIB_OBJECTS) $(MAIN_OBJECT)
#	$(COMPILER) -o $@ \
#	$(LIB_OBJECTS) $(MAIN_OBJECT) \
#	$(INCLUDEFLAGS) $(LINKTIMEFLAGS) \
#	$(FRAMEWORKS)

library: $(LIBTARGET)

relink: cleantarget all
	

$(LIBTARGET): $(LIB_OBJECTS)
	$(AR) $(ARFLAGS) $@ $(LIB_OBJECTS)

%.o: %.cpp
	$(COMPILER) $(COMPILERFLAGS) $(INCLUDEFLAGS) -c $<

%.o: %.c
	$(COMPILER) $(COMPILERFLAGS) $(INCLUDEFLAGS) -c $<

.PHONY : clean tidy
clean: tidy
	rm -f $(TARGET)
	rm -f $(LIBTARGET)
	rm -f $(DEPENDFILE)
	touch $(DEPENDFILE)

tidy:
	rm -f $(LIB_OBJECTS) $(MAIN_OBJECT)

cleantarget:
	rm -f $(TARGET)

depend:
	$(MAKEDEPEND) $(INCLUDEFLAGS) $(SOURCES) > $(DEPENDFILE)

include $(DEPENDFILE)

$(DEPENDFILE):
	touch $(DEPENDFILE)
