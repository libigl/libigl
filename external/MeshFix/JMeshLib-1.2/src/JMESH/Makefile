include ../../makeconf

INCPATH = -I../../include

CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH)

SOURCES.C    =  jmesh.cpp

OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/jmesh.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/jmesh.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/jmesh.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

jmesh.o: ../../include/jmesh.h ../../include/j_mesh.h /usr/include/stdio.h
jmesh.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
jmesh.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
jmesh.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
jmesh.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
jmesh.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
jmesh.o: /usr/include/secure/_common.h /usr/include/math.h
jmesh.o: /usr/include/architecture/i386/math.h /usr/include/limits.h
jmesh.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
jmesh.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
jmesh.o: ../../include/list.h ../../include/point.h ../../include/vertex.h
jmesh.o: ../../include/edge.h ../../include/triangle.h ../../include/matrix.h
jmesh.o: ../../include/tin.h /usr/include/stdlib.h /usr/include/sys/wait.h
jmesh.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
jmesh.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
jmesh.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
jmesh.o: /usr/include/machine/_structs.h /usr/include/sys/resource.h
jmesh.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
jmesh.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
jmesh.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
jmesh.o: /usr/include/machine/types.h /usr/include/i386/types.h
jmesh.o: /usr/include/string.h /usr/include/strings.h
jmesh.o: /usr/include/secure/_string.h
