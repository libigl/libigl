include ../../makeconf

INCPATH = -I../../include -I../THIRD_PARTY
CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH)

SOURCES.C    =  binTree.cpp\
		dijkstraGraph.cpp\
		heap.cpp\
		list.cpp\
		clusterGraph.cpp\
		graph.cpp\
		jqsort.cpp\
		matrix.cpp

OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/primitives.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/primitives.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/primitives.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

binTree.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
binTree.o: /usr/include/sys/_symbol_aliasing.h
binTree.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
binTree.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
binTree.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
binTree.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
binTree.o: /usr/include/secure/_common.h ../../include/binTree.h
dijkstraGraph.o: ../../include/dijkstraGraph.h ../../include/heap.h
dijkstraGraph.o: ../../include/j_mesh.h /usr/include/stdio.h
dijkstraGraph.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
dijkstraGraph.o: /usr/include/sys/_posix_availability.h
dijkstraGraph.o: /usr/include/Availability.h
dijkstraGraph.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
dijkstraGraph.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
dijkstraGraph.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
dijkstraGraph.o: /usr/include/secure/_common.h /usr/include/math.h
dijkstraGraph.o: /usr/include/architecture/i386/math.h /usr/include/limits.h
dijkstraGraph.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
dijkstraGraph.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
dijkstraGraph.o: ../../include/graph.h ../../include/list.h
heap.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
heap.o: /usr/include/sys/_symbol_aliasing.h
heap.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
heap.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
heap.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
heap.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
heap.o: /usr/include/secure/_common.h /usr/include/stdlib.h
heap.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
heap.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
heap.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
heap.o: /usr/include/sys/_structs.h /usr/include/machine/_structs.h
heap.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
heap.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
heap.o: /usr/include/libkern/_OSByteOrder.h
heap.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
heap.o: /usr/include/machine/types.h /usr/include/i386/types.h
heap.o: ../../include/heap.h ../../include/j_mesh.h /usr/include/math.h
heap.o: /usr/include/architecture/i386/math.h /usr/include/limits.h
heap.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
heap.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
list.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
list.o: /usr/include/sys/_symbol_aliasing.h
list.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
list.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
list.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
list.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
list.o: /usr/include/secure/_common.h /usr/include/stdlib.h
list.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
list.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
list.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
list.o: /usr/include/sys/_structs.h /usr/include/machine/_structs.h
list.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
list.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
list.o: /usr/include/libkern/_OSByteOrder.h
list.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
list.o: /usr/include/machine/types.h /usr/include/i386/types.h
list.o: ../../include/list.h ../../include/jqsort.h
clusterGraph.o: ../../include/clusterGraph.h ../../include/heap.h
clusterGraph.o: ../../include/j_mesh.h /usr/include/stdio.h
clusterGraph.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
clusterGraph.o: /usr/include/sys/_posix_availability.h
clusterGraph.o: /usr/include/Availability.h
clusterGraph.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
clusterGraph.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
clusterGraph.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
clusterGraph.o: /usr/include/secure/_common.h /usr/include/math.h
clusterGraph.o: /usr/include/architecture/i386/math.h /usr/include/limits.h
clusterGraph.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
clusterGraph.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
clusterGraph.o: ../../include/graph.h ../../include/list.h
graph.o: ../../include/graph.h ../../include/list.h /usr/include/stdio.h
graph.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
graph.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
graph.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
graph.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
graph.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
graph.o: /usr/include/secure/_common.h
matrix.o: /usr/include/math.h /usr/include/architecture/i386/math.h
matrix.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
matrix.o: /usr/include/sys/_posix_availability.h ../../include/matrix.h
matrix.o: /usr/include/stdio.h /usr/include/Availability.h
matrix.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
matrix.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
matrix.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
matrix.o: /usr/include/secure/_common.h ../../include/list.h
