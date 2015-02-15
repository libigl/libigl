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
binTree.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
binTree.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
binTree.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
binTree.o: ../../include/binTree.h
dijkstraGraph.o: ../../include/dijkstraGraph.h ../../include/heap.h
dijkstraGraph.o: ../../include/j_mesh.h /usr/include/stdio.h
dijkstraGraph.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
dijkstraGraph.o: /usr/include/sys/_posix_availability.h
dijkstraGraph.o: /usr/include/Availability.h
dijkstraGraph.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
dijkstraGraph.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
dijkstraGraph.o: /usr/include/i386/_types.h
dijkstraGraph.o: /usr/include/sys/_types/_va_list.h
dijkstraGraph.o: /usr/include/sys/_types/_size_t.h
dijkstraGraph.o: /usr/include/sys/_types/_null.h /usr/include/secure/_stdio.h
dijkstraGraph.o: /usr/include/secure/_common.h /usr/include/math.h
dijkstraGraph.o: /usr/include/limits.h /usr/include/machine/limits.h
dijkstraGraph.o: /usr/include/i386/limits.h /usr/include/i386/_limits.h
dijkstraGraph.o: /usr/include/sys/syslimits.h ../../include/graph.h
dijkstraGraph.o: ../../include/list.h
heap.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
heap.o: /usr/include/sys/_symbol_aliasing.h
heap.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
heap.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
heap.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
heap.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
heap.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
heap.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
heap.o: /usr/include/stdlib.h /usr/include/sys/wait.h
heap.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
heap.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
heap.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
heap.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
heap.o: /usr/include/mach/i386/_structs.h
heap.o: /usr/include/sys/_types/_sigaltstack.h
heap.o: /usr/include/sys/_types/_ucontext.h
heap.o: /usr/include/sys/_types/_pthread_attr_t.h
heap.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h
heap.o: /usr/include/sys/resource.h /usr/include/sys/_types/_timeval.h
heap.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
heap.o: /usr/include/sys/_endian.h /usr/include/alloca.h
heap.o: /usr/include/sys/_types/_ct_rune_t.h
heap.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
heap.o: /usr/include/machine/types.h /usr/include/i386/types.h
heap.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
heap.o: /usr/include/sys/_types/_int32_t.h /usr/include/sys/_types/_int64_t.h
heap.o: /usr/include/sys/_types/_intptr_t.h
heap.o: /usr/include/sys/_types/_uintptr_t.h
heap.o: /usr/include/sys/_types/___offsetof.h
heap.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
heap.o: ../../include/heap.h ../../include/j_mesh.h /usr/include/math.h
heap.o: /usr/include/limits.h /usr/include/machine/limits.h
heap.o: /usr/include/i386/limits.h /usr/include/i386/_limits.h
heap.o: /usr/include/sys/syslimits.h
list.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
list.o: /usr/include/sys/_symbol_aliasing.h
list.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
list.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
list.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
list.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
list.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
list.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
list.o: /usr/include/stdlib.h /usr/include/sys/wait.h
list.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
list.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
list.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
list.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
list.o: /usr/include/mach/i386/_structs.h
list.o: /usr/include/sys/_types/_sigaltstack.h
list.o: /usr/include/sys/_types/_ucontext.h
list.o: /usr/include/sys/_types/_pthread_attr_t.h
list.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h
list.o: /usr/include/sys/resource.h /usr/include/sys/_types/_timeval.h
list.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
list.o: /usr/include/sys/_endian.h /usr/include/alloca.h
list.o: /usr/include/sys/_types/_ct_rune_t.h
list.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
list.o: /usr/include/machine/types.h /usr/include/i386/types.h
list.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
list.o: /usr/include/sys/_types/_int32_t.h /usr/include/sys/_types/_int64_t.h
list.o: /usr/include/sys/_types/_intptr_t.h
list.o: /usr/include/sys/_types/_uintptr_t.h
list.o: /usr/include/sys/_types/___offsetof.h
list.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
list.o: ../../include/list.h ../../include/jqsort.h
clusterGraph.o: ../../include/clusterGraph.h ../../include/heap.h
clusterGraph.o: ../../include/j_mesh.h /usr/include/stdio.h
clusterGraph.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
clusterGraph.o: /usr/include/sys/_posix_availability.h
clusterGraph.o: /usr/include/Availability.h
clusterGraph.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
clusterGraph.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
clusterGraph.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
clusterGraph.o: /usr/include/sys/_types/_size_t.h
clusterGraph.o: /usr/include/sys/_types/_null.h /usr/include/secure/_stdio.h
clusterGraph.o: /usr/include/secure/_common.h /usr/include/math.h
clusterGraph.o: /usr/include/limits.h /usr/include/machine/limits.h
clusterGraph.o: /usr/include/i386/limits.h /usr/include/i386/_limits.h
clusterGraph.o: /usr/include/sys/syslimits.h ../../include/graph.h
clusterGraph.o: ../../include/list.h
graph.o: ../../include/graph.h ../../include/list.h /usr/include/stdio.h
graph.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
graph.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
graph.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
graph.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
graph.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
graph.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
graph.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
matrix.o: /usr/include/math.h /usr/include/sys/cdefs.h
matrix.o: /usr/include/sys/_symbol_aliasing.h
matrix.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
matrix.o: /usr/include/AvailabilityInternal.h ../../include/matrix.h
matrix.o: /usr/include/stdio.h /usr/include/_types.h
matrix.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
matrix.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
matrix.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
matrix.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
matrix.o: ../../include/list.h
