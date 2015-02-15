include ../../makeconf

INCPATH = -I../../include
CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH) -I../RANGECODER

SOURCES.C    =  point.cpp\
		vertex.cpp\
		edge.cpp\
		triangle.cpp\
		tin.cpp\
		io.cpp\
		checkAndRepair.cpp
		
OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/mesh_structure.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/mesh_structure.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/mesh_structure.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

point.o: ../../include/point.h ../../include/j_mesh.h /usr/include/stdio.h
point.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
point.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
point.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
point.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
point.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
point.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
point.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
point.o: /usr/include/math.h /usr/include/limits.h
point.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
point.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
point.o: /usr/include/stdlib.h /usr/include/sys/wait.h
point.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
point.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
point.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
point.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
point.o: /usr/include/mach/i386/_structs.h
point.o: /usr/include/sys/_types/_sigaltstack.h
point.o: /usr/include/sys/_types/_ucontext.h
point.o: /usr/include/sys/_types/_pthread_attr_t.h
point.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h
point.o: /usr/include/sys/resource.h /usr/include/sys/_types/_timeval.h
point.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
point.o: /usr/include/sys/_endian.h /usr/include/alloca.h
point.o: /usr/include/sys/_types/_ct_rune_t.h
point.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
point.o: /usr/include/machine/types.h /usr/include/i386/types.h
point.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
point.o: /usr/include/sys/_types/_int32_t.h
point.o: /usr/include/sys/_types/_int64_t.h
point.o: /usr/include/sys/_types/_intptr_t.h
point.o: /usr/include/sys/_types/_uintptr_t.h
point.o: /usr/include/sys/_types/___offsetof.h
point.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
point.o: /usr/include/errno.h /usr/include/sys/errno.h
vertex.o: ../../include/vertex.h ../../include/j_mesh.h /usr/include/stdio.h
vertex.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
vertex.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
vertex.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
vertex.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
vertex.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
vertex.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
vertex.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
vertex.o: /usr/include/math.h /usr/include/limits.h
vertex.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
vertex.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
vertex.o: ../../include/list.h ../../include/point.h ../../include/edge.h
vertex.o: ../../include/triangle.h ../../include/jmesh.h
vertex.o: ../../include/matrix.h ../../include/tin.h /usr/include/stdlib.h
vertex.o: /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h
vertex.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
vertex.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
vertex.o: /usr/include/i386/signal.h /usr/include/machine/_mcontext.h
vertex.o: /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h
vertex.o: /usr/include/sys/_types/_sigaltstack.h
vertex.o: /usr/include/sys/_types/_ucontext.h
vertex.o: /usr/include/sys/_types/_pthread_attr_t.h
vertex.o: /usr/include/sys/_types/_sigset_t.h
vertex.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
vertex.o: /usr/include/sys/_types/_timeval.h /usr/include/machine/endian.h
vertex.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
vertex.o: /usr/include/alloca.h /usr/include/sys/_types/_ct_rune_t.h
vertex.o: /usr/include/sys/_types/_rune_t.h
vertex.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
vertex.o: /usr/include/i386/types.h /usr/include/sys/_types/_int8_t.h
vertex.o: /usr/include/sys/_types/_int16_t.h
vertex.o: /usr/include/sys/_types/_int32_t.h
vertex.o: /usr/include/sys/_types/_int64_t.h
vertex.o: /usr/include/sys/_types/_intptr_t.h
vertex.o: /usr/include/sys/_types/_uintptr_t.h
vertex.o: /usr/include/sys/_types/___offsetof.h
vertex.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
vertex.o: /usr/include/errno.h /usr/include/sys/errno.h
edge.o: ../../include/edge.h ../../include/j_mesh.h /usr/include/stdio.h
edge.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
edge.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
edge.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
edge.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
edge.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
edge.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
edge.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
edge.o: /usr/include/math.h /usr/include/limits.h
edge.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
edge.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
edge.o: ../../include/vertex.h ../../include/list.h ../../include/point.h
edge.o: ../../include/triangle.h ../../include/jmesh.h ../../include/matrix.h
edge.o: ../../include/tin.h
triangle.o: ../../include/triangle.h ../../include/jmesh.h
triangle.o: ../../include/j_mesh.h /usr/include/stdio.h
triangle.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
triangle.o: /usr/include/sys/_posix_availability.h
triangle.o: /usr/include/Availability.h /usr/include/AvailabilityInternal.h
triangle.o: /usr/include/_types.h /usr/include/sys/_types.h
triangle.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
triangle.o: /usr/include/sys/_types/_va_list.h
triangle.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
triangle.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
triangle.o: /usr/include/math.h /usr/include/limits.h
triangle.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
triangle.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
triangle.o: ../../include/list.h ../../include/point.h ../../include/vertex.h
triangle.o: ../../include/edge.h ../../include/matrix.h ../../include/tin.h
triangle.o: /usr/include/stdlib.h /usr/include/sys/wait.h
triangle.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
triangle.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
triangle.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
triangle.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
triangle.o: /usr/include/mach/i386/_structs.h
triangle.o: /usr/include/sys/_types/_sigaltstack.h
triangle.o: /usr/include/sys/_types/_ucontext.h
triangle.o: /usr/include/sys/_types/_pthread_attr_t.h
triangle.o: /usr/include/sys/_types/_sigset_t.h
triangle.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
triangle.o: /usr/include/sys/_types/_timeval.h /usr/include/machine/endian.h
triangle.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
triangle.o: /usr/include/alloca.h /usr/include/sys/_types/_ct_rune_t.h
triangle.o: /usr/include/sys/_types/_rune_t.h
triangle.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
triangle.o: /usr/include/i386/types.h /usr/include/sys/_types/_int8_t.h
triangle.o: /usr/include/sys/_types/_int16_t.h
triangle.o: /usr/include/sys/_types/_int32_t.h
triangle.o: /usr/include/sys/_types/_int64_t.h
triangle.o: /usr/include/sys/_types/_intptr_t.h
triangle.o: /usr/include/sys/_types/_uintptr_t.h
triangle.o: /usr/include/sys/_types/___offsetof.h
triangle.o: /usr/include/sys/_types/_dev_t.h
triangle.o: /usr/include/sys/_types/_mode_t.h
tin.o: ../../include/tin.h ../../include/jmesh.h ../../include/j_mesh.h
tin.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
tin.o: /usr/include/sys/_symbol_aliasing.h
tin.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
tin.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
tin.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
tin.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
tin.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
tin.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
tin.o: /usr/include/math.h /usr/include/limits.h
tin.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
tin.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
tin.o: ../../include/list.h ../../include/point.h ../../include/vertex.h
tin.o: ../../include/edge.h ../../include/triangle.h ../../include/matrix.h
tin.o: /usr/include/stdlib.h /usr/include/sys/wait.h
tin.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
tin.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
tin.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
tin.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
tin.o: /usr/include/mach/i386/_structs.h
tin.o: /usr/include/sys/_types/_sigaltstack.h
tin.o: /usr/include/sys/_types/_ucontext.h
tin.o: /usr/include/sys/_types/_pthread_attr_t.h
tin.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h
tin.o: /usr/include/sys/resource.h /usr/include/sys/_types/_timeval.h
tin.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
tin.o: /usr/include/sys/_endian.h /usr/include/alloca.h
tin.o: /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h
tin.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
tin.o: /usr/include/i386/types.h /usr/include/sys/_types/_int8_t.h
tin.o: /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h
tin.o: /usr/include/sys/_types/_int64_t.h /usr/include/sys/_types/_intptr_t.h
tin.o: /usr/include/sys/_types/_uintptr_t.h
tin.o: /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h
tin.o: /usr/include/sys/_types/_mode_t.h /usr/include/string.h
tin.o: /usr/include/secure/_string.h
io.o: ../../include/jmesh.h ../../include/j_mesh.h /usr/include/stdio.h
io.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
io.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
io.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
io.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
io.o: /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h
io.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
io.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
io.o: /usr/include/math.h /usr/include/limits.h /usr/include/machine/limits.h
io.o: /usr/include/i386/limits.h /usr/include/i386/_limits.h
io.o: /usr/include/sys/syslimits.h ../../include/list.h ../../include/point.h
io.o: ../../include/vertex.h ../../include/edge.h ../../include/triangle.h
io.o: ../../include/matrix.h ../../include/tin.h /usr/include/stdlib.h
io.o: /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h
io.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
io.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
io.o: /usr/include/i386/signal.h /usr/include/machine/_mcontext.h
io.o: /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h
io.o: /usr/include/sys/_types/_sigaltstack.h
io.o: /usr/include/sys/_types/_ucontext.h
io.o: /usr/include/sys/_types/_pthread_attr_t.h
io.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h
io.o: /usr/include/sys/resource.h /usr/include/sys/_types/_timeval.h
io.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
io.o: /usr/include/sys/_endian.h /usr/include/alloca.h
io.o: /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h
io.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
io.o: /usr/include/i386/types.h /usr/include/sys/_types/_int8_t.h
io.o: /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h
io.o: /usr/include/sys/_types/_int64_t.h /usr/include/sys/_types/_intptr_t.h
io.o: /usr/include/sys/_types/_uintptr_t.h
io.o: /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h
io.o: /usr/include/sys/_types/_mode_t.h /usr/include/string.h
io.o: /usr/include/secure/_string.h /usr/include/ctype.h
io.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h
checkAndRepair.o: ../../include/jmesh.h ../../include/j_mesh.h
checkAndRepair.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
checkAndRepair.o: /usr/include/sys/_symbol_aliasing.h
checkAndRepair.o: /usr/include/sys/_posix_availability.h
checkAndRepair.o: /usr/include/Availability.h
checkAndRepair.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
checkAndRepair.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
checkAndRepair.o: /usr/include/i386/_types.h
checkAndRepair.o: /usr/include/sys/_types/_va_list.h
checkAndRepair.o: /usr/include/sys/_types/_size_t.h
checkAndRepair.o: /usr/include/sys/_types/_null.h
checkAndRepair.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
checkAndRepair.o: /usr/include/math.h /usr/include/limits.h
checkAndRepair.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
checkAndRepair.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
checkAndRepair.o: ../../include/list.h ../../include/point.h
checkAndRepair.o: ../../include/vertex.h ../../include/edge.h
checkAndRepair.o: ../../include/triangle.h ../../include/matrix.h
checkAndRepair.o: ../../include/tin.h ../../include/jqsort.h
checkAndRepair.o: /usr/include/stdlib.h /usr/include/sys/wait.h
checkAndRepair.o: /usr/include/sys/_types/_pid_t.h
checkAndRepair.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
checkAndRepair.o: /usr/include/sys/appleapiopts.h
checkAndRepair.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
checkAndRepair.o: /usr/include/machine/_mcontext.h
checkAndRepair.o: /usr/include/i386/_mcontext.h
checkAndRepair.o: /usr/include/mach/i386/_structs.h
checkAndRepair.o: /usr/include/sys/_types/_sigaltstack.h
checkAndRepair.o: /usr/include/sys/_types/_ucontext.h
checkAndRepair.o: /usr/include/sys/_types/_pthread_attr_t.h
checkAndRepair.o: /usr/include/sys/_types/_sigset_t.h
checkAndRepair.o: /usr/include/sys/_types/_uid_t.h
checkAndRepair.o: /usr/include/sys/resource.h
checkAndRepair.o: /usr/include/sys/_types/_timeval.h
checkAndRepair.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
checkAndRepair.o: /usr/include/sys/_endian.h /usr/include/alloca.h
checkAndRepair.o: /usr/include/sys/_types/_ct_rune_t.h
checkAndRepair.o: /usr/include/sys/_types/_rune_t.h
checkAndRepair.o: /usr/include/sys/_types/_wchar_t.h
checkAndRepair.o: /usr/include/machine/types.h /usr/include/i386/types.h
checkAndRepair.o: /usr/include/sys/_types/_int8_t.h
checkAndRepair.o: /usr/include/sys/_types/_int16_t.h
checkAndRepair.o: /usr/include/sys/_types/_int32_t.h
checkAndRepair.o: /usr/include/sys/_types/_int64_t.h
checkAndRepair.o: /usr/include/sys/_types/_intptr_t.h
checkAndRepair.o: /usr/include/sys/_types/_uintptr_t.h
checkAndRepair.o: /usr/include/sys/_types/___offsetof.h
checkAndRepair.o: /usr/include/sys/_types/_dev_t.h
checkAndRepair.o: /usr/include/sys/_types/_mode_t.h /usr/include/string.h
checkAndRepair.o: /usr/include/secure/_string.h
