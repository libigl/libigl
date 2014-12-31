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
point.o: /usr/include/i386/_types.h
point.o: /usr/include/sys/_pthread/_pthread_types.h
point.o: /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_size_t.h
point.o: /usr/include/sys/_types/_null.h /usr/include/sys/stdio.h
point.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
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
point.o: /usr/include/sys/_pthread/_pthread_attr_t.h
point.o: /usr/include/sys/_types/_sigaltstack.h
point.o: /usr/include/sys/_types/_ucontext.h
point.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h
point.o: /usr/include/sys/resource.h /usr/include/stdint.h
point.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
point.o: /usr/include/sys/_types/_int32_t.h
point.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
point.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
point.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
point.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
point.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
point.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
point.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
point.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
point.o: /usr/include/sys/_types/_ct_rune_t.h
point.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
point.o: /usr/include/machine/types.h /usr/include/i386/types.h
point.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
point.o: /usr/include/errno.h /usr/include/sys/errno.h
point.o: /usr/include/sys/_types/_errno_t.h
vertex.o: ../../include/vertex.h ../../include/j_mesh.h /usr/include/stdio.h
vertex.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
vertex.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
vertex.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
vertex.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
vertex.o: /usr/include/i386/_types.h
vertex.o: /usr/include/sys/_pthread/_pthread_types.h
vertex.o: /usr/include/sys/_types/_va_list.h
vertex.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
vertex.o: /usr/include/sys/stdio.h /usr/include/sys/_types/_off_t.h
vertex.o: /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h
vertex.o: /usr/include/secure/_common.h /usr/include/math.h
vertex.o: /usr/include/limits.h /usr/include/machine/limits.h
vertex.o: /usr/include/i386/limits.h /usr/include/i386/_limits.h
vertex.o: /usr/include/sys/syslimits.h ../../include/list.h
vertex.o: ../../include/point.h ../../include/edge.h ../../include/triangle.h
vertex.o: ../../include/jmesh.h ../../include/matrix.h ../../include/tin.h
vertex.o: /usr/include/stdlib.h /usr/include/sys/wait.h
vertex.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
vertex.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
vertex.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
vertex.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
vertex.o: /usr/include/mach/i386/_structs.h
vertex.o: /usr/include/sys/_pthread/_pthread_attr_t.h
vertex.o: /usr/include/sys/_types/_sigaltstack.h
vertex.o: /usr/include/sys/_types/_ucontext.h
vertex.o: /usr/include/sys/_types/_sigset_t.h
vertex.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
vertex.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
vertex.o: /usr/include/sys/_types/_int16_t.h
vertex.o: /usr/include/sys/_types/_int32_t.h
vertex.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
vertex.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
vertex.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
vertex.o: /usr/include/sys/_types/_uintptr_t.h
vertex.o: /usr/include/_types/_intmax_t.h /usr/include/_types/_uintmax_t.h
vertex.o: /usr/include/sys/_types/_timeval.h /usr/include/machine/endian.h
vertex.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
vertex.o: /usr/include/libkern/_OSByteOrder.h
vertex.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
vertex.o: /usr/include/sys/_types/_ct_rune_t.h
vertex.o: /usr/include/sys/_types/_rune_t.h
vertex.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
vertex.o: /usr/include/i386/types.h /usr/include/sys/_types/_dev_t.h
vertex.o: /usr/include/sys/_types/_mode_t.h /usr/include/errno.h
vertex.o: /usr/include/sys/errno.h /usr/include/sys/_types/_errno_t.h
edge.o: ../../include/edge.h ../../include/j_mesh.h /usr/include/stdio.h
edge.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
edge.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
edge.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
edge.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
edge.o: /usr/include/i386/_types.h /usr/include/sys/_pthread/_pthread_types.h
edge.o: /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_size_t.h
edge.o: /usr/include/sys/_types/_null.h /usr/include/sys/stdio.h
edge.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
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
triangle.o: /usr/include/sys/_pthread/_pthread_types.h
triangle.o: /usr/include/sys/_types/_va_list.h
triangle.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h
triangle.o: /usr/include/sys/stdio.h /usr/include/sys/_types/_off_t.h
triangle.o: /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h
triangle.o: /usr/include/secure/_common.h /usr/include/math.h
triangle.o: /usr/include/limits.h /usr/include/machine/limits.h
triangle.o: /usr/include/i386/limits.h /usr/include/i386/_limits.h
triangle.o: /usr/include/sys/syslimits.h ../../include/list.h
triangle.o: ../../include/point.h ../../include/vertex.h ../../include/edge.h
triangle.o: ../../include/matrix.h ../../include/tin.h /usr/include/stdlib.h
triangle.o: /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h
triangle.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
triangle.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
triangle.o: /usr/include/i386/signal.h /usr/include/machine/_mcontext.h
triangle.o: /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h
triangle.o: /usr/include/sys/_pthread/_pthread_attr_t.h
triangle.o: /usr/include/sys/_types/_sigaltstack.h
triangle.o: /usr/include/sys/_types/_ucontext.h
triangle.o: /usr/include/sys/_types/_sigset_t.h
triangle.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
triangle.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
triangle.o: /usr/include/sys/_types/_int16_t.h
triangle.o: /usr/include/sys/_types/_int32_t.h
triangle.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
triangle.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
triangle.o: /usr/include/_types/_uint64_t.h
triangle.o: /usr/include/sys/_types/_intptr_t.h
triangle.o: /usr/include/sys/_types/_uintptr_t.h
triangle.o: /usr/include/_types/_intmax_t.h /usr/include/_types/_uintmax_t.h
triangle.o: /usr/include/sys/_types/_timeval.h /usr/include/machine/endian.h
triangle.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
triangle.o: /usr/include/libkern/_OSByteOrder.h
triangle.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
triangle.o: /usr/include/sys/_types/_ct_rune_t.h
triangle.o: /usr/include/sys/_types/_rune_t.h
triangle.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
triangle.o: /usr/include/i386/types.h /usr/include/sys/_types/_dev_t.h
triangle.o: /usr/include/sys/_types/_mode_t.h
tin.o: ../../include/tin.h ../../include/jmesh.h ../../include/j_mesh.h
tin.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
tin.o: /usr/include/sys/_symbol_aliasing.h
tin.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
tin.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
tin.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
tin.o: /usr/include/i386/_types.h /usr/include/sys/_pthread/_pthread_types.h
tin.o: /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_size_t.h
tin.o: /usr/include/sys/_types/_null.h /usr/include/sys/stdio.h
tin.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
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
tin.o: /usr/include/sys/_pthread/_pthread_attr_t.h
tin.o: /usr/include/sys/_types/_sigaltstack.h
tin.o: /usr/include/sys/_types/_ucontext.h
tin.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h
tin.o: /usr/include/sys/resource.h /usr/include/stdint.h
tin.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
tin.o: /usr/include/sys/_types/_int32_t.h /usr/include/sys/_types/_int64_t.h
tin.o: /usr/include/_types/_uint8_t.h /usr/include/_types/_uint16_t.h
tin.o: /usr/include/_types/_uint32_t.h /usr/include/_types/_uint64_t.h
tin.o: /usr/include/sys/_types/_intptr_t.h
tin.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
tin.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
tin.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
tin.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
tin.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
tin.o: /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h
tin.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
tin.o: /usr/include/i386/types.h /usr/include/sys/_types/_dev_t.h
tin.o: /usr/include/sys/_types/_mode_t.h /usr/include/string.h
tin.o: /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h
tin.o: /usr/include/strings.h /usr/include/secure/_string.h
io.o: ../../include/jmesh.h ../../include/j_mesh.h /usr/include/stdio.h
io.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
io.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
io.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
io.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
io.o: /usr/include/i386/_types.h /usr/include/sys/_pthread/_pthread_types.h
io.o: /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_size_t.h
io.o: /usr/include/sys/_types/_null.h /usr/include/sys/stdio.h
io.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
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
io.o: /usr/include/sys/_pthread/_pthread_attr_t.h
io.o: /usr/include/sys/_types/_sigaltstack.h
io.o: /usr/include/sys/_types/_ucontext.h /usr/include/sys/_types/_sigset_t.h
io.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
io.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
io.o: /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h
io.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
io.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
io.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
io.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
io.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
io.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
io.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
io.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
io.o: /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h
io.o: /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h
io.o: /usr/include/i386/types.h /usr/include/sys/_types/_dev_t.h
io.o: /usr/include/sys/_types/_mode_t.h /usr/include/string.h
io.o: /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h
io.o: /usr/include/strings.h /usr/include/secure/_string.h
io.o: /usr/include/ctype.h /usr/include/runetype.h
io.o: /usr/include/sys/_types/_wint_t.h
checkAndRepair.o: ../../include/jmesh.h ../../include/j_mesh.h
checkAndRepair.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
checkAndRepair.o: /usr/include/sys/_symbol_aliasing.h
checkAndRepair.o: /usr/include/sys/_posix_availability.h
checkAndRepair.o: /usr/include/Availability.h
checkAndRepair.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
checkAndRepair.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
checkAndRepair.o: /usr/include/i386/_types.h
checkAndRepair.o: /usr/include/sys/_pthread/_pthread_types.h
checkAndRepair.o: /usr/include/sys/_types/_va_list.h
checkAndRepair.o: /usr/include/sys/_types/_size_t.h
checkAndRepair.o: /usr/include/sys/_types/_null.h /usr/include/sys/stdio.h
checkAndRepair.o: /usr/include/sys/_types/_off_t.h
checkAndRepair.o: /usr/include/sys/_types/_ssize_t.h
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
checkAndRepair.o: /usr/include/sys/_pthread/_pthread_attr_t.h
checkAndRepair.o: /usr/include/sys/_types/_sigaltstack.h
checkAndRepair.o: /usr/include/sys/_types/_ucontext.h
checkAndRepair.o: /usr/include/sys/_types/_sigset_t.h
checkAndRepair.o: /usr/include/sys/_types/_uid_t.h
checkAndRepair.o: /usr/include/sys/resource.h /usr/include/stdint.h
checkAndRepair.o: /usr/include/sys/_types/_int8_t.h
checkAndRepair.o: /usr/include/sys/_types/_int16_t.h
checkAndRepair.o: /usr/include/sys/_types/_int32_t.h
checkAndRepair.o: /usr/include/sys/_types/_int64_t.h
checkAndRepair.o: /usr/include/_types/_uint8_t.h
checkAndRepair.o: /usr/include/_types/_uint16_t.h
checkAndRepair.o: /usr/include/_types/_uint32_t.h
checkAndRepair.o: /usr/include/_types/_uint64_t.h
checkAndRepair.o: /usr/include/sys/_types/_intptr_t.h
checkAndRepair.o: /usr/include/sys/_types/_uintptr_t.h
checkAndRepair.o: /usr/include/_types/_intmax_t.h
checkAndRepair.o: /usr/include/_types/_uintmax_t.h
checkAndRepair.o: /usr/include/sys/_types/_timeval.h
checkAndRepair.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
checkAndRepair.o: /usr/include/sys/_endian.h
checkAndRepair.o: /usr/include/libkern/_OSByteOrder.h
checkAndRepair.o: /usr/include/libkern/i386/_OSByteOrder.h
checkAndRepair.o: /usr/include/alloca.h /usr/include/sys/_types/_ct_rune_t.h
checkAndRepair.o: /usr/include/sys/_types/_rune_t.h
checkAndRepair.o: /usr/include/sys/_types/_wchar_t.h
checkAndRepair.o: /usr/include/machine/types.h /usr/include/i386/types.h
checkAndRepair.o: /usr/include/sys/_types/_dev_t.h
checkAndRepair.o: /usr/include/sys/_types/_mode_t.h /usr/include/string.h
checkAndRepair.o: /usr/include/sys/_types/_rsize_t.h
checkAndRepair.o: /usr/include/sys/_types/_errno_t.h /usr/include/strings.h
checkAndRepair.o: /usr/include/secure/_string.h
