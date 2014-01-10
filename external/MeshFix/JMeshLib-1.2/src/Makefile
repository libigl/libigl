include ../makeconf


#------------------------ PATHS ------------------------

SOURCEDIR    =  .
INCLUDEDIR   =  ../include

#------------------------ OBJECTS ------------------------

PRIMITIVES = 	OBJECTS/primitives.o
MESH_STRUCT = OBJECTS/mesh_structure.o
JMESH = 	OBJECTS/jmesh.o

OBJECTS = $(PRIMITIVES) $(MESH_STRUCT) $(JMESH)

#-------------------- RULES --------------------------------------

lib: $(OBJECTS)
	$(AR) -rs ../lib/libjmesh.a $(OBJECTS)
	$(RM) -f $(OBJECTS)

depend:
	$(MAKE) -C PRIMITIVES depend
	$(MAKE) -C MESH_STRUCTURE depend
	$(MAKE) -C JMESH depend

clean:
	$(MAKE) -C PRIMITIVES clean
	$(MAKE) -C MESH_STRUCTURE clean
	$(MAKE) -C JMESH clean

#----------------------- GLOBAL DEPENDANCIES ---------------------

$(PRIMITIVES):
	$(MAKE) -C PRIMITIVES
$(MESH_STRUCT):
	$(MAKE) -C MESH_STRUCTURE
$(JMESH):
	$(MAKE) -C JMESH
