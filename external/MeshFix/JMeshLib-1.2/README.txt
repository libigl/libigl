----------------------------
JMeshLib - Version 1.2                                                                
----------------------------

by Marco Attene
                                                                    
Consiglio Nazionale delle Ricerche                                        
Istituto di Matematica Applicata e Tecnologie Informatiche                
Sezione di Genova                                                         
IMATI-GE / CNR                                                            
                                                                         
JMeshLib provides a framework to work with manifold triangle meshes.
It implements an edge-based data structure with all its fundamental
functionalities (i.e., file I/O, mesh construction/destruction, traversal).
It is written in C++ and includes support for reading and writing the
following file formats:
OFF (http://shape.cs.princeton.edu/benchmark/documentation/off_format.html)
PLY (http://www.cs.unc.edu/~geom/Powerplant/Ply.doc)
STL (http://www.sdsc.edu/tmf/Stl-specs/stl.html)
VER-TRI (proprietary format used at IMATI-GE / CNR)
and partially:
IV 2.1, VRML 1.0, VRML 2.0, OBJ.

In contrast to other generic libraries dealing with surface meshes,
JMeshLib includes tools to automatically fix the most common problems
present in surface meshes coming from laser scanning (conversion to
oriented manifold, topological noise detection, hole filling, removal
of degenerate faces, ...) through a clear and easy-to-learn C++ API.

This package provides pre-compiled static libraries Windows only (lib/jmesh.lib).
For Linux you must compile the source tree yourself.

See the comments within the source files for details or use doxygen to
produce documentation in a more readable format.
The file "tin.h" is a good starting point.

-------------------
System Rrequirements
--------------------

JMeshLib v1.0 has been extensively tested on 32 bit PCs running either:
 - ELF Linux with standard development tools (gcc/g++)
 - Windows OS with MSVC 8.0 (Visual C++ 2005)

From version 1.1 it should compile and work properly on 64-bit Linux
machines too. Uncomment the proper line in the 'makeconf' to create
64-bit-enabled binaries.

-------------------
Building the tree
-------------------

On WINDOWS: double-click on minJMeshLib.vcproj and press F7.
On Linux: type 'make' on the command line.

On Linux systems you need the command 'makedepend'. If you
don't have it installed on your system, you may need to
download and install the 'imake' package.

-------------------
Using the library
-------------------

On WINDOWS:
Add the path to jmesh's include dir to your project's include path
AND
add the path to jmesh.lib to your linker's command line

On Linux:
Add the path to jmesh's include dir to your project's include path
AND
add the path to jmesh.a to your linker's command line

---------
Copyright
---------

JMeshLib is

Copyright(C) 2006-2011: IMATI-GE / CNR                                         
                                                                          
All rights reserved.                                                      
                                                                          
This program is free software; you can redistribute it and/or modify      
it under the terms of the GNU General Public License as published by      
the Free Software Foundation; either version 2 of the License, or         
(at your option) any later version.                                       
                                                                          
This program is distributed in the hope that it will be useful,           
but WITHOUT ANY WARRANTY; without even the implied warranty of            
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             
GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          
for more details.                                                         
