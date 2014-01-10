---------------
MeshFix V1.0
by Marco Attene
IMATI-GE / CNR
---------------

This software takes as input a polygon mesh and produces a copy of the input
where all the occurrences of a specific set of "defects" are corrected. MeshFix
has been designed to correct typical flaws present in RAW DIGITIZED mesh
models, thus it might fail or produce coarse results if run on other sorts of
input meshes (e.g. tessellated CAD models). When the software fails, it appends
a textual description to the file "meshfix.log".

The input is assumed to represent a single CLOSED SOLID OBJECT, thus the output
will be a SINGLE WATERTIGHT TRIANGLE MESH bounding a polyhedron. All the
singularities, self-intersections and degenerate elements are removed from the
input, while regions of the surface without defects are left unmodified.
Accepted input formats are OFF, PLY and STL. Other formats are supported only
partially. See http://jmeshlib.sourceforge.net for details on supported
formats.

-----------------------------
ALGORITHM AND CITATION POLICY
-----------------------------
To better understand how the algorithm works, please refer to the following paper:

   M. Attene.
   A lightweight approach to repairing digitized polygon meshes.
   The Visual Computer, 2010. (c) Springer. DOI: 10.1007/s00371-010-0416-3

This software is based on ideas published therein. If you use MeshFix for
research purposes you should cite the above paper in your published results.
MeshFix cannot be used for commercial purposes without a written permission by
the author.

----------
PARAMETERS
----------
The user may select a threshold angle to be used when assessing whether a triangle is to be considered degenerate or not. Also, the user may want to force the software to keep only the biggest connected component of the input, while considering the others as "noise". If not specified by the user, the default threshold angle is 0, meaning that only triangles which have null area are considered to be degenerate (exact arithmetic is used for such evaluation).

-------------------
COMMAND LINE SYNTAX
-------------------
Usage: MeshFix meshfile [-a epsilon_angle] [-w] [-n]
  Processes 'meshfile' and saves the result to 'meshfile_fixed.off'
  By default, epsilon_angle is 0. If specified, it must be in the range (0 - 2) degrees.
  With '-w', the output is saved in VRML format instead of OFF.
  With '-n', only the biggest input component is kept, otherwise all of them are used.

---------------------
OTHER LAUNCHING MODES
---------------------
A mesh can also be fixed my simply dragging its icon on MeshFix's icon. In this case, only default parameters can be used.
MeshFix does not require any interaction, thus it can be inserted into a script to automatically repair all the models of a given repository (e.g. a folder). If some failure cases occur, they will be logged to "meshfix.log" and thus can be easily located and possibly processed is a second stage through interactive software tools such as ReMESH (http://remesh.sourceforge.net).


---------------------
SOURCE CODE
---------------------

To compile the source code you need the following libraries correctly installed on your PC:

1)
JMeshLib 1.1 or later
DOWNLOAD AT: http://jmeshlib.sourceforge.net

2)
The fast robust geometric predicates by J.R.Shewchuk
DOWNLOAD AT: http://www.cs.cmu.edu/~quake/robust.html
Rename the two files from predicates.* to jrs_predicates.* and copy them as follows:
jrs_predicates.h to JMeshExt-1.0alpha_src\include\
jrs_predicates.c to JMeshExt-1.0alpha_src\src\JRS_Predicates\

3)
OpenNL
DOWNLOAD AT: http://alice.loria.fr/index.php/software/4-library/23-opennl.html

After that, compile JMeshExt (VC project in "JMeshExt-1.0alpha_src\vc8") and then MeshFix (VC project in "vc8").

---------
Copyright
---------

MeshFix is

Copyright(C) 2010: IMATI-GE / CNR                                       

All rights reserved.                                                      

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.                                       

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License (http://www.gnu.org/licenses/gpl.txt) for more details.                                                         
