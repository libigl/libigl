# License

Libigl is primarily licensed under MPL2:

- http://www.mozilla.org/MPL/2.0/
- http://www.mozilla.org/MPL/2.0/FAQ.html

All `.h` and `.cpp` _files_ directly in `include/igl` (but not necessarily in
sub-directories) are subject only to the terms of the MPL2; they should not
include any code that is covered by other/less-permissive licenses.

The `.h` and `.cpp` _files_ in sub-directories of `include/igl` allow libigl to
integrate with external third-party libraries (e.g., those in `external/`) and
are subject to the MPL2, _**and**_ also the terms of licenses of the
corresponding external library.  The licenses used by these libraries fall under
three categories:

- common "free, non-copyleft licenses" (such as zlib, BSD, MIT, and public
  domain)
  - `include/igl/anttweakbar`
  - `include/igl/embree`
  - `include/igl/opengl`
  - `include/igl/opengl/glfw`
  - `include/igl/opengl2`
  - `include/igl/png`
  - `include/igl/viewer`
  - `include/igl/xml`
- common "copyleft" licences (such as GPL, LGPL, and AGPL)
  - `include/igl/copyleft`
  - `include/igl/copyleft/cgal`
  - `include/igl/copyleft/comiso`
  - `include/igl/copyleft/cork`
  - `include/igl/copyleft/tetgen`
- other "uncommon" licenses or commercial software
  - `include/igl/lim`
  - `include/igl/matlab`
  - `include/igl/mosek`
  - `include/igl/triangle`

The Libigl code that interfaces with "copyleft" libraries is in
`include/igl/copyleft`.  Only include these headers if you are accept the
licensing terms of the corresponding external library.  For example, using
`include/igl/copyleft/tetgen` requires that you accept the terms of the AGPLv3.
