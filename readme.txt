e gl Library - A simple c++ mesh library
Copyright 2011 - Daniele Panozzo, Alec Jacobson, Olga Diamanti
Interactive Geometry Lab - ETH Zurich

This is a *header* library. Each header file should contain a single
function. The function may have multiple prototypes. All functions
should use the igl namespace and should adhere to the conventions and
styles listed below. 

= Example =
The following is a toy example of a function called example() in the
file example.h:

#ifndef IGL_EXAMPLE_H
#define IGL_EXAMPLE_H
namespace igl
{
  // An example function that prints "Hello, world"
  inline void example();
}

// Implementation
#include <iostream>
inline void igl::example()
{
  std::cout<<"Hello, world"<<std::endl;
}

= Standards =

This library is shared by many people. Each function prototype should be
well documented.  Write a summary of what the function does and a
description of each template, input and output in each prototype. 

  - All functions must be inlined, otherwise there is trouble when
    linking included headers used in multiple .cpp files
  - Use a single .h file with the same name as the function
  - Do *not* use any .cpp files (this is *header only* library)
  - At least one version of the function should use references for all
    outputs
  - Use wrappers and additional prototypes for returning single-output
    functions' outputs
  - Use c++ references for inputs and outputs rather than pointers or
    pass-by-copy
  - Take the time to write multiple prototypes if you'd like to have
    optional parameters with default values
  - External dependencies (besides Eigen, OpenGL, etc.) are discouraged
  - External dependencies must be clearly identified at the top of each
    file.
  - Single header external dependencies can go in the external/
    directory
  - Do not use the using namespace directive anywhere. This means never
    write:
    "using namespace std;"
     or 
    "using namespace igl;"
    etc.

  = Specific style conventions =

  Each file should be wrapped in an ifndef compiler directive. If the
  file's (and function's) name is example.h, then the ifndef should
  always begin with IGL_, then the function/file name in all caps then
  _H. As in:
#ifndef IGL_EXAMPLE_H
#define IGL_EXAMPLE_H
...
#endif
  
  Each file should begin with prototypes *without implementations* of
  each version of the function. All defined in the igl namespace. Each
  prototype should have its own comments describing what it is doing,
  and its templates, inputs, outputs.

  Implementations should be at the end of each file, separated from the
  prototypes using:
// Implementation

  Any includes, such as std libraries etc. should come after
  the //Implementation separator (whenever feasibly possible).

  Be generous by templating your inputs and outputs. If you do use
  templates, you must document the template just as you document inputs
  and outputs.

  = Useful checks =
  
  Find files that aren't using "inline"

    grep -L inline *

  Find files that aren't using igl namespace

    grep -L "namespace igl" *

  Find files using [TAB] character

    grep -P '\t' *

  Find files that don't contain // Implementation

    grep -L "^\/\/ Implementation" *

  Find files that don't contain #ifndef IGL_

    grep -L "^#ifndef IGL_" *


