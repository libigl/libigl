# Libigl Style Guidelines

Libigl is used and developed by many people. This document highlights some
style guidelines for _developers_ of the library, but also acts as
best-practices for users.

## One function, one .h/.cpp pair [filefunction]

The structure of libigl is very flat and function-based. For every
function/sub-routine, create a single .h and .cpp file. For example, if you have
a function that determines connected components from a face list `F` you would
create the header `connected_components.h` and `connected_components.cpp` and the only
function defined should be `void connected_components(const ... F, ... C)`. If the
implementation of `connected_components` requires a subroutine to compute an
adjacency matrix then _create another pair_ `adjacency_matrix.h` and
`adjacency_matrix.cpp` with a single function `void adjacency_matrix(const ... F, ... A)`.

### Example
Here is an example function that would be defined in
`include/igl/example_fun.h` and implemented in `include/igl/example_fun.cpp`.

#### `example_fun.h`

```cpp
// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 [Your Name] [your email address]
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#ifndef IGL_EXAMPLE_FUN_H
#define IGL_EXAMPLE_FUN_H

#include "igl_inline.h"

namespace igl
{
  // This is an example of a function, it takes a templated parameter and
  // shovels it into cout
  //
  // Input:
  //   input  some input of a Printable type
  // Returns true for the sake of returning something
  template <typename Printable>
  IGL_INLINE bool example_fun(const Printable & input);
}

#ifndef IGL_STATIC_LIBRARY
#  include "example_fun.cpp"
#endif

#endif
```

#### `example_fun.cpp`

```
// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 [Your Name] [your email address]
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#include "igl/example_fun.h"
#include <iostream>

template <typename Printable>
IGL_INLINE bool igl::example_fun(const Printable & input)
{
  using namespace std;
  cout<<"example_fun: "<<input<<endl;
  return true;
}

#ifdef IGL_STATIC_LIBRARY
template bool igl::example_fun<double>(const double& input);
template bool igl::example_fun<int>(const int& input);
#endif
```


### Avoid static "helper" functions

Strive to encapsulate sub-functions that could possibly be useful outside of
the implementation of your current function. This might mean abstracting the
interface a bit. If it doesn't dramatically effect performance then create a
new pair of .h/.cpp files with this sub-function.

#### Lambda functions

If encapsulation in a separate file is not possible or does not make sense,
then avoid crowding the namespace by creating lambda functions within the
function implementation.

These lambda functions must still be documented with clear [input and output
arguments](#headerdocumentation). Avoid using full capturing of all automatic
variables: do not use `[&]` or `[=]`. Rather specify each captured variable
individually.

### Avoid "helper" classes

Libigl is built around the high-performance paradigm of "struct of arrays"
rather than "array of structs". The way we achieve this is to avoid classes and
pass "basic types" directly. The price we pay is long function interfaces, but
this increases code reuse dramatically. A "basic type" in our context is a
Eigen type, stl type, or basic C type.

## Header Documentation

Each function prototype should be well documented in its corresponding .h
header file. A typical documentation consists of four parts:

```cpp
// [A human readable description of what the function does.]
//
// Inputs:
//   [variable name of first (const) input]   [dimensions and description of
//     this input variable]
//   [variable name of second (const) input]   [dimensions and description of
//     this input variable]
//   ...
// Outputs:
//   [variable name of first output ]   [dimensions and description of this
//     output variable]
//   [variable name of second output ]   [dimensions and description of this
//     output variable]
//   ...
// Returns [description of return value]
```

### Example

For example the header `barycenter.h`

```
// Computes the barycenter of every simplex
//
// Inputs:
//   V  #V by dim matrix of vertex coordinates
//   F  #F by simplex_size  matrix of indices of simplex corners into V
// Output:
//   BC  #F by dim matrix of 3d vertices
//
```

## Const inputs

All input parameters should be demarcated `const`. If an input is also an
output than consider exposing two parameters (one `const`) or be sure to list
the variable under both `// Inputs:` and `// Outputs:` in the header comments.

## Reference parameters

All but simple types should be passed by reference (e.g. `Matrix & mat`) rather
than pointers (e.g. `Matrix * mat`) or value (e.g. `Matrix mat`).


## Returns vs output parameters

All functions should be implemented with at least one overload that has a
`void` or simple return type (e.g. `bool` on success/failure). With this
implementation its then possible to write an overload that returns a single
output. Please see [Templating with Eigen](#templatingwitheigen).

For example:

```cpp
template <typename Atype>
void adjacency_matrix(const ... & F, Eigen::SparseMatrix<AType> & A);

template <typename Atype>
Eigen::SparseMatrix<Atype> adjacency_matrix(const ... & F);
```

## Templating with Eigen

Functions taking Eigen dense matrices/arrays as inputs and outputs (but **not**
return arguments), should template on top of `Eigen::PlainObjectBase`. **Each
parameter** should be derived using its own template.

For example,

```cpp
template <typename DerivedV, typename DerivedF, typename DerivedBC>
void barycenter(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedBC> & BC);
```

The `Derived*` template encodes the scalar type (e.g. `double`, `int`), the
number of rows and cols at compile time, and the data storage (Row-major vs.
column-major). 

Returning Eigen types is discouraged. In cases where the size and scalar type
are a fixed **and matching** function of an input `Derived*` template, then
return that `Derived*` type. **Do not** return
`Eigen::PlainObjectBase<...>` types. For example, this function scales fits a
given set of points to the unit cube. The return is a new set of vertex
positions so its type should _match_ that of the input points:

```cpp
template <typename DerivedV>
void DerivedV fit_to_unit_cube(const Eigen::PlainObjectBase<DerivedV> & V);
```

To implement this function, it is **required** to implement a more generic
output-argument version and call that. So a full implementation looks like:

In `igl/fit_in_unit_cube.h`:

```cpp
template <typename DerivedV, typename DerivedW>
void fit_to_unit_cube(
  const Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedW> & W);
template <typename DerivedV>
void DerivedV fit_to_unit_cube(const Eigen::PlainObjectBase<DerivedV> & V);
```

In `igl/fit_in_unit_cube.cpp`:

```
template <typename DerivedV, typename DerivedW>
void fit_to_unit_cube(
  const Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedW> & W)
{
  W = (V.rowwise()-V.colwise().minCoeff()).array() /
    (V.maxCoeff()-V.minCoeff());
}

template <typename DerivedV>
void DerivedV fit_to_unit_cube(const Eigen::PlainObjectBase<DerivedV> & V)
{
  DerivedV W;
  fit_to_unit_cube(V,W);
  return W;
}
```

Notice that `W` is declared as a `DerivedV` type and **not**
`Eigen::PlainObjectBase<DerivedV>` type.

**Note:** Not all functions are suitable for returning Eigen types. For example
`igl::barycenter` above outputs a #F by dim list of barycenters. Returning a
`DerivedV` type would be inappropriate since the number of rows in `DerivedV`
will be #V and may not match the number of rows in `DerivedF` (#F).

## Function naming conventions 

Functions (and [thus also files](#filefunction)) should have simple,
descriptive names using lowercase letters and underscores between words. Avoid
unnecessary prefaces. For example, instead of `compute_adjacency_matrix`,
`construct_adjacency_matrix`, `extract_adjacency_matrix`,
`get_adjacency_matrix`, or `set_adjacency_matrix` just call the function
`adjacency_matrix`.

## Variable naming conventions

Libigl prefers short (even single character) variable names _with heavy
documentation_ in the comments in the header file or above the declaration of
the function. When possible use `V` to mean a list of vertex positions and `F`
to mean a list of faces/triangles.

## Class naming conventions 

Classes should be avoided. When naming a class use CamelCase (e.g.
SortableRow.h).

## Enum naming convertion

Enums types should be placed in the appropriate `igl::` namespace and should be
named in CamelCase (e.g. `igl::SolverStatus`) and instances should be named in
ALL_CAPS with underscores between words and prefaced with the name of the enum.
For example:

```cpp
namespace igl
{
  enum SolverStatus
  {
    // Good
    SOLVER_STATUS_CONVERGED = 0,
    // OK
    SOLVER_STATUS_MAX_ITER = 1,
    // Bad
    SOLVER_STATUS_ERROR = 2,
    NUM_SOLVER_STATUSES = 3,
  };
};
```

### Exception for file IO

For legacy reasons, file reading and writing functions use a different naming
convention. A functions reading a `.xyz` file should be named `readXYZ` and a
function writing `.xyz` files should be names `writeXYZ`.

## `using namespace ...` in global scope

Writing `using namespace std;`, `using namespace Eigen;` etc. outside of a
global scope is strictly forbidden. Place these lines at the top of each
function instead.


## Namespaces and external dependencies

Functions in the main library (directly in `include/igl`) should only depend on
Eigen and stl. These functions should have the `igl::` namespace.

Functions with other dependencies should be placed into
appropriate sub-directories (e.g. if `myfunction` depends on tetgen then create
`igl/copyleft/tetgen/myfunction.h` and `igl/copyleft/tetgen/myfunction.cpp` and give the function
the namespace `igl::copyleft::tetgen::myfunction`.

### copyleft subdirectory/namespace 

Dependencies that require users of libigl to release their projects open source
(e.g. GPL) are considered aggressively "copyleft" and should be placed in the
`include/igl/copyleft/` sub-directory and `igl::copyleft::` namespace.

## Assertions

Be generous with assertions and always identify the assertion with strings:

```cpp
assert(m < n && "m must be less than n");
```

## ifndef include guard

Every header file should be wrapped in an `#ifndef` compiler directive. The
name of the guard should be in direct correspondence with the path of the .h
file. For example, `include/igl/copyleft/tetgen/tetrahedralize.h` should be

```cpp
#ifndef IGL_COPYLEFT_TETGEN_TETRAHEDRALIZE_H
#define IGL_COPYLEFT_TETGEN_TETRAHEDRALIZE_H
...
#endif
```

## Spaces vs. tabs indentation

Do not use tabs. Use 2 spaces for each indentation level.

## Max line length

Limit lines to 80 characters. Break up long lines into many operations (this
also helps performance).

## Include order

`#include` directives at the top of a .h or .cpp file should be sorted
according to a simple principle: place headers of files most likely to be
edited by you first. This means for
`include/igl/copyleft/tetgen/tetrahedralize.cpp` you might see

```cpp
// [Includes of headers in this directory]
#include "tetrahedralize.h"
#include "mesh_to_tetgenio.h"
#include "tetgenio_to_tetmesh.h"
// [Includes of headers in this project]
#include "../../matrix_to_list.h"
#include "../../list_to_matrix.h"
#include "../../boundary_facets.h"
// [Includes of headers of related projects]
#include <Eigen/Core>
// [Includes of headers of standard libraries]
#include <cassert>
#include <iostream>
```

## Placement of includes

Whenever possible `#include` directives should be placed in the `.cpp`
implementation file rather than the `.h` header file.
