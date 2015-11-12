# Libigl Style Guidelines

This library is shared by many people. This document highlights some style
guidelines for _developers_ of the library, but also act as best-practices for
users.

## One function, one .h/.cpp pair [filefunction]

The structure of libigl is very flat and function-based. For every
function/sub-routine create a single .h and .cpp file. For example, if you have
a function that determines connected components from a face list `F` you would
create the header `connected_components.h` and `connected_components.cpp` and the only
function defined should be `void connected_components(const ... F, ... C)`. If the
implementation of `connected_components` requires a subroutine to compute an
adjacency matrix then _create another pair_ `adjacency_matrix.h` and
`adjacency_matrix.cpp` with a single function `void adjacency_matrix(const ... F, ...
A)`.

### Avoid static "helper" functions

Strive to encapsulate sub-functions that could possibly be useful outside of
the implementation of your current function. This might mean abstracting the
interface a bit. If it doesn't dramatically effect performance then create a
new pair of .h/.cpp files with this sub-function.

#### Lambda functions

If encapsulation in a separate file is not possible or does not make sense,
then avoid crowding the namespace by creating lambda functions within the
function implmentation.

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
//   [variable name of first (const) input]   [dimensions and
//     description of this input variable]
//   [variable name of second (const) input]   [dimensions and
//     description of this input variable]
//   ...
// Outputs:
//   [variable name of first output ]   [dimensions and
//     description of this output variable]
//   [variable name of second output ]   [dimensions and
//     description of this output variable]
//   ...
// Returns [description of return value]
```

### Example

For example the header `barycenter.h`

```
// Computes the barycenter of every simplex
//
// Inputs:
//   V  #V x dim matrix of vertex coordinates
//   F  #F x simplex_size  matrix of indices of simplex corners into V
// Output:
//   BC  #F x dim matrix of 3d vertices
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
output.

For example:

```cpp
template <typename Atype>
void adjacency_matrix(const ... & F, Eigen::SparseMatrix<AType> & A);

template <typename Atype>
Eigen::SparseMatrix<Atype> adjacency_matrix(const ... & F);
```

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
`igl/tetgen/myfunction.h` and `igl/tetgen/myfunction.cpp` and give the function
the namespace `igl::tetgen::myfunction`.

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
file. For example, `include/igl/tetgen/tetrahedralize.h` should be

```cpp
#ifndef IGL_TETGEN_TETRAHEDRALIZE_H
#define IGL_TETGEN_TETRAHEDRALIZE_H
...
#endif
```

## Eigen templates


