# libigl - A simple C++ geometry processing library

![](libigl-teaser.png)

<https://github.com/libigl/libigl/>

libigl is a simple C++ geometry processing library. We have a wide
functionality including construction of sparse discrete differential geometry
operators and finite-elements matrices such as the contangent Laplacian and
diagonalized mass matrix, simple facet and edge-based topology data structures,
mesh-viewing utilities for OpenGL and GLSL, and many core functions for matrix
manipulation which make [Eigen](http://eigen.tuxfamily.org) feel a lot more
like MATLAB.

It is first and foremost a header library. Each header file contains a single
function. Most are tailored to operate on a generic triangle mesh stored in an
n-by-3 matrix of vertex positions V and an m-by-3 matrix of triangle indices F.
The library may also be [compiled](build/) into a statically linked
library, for faster compile times with your projects.

We use the [Eigen](http://eigen.tuxfamily.org) library heavily in our code. Our
group prototypes a lot in MATLAB, and we have a useful [conversion
table](matlab-to-eigen.html) from
MATLAB to libigl/Eigen.

## Tutorial

As of version 1.0, libigl includes an introductory
[tutorial](tutorial/tutorial.html) that covers
its basic functionalities.

## Installation
Libigl is a *header* library. You do **not** need to build anything to install.
Simply add `igl/` to your include path and include relevant headers. Here is a
small "Hello, World" program:

```cpp
#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
int main()
{
  Eigen::MatrixXd V(4,2);
  V<<0,0,
     1,0,
     1,1,
     0,1;
  Eigen::MatrixXi F(2,3);
  F<<0,1,2,
     0,2,3;
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  std::cout<<"Hello, mesh: "<<std::endl<<L*V<<std::endl;
  return 0;
}
```

If you save this in `hello.cpp`, then you could compile this with (assuming
Eigen is installed in /opt/local/include/eigen3):

```bash
gcc -I/opt/local/include/eigen3 -I./igl/ hello.cpp -o hello
```

Running `./hello` would then produce

```
Hello, mesh:
 0.5  0.5
-0.5  0.5
-0.5 -0.5
 0.5 -0.5
```

## Dependencies
Dependencies are on a per-include basis and the majority of the functions in
libigl depends only on the [Eigen](http://eigen.tuxfamily.org) library.

For more information see our [tutorial](tutorial/tutorial.html).

### GCC and the optional CGAL dependency
The `include/igl/cgal/*.h` headers depend on CGAL. It has come to our attention
that CGAL does not work properly with GCC 4.8. To the best of our knowledge,
GCC 4.7 and clang will work correctly.

### OpenMP and Windows
Some of our functions will take advantage of OpenMP if available. However, it
has come to our attention that Visual Studio + Eigen does not work properly
with OpenMP. Since OpenMP only improves performance without affecting
functionality we recommend avoiding OpenMP on Windows or proceeding with
caution.

## Download
You can keep up to date by cloning a read-only copy of our GitHub
[repository](https://github.com/libigl).

## How to contribute

If you are interested in joining development, please fork the repository and
submit a [pull request](https://help.github.com/articles/using-pull-requests/)
with your changes.

## License
libigl is primarily [MPL2](http://www.mozilla.org/MPL/2.0/) licensed
([FAQ](http://www.mozilla.org/MPL/2.0/FAQ.html)). Some files contain
third-party code under other licenses. We're currently in the processes of
identifying these and marking appropriately.

## Attribution
If you use libigl in your academic projects, please cite the papers we
implement as appropriate. To cite the library in general, you could use this
BibTeX entry:

```bibtex
@misc{libigl,
  title = {{libigl}: A simple {C++} geometry processing library},
  author = {Alec Jacobson and Daniele Panozzo and others},
  note = {http://libigl.github.io/libigl/},
  year = {2015},
}
```

## Contact

Libigl is a group endeavor led by [Alec
Jacobson](http://www.cs.columbia.edu/~jacobson/) and [Daniele
Panozzo](http://www.inf.ethz.ch/personal/dpanozzo/). Please [contact
us](mailto:alecjacobson@gmail.com,daniele.panozzo@gmail.com) if you have
questions or comments. We are happy to get feedback!

If you're using libigl in your projects, quickly [drop us a
note](mailto:alecjacobson@gmail.com,daniele.panozzo@gmail.com). Tell us who you
are and what you're using it for. This helps us apply for funding and justify
spending time maintaining this.

If you find bugs or have problems please use our [github issue tracking
page](https://github.com/libigl/libigl/issues).

## Copyright
2015 Alec Jacobson, Daniele Panozzo, Olga Diamanti, Christian Schüller, Kenshi
Takayama, Leo Sacht, Wenzel Jacob, Nico Pietroni, Amir Vaxman

![](tutorial/images/libigl-logo.jpg)
