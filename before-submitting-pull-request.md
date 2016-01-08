# Before submitting a pull request

There are a variety of things you can do before submitting a pull request that
will reduce the effort on the libigl team to merge your code and increase the
likelihood that the merge ever happens.

  1. Test your code and submit a unit test as part of the pull request
  2. Verify that your code matches the [libigl style
  guidelines](style-guidelines.md)
  3. Run the [exhaustive build test](#exhaustivebuildtest) below

## Exhaustive build test

This script will `git clone` libigl to a temporary directory and build 

  1. the static libigl library, 
  2. the tutorial using the default header only libigl, and 
  3. the tutorial using the static library libigl.
  
Eventually this script should also run the unit tests.

```bash
# In scripts/clone_and_build.sh add your email address to the line:
# `recipients="alecjacobson@gmail.com,youremail@domain.com"`
# In your email client (e.g. gmail) create a filter to prevent emails 
# from your local machine from going to spam
scripts/clone_and_build.sh
```

### Direct test of tutorial using static library

This part of the `clone_and_build.sh` script catches 99% of the compilation
issues that _don't_ show up when testing:

```bash
cd tutorial/
mkdir build-use-static
cd build-use-static
cmake -DCMAKE_BUILD_TYPE=Release -DLIBIGL_USE_STATIC_LIBRARY=ON ..
make
```

A typical issue is a missing template instantiation (symbol not found):

    "void igl::cgal::points_inside_component<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&)", referenced from:
        void igl::cgal::outer_hull<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&) in libiglboolean.a(mesh_boolean.cpp.o)
        void igl::cgal::outer_hull<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&) in libiglboolean.a(mesh_boolean.cpp.o)

This looks like a mess, but the solution is very simple. Copy the chunk inside of the quotes, in this case:

    "void igl::cgal::points_inside_component<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&)"

and paste it at the bottom of the relevant .cpp file with the word template in front of it and a semicolon at then. In this case, in include/igl/cgal/points_inside_component.cpp:

```cpp
#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::cgal::points_inside_component<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
```

Then "rinse and repeat".
