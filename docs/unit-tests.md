# Unit Tests for [libigl](https://github.com/libigl/libigl)
[![Build Status](https://travis-ci.org/libigl/libigl-unit-tests.svg?branch=master)](https://travis-ci.org/libigl/libigl-unit-tests)

Get started with

```
git clone --recursive git@github.com:libigl/libigl-unit-tests.git
```

## Dependencies

[googletest](https://github.com/google/googletest) is a submodule


## Build and test

Use `cmake` to generate a `Makefile` that will build _and test_ upon issuing
`make`:

```
mkdir build
cd build
cmake ..
```

Then build and test with

```
make test
```

This will first compile the tests and then immediately run the tests. If tests
are succeeding you should see output similar to:

```
Test project /usr/local/libigl-unit-tests/build
    Start 1: run_igl_mosek_tests
1/4 Test #1: run_igl_mosek_tests ..............***Exception: Other  0.00 sec
    Start 2: run_igl_boolean_tests
2/4 Test #2: run_igl_boolean_tests ............   Passed    1.12 sec
    Start 3: run_igl_cgal_tests
3/4 Test #3: run_igl_cgal_tests ...............   Passed    2.46 sec
    Start 4: run_igl_tests
```

Alternatively, to get more detailed output you can call ctest directly with the
verbose flag:

```
GTEST_COLOR=1 ctest --verbose
```

You'll see outputs for each individual test:

```
UpdateCTestConfiguration  from :/usr/local/libigl-unit-tests/build/DartConfiguration.tcl
UpdateCTestConfiguration  from :/usr/local/libigl-unit-tests/build/DartConfiguration.tcl
Test project /usr/local/libigl-unit-tests/build
Constructing a list of tests
Done constructing a list of tests
Updating test list for fixtures
Added 0 tests to meet fixture requirements
Checking test dependency graph...
Checking test dependency graph end
test 1
    Start 1: run_igl_mosek_tests

1: Test command: /usr/local/libigl-unit-tests/build/include/igl/mosek/igl_mosek_tests
1: Test timeout computed to be: 9.99988e+06
1: [==========] Running 1 test from 1 test case.
1: [----------] Global test environment set-up.
1: [----------] 1 test from mosek_bbw
1: [ RUN      ] mosek_bbw.decimated_knight
1: /usr/local/libigl-unit-tests/include/igl/mosek/bbw.cpp:25: Failure
1: Expected: ((Wmo-W_groundtruth).array().abs().maxCoeff()) < (1e-3), actual: 0.00287895 vs 0.001
1: [  FAILED  ] mosek_bbw.decimated_knight (2126 ms)
1: [----------] 1 test from mosek_bbw (2126 ms total)
1: 
1: [----------] Global test environment tear-down
1: [==========] 1 test from 1 test case ran. (2126 ms total)
1: [  PASSED  ] 0 tests.
1: [  FAILED  ] 1 test, listed below:
1: [  FAILED  ] mosek_bbw.decimated_knight
1: 
1:  1 FAILED TEST
1/4 Test #1: run_igl_mosek_tests ..............***Failed    2.14 sec
test 2
    Start 2: run_igl_boolean_tests

2: Test command: /usr/local/libigl-unit-tests/build/include/igl/copyleft/boolean/igl_boolean_tests
2: Test timeout computed to be: 9.99988e+06
2: [==========] Running 3 tests from 1 test case.
2: [----------] Global test environment set-up.
2: [----------] 3 tests from MeshBoolean
2: [ RUN      ] MeshBoolean.TwoCubes
...
```


## Generating new tests

Many libigl functions act on triangle meshes. To make it easy to add a new test
for a libigl function, we have a script that will output some boilerplate
testing code. So if you issue:

```
cd include/igl/
../../scripts/new.sh myfun
```

This will create a test file `myfun.cpp` containing:

```
#include <test_common.h>
#include <igl/myfun.h>

class myfun : public ::testing::TestWithParam<std::string> {};

TEST_P(myfun, change_to_meaningful_name)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::SparseMatrix<double> L;
  // Load example mesh: GetParam() will be name of mesh file
  test_common::load_mesh(GetParam(), V, F);
  // ASSERT_EQ(a,b);
  // ASSERT_TRUE(a==b);
  // ASSERT_NEAR(a,b,1e-15)
  // ASSERT_LT(a,1e-12);
}

INSTANTIATE_TEST_CASE_P
(
 all_meshes,
 myfun,
 ::testing::ValuesIn(test_common::all_meshes()),
 test_common::string_test_name
);
```

Add a call to `igl::myfun` and an assertion (e.g., `ASSERT_EQ`) and this will
add a test for _all_ meshes in the `data/` folder. (Should also change
`change_to_meaningful_name` to a meaningful name based on what you're testing).

## Conventions

When naming a test for a function `igl::extra::function_name` use:

```cpp
TEST(extra_function_name, meaning_test_name)
{
  ...
}
```

where `meaning_test_name` could identify the type of test or the type of data
being used.

### Example

The test for `igl::copyleft::cgal::order_facets_around_edges` in
`include/igl/copyleft/cgal/order_facets_around_edges.cpp` is:

```cpp
TEST(copyleft_cgal_order_facets_around_edges, TripletFaces)
{
  ...
}
```

which tests this function on example data containing a triplet of faces.

## Guarantees

None.

(Obviously?) The presence of a unit test here for some function (e.g.,
`igl::cotmatrix`) is not a guarantee or even an endorsement of the notion that
the libigl function `igl::cotmatrix` is bug free or "fully tested" or "heavily
tested" or even "adequately tested".

## Need work?

Some of the most used libigl functions

```bash
grep -hr "^#include \"" ../libigl/include/igl | sed -e 's/\(\.\.\/\)//g' | sort | uniq -c | sort
```

still don't have unit tests.
