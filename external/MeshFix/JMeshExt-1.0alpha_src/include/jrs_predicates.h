// Alec: Copied from tetgen.h
// Alec need to extern for calling in C++ program:
#ifdef __cplusplus
extern "C" {
#endif

#ifndef _predicates_h
#define _predicates_h

#define REAL double

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Robust Geometric predicates                                               //
//                                                                           //
// Geometric predicates are simple tests of spatial relations of a set of d- //
// dimensional points, such as the orientation test and the point-in-sphere  //
// test. Each of these tests is performed by evaluating the sign of a deter- //
// minant of a matrix whose entries are the coordinates of these points.  If //
// the computation is performed by using the floating-point numbers, e.g.,   //
// the single or double numbers in C/C++, roundoff error may cause an incor- //
// rect result. This may either lead to a wrong result or eventually lead to //
// a failure of the program.                                                 //
//                                                                           //
// Various techniques are developed to avoid roundoff errors, such as exact  //
// multi-precision computations, interval arthmetics, adaptive exact arthme- //
// tics, and filtered exact arthmetics, etc. Devillers and Pion give a nice  //
// discussion and comparisons of these techniques for robustly computing the //
// Delaunay triangulations [Devillers and Pion 2002].                        //
//                                                                           //
// The following routines implemented the orientation test and the point-in- //
// sphere test use the adaptive exact floating-point arithmetics [Shewchuk   //
// 1997]. They are generously provided by Jonathan Schewchuk in the public   //
// domain, http://www.cs.cmu.edu/~quake/robust.html. The source code are in  //
// file "predicates.cxx".                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
void exactinit();
REAL orient3d( REAL *pa, REAL *pb, REAL *pc, REAL *pd);
REAL orient2d(REAL *pa, REAL *pb, REAL *pc);

#endif

// Alec: see above
#ifdef __cplusplus
}
#endif 
