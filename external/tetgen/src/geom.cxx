#include "../tetgen.h"
//// geom_cxx /////////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

// PI is the ratio of a circle's circumference to its diameter.
REAL tetgenmesh::PI = 3.14159265358979323846264338327950288419716939937510582;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insphere_s()    Insphere test with symbolic perturbation.                 //
//                                                                           //
// Given four points pa, pb, pc, and pd, test if the point pe lies inside or //
// outside the circumscribed sphere of the four points.                      //
//                                                                           //
// Here we assume that the 3d orientation of the point sequence {pa, pb, pc, //
// pd} is positive (NOT zero), i.e., pd lies above the plane passing through //
// points pa, pb, and pc. Otherwise, the returned sign is flipped.           //
//                                                                           //
// Return a positive value (> 0) if pe lies inside, a negative value (< 0)   //
// if pe lies outside the sphere, the returned value will not be zero.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::insphere_s(REAL* pa, REAL* pb, REAL* pc, REAL* pd, REAL* pe)
{
  REAL sign;

  sign = insphere(pa, pb, pc, pd, pe);
  if (sign != 0.0) {
    return sign;
  }

  // Symbolic perturbation.
  point pt[5], swappt;
  REAL oriA, oriB;
  int swaps, count;
  int n, i;

  pt[0] = pa;
  pt[1] = pb;
  pt[2] = pc;
  pt[3] = pd;
  pt[4] = pe;
  
  // Sort the five points such that their indices are in the increasing
  //   order. An optimized bubble sort algorithm is used, i.e., it has
  //   the worst case O(n^2) runtime, but it is usually much faster.
  swaps = 0; // Record the total number of swaps.
  n = 5;
  do {
    count = 0;
    n = n - 1;
    for (i = 0; i < n; i++) {
      if (pointmark(pt[i]) > pointmark(pt[i+1])) {
        swappt = pt[i]; pt[i] = pt[i+1]; pt[i+1] = swappt;
        count++;
      }
    }
    swaps += count;
  } while (count > 0); // Continue if some points are swapped.

  oriA = orient3d(pt[1], pt[2], pt[3], pt[4]);
  if (oriA != 0.0) {
    // Flip the sign if there are odd number of swaps.
    if ((swaps % 2) != 0) oriA = -oriA;
    return oriA;
  }
  
  oriB = -orient3d(pt[0], pt[2], pt[3], pt[4]);
  assert(oriB != 0.0); // SELF_CHECK
  // Flip the sign if there are odd number of swaps.
  if ((swaps % 2) != 0) oriB = -oriB;
  return oriB;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// orient4d_s()    4d orientation test with symbolic perturbation.           //
//                                                                           //
// Given four lifted points pa', pb', pc', and pd' in R^4,test if the lifted //
// point pe' in R^4 lies below or above the hyperplane passing through the   //
// four points pa', pb', pc', and pd'.                                       //
//                                                                           //
// Here we assume that the 3d orientation of the point sequence {pa, pb, pc, //
// pd} is positive (NOT zero), i.e., pd lies above the plane passing through //
// the points pa, pb, and pc. Otherwise, the returned sign is flipped.       //
//                                                                           //
// Return a positive value (> 0) if pe' lies below, a negative value (< 0)   //
// if pe' lies above the hyperplane, the returned value should not be zero.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::orient4d_s(REAL* pa, REAL* pb, REAL* pc, REAL* pd, REAL* pe,
                            REAL aheight, REAL bheight, REAL cheight, 
                            REAL dheight, REAL eheight)
{
  REAL sign;

  sign = orient4d(pa, pb, pc, pd, pe, 
                  aheight, bheight, cheight, dheight, eheight);
  if (sign != 0.0) {
    return sign;
  }

  // Symbolic perturbation.
  point pt[5], swappt;
  REAL oriA, oriB;
  int swaps, count;
  int n, i;

  pt[0] = pa;
  pt[1] = pb;
  pt[2] = pc;
  pt[3] = pd;
  pt[4] = pe;
  
  // Sort the five points such that their indices are in the increasing
  //   order. An optimized bubble sort algorithm is used, i.e., it has
  //   the worst case O(n^2) runtime, but it is usually much faster.
  swaps = 0; // Record the total number of swaps.
  n = 5;
  do {
    count = 0;
    n = n - 1;
    for (i = 0; i < n; i++) {
      if (pointmark(pt[i]) > pointmark(pt[i+1])) {
        swappt = pt[i]; pt[i] = pt[i+1]; pt[i+1] = swappt;
        count++;
      }
    }
    swaps += count;
  } while (count > 0); // Continue if some points are swapped.

  oriA = orient3d(pt[1], pt[2], pt[3], pt[4]);
  if (oriA != 0.0) {
    // Flip the sign if there are odd number of swaps.
    if ((swaps % 2) != 0) oriA = -oriA;
    return oriA;
  }
  
  oriB = -orient3d(pt[0], pt[2], pt[3], pt[4]);
  assert(oriB != 0.0); // SELF_CHECK
  // Flip the sign if there are odd number of swaps.
  if ((swaps % 2) != 0) oriB = -oriB;
  return oriB;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_edge_test()    Triangle-edge intersection test.                       //
//                                                                           //
// This routine takes a triangle T (with vertices A, B, C) and an edge E (P, //
// Q) in 3D, and tests if they intersect each other.                         //
//                                                                           //
// If the point 'R' is not NULL, it lies strictly above the plane defined by //
// A, B, C. It is used in test when T and E are coplanar.                    //
//                                                                           //
// If T and E intersect each other, they may intersect in different ways. If //
// 'level' > 0, their intersection type will be reported 'types' and 'pos'.  //
//                                                                           //
// The return value indicates one of the following cases:                    //
//   - 0, T and E are disjoint.                                              //
//   - 1, T and E intersect each other.                                      //
//   - 2, T and E are not coplanar. They intersect at a single point.        //
//   - 4, T and E are coplanar. They intersect at a single point or a line   //
//        segment (if types[1] != DISJOINT).                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#define SETVECTOR3(V, a0, a1, a2) (V)[0] = (a0); (V)[1] = (a1); (V)[2] = (a2)

#define SWAP2(a0, a1, tmp) (tmp) = (a0); (a0) = (a1); (a1) = (tmp)

int tetgenmesh::tri_edge_2d(point A, point B, point C, point P, point Q, 
                            point R, int level, int *types, int *pos)
{
  point U[3], V[3];  // The permuted vectors of points.
  int pu[3], pv[3];  // The original positions of points.
  REAL abovept[3];
  REAL sA, sB, sC;
  REAL s1, s2, s3, s4;
  int z1;

  if (R == NULL) {
    // Calculate a lift point.
    if (1) {
      REAL n[3], len;
      // Calculate a lift point, saved in dummypoint.
      facenormal(A, B, C, n, 1, NULL);
      len = sqrt(dot(n, n));
      if (len != 0) {
        n[0] /= len;
        n[1] /= len;
        n[2] /= len;
        len = distance(A, B);
        len += distance(B, C);
        len += distance(C, A);
        len /= 3.0;
        R = abovept; //dummypoint;
        R[0] = A[0] + len * n[0];
        R[1] = A[1] + len * n[1];
        R[2] = A[2] + len * n[2];
      } else {
        // The triangle [A,B,C] is (nearly) degenerate, i.e., it is (close)
        //   to a line.  We need a line-line intersection test.
        //assert(0);
        // !!! A non-save return value.!!!
        return 0;  // DISJOINT
      }
    }
  }

  // Test A's, B's, and C's orientations wrt plane PQR. 
  sA = orient3d(P, Q, R, A);
  sB = orient3d(P, Q, R, B);
  sC = orient3d(P, Q, R, C);


  if (sA < 0) {
    if (sB < 0) {
      if (sC < 0) { // (---).
        return 0; 
      } else {
        if (sC > 0) { // (--+).
          // All points are in the right positions.
          SETVECTOR3(U, A, B, C);  // I3
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 0, 1, 2);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else { // (--0).
          SETVECTOR3(U, A, B, C);  // I3
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 0, 1, 2);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        }
      }
    } else { 
      if (sB > 0) {
        if (sC < 0) { // (-+-).
          SETVECTOR3(U, C, A, B);  // PT = ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 2, 0, 1);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else {
          if (sC > 0) { // (-++).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 0;
          } else { // (-+0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          }
        }
      } else {
        if (sC < 0) { // (-0-).
          SETVECTOR3(U, C, A, B);  // PT = ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 2, 0, 1);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        } else {
          if (sC > 0) { // (-0+).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 2;
          } else { // (-00).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 3; 
          }
        }
      }
    }
  } else {
    if (sA > 0) {
      if (sB < 0) {
        if (sC < 0) { // (+--).
          SETVECTOR3(U, B, C, A);  // PT = ST x ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 1, 2, 0);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else {
          if (sC > 0) { // (+-+).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 0;
          } else { // (+-0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 2;
          }
        }
      } else { 
        if (sB > 0) {
          if (sC < 0) { // (++-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 0;
          } else {
            if (sC > 0) { // (+++).
              return 0; 
            } else { // (++0).
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 0, 1, 2);
              SETVECTOR3(pv, 1, 0, 2);
              z1 = 1; 
            }
          }
        } else { // (+0#)
          if (sC < 0) { // (+0-).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          } else {
            if (sC > 0) { // (+0+).
              SETVECTOR3(U, C, A, B);  // PT = ST
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 2, 0, 1);
              SETVECTOR3(pv, 1, 0, 2);
              z1 = 1;
            } else { // (+00).
              SETVECTOR3(U, B, C, A);  // PT = ST x ST
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 1, 2, 0);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            }
          }
        }
      }
    } else { 
      if (sB < 0) {
        if (sC < 0) { // (0--).
          SETVECTOR3(U, B, C, A);  // PT = ST x ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 1, 2, 0);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        } else {
          if (sC > 0) { // (0-+).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          } else { // (0-0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 3; 
          }
        }
      } else { 
        if (sB > 0) {
          if (sC < 0) { // (0+-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 2;
          } else {
            if (sC > 0) { // (0++).
              SETVECTOR3(U, B, C, A);  // PT = ST x ST
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 1, 2, 0);
              SETVECTOR3(pv, 1, 0, 2);
              z1 = 1;
            } else { // (0+0).
              SETVECTOR3(U, C, A, B);  // PT = ST
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 2, 0, 1);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            }
          }
        } else { // (00#)
          if (sC < 0) { // (00-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 3; 
          } else {
            if (sC > 0) { // (00+).
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 0, 1, 2);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            } else { // (000)
              // Not possible unless ABC is degenerate.
              // Avoiding compiler warnings.
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 0, 1, 2);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 4;
            }
          }
        }
      }
    }
  }

  s1 = orient3d(U[0], U[2], R, V[1]);  // A, C, R, Q
  s2 = orient3d(U[1], U[2], R, V[0]);  // B, C, R, P

  if (s1 > 0) {
    return 0;
  }
  if (s2 < 0) {
    return 0;
  }

  if (level == 0) {
    return 1;  // They are intersected.
  }

  assert(z1 != 4); // SELF_CHECK

  if (z1 == 1) {
    if (s1 == 0) {  // (0###)
      // C = Q.
      types[0] = (int) SHAREVERT;
      pos[0] = pu[2]; // C
      pos[1] = pv[1]; // Q
      types[1] = (int) DISJOINT;
    } else {
      if (s2 == 0) { // (#0##)
        // C = P.
        types[0] = (int) SHAREVERT;
        pos[0] = pu[2]; // C
        pos[1] = pv[0]; // P
        types[1] = (int) DISJOINT;
      } else { // (-+##)
        // C in [P, Q].
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[2]; // C
        pos[1] = pv[0]; // [P, Q]
        types[1] = (int) DISJOINT;
      }
    }
    return 4;
  }

  s3 = orient3d(U[0], U[2], R, V[0]);  // A, C, R, P
  s4 = orient3d(U[1], U[2], R, V[1]);  // B, C, R, Q
      
  if (z1 == 0) {  // (tritri-03)
    if (s1 < 0) {
      if (s3 > 0) {
        assert(s2 > 0); // SELF_CHECK
        if (s4 > 0) {
          // [P, Q] overlaps [k, l] (-+++).
          types[0] = (int) ACROSSEDGE;
          pos[0] = pu[2]; // [C, A]
          pos[1] = pv[0]; // [P, Q]
          types[1] = (int) TOUCHFACE;
          pos[2] = 3;     // [A, B, C]
          pos[3] = pv[1]; // Q
        } else {
          if (s4 == 0) {
            // Q = l, [P, Q] contains [k, l] (-++0).
            types[0] = (int) ACROSSEDGE;
            pos[0] = pu[2]; // [C, A]
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) TOUCHEDGE;
            pos[2] = pu[1]; // [B, C]
            pos[3] = pv[1]; // Q
          } else { // s4 < 0
            // [P, Q] contains [k, l] (-++-).
            types[0] = (int) ACROSSEDGE;
            pos[0] = pu[2]; // [C, A]
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) ACROSSEDGE;
            pos[2] = pu[1]; // [B, C]
            pos[3] = pv[0]; // [P, Q]
          }
        }
      } else {
        if (s3 == 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // P = k, [P, Q] in [k, l] (-+0+).
            types[0] = (int) TOUCHEDGE;
            pos[0] = pu[2]; // [C, A]
            pos[1] = pv[0]; // P
            types[1] = (int) TOUCHFACE;
            pos[2] = 3;     // [A, B, C]
            pos[3] = pv[1]; // Q
          } else {
            if (s4 == 0) {
              // [P, Q] = [k, l] (-+00).
              types[0] = (int) TOUCHEDGE;
              pos[0] = pu[2]; // [C, A]
              pos[1] = pv[0]; // P
              types[1] = (int) TOUCHEDGE;
              pos[2] = pu[1]; // [B, C]
              pos[3] = pv[1]; // Q
            } else {
              // P = k, [P, Q] contains [k, l] (-+0-).
              types[0] = (int) TOUCHEDGE;
              pos[0] = pu[2]; // [C, A]
              pos[1] = pv[0]; // P
              types[1] = (int) ACROSSEDGE;
              pos[2] = pu[1]; // [B, C]
              pos[3] = pv[0]; // [P, Q]
            }
          }
        } else { // s3 < 0
          if (s2 > 0) {
            if (s4 > 0) {
              // [P, Q] in [k, l] (-+-+).
              types[0] = (int) TOUCHFACE;
              pos[0] = 3;     // [A, B, C]
              pos[1] = pv[0]; // P
              types[1] = (int) TOUCHFACE;
              pos[2] = 3;     // [A, B, C]
              pos[3] = pv[1]; // Q
            } else {
              if (s4 == 0) {
                // Q = l, [P, Q] in [k, l] (-+-0).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // [A, B, C]
                pos[1] = pv[0]; // P
                types[1] = (int) TOUCHEDGE;
                pos[2] = pu[1]; // [B, C]
                pos[3] = pv[1]; // Q
              } else { // s4 < 0
                // [P, Q] overlaps [k, l] (-+--).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // [A, B, C]
                pos[1] = pv[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = pu[1]; // [B, C]
                pos[3] = pv[0]; // [P, Q]
              }
            }
          } else { // s2 == 0
            // P = l (#0##).
            types[0] = (int) TOUCHEDGE;
            pos[0] = pu[1]; // [B, C]
            pos[1] = pv[0]; // P
            types[1] = (int) DISJOINT;
          }
        }
      }
    } else { // s1 == 0
      // Q = k (0####)
      types[0] = (int) TOUCHEDGE;
      pos[0] = pu[2]; // [C, A]
      pos[1] = pv[1]; // Q
      types[1] = (int) DISJOINT;
    }
  } else if (z1 == 2) {  // (tritri-23)
    if (s1 < 0) {
      if (s3 > 0) {
        assert(s2 > 0); // SELF_CHECK
        if (s4 > 0) {
          // [P, Q] overlaps [A, l] (-+++).
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[0]; // A
          pos[1] = pv[0]; // [P, Q]
          types[1] = (int) TOUCHFACE;
          pos[2] = 3;     // [A, B, C]
          pos[3] = pv[1]; // Q
        } else {
          if (s4 == 0) {
            // Q = l, [P, Q] contains [A, l] (-++0).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) TOUCHEDGE;
            pos[2] = pu[1]; // [B, C]
            pos[3] = pv[1]; // Q
          } else { // s4 < 0
            // [P, Q] contains [A, l] (-++-).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) ACROSSEDGE;
            pos[2] = pu[1]; // [B, C]
            pos[3] = pv[0]; // [P, Q]
          }
        }
      } else {
        if (s3 == 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // P = A, [P, Q] in [A, l] (-+0+).
            types[0] = (int) SHAREVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // P
            types[1] = (int) TOUCHFACE;
            pos[2] = 3;     // [A, B, C]
            pos[3] = pv[1]; // Q
          } else {
            if (s4 == 0) {
              // [P, Q] = [A, l] (-+00).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pv[0]; // P
              types[1] = (int) TOUCHEDGE;
              pos[2] = pu[1]; // [B, C]
              pos[3] = pv[1]; // Q
            } else { // s4 < 0
              // Q = l, [P, Q] in [A, l] (-+0-).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pv[0]; // P
              types[1] = (int) ACROSSEDGE;
              pos[2] = pu[1]; // [B, C]
              pos[3] = pv[0]; // [P, Q]
            }
          }
        } else { // s3 < 0
          if (s2 > 0) {
            if (s4 > 0) {
              // [P, Q] in [A, l] (-+-+).
              types[0] = (int) TOUCHFACE;
              pos[0] = 3;     // [A, B, C]
              pos[1] = pv[0]; // P
              types[0] = (int) TOUCHFACE;
              pos[0] = 3;     // [A, B, C]
              pos[1] = pv[1]; // Q
            } else {
              if (s4 == 0) {
                // Q = l, [P, Q] in [A, l] (-+-0).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // [A, B, C]
                pos[1] = pv[0]; // P
                types[0] = (int) TOUCHEDGE;
                pos[0] = pu[1]; // [B, C]
                pos[1] = pv[1]; // Q
              } else { // s4 < 0
                // [P, Q] overlaps [A, l] (-+--).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // [A, B, C]
                pos[1] = pv[0]; // P
                types[0] = (int) ACROSSEDGE;
                pos[0] = pu[1]; // [B, C]
                pos[1] = pv[0]; // [P, Q]
              }
            }
          } else { // s2 == 0
            // P = l (#0##).
            types[0] = (int) TOUCHEDGE;
            pos[0] = pu[1]; // [B, C]
            pos[1] = pv[0]; // P
            types[1] = (int) DISJOINT;
          }
        }
      }
    } else { // s1 == 0
      // Q = A (0###).
      types[0] = (int) SHAREVERT;
      pos[0] = pu[0]; // A
      pos[1] = pv[1]; // Q
      types[1] = (int) DISJOINT;
    }
  } else if (z1 == 3) {  // (tritri-33)
    if (s1 < 0) {
      if (s3 > 0) {
        assert(s2 > 0); // SELF_CHECK
        if (s4 > 0) {
          // [P, Q] overlaps [A, B] (-+++).
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[0]; // A
          pos[1] = pv[0]; // [P, Q]
          types[1] = (int) TOUCHEDGE;
          pos[2] = pu[0]; // [A, B]
          pos[3] = pv[1]; // Q
        } else {
          if (s4 == 0) {
            // Q = B, [P, Q] contains [A, B] (-++0).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) SHAREVERT;
            pos[2] = pu[1]; // B
            pos[3] = pv[1]; // Q
          } else { // s4 < 0
            // [P, Q] contains [A, B] (-++-).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) ACROSSVERT;
            pos[2] = pu[1]; // B
            pos[3] = pv[0]; // [P, Q]
          }
        }
      } else {
        if (s3 == 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // P = A, [P, Q] in [A, B] (-+0+).
            types[0] = (int) SHAREVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // P
            types[1] = (int) TOUCHEDGE;
            pos[2] = pu[0]; // [A, B]
            pos[3] = pv[1]; // Q
          } else {
            if (s4 == 0) {
              // [P, Q] = [A, B] (-+00).
              types[0] = (int) SHAREEDGE;
              pos[0] = pu[0]; // [A, B]
              pos[1] = pv[0]; // [P, Q]
              types[1] = (int) DISJOINT;
            } else { // s4 < 0
              // P= A, [P, Q] in [A, B] (-+0-).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pv[0]; // P
              types[1] = (int) ACROSSVERT;
              pos[2] = pu[1]; // B
              pos[3] = pv[0]; // [P, Q]
            }
          }
        } else { // s3 < 0
          if (s2 > 0) {
            if (s4 > 0) {
              // [P, Q] in [A, B] (-+-+).
              types[0] = (int) TOUCHEDGE;
              pos[0] = pu[0]; // [A, B]
              pos[1] = pv[0]; // P
              types[1] = (int) TOUCHEDGE;
              pos[2] = pu[0]; // [A, B]
              pos[3] = pv[1]; // Q
            } else {
              if (s4 == 0) {
                // Q = B, [P, Q] in [A, B] (-+-0).
                types[0] = (int) TOUCHEDGE;
                pos[0] = pu[0]; // [A, B]
                pos[1] = pv[0]; // P
                types[1] = (int) SHAREVERT;
                pos[2] = pu[1]; // B
                pos[3] = pv[1]; // Q
              } else { // s4 < 0
                // [P, Q] overlaps [A, B] (-+--).
                types[0] = (int) TOUCHEDGE;
                pos[0] = pu[0]; // [A, B]
                pos[1] = pv[0]; // P
                types[1] = (int) ACROSSVERT;
                pos[2] = pu[1]; // B
                pos[3] = pv[0]; // [P, Q]
              }
            }
          } else { // s2 == 0
            // P = B (#0##).
            types[0] = (int) SHAREVERT;
            pos[0] = pu[1]; // B
            pos[1] = pv[0]; // P
            types[1] = (int) DISJOINT;
          }
        }
      }
    } else { // s1 == 0
      // Q = A (0###).
      types[0] = (int) SHAREVERT;
      pos[0] = pu[0]; // A
      pos[1] = pv[1]; // Q
      types[1] = (int) DISJOINT;
    }
  }

  return 4;
}

int tetgenmesh::tri_edge_tail(point A,point B,point C,point P,point Q,point R,
                              REAL sP,REAL sQ,int level,int *types,int *pos)
{
  point U[3], V[3]; //, Ptmp;
  int pu[3], pv[3]; //, itmp;
  REAL s1, s2, s3;
  int z1;


  if (sP < 0) {
    if (sQ < 0) { // (--) disjoint
      return 0;
    } else {
      if (sQ > 0) { // (-+)
        SETVECTOR3(U, A, B, C);
        SETVECTOR3(V, P, Q, R);
        SETVECTOR3(pu, 0, 1, 2);
        SETVECTOR3(pv, 0, 1, 2);
        z1 = 0;
      } else { // (-0)
        SETVECTOR3(U, A, B, C);
        SETVECTOR3(V, P, Q, R);
        SETVECTOR3(pu, 0, 1, 2);
        SETVECTOR3(pv, 0, 1, 2);
        z1 = 1;
      }
    }
  } else {
    if (sP > 0) { // (+-)
      if (sQ < 0) {
        SETVECTOR3(U, A, B, C);
        SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
        SETVECTOR3(pu, 0, 1, 2);
        SETVECTOR3(pv, 1, 0, 2);
        z1 = 0;
      } else {
        if (sQ > 0) { // (++) disjoint
          return 0;
        } else { // (+0)
          SETVECTOR3(U, B, A, C); // A and B are flipped.
          SETVECTOR3(V, P, Q, R);
          SETVECTOR3(pu, 1, 0, 2);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        }
      }
    } else { // sP == 0
      if (sQ < 0) { // (0-)
        SETVECTOR3(U, A, B, C);
        SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
        SETVECTOR3(pu, 0, 1, 2);
        SETVECTOR3(pv, 1, 0, 2);
        z1 = 1;
      } else {
        if (sQ > 0) { // (0+)
          SETVECTOR3(U, B, A, C);  // A and B are flipped.
          SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
          SETVECTOR3(pu, 1, 0, 2);
          SETVECTOR3(pv, 1, 0, 2);
          z1 = 1;
        } else { // (00)
          // A, B, C, P, and Q are coplanar.
          z1 = 2;
        }
      }
    }
  }

  if (z1 == 2) {
    // The triangle and the edge are coplanar.
    return tri_edge_2d(A, B, C, P, Q, R, level, types, pos);
  }

  s1 = orient3d(U[0], U[1], V[0], V[1]);
  if (s1 < 0) {
    return 0;
  }

  s2 = orient3d(U[1], U[2], V[0], V[1]);
  if (s2 < 0) {
    return 0;
  }

  s3 = orient3d(U[2], U[0], V[0], V[1]);
  if (s3 < 0) {
    return 0;
  }

  if (level == 0) {
    return 1;  // The are intersected.
  }

  types[1] = (int) DISJOINT; // No second intersection point.

  if (z1 == 0) {
    if (s1 > 0) {
      if (s2 > 0) {
        if (s3 > 0) { // (+++)
          // [P, Q] passes interior of [A, B, C].
          types[0] = (int) ACROSSFACE;
          pos[0] = 3;  // interior of [A, B, C]
          pos[1] = 0;  // [P, Q]
        } else { // s3 == 0 (++0)
          // [P, Q] intersects [C, A].
          types[0] = (int) ACROSSEDGE;
          pos[0] = pu[2];  // [C, A]
          pos[1] = 0;  // [P, Q]
        }
      } else { // s2 == 0
        if (s3 > 0) { // (+0+)
          // [P, Q] intersects [B, C].
          types[0] = (int) ACROSSEDGE;
          pos[0] = pu[1];  // [B, C]
          pos[1] = 0;  // [P, Q]
        } else { // s3 == 0 (+00)
          // [P, Q] passes C.
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[2];  // C
          pos[1] = 0;  // [P, Q]
        }
      }
    } else { // s1 == 0
      if (s2 > 0) {
        if (s3 > 0) { // (0++)
          // [P, Q] intersects [A, B].
          types[0] = (int) ACROSSEDGE;
          pos[0] = pu[0];  // [A, B]
          pos[1] = 0;  // [P, Q]
        } else { // s3 == 0 (0+0)
          // [P, Q] passes A.
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[0];  // A
          pos[1] = 0;  // [P, Q]
        }
      } else { // s2 == 0
        if (s3 > 0) { // (00+)
          // [P, Q] passes B.
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[1];  // B
          pos[1] = 0;  // [P, Q]
        } else { // s3 == 0 (000)
          // Impossible.
          assert(0);
        }
      }
    }
  } else { // z1 == 1
    if (s1 > 0) {
      if (s2 > 0) {
        if (s3 > 0) { // (+++)
          // Q lies in [A, B, C].
          types[0] = (int) TOUCHFACE;
          pos[0] = 0; // [A, B, C]
          pos[1] = pv[1]; // Q
        } else { // s3 == 0 (++0)
          // Q lies on [C, A].
          types[0] = (int) TOUCHEDGE;
          pos[0] = pu[2]; // [C, A]
          pos[1] = pv[1]; // Q
        }
      } else { // s2 == 0
        if (s3 > 0) { // (+0+)
          // Q lies on [B, C].
          types[0] = (int) TOUCHEDGE;
          pos[0] = pu[1]; // [B, C]
          pos[1] = pv[1]; // Q
        } else { // s3 == 0 (+00)
          // Q = C.
          types[0] = (int) SHAREVERT;
          pos[0] = pu[2]; // C
          pos[1] = pv[1]; // Q
        }
      }
    } else { // s1 == 0
      if (s2 > 0) {
        if (s3 > 0) { // (0++)
          // Q lies on [A, B].
          types[0] = (int) TOUCHEDGE;
          pos[0] = pu[0]; // [A, B]
          pos[1] = pv[1]; // Q
        } else { // s3 == 0 (0+0)
          // Q = A.
          types[0] = (int) SHAREVERT;
          pos[0] = pu[0]; // A
          pos[1] = pv[1]; // Q
        }
      } else { // s2 == 0
        if (s3 > 0) { // (00+)
          // Q = B.
          types[0] = (int) SHAREVERT;
          pos[0] = pu[1]; // B
          pos[1] = pv[1]; // Q
        } else { // s3 == 0 (000)
          // Impossible.
          assert(0);
        }
      }
    }
  }

  // T and E intersect in a single point.
  return 2;
}

int tetgenmesh::tri_edge_test(point A, point B, point C, point P, point Q, 
                              point R, int level, int *types, int *pos)
{
  REAL sP, sQ;

  // Test the locations of P and Q with respect to ABC.
  sP = orient3d(A, B, C, P);
  sQ = orient3d(A, B, C, Q);

  return tri_edge_tail(A, B, C, P, Q, R, sP, sQ, level, types, pos);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_tri_inter()    Test whether two triangle (abc) and (opq) are          //
//                    intersecting or not.                                   //
//                                                                           //
// Return 0 if they are disjoint. Otherwise, return 1. 'type' returns one of //
// the four cases: SHAREVERTEX, SHAREEDGE, SHAREFACE, and INTERSECT.         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::tri_edge_inter_tail(REAL* A, REAL* B, REAL* C,  REAL* P, 
                                    REAL* Q, REAL s_p, REAL s_q)
{
  int types[2], pos[4];
  int ni;  // =0, 2, 4

  ni = tri_edge_tail(A, B, C, P, Q, NULL, s_p, s_q, 1, types, pos);

  if (ni > 0) {
    if (ni == 2) {
      // Get the intersection type.
      if (types[0] == (int) SHAREVERT) {
        return (int) SHAREVERT;
      } else {
        return (int) INTERSECT;
      }
    } else if (ni == 4) { 
      // There may be two intersections.
      if (types[0] == (int) SHAREVERT) {
        if (types[1] == (int) DISJOINT) {
          return (int) SHAREVERT;
        } else {
          assert(types[1] != (int) SHAREVERT);
          return (int) INTERSECT;
        }
      } else {
        if (types[0] == (int) SHAREEDGE) {
          return (int) SHAREEDGE;
        } else {
          return (int) INTERSECT;
        }
      }
    } else {
      assert(0);
    }
  }

  return (int) DISJOINT;
}

int tetgenmesh::tri_tri_inter(REAL* A,REAL* B,REAL* C,REAL* O,REAL* P,REAL* Q)
{
  REAL s_o, s_p, s_q;
  REAL s_a, s_b, s_c;

  s_o = orient3d(A, B, C, O);
  s_p = orient3d(A, B, C, P);
  s_q = orient3d(A, B, C, Q);
  if ((s_o * s_p > 0.0) && (s_o * s_q > 0.0)) {
    // o, p, q are all in the same halfspace of ABC.
    return 0; // DISJOINT;
  }

  s_a = orient3d(O, P, Q, A);
  s_b = orient3d(O, P, Q, B);
  s_c = orient3d(O, P, Q, C);
  if ((s_a * s_b > 0.0) && (s_a * s_c > 0.0)) {
    // a, b, c are all in the same halfspace of OPQ.
    return 0; // DISJOINT;
  }

  int abcop, abcpq, abcqo;
  int shareedge = 0;

  abcop = tri_edge_inter_tail(A, B, C, O, P, s_o, s_p);
  if (abcop == (int) INTERSECT) {
    return (int) INTERSECT;
  } else if (abcop == (int) SHAREEDGE) {
    shareedge++;
  }
  abcpq = tri_edge_inter_tail(A, B, C, P, Q, s_p, s_q);
  if (abcpq == (int) INTERSECT) {
    return (int) INTERSECT;
  } else if (abcpq == (int) SHAREEDGE) {
    shareedge++;
  }
  abcqo = tri_edge_inter_tail(A, B, C, Q, O, s_q, s_o);
  if (abcqo == (int) INTERSECT) {
    return (int) INTERSECT;
  } else if (abcqo == (int) SHAREEDGE) {
    shareedge++;
  }
  if (shareedge == 3) {
    // opq are coincident with abc.
    return (int) SHAREFACE;
  }

  // It is only possible either no share edge or one.
  assert(shareedge == 0 || shareedge == 1);

  // Continue to detect whether opq and abc are intersecting or not.
  int opqab, opqbc, opqca;

  opqab = tri_edge_inter_tail(O, P, Q, A, B, s_a, s_b);
  if (opqab == (int) INTERSECT) {
    return (int) INTERSECT;
  }
  opqbc = tri_edge_inter_tail(O, P, Q, B, C, s_b, s_c);
  if (opqbc == (int) INTERSECT) {
    return (int) INTERSECT;
  }
  opqca = tri_edge_inter_tail(O, P, Q, C, A, s_c, s_a);
  if (opqca == (int) INTERSECT) {
    return (int) INTERSECT;
  }

  // At this point, two triangles are not intersecting and not coincident.
  //   They may be share an edge, or share a vertex, or disjoint.
  if (abcop == (int) SHAREEDGE) {
    assert((abcpq == (int) SHAREVERT) && (abcqo == (int) SHAREVERT));
    // op is coincident with an edge of abc.
    return (int) SHAREEDGE;
  }
  if (abcpq == (int) SHAREEDGE) {
    assert((abcop == (int) SHAREVERT) && (abcqo == (int) SHAREVERT));
    // pq is coincident with an edge of abc.
    return (int) SHAREEDGE;
  }
  if (abcqo == (int) SHAREEDGE) {
    assert((abcop == (int) SHAREVERT) && (abcpq == (int) SHAREVERT));
    // qo is coincident with an edge of abc.
    return (int) SHAREEDGE;
  }

  // They may share a vertex or disjoint.
  if (abcop == (int) SHAREVERT) {
    // o or p is coincident with a vertex of abc.
    if (abcpq == (int) SHAREVERT) {
      // p is the coincident vertex.
      assert(abcqo != (int) SHAREVERT);
    } else {
      // o is the coincident vertex.
      assert(abcqo == (int) SHAREVERT);
    }
    return (int) SHAREVERT;
  }
  if (abcpq == (int) SHAREVERT) {
    // q is the coincident vertex.
    assert(abcqo == (int) SHAREVERT);
    return (int) SHAREVERT;
  }

  // They are disjoint.
  return (int) DISJOINT;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lu_decmp()    Compute the LU decomposition of a matrix.                   //
//                                                                           //
// Compute the LU decomposition of a (non-singular) square matrix A using    //
// partial pivoting and implicit row exchanges.  The result is:              //
//     A = P * L * U,                                                        //
// where P is a permutation matrix, L is unit lower triangular, and U is     //
// upper triangular.  The factored form of A is used in combination with     //
// 'lu_solve()' to solve linear equations: Ax = b, or invert a matrix.       //
//                                                                           //
// The inputs are a square matrix 'lu[N..n+N-1][N..n+N-1]', it's size is 'n'.//
// On output, 'lu' is replaced by the LU decomposition of a rowwise permuta- //
// tion of itself, 'ps[N..n+N-1]' is an output vector that records the row   //
// permutation effected by the partial pivoting, effectively,  'ps' array    //
// tells the user what the permutation matrix P is; 'd' is output as +1/-1   //
// depending on whether the number of row interchanges was even or odd,      //
// respectively.                                                             //
//                                                                           //
// Return true if the LU decomposition is successfully computed, otherwise,  //
// return false in case that A is a singular matrix.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::lu_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N)
{
  REAL scales[4];
  REAL pivot, biggest, mult, tempf;
  int pivotindex = 0;
  int i, j, k;

  *d = 1.0;                                      // No row interchanges yet.

  for (i = N; i < n + N; i++) {                             // For each row.
    // Find the largest element in each row for row equilibration
    biggest = 0.0;
    for (j = N; j < n + N; j++)
      if (biggest < (tempf = fabs(lu[i][j])))
        biggest  = tempf;
    if (biggest != 0.0)
      scales[i] = 1.0 / biggest;
    else {
      scales[i] = 0.0;
      return false;                            // Zero row: singular matrix.
    }
    ps[i] = i;                                 // Initialize pivot sequence.
  }

  for (k = N; k < n + N - 1; k++) {                      // For each column.
    // Find the largest element in each column to pivot around.
    biggest = 0.0;
    for (i = k; i < n + N; i++) {
      if (biggest < (tempf = fabs(lu[ps[i]][k]) * scales[ps[i]])) {
        biggest = tempf;
        pivotindex = i;
      }
    }
    if (biggest == 0.0) {
      return false;                         // Zero column: singular matrix.
    }
    if (pivotindex != k) {                         // Update pivot sequence.
      j = ps[k];
      ps[k] = ps[pivotindex];
      ps[pivotindex] = j;
      *d = -(*d);                          // ...and change the parity of d.
    }

    // Pivot, eliminating an extra variable  each time
    pivot = lu[ps[k]][k];
    for (i = k + 1; i < n + N; i++) {
      lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;
      if (mult != 0.0) {
        for (j = k + 1; j < n + N; j++)
          lu[ps[i]][j] -= mult * lu[ps[k]][j];
      }
    }
  }

  // (lu[ps[n + N - 1]][n + N - 1] == 0.0) ==> A is singular.
  return lu[ps[n + N - 1]][n + N - 1] != 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lu_solve()    Solves the linear equation:  Ax = b,  after the matrix A    //
//               has been decomposed into the lower and upper triangular     //
//               matrices L and U, where A = LU.                             //
//                                                                           //
// 'lu[N..n+N-1][N..n+N-1]' is input, not as the matrix 'A' but rather as    //
// its LU decomposition, computed by the routine 'lu_decmp'; 'ps[N..n+N-1]'  //
// is input as the permutation vector returned by 'lu_decmp';  'b[N..n+N-1]' //
// is input as the right-hand side vector, and returns with the solution     //
// vector. 'lu', 'n', and 'ps' are not modified by this routine and can be   //
// left in place for successive calls with different right-hand sides 'b'.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::lu_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N)
{
  int i, j;
  REAL X[4], dot;

  for (i = N; i < n + N; i++) X[i] = 0.0;

  // Vector reduction using U triangular matrix.
  for (i = N; i < n + N; i++) {
    dot = 0.0;
    for (j = N; j < i + N; j++)
      dot += lu[ps[i]][j] * X[j];
    X[i] = b[ps[i]] - dot;
  }

  // Back substitution, in L triangular matrix.
  for (i = n + N - 1; i >= N; i--) {
    dot = 0.0;
    for (j = i + 1; j < n + N; j++)
      dot += lu[ps[i]][j] * X[j];
    X[i] = (X[i] - dot) / lu[ps[i]][i];
  }

  for (i = N; i < n + N; i++) b[i] = X[i];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// incircle3d()    3D in-circle test.                                        //
//                                                                           //
// Return a negative value if pd is inside the circumcircle of the triangle  //
// pa, pb, and pc.                                                           //
//                                                                           //
// IMPORTANT: It assumes that [a,b] is the common edge, i.e., the two input  //
// triangles are [a,b,c] and [b,a,d].                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::incircle3d(point pa, point pb, point pc, point pd)
{
  REAL area2[2], n1[3], n2[3], c[3];
  REAL sign, r, d;

  // Calculate the areas of the two triangles [a, b, c] and [b, a, d].
  facenormal(pa, pb, pc, n1, 1, NULL);
  area2[0] = dot(n1, n1);
  facenormal(pb, pa, pd, n2, 1, NULL);
  area2[1] = dot(n2, n2);

  if (area2[0] > area2[1]) {
    // Choose [a, b, c] as the base triangle.
    circumsphere(pa, pb, pc, NULL, c, &r);
    d = distance(c, pd);
  } else {
    // Choose [b, a, d] as the base triangle.
    if (area2[1] > 0) {
      circumsphere(pb, pa, pd, NULL, c, &r);
      d = distance(c, pc);
    } else {
      // The four points are collinear. This case only happens on the boundary.
      return 0; // Return "not inside".
    }
  }

  sign = d - r;
  if (fabs(sign) / r < b->epsilon) {
    sign = 0;
  }

  return sign;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// facenormal()    Calculate the normal of the face.                         //
//                                                                           //
// The normal of the face abc can be calculated by the cross product of 2 of //
// its 3 edge vectors.  A better choice of two edge vectors will reduce the  //
// numerical error during the calculation.  Burdakov proved that the optimal //
// basis problem is equivalent to the minimum spanning tree problem with the //
// edge length be the functional, see Burdakov, "A greedy algorithm for the  //
// optimal basis problem", BIT 37:3 (1997), 591-599. If 'pivot' > 0, the two //
// short edges in abc are chosen for the calculation.                        //
//                                                                           //
// If 'lav' is not NULL and if 'pivot' is set, the average edge length of    //
// the edges of the face [a,b,c] is returned.                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::facenormal(point pa, point pb, point pc, REAL *n, int pivot,
                            REAL* lav)
{
  REAL v1[3], v2[3], v3[3], *pv1, *pv2;
  REAL L1, L2, L3;

  v1[0] = pb[0] - pa[0];  // edge vector v1: a->b
  v1[1] = pb[1] - pa[1];
  v1[2] = pb[2] - pa[2];
  v2[0] = pa[0] - pc[0];  // edge vector v2: c->a
  v2[1] = pa[1] - pc[1];
  v2[2] = pa[2] - pc[2];

  // Default, normal is calculated by: v1 x (-v2) (see Fig. fnormal).
  if (pivot > 0) {
    // Choose edge vectors by Burdakov's algorithm.
    v3[0] = pc[0] - pb[0];  // edge vector v3: b->c
    v3[1] = pc[1] - pb[1];
    v3[2] = pc[2] - pb[2];
    L1 = dot(v1, v1);
    L2 = dot(v2, v2);
    L3 = dot(v3, v3);
    // Sort the three edge lengths.
    if (L1 < L2) {
      if (L2 < L3) {
        pv1 = v1; pv2 = v2; // n = v1 x (-v2).
      } else {
        pv1 = v3; pv2 = v1; // n = v3 x (-v1).
      }
    } else {
      if (L1 < L3) {
        pv1 = v1; pv2 = v2; // n = v1 x (-v2).
      } else {
        pv1 = v2; pv2 = v3; // n = v2 x (-v3).
      }
    }
    if (lav) {
      // return the average edge length.
      *lav = (sqrt(L1) + sqrt(L2) + sqrt(L3)) / 3.0;
    }
  } else {
    pv1 = v1; pv2 = v2; // n = v1 x (-v2).
  }

  // Calculate the face normal.
  cross(pv1, pv2, n);
  // Inverse the direction;
  n[0] = -n[0];
  n[1] = -n[1];
  n[2] = -n[2];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// shortdistance()    Returns the shortest distance from point p to a line   //
//                    defined by two points e1 and e2.                       //
//                                                                           //
// First compute the projection length l_p of the vector v1 = p - e1 along   //
// the vector v2 = e2 - e1. Then Pythagoras' Theorem is used to compute the  //
// shortest distance.                                                        //
//                                                                           //
// This routine allows that p is collinear with the line. In this case, the  //
// return value is zero. The two points e1 and e2 should not be identical.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::shortdistance(REAL* p, REAL* e1, REAL* e2)
{
  REAL v1[3], v2[3];
  REAL len, l_p;

  v1[0] = e2[0] - e1[0];
  v1[1] = e2[1] - e1[1];
  v1[2] = e2[2] - e1[2];
  v2[0] = p[0] - e1[0];
  v2[1] = p[1] - e1[1];
  v2[2] = p[2] - e1[2];

  len = sqrt(dot(v1, v1));
  assert(len != 0.0);

  v1[0] /= len;
  v1[1] /= len;
  v1[2] /= len;
  l_p = dot(v1, v2);

  return sqrt(dot(v2, v2) - l_p * l_p);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// triarea()    Return the area of a triangle.                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::triarea(REAL* pa, REAL* pb, REAL* pc)
{
  REAL A[4][4];  

  // Compute the coefficient matrix A (3x3).
  A[0][0] = pb[0] - pa[0];
  A[0][1] = pb[1] - pa[1];
  A[0][2] = pb[2] - pa[2]; // vector V1 (pa->pb)
  A[1][0] = pc[0] - pa[0];
  A[1][1] = pc[1] - pa[1];
  A[1][2] = pc[2] - pa[2]; // vector V2 (pa->pc)

  cross(A[0], A[1], A[2]); // vector V3 (V1 X V2)

  return 0.5 * sqrt(dot(A[2], A[2])); // The area of [a,b,c].
}

REAL tetgenmesh::orient3dfast(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
{
  REAL adx, bdx, cdx;
  REAL ady, bdy, cdy;
  REAL adz, bdz, cdz;

  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];
  adz = pa[2] - pd[2];
  bdz = pb[2] - pd[2];
  cdz = pc[2] - pd[2];

  return adx * (bdy * cdz - bdz * cdy)
       + bdx * (cdy * adz - cdz * ady)
       + cdx * (ady * bdz - adz * bdy);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// interiorangle()    Return the interior angle (0 - 2 * PI) between vectors //
//                    o->p1 and o->p2.                                       //
//                                                                           //
// 'n' is the normal of the plane containing face (o, p1, p2).  The interior //
// angle is the total angle rotating from o->p1 around n to o->p2.  Exchange //
// the position of p1 and p2 will get the complement angle of the other one. //
// i.e., interiorangle(o, p1, p2) = 2 * PI - interiorangle(o, p2, p1).  Set  //
// 'n' be NULL if you only want the interior angle between 0 - PI.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::interiorangle(REAL* o, REAL* p1, REAL* p2, REAL* n)
{
  REAL v1[3], v2[3], np[3];
  REAL theta, costheta, lenlen;
  REAL ori, len1, len2;

  // Get the interior angle (0 - PI) between o->p1, and o->p2.
  v1[0] = p1[0] - o[0];
  v1[1] = p1[1] - o[1];
  v1[2] = p1[2] - o[2];
  v2[0] = p2[0] - o[0];
  v2[1] = p2[1] - o[1];
  v2[2] = p2[2] - o[2];
  len1 = sqrt(dot(v1, v1));
  len2 = sqrt(dot(v2, v2));
  lenlen = len1 * len2;
  assert(lenlen != 0.0);

  costheta = dot(v1, v2) / lenlen;
  if (costheta > 1.0) {
    costheta = 1.0; // Roundoff. 
  } else if (costheta < -1.0) {
    costheta = -1.0; // Roundoff. 
  }
  theta = acos(costheta);
  if (n != NULL) {
    // Get a point above the face (o, p1, p2);
    np[0] = o[0] + n[0];
    np[1] = o[1] + n[1];
    np[2] = o[2] + n[2];
    // Adjust theta (0 - 2 * PI).
    ori = orient3d(p1, o, np, p2);
    if (ori > 0.0) {
      theta = 2 * PI - theta;
    }
  }

  return theta;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// projpt2edge()    Return the projection point from a point to an edge.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::projpt2edge(REAL* p, REAL* e1, REAL* e2, REAL* prj)
{
  REAL v1[3], v2[3];
  REAL len, l_p;

  v1[0] = e2[0] - e1[0];
  v1[1] = e2[1] - e1[1];
  v1[2] = e2[2] - e1[2];
  v2[0] = p[0] - e1[0];
  v2[1] = p[1] - e1[1];
  v2[2] = p[2] - e1[2];

  len = sqrt(dot(v1, v1));
  assert(len != 0.0);
  v1[0] /= len;
  v1[1] /= len;
  v1[2] /= len;
  l_p = dot(v1, v2);

  prj[0] = e1[0] + l_p * v1[0];
  prj[1] = e1[1] + l_p * v1[1];
  prj[2] = e1[2] + l_p * v1[2];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// projpt2face()    Return the projection point from a point to a face.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::projpt2face(REAL* p, REAL* f1, REAL* f2, REAL* f3, REAL* prj)
{
  REAL fnormal[3], v1[3];
  REAL len, dist;

  // Get the unit face normal.
  facenormal(f1, f2, f3, fnormal, 1, NULL);
  len = sqrt(fnormal[0]*fnormal[0] + fnormal[1]*fnormal[1] + 
             fnormal[2]*fnormal[2]);
  fnormal[0] /= len;
  fnormal[1] /= len;
  fnormal[2] /= len;
  // Get the vector v1 = |p - f1|.
  v1[0] = p[0] - f1[0];
  v1[1] = p[1] - f1[1];
  v1[2] = p[2] - f1[2];
  // Get the project distance.
  dist = dot(fnormal, v1);
  
  // Get the project point.
  prj[0] = p[0] - dist * fnormal[0];
  prj[1] = p[1] - dist * fnormal[1];
  prj[2] = p[2] - dist * fnormal[2];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// facedihedral()    Return the dihedral angle (in radian) between two       //
//                   adjoining faces.                                        //
//                                                                           //
// 'pa', 'pb' are the shared edge of these two faces, 'pc1', and 'pc2' are   //
// apexes of these two faces.  Return the angle (between 0 to 2*pi) between  //
// the normal of face (pa, pb, pc1) and normal of face (pa, pb, pc2).        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::facedihedral(REAL* pa, REAL* pb, REAL* pc1, REAL* pc2)
{
  REAL n1[3], n2[3];
  REAL n1len, n2len;
  REAL costheta, ori;
  REAL theta;

  facenormal(pa, pb, pc1, n1, 1, NULL);
  facenormal(pa, pb, pc2, n2, 1, NULL);
  n1len = sqrt(dot(n1, n1));
  n2len = sqrt(dot(n2, n2));
  costheta = dot(n1, n2) / (n1len * n2len);
  // Be careful rounding error!
  if (costheta > 1.0) {
    costheta = 1.0;
  } else if (costheta < -1.0) {
    costheta = -1.0;
  }
  theta = acos(costheta);
  ori = orient3d(pa, pb, pc1, pc2);
  if (ori > 0.0) {
    theta = 2 * PI - theta;
  }

  return theta;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetalldihedral()    Get all (six) dihedral angles of a tet.               //
//                                                                           //
// If 'cosdd' is not NULL, it returns the cosines of the 6 dihedral angles,  //
// the edge indices are given in the global array 'edge2ver'. If 'cosmaxd'   //
// (or 'cosmind') is not NULL, it returns the cosine of the maximal (or      //
// minimal) dihedral angle.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::tetalldihedral(point pa, point pb, point pc, point pd,
                                REAL* cosdd, REAL* cosmaxd, REAL* cosmind)
{
  REAL N[4][3], vol, cosd, len;
  int f1 = 0, f2 = 0, i, j;

  vol = 0; // Check if the tet is valid or not.

  // Get four normals of faces of the tet.
  tetallnormal(pa, pb, pc, pd, N, &vol);

  if (vol > 0) {
    // Normalize the normals.
    for (i = 0; i < 4; i++) {
      len = sqrt(dot(N[i], N[i]));
      if (len != 0.0) {
        for (j = 0; j < 3; j++) N[i][j] /= len;
      } else {
        // There are degeneracies, such as duplicated vertices.
        vol = 0; //assert(0);
      }
    }
  }

  if (vol <= 0) { // if (vol == 0.0) {
    // A degenerated tet or an inverted tet.
    facenormal(pc, pb, pd, N[0], 1, NULL);
    facenormal(pa, pc, pd, N[1], 1, NULL);
    facenormal(pb, pa, pd, N[2], 1, NULL);
    facenormal(pa, pb, pc, N[3], 1, NULL);
    // Normalize the normals.
    for (i = 0; i < 4; i++) {
      len = sqrt(dot(N[i], N[i]));
      if (len != 0.0) {
        for (j = 0; j < 3; j++) N[i][j] /= len;
      } else {
        // There are degeneracies, such as duplicated vertices.
        break; // Not a valid normal.
      }
    }
    if (i < 4) {
      // Do not calculate dihedral angles.
      // Set all angles be 0 degree. There will be no quality optimization for
      //   this tet! Use volume optimization to correct it.
      if (cosdd != NULL) {
        for (i = 0; i < 6; i++) {
          cosdd[i] = -1.0; // 180 degree.
        }
      }
      // This tet has zero volume.
      if (cosmaxd != NULL) {
        *cosmaxd = -1.0; // 180 degree.
      }
      if (cosmind != NULL) {
        *cosmind = -1.0; // 180 degree.
      }
      return false;
    }
  }

  // Calculate the cosine of the dihedral angles of the edges.
  for (i = 0; i < 6; i++) {
    switch (i) {
    case 0: f1 = 0; f2 = 1; break; // [c,d].
    case 1: f1 = 1; f2 = 2; break; // [a,d].
    case 2: f1 = 2; f2 = 3; break; // [a,b].
    case 3: f1 = 0; f2 = 3; break; // [b,c].
    case 4: f1 = 2; f2 = 0; break; // [b,d].
    case 5: f1 = 1; f2 = 3; break; // [a,c].
    }
    cosd = -dot(N[f1], N[f2]);
    if (cosd < -1.0) cosd = -1.0; // Rounding.
    if (cosd >  1.0) cosd =  1.0; // Rounding.
    if (cosdd) cosdd[i] = cosd;
    if (cosmaxd || cosmind) {
      if (i == 0) {
        if (cosmaxd) *cosmaxd = cosd;
        if (cosmind) *cosmind = cosd;
      } else {
        if (cosmaxd) *cosmaxd = cosd < *cosmaxd ? cosd : *cosmaxd;
        if (cosmind) *cosmind = cosd > *cosmind ? cosd : *cosmind;
      }
    }
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetallnormal()    Get the in-normals of the four faces of a given tet.    //
//                                                                           //
// Let tet be abcd. N[4][3] returns the four normals, which are: N[0] cbd,   //
// N[1] acd, N[2] bad, N[3] abc (exactly corresponding to the face indices   //
// of the mesh data structure). These normals are unnormalized.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::tetallnormal(point pa, point pb, point pc, point pd,
                              REAL N[4][3], REAL* volume)
{
  REAL A[4][4], rhs[4], D;
  int indx[4];
  int i, j;

  // get the entries of A[3][3].
  for (i = 0; i < 3; i++) A[0][i] = pa[i] - pd[i];  // d->a vec
  for (i = 0; i < 3; i++) A[1][i] = pb[i] - pd[i];  // d->b vec
  for (i = 0; i < 3; i++) A[2][i] = pc[i] - pd[i];  // d->c vec

  // Compute the inverse of matrix A, to get 3 normals of the 4 faces.
  if (lu_decmp(A, 3, indx, &D, 0)) { // Decompose the matrix just once.
    if (volume != NULL) {
      // Get the volume of the tet.
      *volume = fabs((A[indx[0]][0] * A[indx[1]][1] * A[indx[2]][2])) / 6.0;
    }
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) rhs[i] = 0.0;
      rhs[j] = 1.0;  // Positive means the inside direction
      lu_solve(A, 3, indx, rhs, 0);
      for (i = 0; i < 3; i++) N[j][i] = rhs[i];
    }
    // Get the fourth normal by summing up the first three.
    for (i = 0; i < 3; i++) N[3][i] = - N[0][i] - N[1][i] - N[2][i];
  } else {
    // The tet is degenerated.
    if (volume != NULL) {
      *volume = 0;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetaspectratio()    Calculate the aspect ratio of the tetrahedron.        //
//                                                                           //
// The aspect ratio of a tet is R/h, where R is the circumradius and h is    //
// the shortest height of the tet.                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::tetaspectratio(point pa, point pb, point pc, point pd)
{
  REAL vda[3], vdb[3], vdc[3];
  REAL N[4][3], A[4][4], rhs[4], D;
  REAL H[4], volume, radius2, minheightinv;
  int indx[4];
  int i, j; 

  // Set the matrix A = [vda, vdb, vdc]^T.
  for (i = 0; i < 3; i++) A[0][i] = vda[i] = pa[i] - pd[i];
  for (i = 0; i < 3; i++) A[1][i] = vdb[i] = pb[i] - pd[i];
  for (i = 0; i < 3; i++) A[2][i] = vdc[i] = pc[i] - pd[i];
  // Lu-decompose the matrix A.
  lu_decmp(A, 3, indx, &D, 0);
  // Get the volume of abcd.
  volume = (A[indx[0]][0] * A[indx[1]][1] * A[indx[2]][2]) / 6.0;
  // Check if it is zero.
  if (volume == 0.0) return 1.0e+200; // A degenerate tet.
  // if (volume < 0.0) volume = -volume;
  // Check the radiu-edge ratio of the tet.
  rhs[0] = 0.5 * dot(vda, vda);
  rhs[1] = 0.5 * dot(vdb, vdb);
  rhs[2] = 0.5 * dot(vdc, vdc);
  lu_solve(A, 3, indx, rhs, 0);
  // Get the circumcenter.
  // for (i = 0; i < 3; i++) circumcent[i] = pd[i] + rhs[i];
  // Get the square of the circumradius.
  radius2 = dot(rhs, rhs);

  // Compute the 4 face normals (N[0], ..., N[3]).
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) rhs[i] = 0.0;
    rhs[j] = 1.0;  // Positive means the inside direction
    lu_solve(A, 3, indx, rhs, 0);
    for (i = 0; i < 3; i++) N[j][i] = rhs[i];
  }
  // Get the fourth normal by summing up the first three.
  for (i = 0; i < 3; i++) N[3][i] = - N[0][i] - N[1][i] - N[2][i];
  // Normalized the normals.
  for (i = 0; i < 4; i++) {
    // H[i] is the inverse of the height of its corresponding face.
    H[i] = sqrt(dot(N[i], N[i]));
    // if (H[i] > 0.0) {
    //   for (j = 0; j < 3; j++) N[i][j] /= H[i];
    // }
  }
  // Get the radius of the inscribed sphere.
  // insradius = 1.0 / (H[0] + H[1] + H[2] + H[3]);
  // Get the biggest H[i] (corresponding to the smallest height).
  minheightinv = H[0];
  for (i = 1; i < 3; i++) {
    if (H[i] > minheightinv) minheightinv = H[i];
  }

  return sqrt(radius2) * minheightinv;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// circumsphere()    Calculate the smallest circumsphere (center and radius) //
//                   of the given three or four points.                      //
//                                                                           //
// The circumsphere of four points (a tetrahedron) is unique if they are not //
// degenerate. If 'pd = NULL', the smallest circumsphere of three points is  //
// the diametral sphere of the triangle if they are not degenerate.          //
//                                                                           //
// Return TRUE if the input points are not degenerate and the circumcenter   //
// and circumradius are returned in 'cent' and 'radius' respectively if they //
// are not NULLs.  Otherwise, return FALSE, the four points are co-planar.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::circumsphere(REAL* pa, REAL* pb, REAL* pc, REAL* pd, 
                              REAL* cent, REAL* radius)
{
  REAL A[4][4], rhs[4], D;
  int indx[4];

  // Compute the coefficient matrix A (3x3).
  A[0][0] = pb[0] - pa[0];
  A[0][1] = pb[1] - pa[1];
  A[0][2] = pb[2] - pa[2];
  A[1][0] = pc[0] - pa[0];
  A[1][1] = pc[1] - pa[1];
  A[1][2] = pc[2] - pa[2];
  if (pd != NULL) {
    A[2][0] = pd[0] - pa[0];
    A[2][1] = pd[1] - pa[1]; 
    A[2][2] = pd[2] - pa[2];
  } else {
    cross(A[0], A[1], A[2]);
  }

  // Compute the right hand side vector b (3x1).
  rhs[0] = 0.5 * dot(A[0], A[0]);
  rhs[1] = 0.5 * dot(A[1], A[1]);
  if (pd != NULL) {
    rhs[2] = 0.5 * dot(A[2], A[2]);
  } else {
    rhs[2] = 0.0;
  }

  // Solve the 3 by 3 equations use LU decomposition with partial pivoting
  //   and backward and forward substitute..
  if (!lu_decmp(A, 3, indx, &D, 0)) {
    if (radius != (REAL *) NULL) *radius = 0.0;
    return false;
  }    
  lu_solve(A, 3, indx, rhs, 0);
  if (cent != (REAL *) NULL) {
    cent[0] = pa[0] + rhs[0];
    cent[1] = pa[1] + rhs[1];
    cent[2] = pa[2] + rhs[2];
  }
  if (radius != (REAL *) NULL) {
    *radius = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// orthosphere()    Calulcate the orthosphere of four weighted points.       //
//                                                                           //
// A weighted point (p, P^2) can be interpreted as a sphere centered at the  //
// point 'p' with a radius 'P'. The 'height' of 'p' is pheight = p[0]^2 +    //
// p[1]^2 + p[2]^2 - P^2.                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::orthosphere(REAL* pa, REAL* pb, REAL* pc, REAL* pd,
                             REAL  aheight, REAL bheight, REAL cheight, 
                             REAL  dheight, REAL* orthocent, REAL* radius)
{
  REAL A[4][4], rhs[4], D;
  int indx[4];

  // Set the coefficient matrix A (4 x 4).
  A[0][0] = 1.0; A[0][1] = pa[0]; A[0][2] = pa[1]; A[0][3] = pa[2];
  A[1][0] = 1.0; A[1][1] = pb[0]; A[1][2] = pb[1]; A[1][3] = pb[2];
  A[2][0] = 1.0; A[2][1] = pc[0]; A[2][2] = pc[1]; A[2][3] = pc[2];
  A[3][0] = 1.0; A[3][1] = pd[0]; A[3][2] = pd[1]; A[3][3] = pd[2];

  // Set the right hand side vector (4 x 1).
  rhs[0] = 0.5 * aheight;
  rhs[1] = 0.5 * bheight;
  rhs[2] = 0.5 * cheight;
  rhs[3] = 0.5 * dheight;

  // Solve the 4 by 4 equations use LU decomposition with partial pivoting
  //   and backward and forward substitute..
  if (!lu_decmp(A, 4, indx, &D, 0)) {
    if (radius != (REAL *) NULL) *radius = 0.0;
    return false;
  }
  lu_solve(A, 4, indx, rhs, 0);

  if (orthocent != (REAL *) NULL) {
    orthocent[0] = rhs[1];
    orthocent[1] = rhs[2];
    orthocent[2] = rhs[3];
  }
  if (radius != (REAL *) NULL) {
    // rhs[0] = - rheight / 2;
    // rheight  = - 2 * rhs[0];
    //          =  r[0]^2 + r[1]^2 + r[2]^2 - radius^2
    // radius^2 = r[0]^2 + r[1]^2 + r[2]^2 -rheight
    //          = r[0]^2 + r[1]^2 + r[2]^2 + 2 * rhs[0]
    *radius = sqrt(rhs[1] * rhs[1] + rhs[2] * rhs[2] + rhs[3] * rhs[3]
                   + 2.0 * rhs[0]);
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// planelineint()    Calculate the intersection of a line and a plane.       //
//                                                                           //
// The equation of a plane (points P are on the plane with normal N and P3   //
// on the plane) can be written as: N dot (P - P3) = 0. The equation of the  //
// line (points P on the line passing through P1 and P2) can be written as:  //
// P = P1 + u (P2 - P1). The intersection of these two occurs when:          //
//   N dot (P1 + u (P2 - P1)) = N dot P3.                                    //
// Solving for u gives:                                                      //
//         N dot (P3 - P1)                                                   //
//   u = ------------------.                                                 //
//         N dot (P2 - P1)                                                   //
// If the denominator is 0 then N (the normal to the plane) is perpendicular //
// to the line.  Thus the line is either parallel to the plane and there are //
// no solutions or the line is on the plane in which case there are an infi- //
// nite number of solutions.                                                 //
//                                                                           //
// The plane is given by three points pa, pb, and pc, e1 and e2 defines the  //
// line. If u is non-zero, The intersection point (if exists) returns in ip. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::planelineint(REAL* pa, REAL* pb, REAL* pc, REAL* e1, REAL* e2,
                              REAL* ip, REAL* u)
{
  REAL n[3], det, det1;

  // Calculate N.
  facenormal(pa, pb, pc, n, 1, NULL);
  // Calculate N dot (e2 - e1).
  det = n[0] * (e2[0] - e1[0]) + n[1] * (e2[1] - e1[1])
      + n[2] * (e2[2] - e1[2]);
  if (det != 0.0) {
    // Calculate N dot (pa - e1)
    det1 = n[0] * (pa[0] - e1[0]) + n[1] * (pa[1] - e1[1])
         + n[2] * (pa[2] - e1[2]);
    *u = det1 / det;
    ip[0] = e1[0] + *u * (e2[0] - e1[0]);
    ip[1] = e1[1] + *u * (e2[1] - e1[1]);
    ip[2] = e1[2] + *u * (e2[2] - e1[2]);
  } else {
    *u = 0.0;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// linelineint()    Calculate the intersection(s) of two line segments.      //
//                                                                           //
// Calculate the line segment [P, Q] that is the shortest route between two  //
// lines from A to B and C to D. Calculate also the values of tp and tq      //
// where: P = A + tp (B - A), and Q = C + tq (D - C).                        //
//                                                                           //
// Return 1 if the line segment exists. Otherwise, return 0.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::linelineint(REAL* A, REAL* B, REAL* C, REAL* D, REAL* P, 
		            REAL* Q, REAL* tp, REAL* tq)
{
  REAL vab[3], vcd[3], vca[3];
  REAL vab_vab, vcd_vcd, vab_vcd;
  REAL vca_vab, vca_vcd;
  REAL det, eps;
  int i;

  for (i = 0; i < 3; i++) {
    vab[i] = B[i] - A[i];
    vcd[i] = D[i] - C[i];
    vca[i] = A[i] - C[i];
  }

  vab_vab = dot(vab, vab);
  vcd_vcd = dot(vcd, vcd);
  vab_vcd = dot(vab, vcd);

  det = vab_vab * vcd_vcd - vab_vcd * vab_vcd;
  // Round the result.
  eps = det / (fabs(vab_vab * vcd_vcd) + fabs(vab_vcd * vab_vcd));
  if (eps < b->epsilon) {
    return 0;
  }

  vca_vab = dot(vca, vab);
  vca_vcd = dot(vca, vcd);

  *tp = (vcd_vcd * (- vca_vab) + vab_vcd * vca_vcd) / det;
  *tq = (vab_vcd * (- vca_vab) + vab_vab * vca_vcd) / det;

  for (i = 0; i < 3; i++) P[i] = A[i] + (*tp) * vab[i];
  for (i = 0; i < 3; i++) Q[i] = C[i] + (*tq) * vcd[i];

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetprismvol()    Calculate the volume of a tetrahedral prism in 4D.       //
//                                                                           //
// A tetrahedral prism is a convex uniform polychoron (four dimensional poly-//
// tope). It has 6 polyhedral cells: 2 tetrahedra connected by 4 triangular  //
// prisms. It has 14 faces: 8 triangular and 6 square. It has 16 edges and 8 //
// vertices. (Wikipedia).                                                    //
//                                                                           //
// Let 'p0', ..., 'p3' be four affinely independent points in R^3. They form //
// the lower tetrahedral facet of the prism.  The top tetrahedral facet is   //
// formed by four vertices, 'p4', ..., 'p7' in R^4, which is obtained by     //
// lifting each vertex of the lower facet into R^4 by a weight (height).  A  //
// canonical choice of the weights is the square of Euclidean norm of of the //
// points (vectors).                                                         //
//                                                                           //
//                                                                           //
// The return value is (4!) 24 times of the volume of the tetrahedral prism. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::tetprismvol(REAL* p0, REAL* p1, REAL* p2, REAL* p3)
{
  REAL *p4, *p5, *p6, *p7;
  REAL w4, w5, w6, w7;
  REAL vol[4];

  p4 = p0;
  p5 = p1;
  p6 = p2;
  p7 = p3;

  // TO DO: these weights can be pre-calculated!
  w4 = dot(p0, p0);
  w5 = dot(p1, p1);
  w6 = dot(p2, p2);
  w7 = dot(p3, p3);

  // Calculate the volume of the tet-prism.
  vol[0] = orient4d(p5, p6, p4, p3, p7, w5, w6, w4, 0, w7);
  vol[1] = orient4d(p3, p6, p2, p0, p1,  0, w6,  0, 0,  0);
  vol[2] = orient4d(p4, p6, p3, p0, p1, w4, w6,  0, 0,  0);
  vol[3] = orient4d(p6, p5, p4, p3, p1, w6, w5, w4, 0,  0);

  return fabs(vol[0]) + fabs(vol[1]) + fabs(vol[2]) + fabs(vol[3]);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// calculateabovepoint()    Calculate a point above a facet in 'dummypoint'. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::calculateabovepoint(arraypool *facpoints, point *ppa,
                                     point *ppb, point *ppc)
{
  point *ppt, pa, pb, pc;
  REAL v1[3], v2[3], n[3];
  REAL lab, len, A, area;
  REAL x, y, z;
  int i;

  ppt = (point *) fastlookup(facpoints, 0);
  pa = *ppt; // a is the first point.
  pb = pc = NULL; // Avoid compiler warnings.

  // Get a point b s.t. the length of [a, b] is maximal.
  lab = 0;
  for (i = 1; i < facpoints->objects; i++) {
    ppt = (point *) fastlookup(facpoints, i);
    x = (*ppt)[0] - pa[0];
    y = (*ppt)[1] - pa[1];
    z = (*ppt)[2] - pa[2];
    len = x * x + y * y + z * z;
    if (len > lab) {
      lab = len;
      pb = *ppt;
    }
  }
  lab = sqrt(lab);
  if (lab == 0) {
    if (!b->quiet) {
      printf("Warning:  All points of a facet are coincident with %d.\n",
        pointmark(pa));
    }
    return false;
  }

  // Get a point c s.t. the area of [a, b, c] is maximal.
  v1[0] = pb[0] - pa[0];
  v1[1] = pb[1] - pa[1];
  v1[2] = pb[2] - pa[2];
  A = 0;
  for (i = 1; i < facpoints->objects; i++) {
    ppt = (point *) fastlookup(facpoints, i);
    v2[0] = (*ppt)[0] - pa[0];
    v2[1] = (*ppt)[1] - pa[1];
    v2[2] = (*ppt)[2] - pa[2];
    cross(v1, v2, n);
    area = dot(n, n);
    if (area > A) {
      A = area;
      pc = *ppt;
    }
  }
  if (A == 0) {
    // All points are collinear. No above point.
    if (!b->quiet) {
      printf("Warning:  All points of a facet are collinaer with [%d, %d].\n",
        pointmark(pa), pointmark(pb));
    }
    return false;
  }

  // Calculate an above point of this facet.
  facenormal(pa, pb, pc, n, 1, NULL);
  len = sqrt(dot(n, n));
  n[0] /= len;
  n[1] /= len;
  n[2] /= len;
  lab /= 2.0; // Half the maximal length.
  dummypoint[0] = pa[0] + lab * n[0];
  dummypoint[1] = pa[1] + lab * n[1];
  dummypoint[2] = pa[2] + lab * n[2];

  if (ppa != NULL) {
    // Return the three points.
    *ppa = pa;
    *ppb = pb;
    *ppc = pc;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Calculate an above point. It lies above the plane containing  the subface //
//   [a,b,c], and save it in dummypoint. Moreover, the vector pa->dummypoint //
//   is the normal of the plane.                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::calculateabovepoint4(point pa, point pb, point pc, point pd)
{
  REAL n1[3], n2[3], *norm;
  REAL len, len1, len2;

  // Select a base.
  facenormal(pa, pb, pc, n1, 1, NULL);
  len1 = sqrt(dot(n1, n1));
  facenormal(pa, pb, pd, n2, 1, NULL);
  len2 = sqrt(dot(n2, n2));
  if (len1 > len2) {
    norm = n1;
    len = len1;
  } else {
    norm = n2;
    len = len2;
  }
  assert(len > 0);
  norm[0] /= len;
  norm[1] /= len;
  norm[2] /= len;
  len = distance(pa, pb);
  dummypoint[0] = pa[0] + len * norm[0];
  dummypoint[1] = pa[1] + len * norm[1];
  dummypoint[2] = pa[2] + len * norm[2];
}

////                                                                       ////
////                                                                       ////
//// geom_cxx /////////////////////////////////////////////////////////////////

