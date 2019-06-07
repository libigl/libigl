// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
const char *__doc_igl_active_set = R"igl_Qu8mg5v7(// Known Bugs: rows of [Aeq;Aieq] **must** be linearly independent. Should be
  // using QR decomposition otherwise:
  //   http://www.okstate.edu/sas/v8/sashtml/ormp/chap5/sect32.htm
  //
  // ACTIVE_SET Minimize quadratic energy
  //
  // 0.5*Z'*A*Z + Z'*B + C with constraints
  //
  // that Z(known) = Y, optionally also subject to the constraints Aeq*Z = Beq,
  // and further optionally subject to the linear inequality constraints that
  // Aieq*Z <= Bieq and constant inequality constraints lx <= x <= ux
  //
  // Inputs:
  //   A  n by n matrix of quadratic coefficients
  //   B  n by 1 column of linear coefficients
  //   known  list of indices to known rows in Z
  //   Y  list of fixed values corresponding to known rows in Z
  //   Aeq  meq by n list of linear equality constraint coefficients
  //   Beq  meq by 1 list of linear equality constraint constant values
  //   Aieq  mieq by n list of linear inequality constraint coefficients
  //   Bieq  mieq by 1 list of linear inequality constraint constant values
  //   lx  n by 1 list of lower bounds [] implies -Inf
  //   ux  n by 1 list of upper bounds [] implies Inf
  //   params  struct of additional parameters (see below)
  //   Z  if not empty, is taken to be an n by 1 list of initial guess values
  //     (see output)
  // Outputs:
  //   Z  n by 1 list of solution values
  // Returns true on success, false on error
  //
  // Benchmark: For a harmonic solve on a mesh with 325K facets, matlab 2.2
  // secs, igl/min_quad_with_fixed.h 7.1 secs
  //)igl_Qu8mg5v7";
const char *__doc_igl_adjacency_list = R"igl_Qu8mg5v7(// Constructs the graph adjacency list of a given mesh (V,F)
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   F       #F by dim list of mesh faces (must be triangles)
  //   sorted  flag that indicates if the list should be sorted counter-clockwise
  // Outputs:
  //   A  vector<vector<T> > containing at row i the adjacent vertices of vertex i
  //
  // Example:
  //   // Mesh in (V,F)
  //   vector<vector<double> > A;
  //   adjacency_list(F,A);
  //
  // See also: edges, cotmatrix, diag)igl_Qu8mg5v7";
const char *__doc_igl_arap_precomputation = R"igl_Qu8mg5v7(// Compute necessary information to start using an ARAP deformation
  //
  // Inputs:
  //   V  #V by dim list of mesh positions
  //   F  #F by simplex-size list of triangle|tet indices into V
  //   dim  dimension being used at solve time. For deformation usually dim =
  //     V.cols(), for surface parameterization V.cols() = 3 and dim = 2
  //   b  #b list of "boundary" fixed vertex indices into V
  // Outputs:
  //   data  struct containing necessary precomputation)igl_Qu8mg5v7";
const char *__doc_igl_arap_solve = R"igl_Qu8mg5v7(// Inputs:
  //   bc  #b by dim list of boundary conditions
  //   data  struct containing necessary precomputation and parameters
  //   U  #V by dim initial guess)igl_Qu8mg5v7";
const char *__doc_igl_avg_edge_length = R"igl_Qu8mg5v7(// Compute the average edge length for the given triangle mesh
  // Templates:
  //   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
  //   DerivedF derived from face indices matrix type: i.e. MatrixXi
  //   DerivedL derived from edge lengths matrix type: i.e. MatrixXd
  // Inputs:
  //   V  eigen matrix #V by 3
  //   F  #F by simplex-size list of mesh faces (must be simplex)
  // Outputs:
  //   l  average edge length
  //
  // See also: adjacency_matrix)igl_Qu8mg5v7";
const char *__doc_igl_barycenter = R"igl_Qu8mg5v7(// Computes the barycenter of every simplex
  //
  // Inputs:
  //   V  #V x dim matrix of vertex coordinates
  //   F  #F x simplex_size  matrix of indices of simplex corners into V
  // Output:
  //   BC  #F x dim matrix of 3d vertices
  //)igl_Qu8mg5v7";
const char *__doc_igl_barycentric_coordinates = R"igl_Qu8mg5v7(// Compute barycentric coordinates in a tet
  //
  // Inputs:
  //   P  #P by 3 Query points in 3d
  //   A  #P by 3 Tet corners in 3d
  //   B  #P by 3 Tet corners in 3d
  //   C  #P by 3 Tet corners in 3d
  //   D  #P by 3 Tet corners in 3d
  // Outputs:
  //   L  #P by 4 list of barycentric coordinates
  //   )igl_Qu8mg5v7";
const char *__doc_igl_barycentric_to_global = R"igl_Qu8mg5v7(// Converts barycentric coordinates in the embree form to 3D coordinates
  // Embree stores barycentric coordinates as triples: fid, bc1, bc2
  // fid is the id of a face, bc1 is the displacement of the point wrt the
  // first vertex v0 and the edge v1-v0. Similarly, bc2 is the displacement
  // wrt v2-v0.
  //
  // Input:
  // V:  #Vx3 Vertices of the mesh
  // F:  #Fxe Faces of the mesh
  // bc: #Xx3 Barycentric coordinates, one row per point
  //
  // Output:
  // #X: #Xx3 3D coordinates of all points in bc)igl_Qu8mg5v7";

const char *__doc_igl_bbw = R"igl_Qu8mg5v7(// Compute Bounded Biharmonic Weights on a given domain (V,Ele) with a given
  // set of boundary conditions
  //
  // Templates
  //   DerivedV  derived type of eigen matrix for V (e.g. MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. MatrixXi)
  //   Derivedb  derived type of eigen matrix for b (e.g. VectorXi)
  //   Derivedbc  derived type of eigen matrix for bc (e.g. MatrixXd)
  //   DerivedW  derived type of eigen matrix for W (e.g. MatrixXd)
  // Inputs:
  //   V  #V by dim vertex positions
  //   Ele  #Elements by simplex-size list of element indices
  //   b  #b boundary indices into V
  //   bc #b by #W list of boundary values
  //   data  object containing options, initial guess --> solution and results
  // Outputs:
  //   W  #V by #W list of *unnormalized* weights to normalize use
  //    igl::normalize_row_sums(W,W);
  // Returns true on success, false on failure)igl_Qu8mg5v7";
const char *__doc_igl_boundary_conditions = R"igl_Qu8mg5v7(// Compute boundary conditions for automatic weights computation. This
  // function expects that the given mesh (V,Ele) has sufficient samples
  // (vertices) exactly at point handle locations and exactly along bone and
  // cage edges.
  //
  // Inputs:
  //   V  #V by dim list of domain vertices
  //   Ele  #Ele by simplex-size list of simplex indices
  //   C  #C by dim list of handle positions
  //   P  #P by 1 list of point handle indices into C
  //   BE  #BE by 2 list of bone edge indices into C
  //   CE  #CE by 2 list of cage edge indices into *P*
  // Outputs:
  //   b  #b list of boundary indices (indices into V of vertices which have
  //     known, fixed values)
  //   bc #b by #weights list of known/fixed values for boundary vertices
  //     (notice the #b != #weights in general because #b will include all the
  //     intermediary samples along each bone, etc.. The ordering of the
  //     weights corresponds to [P;BE]
  // Returns false if boundary conditions are suspicious:
  //   P and BE are empty
  //   bc is empty
  //   some column of bc doesn't have a 0 (assuming bc has >1 columns)
  //   some column of bc doesn't have a 1 (assuming bc has >1 columns))igl_Qu8mg5v7";
const char *__doc_igl_boundary_facets = R"igl_Qu8mg5v7(// BOUNDARY_FACETS Determine boundary faces (edges) of tetrahedra (triangles)
  // stored in T (analogous to qptoolbox's `outline` and `boundary_faces`).
  //
  // Templates:
  //   IntegerT  integer-value: e.g. int
  //   IntegerF  integer-value: e.g. int
  // Input:
  //  T  tetrahedron (triangle) index list, m by 4 (3), where m is the number of tetrahedra
  // Output:
  //  F  list of boundary faces, n by 3 (2), where n is the number of boundary faces
  //
  //)igl_Qu8mg5v7";
const char *__doc_igl_boundary_loop = R"igl_Qu8mg5v7(// Compute list of ordered boundary loops for a manifold mesh.
  //
  // Templates:
  //  Index  index type
  // Inputs:
  //   F  #V by dim list of mesh faces
  // Outputs:
  //   L  list of loops where L[i] = ordered list of boundary vertices in loop i
  //)igl_Qu8mg5v7";
const char *__doc_igl_cat = R"igl_Qu8mg5v7(// Perform concatenation of a two matrices along a single dimension
  // If dim == 1, then C = [A;B]. If dim == 2 then C = [A B]
  //
  // Template:
  //   Scalar  scalar data type for sparse matrices like double or int
  //   Mat  matrix type for all matrices (e.g. MatrixXd, SparseMatrix)
  //   MatC  matrix type for output matrix (e.g. MatrixXd) needs to support
  //     resize
  // Inputs:
  //   A  first input matrix
  //   B  second input matrix
  //   dim  dimension along which to concatenate, 1 or 2
  // Outputs:
  //   C  output matrix
  //   )igl_Qu8mg5v7";
const char *__doc_igl_collapse_edge = R"igl_Qu8mg5v7(See collapse_edge for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_colon = R"igl_Qu8mg5v7(// Colon operator like matlab's colon operator. Enumerats values between low
  // and hi with step step.
  // Templates:
  //   L  should be a eigen matrix primitive type like int or double
  //   S  should be a eigen matrix primitive type like int or double
  //   H  should be a eigen matrix primitive type like int or double
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   low  starting value if step is valid then this is *always* the first
  //     element of I
  //   step  step difference between sequential elements returned in I,
  //     remember this will be cast to template T at compile time. If low<hi
  //     then step must be positive. If low>hi then step must be negative.
  //     Otherwise I will be set to empty.
  //   hi  ending value, if (hi-low)%step is zero then this will be the last
  //     element in I. If step is positive there will be no elements greater
  //     than hi, vice versa if hi<low
  // Output:
  //   I  list of values from low to hi with step size step)igl_Qu8mg5v7";
const char *__doc_igl_column_to_quats = R"igl_Qu8mg5v7(// "Columnize" a list of quaternions (q1x,q1y,q1z,q1w,q2x,q2y,q2z,q2w,...)
  //
  // Inputs:
  //   Q  n*4-long list of coefficients
  // Outputs:
  //   vQ  n-long list of quaternions
  // Returns false if n%4!=0)igl_Qu8mg5v7";
const char *__doc_igl_comb_cross_field = R"igl_Qu8mg5v7(// Inputs:
  //   V          #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F          #F by 4 eigen Matrix of face (quad) indices
  //   PD1in      #F by 3 eigen Matrix of the first per face cross field vector
  //   PD2in      #F by 3 eigen Matrix of the second per face cross field vector
  // Output:
  //   PD1out      #F by 3 eigen Matrix of the first combed cross field vector
  //   PD2out      #F by 3 eigen Matrix of the second combed cross field vector
  //)igl_Qu8mg5v7";
const char *__doc_igl_comb_frame_field = R"igl_Qu8mg5v7(// Inputs:
  //   V            #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F            #F by 4 eigen Matrix of face (quad) indices
  //   PD1          #F by 3 eigen Matrix of the first per face cross field vector
  //   PD2          #F by 3 eigen Matrix of the second per face cross field vector
  //   BIS1_combed  #F by 3 eigen Matrix of the first combed bisector field vector
  //   BIS2_combed  #F by 3 eigen Matrix of the second combed bisector field vector
  // Output:
  //   PD1_combed  #F by 3 eigen Matrix of the first combed cross field vector
  //   PD2_combed  #F by 3 eigen Matrix of the second combed cross field vector
  //)igl_Qu8mg5v7";
const char *__doc_igl_compute_frame_field_bisectors = R"igl_Qu8mg5v7(// Compute bisectors of a frame field defined on mesh faces
  // Inputs:
  //   V     #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F     #F by 3 eigen Matrix of face (triangle) indices
  //   B1    #F by 3 eigen Matrix of face (triangle) base vector 1
  //   B2    #F by 3 eigen Matrix of face (triangle) base vector 2
  //   PD1   #F by 3 eigen Matrix of the first per face frame field vector
  //   PD2   #F by 3 eigen Matrix of the second per face frame field vector
  // Output:
  //   BIS1  #F by 3 eigen Matrix of the first per face frame field bisector
  //   BIS2  #F by 3 eigen Matrix of the second per face frame field bisector
  //)igl_Qu8mg5v7";
const char *__doc_igl_copyleft_cgal_mesh_boolean = R"igl_Qu8mg5v7(//  MESH_BOOLEAN Compute boolean csg operations on "solid", consistently
      //  oriented meshes.
      //
      //  Inputs:
      //    VA  #VA by 3 list of vertex positions of first mesh
      //    FA  #FA by 3 list of triangle indices into VA
      //    VB  #VB by 3 list of vertex positions of second mesh
      //    FB  #FB by 3 list of triangle indices into VB
      //    type  type of boolean operation
      //  Outputs:
      //    VC  #VC by 3 list of vertex positions of boolean result mesh
      //    FC  #FC by 3 list of triangle indices into VC
      //    J  #FC list of indices into [FA;FA.rows()+FB] revealing "birth" facet
      //  Returns true if inputs induce a piecewise constant winding number
      //  field and type is valid
      //
      //  See also: mesh_boolean_cork, intersect_other,
      //  remesh_self_intersections)igl_Qu8mg5v7";
const char *__doc_igl_copyleft_cgal_remesh_self_intersections = R"igl_Qu8mg5v7(// Given a triangle mesh (V,F) compute a new mesh (VV,FF) which is the same
      // as (V,F) except that any self-intersecting triangles in (V,F) have been
      // subdivided (new vertices and face created) so that the self-intersection
      // contour lies exactly on edges in (VV,FF). New vertices will appear in
      // original faces or on original edges. New vertices on edges are "merged"
      // only across original faces sharing that edge. This means that if the input
      // triangle mesh is a closed manifold the output will be too.
      //
      // Inputs:
      //   V  #V by 3 list of vertex positions
      //   F  #F by 3 list of triangle indices into V
      //   params  struct of optional parameters
      // Outputs:
      //   VV  #VV by 3 list of vertex positions
      //   FF  #FF by 3 list of triangle indices into VV
      //   IF  #intersecting face pairs by 2  list of intersecting face pairs,
      //     indexing F
      //   J  #FF list of indices into F denoting birth triangle
      //   IM  #VV list of indices into VV of unique vertices.
      //
      // Known bugs: If an existing edge in (V,F) lies exactly on another face then
      // any resulting additional vertices along that edge may not get properly
      // connected so that the output mesh has the same global topology. This is
      // because
      //
      // Example:
      //     // resolve intersections
      //     igl::copyleft::cgal::remesh_self_intersections(V,F,params,VV,FF,IF,J,IM);
      //     // _apply_ duplicate vertex mapping IM to FF
      //     for_each(FF.data(),FF.data()+FF.size(),[&IM](int & a){a=IM(a);});
      //     // remove any vertices now unreferenced after duplicate mapping.
      //     igl::remove_unreferenced(VV,FF,SV,SF,UIM);
      //     // Now (SV,SF) is ready to extract outer hull
      //     igl::copyleft::cgal::outer_hull(SV,SF,G,J,flip);
      //)igl_Qu8mg5v7";
const char *__doc_igl_copyleft_comiso_miq = R"igl_Qu8mg5v7(// Inputs:
    //   V              #V by 3 list of mesh vertex 3D positions
    //   F              #F by 3 list of faces indices in V
    //   PD1            #V by 3 first line of the Jacobian per triangle
    //   PD2            #V by 3 second line of the Jacobian per triangle
    //                  (optional, if empty it will be a vector in the tangent plane orthogonal to PD1)
    //   scale          global scaling for the gradient (controls the quads resolution)
    //   stiffness      weight for the stiffness iterations
    //   direct_round   greedily round all integer variables at once (greatly improves optimization speed but lowers quality)
    //   iter           stiffness iterations (0 = no stiffness)
    //   local_iter     number of local iterations for the integer rounding
    //   do_round       enables the integer rounding (disabling it could be useful for debugging)
    //   round_vertices id of additional vertices that should be snapped to integer coordinates
    //   hard_features  #H by 2 list of pairs of vertices that belongs to edges that should be snapped to integer coordinates
    //
    // Output:
    //   UV             #UV by 2 list of vertices in 2D
    //   FUV            #FUV by 3 list of face indices in UV
    //
    // TODO: rename the parameters name in the cpp consistently
    //       improve the handling of hard_features, right now it might fail in difficult cases)igl_Qu8mg5v7";
const char *__doc_igl_copyleft_comiso_nrosy = R"igl_Qu8mg5v7(// Generate a N-RoSy field from a sparse set of constraints
    //
    // Inputs:
    //   V       #V by 3 list of mesh vertex coordinates
    //   F       #F by 3 list of mesh faces (must be triangles)
    //   b       #B by 1 list of constrained face indices
    //   bc      #B by 3 list of representative vectors for the constrained
    //     faces
    //   b_soft  #S by 1 b for soft constraints
    //   w_soft  #S by 1 weight for the soft constraints (0-1)
    //   bc_soft #S by 3 bc for soft constraints
    //   N       the degree of the N-RoSy vector field
    //   soft    the strength of the soft constraints w.r.t. smoothness
    //           (0 -> smoothness only, 1->constraints only)
    // Outputs:
    //   R       #F by 3 the representative vectors of the interpolated field
    //   S       #V by 1 the singularity index for each vertex (0 = regular))igl_Qu8mg5v7";
const char *__doc_igl_copyleft_marching_cubes = R"igl_Qu8mg5v7(// marching_cubes( values, points, x_res, y_res, z_res, vertices, faces )
    //
    // performs marching cubes reconstruction on the grid defined by values, and
    // points, and generates vertices and faces
    //
    // Input:
    //  values  #number_of_grid_points x 1 array -- the scalar values of an
    //    implicit function defined on the grid points (<0 in the inside of the
    //    surface, 0 on the border, >0 outside)
    //  points  #number_of_grid_points x 3 array -- 3-D positions of the grid
    //    points, ordered in x,y,z order:
    //      points[index] = the point at (x,y,z) where :
    //      x = (index % (xres -1),
    //      y = (index / (xres-1)) %(yres-1),
    //      z = index / (xres -1) / (yres -1) ).
    //      where x,y,z index x, y, z dimensions
    //      i.e. index = x + y*xres + z*xres*yres
    //  xres  resolutions of the grid in x dimension
    //  yres  resolutions of the grid in y dimension
    //  zres  resolutions of the grid in z dimension
    // Output:
    //   vertices  #V by 3 list of mesh vertex positions
    //   faces  #F by 3 list of mesh triangle indices
    //)igl_Qu8mg5v7";
const char *__doc_igl_copyleft_swept_volume = R"igl_Qu8mg5v7(// Compute the surface of the swept volume of a solid object with surface
    // (V,F) mesh under going rigid motion.
    //
    // Inputs:
    //   V  #V by 3 list of mesh positions in reference pose
    //   F  #F by 3 list of mesh indices into V
    //   transform  function handle so that transform(t) returns the rigid
    //     transformation at time t∈[0,1]
    //   steps  number of time steps: steps=3 --> t∈{0,0.5,1}
    //   grid_res  number of grid cells on the longest side containing the
    //     motion (isolevel+1 cells will also be added on each side as padding)
    //   isolevel  distance level to be contoured as swept volume
    // Outputs:
    //   SV  #SV by 3 list of mesh positions of the swept surface
    //   SF  #SF by 3 list of mesh faces into SV)igl_Qu8mg5v7";
const char *__doc_igl_copyleft_tetgen_tetrahedralize = R"igl_Qu8mg5v7(// Mesh the interior of a surface mesh (V,F) using tetgen
      //
      // Inputs:
      //   V  #V by 3 vertex position list
      //   F  #F list of polygon face indices into V (0-indexed)
      //   switches  string of tetgen options (See tetgen documentation) e.g.
      //     "pq1.414a0.01" tries to mesh the interior of a given surface with
      //       quality and area constraints
      //     "" will mesh the convex hull constrained to pass through V (ignores F)
      // Outputs:
      //   TV  #V by 3 vertex position list
      //   TT  #T by 4 list of tet face indices
      //   TF  #F by 3 list of triangle face indices
      // Returns status:
      //   0 success
      //   1 tetgen threw exception
      //   2 tetgen did not crash but could not create any tets (probably there are
      //     holes, duplicate faces etc.)
      //   -1 other error)igl_Qu8mg5v7";
const char *__doc_igl_cotmatrix = R"igl_Qu8mg5v7(// Constructs the cotangent stiffness matrix (discrete laplacian) for a given
  // mesh (V,F).
  //
  // Templates:
  //   DerivedV  derived type of eigen matrix for V (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. derived from
  //     MatrixXi)
  //   Scalar  scalar type for eigen sparse matrix (e.g. double)
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by simplex_size list of mesh faces (must be triangles)
  // Outputs:
  //   L  #V by #V cotangent matrix, each row i corresponding to V(i,:)
  //
  // See also: adjacency_matrix
  //
  // Note: This Laplacian uses the convention that diagonal entries are
  // **minus** the sum of off-diagonal entries. The diagonal entries are
  // therefore in general negative and the matrix is **negative** semi-definite
  // (immediately, -L is **positive** semi-definite)
  //)igl_Qu8mg5v7";
const char *__doc_igl_covariance_scatter_matrix = R"igl_Qu8mg5v7(// Construct the covariance scatter matrix for a given arap energy
  // Inputs:
  //   V  #V by Vdim list of initial domain positions
  //   F  #F by 3 list of triangle indices into V
  //   energy  ARAPEnergyType enum value defining which energy is being used.
  //     See ARAPEnergyType.h for valid options and explanations.
  // Outputs:
  //   CSM dim*#V/#F by dim*#V sparse matrix containing special laplacians along
  //     the diagonal so that when multiplied by V gives covariance matrix
  //     elements, can be used to speed up covariance matrix computation)igl_Qu8mg5v7";
const char *__doc_igl_cross_field_mismatch = R"igl_Qu8mg5v7(// Inputs:
  //   V         #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F         #F by 3 eigen Matrix of face (quad) indices
  //   PD1       #F by 3 eigen Matrix of the first per face cross field vector
  //   PD2       #F by 3 eigen Matrix of the second per face cross field vector
  //   isCombed  boolean, specifying whether the field is combed (i.e. matching has been precomputed.
  //             If not, the field is combed first.
  // Output:
  //   Handle_MMatch    #F by 3 eigen Matrix containing the integer mismatch of the cross field
  //                    across all face edges
  //)igl_Qu8mg5v7";
const char *__doc_igl_cut_mesh_from_singularities = R"igl_Qu8mg5v7(// Given a mesh (V,F) and the integer mismatch of a cross field per edge
  // (mismatch), finds the cut_graph connecting the singularities (seams) and the
  // degree of the singularities singularity_index
  //
  // Input:
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of faces
  //   mismatch  #F by 3 list of per corner integer mismatch
  // Outputs:
  //   seams  #F by 3 list of per corner booleans that denotes if an edge is a
  //     seam or not
  //)igl_Qu8mg5v7";
const char *__doc_igl_deform_skeleton = R"igl_Qu8mg5v7(// Deform a skeleton.
  //
  // Inputs:
  //   C  #C by 3 list of joint positions
  //   BE  #BE by 2 list of bone edge indices
  //   vA  #BE list of bone transformations
  // Outputs
  //   CT  #BE*2 by 3 list of deformed joint positions
  //   BET  #BE by 2 list of bone edge indices (maintains order)
  //)igl_Qu8mg5v7";
const char *__doc_igl_directed_edge_orientations = R"igl_Qu8mg5v7(// Determine rotations that take each edge from the x-axis to its given rest
  // orientation.
  //
  // Inputs:
  //   C  #C by 3 list of edge vertex positions
  //   E  #E by 2 list of directed edges
  // Outputs:
  //   Q  #E list of quaternions
  //)igl_Qu8mg5v7";
const char *__doc_igl_directed_edge_parents = R"igl_Qu8mg5v7(// Recover "parents" (preceding edges) in a tree given just directed edges.
  //
  // Inputs:
  //   E  #E by 2 list of directed edges
  // Outputs:
  //   P  #E list of parent indices into E (-1) means root
  //)igl_Qu8mg5v7";
const char *__doc_igl_doublearea = R"igl_Qu8mg5v7(// DOUBLEAREA computes twice the area for each input triangle[quad]
  //
  // Templates:
  //   DerivedV  derived type of eigen matrix for V (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. derived from
  //     MatrixXi)
  //   DeriveddblA  derived type of eigen matrix for dblA (e.g. derived from
  //     MatrixXd)
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by simplex_size list of mesh faces (must be triangles or quads)
  // Outputs:
  //   dblA  #F list of triangle[quad] double areas (SIGNED only for 2D input)
  //
  // Known bug: For dim==3 complexity is O(#V + #F)!! Not just O(#F). This is a big deal
  // if you have 1million unreferenced vertices and 1 face)igl_Qu8mg5v7";
const char *__doc_igl_doublearea_single = R"igl_Qu8mg5v7(// Single triangle in 2D!
  //
  // This should handle streams of corners not just single corners)igl_Qu8mg5v7";
const char *__doc_igl_doublearea_quad = R"igl_Qu8mg5v7(// DOUBLEAREA_QUAD computes twice the area for each input quadrilateral
  //
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by simplex_size list of mesh faces (must be quadrilaterals)
  // Outputs:
  //   dblA  #F list of quadrilateral double areas
  //)igl_Qu8mg5v7";
const char *__doc_igl_dqs = R"igl_Qu8mg5v7(// Dual quaternion skinning
  //
  // Inputs:
  //   V  #V by 3 list of rest positions
  //   W  #W by #C list of weights
  //   vQ  #C list of rotation quaternions
  //   vT  #C list of translation vectors
  // Outputs:
  //   U  #V by 3 list of new positions)igl_Qu8mg5v7";
const char *__doc_igl_edge_lengths = R"igl_Qu8mg5v7(// Constructs a list of lengths of edges opposite each index in a face
  // (triangle/tet) list
  //
  // Templates:
  //   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
  //   DerivedF derived from face indices matrix type: i.e. MatrixXi
  //   DerivedL derived from edge lengths matrix type: i.e. MatrixXd
  // Inputs:
  //   V  eigen matrix #V by 3
  //   F  #F by 2 list of mesh edges
  //    or
  //   F  #F by 3 list of mesh faces (must be triangles)
  //    or
  //   T  #T by 4 list of mesh elements (must be tets)
  // Outputs:
  //   L  #F by {1|3|6} list of edge lengths
  //     for edges, column of lengths
  //     for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  //     for tets, columns correspond to edges
  //     [3 0],[3 1],[3 2],[1 2],[2 0],[0 1]
  //)igl_Qu8mg5v7";
const char *__doc_igl_edge_topology = R"igl_Qu8mg5v7(// Initialize Edges and their topological relations (assumes an edge-manifold
  // mesh)
  //
  // Output:
  // EV  : #Ex2, Stores the edge description as pair of indices to vertices
  // FE : #Fx3, Stores the Triangle-Edge relation
  // EF : #Ex2: Stores the Edge-Triangle relation
  //
  // TODO: This seems to be a inferior duplicate of edge_flaps.h:
  //   - unused input parameter V
  //   - roughly 2x slower than edge_flaps
  //   - outputs less information: edge_flaps reveals corner opposite edge
  //   - FE uses non-standard and ambiguous order: FE(f,c) is merely an edge
  //     incident on corner c of face f. In contrast, edge_flaps's EMAP(f,c) reveals
  //     the edge _opposite_ corner c of face f)igl_Qu8mg5v7";
const char *__doc_igl_eigs = R"igl_Qu8mg5v7(See eigs for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_embree_ambient_occlusion = R"igl_Qu8mg5v7(// Compute ambient occlusion per given point
    //
    // Inputs:
    //    ei  EmbreeIntersector containing (V,F)
    //    P  #P by 3 list of origin points
    //    N  #P by 3 list of origin normals
    // Outputs:
    //    S  #P list of ambient occlusion values between 1 (fully occluded) and
    //      0 (not occluded)
    //)igl_Qu8mg5v7";
const char *__doc_igl_embree_line_mesh_intersection = R"igl_Qu8mg5v7(// Project the point cloud V_source onto the triangle mesh
    // V_target,F_target.
    // A ray is casted for every vertex in the direction specified by
    // N_source and its opposite.
    //
    // Input:
    // V_source: #Vx3 Vertices of the source mesh
    // N_source: #Vx3 Normals of the point cloud
    // V_target: #V2x3 Vertices of the target mesh
    // F_target: #F2x3 Faces of the target mesh
    //
    // Output:
    // #Vx3 matrix of baricentric coordinate. Each row corresponds to
    // a vertex of the projected mesh and it has the following format:
    // id b1 b2. id is the id of a face of the source mesh. b1 and b2 are
    // the barycentric coordinates wrt the first two edges of the triangle
    // To convert to standard global coordinates, see barycentric_to_global.h)igl_Qu8mg5v7";
const char *__doc_igl_embree_reorient_facets_raycast = R"igl_Qu8mg5v7(// Orient each component (identified by C) of a mesh (V,F) using ambient
    // occlusion such that the front side is less occluded than back side, as
    // described in "A Simple Method for Correcting Facet Orientations in
    // Polygon Meshes Based on Ray Casting" [Takayama et al. 2014].
    //
    // Inputs:
    //   V  #V by 3 list of vertex positions
    //   F  #F by 3 list of triangle indices
    //   rays_total  Total number of rays that will be shot
    //   rays_minimum  Minimum number of rays that each patch should receive
    //   facet_wise  Decision made for each face independently, no use of patches
    //     (i.e., each face is treated as a patch)
    //   use_parity  Use parity mode
    //   is_verbose  Verbose output to cout
    // Outputs:
    //   I  #F list of whether face has been flipped
    //   C  #F list of patch ID (output of bfs_orient > manifold patches))igl_Qu8mg5v7";
const char *__doc_igl_exact_geodesic = R"igl_Qu8mg5v7( 
    // Exact geodesic algorithm for triangular mesh with the implementation from https://code.google.com/archive/p/geodesic/,  
    // and the algorithm first described by Mitchell, Mount and Papadimitriou in 1987 
    //  
    // Inputs: 
    //   V  #V by 3 list of 3D vertex positions 
    //   F  #F by 3 list of mesh faces 
    //   VS #VS by 1 vector specifying indices of source vertices 
    //   FS #FS by 1 vector specifying indices of source faces 
    //   VT #VT by 1 vector specifying indices of target vertices 
    //   FT #FT by 1 vector specifying indices of target faces 
    // Output: 
    //   D  #VT+#FT by 1 vector of geodesic distances of each target w.r.t. the nearest one in the source set 
    // 
    // Note:  
    //      Specifying a face as target/source means its center.  
    //)igl_Qu8mg5v7";
const char *__doc_igl_heat_geodesics_precompute = R"igl_Qu8mg5v7(
    // Precompute factorized solvers for computing a fast approximation of
    // geodesic distances on a mesh (V,F). [Crane et al. 2013]
    //
    // Inputs:
    //   V  #V by dim list of mesh vertex positions
    //   F  #F by 3 list of mesh face indices into V
    // Outputs:
    //   data  precomputation data (see heat_geodesics_solve)
    //)igl_Qu8mg5v7";
const char *__doc_igl_heat_geodesics_solve = R"igl_Qu8mg5v7(
    // Compute fast approximate geodesic distances using precomputed data from a
    // set of selected source vertices (gamma)
    //
    // Inputs:
    //   data  precomputation data (see heat_geodesics_precompute)
    //   gamma  #gamma list of indices into V of source vertices
    // Outputs:
    //   D  #V list of distances to gamma
    //)igl_Qu8mg5v7";
const char *__doc_igl_find_cross_field_singularities = R"igl_Qu8mg5v7(// Inputs:
  //   V                #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F                #F by 3 eigen Matrix of face (quad) indices
  //   Handle_MMatch    #F by 3 eigen Matrix containing the integer mismatch of the cross field
  //                    across all face edges
  // Output:
  //   isSingularity    #V by 1 boolean eigen Vector indicating the presence of a singularity on a vertex
  //   singularityIndex #V by 1 integer eigen Vector containing the singularity indices
  //)igl_Qu8mg5v7";
const char *__doc_igl_fit_rotations = R"igl_Qu8mg5v7(// Known issues: This seems to be implemented in Eigen/Geometry:
  // Eigen::umeyama
  //
  // FIT_ROTATIONS Given an input mesh and new positions find rotations for
  // every covariance matrix in a stack of covariance matrices
  //
  // Inputs:
  //   S  nr*dim by dim stack of covariance matrices
  //   single_precision  whether to use single precision (faster)
  // Outputs:
  //   R  dim by dim * nr list of rotations
  //)igl_Qu8mg5v7";
const char *__doc_igl_fit_rotations_planar = R"igl_Qu8mg5v7(// FIT_ROTATIONS Given an input mesh and new positions find 2D rotations for
  // every vertex that best maps its one ring to the new one ring
  //
  // Inputs:
  //   S  nr*dim by dim stack of covariance matrices, third column and every
  //   third row will be ignored
  // Outputs:
  //   R  dim by dim * nr list of rotations, third row and third column of each
  //   rotation will just be identity
  //)igl_Qu8mg5v7";
const char *__doc_igl_fit_rotations_SSE = R"igl_Qu8mg5v7(See fit_rotations_SSE for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_floor = R"igl_Qu8mg5v7(// Floor a given matrix to nearest integers
  //
  // Inputs:
  //   X  m by n matrix of scalars
  // Outputs:
  //   Y  m by n matrix of floored integers)igl_Qu8mg5v7";
const char *__doc_igl_forward_kinematics = R"igl_Qu8mg5v7(// Given a skeleton and a set of relative bone rotations compute absolute
  // rigid transformations for each bone.
  //
  // Inputs:
  //   C  #C by dim list of joint positions
  //   BE  #BE by 2 list of bone edge indices
  //   P  #BE list of parent indices into BE
  //   dQ  #BE list of relative rotations
  //   dT  #BE list of relative translations
  // Outputs:
  //   vQ  #BE list of absolute rotations
  //   vT  #BE list of absolute translations)igl_Qu8mg5v7";
const char *__doc_igl_gaussian_curvature = R"igl_Qu8mg5v7(// Compute discrete local integral gaussian curvature (angle deficit, without
  // averaging by local area).
  //
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigen Matrix of face (triangle) indices
  // Output:
  //   K  #V by 1 eigen Matrix of discrete gaussian curvature values
  //)igl_Qu8mg5v7";
const char *__doc_igl_get_seconds = R"igl_Qu8mg5v7(// Return the current time in seconds since program start
  //
  // Example:
  //    const auto & tictoc = []()
  //    {
  //      static double t_start = igl::get_seconds();
  //      double diff = igl::get_seconds()-t_start;
  //      t_start += diff;
  //      return diff;
  //    };
  //    tictoc();
  //    ... // part 1
  //    cout<<"part 1: "<<tictoc()<<endl;
  //    ... // part 2
  //    cout<<"part 2: "<<tictoc()<<endl;
  //    ... // etc)igl_Qu8mg5v7";
const char *__doc_igl_grad = R"igl_Qu8mg5v7(// Gradient of a scalar function defined on piecewise linear elements (mesh)
  // is constant on each triangle [tetrahedron] i,j,k:
  // grad(Xijk) = (Xj-Xi) * (Vi - Vk)^R90 / 2A + (Xk-Xi) * (Vj - Vi)^R90 / 2A
  // where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
  // i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of
  // 90 degrees
  //)igl_Qu8mg5v7";
const char *__doc_igl_harmonic = R"igl_Qu8mg5v7(// Compute k-harmonic weight functions "coordinates".
  //
  //
  // Inputs:
  //   V  #V by dim vertex positions
  //   F  #F by simplex-size list of element indices
  //   b  #b boundary indices into V
  //   bc #b by #W list of boundary values
  //   k  power of harmonic operation (1: harmonic, 2: biharmonic, etc)
  // Outputs:
  //   W  #V by #W list of weights
  //)igl_Qu8mg5v7";
const char *__doc_igl_hsv_to_rgb = R"igl_Qu8mg5v7(// Convert RGB to HSV
  //
  // Inputs:
  //   h  hue value (degrees: [0,360])
  //   s  saturation value ([0,1])
  //   v  value value ([0,1])
  // Outputs:
  //   r  red value ([0,1])
  //   g  green value ([0,1])
  //   b  blue value ([0,1]))igl_Qu8mg5v7";
const char *__doc_igl_internal_angles = R"igl_Qu8mg5v7(// Compute internal angles for a triangle mesh
  //
  // Inputs:
  //   V  #V by dim eigen Matrix of mesh vertex nD positions
  //   F  #F by poly-size eigen Matrix of face (triangle) indices
  // Output:
  //   K  #F by poly-size eigen Matrix of internal angles
  //     for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  //
  // Known Issues:
  //   if poly-size ≠ 3 then dim must equal 3.)igl_Qu8mg5v7";
const char *__doc_igl_internal_angles_using_squared_edge_lengths = R"igl_Qu8mg5v7(// Inputs:
  //   L_sq  #F by 3 list of squared edge lengths
  // Output:
  //   K  #F by poly-size eigen Matrix of internal angles
  //     for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  //
  // Note:
  //   Usage of internal_angles_using_squared_edge_lengths is preferred to internal_angles_using_squared_edge_lengths)igl_Qu8mg5v7";
const char *__doc_igl_internal_angles_using_edge_lengths = R"igl_Qu8mg5v7(// Inputs:
  //   L  #F by 3 list of edge lengths
  // Output:
  //   K  #F by poly-size eigen Matrix of internal angles
  //     for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  //
  // Note:
  //   Usage of internal_angles_using_squared_edge_lengths is preferred to internal_angles_using_squared_edge_lengths
  //   This function is deprecated and probably will be removed in future versions)igl_Qu8mg5v7";
const char *__doc_igl_invert_diag = R"igl_Qu8mg5v7(// Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   X  an m by n sparse matrix
  // Outputs:
  //   Y  an m by n sparse matrix)igl_Qu8mg5v7";
const char *__doc_igl_is_irregular_vertex = R"igl_Qu8mg5v7(// Determine if a vertex is irregular, i.e. it has more than 6 (triangles)
  // or 4 (quads) incident edges. Vertices on the boundary are ignored.
  //
  // Inputs:
  //   V  #V by dim list of vertex positions
  //   F  #F by 3[4] list of triangle[quads] indices
  // Returns #V vector of bools revealing whether vertices are singular
  //)igl_Qu8mg5v7";
const char *__doc_igl_jet = R"igl_Qu8mg5v7(// JET like MATLAB's jet
  //
  // Inputs:
  //   m  number of colors
  // Outputs:
  //   J  m by list of RGB colors between 0 and 1
  //
//#ifndef IGL_NO_EIGEN
//  void jet(const int m, Eigen::MatrixXd & J);
//#endif
  // Wrapper for directly computing [r,g,b] values for a given factor f between
  // 0 and 1
  //
  // Inputs:
  //   f  factor determining color value as if 0 was min and 1 was max
  // Outputs:
  //   r  red value
  //   g  green value
  //   b  blue value)igl_Qu8mg5v7";
const char *__doc_igl_lbs_matrix = R"igl_Qu8mg5v7(// LBS_MATRIX Linear blend skinning can be expressed by V' = M * T where V' is
  // a #V by dim matrix of deformed vertex positions (one vertex per row), M is a
  // #V by (dim+1)*#T (composed of weights and rest positions) and T is a
  // #T*(dim+1) by dim matrix of #T stacked transposed transformation matrices.
  // See equations (1) and (2) in "Fast Automatic Skinning Transformations"
  // [Jacobson et al 2012]
  //
  // Inputs:
  //   V  #V by dim list of rest positions
  //   W  #V+ by #T  list of weights
  // Outputs:
  //   M  #V by #T*(dim+1)
  //
  // In MATLAB:
  //   kron(ones(1,size(W,2)),[V ones(size(V,1),1)]).*kron(W,ones(1,size(V,2)+1)))igl_Qu8mg5v7";
const char *__doc_igl_lbs_matrix_column = R"igl_Qu8mg5v7(// LBS_MATRIX  construct a matrix that when multiplied against a column of
  // affine transformation entries computes new coordinates of the vertices
  //
  // I'm not sure it makes since that the result is stored as a sparse matrix.
  // The number of non-zeros per row *is* dependent on the number of mesh
  // vertices and handles.
  //
  // Inputs:
  //   V  #V by dim list of vertex rest positions
  //   W  #V by #handles list of correspondence weights
  // Output:
  //   M  #V * dim by #handles * dim * (dim+1) matrix such that
  //     new_V(:) = LBS(V,W,A) = reshape(M * A,size(V)), where A is a column
  //     vectors formed by the entries in each handle's dim by dim+1
  //     transformation matrix. Specifcally, A =
  //       reshape(permute(Astack,[3 1 2]),n*dim*(dim+1),1)
  //     or A = [Lxx;Lyx;Lxy;Lyy;tx;ty], and likewise for other dim
  //     if Astack(:,:,i) is the dim by (dim+1) transformation at handle i)igl_Qu8mg5v7";
const char *__doc_igl_local_basis = R"igl_Qu8mg5v7(// Compute a local orthogonal reference system for each triangle in the given mesh
  // Templates:
  //   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
  //   DerivedF derived from face indices matrix type: i.e. MatrixXi
  // Inputs:
  //   V  eigen matrix #V by 3
  //   F  #F by 3 list of mesh faces (must be triangles)
  // Outputs:
  //   B1 eigen matrix #F by 3, each vector is tangent to the triangle
  //   B2 eigen matrix #F by 3, each vector is tangent to the triangle and perpendicular to B1
  //   B3 eigen matrix #F by 3, normal of the triangle
  //
  // See also: adjacency_matrix)igl_Qu8mg5v7";
const char *__doc_igl_lscm = R"igl_Qu8mg5v7(// Compute a Least-squares conformal map parametrization (equivalently
  // derived in "Intrinsic Parameterizations of Surface Meshes" [Desbrun et al.
  // 2002] and "Least Squares Conformal Maps for Automatic Texture Atlas
  // Generation" [Lévy et al. 2002]), though this implementation follows the
  // derivation in: "Spectral Conformal Parameterization" [Mullen et al. 2008]
  // (note, this does **not** implement the Eigen-decomposition based method in
  // [Mullen et al. 2008], which is not equivalent). Input should be a manifold
  // mesh (also no unreferenced vertices) and "boundary" (fixed vertices) `b`
  // should contain at least two vertices per connected component.
  //
  // Inputs:
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of mesh faces (must be triangles)
  //   b  #b boundary indices into V
  //   bc #b by 3 list of boundary values
  // Outputs:
  //   UV #V by 2 list of 2D mesh vertex positions in UV space
  // Returns true only on solver success.
  //)igl_Qu8mg5v7";
const char *__doc_igl_map_vertices_to_circle = R"igl_Qu8mg5v7(// Map the vertices whose indices are in a given boundary loop (bnd) on the
  // unit circle with spacing proportional to the original boundary edge
  // lengths.
  //
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   b  #W list of vertex ids
  // Outputs:
  //   UV   #W by 2 list of 2D position on the unit circle for the vertices in b)igl_Qu8mg5v7";
const char *__doc_igl_massmatrix = R"igl_Qu8mg5v7(// Constructs the mass (area) matrix for a given mesh (V,F).
  //
  // Templates:
  //   DerivedV  derived type of eigen matrix for V (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. derived from
  //     MatrixXi)
  //   Scalar  scalar type for eigen sparse matrix (e.g. double)
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by simplex_size list of mesh faces (must be triangles)
  //   type  one of the following ints:
  //     MASSMATRIX_TYPE_BARYCENTRIC  barycentric
  //     MASSMATRIX_TYPE_VORONOI voronoi-hybrid {default}
  //     MASSMATRIX_TYPE_FULL full {not implemented}
  // Outputs:
  //   M  #V by #V mass matrix
  //
  // See also: adjacency_matrix
  //)igl_Qu8mg5v7";
const char *__doc_igl_min_quad_with_fixed_precompute = R"igl_Qu8mg5v7(// Known Bugs: rows of Aeq **should probably** be linearly independent.
  // During precomputation, the rows of a Aeq are checked via QR. But in case
  // they're not then resulting probably will no longer be sparse: it will be
  // slow.
  //
  // MIN_QUAD_WITH_FIXED Minimize a quadratic energy of the form
  //
  // trace( 0.5*Z'*A*Z + Z'*B + constant )
  //
  // subject to
  //
  //   Z(known,:) = Y, and
  //   Aeq*Z = Beq
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   A  n by n matrix of quadratic coefficients
  //   known list of indices to known rows in Z
  //   Y  list of fixed values corresponding to known rows in Z
  //   Aeq  m by n list of linear equality constraint coefficients
  //   pd flag specifying whether A(unknown,unknown) is positive definite
  // Outputs:
  //   data  factorization struct with all necessary information to solve
  //     using min_quad_with_fixed_solve
  // Returns true on success, false on error
  //
  // Benchmark: For a harmonic solve on a mesh with 325K facets, matlab 2.2
  // secs, igl/min_quad_with_fixed.h 7.1 secs
  //)igl_Qu8mg5v7";
const char *__doc_igl_min_quad_with_fixed_solve = R"igl_Qu8mg5v7(// Solves a system previously factored using min_quad_with_fixed_precompute
  //
  // Template:
  //   T  type of sparse matrix (e.g. double)
  //   DerivedY  type of Y (e.g. derived from VectorXd or MatrixXd)
  //   DerivedZ  type of Z (e.g. derived from VectorXd or MatrixXd)
  // Inputs:
  //   data  factorization struct with all necessary precomputation to solve
  //   B  n by k column of linear coefficients
  //   Y  b by k list of constant fixed values
  //   Beq  m by k list of linear equality constraint constant values
  // Outputs:
  //   Z  n by k solution
  //   sol  #unknowns+#lagrange by k solution to linear system
  // Returns true on success, false on error)igl_Qu8mg5v7";
const char *__doc_igl_min_quad_with_fixed = R"igl_Qu8mg5v7(See min_quad_with_fixed for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_normalize_row_lengths = R"igl_Qu8mg5v7(// Obsolete: just use A.rowwise().normalize() or B=A.rowwise().normalized();
  //
  // Normalize the rows in A so that their lengths are each 1 and place the new
  // entries in B
  // Inputs:
  //   A  #rows by k input matrix
  // Outputs:
  //   B  #rows by k input matrix, can be the same as A)igl_Qu8mg5v7";
const char *__doc_igl_normalize_row_sums = R"igl_Qu8mg5v7(// Normalize the rows in A so that their sums are each 1 and place the new
  // entries in B
  // Inputs:
  //   A  #rows by k input matrix
  // Outputs:
  //   B  #rows by k input matrix, can be the same as A
  //
  // Note: This is just calling an Eigen one-liner.)igl_Qu8mg5v7";
const char *__doc_igl_parula = R"igl_Qu8mg5v7(// PARULA like MATLAB's parula
  //
  // Inputs:
  //   m  number of colors
  // Outputs:
  //   J  m by list of RGB colors between 0 and 1
  //
  // Wrapper for directly computing [r,g,b] values for a given factor f between
  // 0 and 1
  //
  // Inputs:
  //   f  factor determining color value as if 0 was min and 1 was max
  // Outputs:
  //   r  red value
  //   g  green value
  //   b  blue value)igl_Qu8mg5v7";
const char *__doc_igl_per_corner_normals = R"igl_Qu8mg5v7(// Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  //   corner_threshold  threshold in degrees on sharp angles
  // Output:
  //   CN  #F*3 by 3 eigen Matrix of mesh vertex 3D normals, where the normal
  //     for corner F(i,j) is at CN(i*3+j,:) )igl_Qu8mg5v7";
const char *__doc_igl_per_edge_normals = R"igl_Qu8mg5v7(// Compute face normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigen Matrix of face (triangle) indices
  //   weight  weighting type
  //   FN  #F by 3 matrix of 3D face normals per face
  // Output:
  //   N  #2 by 3 matrix of mesh edge 3D normals per row
  //   E  #E by 2 matrix of edge indices per row
  //   EMAP  #E by 1 matrix of indices from all edges to E
  //)igl_Qu8mg5v7";
const char *__doc_igl_per_face_normals = R"igl_Qu8mg5v7(// Compute face normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigen Matrix of face (triangle) indices
  //   Z  3 vector normal given to faces with degenerate normal.
  // Output:
  //   N  #F by 3 eigen Matrix of mesh face (triangle) 3D normals
  //
  // Example:
  //   // Give degenerate faces (1/3,1/3,1/3)^0.5
  //   per_face_normals(V,F,Vector3d(1,1,1).normalized(),N);)igl_Qu8mg5v7";
const char *__doc_igl_per_face_normals_stable = R"igl_Qu8mg5v7(// Special version where order of face indices is guaranteed not to effect
  // output.)igl_Qu8mg5v7";
const char *__doc_igl_per_vertex_normals = R"igl_Qu8mg5v7(// Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  //   weighting  Weighting type
  // Output:
  //   N  #V by 3 eigen Matrix of mesh vertex 3D normals)igl_Qu8mg5v7";
const char *__doc_igl_planarize_quad_mesh = R"igl_Qu8mg5v7(// Inputs:
  //   Vin        #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F          #F by 4 eigen Matrix of face (quad) indices
  //   maxIter    maximum numbers of iterations
  //   threshold  minimum allowed threshold for non-planarity
  // Output:
  //   Vout       #V by 3 eigen Matrix of planar mesh vertex 3D positions
  //)igl_Qu8mg5v7";
const char *__doc_igl_png_readPNG = R"igl_Qu8mg5v7(// Read an image from a .png file into 4 memory buffers
    //
    // Input:
    //  png_file  path to .png file
    // Output:
    //  R,G,B,A texture channels
    // Returns true on success, false on failure
    //)igl_Qu8mg5v7";
const char *__doc_igl_png_writePNG = R"igl_Qu8mg5v7(// Writes an image to a png file
    //
    // Input:
    //  R,G,B,A texture channels
    // Output:
    //  png_file  path to .png file
    // Returns true on success, false on failure
    //)igl_Qu8mg5v7";
const char *__doc_igl_point_mesh_squared_distance = R"igl_Qu8mg5v7(// Compute distances from a set of points P to a triangle mesh (V,F)
  //
  // Inputs:
  //   P  #P by 3 list of query point positions
  //   V  #V by 3 list of vertex positions
  //   Ele  #Ele by (3|2|1) list of (triangle|edge|point) indices
  // Outputs:
  //   sqrD  #P list of smallest squared distances
  //   I  #P list of primitive indices corresponding to smallest distances
  //   C  #P by 3 list of closest points
  //
  // Known bugs: This only computes distances to given primitivess. So
  // unreferenced vertices are ignored. However, degenerate primitives are
  // handled correctly: triangle [1 2 2] is treated as a segment [1 2], and
  // triangle [1 1 1] is treated as a point. So one _could_ add extra
  // combinatorially degenerate rows to Ele for all unreferenced vertices to
  // also get distances to points.)igl_Qu8mg5v7";
const char *__doc_igl_polar_svd = R"igl_Qu8mg5v7(// Computes the polar decomposition (R,T) of a matrix A using SVD singular
  // value decomposition
  //
  // Inputs:
  //   A  3 by 3 matrix to be decomposed
  // Outputs:
  //   R  3 by 3 rotation matrix part of decomposition (**always rotataion**)
  //   T  3 by 3 stretch matrix part of decomposition
  //   U  3 by 3 left-singular vectors
  //   S  3 by 1 singular values
  //   V  3 by 3 right-singular vectors
  //
  //)igl_Qu8mg5v7";
const char *__doc_igl_principal_curvature = R"igl_Qu8mg5v7(// Compute the principal curvature directions and magnitude of the given triangle mesh
  //   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
  //   DerivedF derived from face indices matrix type: i.e. MatrixXi
  // Inputs:
  //   V       eigen matrix #V by 3
  //   F       #F by 3 list of mesh faces (must be triangles)
  //   radius  controls the size of the neighbourhood used, 1 = average edge length
  //
  // Outputs:
  //   PD1 #V by 3 maximal curvature direction for each vertex.
  //   PD2 #V by 3 minimal curvature direction for each vertex.
  //   PV1 #V by 1 maximal curvature value for each vertex.
  //   PV2 #V by 1 minimal curvature value for each vertex.
  //
  // See also: average_onto_faces, average_onto_vertices
  //
  // This function has been developed by: Nikolas De Giorgis, Luigi Rocca and Enrico Puppo.
  // The algorithm is based on:
  // Efficient Multi-scale Curvature and Crease Estimation
  // Daniele Panozzo, Enrico Puppo, Luigi Rocca
  // GraVisMa, 2010)igl_Qu8mg5v7";
const char *__doc_igl_quad_planarity = R"igl_Qu8mg5v7(// Compute planarity of the faces of a quad mesh
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 4 eigen Matrix of face (quad) indices
  // Output:
  //   P  #F by 1 eigen Matrix of mesh face (quad) planarities
  //)igl_Qu8mg5v7";
const char *__doc_igl_randperm = R"igl_Qu8mg5v7(// Like matlab's randperm(n) but minus 1
  //
  // Inputs:
  //   n  number of elements
  // Outputs:
  //   I  n list of rand permutation of 0:n-1)igl_Qu8mg5v7";
const char *__doc_igl_readDMAT = R"igl_Qu8mg5v7(See readDMAT for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_readMESH = R"igl_Qu8mg5v7(// load a tetrahedral volume mesh from a .mesh file
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Input:
  //   mesh_file_name  path of .mesh file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   T  #T list of tet indices into vertex positions
  //   F  #F list of face indices into vertex positions
  //
  // Known bugs: Holes and regions are not supported)igl_Qu8mg5v7";
const char *__doc_igl_readOBJ = R"igl_Qu8mg5v7(// Read a mesh from an ascii obj file, filling in vertex positions, normals
  // and texture coordinates. Mesh may have faces of any number of degree
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Inputs:
  //  str  path to .obj file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   TC  double matrix of texture coordinats #TC by 2
  //   N  double matrix of corner normals #N by 3
  //   F  #F list of face indices into vertex positions
  //   FTC  #F list of face indices into vertex texture coordinates
  //   FN  #F list of face indices into vertex normals
  // Returns true on success, false on errors)igl_Qu8mg5v7";
const char *__doc_igl_readOFF = R"igl_Qu8mg5v7(// Read a mesh from an ascii OFF file, filling in vertex positions, normals
  // and texture coordinates. Mesh may have faces of any number of degree
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Inputs:
  //  str  path to .obj file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   F  #F list of face indices into vertex positions
  //   N  list of vertex normals #V by 3
  //   C  list of rgb color values per vertex #V by 3
  // Returns true on success, false on errors)igl_Qu8mg5v7";
const char *__doc_igl_readTGF = R"igl_Qu8mg5v7(// READTGF
  //
  // [V,E,P,BE,CE,PE] = readTGF(filename)
  //
  // Read a graph from a .tgf file
  //
  // Input:
  //  filename  .tgf file name
  // Output:
  //  V  # vertices by 3 list of vertex positions
  //  E  # edges by 2 list of edge indices
  //  P  # point-handles list of point handle indices
  //  BE # bone-edges by 2 list of bone-edge indices
  //  CE # cage-edges by 2 list of cage-edge indices
  //  PE # pseudo-edges by 2 list of pseudo-edge indices
  //
  // Assumes that graph vertices are 3 dimensional)igl_Qu8mg5v7";
const char *__doc_igl_read_triangle_mesh = R"igl_Qu8mg5v7(// read mesh from an ascii file with automatic detection of file format.
  // supported: obj, off, stl, wrl, ply, mesh)
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Inputs:
  //   str  path to file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   F  eigen int matrix #F by 3
  // Returns true iff success)igl_Qu8mg5v7";
const char *__doc_igl_remove_duplicate_vertices = R"igl_Qu8mg5v7(// REMOVE_DUPLICATE_VERTICES Remove duplicate vertices upto a uniqueness
  // tolerance (epsilon)
  //
  // Inputs:
  //   V  #V by dim list of vertex positions
  //   epsilon  uniqueness tolerance (significant digit), can probably think of
  //     this as a tolerance on L1 distance
  // Outputs:
  //   SV  #SV by dim new list of vertex positions
  //   SVI #V by 1 list of indices so SV = V(SVI,:)
  //   SVJ #SV by 1 list of indices so V = SV(SVJ,:)
  //
  // Example:
  //   % Mesh in (V,F)
  //   [SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-7);
  //   % remap faces
  //   SF = SVJ(F);
  //)igl_Qu8mg5v7";
const char *__doc_igl_rotate_vectors = R"igl_Qu8mg5v7(// Rotate the vectors V by A radiants on the tangent plane spanned by B1 and
  // B2
  //
  // Inputs:
  //   V     #V by 3 eigen Matrix of vectors
  //   A     #V eigen vector of rotation angles or a single angle to be applied
  //     to all vectors
  //   B1    #V by 3 eigen Matrix of base vector 1
  //   B2    #V by 3 eigen Matrix of base vector 2
  //
  // Output:
  //   Returns the rotated vectors
  //)igl_Qu8mg5v7";
const char *__doc_igl_setdiff = R"igl_Qu8mg5v7(// Set difference of elements of matrices
  //
  // Inputs:
  //   A  m-long vector of indices
  //   B  n-long vector of indices
  // Outputs:
  //   C  (k<=m)-long vector of unique elements appearing in A but not in B
  //   IA  (k<=m)-long list of indices into A so that C = A(IA)
  //)igl_Qu8mg5v7";
const char *__doc_igl_shape_diameter_function = R"igl_Qu8mg5v7(// Compute shape diamater function per given point. In the parlence of the
  // paper "Consistent Mesh Partitioning and Skeletonisation using the Shape
  // Diameter Function" [Shapiro et al. 2008], this implementation uses a 180°
  // cone and a _uniform_ average (_not_ a average weighted by inverse angles).
  //
  // Inputs:
  //    shoot_ray  function handle that outputs hits of a given ray against a
  //      mesh (embedded in function handles as captured variable/data)
  //    P  #P by 3 list of origin points
  //    N  #P by 3 list of origin normals
  // Outputs:
  //    S  #P list of shape diamater function values between bounding box
  //    diagonal (perfect sphere) and 0 (perfect needle hook)
  //)igl_Qu8mg5v7";
const char *__doc_igl_signed_distance = R"igl_Qu8mg5v7(// Computes signed distance to a mesh
  //
  // Inputs:
  //   P  #P by 3 list of query point positions
  //   V  #V by 3 list of vertex positions
  //   F  #F by ss list of triangle indices, ss should be 3 unless sign_type ==
  //     SIGNED_DISTANCE_TYPE_UNSIGNED
  //   sign_type  method for computing distance _sign_ S
  // Outputs:
  //   S  #P list of smallest signed distances
  //   I  #P list of facet indices corresponding to smallest distances
  //   C  #P by 3 list of closest points
  //   N  #P by 3 list of closest normals (only set if
  //     sign_type=SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
  //
  // Known bugs: This only computes distances to triangles. So unreferenced
  // vertices and degenerate triangles are ignored.)igl_Qu8mg5v7";
const char *__doc_igl_signed_distance_pseudonormal = R"igl_Qu8mg5v7(// Computes signed distance to mesh
  //
  // Inputs:
  //   tree  AABB acceleration tree (see AABB.h)
  //   F  #F by 3 list of triangle indices
  //   FN  #F by 3 list of triangle normals
  //   VN  #V by 3 list of vertex normals (ANGLE WEIGHTING)
  //   EN  #E by 3 list of edge normals (UNIFORM WEIGHTING)
  //   EMAP  #F*3 mapping edges in F to E
  //   q  Query point
  // Returns signed distance to mesh
  //)igl_Qu8mg5v7";
const char *__doc_igl_signed_distance_winding_number = R"igl_Qu8mg5v7(// Inputs:
  //   tree  AABB acceleration tree (see cgal/point_mesh_squared_distance.h)
  //   hier  Winding number evaluation hierarchy
  //   q  Query point
  // Returns signed distance to mesh)igl_Qu8mg5v7";
const char *__doc_igl_slice = R"igl_Qu8mg5v7(// Act like the matlab X(row_indices,col_indices) operator, where
  // row_indices, col_indices are non-negative integer indices.
  //
  // Inputs:
  //   X  m by n matrix
  //   R  list of row indices
  //   C  list of column indices
  // Output:
  //   Y  #R by #C matrix
  //
  // See also: slice_mask)igl_Qu8mg5v7";
const char *__doc_igl_slice_into = R"igl_Qu8mg5v7(// Act like the matlab Y(row_indices,col_indices) = X
  //
  // Inputs:
  //   X  xm by xn rhs matrix
  //   R  list of row indices
  //   C  list of column indices
  //   Y  ym by yn lhs matrix
  // Output:
  //   Y  ym by yn lhs matrix, same as input but Y(R,C) = X)igl_Qu8mg5v7";
const char *__doc_igl_slice_mask = R"igl_Qu8mg5v7(// Act like the matlab X(row_mask,col_mask) operator, where
  // row_mask, col_mask are non-negative integer indices.
  //
  // Inputs:
  //   X  m by n matrix
  //   R  m list of row bools
  //   C  n list of column bools
  // Output:
  //   Y  #trues-in-R by #trues-in-C matrix
  //
  // See also: slice_mask)igl_Qu8mg5v7";
const char *__doc_igl_marching_tets = R"igl_Qu8mg5v7(// SLICE_TETS Slice through a tet mesh (V,T) along a given plane (via its
  // implicit equation).
  //
  // Inputs:
  //   V  #V by 3 list of tet mesh vertices
  //   T  #T by 4 list of tet indices into V
  //   plane  list of 4 coefficients in the plane equation: [x y z 1]'*plane = 0
  //   Optional:
  //     'Manifold' followed by whether to stitch together triangles into a
  //       manifold mesh {true}: results in more compact U but slightly slower.
  // Outputs:
  //   U  #U by 3 list of triangle mesh vertices along slice
  //   G  #G by 3 list of triangles indices into U
  //   J  #G list of indices into T revealing from which tet each faces comes
  //   BC  #U by #V list of barycentric coordinates (or more generally: linear
  //     interpolation coordinates) so that U = BC*V
  // )igl_Qu8mg5v7";
const char *__doc_igl_sortrows = R"igl_Qu8mg5v7(// Act like matlab's [Y,I] = sortrows(X)
  //
  // Templates:
  //   DerivedX derived scalar type, e.g. MatrixXi or MatrixXd
  //   DerivedI derived integer type, e.g. MatrixXi
  // Inputs:
  //   X  m by n matrix whose entries are to be sorted
  //   ascending  sort ascending (true, matlab default) or descending (false)
  // Outputs:
  //   Y  m by n matrix whose entries are sorted (**should not** be same
  //     reference as X)
  //   I  m list of indices so that
  //     Y = X(I,:);)igl_Qu8mg5v7";
const char *__doc_igl_streamlines_init = R"igl_Qu8mg5v7(// Given a mesh and a field the function computes the /data/ necessary for tracing the field'
    // streamlines, and creates the initial /state/ for the tracing.
    // Inputs:
    //   V             #V by 3 list of mesh vertex coordinates
    //   F             #F by 3 list of mesh faces
    //   temp_field    #F by 3n list of the 3D coordinates of the per-face vectors
    //                    (n-degrees stacked horizontally for each triangle)
    //   treat_as_symmetric
    //              if true, adds n symmetry directions to the field (N = 2n). Else N = n
    //   percentage    [0-1] percentage of faces sampled
    // Outputs:
    //   data          struct containing topology information of the mesh and field
    //   state         struct containing the state of the tracing)igl_Qu8mg5v7";
const char *__doc_igl_streamlines_next = R"igl_Qu8mg5v7(// The function computes the next state for each point in the sample
    //   V             #V by 3 list of mesh vertex coordinates
    //   F             #F by 3 list of mesh faces
    //   data          struct containing topology information
    //   state         struct containing the state of the tracing)igl_Qu8mg5v7";
const char *__doc_igl_triangle_triangle_adjacency = R"igl_Qu8mg5v7(// Constructs the triangle-triangle adjacency matrix for a given
  // mesh (V,F).
  //
  // Templates:
  //   Scalar derived type of eigen matrix for V (e.g. derived from
  //     MatrixXd)
  //   Index  derived type of eigen matrix for F (e.g. derived from
  //     MatrixXi)
  // Inputs:
  //   F  #F by simplex_size list of mesh faces (must be triangles)
  // Outputs:
  //   TT   #F by #3 adjacent matrix, the element i,j is the id of the triangle adjacent to the j edge of triangle i
  //   TTi  #F by #3 adjacent matrix, the element i,j is the id of edge of the triangle TT(i,j) that is adjacent with triangle i
  // NOTE: the first edge of a triangle is [0,1] the second [1,2] and the third [2,3].
  //       this convention is DIFFERENT from cotmatrix_entries.h
  // Known bug: this should not need to take V as input.)igl_Qu8mg5v7";
const char *__doc_igl_triangle_triangle_adjacency_preprocess = R"igl_Qu8mg5v7(// Preprocessing)igl_Qu8mg5v7";
const char *__doc_igl_triangle_triangle_adjacency_extractTT = R"igl_Qu8mg5v7(// Extract the face adjacencies)igl_Qu8mg5v7";
const char *__doc_igl_triangle_triangle_adjacency_extractTTi = R"igl_Qu8mg5v7(// Extract the face adjacencies indices (needed for fast traversal))igl_Qu8mg5v7";
const char *__doc_igl_triangle_triangulate = R"igl_Qu8mg5v7(// Triangulate the interior of a polygon using the triangle library.
    //
    // Inputs:
    //   V #V by 2 list of 2D vertex positions
    //   E #E by 2 list of vertex ids forming unoriented edges of the boundary of the polygon
    //   H #H by 2 coordinates of points contained inside holes of the polygon
    //   flags  string of options pass to triangle (see triangle documentation)
    // Outputs:
    //   V2  #V2 by 2  coordinates of the vertives of the generated triangulation
    //   F2  #F2 by 3  list of indices forming the faces of the generated triangulation
    //)igl_Qu8mg5v7";
const char *__doc_igl_unique = R"igl_Qu8mg5v7(// Act like matlab's [C,IA,IC] = unique(X)
  //
  // Templates:
  //   T  comparable type T
  // Inputs:
  //   A  #A vector of type T
  // Outputs:
  //   C  #C vector of unique entries in A
  //   IA  #C index vector so that C = A(IA);
  //   IC  #A index vector so that A = C(IC);)igl_Qu8mg5v7";
const char *__doc_igl_unique_rows = R"igl_Qu8mg5v7(// Act like matlab's [C,IA,IC] = unique(X,'rows')
  //
  // Templates:
  //   DerivedA derived scalar type, e.g. MatrixXi or MatrixXd
  //   DerivedIA derived integer type, e.g. MatrixXi
  //   DerivedIC derived integer type, e.g. MatrixXi
  // Inputs:
  //   A  m by n matrix whose entries are to unique'd according to rows
  // Outputs:
  //   C  #C vector of unique rows in A
  //   IA  #C index vector so that C = A(IA,:);
  //   IC  #A index vector so that A = C(IC,:);)igl_Qu8mg5v7";
const char *__doc_igl_unproject_onto_mesh = R"igl_Qu8mg5v7(// Unproject a screen location (using current opengl viewport, projection, and
  // model view) to a 3D position _onto_ a given mesh, if the ray through the
  // given screen location (x,y) _hits_ the mesh.
  //
  // Inputs:
  //    pos        screen space coordinates
  //    model      model matrix
  //    proj       projection matrix
  //    viewport   vieweport vector
  //    V   #V by 3 list of mesh vertex positions
  //    F   #F by 3 list of mesh triangle indices into V
  // Outputs:
  //    fid  id of the first face hit
  //    bc  barycentric coordinates of hit
  // Returns true if there's a hit)igl_Qu8mg5v7";
const char *__doc_igl_upsample = R"igl_Qu8mg5v7(// Subdivide without moving vertices: Given the triangle mesh [V, F],
  // where n_verts = V.rows(), computes newV and a sparse matrix S s.t.
  // [newV, newF] is the subdivided mesh where newV = S*V.
  //
  // Inputs:
  //   n_verts  an integer (number of mesh vertices)
  //   F  an m by 3 matrix of integers of triangle faces
  // Outputs:
  //   S  a sparse matrix (will become the subdivision matrix)
  //   newF  a matrix containing the new faces)igl_Qu8mg5v7";
const char *__doc_igl_winding_number = R"igl_Qu8mg5v7(// WINDING_NUMBER Compute the sum of solid angles of a triangle/tetrahedron
  // described by points (vectors) V
  //
  // Templates:
  //   dim  dimension of input
  // Inputs:
  //  V  n by 3 list of vertex positions
  //  F  #F by 3 list of triangle indices, minimum index is 0
  //  O  no by 3 list of origin positions
  // Outputs:
  //  S  no by 1 list of winding numbers
  //
  // 3d)igl_Qu8mg5v7";
const char *__doc_igl_winding_number_3 = R"igl_Qu8mg5v7(// Inputs:
  //   V  pointer to array containing #V by 3 vertex positions along rows,
  //     given in column major order
  //   n  number of mesh vertices
  //   F  pointer to array containing #F by 3 face indices along rows,
  //     given in column major order
  //   m  number of faces
  //   O  pointer to array containing #O by 3 query positions along rows,
  //     given in column major order
  //   no  number of origins
  // Outputs:
  //   S  no by 1 list of winding numbers)igl_Qu8mg5v7";
const char *__doc_igl_winding_number_2 = R"igl_Qu8mg5v7(//// Only one evaluation origin
  //template <typename DerivedF>
  //IGL_INLINE void winding_number_3(
  //  const double * V,
  //  const int n,
  //  const DerivedF * F,
  //  const int m,
  //  const double * O,
  //  double * S);
  // 2d)igl_Qu8mg5v7";
const char *__doc_igl_writeMESH = R"igl_Qu8mg5v7(// save a tetrahedral volume mesh to a .mesh file
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be cast as double)
  //   Index  type for indices (will be cast to int)
  // Input:
  //   mesh_file_name  path of .mesh file
  //   V  double matrix of vertex positions  #V by 3
  //   T  #T list of tet indices into vertex positions
  //   F  #F list of face indices into vertex positions
  //
  // Known bugs: Holes and regions are not supported)igl_Qu8mg5v7";
const char *__doc_igl_writeOBJ = R"igl_Qu8mg5v7(// Write a mesh in an ascii obj file
  // Inputs:
  //   str  path to outputfile
  //   V  #V by 3 mesh vertex positions
  //   F  #F by 3|4 mesh indices into V
  //   CN #CN by 3 normal vectors
  //   FN  #F by 3|4 corner normal indices into CN
  //   TC  #TC by 2|3 texture coordinates
  //   FTC #F by 3|4 corner texture coord indices into TC
  // Returns true on success, false on error
  //
  // Known issues: Horrifyingly, this does not have the same order of
  // parameters as readOBJ.)igl_Qu8mg5v7";
const char *__doc_igl_writePLY = R"igl_Qu8mg5v7(// Write a mesh in an ascii ply file
  // Inputs:
  //   str  path to outputfile
  //   V  #V by 3 mesh vertex positions
  //   F  #F by 3 mesh indices into V
  //   N  #V by 3 normal vectors
  //   UV #V by 2 texture coordinates
  // Returns true on success, false on error)igl_Qu8mg5v7";
const char *__doc_igl_readPLY= R"igl_Qu8mg5v7(// Read a mesh from an ascii ply file, filling in vertex positions,
  // mesh indices, normals and texture coordinates
  // Inputs:
  //  str path to .obj file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   F  #F list of face indices into vertex positions
  //   N  double matrix of corner normals #N by 3
  //   UV #V by 2 texture coordinates
  // Returns true on success, false on errors)igl_Qu8mg5v7";
const char *__doc_igl_seam_edges=R"igl_Qu8mg5v7(// Finds all UV-space boundaries of a mesh.
  //
  // Inputs:
  //   V  #V by dim list of positions of the input mesh.
  //   TC  #TC by 2 list of 2D texture coordinates of the input mesh
  //   F  #F by 3 list of triange indices into V representing a
  //     manifold-with-boundary triangle mesh
  //   FTC  #F by 3 list of indices into TC for each corner
  // Outputs:
  //   seams  Edges where the forwards and backwards directions have different
  //     texture coordinates, as a #seams-by-4 matrix of indices. Each row is
  //     organized as [ forward_face_index, forward_face_vertex_index,
  //     backwards_face_index, backwards_face_vertex_index ] such that one side
  //     of the seam is the edge:
  //         F[ seams( i, 0 ), seams( i, 1 ) ], F[ seams( i, 0 ), (seams( i, 1 ) + 1) % 3 ]
  //     and the other side is the edge:
  //         F[ seams( i, 2 ), seams( i, 3 ) ], F[ seams( i, 2 ), (seams( i, 3 ) + 1) % 3 ]
  //   boundaries  Edges with only one incident triangle, as a #boundaries-by-2
  //     matrix of indices. Each row is organized as 
  //         [ face_index, face_vertex_index ]
  //     such that the edge is:
  //         F[ boundaries( i, 0 ), boundaries( i, 1 ) ], F[ boundaries( i, 0 ), (boundaries( i, 1 ) + 1) % 3 ]
  //   foldovers  Edges where the two incident triangles fold over each other
  //     in UV-space, as a #foldovers-by-4 matrix of indices.
  //     Each row is organized as [ forward_face_index, forward_face_vertex_index,
  //     backwards_face_index, backwards_face_vertex_index ]
  //     such that one side of the foldover is the edge:
  //       F[ foldovers( i, 0 ), foldovers( i, 1 ) ], F[ foldovers( i, 0 ), (foldovers( i, 1 ) + 1) % 3 ]
  //     and the other side is the edge:
  //       F[ foldovers( i, 2 ), foldovers( i, 3 ) ], F[ foldovers( i, 2 ), (foldovers( i, 3 ) + 1) % 3 ])igl_Qu8mg5v7";
