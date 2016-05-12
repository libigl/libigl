const char *__doc_igl_principal_curvature = R"igl_Qu8mg5v7(// Compute the principal curvature directions and magnitude of the given triangle mesh
  //   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
  //   DerivedF derived from face indices matrix type: i.e. MatrixXi
  // Inputs:
  //   V       eigen matrix #V by 3
  //   F       #F by 3 list of mesh faces (must be triangles)
  //   radius  controls the size of the neighbourhood used, 1 = average edge lenght
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
  //
  // Known bugs: off by 1e-16 on regular grid. I think its a problem of
  // arithmetic order in cotmatrix_entries.h: C(i,e) = (arithmetic)/dblA/4)igl_Qu8mg5v7";
const char *__doc_igl_floor = R"igl_Qu8mg5v7(// Floor a given matrix to nearest integers 
  //
  // Inputs:
  //   X  m by n matrix of scalars
  // Outputs:
  //   Y  m by n matrix of floored integers)igl_Qu8mg5v7";
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
const char *__doc_igl_readOFF = R"igl_Qu8mg5v7(// Read a mesh from an ascii obj file, filling in vertex positions, normals
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
  //   TC  double matrix of texture coordinats #TC by 2
  //   FTC  #F list of face indices into vertex texture coordinates
  //   N  double matrix of corner normals #N by 3
  //   FN  #F list of face indices into vertex normals
  // Returns true on success, false on errors)igl_Qu8mg5v7";
const char *__doc_igl_per_vertex_normals = R"igl_Qu8mg5v7(// Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  //   weighting  Weighting type
  // Output:
  //   N  #V by 3 eigen Matrix of mesh vertex 3D normals)igl_Qu8mg5v7";
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
const char *__doc_igl_barycenter = R"igl_Qu8mg5v7(// Computes the barycenter of every simplex
  //
  // Inputs:
  //   V  #V x dim matrix of vertex coordinates
  //   F  #F x simplex_size  matrix of indices of simplex corners into V
  // Output:
  //   BC  #F x dim matrix of 3d vertices
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
const char *__doc_igl_eigs = R"igl_Qu8mg5v7(See eigs for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_per_corner_normals = R"igl_Qu8mg5v7(// Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  //   corner_threshold  threshold in degrees on sharp angles
  // Output:
  //   CN  #F*3 by 3 eigen Matrix of mesh vertex 3D normals, where the normal
  //     for corner F(i,j) is at CN(i*3+j,:) )igl_Qu8mg5v7";
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
const char *__doc_igl_gaussian_curvature = R"igl_Qu8mg5v7(// Compute discrete local integral gaussian curvature (angle deficit, without
  // averaging by local area).
  //
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigen Matrix of face (triangle) indices
  // Output:
  //   K  #V by 1 eigen Matrix of discrete gaussian curvature values
  //)igl_Qu8mg5v7";
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
const char *__doc_igl_find_cross_field_singularities = R"igl_Qu8mg5v7(// Inputs:
  //   V                #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F                #F by 3 eigen Matrix of face (quad) indices
  //   Handle_MMatch    #F by 3 eigen Matrix containing the integer missmatch of the cross field
  //                    across all face edges
  // Output:
  //   isSingularity    #V by 1 boolean eigen Vector indicating the presence of a singularity on a vertex
  //   singularityIndex #V by 1 integer eigen Vector containing the singularity indices
  //)igl_Qu8mg5v7";
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
const char *__doc_igl_setdiff = R"igl_Qu8mg5v7(// Set difference of elements of matrices
  //
  // Inputs:
  //   A  m-long vector of indices
  //   B  n-long vector of indices
  // Outputs:
  //   C  (k<=m)-long vector of unique elements appearing in A but not in B
  //   IA  (k<=m)-long list of indices into A so that C = A(IA)
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
const char *__doc_igl_map_vertices_to_circle = R"igl_Qu8mg5v7(// Map the vertices whose indices are in a given boundary loop (bnd) on the
  // unit circle with spacing proportional to the original boundary edge
  // lengths.
  //
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   b  #W list of vertex ids
  // Outputs:
  //   UV   #W by 2 list of 2D position on the unit circle for the vertices in b)igl_Qu8mg5v7";
const char *__doc_igl_writeOBJ = R"igl_Qu8mg5v7(// Write a mesh in an ascii obj file
  // Inputs:
  //   str  path to outputfile
  //   V  #V by 3 mesh vertex positions
  //   F  #F by 3|4 mesh indices into V
  //   CN #CN by 3 normal vectors
  //   FN  #F by 3|4 corner normal indices into CN
  //   TC  #TC by 2|3 texture coordinates
  //   FTC #F by 3|4 corner texture coord indices into TC
  // Returns true on success, false on error)igl_Qu8mg5v7";
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
const char *__doc_igl_cut_mesh_from_singularities = R"igl_Qu8mg5v7(// Given a mesh (V,F) and the integer mismatch of a cross field per edge
  // (MMatch), finds the cut_graph connecting the singularities (seams) and the
  // degree of the singularities singularity_index
  //
  // Input:
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of faces
  //   MMatch  #F by 3 list of per corner integer mismatch
  // Outputs:
  //   seams  #F by 3 list of per corner booleans that denotes if an edge is a
  //     seam or not
  //)igl_Qu8mg5v7";
const char *__doc_igl_readDMAT = R"igl_Qu8mg5v7(See readDMAT for the documentation.)igl_Qu8mg5v7";
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
const char *__doc_igl_min_quad_with_fixed_precompute = R"igl_Qu8mg5v7(// Known Bugs: rows of Aeq **should probably** be linearly independent.
  // During precomputation, the rows of a Aeq are checked via QR. But in case
  // they're not then resulting probably will no longer be sparse: it will be
  // slow.
  //
  // MIN_QUAD_WITH_FIXED Minimize quadratic energy 
  //
  // 0.5*Z'*A*Z + Z'*B + C with
  //
  // constraints that Z(known) = Y, optionally also subject to the constraints
  // Aeq*Z = Beq
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
  //   B  n by 1 column of linear coefficients
  //   Y  b by 1 list of constant fixed values
  //   Beq  m by 1 list of linear equality constraint constant values
  // Outputs:
  //   Z  n by cols solution
  //   sol  #unknowns+#lagrange by cols solution to linear system
  // Returns true on success, false on error)igl_Qu8mg5v7";
const char *__doc_igl_min_quad_with_fixed = R"igl_Qu8mg5v7(See min_quad_with_fixed for the documentation.)igl_Qu8mg5v7";
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
const char *__doc_igl_cross_field_missmatch = R"igl_Qu8mg5v7(// Inputs:
  //   V         #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F         #F by 3 eigen Matrix of face (quad) indices
  //   PD1       #F by 3 eigen Matrix of the first per face cross field vector
  //   PD2       #F by 3 eigen Matrix of the second per face cross field vector
  //   isCombed  boolean, specifying whether the field is combed (i.e. matching has been precomputed.
  //             If not, the field is combed first.
  // Output:
  //   Handle_MMatch    #F by 3 eigen Matrix containing the integer missmatch of the cross field
  //                    across all face edges
  //)igl_Qu8mg5v7";
const char *__doc_igl_grad = R"igl_Qu8mg5v7(// Gradient of a scalar function defined on piecewise linear elements (mesh)
  // is constant on each triangle i,j,k:
  // grad(Xijk) = (Xj-Xi) * (Vi - Vk)^R90 / 2A + (Xk-Xi) * (Vj - Vi)^R90 / 2A
  // where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
  // i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of
  // 90 degrees
  //)igl_Qu8mg5v7";
const char *__doc_igl_slice_into = R"igl_Qu8mg5v7(// Act like the matlab Y(row_indices,col_indices) = X
  // 
  // Inputs:
  //   X  xm by xn rhs matrix
  //   R  list of row indices
  //   C  list of column indices
  //   Y  ym by yn lhs matrix
  // Output:
  //   Y  ym by yn lhs matrix, same as input but Y(R,C) = X)igl_Qu8mg5v7";
const char *__doc_igl_n_polyvector = R"igl_Qu8mg5v7(// Inputs:
  //   v0, v1         the two #3 by 1 vectors
  //   normalized     boolean, if false, then the vectors are normalized prior to the calculation
  // Output:
  //                  3 by 3 rotation matrix that takes v0 to v1
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
const char *__doc_igl_boundary_loop = R"igl_Qu8mg5v7(// Compute list of ordered boundary loops for a manifold mesh.
  //
  // Templates:
  //  Index  index type
  // Inputs:
  //   F  #V by dim list of mesh faces
  // Outputs:
  //   L  list of loops where L[i] = ordered list of boundary vertices in loop i
  //)igl_Qu8mg5v7";
const char *__doc_igl_comb_cross_field = R"igl_Qu8mg5v7(// Inputs:
  //   V          #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F          #F by 4 eigen Matrix of face (quad) indices
  //   PD1in      #F by 3 eigen Matrix of the first per face cross field vector
  //   PD2in      #F by 3 eigen Matrix of the second per face cross field vector
  // Output:
  //   PD1out      #F by 3 eigen Matrix of the first combed cross field vector
  //   PD2out      #F by 3 eigen Matrix of the second combed cross field vector
  //)igl_Qu8mg5v7";
const char *__doc_igl_invert_diag = R"igl_Qu8mg5v7(// Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   X  an m by n sparse matrix
  // Outputs:
  //   Y  an m by n sparse matrix)igl_Qu8mg5v7";
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
    // TODO: rename the parameters name in the cpp consistenly
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
    //   soft    the strenght of the soft contraints w.r.t. smoothness
    //           (0 -> smoothness only, 1->constraints only)
    // Outputs:
    //   R       #F by 3 the representative vectors of the interpolated field
    //   S       #V by 1 the singularity index for each vertex (0 = regular))igl_Qu8mg5v7";
