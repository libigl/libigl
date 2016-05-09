const char *__doc_igl_principal_curvature = R"igl_Qu8mg5v7(See principal_curvature for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_local_basis = R"igl_Qu8mg5v7(See local_basis for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_cotmatrix = R"igl_Qu8mg5v7(See cotmatrix for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_floor = R"igl_Qu8mg5v7(See floor for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_slice = R"igl_Qu8mg5v7(See slice for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_per_face_normals = R"igl_Qu8mg5v7(See per_face_normals for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_per_face_normals_stable = R"igl_Qu8mg5v7(See per_face_normals_stable for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_readOFF = R"igl_Qu8mg5v7(See readOFF for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_per_vertex_normals = R"igl_Qu8mg5v7(See per_vertex_normals for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_sortrows = R"igl_Qu8mg5v7(See sortrows for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_barycenter = R"igl_Qu8mg5v7(See barycenter for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_jet = R"igl_Qu8mg5v7(See jet for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_eigs = R"igl_Qu8mg5v7(See eigs for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_per_corner_normals = R"igl_Qu8mg5v7(See per_corner_normals for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_massmatrix = R"igl_Qu8mg5v7(See massmatrix for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_colon = R"igl_Qu8mg5v7(See colon for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_rotate_vectors = R"igl_Qu8mg5v7(See rotate_vectors for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_read_triangle_mesh = R"igl_Qu8mg5v7(See read_triangle_mesh for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_gaussian_curvature = R"igl_Qu8mg5v7(See gaussian_curvature for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_avg_edge_length = R"igl_Qu8mg5v7(See avg_edge_length for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_lscm = R"igl_Qu8mg5v7(See lscm for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_find_cross_field_singularities = R"igl_Qu8mg5v7(See find_cross_field_singularities for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_parula = R"igl_Qu8mg5v7(See parula for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_setdiff = R"igl_Qu8mg5v7(See setdiff for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_comb_frame_field = R"igl_Qu8mg5v7(See comb_frame_field for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_map_vertices_to_circle = R"igl_Qu8mg5v7(See map_vertices_to_circle for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_writeOBJ = R"igl_Qu8mg5v7(See writeOBJ for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_active_set = R"igl_Qu8mg5v7(See active_set for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_boundary_facets = R"igl_Qu8mg5v7(See boundary_facets for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_compute_frame_field_bisectors = R"igl_Qu8mg5v7(See compute_frame_field_bisectors for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_readOBJ = R"igl_Qu8mg5v7(See readOBJ for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_cut_mesh_from_singularities = R"igl_Qu8mg5v7(See cut_mesh_from_singularities for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_readDMAT = R"igl_Qu8mg5v7(See readDMAT for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_doublearea = R"igl_Qu8mg5v7(/**
     DOUBLEAREA computes twice the area for each input triangle[quad]
    
     Templates:
       DerivedV  derived type of eigen matrix for V (e.g. derived from
         MatrixXd)
       DerivedF  derived type of eigen matrix for F (e.g. derived from
         MatrixXi)
       DeriveddblA  derived type of eigen matrix for dblA (e.g. derived from
         MatrixXd)
     Inputs:
       V  #V by dim list of mesh vertex positions
       F  #F by simplex_size list of mesh faces (must be triangles or quads)
     Outputs:
       dblA  #F list of triangle[quad] double areas (SIGNED only for 2D input)
    
     Known bug: For dim==3 complexity is O(#V + #F)!! Not just O(#F). This is a big deal
     if you have 1million unreferenced vertices and 1 face
    **/)igl_Qu8mg5v7";
const char *__doc_igl_doublearea_single = R"igl_Qu8mg5v7(/**
     Single triangle in 2D!
     This should handle streams of corners not just single corners
    **/)igl_Qu8mg5v7";
const char *__doc_igl_doublearea_quad = R"igl_Qu8mg5v7(/**
     DOUBLEAREA_QUAD computes twice the area for each input quadrilateral
    
     Inputs:
       V  #V by dim list of mesh vertex positions
       F  #F by simplex_size list of mesh faces (must be quadrilaterals)
     Outputs:
       dblA  #F list of quadrilateral double areas
    
    **/)igl_Qu8mg5v7";
const char *__doc_igl_min_quad_with_fixed_precompute = R"igl_Qu8mg5v7(See min_quad_with_fixed_precompute for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_min_quad_with_fixed_solve = R"igl_Qu8mg5v7(See min_quad_with_fixed_solve for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_min_quad_with_fixed = R"igl_Qu8mg5v7(See min_quad_with_fixed for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_unique = R"igl_Qu8mg5v7(See unique for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_unique_rows = R"igl_Qu8mg5v7(See unique_rows for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_arap_precomputation = R"igl_Qu8mg5v7(See arap_precomputation for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_arap_solve = R"igl_Qu8mg5v7(See arap_solve for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_cross_field_missmatch = R"igl_Qu8mg5v7(See cross_field_missmatch for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_grad = R"igl_Qu8mg5v7(See grad for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_slice_into = R"igl_Qu8mg5v7(See slice_into for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_n_polyvector = R"igl_Qu8mg5v7(See n_polyvector for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_harmonic = R"igl_Qu8mg5v7(See harmonic for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_boundary_loop = R"igl_Qu8mg5v7(See boundary_loop for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_comb_cross_field = R"igl_Qu8mg5v7(See comb_cross_field for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_invert_diag = R"igl_Qu8mg5v7(See invert_diag for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_copyleft_comiso_miq = R"igl_Qu8mg5v7(See miq for the documentation.)igl_Qu8mg5v7";
const char *__doc_igl_copyleft_comiso_nrosy = R"igl_Qu8mg5v7(See nrosy for the documentation.)igl_Qu8mg5v7";
