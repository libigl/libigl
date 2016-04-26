#include "python_shared.h"
#include <sstream>
#include <string>
#include <fstream>

extern void python_export_vector(py::module &);
extern void python_export_igl(py::module &);

#ifdef PY_VIEWER
extern void python_export_igl_viewer(py::module &);
#endif

#ifdef PY_COMISO
extern void python_export_igl_comiso(py::module &);
#endif

PYTHON_PLUGIN(pyigl) {
    py::init_threading();
    py::module m("pyigl", R"pyigldoc(
        Python wrappers for libigl
        --------------------------

        .. currentmodule:: pyigl

        .. autosummary::
           :toctree: _generate

           principal_curvature
           local_basis
           cotmatrix
           floor
           slice
           per_face_normals
           ARAPEnergyType
           readOFF
           per_vertex_normals
           sortrows
           barycenter
           jet
           eigs
           per_corner_normals
           massmatrix
           colon
           rotate_vectors
           read_triangle_mesh
           SolverStatus
           gaussian_curvature
           avg_edge_length
           lscm
           find_cross_field_singularities
           parula
           setdiff
           comb_frame_field
           map_vertices_to_circle
           writeOBJ
           active_set
           boundary_facets
           compute_frame_field_bisectors
           readOBJ
           cut_mesh_from_singularities
           readDMAT
           doublearea
           min_quad_with_fixed
           unique
           arap
           cross_field_missmatch
           grad
           slice_into
           n_polyvector
           harmonic
           boundary_loop
           comb_cross_field
           invert_diag
           copyleft_comiso_miq
           copyleft_comiso_nrosy

    )pyigldoc");

    python_export_vector(m);
    python_export_igl(m);


    #ifdef PY_VIEWER
    python_export_igl_viewer(m);
    #endif

    #ifdef PY_COMISO
    python_export_igl_comiso(m);
    #endif

    return m.ptr();
}
