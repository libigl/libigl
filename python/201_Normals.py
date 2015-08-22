import igl

# Load a mesh in OFF format
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF("../tutorial/shared/fandisk.off", V, F);

# Compute per-face normals
N_faces = igl.eigen.MatrixXd()
igl.per_face_normals(V,F,N_faces);
print("igl::per_face_normals: \n", N_faces, sep='')

# Compute per-vertex normals
N_vertices = igl.eigen.MatrixXd()
igl.per_vertex_normals(V,F,igl.PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA,N_vertices);
print("igl::per_vertex_normals: \n", N_vertices, sep='')

# Compute per-corner normals, |dihedral angle| > 20 degrees --> crease
N_corners = igl.eigen.MatrixXd()
igl.per_corner_normals(V,F,20,N_corners);
print("igl::per_corner_normals: \n", N_corners, sep='')
