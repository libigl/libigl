import igl

# Load mesh
V = igl.eigen.MatrixXd();
F = igl.eigen.MatrixXi();
igl.readOFF("../tutorial/shared/bumpy.off",V,F);

# Compute Gaussian curvature
K = igl.eigen.VectorXd();
igl.gaussian_curvature(V,F,K);

print("igl::gaussian_curvature: \n", K, sep='')

# Compute pseudocolor
C = igl.eigen.MatrixXd();
igl.jet(K,True,C);
print("igl::jet: \n", C, sep='')
