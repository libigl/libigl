%% Compile the wrapper for readOBJ
mex readOBJ_mex.cpp ...
    -I../../include ...
    -I/opt/local/include/eigen3 %% Change this path to point to your copy of Eigen

%% Load an OBJ mesh
[V,F] = readOBJ_mex('../shared/bumpy-cube.obj');

%% Plot the mesh
trimesh(F,V(:,1),V(:,2),V(:,3));
axis vis3d;
