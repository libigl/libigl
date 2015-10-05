# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import igl
import random
from math import cos,sin,pi

# Input mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

# Per face bases
B1 = igl.eigen.MatrixXd()
B2 = igl.eigen.MatrixXd()
B3 = igl.eigen.MatrixXd()

# Face barycenters
B = igl.eigen.MatrixXd()

# Scale for visualizing the fields
global_scale = 1

# Random length factor
rand_factor = 5

samples = igl.eigen.MatrixXi()

def readSamples(fname):
    samples = igl.eigen.MatrixXi()
    numSamples = 0

    fp = open(fname, 'r')

    numSamples = int(fp.readline())

    samples.resize(numSamples,1)

    for i in range(0,numSamples):
        samples[i] = int(fp.readline())

    fp.close()

    return samples

# Create a random set of tangent vectors
def random_constraints(b1, b2, n):

    r = igl.eigen.MatrixXd(1,n*3)

    for i in range(0,n):
        a = random.random()*2*pi
        s = 1 + random.random() * rand_factor
        t = s * (cos(a) * b1 + sin(a) * b2)
        r.setBlock(0,i*3,1,3,t)

    return r

def key_down(viewer, key, modifier):
    if key < ord('1') or key > ord('8'):
        return False

    viewer.data.lines.resize(0,9)

    num = key  - ord('0')

    # Interpolate
    print("Interpolating " + repr(num * 2) + "-PolyVector field")

    b = igl.eigen.MatrixXi([[4550, 2321, 5413, 5350]]).transpose()

    bc = igl.eigen.MatrixXd(b.size(),num*3)

    for i in range(0,b.size()):
        t = random_constraints(B1.row(b[i]),B2.row(b[i]),num)
        bc.setRow(i,t)

    # Interpolated PolyVector field
    pvf = igl.eigen.MatrixXd()
    igl.n_polyvector(V, F, b, bc, pvf)

    # Highlight in red the constrained faces
    C = igl.eigen.MatrixXd.Constant(F.rows(),3,1)

    for i in range(0,b.size()):
        C.setRow(b[i],igl.eigen.MatrixXd([[1, 0, 0]]))
    viewer.data.set_colors(C)

    for n in range(0,num):
        VF = igl.eigen.MatrixXd.Zero(F.rows(),3)

        for i in range(0,b.size()):
            VF.setRow(b[i],bc.block(i,n*3,1,3))

        for i in range(0,samples.rows()):
            VF.setRow(samples[i],pvf.block(samples[i],n*3,1,3))

        c = VF.rowwiseNorm()

        C2 = igl.eigen.MatrixXd()
        igl.jet(c,1,1+rand_factor,C2)
        viewer.data.add_edges(B - global_scale*VF, B + global_scale*VF , C2)

    return False


# Load a mesh in OBJ format
igl.readOBJ("../../tutorial/shared/lilium.obj", V, F)
samples = readSamples("../../tutorial/shared/lilium.samples.0.2")

# Compute local basis for faces
igl.local_basis(V,F,B1,B2,B3)

# Compute face barycenters
igl.barycenter(V, F, B)

# Compute scale for visualizing fields
global_scale = 0.2*igl.avg_edge_length(V, F)

# Make the example deterministic
random.seed(0)

viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)
viewer.callback_key_down = key_down
viewer.core.show_lines = False

key_down(viewer,ord('2'),0)

viewer.launch()
