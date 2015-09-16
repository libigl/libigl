import igl
import tcpviewer
import time

# Load a mesh in OFF format
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
time1 = time.time()
igl.read_triangle_mesh("../tutorial/shared/armadillo.obj", V, F)
time2 = time.time()


print 'Loading mesh (%d vertices) %0.3f ms' % (V.rows(),(time2-time1)*1000.0)


# Plot the mesh
viewer = tcpviewer.TCPViewer();
viewer.data.set_mesh(V, F);
viewer.launch();
time3 = time.time()

print 'Sending to TCP viewer took %0.3f ms' % ((time3-time2)*1000.0)
