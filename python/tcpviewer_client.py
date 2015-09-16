# Echo client program
import socket
import igl
import array

HOST = ''    # The remote host
PORT = 50008              # The same port as used by the server
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((HOST, PORT))

V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF('../tutorial/shared/beetle.off', V, F)
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V,F)

data = viewer.serialize()
print(type(data[0]).__name__)
viewer.deserialize(data)
a = array.array('u',data)
s.sendall(a)
s.close()
