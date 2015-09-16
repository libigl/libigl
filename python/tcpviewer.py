import socket
import multiprocessing
import igl
import array

HOST = ''                 # Symbolic name meaning all available interfaces
PORT = 50008              # Arbitrary non-privileged port

class TCPViewer(igl.viewer.Viewer):
    def launch(self):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((HOST, PORT))
        a = array.array('u',self.serialize())
        s.sendall(a)
        s.close()

if __name__ == "__main__": # The main script is a server
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind((HOST, PORT))
    s.listen(1)
    print("TCP iglviewer server listening on port " + str(PORT))
    try:
        while True:
            conn, addr = s.accept()
            data = ''
            while True:
                datanew = conn.recv(1024)
                if not datanew:
                    break
                data = data+datanew
            conn.close()

            def worker(data):
                viewer = igl.viewer.Viewer()
                temp = list(data.decode('unicode_internal','ignore'))
                viewer.deserialize(temp)
                viewer.launch(True,False)
                return

            t = multiprocessing.Process(target=worker, args=(data,))
            t.start()
    except:
        s.close()
