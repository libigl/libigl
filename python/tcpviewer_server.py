# Echo server program
import socket
import multiprocessing
import igl

HOST = ''                 # Symbolic name meaning all available interfaces
PORT = 50008              # Arbitrary non-privileged port
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((HOST, PORT))
s.listen(1)

try:
    while True:
        conn, addr = s.accept()
        print 'Connected by', addr
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
    print "Closing socket:"
    s.close()
