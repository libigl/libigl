# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
import socket
import threading
import pyigl as igl
import array
import time

HOST = 'localhost'                 # Symbolic name meaning all available interfaces
PORT = 50008              # Arbitrary non-privileged port

def worker(viewer,lock,s):

    print("TCP iglviewer server listening on port " + str(PORT))
    try:
        while True:
            conn, addr = s.accept()
            lock.acquire()
            slist = []
            while True:
                buf = conn.recv(4096)
                if not buf:
                    break
                slist.append(buf.decode('unicode_internal','ignore'))
            conn.close()

            data = ''.join(slist)
            temp = list(data)

            isempty = viewer.data().V.rows() == 0
            viewer.data().deserialize(temp)
            if isempty and viewer.data().V.rows() != 0:
                viewer.core.align_camera_center(viewer.data().V,viewer.data().F)

            lock.release()

    except:
        s.close()
    return

class TCPViewer(igl.glfw.Viewer):
    def launch(self):
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.connect((HOST, PORT))
            ser = self.data().serialize()
            a = array.array('u', ser)
            s.sendall(a)
            s.close()
        except:
            print("Failed to open socket, is tcpviewer running?")

if __name__ == "__main__": # The main script is a server

    ## Try to open the socket first
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.bind((HOST, PORT))
    except:
        print("Failed to bind, port already used.")
        exit(1)
    s.listen(1)

    viewer = igl.glfw.Viewer()

    lock = threading.Lock()
    t = threading.Thread(target=worker, args=(viewer,lock,s,))
    t.setDaemon(True)
    t.start()

    viewer.core.is_animating = True
    # viewer.data().dirty = int(0x03FF)

    viewer.launch_init(True,False)
    done = False
    while not done:
        lock.acquire()
        done = not viewer.launch_rendering(False)
        lock.release()

        time.sleep(0.000001) # DO NOT REMOVE ME

    viewer.launch_shut()
