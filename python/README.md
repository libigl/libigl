# Python wrappers for libigl


## Work in progress
<span style="color:#F62217">
Everything in this folder is currently being developed and it is likely to be
changed radically in the next couple of months, breaking compatibility between
different version. We plan to stabilize the python API by the end of 2016.
</span>

## Introduction

libigl functions can be called natively from python by compiling the wrappers
in this folder. The wrappers supports both python 2.7 and python 3.5 and are
generated using [pybind11](https://github.com/wjakob/pybind11).

The generated library will statically link against all dependencies producing a single,
self-contained binary.

## Installation

The python bindings can be compiled with the following instructions, assuming
that your terminal is pointing to the root of libigl:

```bash
cd python
mkdir build
cd build; cmake ..; make; cd ..
```

The cmake script will complain if it is not able to find python. In that case
you can specify the location of the interpreter by specifying the following
cmake variables.

MacOSX/Linux:

```cmake
SET(PYTHON_LIBRARIES "/usr/local/Cellar/python3/3.5.0/Frameworks/Python.framework/Versions/3.5/lib/libpython3.5m.dylib")
SET(PYTHON_INCLUDE_DIR "/usr/local/Cellar/python3/3.5.0/Frameworks/Python.framework/Versions/3.5/include/python3.5m")
```

Windows:

```cmake
SET(PYTHON_LIBRARIES "C:/Python35/libs/python35.lib")
SET(PYTHON_INCLUDE_DIR "C:/Python35/include")
```

## Tutorial

All libigl tutorials will be ported to python and will use the same naming
scheme. You can find the tutorials in the folder python/tutorials and you can
launch them with the following commands:

```bash
cd python
python 102_DrawMesh.py
```

## Matrix Representation

TODO: describe in detail the wrapped eigen classes and how to convert them to
numpy.

## Viewer and callbacks

The igl viewer provides a convenient and efficient way of visualizing 3D
surfaces in python. It behaves in the same way as the C++ viewer and supports
native python functions as callbacks. This is a simple example that loads
two meshes and switches between the two when a key is pressed:

```python
import pyigl as igl

V1 = igl.eigen.MatrixXd()
F1 = igl.eigen.MatrixXi()

V2 = igl.eigen.MatrixXd()
F2 = igl.eigen.MatrixXi()

def key_pressed(viewer, key, modifier):
    print("Key: ", chr(key))

    if key == ord('1'):
        # # Clear should be called before drawing the mesh
        viewer.data.clear();
        # # Draw_mesh creates or updates the vertices and faces of the displayed mesh.
        # # If a mesh is already displayed, draw_mesh returns an error if the given V and
        # # F have size different than the current ones
        viewer.data.set_mesh(V1, F1);
        viewer.core.align_camera_center(V1,F1);
    elif key == ord('2'):
        viewer.data.clear();
        viewer.data.set_mesh(V2, F2);
        viewer.core.align_camera_center(V2,F2);
    return False


#  Load two meshes
igl.readOFF("../tutorial/shared/bumpy.off", V1, F1);
igl.readOFF("../tutorial/shared/fertility.off", V2, F2);

print("1 Switch to bump mesh")
print("2 Switch to fertility mesh")

viewer = igl.viewer.Viewer()

# Register a keyboard callback that allows to switch between
# the two loaded meshes
viewer.callback_key_pressed = key_pressed
viewer.data.set_mesh(V1, F1)
viewer.launch()
```

### Remote viewer

Whe using the viewer from an interactive python shell (iPython), it is
inconvenient to let the viewer take control of the main thread for rendering
purposes. We provide a simple wrapper for the viewer that allows to launch
a remote process and send meshes to it via a TCP/IP socket. For more
informations on how to use it see the documentation in [tcpviewer.py](tcpviewer.py)

## Matlab

The python wrappers can be natively being used from MATLAB.
We provide a few examples in the folder python/matlab.

## Documentation

The python functions have exactly the same prototypes as their C++ counterpart.
Docstrings for all available python functions are extracted from the C++ header files and compiled into the python module. To get help for a certain function, you can run `help(pyigl.<function_name>)` in the python console.

In the scripts folder there is the script `generate_docstrings.py` that automatically generates python docstrings for a new function. You can run it with `generate_docstrings.py <path_to_cpp_header_files> <path_to_python_files>`. 
The script depends on additional libraries (joblib, mako, clang), make sure to install them (e.g. through pip. python-clang is included in external/nanogui/ext/pybind11/tools/clang).


## Known Issues

## Contact

Libigl is a group endeavor led by [Alec
Jacobson](http://www.cs.columbia.edu/~jacobson/) and [Daniele
Panozzo](http://cs.nyu.edu/~panozzo/). Please [contact
us](mailto:alecjacobson@gmail.com,daniele.panozzo@gmail.com) if you have
questions or comments. For troubleshooting, please post an
[issue](https://github.com/libigl/libigl/issues) on github.

If you're using libigl in your projects, quickly [drop us a
note](mailto:alecjacobson@gmail.com,daniele.panozzo@gmail.com). Tell us who you
are and what you're using it for. This helps us apply for funding and justify
spending time maintaining this.

If you find bugs or have problems please use our [github issue tracking
page](https://github.com/libigl/libigl/issues).

## Copyright
2015 Alec Jacobson, Daniele Panozzo, Christian Schüller, Olga Diamanti, Qingnan
Zhou, Sebastian Koch, Nico Pietroni, Stefan Brugger, Kenshi Takayama, Wenzel Jakob, Nikolas De
Giorgis, Luigi Rocca, Leonardo Sacht, Olga Sorkine-Hornung, and others.

Please see individual files for appropriate copyright notices.
