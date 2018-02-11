glad
====

The files in this folder were generated using [glad](https://github.com/Dav1dde/glad).

### Using a different version of OpenGL

The files in this folder correspond to the OpenGL 3.3 Core Profile. If you need to use libigl with a different version of OpenGL, you can use this [webservice](http://glad.dav1d.de/) to generate the loader for your desired version of OpenGL. Then, in your CMake build system, simply define the target `glad` to point to your own version of glad before including libigl.
