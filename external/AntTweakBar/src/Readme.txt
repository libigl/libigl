Under Windows, it is not necessary to rebuild AntTweakBar since a precompiled
version is provided in the lib directory. But if you want to recompile it,
you can use the provided Visual Studio solution. You'd also need the DirectX 
SDK (http://msdn.microsoft.com/directx) and the path to the DirectX shader 
compiler fxc.exe (included in the DirectX SDK) must be listed in the VC++ 
directories.

To build the library on Linux, open a terminal, go in the src directory and
type make

To build the library on MacOSX, open a terminal, go in the src directory and
type make -f Makefile.osx
