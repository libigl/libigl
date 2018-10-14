# FAQ

> **Q:** I'd like to merge two 3D meshes into one. How to do it with igl?

It sounds like you're trying to compute a union. Are your two meshes closed, watertight manifolds? Then you could call `igl/boolean/mesh_boolean` with the union option. If not, then there's still hope with something else. _[Alec]_

> **Q:** In other apps I have seen the the user is asked to specify singularities , and the then the rosy field is generated. But in libigl it seems like you have to specify faces and direction vectors to design a field. Is it possible to specify singularities?

It is not possible to specify singularities right now. To specify the directions, the vectors should be in global coordinates (the vectors are 3D vectors, the libigl function takes care of projecting them onto the corresponding face), you can take a look here for a basic example that fixes only one face: [Example 505](http://libigl.github.io/libigl/tutorial/#global-seamless-integer-grid-parametrization) _[Daniele]_

> **Q:** Does Libigl use the same 2D Triangle code (my search in the Libigl source code indicates NO, but a confirmation would be reassuring)?

No, it uses CGAL for triangulation. _[Alec]_

> **Q:** Libigl's Boolean depends on some GPL-licensed header files from CGAL. Is it possible to remove this dependency?

No, the dependency on CGAL would require severely rewriting a core function. It is possible to do, but I will not do it. _[Alec]_

> **Q:** Do you have a ready to run command line program so that I can run a test with a few of my sample data sets?

No, but it would be very easy to alter the [boolean tutorial example](http://libigl.github.io/libigl/tutorial/#boolean-operations-on-meshes) to do that. Basically drop the viewer and change the hardcoded paths to command line arguments and write out the result to an obj. _[Alec]_

> **Q:** I see that it can generate N-rosy fields, but is it possible to remesh based on the rosy field?

Here is an example that uses libigl to produce a seamless parametrization:
[Example 505](http://libigl.github.io/libigl/tutorial/#global-seamless-integer-grid-parametrization)

If you want a mesh, you can pass this parametrization to libQEX ([https://github.com/hcebke/libQEx](https://github.com/hcebke/libQEx)) to extract it. We do not have it built in in the tutorial due to a more restrictive licence used by libQEx. _[Daniele]_

> **Q:** I am having issues with parameterization (igl::miq). Even at 100 iterations, there are still distortions. What is the cause?

This is unfortunately the expected behaviour, the MIQ parametrization tends to concentrate the distortion around singularities. _[Daniele]_

> **Q:** I am receiving compilation errors along the lines of "ISO C++ forbids declaration of ... with no type" when compiling under Windows using gcc.

We never tried to compile libigl on Windows with gcc, but we did test our library on:

* windows using visual studio
* linux with gcc
* macosx with clang

You might have a version of gcc that doesn't support (enough of) c++11. Try using Cygwin and g++ 4.9.2. _[Alec, Daniele]_
