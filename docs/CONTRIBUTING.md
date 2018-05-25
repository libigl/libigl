Before opening an issue on creating a pull request, please check the following:

## Compilation Issues

- If you are on Windows, did you select the **x64** version of the Visual Studio compiler?

- If you have a **CMake issue**, make sure you follow the same approach as the  [libigl-example-project](https://github.com/libigl/libigl-example-project) to build libigl with your project, and make sure that you can compile the example project.

- If you have an issue with a **submodule**, check if your submodules are up to date. If you have a doubt about a submodule, delete its folder and run `git submodule update --init --recursive` in the libigl directory.

- If you have an issue with a missing **template issue**, check if your code compile with the *header-only* option of libigl activated. Turn **`OFF`** the CMake option `LIBIGL_USE_STATIC_LIBRARY`: either modify your `CMakeCache.txt` via CMake GUI or ccmake, or delete your `CMakeCache.txt` and re-run `cmake -DLIBIGL_USE_STATIC_LIBRARY=OFF ..` in your build folder.

- Make sure your read the [**FAQ**](https://github.com/libigl/libigl/wiki/FAQ) before asking a new question, and search [**existing issues**](https://github.com/libigl/libigl/issues?q=is%3Aissue+is%3Aclosed) for a problem similar to yours.

- Make sure you read the informations contained in the libigl [homepage](https://github.com/libigl/libigl) as well as the [tutorials](http://libigl.github.io/libigl/tutorial/tutorial.html).

- If none of these solve your problem, then please report your issue in the bug tracker!
