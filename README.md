# libigl - A simple C++ geometry processing library

[![](https://github.com/libigl/libigl/workflows/Build/badge.svg?event=push)](https://github.com/libigl/libigl/actions?query=workflow%3ABuild+branch%3Amain+event%3Apush)
[![](https://anaconda.org/conda-forge/igl/badges/installer/conda.svg)](https://anaconda.org/conda-forge/igl)

![](https://libigl.github.io/libigl-teaser.png)

Documentation, tutorial, and instructions at <https://libigl.github.io>.

## Upstream CMake Difficulties

The upstream libigl repository acquires all of its dependencies as sub-projects through CMake's FetchContent module.
This means that as part of configuring the CMake project, all of the source code for all of the dependencies are
downloaded, built from source, and added as subdirectories of libigl to be built along with its build process. This 
not only means that for each build directory or each project, all dependencies must be rebuilt, but it also means 
that all dependency build messages are now your project's build messages, flooding the console.
Furthermore, there's no option to turn this feature off. To avoid downloading, consuming projects must first call *
find_package* for all of libigl's
dependencies before finding/adding libigl. Although this pattern is simple across a single link, if the consuming
project is also a library, and one of its consuming projects happen to also use libigl, violations of ODR are very
likely.

Beyond requiring downloading and rebuilding for each use, this pattern is additionally problematic because it doesn't 
follow's CMake's standard *find_package* flow and makes it very difficult for consumers to use libigl through 
*find_package*. The latter of which is primarily driven through libigl's improper installation, despite offering the 
option *LIBIGL_INSTALL*.

Issues with existing installation include:

- installing only the igl::core target
- installing package config files in an improperly named directory. These files need to be installed into a 
directory that starts with the project name, yet, they are installed into ...igl/ instead of ...libigl*/
- manually installing artifacts of libigl's built dependencies, but doing it improperly. Consider Eigen3 - a 
required dependency of every target. This is built, and when installing, is installed into another imroperly named 
directory, likely to avoid conflicting with an existing Eigen3 installation.
- CMake package configuration files (libigl-config.cmake) don't actually include any of the targets that may be have been 
installed & exported. So, even when a consuming project finds libigl-config.cmake, it doesn't do anything.
- header files with *.hpp* extention are disregarded

All of this makes it difficult to package libigl for package managers, like vcpkg or conan. Vcpkg applies patches to 
libigl to ignore FetchConent and replace with *find_package*.

## Modifications

All changes included in this fork are to allow libigl to find its dependencies and consuming projects to find libigl 
through CMake's standard patterns. This effectively means

1. stop forcefully downloading dependencies
2. find dependencies via CMake's *find_package* command
3. provide correct package config files
4. install all targets provided by libigl (but not its dependencies)
5. Consider headers with *.hpp* extension, such that they are installed
 
All of the existing libigl target names are maintained for compatibility.

### Specific Changes

1 & 2. The option *LIBIGL_FIND_PACKAGES* is available to selectively find packages instead of download them.
3. The function *igl_install* has been updated to install the appropriate *-config.cmake and *-targets.cmake for the 
given "module".
4. *igl_install* is used across all igl targets, not just igl::core
5. A new CMake module, *igl_glob_sources*, has been added to perform the globbing that was previously being done 
explicitly for each "module", and it includes *.hpp* files.

## State

wip & untested

## Notes

Some people say that all of a project's dependencies should be built with the same compiler invocation with the same 
set of options at the same time. This just isn't practical, overly beneficial, nor does this occur in reality.

- igl::predicates is still being installed as a sub-project.
- igl_restricted::* are not handled at all


| 🚨 Important |
|:---|
| The latest version of libigl (v2.4.0) introduces some **breaking changes** to its CMake build system. Please read our [changelog](https://libigl.github.io/changelog/) page for instructions on how to update your project accordingly. |
