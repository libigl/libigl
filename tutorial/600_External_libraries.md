title: libigl Tutorial
author: Daniele Panozzo, Alec Jacobson and others
date: 20 June 2014
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

* [Chapter 6: External libraries][600]
    * [601 State serialization][601]
    * [602 Mixing matlab code][602]
    * [603 Calling igl functions from matlab][603]
    * [604 Triangulation of closed polygons][604]
    * [605 Tetrahedralization of closed surfaces][605]
    * [606 Baking ambient occlusion][606]

# Chapter 6: External libraries [600]

An additional positive side effect of using matrices as basic types is that it is easy to exchange data between libigl and other softwares and libraries.

## State serialization [601]

Geometry processing applications often require a considerable amount of computational time and/or manual input. In order to make the development efficient it must be possible to serialize and deserialize the state of the application.

Having a good serialization framework allows to quickly start debugging just before the crash happens, avoiding to wait for the precomputation to take place every time. It also makes it easier to define unit testing that can be used to find bugs in interactive applications: if the input is slightly different every time the algorithm is executed, it is very difficult to find bugs.

Unfortunately, serialization is often not considered in geoemtry processing due to the extreme difficulty in serializing pointer-based data structures (like an helf-edge).

In libigl, serialization is simpler, since the majority of the functions use basic types, and pointers are used in very rare cases (usually to interface with external libraries). libigl provides an extremely easy to use XML serialization framework, that drastically reduces the overhead required to add serialization to your applications.

Assume that the state of your application is composed of a mesh and set of integer ids:
``` cpp
class State : public ::igl::XMLSerialization
{
public:
  State() : XMLSerialization("dummy") {}

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  std::vector<int> ids;

  void InitSerialization()
  {
    xmlSerializer->Add(V  , "V");
    xmlSerializer->Add(F  , "F");
    xmlSerializer->Add(ids, "ids");
  }
};
```

A class can be made serializable by inheriting from ::igl::XMLSerialization and trivially implementing the InitSerialization method. Note that you don't have to care the types, Add is able to serialize all basic stl types, all Eigen types and any class inheriting from ::igl::XMLSerialization.

It is then possible to save the state to an xml file:

``` cpp
::igl::XMLSerializer serializer_save("601_Serialization");
serializer_save.Add(state,"State");
serializer_save.Save("temp.xml",true);
```

This code generates the following xml file (assuming V and F contains a simple mesh with two triangles, and ids contains the numbers 6 and 7):
``` xml
<:::601_Serialization>
    <State>
        <V rows="4" cols="3" matrix="
0,0,0,
1,0,0,
1,1,1,
2,1,0"/>
        <F rows="2" cols="3" matrix="
0,1,2,
1,3,2"/>
        <ids size="2" vector_int="
6,7"/>
    </State>
</:::601_Serialization>
```

The xml file can then be loaded in a similar way:

``` cpp
State loaded_state;
::igl::XMLSerializer serializer_load("601_Serialization");
serializer_load.Add(loaded_state,"State");
serializer_load.Load("temp.xml");
```

This can also be used as a convenient interface to provide parameters to command line applications, since the xml files can be directly edited with a standard text editor.

We demonstrate the serialization framework in [Example 601](601_Serialization/main.cpp). We strongly suggest that you make the entire state of your application always serializable: this will save you a lot of troubles when you'll be making figures for a scientific publication. It is very common to have to do small changes to figures during the production of a paper, and being able to serialize the entire state just before you take screenshots will save you many painful hours before a submission deadline.

## Mixing matlab code [602]

libigl can be interfaced matlab, to offload some of the numerically heavy computation to a matlab script. This has the major advantage of allowing to develop efficient and complex UI in C++, while keeping the advantage of fast protototyping of matlab. In particular, using an external matlab script in a libigl application allows to change the algorithm in the matlab script without having to recompile the C++ part.

We demonstrate how to integrate matlab in a libigl application in [Example 602](602_Matlab/main.cpp). The example uses matlab to compute the Eigenfunctions of the discrete Laplacian operator, relying on libigl for mesh IO, visualization and for computing the Laplacian operator.

libigl can connect to an existing instance of matlab (or launching a new one on Linux/MacOSX) using:

``` cpp
igl::mlinit(&engine);
```

The cotangent laplacian is computed using igl::cotmatrix and uploaded to the matlab workspace:

``` cpp
igl::cotmatrix(V,F,L);
igl::mlsetmatrix(&engine,"L",L);
```

It is now possible to use any matlab function on the data. For example, we can see the sparsity pattern of L using spy:

``` cpp
igl::mleval(&engine,"spy(L)");
```

![The matlab spy function is called from a libigl-based application.](images/602_Matlab_1.png)

You can also do some computation and then return it back to the C++ application

``` cpp
igl::mleval(&engine,"[EV,~] = eigs(-L,10,'sm')");
igl::mlgetmatrix(&engine,"EV",EV);
```

and then use libigl functions to plot the eigenfunctions.

![4 Eigenfunctions of the Laplacian plotted in the libigl viewer.](images/602_Matlab_2.png)

## Calling igl functions from matlab [603]

It is also possible to call libigl functions from matlab, compiling them as MEX functions. This can be very useful to offload to C++ code the computationally intensive parts of a matlab application.

We provide a wrapper for igl::readOBJ in [Example 603](603_MEX/compileMEX.m). We plan to provide wrappers for all our functions in the future, if you are interested in this feature (or if you want to help implementing it) please let us know.

## Triangulation of closed polygons [604]

The generation of high-quality triangle and tetrahedral meshes is a very common task in geometry processing. We provide wrappers in libigl to triangle and tetegen.

A triangle mesh canb e cerated starting from a set of boundary edges using igl::triangulate.

``` cpp
igl::triangulate(V,E,H,V2,F2,"a0.005q");
```

where E is a set of boundary edges, H a set of 2D positions of points contained in holes of the triangulation and (V2,F2) is the generate triangulation. Additional parameters can be given to triangles, to control the quality: "a0.005q" puts a bound on the maximal area of the triangles and a minimal angle of 20 degrees. In Example [Example 604](604_Triangle/main.m), the interior of a square (excluded a smaller square in its interior) is triangulated.

![Triangulation of the interior of a polygon.](images/604_Triangle.png)

## Tetrahedralization of closed surfaces [605]

Similarly, the interior of a closed manifold surface can be tetrahedralized using the function igl::tetrahedralize which wraps the tetgen library ([Example 605](605_Tetgen/main.c)):

``` cpp
igl::tetrahedralize(V,F,"pq1.414", TV,TT,TF);
```

![Tetrahedralization of the interior of a surface mesh.](images/605_Tetgen.png)

## Baking ambient occlusion [606]

[Ambient occlusion](http://en.wikipedia.org/wiki/Ambient_occlusion) is a rendering technique used to calculate the exposure of each point in a surface to ambient lighting. It is usually encoded as a scalar (normalized between 0 and 1) associated with the vertice of a mesh.

Formally, ambient occlusion is defined as:

\\[ A_p = \frac{1}{\pi} \int_\omega V_{p,\omega}(n \cdot \omega) d\omega \\]

where \\( V_{p,\omega} \\) is the visibility function at  p, defined to be zero if p is occluded in the direction \\( \omega \\) and one otherwise, and \\( d\omega \\) is the infinitesimal solid angle step of the integration variable \\( \omega \\).

The integral is usually approximate by casting rays in random directions around each vertex. This approximation can be computed using the function:

``` cpp
igl::ambient_occlusion(V,F,V_samples,N_samples,500,AO);
```

that given a scene described in V,F, computes the ambient occlusion of the points in V_samples whose associated normals are N_samples. The number of casted rays can be controlled (usually at least 400-500 rays are required to get a smooth result) and the result is return in AO, as a single scalar for each sample.

Ambient occlusion can be used to darken the surface colors, as shown in [Example 606](606_AmbientOcclusion/main.c)

![A mesh rendered without (left) and with (right) ambient occlusion.](images/606_AmbientOcclusion.png)

# Outlook for continuing development

* better documentation
* your contributions are welcome, using pull request
* open things to do
  * isotropic remeshing
  * matlab wrappers
  * mixed integer solvers
  * fast spatial indices
