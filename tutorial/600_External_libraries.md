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



* discretization are useful to solve PDE
* to create a 2D triangulation we provide a convenient wrapper for triangle

#Tetrahedralization of closed surfaces

* the same in 3D, similar interface

#Baking ambient occlusion

* intro to ambinet occlusion
* can be easily computed with raytracing
* again a simple wrapper
* then you multiply to the colors


#Outlook for continuing development

* better documentation
* your contributions are welcome, using pull request
* open things to do
  * isotropic remeshing
  * matlab wrappers
  * mixed integer solvers
  * fast spatial indices
