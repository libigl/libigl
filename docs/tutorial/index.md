title: libigl Tutorial
author: Daniele Panozzo and Alec Jacobson
date: 07 November 2015
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

# libigl tutorial notes

#### originally presented by Daniele Panozzo and Alec Jacobson at SGP Graduate School 2014

![](images/libigl-logo.jpg)

Libigl is an open source C++ library for geometry processing research and development.  Dropping the heavy data structures of tradition geometry libraries, libigl is a simple header-only library of encapsulated functions. This combines the rapid prototyping familiar to Matlab or Python programmers with the performance and versatility of C++.  The tutorial is a self-contained, hands-on introduction to libigl.  Via interactive, step-by-step examples, we demonstrate how to accomplish common geometry processing tasks such as computation of differential quantities and operators, real-time deformation, parametrization, numerical optimization and remeshing. Each section of the lecture notes links to a cross-platform example application.


[#attene_2014]: Marco Attene. [Direct repair of self-intersecting meshes](https://www.google.com/search?q=Direct+repair+of+self-intersecting+meshes), 2014.
[#baerentzen_2005]: J Andreas Baerentzen and Henrik Aanaes. [Signed distance computation using the angle weighted pseudonormal](https://www.google.com/search?q=Signed+distance+computation+using+the+angle+weighted+pseudonormal), 2005.
[#barbic_2005]: Jernej Barbic and Doug James. [Real-Time Subspace Integration for St.Venant-Kirchhoff Deformable Models](https://www.google.com/search?q=Real-Time+Subspace+Integration+for+St.Venant-Kirchhoff+Deformable+Models), 2005.
[#bommes_2009]: David Bommes, Henrik Zimmer, Leif Kobbelt. [Mixed-integer quadrangulation](http://www-sop.inria.fr/members/David.Bommes/publications/miq.pdf), 2009.
[#botsch_2004]: Matrio Botsch and Leif Kobbelt. [An Intuitive Framework for Real-Time Freeform Modeling](https://www.google.com/search?q=An+Intuitive+Framework+for+Real-Time+Freeform+Modeling), 2004.
[#bouaziz_2012]: Sofien Bouaziz, Mario Deuss, Yuliy Schwartzburg, Thibaut Weise, Mark Pauly [Shape-Up: Shaping Discrete Geometry with Projections](http://lgg.epfl.ch/publications/2012/shapeup.pdf), 2012
[#chao_2010]: Isaac Chao, Ulrich Pinkall, Patrick Sanan, Peter Schröder. [A Simple Geometric Model for Elastic Deformations](https://www.google.com/search?q=A+Simple+Geometric+Model+for+Elastic+Deformations), 2010.
[#diamanti_2014]: Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga Sorkine-Hornung. [Designing N-PolyVector Fields with Complex Polynomials](http://igl.ethz.ch/projects/complex-roots/), 2014
[#diamanti_2015]: Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga Sorkine-Hornung. [Integrable PolyVector Fields](http://igl.ethz.ch/projects/integrable/), 2015
[#eck_2005]: Matthias Eck, Tony DeRose, Tom Duchamp, Hugues Hoppe, Michael Lounsbery, Werner Stuetzle.  [Multiresolution Analysis of Arbitrary Meshes](http://research.microsoft.com/en-us/um/people/hoppe/mra.pdf), 2005.
[#garg_2016]: Akash Garg, Alec Jacobson, Eitan Grinspun. [Computational Design of Reconfigurables](https://www.google.com/search?q=Computational+Design+of+Reconfigurables), 2016
[#hildebrandt_2011]: Klaus Hildebrandt, Christian Schulz, Christoph von Tycowicz, and Konrad Polthier. [Interactive Surface Modeling using Modal Analysis](https://www.google.com/search?q=Interactive+Surface+Modeling+using+Modal+Analysis), 2011.
[#hoppe_1996]: Hugues Hoppe. [Progressive Meshes](https://www.google.com/search?q=Progressive+meshes), 1996
[#jacobson_skinning_course_2014]: Alec Jacobson, Zhigang Deng, Ladislav Kavan, J.P. Lewis. [_Skinning: Real-Time Shape Deformation_](https://www.google.com/search?q=Skinning+Real-Time+Shape+Deformation), 2014.
[#jacobson_thesis_2013]: Alec Jacobson, [_Algorithms and Interfaces for Real-Time Deformation of 2D and 3D Shapes_](https://www.google.com/search?q=Algorithms+and+Interfaces+for+Real-Time+Deformation+of+2D+and+3D+Shapes), 2013.
[#jacobson_2013]: Alec Jacobson, Ladislav Kavan, and Olga Sorkine. [Robust Inside-Outside Segmentation using Generalized Winding Numbers](https://www.google.com/search?q=Robust+Inside-Outside+Segmentation+using+Generalized+Winding+Numbers), 2013.
[#jacobson_2012]: Alec Jacobson, Ilya Baran, Ladislav Kavan, Jovan Popović, and Olga Sorkine. [Fast Automatic Skinning Transformations](https://www.google.com/search?q=Fast+Automatic+Skinning+Transformations), 2012.
[#jacobson_2011]: Alec Jacobson, Ilya Baran, Jovan Popović, and Olga Sorkine. [Bounded Biharmonic Weights for Real-Time Deformation](https://www.google.com/search?q=Bounded+biharmonic+weights+for+real-time+deformation), 2011.
[#jacobson_mixed_2010]: Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin. [Mixed Finite Elements for Variational Surface Modeling](https://www.google.com/search?q=Mixed+Finite+Elements+for+Variational+Surface+Modeling), 2010.
[#kavan_2008]: Ladislav Kavan, Steven Collins, Jiri Zara, and Carol O'Sullivan. [Geometric Skinning with Approximate Dual Quaternion Blending](https://www.google.com/search?q=Geometric+Skinning+with+Approximate+Dual+Quaternion+Blending), 2008.
[#kazhdan_2012]: Michael Kazhdan, Jake Solomon, Mirela Ben-Chen, [Can Mean-Curvature Flow Be Made Non-Singular](https://www.google.com/search?q=Can+Mean-Curvature+Flow+Be+Made+Non-Singular), 2012.
[#knoppel_2013]: Felix Knöppel, Keenan Crane, Ulrich Pinkall, and Peter Schröder. [Globally Optimal Direction Fields](http://www.cs.columbia.edu/~keenan/Projects/GloballyOptimalDirectionFields/paper.pdf), 2013.
[#levy_2002]: Bruno Lévy, Sylvain Petitjean, Nicolas Ray, Jérome Maillot. [Least Squares Conformal Maps, for Automatic Texture Atlas Generation,](http://www.cs.jhu.edu/~misha/Fall09/Levy02.pdf), 2002.
[#levy_2008]: Nicolas Ray, Bruno Vallet, Wan Chiu Li, Bruno Lévy. [N-Symmetry Direction Field Design](http://alice.loria.fr/publications/papers/2008/DGF/NSDFD-TOG.pdf), 2008.
[#liu_2008]: Ligang Liu, Lei Zhang, Yin Xu, Craig Gotsman, Steven J. Gortler. [A Local/Global Approach to Mesh Parameterization](http://cs.harvard.edu/~sjg/papers/arap.pdf), 2008.
[#liu_2011]: Yang Liu, Weiwei Xu, Jun Wang, Lifeng Zhu, Baining Guo, Falai Chen, Guoping Wang.  [General Planar Quadrilateral Mesh Design Using Conjugate Direction Field](http://research.microsoft.com/en-us/um/people/yangliu/publication/cdf.pdf), 2008.
[#loop_1987]: Charles Loop. [Smooth Subdivision Surfaces Based on Triangles](https://www.google.com/search?q=smooth+subdivision+surfaces+based+on+triangles), 1987.
[#lorensen_1987]: W.E. Lorensen and Harvey E. Cline. [Marching cubes: A high resolution 3d surface construction algorithm](https://www.google.com/search?q=Marching+cubes:+A+high+resolution+3d+surface+construction+algorithm), 1987.
[#mcadams_2011]: Alexa McAdams, Andrew Selle, Rasmus Tamstorf, Joseph Teran, Eftychios Sifakis. [Computing the Singular Value Decomposition of 3x3 matrices with minimal branching and elementary floating point operations](https://www.google.com/search?q=Computing+the+Singular+Value+Decomposition+of+3x3+matrices+with+minimal+branching+and+elementary+floating+point+operations), 2011.
[#meyer_2003]: Mark Meyer, Mathieu Desbrun, Peter Schröder and Alan H.  Barr, [Discrete Differential-Geometry Operators for Triangulated 2-Manifolds](https://www.google.com/search?q=Discrete+Differential-Geometry+Operators+for+Triangulated+2-Manifolds), 2003.
[#mullen_2008]: Patrick Mullen, Yiying Tong, Pierre Alliez, Mathieu Desbrun. [Spectral Conformal Parameterization](http://www.geometry.caltech.edu/pubs/MTAD08.pdf), 2008.
[#panozzo_2010]: Daniele Panozzo, Enrico Puppo, Luigi Rocca, [Efficient Multi-scale Curvature and Crease Estimation](https://www.google.com/search?q=Efficient+Multi-scale+Curvature+and+Crease+Estimation), 2010.
[#panozzo_2014]: Daniele Panozzo, Enrico Puppo, Marco Tarini, Olga Sorkine-Hornung.  [Frame Fields: Anisotropic and Non-Orthogonal Cross Fields](http://cs.nyu.edu/~panozzo/papers/frame-fields-2014.pdf), 2014.
[#rabinovich_2016]: Michael Rabinovich, Roi Poranne, Daniele Panozzo, Olga Sorkine-Hornung. [Scalable Locally Injective Mappings](http://cs.nyu.edu/~panozzo/papers/SLIM-2016.pdf), 2016.
[#rustamov_2011]: Raid M. Rustamov, [Multiscale Biharmonic Kernels](https://www.google.com/search?q=Multiscale+Biharmonic+Kernels), 2011.
[#schroeder_1994]: William J. Schroeder, William E. Lorensen, and Steve Linthicum. [Implicit Modeling of Swept Surfaces and Volumes](https://www.google.com/search?q=implicit+modeling+of+swept+surfaces+and+volumes), 1994.
[#schuller_2013]: Christian Schüller, Ladislav Kavan, Daniele Panozzo, Olga Sorkine-Hornung.  [Locally Injective Mappings](http://igl.ethz.ch/projects/LIM/), 2013.
[#sharf_2007]: Andrei Sharf, Thomas Lewiner, Gil Shklarski, Sivan Toledo, and Daniel Cohen-Or. [Interactive topology-aware surface reconstruction](https://www.google.com/search?q=Interactive+topology-aware+surface+reconstruction), 2007.
[#sorkine_2004]: Olga Sorkine, Yaron Lipman, Daniel Cohen-Or, Marc Alexa, Christian Rössl and Hans-Peter Seidel. [Laplacian Surface Editing](https://www.google.com/search?q=Laplacian+Surface+Editing), 2004.
[#sorkine_2007]: Olga Sorkine and Marc Alexa. [As-rigid-as-possible Surface Modeling](https://www.google.com/search?q=As-rigid-as-possible+Surface+Modeling), 2007.
[#stein_2017]: Oded Stein, Eitan Grinspun, Max Wardetzky, Alec Jacobson. [Natural Boundary Conditions for Smoothing in Geometry Processing](https://arxiv.org/abs/1707.04348), 2017.
[#takayama14]: Kenshi Takayama, Alec Jacobson, Ladislav Kavan, Olga Sorkine-Hornung. [A Simple Method for Correcting Facet Orientations in Polygon Meshes Based on Ray Casting](https://www.google.com/search?q=A+Simple+Method+for+Correcting+Facet+Orientations+in+Polygon+Meshes+Based+on+Ray+Casting), 2014.
[#vallet_2008]: Bruno Vallet and Bruno Lévy. [Spectral Geometry Processing with Manifold Harmonics](https://www.google.com/search?q=Spectral+Geometry+Processing+with+Manifold+Harmonics), 2008.
[#vaxman_2016]: Amir Vaxman, Marcel Campen, Olga Diamanti, Daniele Panozzo, David Bommes, Klaus Hildebrandt, Mirela Ben-Chen. [Directional Field Synthesis, Design, and Processing](https://www.google.com/search?q=Directional+Field+Synthesis+Design+and+Processing), 2016
[#wang_bc_2015]: Yu Wang, Alec Jacobson, Jernej Barbic, Ladislav Kavan. [Linear Subspace Design for Real-Time Shape Deformation](https://www.google.com/search?q=Linear+Subspace+Design+for+Real-Time+Shape+Deformation), 2015
[#zhou_2016]: Qingnan Zhou, Eitan Grinspun, Denis Zorin. [Mesh Arrangements for Solid Geometry](https://www.google.com/search?q=Mesh+Arrangements+for+Solid+Geometry), 2016
[#mitchell_1987]: Joseph S. B. Mitchell, David M. Mount, Christos H. Papadimitriou. [The Discrete Geodesic Problem](https://www.google.com/search?q=The+Discrete+Geodesic+Problem), 1987
