This directory contains a slew of models and data needed to run the tutorial
examples. This README.md contains an attempt to track down original sources of
this data for purposes of attribution. In some cases, only the oldest known
use is listed.

| Filename                           | Source             |
|------------------------------------|--------------------|
| 2triangles.off                     | [#libigl][]        |
| 3holes.off                         | ?                  |
| arm-weights.dmat                   | [#libigl][]        |
| arm.obj                            | ?                  |
| arm.tgf                            | [#libigl][]        |
| armadillo-weights.dmat             | [#libigl][]        |
| armadillo.obj                      | [#stanford][]      |
| beetle.off                         | [#ivan][]          |
| big-sigcat.mesh                    | [#jacobson_2013][] |
| bump-domain.obj                    | ?                  |
| bumpy-cube.dmat                    | ?                  |
| bumpy-cube.obj                     | ?                  |
| bumpy.off                          | ?                  |
| bunny.mesh                         | [#stanford][]      |
| bunny.off                          | [#stanford][]      |
| camelhead.off                      | [#sorkine_2004][]  |
| cheburashka-scalar.dmat            | [#libigl][]        |
| cheburashka.off                    | [#cosmic_blobs][]  |
| cow.off                            | ?                  |
| cube.obj                           | [#libigl][]        |
| cube.off                           | [#libigl][]        |
| decimated-knight-selection.dmat    | [#libigl][]        |
| decimated-knight.off               | [#cosmic_blobs][]  |
| decimated-max-selection.dmat       | [#libigl][]        |
| decimated-max.obj                  | [#mpi][]           |
| fandisk.off                        | [#aim_at_shape][]  |
| fertility.off                      | [#aim_at_shape][]  |
| grid.off                           | [#libigl][]        |
| hand-pose.dmat                     | [#libigl][]        |
| hand.mesh                          | ?                  |
| hand.tgf                           | [#libigl][]        |
| horse_quad.obj                     | ?                  |
| inspired_mesh.dmat                 | [#libigl][]        |
| inspired_mesh.obj                  | ?                  |
| inspired_mesh_b.dmat               | [#libigl][]        |
| inspired_mesh_bc.dmat              | [#libigl][]        |
| inspired_mesh_quads_Conjugate.off  | ?                  |
| inspired_mesh_quads_Smooth.off     | ?                  |
| lilium.crossfield                  | [#libigl][]        |
| lilium.obj                         | ?                  |
| lilium.samples.0.2                 | [#libigl][]        |
| lilium_b.dmat                      | [#libigl][]        |
| lilium_bc.dmat                     | [#libigl][]        |
| lion.off                           | ?                  |
| octopus-high.mesh                  | [#pauly][]         |
| octopus-low.mesh                   | [#pauly][]         |
| planexy.off                        | ?                  |
| screwdriver.off                    | ?                  |
| snail.dmat                         | [#libigl][]        |
| snail.obj                          | ?                  |
| snail.samples.0.2                  | [#libigl][]        |
| snail1.dmat                        | [#libigl][]        |
| snail2.dmat                        | [#libigl][]        |
| snail3.dmat                        | [#libigl][]        |
| snail4.dmat                        | [#libigl][]        |
| sphere.obj                         | ?                  |
| truck.obj                          | [#shrec][]         |
| xcylinder.obj                      | [#libigl][]        |
| ycylinder.obj                      | [#libigl][]        |
| zcylinder.obj                      | [#libigl][]        |

[#aim_at_shape]: AIM@SHAPE,IAL_SHARED_PATH "/cow.off".cnr.it:8080/ontologies/shapes/viewgroup.jsp?id=225-Fandisk_MC
[#cosmic_blobs]: Cosmic blobs, http://www.mit.edu/~ibaran/autorig/
[#ivan]: Hand measured from physical object, http://www.cs.utah.edu/docs/misc/Uteapot03.pdf
[#jacobson_2013]: Alec Jacobson, Ladislav Kavan, and Olga Sorkine.
  [Robust Inside-Outside Segmentation using Generalized Winding
  Numbers](https://www.google.com/search?q=Robust+Inside-Outside+Segmentation+using+Generalized+Winding+Numbers),
  2013.
[#libigl]: Original data produced by libigl authors, http://libigl.github.io/libigl/
[#mpi]: Max Planck Institute at Saarbrucken.
[#pauly]: Modeled by Mark Pauly, "Shape Modeling with Point-Sampled Geometry"
[#shrec]: SHREC 2009 Dataset,
  http://www.itl.nist.gov/iad/vug/sharp/benchmark/shrecPartial/data.html
[#sorkine_2004]: Olga Sorkine, Daniel Cohen-Or, [Least-squares
  Meshes](https://www.google.com/search?q=Least+squares+meshes), 2004.
[#stanford]: The Stanford 3D Scanning Repository, http://graphics.stanford.edu/data/3Dscanrep/
