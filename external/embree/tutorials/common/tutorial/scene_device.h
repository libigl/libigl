// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

struct ISPCTriangle 
{
  int v0;                /*< first triangle vertex */
  int v1;                /*< second triangle vertex */
  int v2;                /*< third triangle vertex */
  int materialID;        /*< material of triangle */
};

struct ISPCQuad
{
  int v0;                /*< first triangle vertex */
  int v1;                /*< second triangle vertex */
  int v2;                /*< third triangle vertex */
  int v4;                /*< fourth triangle vertex */
};

struct ISPCHair
{
  int vertex;
  int id;
};

struct ISPCHairSet
{
  Vec3fa* v;       //!< hair control points (x,y,z,r)
  Vec3fa* v2;       //!< hair control points (x,y,z,r)
  ISPCHair* hairs; //!< for each hair, index to first control point
  int numVertices;
  int numHairs;
};


struct ISPCMesh
{
  Vec3fa* positions;    //!< vertex position array
  Vec3fa* positions2;    //!< vertex position array
  Vec3fa* normals;       //!< vertex normal array
  Vec2f* texcoords;     //!< vertex texcoord array
  ISPCTriangle* triangles;  //!< list of triangles
  ISPCQuad* quads;  //!< list of triangles
  float* edge_level;
  int numVertices;
  int numTriangles;
  int numQuads;
  int geomID;
};

struct ISPCSubdivMesh
{
  Vec3fa* positions;       //!< vertex positions
  Vec3fa* normals;         //!< face vertex normals
  Vec2f* texcoords;        //!< face texture coordinates
  int* position_indices;   //!< position indices for all faces
  int* normal_indices;     //!< normal indices for all faces
  int* texcoord_indices;   //!< texcoord indices for all faces
  int* verticesPerFace;    //!< number of indices of each face
  int* holes;              //!< face ID of holes
  float* subdivlevel;      //!< subdivision level
  Vec2i* edge_creases;          //!< crease index pairs
  float* edge_crease_weights;   //!< weight for each crease
  int* vertex_creases;          //!< indices of vertex creases
  float* vertex_crease_weights; //!< weight for each vertex crease
  int numVertices;
  int numFaces;
  int numEdges;
  int numEdgeCreases;
  int numVertexCreases;
  int numHoles;
  int materialID;
  int geomID;
};

struct ISPCAmbientLight
{
  Vec3fa L;                  //!< radiance of ambient light
};

struct ISPCPointLight
{
  Vec3fa P;                  //!< position of point light
  Vec3fa I;                  //!< radiant intensity of point light
};

struct ISPCDirectionalLight
{
  Vec3fa D;                  //!< Light direction
  Vec3fa E;                  //!< Irradiance (W/m^2)
};

struct ISPCDistantLight
{
  Vec3fa D;             //!< Light direction
  Vec3fa L;             //!< Radiance (W/(m^2*sr))
  float halfAngle;     //!< Half illumination angle
  float radHalfAngle;  //!< Half illumination angle in radians
  float cosHalfAngle;  //!< Cosine of half illumination angle
};

enum MaterialTy { MATERIAL_OBJ, MATERIAL_THIN_DIELECTRIC, MATERIAL_METAL, MATERIAL_VELVET, MATERIAL_DIELECTRIC, MATERIAL_METALLIC_PAINT, MATERIAL_MATTE, MATERIAL_MIRROR, MATERIAL_REFLECTIVE_METAL };

struct ISPCMaterial
{
  int ty;
  int align0,align1,align2;
  Vec3fa v[7];
};

struct MatteMaterial
{
  int ty;
  int align[3];
  Vec3fa reflectance;
};

struct MirrorMaterial
{
  int ty;
  int align[3];
  Vec3fa reflectance;
};

struct OBJMaterial
{
  int ty;
  int align[3];

  int illum;             /*< illumination model */
  float d;               /*< dissolve factor, 1=opaque, 0=transparent */
  float Ns;              /*< specular exponent */
  float Ni;              /*< optical density for the surface (index of refraction) */
  
  Vec3fa Ka;              /*< ambient reflectivity */
  Vec3fa Kd;              /*< diffuse reflectivity */
  Vec3fa Ks;              /*< specular reflectivity */
  Vec3fa Kt;              /*< transmission filter */
  Vec3fa v[2];
};

struct MetalMaterial
{
  int ty;
  int align[3];

  Vec3fa reflectance;
  Vec3fa eta;
  Vec3fa k;
  float roughness;
};

struct ReflectiveMetalMaterial
{
  int ty;
  int align[3];

  Vec3fa reflectance;
  Vec3fa eta;
  Vec3fa k;
  float roughness;
};

struct VelvetMaterial
{
  int ty;
  int align[3];

  Vec3fa reflectance;
  Vec3fa horizonScatteringColor;
  float backScattering;
  float horizonScatteringFallOff;
};

struct DielectricMaterial
{
  int ty;
  int align[3];
  Vec3fa transmissionOutside;
  Vec3fa transmissionInside;
  float etaOutside;
  float etaInside;
};


struct ThinDielectricMaterial
{
  int ty;
  int align[3];
  Vec3fa transmission;
  float eta;
};

struct MetallicPaintMaterial
{
  int ty;
  int align[3];
  Vec3fa shadeColor;
  Vec3fa glitterColor;
  float glitterSpread;
  float eta;
};

struct ISPCScene {

  ISPCMesh** meshes;   //!< list of meshes
  ISPCMaterial* materials;     //!< material list
  int numMeshes;                       //!< number of meshes
  int numMaterials;                    //!< number of materials

  ISPCHairSet** hairs;
  int numHairSets;

  ISPCAmbientLight* ambientLights; //!< list of ambient lights
  int numAmbientLights;                    //!< number of ambient lights
  
  ISPCPointLight* pointLights;     //!< list of point lights
  int numPointLights;                      //!< number of point lights
  
  ISPCDirectionalLight* dirLights; //!< list of directional lights
  int numDirectionalLights;                //!< number of directional lights

  ISPCDistantLight* distantLights; //!< list of distant lights
  int numDistantLights;                    //!< number of distant lights

  ISPCSubdivMesh** subdiv;                   //!< list of subdiv meshes
  int numSubdivMeshes;                       //!< number of subdiv meshes

}; // ISPCScene



