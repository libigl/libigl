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

#include "tutorial/tutorial.h"
#include "tutorial/obj_loader.h"
#include "tutorial/hair_loader.h"
#include "tutorial/cy_hair_loader.h"
#include "sys/taskscheduler.h"
#include "image/image.h"

extern "C" embree::Vec3fa g_dirlight_direction = embree::normalize(embree::Vec3fa(1,-1,1));
extern "C" embree::Vec3fa g_dirlight_intensity = embree::Vec3fa(4.0f);
extern "C" embree::Vec3fa g_ambient_intensity = embree::Vec3fa(1.0f);

namespace embree
{
  /* name of the tutorial */
  const char* tutorialName = "tutorial07";

  /* configuration */
  static std::string g_rtcore = "";

  /* output settings */
  static size_t g_width = 512;
  static size_t g_height = 512;
  static bool g_fullscreen = false;
  static size_t g_numThreads = 0;

  static int tessellate_subdivisions = -1;
  static int tessellate_strips = -1;
  extern float g_reduce_hair_segment_error;

  /* scene */
  OBJScene g_obj_scene;
  OBJScene g_obj_scene2;
  static FileName objFilename = "";
  static FileName objFilename2 = "";
  static FileName hairFilename = "";
  static FileName hairFilename2 = "";
  static FileName cy_hairFilename = "";
  static FileName outFilename = "";
  static int g_skipBenchmarkFrames = 0;
  static int g_numBenchmarkFrames = 0;
  static bool g_interactive = true;

  static bool hairy_triangles = false;
  static float hairy_triangles_length = 1.0f;
  static float hairy_triangles_thickness = 0.1f;
  static int hairy_triangles_strands_per_triangle = 1;

  Vec3fa offset = 0.0f;
  Vec3fa offset_mb = 0.0f;

  void addHairSegment(OBJScene::Mesh& mesh, int materialID, const Vec3fa& p0, const Vec3fa& p1)
  {
    const size_t N = tessellate_strips;
    if (length(p1-p0) <= 1E-6f) return;
    const LinearSpace3fa xfm = frame(p1-p0);
    const int base = mesh.v.size();
    for (size_t i=0; i<N; i++) {
      const float a = 2.0f*float(pi)*float(i)/float(N);
      const float u = sinf(a), v = cosf(a);
      mesh.v.push_back(p0+u*p0.w*xfm.vx+v*p0.w*xfm.vy);
      mesh.v.push_back(p1+u*p1.w*xfm.vx+v*p1.w*xfm.vy);
      mesh.vn.push_back(p1-p0);
      mesh.vn.push_back(p1-p0);
      mesh.vt.push_back(Vec2f(max(p0.w,p1.w)));
      mesh.vt.push_back(Vec2f(max(p0.w,p1.w)));
    }
    for (size_t i=0; i<N; i++) {
      const int v0 = base + (2*i+0)%(2*N);
      const int v1 = base + (2*i+1)%(2*N);
      const int v2 = base + (2*i+2)%(2*N);
      const int v3 = base + (2*i+3)%(2*N);
      mesh.triangles.push_back(OBJScene::Triangle(v0,v1,v2,materialID));
      mesh.triangles.push_back(OBJScene::Triangle(v1,v3,v2,materialID));
    }
  }

  void tessellateHair(OBJScene::Mesh& mesh,
                      int materialID, 
                      const Vec3fa& p00,
                      const Vec3fa& p01,
                      const Vec3fa& p02,
                      const Vec3fa& p03,
                      const int depth)
  {
    if (depth > 0) 
    {
      const Vec3fa p10 = 0.5f*(p00+p01);
      const Vec3fa p11 = 0.5f*(p01+p02);
      const Vec3fa p12 = 0.5f*(p02+p03);
      const Vec3fa p20 = 0.5f*(p10+p11);
      const Vec3fa p21 = 0.5f*(p11+p12);
      const Vec3fa p30 = 0.5f*(p20+p21);
      tessellateHair(mesh,materialID,p00,p10,p20,p30,depth-1);
      tessellateHair(mesh,materialID,p30,p21,p12,p03,depth-1);
    } else {
      addHairSegment(mesh,materialID,p00,p03);
    }
  }

  void tessellateHair(OBJScene::Mesh& mesh, int materialID, const OBJScene::HairSet& hairSet)
  {
    for (size_t i=0; i<hairSet.hairs.size(); i++) 
    {
      OBJScene::Hair hair = hairSet.hairs[i];
      const Vec3fa p00 = hairSet.v[hair.vertex+0];
      const Vec3fa p01 = hairSet.v[hair.vertex+1];
      const Vec3fa p02 = hairSet.v[hair.vertex+2];
      const Vec3fa p03 = hairSet.v[hair.vertex+3];
      tessellateHair(mesh,materialID,p00,p01,p02,p03,tessellate_subdivisions);
    }
  }

  void tessellateHair(OBJScene& scene)
  {
    OBJScene::OBJMaterial material; material.illum = 1;
    int materialID = scene.materials.size();
    scene.materials.push_back(material);

    for (int i=0; i<scene.hairsets.size(); i++) 
    {
      OBJScene::Mesh* mesh = new OBJScene::Mesh;
      OBJScene::HairSet* hairset = scene.hairsets[i];
      tessellateHair(*mesh,materialID,*hairset);
      scene.meshes.push_back(mesh);
      delete hairset;
    }
    scene.hairsets.clear();
  }

  void generateHairOnTriangleMesh(OBJScene& scene, OBJScene::Mesh* mesh, float hair_length, float thickness, size_t strands_per_triangle)
  {
    OBJScene::HairSet* hairset = new OBJScene::HairSet;

    for (size_t t=0; t<mesh->triangles.size(); t++) 
    {
      const size_t i0 = mesh->triangles[t].v0;
      const size_t i1 = mesh->triangles[t].v1;
      const size_t i2 = mesh->triangles[t].v2;
      
      const Vec3fa& v0 = mesh->v[ i0 ];
      const Vec3fa& v1 = mesh->v[ i1 ];
      const Vec3fa& v2 = mesh->v[ i2 ];
      
      Vec3fa n0,n1,n2;
      if (mesh->vn.size() != 0)
      {
        n0 = mesh->vn[ i0 ];
        n1 = mesh->vn[ i1 ];
        n2 = mesh->vn[ i2 ];
      }
      else
      {
        const Vec3fa edge0 = v1-v0;
        const Vec3fa edge1 = v2-v0;
        Vec3fa normal = cross(edge0,edge1);
        n0 = n1 = n2 = normal;
      }
      
      //if (length(normal) < 1E-30) continue;
      
      //const float thickness = thickness;
      
      for (size_t i=0; i<strands_per_triangle; i++)
      {
        float ru = drand48();
        float rv = drand48();
        float delta = 0.1f*drand48();
        float delta2 = 0.1f*drand48();
	
        const float w0 = 1.0f - sqrtf(ru);
        const float w1 = sqrtf(ru) * (1.0f - rv);
        const float w2 = sqrtf(ru) * rv;
        const Vec3fa position = w0*v0 + w1*v1 + w2*v2;
        const Vec3fa normal = w0*n0 + w1*n1 + w2*n2;
        if (length(normal) < 1E-30) continue;
        
        const LinearSpace3fa local_frame = frame(normalize(normal));
        const float length = hair_length;
        const Vec3fa dx = local_frame.vx * length;
        const Vec3fa dy = local_frame.vy * length;
        const Vec3fa dz = local_frame.vz * length;
	
        const Vec3fa p0(   0, 0,0);
        //const Vec3fa p1(0.25, 0,0);
        //const Vec3fa p2(0.74, 0,0);
        const Vec3fa p1(0.25, delta, delta2);
        const Vec3fa p2(0.75,-delta,-delta2);
        const Vec3fa p3(   1, 0,0);
	
        Vec3fa l0 = position + p0.x * dz + p0.y * dx + p0.z * dy; l0.w = thickness;
        Vec3fa l1 = position + p1.x * dz + p1.y * dx + p1.z * dy; l1.w = thickness;
        Vec3fa l2 = position + p2.x * dz + p2.y * dx + p2.z * dy; l2.w = thickness;
        Vec3fa l3 = position + p3.x * dz + p3.y * dx + p3.z * dy; l3.w = thickness;
        
        const unsigned int v_index = hairset->v.size();
        hairset->v.push_back(l0);
        hairset->v.push_back(l1);
        hairset->v.push_back(l2);
        hairset->v.push_back(l3);
	
        hairset->hairs.push_back( OBJScene::Hair(v_index,hairset->hairs.size()) );
      }
    }
    scene.hairsets.push_back(hairset);
  }

  void generateHairOnTriangles(OBJScene& scene)
  {
    for (size_t m=0; m<scene.meshes.size(); m++) 
      generateHairOnTriangleMesh(scene,scene.meshes[m],hairy_triangles_length,hairy_triangles_thickness,hairy_triangles_strands_per_triangle);
  }

  Vec3fa uniformSampleSphere(const float& u, const float& v) 
  {
    const float phi = float(two_pi) * u;
    const float cosTheta = 1.0f - 2.0f * v, sinTheta = 2.0f * sqrt(v * (1.0f - v));
    return Vec3fa(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
  }

static int p[513] = { 
  151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
  190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
  88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
  77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
  102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
  135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
  5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
  223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
  129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
  251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
  49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
  138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
  151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
  190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
  88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
  77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
  102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
  135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
  5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
  223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
  129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
  251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
  49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
  138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
  151
};

static float g1[128] = {
  -0.20707, 0.680971, -0.293328, -0.106833, -0.362614, 0.772857, -0.968834, 0.16818, -0.681263, -0.232568, 0.382009, -0.882282, 0.799709, -0.672908, -0.681857, 0.0661294, 0.208288, 0.165398, -0.460058, -0.219044, -0.413199, 0.484755, -0.402949, -0.848924, -0.190035, 0.714756, 0.883937, 0.325661, 0.692952, -0.99449, -0.0752415, 0.065192, 0.575753, -0.468776, 0.965505, -0.38643, 0.20171, 0.217431, -0.575122, 0.77179, -0.390686, -0.69628, -0.324676, -0.225046, 0.28722, 0.507107, 0.207232, 0.0632565, -0.0812794, 0.304977, -0.345638, 0.892741, -0.26392, 0.887781, -0.985143, 0.0331999, -0.454458, -0.951402, 0.183909, -0.590073, 0.755387, -0.881263, -0.478315, -0.394342, 0.78299, -0.00360388, 0.420051, -0.427172, 0.729847, 0.351081, -0.0830201, 0.919271, 0.549351, -0.246897, -0.542722, -0.290932, -0.399364, 0.339532, 0.437933, 0.131909, 0.648931, -0.218776, 0.637533, 0.688017, -0.639064, 0.886792, -0.150226, 0.0413316, -0.868712, 0.827016, 0.765169, 0.522728, -0.202155, 0.376514, 0.523097, -0.189982, -0.749498, -0.0307322, -0.555075, 0.746242, 0.0576438, -0.997172, 0.721028, -0.962605, 0.629784, -0.514232, -0.370856, 0.931465, 0.87112, 0.618863, -0.0157817, -0.559729, 0.152707, -0.421942, -0.357866, -0.477353, -0.652024, -0.996365, -0.910432, -0.517651, -0.169098, 0.403249, -0.556309, 0.00782069, -0.86594, -0.213873, -0.0410469, -0.563716
};

static float g2[128*2] = {
0.605609, 0.538399, 0.796519, -0.944204, 0.908294, 0.756016, 0.0977536, -0.863638, 0.842196, -0.744751, -0.932081, 0.932392, -0.588525, 0.516884, 0.841188, -0.978497, -0.608649, -0.868011, 0.992137, -0.772425, 0.963049, -0.0478757, 0.953878, 0.889467, 0.562174, 0.624644, -0.356598, -0.520726, -0.821833, 0.99985, 0.234183, -0.9791, -0.971815, -0.0979374, -0.108159, -0.34927, -0.592124, -0.775632, 0.97228, 0.753819, 0.941608, 0.578291, 0.852108, -0.760312, -0.784772, 0.0223242, -0.606013, -0.980319, 0.252581, -0.575064, 0.884701, 0.943763, 0.737344, 0.938496, 0.0466562, -0.994566, 0.989782, 0.988368, -0.546155, 0.279211, -0.69504, 0.931229, 0.99768, -0.325874, -0.630157, -0.999936, -0.968623, -0.226805, -0.750428, -0.450961, 0.257868, 0.968011, -0.988005, -0.713965, 0.991007, -0.61059, 0.950437, -0.483042, -0.98105, -0.915356, -0.892527, -0.772958, -0.9081, 0.55692, 0.906075, 0.937419, 0.454624, -0.991582, 0.400857, 0.855933, -0.672619, 0.0713424, 0.593249, -0.378286, -0.997369, -0.827112, 0.708222, -0.995343, 0.985069, 0.698711, -0.180105, 0.999961, -0.768451, 0.993107, -0.918024, 0.0446961, 0.91882, 0.97691, -0.393915, 0.364803, 0.0495592, 0.186545, -0.461553, -0.242776, 0.901952, -0.0710866, 0.888101, 0.999935, 0.277688, 0.0554235, 0.506599, -0.299293, 0.984394, -0.999698, 0.408822, -0.782639, 0.128596, 0.198834, 0.981707, 0.864566, 0.808197, 0.352335, 0.970484, -0.667503, 0.330243, 0.208392, 0.191539, -0.938943, 0.895002, 0.910575, -0.537691, -0.98548, -0.721635, -0.335382, -0.424701, -0.960452, 0.595047, 0.783579, -0.937749, 0.529096, -0.997906, -0.581313, -0.899828, -0.88461, 0.989469, 0.91872, -0.850793, 0.955954, 0.715768, -0.736686, 0.80392, -0.717276, -0.788579, 0.987003, -0.839648, 0.885176, -0.998929, -0.0376033, -0.578371, -0.718771, 0.906081, 0.239947, -0.803563, -0.00826282, 0.991011, -0.0057943, -0.349232, 0.65319, 0.992067, -0.953535, 0.893781, 0.661689, 0.957253, -0.425442, -0.866609, 0.712892, -0.807777, 0.89632, -0.595147, -0.0224999, -0.643786, 0.545815, -0.870124, -0.696306, -0.99902, 0.773648, -0.806008, -0.931319, -0.780114, -0.552154, -0.933812, -0.563108, -0.619909, 0.966532, 0.692454, 0.993284, 0.338885, -0.75104, 0.237272, -0.713619, -0.160187, -0.199242, -0.371265, -0.781439, -0.914125, -0.944104, 0.169525, -0.984403, 0.976056, -0.265228, 0.94232, 0.993906, -0.877517, -0.89618, 0.611817, -0.106758, 0.680403, 0.163329, -0.325386, -0.0687362, -0.901164, 0.460314, 0.999981, -0.0408026, 0.850356, -0.763343, -0.170806, -0.102919, 0.581564, 0.688634, 0.284368, -0.276419, 0.616641, -0.929771, 0.927865, 0.440373, 0.153446, 0.840456, 0.996966, 0.867209, -0.135077, -0.493238, -0.577193, 0.0588088, 0.715215, 0.0143633
};

static float g3[128*4] = {
  -0.582745, 0.443494, -0.680971, 0, -0.601153, 0.791961, 0.106833, 0, -0.265466, 0.576385, -0.772857, 0, 0.981035, 0.0963612, -0.16818, 0, 0.524388, 0.819103, 0.232568, 0, -0.170518, -0.43875, 0.882282, 0, 0.598053, -0.435348, 0.672908, 0, 0.53956, 0.839346, -0.0661294, 0, -0.782511, -0.600267, -0.165398, 0, -0.122114, 0.968043, 0.219044, 0, -0.235567, 0.842331, -0.484755, 0, -0.158657, 0.504139, 0.848924, 0, -0.578396, 0.39317, -0.714756, 0, 0.883328, -0.337159, -0.325661, 0, 0.0597264, -0.0861552, 0.99449, 0, -0.970124, 0.233685, -0.0651921, 0, 0.208238, -0.858421, 0.468776, 0, 0.916908, -0.0997567, 0.38643, 0, -0.786568, -0.577957, -0.217431, 0, 0.14868, 0.618251, -0.77179, 0, -0.24168, 0.675858, 0.69628, 0, -0.50994, 0.83025, 0.225046, 0, -0.534183, -0.676382, -0.507107, 0, -0.793861, -0.6048, -0.0632565, 0, -0.92148, 0.240548, -0.304977, 0, -0.210037, 0.39862, -0.892741, 0, -0.310918, 0.339375, -0.887781, 0, 0.99836, 0.0466305, -0.0331999, 0, -0.0439099, 0.304806, 0.951402, 0, -0.676304, -0.440938, 0.590073, 0, 0.339805, -0.328495, 0.881263, 0, -0.0625568, 0.916832, 0.394342, 0, 0.776463, -0.630153, 0.00360388, 0, -0.224717, -0.8758, 0.427172, 0, 0.618879, -0.70266, -0.351081, 0, -0.380313, 0.101503, -0.919271, 0, 0.149639, -0.957418, 0.246897, 0, 0.128024, 0.948139, 0.290932, 0, -0.292448, 0.893976, -0.339532, 0, -0.192062, -0.972477, -0.131909, 0, 0.44007, -0.870905, 0.218776, 0, 0.303887, -0.659003, -0.688017, 0, 0.195552, 0.41876, -0.886792, 0, -0.889922, 0.454236, -0.0413315, 0, 0.515034, 0.225353, -0.827016, 0, 0.63084, -0.573408, -0.522728, 0, -0.745779, 0.549592, -0.376514, 0, 0.0711763, -0.979204, 0.189982, 0, 0.705657, 0.707887, 0.0307322, 0, 0.114603, 0.655735, -0.746242, 0, -0.0739232, -0.0135353, 0.997172, 0, 0.173356, -0.20818, 0.962605, 0, 0.34008, -0.787344, 0.514232, 0, -0.143596, 0.334295, -0.931465, 0, 0.721989, -0.30942, -0.618863, 0, -0.827657, 0.0410685, 0.559729, 0, -0.804277, -0.418454, 0.421942, 0, -0.379459, 0.792556, 0.477353, 0, 0.0391537, 0.0756503, 0.996365, 0, 0.821943, 0.237588, 0.517651, 0, -0.788974, 0.463584, -0.403249, 0, 0.175972, 0.984364, -0.00782073, 0, 0.891497, 0.399363, 0.213873, 0, -0.819111, 0.106216, 0.563716, 0, 0.105511, 0.544028, -0.832406, 0, -0.464551, 0.63753, 0.614612, 0, 0.232387, 0.935154, -0.267363, 0, 0.777619, 0.272068, -0.566823, 0, 0.975331, 0.190338, 0.111807, 0, 0.224313, 0.450072, -0.86436, 0, 0.841897, -0.536898, 0.0543103, 0, 0.637123, -0.664145, -0.391135, 0, 0.901675, -0.422984, 0.0898189, 0, -0.496241, 0.367413, -0.786608, 0, -0.255468, -0.689763, -0.677469, 0, -0.0616459, -0.951141, -0.302539, 0, -0.431011, -0.889035, -0.154425, 0, -0.0711688, 0.486502, -0.870776, 0, -0.223359, -0.36162, 0.905175, 0, -0.678546, 0.695482, -0.23639, 0, 0.576553, 0.77934, 0.245389, 0, -0.194568, -0.24951, 0.948624, 0, 0.28962, -0.447736, 0.845962, 0, -0.0403821, -0.871893, 0.488028, 0, 0.790972, -0.560788, 0.244705, 0, -0.34553, 0.739953, 0.57713, 0, -0.516376, -0.697122, 0.49737, 0, 0.115998, 0.859293, 0.498156, 0, 0.643831, -0.239955, 0.72657, 0, -0.125114, 0.987348, -0.0974144, 0, -0.306452, 0.610699, -0.73016, 0, -0.269845, 0.893027, -0.360119, 0, 0.328563, -0.570628, -0.752615, 0, -0.306918, -0.42057, 0.853769, 0, 0.699245, -0.51785, 0.492837, 0, -0.558362, -0.469763, -0.68378, 0, 0.476563, -0.841398, 0.254826, 0, 0.0276172, -0.623206, 0.78157, 0, 0.587723, -0.800313, -0.118659, 0, 0.594035, -0.740708, 0.313806, 0, -0.340185, -0.887929, 0.309605, 0, 0.312245, -0.246681, -0.917416, 0, 0.194206, 0.186398, -0.963089, 0, 0.915704, 0.329835, -0.229553, 0, 0.94133, 0.229917, 0.247055, 0, -0.888253, -0.144148, 0.436152, 0, -0.906917, -0.362625, -0.214486, 0, 0.403108, -0.908884, 0.10693, 0, 0.983963, 0.169256, 0.056292, 0, -0.197949, 0.888236, 0.414553, 0, 0.0879741, 0.247673, 0.964841, 0, 0.474384, -0.868071, -0.146331, 0, 0.699884, 0.541342, -0.465953, 0, 0.610965, 0.567249, 0.552223, 0, 0.830508, -0.285788, -0.478103, 0, 0.328573, -0.683076, -0.652263, 0, -0.00537775, 0.873381, 0.487009, 0, -0.51289, 0.828835, 0.223557, 0, -0.871168, -0.15102, 0.467182, 0, -0.545561, 0.390016, -0.741789, 0, 0.874063, 0.259258, 0.410852, 0, -0.781555, 0.612184, -0.120005, 0, -0.284928, 0.708938, -0.645154, 0, -0.568809, 0.0883274, 0.817713, 0, -0.0429388, 0.549957, -0.834088, 0, 0.933296, -0.127233, 0.335813, 0, 0.698149, -0.493464, 0.51873, 0, -0.603413, 0.617495, -0.504572, 0 
};

__forceinline float fade(float t) {
  return (t * t * t) * (t * (t * 6 - 15) + 10); 
}

/*__forceinline float fade(float t) {
  return t * t * (3 - t * 2); 
  }*/

__forceinline float lerp(float t, float a, float b) { 
  return a + t * (b - a); 
}

__forceinline float grad(int hash, float x) {
  return x*g1[hash&127];
}

__forceinline float grad(int hash, float x, float y) {
  int h = hash&127;
  return x*g2[2*h+0]+y*g2[2*h+1];
}

__forceinline float grad(int hash, float x, float y, float z) {
  int h = hash&127;
  return x*g3[4*h+0]+y*g3[4*h+1]+z*g3[4*h+2];
}

float noise(float x) 
{
  float fx = floorf(x);
  int X = (int)fx & 255;
  x -= fx;
  float u = fade(x);
  float g00 = grad(p[X  ],x  );
  float g10 = grad(p[X+1],x-1);
  return lerp(u,g00,g10);
}

float noise(float x, float y) 
{
  float fx = floorf(x);
  float fy = floorf(y);
  
  int X = (int)fx & 255;
  int Y = (int)fy & 255;
  
  x -= fx;
  y -= fy;
  
  float u = fade(x);
  float v = fade(y);
  
  int p00  = p[X  ]+Y;
  int p10  = p[X+1]+Y;
  int p01  = p[X  ]+Y+1;
  int p11  = p[X+1]+Y+1;
  
  float g00 = grad(p[p00],x  ,y  );
  float g10 = grad(p[p10],x-1,y  );
  float g01 = grad(p[p01],x  ,y-1);
  float g11 = grad(p[p11],x-1,y-1);
  
  return lerp(v,lerp(u,g00,g10),lerp(u,g01,g11));
}

float noise(float x, float y, float z) 
{
  float fx = floorf(x);
  float fy = floorf(y);
  float fz = floorf(z);
  
  int X = (int)fx & 255;
  int Y = (int)fy & 255;
  int Z = (int)fz & 255;
  
  x -= fx;
  y -= fy;
  z -= fz;
  
  float u = fade(x);
  float v = fade(y);
  float w = fade(z);
  
  int p00  = p[X]+Y;
  int p000 = p[p00]+Z;
  int p010 = p[p00+1]+Z;
  int p001 = p000+1; 
  int p011 = p010+1;
  int p10  = p[X+1]+Y;
  int p100 = p[p10]+Z;
  int p110 = p[p10+1]+Z;
  int p101 = p100+1;
  int p111 = p110+1;
  
  float g000 = grad(p[p000],x  ,y  ,z  );
  float g100 = grad(p[p100],x-1,y  ,z  );
  float g010 = grad(p[p010],x  ,y-1,z  );
  float g110 = grad(p[p110],x-1,y-1,z  );
  float g001 = grad(p[p001],x  ,y  ,z-1);
  float g101 = grad(p[p101],x-1,y  ,z-1);
  float g011 = grad(p[p011],x  ,y-1,z-1);
  float g111 = grad(p[p111],x-1,y-1,z-1);
  
  return lerp(w,
              lerp(v,lerp(u,g000,g100),lerp(u,g010,g110)),
              lerp(v,lerp(u,g001,g101),lerp(u,g011,g111)));
}

  Vec3fa noise3D(const Vec3fa& p)
  {
    float x = noise(4.0f*p.x);
    float y = noise(4.0f*p.y);
    float z = noise(4.0f*p.z);
    return p+0.2f*Vec3fa(x,y,z);
  }

  void addHairySphere (OBJScene& scene, const Vec3fa& p, float r)
  {
    const size_t numPhi   = 20;
    const size_t numTheta = 2*numPhi;
    OBJScene::Mesh* mesh = new OBJScene::Mesh;

    OBJScene::Material material;
    int materialID = scene.materials.size();
    scene.materials.push_back(material);
    
    const float rcpNumTheta = rcp((float)numTheta);
    const float rcpNumPhi   = rcp((float)numPhi);
    for (int phi=0; phi<=numPhi; phi++)
    {
      for (int theta=0; theta<numTheta; theta++)
      {
        const float phif   = phi*float(pi)*rcpNumPhi;
        const float thetaf = theta*2.0f*float(pi)*rcpNumTheta;
        const Vec3fa dp(sin(phif)*sin(thetaf),cos(phif),sin(phif)*cos(thetaf));
        mesh->v. push_back(p+r*dp);
        mesh->vn.push_back(dp);
      }
      if (phi == 0) continue;
      
      for (int theta=1; theta<=numTheta; theta++) 
      {
        int p00 = (phi-1)*numTheta+theta-1;
        int p01 = (phi-1)*numTheta+theta%numTheta;
        int p10 = phi*numTheta+theta-1;
        int p11 = phi*numTheta+theta%numTheta;

        if (phi > 1)
          mesh->triangles.push_back(OBJScene::Triangle(p10,p00,p01,materialID));
        
        if (phi < numPhi) 
          mesh->triangles.push_back(OBJScene::Triangle(p11,p10,p01,materialID));
      }
    }
    scene.meshes.push_back(mesh);
    //generateHairOnTriangleMesh(scene,mesh,0.5f*r,0.001f*r,80);

#if 0
    const float thickness = 0.01f*r;
    OBJScene::HairSet* hairset = new OBJScene::HairSet;
    srand48(123456789);
    for (size_t t=0; t<16; t++) 
      {
	Vec3fa dp = uniformSampleSphere(drand48(),drand48());

	Vec3fa l0 = p + r*       (dp + 0.00f*dp); l0.w = thickness;
	Vec3fa l1 = p + r*       (dp + 0.25f*dp); l1.w = thickness;
	Vec3fa l2 = p + r*noise3D(dp + 0.50f*dp); l2.w = thickness;
	Vec3fa l3 = p + r*noise3D(dp + 0.75f*dp); l3.w = thickness;

	const unsigned int v_index = hairset->v.size();
	hairset->v.push_back(l0);
	hairset->v.push_back(l1);
	hairset->v.push_back(l2);
	hairset->v.push_back(l3);

	hairset->hairs.push_back( OBJScene::Hair(v_index,hairset->hairs.size()) );
      }
    scene.hairsets.push_back(hairset);
#else
    const float thickness = 0.001f*r;
    OBJScene::HairSet* hairset = new OBJScene::HairSet;

    for (size_t t=0; t<100000; t++) 
    {
      Vec3fa dp = uniformSampleSphere(drand48(),drand48());

      Vec3fa l0 = p + r*       (dp + 0.00f*dp); l0.w = thickness;
      Vec3fa l1 = p + r*       (dp + 0.25f*dp); l1.w = thickness;
      Vec3fa l2 = p + r*noise3D(dp + 0.50f*dp); l2.w = thickness;
      Vec3fa l3 = p + r*noise3D(dp + 0.75f*dp); l3.w = thickness;
        
      const unsigned int v_index = hairset->v.size();
      hairset->v.push_back(l0);
      hairset->v.push_back(l1);
      hairset->v.push_back(l2);
      hairset->v.push_back(l3);
	
      hairset->hairs.push_back( OBJScene::Hair(v_index,hairset->hairs.size()) );
    }
    scene.hairsets.push_back(hairset);
#endif
  }

  void addGroundPlane (OBJScene& scene, const Vec3fa& p00, const Vec3fa& p01, const Vec3fa& p10, const Vec3fa& p11)
  {
    OBJScene::Mesh* mesh = new OBJScene::Mesh;

    OBJScene::Material material;
    int materialID = scene.materials.size();
    scene.materials.push_back(material);

    mesh->v.push_back(p00);
    mesh->v.push_back(p01);
    mesh->v.push_back(p10);
    mesh->v.push_back(p11);

    mesh->triangles.push_back(OBJScene::Triangle(0,1,2,materialID));
    mesh->triangles.push_back(OBJScene::Triangle(2,1,3,materialID));

    scene.meshes.push_back(mesh);
  }

  static void parseCommandLine(Ref<ParseStream> cin, const FileName& path)
  {
    while (true)
    {
      std::string tag = cin->getString();
      if (tag == "") return;

      /* parse command line parameters from a file */
      if (tag == "-c") {
        FileName file = path + cin->getFileName();
        parseCommandLine(new ParseStream(new LineCommentFilter(file, "#")), file.path());
      }

      /* load OBJ model */
      else if (tag == "-i") {
        objFilename = path + cin->getFileName();
      }

      /* load motion blur OBJ model */
      else if (tag == "-i_mb") {
        objFilename = path + cin->getFileName();
        objFilename2 = path + cin->getFileName();
      }

      /* load hair model */
      else if (tag == "--hair") {
        hairFilename = path + cin->getFileName();
      }

      /* motion blur hair model */
      else if (tag == "--hair_mb") {
        hairFilename = path + cin->getFileName();
        hairFilename2 = path + cin->getFileName();
      }

      /* load hair model */
      else if (tag == "--cy_hair") {
        cy_hairFilename = path + cin->getFileName();
      }

      /* scene offset */
      else if (tag == "--offset") {
        offset = cin->getVec3fa();
      }

      /* scene offset */
      else if (tag == "--offset_mb") {
        offset_mb = cin->getVec3fa();
      }

      /* directional light */
      else if (tag == "--dirlight") {
        g_dirlight_direction = normalize(cin->getVec3fa());
        g_dirlight_intensity = cin->getVec3fa();
      }

      /* ambient light */
      else if (tag == "--ambient") {
        g_ambient_intensity = cin->getVec3fa();
      }

      /* tessellation flags */
      else if (tag == "--tessellate-hair") {
        tessellate_subdivisions = cin->getInt();
        tessellate_strips = cin->getInt();
      }

      /* create hairy spheres */
      else if (tag == "--hairy-sphere") {
        const Vec3fa p = cin->getVec3fa();
        const float  r = cin->getFloat();
        addHairySphere(g_obj_scene,p,r);
      }

      /* create hair on geometry */
      else if (tag == "--hairy-triangles") {
        hairy_triangles_strands_per_triangle = cin->getInt();
        hairy_triangles_length = cin->getFloat();
        hairy_triangles_thickness = cin->getFloat();
        hairy_triangles = true;
      }

      /* reduce number of hair segments */
      else if (tag == "--reduce-hair-segment-error") {
        g_reduce_hair_segment_error = cin->getFloat();
      }

      /* output filename */
      else if (tag == "-o") {
        outFilename = cin->getFileName();
	g_interactive = false;
      }

      /* number of frames to render in benchmark mode */
      else if (tag == "-benchmark") {
        g_skipBenchmarkFrames = cin->getInt();
        g_numBenchmarkFrames  = cin->getInt();
	g_interactive = false;
      }

      /* parse camera parameters */
      else if (tag == "-vp") {
	g_camera.from = cin->getVec3fa();
      }
      else if (tag == "-vi") g_camera.to = cin->getVec3fa();
      else if (tag == "-vd") g_camera.to = g_camera.from + cin->getVec3fa();
      else if (tag == "-vu") g_camera.up = cin->getVec3fa();
      else if (tag == "-fov") g_camera.fov = cin->getFloat();

      /* frame buffer size */
      else if (tag == "-size") {
        g_width = cin->getInt();
        g_height = cin->getInt();
      }

      /* full screen mode */
      else if (tag == "-fullscreen") 
        g_fullscreen = true;
      
      /* rtcore configuration */
      else if (tag == "-rtcore")
        g_rtcore = cin->getString();

      /* number of threads to use */
      else if (tag == "-threads")
        g_numThreads = cin->getInt();

      /* skip unknown command line parameter */
      else {
        std::cerr << "unknown command line parameter: " << tag << " ";
        while (cin->peek() != "" && cin->peek()[0] != '-') std::cerr << cin->getString() << " ";
        std::cerr << std::endl;
      }
    }
  }

  void renderBenchmark(const FileName& fileName)
  {
    resize(g_width,g_height);
    AffineSpace3fa pixel2world = g_camera.pixel2world(g_width,g_height);

    double dt = 0.0f;
    size_t numTotalFrames = g_skipBenchmarkFrames + g_numBenchmarkFrames;
    for (size_t i=0; i<numTotalFrames; i++) 
    {
      double t0 = getSeconds();
      render(0.0f,pixel2world.l.vx,pixel2world.l.vy,pixel2world.l.vz,pixel2world.p);
      double t1 = getSeconds();
      std::cout << "frame [" << i << " / " << numTotalFrames << "] ";
      std::cout << 1.0/(t1-t0) << "fps ";
      if (i < g_skipBenchmarkFrames) std::cout << "(skipped)";
      std::cout << std::endl;
      if (i >= g_skipBenchmarkFrames) dt += t1-t0;
    }
    std::cout << "frame [" << g_skipBenchmarkFrames << " - " << numTotalFrames << "] " << std::flush;
    std::cout << double(g_numBenchmarkFrames)/dt << "fps " << std::endl;
    std::cout << "BENCHMARK_RENDER " << double(g_numBenchmarkFrames)/dt << std::endl;
  }

  void renderToFile(const FileName& fileName)
  {
    resize(g_width,g_height);
    AffineSpace3fa pixel2world = g_camera.pixel2world(g_width,g_height);

    render(0.0f,
	   pixel2world.l.vx,
	   pixel2world.l.vy,
	   pixel2world.l.vz,
	   pixel2world.p);
    
    void* ptr = map();
    Ref<Image> image = new Image4c(g_width, g_height, (Col4c*)ptr);
    storeImage(image, fileName);
    unmap();
    cleanup();
  }

  /* main function in embree namespace */
  int main(int argc, char** argv) 
  {
    g_camera.from = Vec3fa(3.21034f,0.320831f,-0.162478f);
    g_camera.to   = Vec3fa(2.57003f,0.524887f, 0.163145f);

    g_camera.from = Vec3fa(-3,3,3);
    g_camera.to = Vec3fa(0,1,0);
    g_camera.up = Vec3fa(0,1,0);

    /* create stream for parsing */
    Ref<ParseStream> stream = new ParseStream(new CommandLineStream(argc, argv));

    /* parse command line */  
    parseCommandLine(stream, FileName());
    //if (g_numThreads) // FIXME
    //g_rtcore += ",threads=" + std::stringOf(g_numThreads);

    /* initialize task scheduler */
#if !defined(RTCORE_EXPORT_ALL_SYMBOLS)
    TaskScheduler::create(g_numThreads);
#endif

    /* initialize ray tracing core */
    init(g_rtcore.c_str());

    /* load scene */
    if (objFilename.str() != "" && objFilename.str() != "none") {
      loadOBJ(objFilename,AffineSpace3f::translate(-offset),g_obj_scene);
      if (objFilename2.str() != "")
        loadOBJ(objFilename2,AffineSpace3f::translate(-offset_mb),g_obj_scene2);
    }

    /* load hair */
    if (hairFilename.str() != "" && hairFilename.str() != "none") {
      loadHair(hairFilename,g_obj_scene,offset);
      if (hairFilename2.str() != "")
        loadHair(hairFilename2,g_obj_scene2,offset_mb);
    }

    if (!g_obj_scene2.empty()) {
      g_obj_scene.set_motion_blur(g_obj_scene2);
    }

   /* load cy_hair */
    if (cy_hairFilename.str() != "") {
      loadCYHair(cy_hairFilename,g_obj_scene,offset);
    }

    /* generate hair on triangles */
    if (hairy_triangles)
      generateHairOnTriangles(g_obj_scene);

    /* tessellate hair */
    if (tessellate_strips > 0) 
      tessellateHair(g_obj_scene);

    /* if scene is empty, create default scene */
    if (g_obj_scene.meshes.size() + g_obj_scene.hairsets.size() == 0) 
    {
      addHairySphere(g_obj_scene,Vec3fa(0,1.5f,0),1.5f);
      addGroundPlane(g_obj_scene,Vec3fa(-10,0,-10),Vec3fa(-10,0,+10),Vec3fa(+10,0,-10),Vec3fa(+10,0,+10));
    }

    /* send model */
    set_scene(&g_obj_scene);

    /* benchmark mode */
    if (g_numBenchmarkFrames)
      renderBenchmark(outFilename);
    
    /* render to disk */
    if (outFilename.str() != "") 
      renderToFile(outFilename);

    /* interactive mode */
    if (g_interactive) {
      initWindowState(argc,argv,tutorialName, g_width, g_height, g_fullscreen);
      enterWindowRunLoop();
    }
    
    return 0;
  }
}

int main(int argc, char** argv)
{
  try {
    return embree::main(argc, argv);
  }
  catch (const std::exception& e) {
    std::cout << "Error: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Error: unknown exception caught." << std::endl;
    return 1;
  }
}
