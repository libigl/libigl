// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
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

#ifndef __EMBREE_SCENE_H__
#define __EMBREE_SCENE_H__

#include "common/default.h"

#include "scene_triangle_mesh.h"
#include "scene_user_geometry.h"
#include "scene_quadratic_bezier_curves.h"

#include "common/acceln.h"
#include "geometry.h"
#include "common/buildsource.h"

namespace embree
{
  /*! Base class all scenes are derived from */
  class Scene : public Accel
  {
    ALIGNED_CLASS;
  public:

    typedef TriangleMeshScene::TriangleMesh TriangleMesh;
    
    /*! Scene construction */
    Scene (RTCSceneFlags flags, RTCAlgorithmFlags aflags);

    /*! Scene destruction */
    ~Scene ();

    /*! Creates new user geometry. */
    unsigned int newUserGeometry (size_t items);

    /*! Creates a new scene instance. */
    unsigned int newInstance (Scene* scene);

    /*! Creates a new triangle mesh. */
    unsigned int newTriangleMesh (RTCGeometryFlags flags, size_t maxTriangles, size_t maxVertices, size_t numTimeSteps);

    /*! Creates a new collection of quadratic bezier curves. */
    unsigned int newQuadraticBezierCurves (RTCGeometryFlags flags, size_t maxCurves, size_t maxVertices, size_t numTimeSteps);

    /*! Builds acceleration structure for the scene. */
    void build ();

    void build (size_t threadIndex, size_t threadCount);

    /*! build task */
    TASK_COMPLETE_FUNCTION(Scene,task_build);
    TaskScheduler::Task task;

    /* return number of geometries */
    __forceinline size_t size() const { return geometries.size(); }
    
    /* add user geometry to scene */
    unsigned int add (Geometry* geometry);
    
    /* removes user geometry from scene again */
    void remove(Geometry* geometry);

    /* determines of the scene is ready to get build */
    bool ready() { return numMappedBuffers == 0; }

    /* get mesh by ID */
    __forceinline       Geometry* get(size_t i)       { assert(i < geometries.size()); return geometries[i]; }
    __forceinline const Geometry* get(size_t i) const { assert(i < geometries.size()); return geometries[i]; }
    __forceinline       Geometry* get_locked(size_t i)  { 
      Lock<AtomicMutex> lock(geometriesMutex);
      Geometry *g = geometries[i]; 
      assert(i < geometries.size()); 
      return g; 
    }

    /* get triangle mesh by ID */
    __forceinline TriangleMesh* getTriangleMesh(size_t i) { 
      assert(i < geometries.size()); 
      assert(geometries[i]);
      assert(geometries[i]->type == TRIANGLE_MESH);
      return (TriangleMesh*) geometries[i]; 
    }
    __forceinline const TriangleMesh* getTriangleMesh(size_t i) const { 
      assert(i < geometries.size()); 
      assert(geometries[i]);
      assert(geometries[i]->type == TRIANGLE_MESH);
      return (TriangleMesh*) geometries[i]; 
    }
    __forceinline TriangleMesh* getTriangleMeshSafe(size_t i) { 
      assert(i < geometries.size()); 
      if (geometries[i] == NULL) return NULL;
      if (geometries[i]->type != TRIANGLE_MESH) return NULL;
      else return (TriangleMesh*) geometries[i]; 
    }
    __forceinline UserGeometryScene::Base* getUserGeometrySafe(size_t i) { 
      assert(i < geometries.size()); 
      if (geometries[i] == NULL) return NULL;
      if (geometries[i]->type != USER_GEOMETRY && geometries[i]->type != INSTANCES) return NULL;
      else return (UserGeometryScene::Base*) geometries[i]; 
    }


    /* test if this is a static scene */
    __forceinline bool isStatic() const { return embree::isStatic(flags); }

    /* test if this is a dynamic scene */
    __forceinline bool isDynamic() const { return embree::isDynamic(flags); }

    __forceinline bool isCompact() const { return embree::isCompact(flags); }
    __forceinline bool isCoherent() const { return embree::isCoherent(flags); }
    __forceinline bool isRobust() const { return embree::isRobust(flags); }
    __forceinline bool isHighQuality() const { return embree::isHighQuality(flags); }

    /* test if scene got already build */
    __forceinline bool isBuild() const { return is_build; }

    struct FlatTriangleAccelBuildSource : public BuildSource
    {
      FlatTriangleAccelBuildSource (Scene* scene, size_t numTimeSteps = 1)
        : scene(scene), numTimeSteps(numTimeSteps) {}

      bool isEmpty () const { 
        if (numTimeSteps == 1) return scene->numTriangleMeshes  == 0;
        else                   return scene->numTriangleMeshes2 == 0;
      }
      
      size_t groups () const { 
        return scene->geometries.size();
      }
      
      size_t prims (size_t group, size_t* numVertices) const 
      {
        if (scene->get(group) == NULL || scene->get(group)->type != TRIANGLE_MESH) return 0;
        TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(group);
        if (mesh == NULL || !mesh->isEnabled() || mesh->numTimeSteps != numTimeSteps) return 0;
        if (numVertices) *numVertices = mesh->numVertices;
        return mesh->numTriangles;
      }

      const BBox3f bounds(size_t group, size_t prim) const 
      {
	assert(scene->get(group) != NULL);
	assert(scene->get(group)->type == TRIANGLE_MESH);

        TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(group);
        if (mesh == NULL) return empty;
        return mesh->bounds(prim);
      }

      void bounds(size_t group, size_t begin, size_t end, BBox3f* bounds_o) const 
      {
	assert(scene->get(group) != NULL);
	assert(scene->get(group)->type == TRIANGLE_MESH);

        TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(group);
        if (mesh == NULL) { 
          for (size_t i=0; i<end-begin; i++)
            bounds_o[i] = empty;
        } else {
          for (size_t i=begin; i<end; i++)
            bounds_o[i-begin] = mesh->bounds(i);
        }
      }

      const Vec3fa vertex(size_t group, size_t prim, size_t vtxID) const 
      {
	assert(scene->get(group) != NULL);
	assert(scene->get(group)->type == TRIANGLE_MESH);

        const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(group);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(prim);
	return mesh->vertex(tri.v[vtxID]);
      }
      
      void split (const PrimRef& prim, int dim, float pos, PrimRef& left_o, PrimRef& right_o) const 
      {
	assert(scene->get(prim.geomID()));
	assert(scene->get(prim.geomID())->type == TRIANGLE_MESH);

        const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(prim.geomID());
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(prim.primID());
        const Vec3fa& v0 = mesh->vertex(tri.v[0]);
        const Vec3fa& v1 = mesh->vertex(tri.v[1]);
        const Vec3fa& v2 = mesh->vertex(tri.v[2]);
        splitTriangle(prim,dim,pos,v0,v1,v2,left_o,right_o);
      }
      
    public:
      Scene* scene;
      size_t numTimeSteps;
    };

    
  public:
    std::vector<int> usedIDs;
    std::vector<Geometry*> geometries; //!< list of all user geometries
    
  public:
    AccelN accels;
    atomic_t numMappedBuffers;         //!< number of mapped buffers
    RTCSceneFlags flags;
    RTCAlgorithmFlags aflags;
    bool needTriangles;
    bool needVertices;
    bool is_build;
    MutexSys mutex;
    AtomicMutex geometriesMutex;

  public:
    atomic_t numTriangleMeshes;        //!< number of enabled triangle meshes
    atomic_t numTriangleMeshes2;       //!< number of enabled motion blur triangle meshes
    atomic_t numCurveSets;             //!< number of enabled curve sets
    atomic_t numCurveSets2;            //!< number of enabled motion blur curve sets
    atomic_t numUserGeometries;        //!< number of enabled user geometries
    
  public:
    FlatTriangleAccelBuildSource flat_triangle_source_1;
    FlatTriangleAccelBuildSource flat_triangle_source_2;
  };

  typedef Builder* (*TriangleMeshBuilderFunc)(void* accel, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize);
  typedef Builder* (*BuilderFunc)            (void* accel, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize);

}

#endif

