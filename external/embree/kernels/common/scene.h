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

#include "common/default.h"

#include "scene_triangle_mesh.h"
#include "scene_user_geometry.h"
#include "scene_bezier_curves.h"
#include "scene_subdiv_mesh.h"

#include "common/acceln.h"
#include "geometry.h"

namespace embree
{
  /*! Base class all scenes are derived from */
  class Scene : public Accel
  {
    ALIGNED_CLASS;

  public:
    template<typename Ty> // FIXME: give motion blur meshes different type to iterate over them here
    class Iterator
    {
    public:
      Iterator ()  {}

      Iterator (Scene* scene) 
        : scene(scene) {}

      __forceinline Ty* operator[] (const size_t i) 
      {
        Geometry* geom = scene->geometries[i];
        if (geom == NULL) return NULL;
        if (!geom->isEnabled()) return NULL;
        if (geom->type != Ty::geom_type) return NULL;
        return (Ty*) geom;
      }

      __forceinline size_t size() const {
        return scene->size();
      }
      
    private:
      Scene* scene;
    };

  public:
    
    /*! Scene construction */
    Scene (RTCSceneFlags flags, RTCAlgorithmFlags aflags);

    void createTriangleAccel();
    void createHairAccel();
    void createSubdivAccel();

    /*! Scene destruction */
    ~Scene ();

    /*! Creates new user geometry. */
    unsigned int newUserGeometry (size_t items);

    /*! Creates a new scene instance. */
    unsigned int newInstance (Scene* scene);

    /*! Creates a new triangle mesh. */
    unsigned int newTriangleMesh (RTCGeometryFlags flags, size_t maxTriangles, size_t maxVertices, size_t numTimeSteps);

    /*! Creates a new collection of quadratic bezier curves. */
    unsigned int newBezierCurves (RTCGeometryFlags flags, size_t maxCurves, size_t maxVertices, size_t numTimeSteps);

    /*! Creates a new subdivision mesh. */
    unsigned int newSubdivisionMesh (RTCGeometryFlags flags, size_t numFaces, size_t numEdges, size_t numVertices, size_t numEdgeCreases, size_t numVertexCreases, size_t numHoles, size_t numTimeSteps);

    /*! Builds acceleration structure for the scene. */
    void build (size_t threadIndex, size_t threadCount);

    /*! stores scene into binary file */
    void write(std::ofstream& file);

    /*! build task */
    TASK_RUN_FUNCTION(Scene,task_build_parallel);
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
    __forceinline SubdivMesh* getSubdivMesh(size_t i) { 
      assert(i < geometries.size()); 
      assert(geometries[i]);
      assert(geometries[i]->type == SUBDIV_MESH);
      return (SubdivMesh*) geometries[i]; 
    }
    __forceinline const SubdivMesh* getSubdivMesh(size_t i) const { 
      assert(i < geometries.size()); 
      assert(geometries[i]);
      assert(geometries[i]->type == SUBDIV_MESH);
      return (SubdivMesh*) geometries[i]; 
    }
    __forceinline UserGeometryBase* getUserGeometrySafe(size_t i) { 
      assert(i < geometries.size()); 
      if (geometries[i] == NULL) return NULL;
      if (geometries[i]->type != USER_GEOMETRY) return NULL;
      else return (UserGeometryBase*) geometries[i]; 
    }
    __forceinline BezierCurves* getBezierCurves(size_t i) { 
      assert(i < geometries.size()); 
      assert(geometries[i]);
      assert(geometries[i]->type == BEZIER_CURVES);
      return (BezierCurves*) geometries[i]; 
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

  public:
    std::vector<int> usedIDs;
    std::vector<Geometry*> geometries; //!< list of all user geometries
    
  public:
    AccelN accels;
    unsigned int commitCounter;
    atomic_t numMappedBuffers;         //!< number of mapped buffers
    RTCSceneFlags flags;
    RTCAlgorithmFlags aflags;
    bool needTriangles; 
    bool needVertices; // FIXME: this flag is also used for hair geometry, but there should be a second flag
    bool is_build;
    MutexSys mutex;
    AtomicMutex geometriesMutex;
    
    /*! global lock step task scheduler */
    __aligned(64) LockStepTaskScheduler lockstep_scheduler;

  public:
    atomic_t numTriangles;             //!< number of enabled triangles
    atomic_t numTriangles2;            //!< number of enabled motion blur triangles
    atomic_t numBezierCurves;          //!< number of enabled curves
    atomic_t numBezierCurves2;         //!< number of enabled motion blur curves
    atomic_t numSubdivPatches;         //!< number of enabled subdivision patches
    atomic_t numSubdivPatches2;        //!< number of enabled motion blur subdivision patches
    atomic_t numUserGeometries1;       //!< number of enabled user geometries

    atomic_t numIntersectionFilters4;   //!< number of enabled intersection/occlusion filters for 4-wide ray packets
    atomic_t numIntersectionFilters8;   //!< number of enabled intersection/occlusion filters for 8-wide ray packets
    atomic_t numIntersectionFilters16;  //!< number of enabled intersection/occlusion filters for 16-wide ray packets
  };
}
