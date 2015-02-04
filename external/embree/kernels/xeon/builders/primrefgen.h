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

#include "common/scene.h"
#include "primrefalloc.h"
#include "primrefblock.h"
#include "heuristic_fallback.h"

namespace embree
{
  namespace isa
  {
    /*! Generates a list of triangle build primitives from the scene. */
    class PrimRefListGen
    {
      typedef atomic_set<PrimRefBlockT<PrimRef> > PrimRefList;

    public:      
      static void generate(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>* alloc, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRefList& prims, PrimInfo& pinfo);
      
    private:
      
      /*! standard constructor that schedules the task */
      PrimRefListGen (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>* alloc, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRefList& prims, PrimInfo& pinfo);
            
      /*! parallel task to iterate over the primitives */
      TASK_SET_FUNCTION(PrimRefListGen,task_gen_parallel);
      
    private:
      const Scene* scene;                  //!< input geometry
      GeometryTy ty;                       //!< types of geometry to generate
      size_t numTimeSteps;                 //!< number of timesteps to generate
      PrimRefBlockAlloc<PrimRef>* alloc;   //!< allocator for build primitive blocks
      size_t numPrimitives;                //!< number of generated primitives
      TaskScheduler::Task task;
      PrimRefList& prims_o;                 //!< list of build primitives
      PrimInfo& pinfo_o;                   //!< bounding information of primitives
    };

    /*! Generates a list of triangle build primitives from some geometry. */
    template<typename Ty>
    class PrimRefListGenFromGeometry
    {
      typedef atomic_set<PrimRefBlockT<PrimRef> > PrimRefList;

    public:      
      static void generate(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>* alloc, const Ty* geom, PrimRefList& prims, PrimInfo& pinfo);
      
    private:
      
      /*! standard constructor that schedules the task */
      PrimRefListGenFromGeometry (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>* alloc, const Ty* geom, PrimRefList& prims, PrimInfo& pinfo);
            
      /*! parallel task to iterate over the primitives */
      TASK_SET_FUNCTION(PrimRefListGenFromGeometry,task_gen_parallel);
      
      /* input data */
    private:
      const Ty* geom;                      //!< input geometry
      PrimRefBlockAlloc<PrimRef>* alloc;   //!< allocator for build primitive blocks
      TaskScheduler::Task task;
      PrimRefList& prims_o;                 //!< list of build primitives
      PrimInfo& pinfo_o;                   //!< bounding information of primitives
    };

    /*! Generates an array of triangle build primitives from the scene. */
    class PrimRefArrayGen
    {
    public:   
      static void generate_sequential(size_t threadIndex, size_t threadCount, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRef* prims, PrimInfo& pinfo);
      static void generate_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRef* prims, PrimInfo& pinfo);
      
    private:
      
      /*! standard constructor that schedules the task */
      PrimRefArrayGen (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRef* prims_o, PrimInfo& pinfo_o, bool parallel);
            
      /*! parallel task to iterate over the primitives */
      TASK_FUNCTION(PrimRefArrayGen,task_gen_parallel);
      
      /* input data */
    private:
      const Scene* scene;           //!< input geometry
      GeometryTy ty;                //!< types of geometry to generate
      size_t numTimeSteps;          //!< number of timesteps to generate
      size_t numPrimitives;         //!< number of generated primitives
      PrimRef* prims_o;             //!< list of build primitives
      PrimInfo& pinfo_o;            //!< bounding information of primitives
      size_t* dst;                  //!< write-start offset for each thread
    };

    /*! Generates an array of triangle build primitives from some geometry. */
    template<typename Ty>
    class PrimRefArrayGenFromGeometry
    {
    public:   
      static void generate_sequential(size_t threadIndex, size_t threadCount, const Ty* geom, PrimRef* prims, PrimInfo& pinfo);
      static void generate_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Ty* geom, PrimRef* prims, PrimInfo& pinfo);
      
    private:
      
      /*! standard constructor */
      PrimRefArrayGenFromGeometry (const Ty* geom, PrimRef* prims_o, PrimInfo& pinfo_o);
            
      /*! parallel task to iterate over the primitives */
      TASK_FUNCTION(PrimRefArrayGenFromGeometry,task_gen_parallel);
      
      /* input data */
    private:
      const Ty* geom;               //!< input geometry
      PrimRef* prims_o;             //!< list of build primitives
      PrimInfo& pinfo_o;            //!< bounding information of primitives
      size_t* dst;                  //!< write-start offset for each thread
    };
  }
}
