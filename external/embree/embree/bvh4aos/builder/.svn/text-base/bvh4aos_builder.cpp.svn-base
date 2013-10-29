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

#include "bvh4aos_builder.h"
#include "geometry/triangles.h"
#include "geometry/triangle_mesh.h"
#include "simd/mic_m.h"
#include "simd/mic_f.h"
#include "simd/mic_i.h"

#include "bvh4aos_globals.h"
#include "bvh4aos_builder_common.h"

#define TIMER(x) x
#define PARALLEL_BUILD 1
#define CHECK_TREE 0
#define TOUCH_PAGES 1

#define PRE_ALLOCATION_FACTOR_QBVH_HIGH_QUALITY 0.77f
#define PRE_ALLOCATION_FACTOR_QBVH_FAST         0.86f

namespace embree
{
  enum {
    QBVH_HIGH_QUALITY  = 1,
    QBVH_FAST          = 2,
  };

  MIC_ALIGN unsigned int BVH4AOSBuilder::buildMode = QBVH_HIGH_QUALITY;

  _INLINE size_t align2MB(size_t size)
  {
    if (size == 0) return (1LL << 21);
    const size_t pages = size >> 21;
    const size_t mod   = size % (1LL << 21) ? 1 : 0;
    return (pages+mod) * (1LL << 21);
  }

  _INLINE void touch2MPages(void *addr, size_t size)
  {
    const size_t pageSize = (1LL << 21);
    for (size_t i=0;i<size;i+=pageSize)
      ((char*)addr)[i] = 0;
  }
  
  BVH4AOSBuilder::BVH4AOSBuilder(TaskScheduler::Event* event_i, BVH4AOS* bvh, const size_t minLeafSize, const size_t maxLeafSize,const bool highQuality)
    : bvh(bvh)
  {

    TIMER(Timer timer; unsigned long cycles = 0); 
    TIMER(if (unlikely(g_verbose)) timer.start());

    /*! verify correct input parameters */
    RTCGeometry* geom = bvh->geom;
    mesh = dynamic_cast<TriangleMesh*>(geom);
    if (mesh == NULL    ) throw std::runtime_error("BVH4AOSBuilder: triangle mesh expected");
    if (minLeafSize != 4) throw std::runtime_error("BVH4AOSBuilder: setting minLeafSize not supported.");
    if (maxLeafSize != 4) throw std::runtime_error("BVH4AOSBuilder: setting maxLeafSize not supported.");
    if (bvh->name() != "bvh4aos.triangle1") throw std::runtime_error("BVH4AOSBuilder: only bvh4aos.triangle1 is supported.");

    /*! input triangle mesh */
    RTCTriangle* triangles = mesh->triangles;
    size_t numTriangles    = mesh->numTriangles;
    RTCVertex* vertices    = mesh->vertices;
    size_t numVertices     = mesh->numVertices;


    buildMode = highQuality ? QBVH_HIGH_QUALITY : QBVH_FAST;

    float pre_allocation_factor = 1.0f;

    if (buildMode == QBVH_HIGH_QUALITY)
      pre_allocation_factor = PRE_ALLOCATION_FACTOR_QBVH_HIGH_QUALITY;
    else if (buildMode == QBVH_FAST)      
      pre_allocation_factor = PRE_ALLOCATION_FACTOR_QBVH_FAST;
    else
      FATAL("unknown builder");

    /*! output BVH */
    const size_t pre_alloc_nodes    = (size_t)(pre_allocation_factor*(float)numTriangles);
    const size_t numEmbreeThreads   = TaskScheduler::getNumThreads();
    const size_t bytesAABBArray     = align2MB(numTriangles * sizeof(AABB) + 16 * CACHELINE_SIZE);
    const size_t bytesNodeArray     = align2MB(max(pre_alloc_nodes,(size_t)numEmbreeThreads*2*64) * sizeof(BVHNode) +  2 * CACHELINE_SIZE);
    const size_t bytesTriangleArray = align2MB(sizeof(TriangleAccel)*numTriangles); 

#if 0
    PRINT(numTriangles);
    PRINT(bytesAABBArray);
    PRINT(bytesNodeArray);
    PRINT(bytesTriangleArray);
#endif

    //ParallelQBVHBuilder::max_nodes = max(pre_alloc_nodes,(unsigned int)numEmbreeThreads*2*64);
    bvh->init(bytesNodeArray,bytesTriangleArray);

    ParallelQBVHBuilder::max_nodes = bytesNodeArray / sizeof(AABB);

    BVH4AOS::Node* nodes_o     = (BVH4AOS::Node*) bvh->alloc_nodes->base();
    Triangle1v*    triangles_o = (Triangle1v*   ) bvh->alloc_tris ->base();

    if (unlikely(g_verbose)) 
      std::cout << "Allocating temporary array for AABBs (" << (float)bytesAABBArray / (1024.0f*1024.0f) << " MB)" << std::endl;

    if (ParallelQBVHBuilder::aabb == NULL || ParallelQBVHBuilder::triangles != numTriangles) 
      {
	if (ParallelQBVHBuilder::aabb != NULL)
	  {
	    os_free(ParallelQBVHBuilder::aabb,bytesAABBArray);
	  }
	ParallelQBVHBuilder::aabb = (AABB*)os_reserve(bytesAABBArray);	
      }


    ParallelQBVHBuilder::atomicID.reset();

    ParallelQBVHBuilder::globalBuildTrianglePtr  = (BuilderTriangle*)triangles;
    ParallelQBVHBuilder::globalVertexPtr         = (Vec3fa*)vertices;
    ParallelQBVHBuilder::globalAccelPtr          = (TriangleAccel*)triangles_o;
    ParallelQBVHBuilder::globalNodePtr           = (BVHNode*)nodes_o;
    ParallelQBVHBuilder::triangles    = numTriangles;
    ParallelQBVHBuilder::vertices     = numVertices;
    ParallelQBVHBuilder::tmp_aabb     = (AABB*)triangles_o;    
    ParallelQBVHBuilder::initQBVHNode[0].setToInfPoint();
    ParallelQBVHBuilder::initQBVHNode[0].ext_min.t = (unsigned int)(1 << 31);
    ParallelQBVHBuilder::initQBVHNode[1].setToInfPoint();
    ParallelQBVHBuilder::initQBVHNode[1].ext_min.t = (unsigned int)(1 << 31);
    ParallelQBVHBuilder::CONTROL_THREAD_ID  = 0;	  
    ParallelQBVHBuilder::numActiveLocalWorkQueues = 1; 
    

#if TOUCH_PAGES == 1
    touch2MPages(ParallelQBVHBuilder::globalNodePtr,bytesNodeArray);
    touch2MPages(ParallelQBVHBuilder::globalAccelPtr,bytesTriangleArray);
    touch2MPages(ParallelQBVHBuilder::aabb,bytesAABBArray);
#endif

#if PARALLEL_BUILD == 1
    bool parallelBuild = true;
#else
    bool parallelBuild = false;
#endif
    const unsigned int minTrisForParallelBuild = MAX_MIC_THREADS * 64;


    TIMER(
	  if (unlikely(g_verbose)) 
	    {
	      cycles = (unsigned long)timer.stop();        
	      std::cout << "==> Allocate data arrays for QBVH builder " << cycles << " cycles " << cycles/1000000.0f << " mcycles" << std::endl;
	    }
	  );

    if (parallelBuild && numTriangles >= minTrisForParallelBuild)
      {
	size_t numThreads = TaskScheduler::getNumThreads();
	ParallelQBVHBuilder::NUM_TOTAL_THREADS        = numThreads; 
	ParallelQBVHBuilder::NUM_TOTAL_CORES          = numThreads/4;
	ParallelQBVHBuilder::numActiveLocalWorkQueues = numThreads/4;
        new (&task) TaskScheduler::Task(event_i,_task_build_parallel,this,numThreads,_finish,this,"BVH4AOS::Build");
	TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
      }
    else
      {
	// =====================
	// ==== SERIAL PATH ====
	// =====================

	TIMER(Timer timer; unsigned long cycles = 0); 
	TIMER(timer.start());

	if (unlikely(g_verbose && numTriangles < minTrisForParallelBuild))
	  std::cout << "WARNING: not enough triangles for parallel build, falling back to serial mode" << std::endl << std::flush;
	    
	ParallelQBVHBuilder::NUM_TOTAL_THREADS  = 1; 
	ParallelQBVHBuilder::NUM_TOTAL_CORES    = 1;
   
	if (buildMode == QBVH_HIGH_QUALITY)
	  {
	    if (unlikely(g_verbose)) std::cout << "SAH QBVH Builder" << std::endl;

	    ParallelQBVHBuilder::numActiveLocalWorkQueues = 1; 
	    ParallelQBVHBuilder::build();    
	  }
	else
	  FATAL("unknown builder");

#if CHECK_TREE == 1
	if (buildMode == QBVH_HIGH_QUALITY)
	  ParallelQBVHBuilder::check_tree( true );
	else
	  FATAL("unknown builder");
#endif

	ParallelQBVHBuilder::convertToOptimizedQBVH();

        /* set root node and bounds */
        BVH4AOS::Node* root = (BVH4AOS::Node*) bvh->nodePtr();
        bvh->root = root->child(0);
        bvh->bounds = root->bounds(0);


        if (unlikely(g_verbose)) 
        {
          TIMER(
                cycles = (unsigned long)timer.stop();    
                const double mtris_per_sec = (double)(Timer::freqInMHz*1E6) / (float)cycles * (float)ParallelQBVHBuilder::triangles / 1000000.0f;
                std::cout << "==> Build QBVH for " << ParallelQBVHBuilder::triangles << " triangles in " << cycles/1000000.0f << " Mcycles => "  << mtris_per_sec << " mtris/s (@ " << Timer::freqInMHz << " MHz)" << std::endl;
                std::cout << "build time = " << (double)cycles / (double)(Timer::freqInMHz*1E3) << " ms"  << std::endl;
                );
          
          
          const unsigned int used_nodes = ParallelQBVHBuilder::getQBVHNodes() << 2;
          const unsigned int max_nodes  = ParallelQBVHBuilder::getMaxAllocatedNodes();
          std::cout << "max_nodes " << max_nodes << " <-> used nodes " << used_nodes << " -> " << (float)used_nodes * 100.0f / max_nodes << " %" << std::endl;
          
          
          unsigned long accelSizeInBytes = 0;
          accelSizeInBytes += ParallelQBVHBuilder::getQBVHNodes() * 4 * sizeof(BVHNode) + CACHELINE_SIZE;
          accelSizeInBytes += sizeof(TriangleAccel)*ParallelQBVHBuilder::triangles;
          std::cout << "accel size = " << (float)accelSizeInBytes / 1E6 << " MB"  << std::endl;
        }

	if (buildMode == QBVH_HIGH_QUALITY /* || buildMode == QBVH_FAST */)
	  {
	    os_free(ParallelQBVHBuilder::aabb,bytesAABBArray); 
	    ParallelQBVHBuilder::aabb = NULL;
	  }

      }
  }


  void BVH4AOSBuilder::task_build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* group)
  {
    size_t numThreads = TaskScheduler::getNumThreads();

    if (threadIndex == 0)
      {
	TIMER(Timer timer; unsigned long cycles = 0); 
	TIMER(timer.start());


	if (buildMode == QBVH_HIGH_QUALITY)
	  {
	    if (unlikely(g_verbose)) 
	      std::cout << "SAH QBVH Builder" << std::endl;
	    ParallelQBVHBuilder::build();    
	  }
	else if (buildMode == QBVH_FAST)      
	  {
	    if (unlikely(g_verbose)) std::cout << "Morton-Code QBVH Builder" << std::endl;

	    ParallelQBVHBuilder::numActiveLocalWorkQueues =  MAX_MIC_THREADS / NUM_LOCAL_ATOMIC_BUILD_RECORDS;
	    ParallelQBVHBuilderMortonCode::build();
	  }
	else
	  FATAL("unknown builder");


#if CHECK_TREE == 1
	if (buildMode == QBVH_HIGH_QUALITY)
	  ParallelQBVHBuilder::check_tree( true );
	else if (buildMode == QBVH_FAST)      
	  ParallelQBVHBuilderMortonCode::check_tree();
	else
	  FATAL("unknown builder");
#endif

	ParallelQBVHBuilder::convertToOptimizedQBVH();


        if (unlikely(g_verbose)) 
        {
          TIMER(
                cycles = (unsigned long)timer.stop();    
                const double mtris_per_sec = (double)(Timer::freqInMHz*1E6) / (float)cycles * (float)ParallelQBVHBuilder::triangles / 1000000.0f;

                std::cout << "==> Build QBVH for " << ParallelQBVHBuilder::triangles << " triangles in " << cycles/1000000.0f << " Mcycles => "  << mtris_per_sec << " mtris/s (@ " << Timer::freqInMHz << " MHz)" << std::endl;
                std::cout << "build time = " <<  (float)cycles / (Timer::freqInMHz*1E3) << " ms"  << std::endl;
                );
          
          
          const unsigned int used_nodes = ParallelQBVHBuilder::getQBVHNodes() << 2;
          const unsigned int max_nodes  = ParallelQBVHBuilder::getMaxAllocatedNodes();
          std::cout << "max_nodes " << max_nodes << " <-> used nodes " << used_nodes << " -> " << (float)used_nodes * 100.0f / max_nodes << " %" << std::endl;
          
          
          unsigned long accelSizeInBytes = 0;
          accelSizeInBytes += ParallelQBVHBuilder::getQBVHNodes() * 4 * sizeof(BVHNode) + CACHELINE_SIZE;
          accelSizeInBytes += sizeof(TriangleAccel)*ParallelQBVHBuilder::triangles;
          std::cout << "accel size = " << (float)accelSizeInBytes / 1E6 << " MB"  << std::endl;
        }

        size_t bytesAABBArray = align2MB(mesh->numTriangles * sizeof(AABB) + 16 * CACHELINE_SIZE);

	if (buildMode == QBVH_HIGH_QUALITY /* || buildMode == QBVH_FAST */)
	  {
	    os_free(ParallelQBVHBuilder::aabb,bytesAABBArray); 
	    ParallelQBVHBuilder::aabb = NULL;
	  }

	QBVHTaskScheduler::releaseThreads();

      }
    else
      {
	QBVHTaskScheduler::taskDispatchMainLoop(threadIndex);
      }

  }

  void BVH4AOSBuilder::finish(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event) 
  {
    /* set root node and bounds */
    BVH4AOS::Node* root = (BVH4AOS::Node*) bvh->nodePtr();
    bvh->root = root->child(0);
    bvh->bounds = root->bounds(0);

    /* delete has to be called in finish function */
    delete this;
  }
}

