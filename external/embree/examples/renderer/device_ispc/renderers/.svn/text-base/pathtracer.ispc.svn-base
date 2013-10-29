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

#define PDF_CULLING 0.0f 

#include "renderer.isph"
#include "materials/medium.isph"
#include "materials/material.isph"
#include "lights/light.isph"
#include "brdfs/compositedbrdf.isph"
#include "textures/image.isph"
#include "samplers/precomputed_sampler.isph"
#include "framebuffers/framebuffer.isph"

#if defined(__MIC__) 
#  define PACKET_WIDTH 4
#  define PACKET_HEIGHT 4
#elif defined (__AVX__)
#  define PACKET_WIDTH 4
#  define PACKET_HEIGHT 2
#elif defined (__SSE4__)
#  define PACKET_WIDTH 2
#  define PACKET_HEIGHT 2
#else
#  error "unsupported configuration"
#endif

#define TILE_SIZE_X 8
#define TILE_SIZE_Y 8

//////////////////////////////////////////////////////////////////
// LightPath

struct LightPath 
{
  Ray    ray;                /*! Last ray in the path. */
  Medium lastMedium;             /*! Medium the last ray travels inside. */
  uniform uint depth;            /*! Recursion depth of path. */
  vec3f  throughput;             /*! Determines the fraction of
                                     radiance reaches the pixel along
                                     the path. */
  bool   ignoreVisibleLights;    /*! If the previous shade point used
                                     shadow rays we have to ignore the
                                     emission of geometrical lights to
                                     not double count them. */
  bool   unbent;                    /*! True of the ray path is a straight line. */
};

inline void init_LightPath(LightPath& lp, const Ray &ray)
{
  lp.ray = ray;
  lp.lastMedium = make_Medium_Vacuum();
  lp.depth = 0;
  lp.throughput = make_vec3f(1.f);
  lp.ignoreVisibleLights = false;
  lp.unbent = true;
}

inline void extend_fast(LightPath& lp,
                   const vec3f nextRay_org, 
                   const vec3f nextRay_dir, 
                   const float nextRay_near, 
                   const float nextRay_far, 
                   const vec3f weight, 
                   const bool ignoreVL)
{
  lp.unbent = lp.unbent & eq(nextRay_dir,lp.ray.dir);
  init_Ray(lp.ray,nextRay_org,nextRay_dir,nextRay_near,nextRay_far);


  lp.depth = lp.depth+1;
  lp.throughput = mul(lp.throughput,weight);
  lp.ignoreVisibleLights = ignoreVL;
}

//////////////////////////////////////////////////////////////////
// PathTracer

#define MAX_LIGHTS 1000

struct PathTracer 
{
  uniform Renderer base;
  uniform int32 numRays;
  uniform int maxDepth;
  uniform float minContribution;
  uniform float epsilon;
  uniform int spp;
  uniform int iteration;
  uniform Image* uniform backplate;

  /*! Random variables. */
  uniform int lightSampleID;            //!< 2D random variable to sample the light source.
  uniform int firstScatterSampleID;     //!< 2D random variable to sample the BRDF.
  uniform int firstScatterTypeSampleID; //!< 1D random variable to sample the BRDF type to choose.
  uniform int precomputedLightSampleID[MAX_LIGHTS]; //!< ID of precomputed light samples for lights that need precomputations.
  uniform PrecomputedSampler sampler;
};

vec3f PathTraceIntegrator_Li(const uniform PathTracer* uniform this,
                             const vec2f &pixel,
                             LightPath &lightPath, 
                             const uniform Scene *uniform scene,
                             const uniform PrecomputedSample* uniform sample_,
                             uint &numRays)
{
  //uniform uint/*BRDFType*/ directLightingBRDFTypes = (uniform uint)(DIFFUSE|GLOSSY);
  uniform uint/*BRDFType*/ directLightingBRDFTypes = (uniform uint)(DIFFUSE);
  uniform uint/*BRDFType*/ giBRDFTypes = (uniform uint)(ALL);
  
  vec3f L = make_vec3f(0.f);
  vec3f Lw = make_vec3f(1.f);

  while (true) 
  {
    /*! Terminate path if too long or contribution too low. */
    if (lightPath.depth >= this->maxDepth)
      return L;
	if (reduce_max(lightPath.throughput) <= this->minContribution)
	  return L;

    /*! Traverse ray. */
    DifferentialGeometry dg;
    scene->intersector->intersect(scene->intersector,lightPath.ray); 
    postIntersect(scene,lightPath.ray,dg);
    numRays++;  

    const vec3f wo = neg(lightPath.ray.dir);

    /*! Environment shading when nothing hit. */
    if (noHit(lightPath.ray)) 
    {
      if ((bool)this->backplate & lightPath.unbent) {
        const int x = clamp((int)(pixel.x * this->backplate->size.x), 0, (int)this->backplate->size.x-1);
        const int y = clamp((int)(pixel.y * this->backplate->size.y), 0, (int)this->backplate->size.y-1);
        L = this->backplate->get_nearest_varying(this->backplate,x,y);
      }
      else {
        if (!lightPath.ignoreVisibleLights) {
          for (uniform int i=0; i<scene->num_envLights; i++) {
            uniform const EnvironmentLight *uniform l = scene->envLights[i]; 
            L = add(L, mul(Lw,l->Le(l,wo)));
          }      
        }
      }
      return L;
    }

    /*! Shade surface. */
    uniform CompositedBRDF brdfs;
    CompositedBRDF__Constructor(&brdfs);
#if 0
    uniform Material* material = scene->geometry[lightPath.ray.id0]->material;
    foreach_unique(m in material) 
      if (m != NULL) m->shade(m,lightPath.ray, lightPath.lastMedium, dg, brdfs);
#else
    foreach_unique(geomID in lightPath.ray.id0) {
      uniform Material* uniform m = scene->geometry[geomID]->material;
      if (m != NULL) m->shade(m,lightPath.ray, lightPath.lastMedium, dg, brdfs);
    }
#endif

    /*! Add light emitted by hit area light source. */
    if (!lightPath.ignoreVisibleLights) {
      foreach_unique(geomID in lightPath.ray.id0) {
        const uniform AreaLight* uniform l = scene->geometry[geomID]->light;
        if (l != NULL) L = add(L, mul(Lw, l->Le(l,dg,wo)));
      }
    }

    /*! Check if any BRDF component uses direct lighting. */
    bool useDirectLighting = brdfs.brdfTypes & directLightingBRDFTypes;

    /*! Direct lighting. Shoot shadow rays to all light sources. */
    if (useDirectLighting) 
    {
      uniform int numAllLights = min(MAX_LIGHTS,scene->num_allLights);
      for (uniform int i=0; i<numAllLights; i++) 
      {
        uniform Light* uniform light = scene->allLights[i];

        /*! Either use precomputed samples for the light or sample light now. */
        LightSample ls; 
	ls.wi.v = make_vec3f(0.0f,0.0f,0.0f); ls.wi.pdf = 0.0f;
        if (light->type & TY_PRECOMPUTE_LIGHT_SAMPLES) {
          ls = PrecomputedSample__getLightSample(sample_,this->precomputedLightSampleID[i]);
        }
        else {
          ls.L = light->sample(light, dg, ls.wi, ls.tMax, PrecomputedSample__getVec2f(sample_,this->lightSampleID));
        }

        /*! Ignore zero radiance or illumination from the back. */
        //if (reduce_max(ls.L) <= 0.0f | ls.wi.pdf <= PDF_CULLING | dot(dg.Ns,ls.wi.v) <= 1e-8f) 
        if (reduce_max(ls.L) <= 0.0f | ls.wi.pdf <= PDF_CULLING) 
          continue;

        /*! Evaluate BRDF */
        vec3f brdf = CompositedBRDF__eval(&brdfs,wo,dg,ls.wi.v,directLightingBRDFTypes);
        if (reduce_max(brdf) <= 0.0f)
          continue;
        
        /*! Test for shadows. */
        numRays++;
        Ray shadow_ray; 
        init_Ray(shadow_ray,dg.P,ls.wi.v,dg.error*this->epsilon,ls.tMax-dg.error*this->epsilon);
        if (scene->intersector->occluded(scene->intersector,shadow_ray)) continue; 

        L = add(L,mul(mul(Lw,ls.L),mul(brdf,rcp(ls.wi.pdf))));
      }
    }

    /*! Global illumination. Pick one BRDF component and sample it. */
    if (lightPath.depth >= this->maxDepth) 
      return L;
    
    /*! sample brdf */
    Sample3f wi = make_Sample3f(make_vec3f(0.0f),0.0f); uint type = 0;
    vec2f s  = PrecomputedSample__getVec2f(sample_,this->firstScatterSampleID     + lightPath.depth);
    float ss = PrecomputedSample__getFloat(sample_,this->firstScatterTypeSampleID + lightPath.depth);
    vec3f c = CompositedBRDF__sample(&brdfs,wo,dg,wi,type,s,ss,giBRDFTypes);
    
    /*! Continue only if we hit something valid. */
    if (reduce_max(c) <= 0.0f | wi.pdf <= PDF_CULLING) 
      return L;

    /*! Compute  simple volumetric effect. */
    const vec3f transmission = lightPath.lastMedium.transmission;
    if (ne(transmission,make_vec3f(1.f)))
      c = mul(c, pow(transmission,lightPath.ray.tfar));
    
    /*! Tracking medium if we hit a medium interface. */
    if (type & TRANSMISSION) {
#if 0
      foreach_unique(m in material) 
        if (m != NULL)  m->selectNextMedium(m,lightPath.lastMedium);
#else
      foreach_unique(geomID in lightPath.ray.id0) {
        uniform Material* uniform m = scene->geometry[geomID]->material;
        if (m != NULL)  m->selectNextMedium(m,lightPath.lastMedium);
      }
#endif
    }
    
    /*! Continue the path. */
    extend_fast(lightPath, 
           dg.P,wi.v,dg.error*this->epsilon,inf,
           c,(type & directLightingBRDFTypes) != NONE);

    Lw = mul(mul(Lw,c),rcp(wi.pdf));
  }
  return L;
}

inline vec3f PathTracer__renderPixel(const uniform PathTracer* uniform this,
                                     const uniform Camera *uniform camera,
                                     const uniform Scene  *uniform scene,
                                     const uniform FrameBuffer *uniform fb,
                                     uniform Random& rnd,
                                     const uint ix, const uint iy, 
                                     uint &numRays)
{
  vec3f L = make_vec3f(0.f);
  uniform int set = Random__getInt(&rnd);
  for (uniform int s=0; s<this->spp; s++) 
  {
    uniform PrecomputedSample* uniform sample = 
      PrecomputedSampler__get(&this->sampler,set,this->iteration*this->spp+s);
    
    const vec2f pixelSample = PrecomputedSample__getPixel(sample);
    const vec2f lensSample  = PrecomputedSample__getLens(sample);
    const vec2f screenSample = mul(add(make_vec2f(ix,iy),pixelSample),fb->invSize);
    
    Ray ray;
    camera->initRay(camera,ray,screenSample,lensSample);
    LightPath lightPath; init_LightPath(lightPath,ray);
    L = add(L, PathTraceIntegrator_Li(this,screenSample,lightPath,scene,sample,numRays));
  }
  return mul(L, rcp(this->spp));
} 

task void PathTracer__renderTile(uniform PathTracer* uniform this,
                                 const uniform Camera *uniform camera,
                                 const uniform Scene  *uniform scene,
                                 const uniform ToneMapper* uniform toneMapper,
                                 uniform FrameBuffer *uniform fb,
                                 uniform AccuBuffer *uniform accu,
                                 const uniform int accuMode,
                                 const uniform uint numTiles_x) 
{
  uint numRays = 0;
  const uniform uint tile_y = taskIndex / numTiles_x;
  const uniform uint tile_x = taskIndex - tile_y * numTiles_x;
  const uint sample_y = programIndex / PACKET_WIDTH; 
  const uint sample_x = programIndex - sample_y * PACKET_WIDTH;

  uniform Random rnd; 
  uniform int uniqueID = tile_x * 917 + tile_y * 81551 + 3433*g_serverID;
  Random__setSeed(&rnd,uniqueID); // expensive
  
  const uniform uint tile_y0 = tile_y * TILE_SIZE_Y;
  const uniform uint tile_x0 = tile_x * TILE_SIZE_X;

  for (uniform uint iy=0; iy<TILE_SIZE_Y; iy+=PACKET_HEIGHT)
  {
    const uint y = (tile_y0 + iy) + sample_y;
    if (y >= fb->size.y) continue;
    
    if (!activeLine(y)) continue;
    size_t _y = raster2buffer(y);

    for (uniform unsigned int ix=0; ix<TILE_SIZE_X; ix+=PACKET_WIDTH) 
    { 
      const uint x = (tile_x0 + ix) + sample_x;
      if (x >= fb->size.x) continue;

      vec3f R = PathTracer__renderPixel(this,camera,scene,fb,rnd,x,y,numRays);
      vec3f d = AccuBuffer__update(accu,x,_y,R,accuMode);
      if (toneMapper) d = toneMapper->toneMap(toneMapper,d,x,y,fb->size);
      fb->set(fb,x,_y,d);
    }
  }

  /* count number of rays */
  uniform int num = 0;
  foreach_active(i) {
    num += extract(numRays,i);
  }
  atomic_add_global(&this->numRays,num);
}

void PathTracer__initSampler(uniform PathTracer* uniform this, const uniform Scene* uniform scene)
{
  PrecomputedSampler__reset(&this->sampler);
  PrecomputedSampler__Constructor(&this->sampler,16*16,64); // first argument HAS to be a square number, both have to be a power of two

  this->lightSampleID = PrecomputedSampler__request2D(&this->sampler,1);
  uniform int numAllLights = min(MAX_LIGHTS,scene->num_allLights);

  for (uniform int i=0; i<numAllLights; i++) {
    this->precomputedLightSampleID[i] = -1;
    if (scene->allLights[i]->type & TY_PRECOMPUTE_LIGHT_SAMPLES) 
      this->precomputedLightSampleID[i] = PrecomputedSampler__requestLightSample(&this->sampler,this->lightSampleID,scene->allLights[i]);
  }

  this->firstScatterSampleID = PrecomputedSampler__request2D(&this->sampler,this->maxDepth);
  this->firstScatterTypeSampleID = PrecomputedSampler__request1D(&this->sampler,this->maxDepth);
  PrecomputedSampler__init(&this->sampler);
}
 
void PathTracer_renderFrameInit(uniform Renderer* uniform _this,
                                const uniform Scene* uniform scene)
{
  uniform PathTracer* uniform this = (uniform PathTracer* uniform) _this;
  if (PrecomputedSampler__empty(&this->sampler)) PathTracer__initSampler(this,scene);
}

uniform int PathTracer_renderFrame(uniform Renderer* uniform _this,
                                   const uniform Camera* uniform camera,
                                   const uniform Scene* uniform scene,
                                   const uniform ToneMapper* uniform toneMapper,
                                   uniform SwapChain *uniform swapchain,
                                   const uniform int accuMode)
{
  uniform PathTracer* uniform this = (uniform PathTracer* uniform) _this;
  this->numRays = 0;
  if (accuMode == 0) this->iteration = 0;
  uniform int numTiles_x = (swapchain->width+ (TILE_SIZE_X-1))/TILE_SIZE_X;
  uniform int numTiles_y = (swapchain->height+(TILE_SIZE_Y-1))/TILE_SIZE_Y;
  uniform int numTiles = numTiles_x * numTiles_y;
  uniform FrameBuffer* uniform fb = SwapChain__get_buffer(swapchain);
  uniform AccuBuffer* uniform accu = SwapChain__get_accu(swapchain);
  launch[numTiles] PathTracer__renderTile(this,camera,scene,toneMapper,fb,accu,accuMode,numTiles_x);
  sync;
  rtcDebug();
  this->iteration++;
  return this->numRays;
}

void PathTracer__Destructor(uniform RefCount* uniform _this)
{
  uniform PathTracer* uniform this = (uniform PathTracer* uniform) _this;
  PrecomputedSampler__Destructor(&this->sampler);
  RefCount__DecRef(&this->backplate->base);
  Renderer__Destructor(_this);
}

void PathTracer__Constructor(uniform PathTracer* uniform this,
                             const uniform int& maxDepth, 
                             const uniform float& minContribution,
                             const uniform float& epsilon,
                             const uniform int& spp,
                             uniform Image* uniform backplate)
{
  Renderer__Constructor(&this->base,PathTracer__Destructor,PathTracer_renderFrameInit,PathTracer_renderFrame);

  this->maxDepth = maxDepth;
  uniform int maxSupportedDepth = min(MAX_SAMPLES_2D,MAX_SAMPLES_2D-1);
  if (this->maxDepth > maxSupportedDepth) {
    print("WARNING: setting recursion depth to maximum of %\n",maxSupportedDepth);
    this->maxDepth = maxSupportedDepth;
  }

  this->minContribution = minContribution;
  this->epsilon = epsilon;
  this->spp = spp;
  this->iteration = 0;
  RefCount__IncRef(&backplate->base);
  this->backplate = backplate;

  this->lightSampleID = 0;
  this->firstScatterSampleID = 0;
  this->firstScatterTypeSampleID = 0;
  PrecomputedSampler__Constructor(&this->sampler,0,0);
}

export void* uniform PathTracer__new(const uniform int& maxDepth, 
                                     const uniform float& minContribution,
                                     const uniform float& epsilon,
                                     const uniform int& spp,
                                     void* uniform backplate)
{
  uniform PathTracer *uniform this = uniform new uniform PathTracer;
  PathTracer__Constructor(this,maxDepth,minContribution,epsilon,spp,(uniform Image* uniform) backplate);
  return this;
}
