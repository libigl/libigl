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

#ifndef __EMBREE_SAMPLER_H__
#define __EMBREE_SAMPLER_H__

#include "default.h"
#include "math/random.h"
#include "lights/light.h"
#include "filters/filter.h"

namespace embree
{
  /*! A precomputed sample on light source. */
  struct LightSample 
  {
    ALIGNED_CLASS;
  public:
    Sample3f wi; //!< The direction towards the light source.
    float tMax;  //!< The largest valid parameter value for a shadow ray.
    Color L;     //!< The importance weighted radiance for this sample.
  };

  /*! A complete high-dimensional sample, attached to one sample in the image plane. */
  struct PrecomputedSample 
  {
    /*! Get the raster position of the current sample. */
    Vec2f getPrimary() const { return raster; }

    /*! Get the integer raster coordinates of the current sample. */
    Vec2i getIntegerRaster() const { return integerRaster; }

    /*! Get the integer raster coordinates of the current sample. */
    Vec2i getPixel() const { return integerRaster; }

    /*! Get the current lens sample. */
    Vec2f getLens() const { return lens; }

    /*! Get the current time sample. */
    float getTime() const { return time; }

    /*! Get the specified additional 1D sample for the current sample. */
    float getFloat (int dim) const { return samples1D[dim]; }

    /*! Get the specified additional 2D sample for the current sample. */
    Vec2f getVec2f(int dim) const { return samples2D[dim]; }

    /*! Get the specified precomputed light sample for the current sample. */
    LightSample getLightSample(int lightSampleId) const { return lightSamples[lightSampleId]; }

    Vec2f pixel;               //!< Sample location inside the pixel. [0.5;0.5] is the pixel center.
    Vec2f raster;              //!< Floating point coordinates of the sample in [0;width]x[0;height].
    Vec2i integerRaster;       //!< Integer coordinates of the associated pixel in [0;width]x[0;height].
    float time;                //!< time sample for motion blur.
    Vec2f lens;                //!< 2D lens sample for depth of field.
    float* samples1D;          //!< Additional 1D samples requested by the integrator.
    Vec2f* samples2D;          //!< Additional 2D samples requested by the integrator.
    LightSample* lightSamples; //!< Precomputed light samples.
    Vec2i imageSize;           //!< Dimesions of the entire image.
  };

  class SamplerFactory;

  /*! A sampler thread responsible for sampling a complete tile.
   *  The sampler uses precomputed samlpes from a sampler factory. */
  class Sampler 
  {
  public:
    /*! Create a sampler using the specified sampler factory. */
    Sampler(const Ref<SamplerFactory>& factory) : factory(factory) {}

    /*! Initialize the sampler for a given tile. */
    void init(const Vec2i& imageSize, const Vec2i& tileBegin,
              const Vec2i& tileEnd, int iteration = 0);

    /*! Proceed to the next sample. */
    void proceed(PrecomputedSample& sample);

    /*! Has the tile been sampled completely? */
    bool finished() const { return done; }

    /*! Get dimensions of entire image. */
    const Vec2i& getImageSize() const { return imageSize; }

    /*! Get coordinates of first pixel in tile. */
    const Vec2i& getTileBegin() const { return tileBegin; }

    /*! Get coordinates of last pixel in tile. */
    const Vec2i& getTileEnd() const { return tileEnd; }

  protected:
    Vec2i imageSize; //!< Dimesions of the entire image.
    Vec2i tileBegin; //!< Coordinates of first pixel in current tile.
    Vec2i tileEnd;   //!< Coordinates of last pixel in current tile.
    int iteration;   //!< Current iteration.

    Random randomNumberGenerator; //!< Random number generated used by this sampler thread.
    Ref<SamplerFactory> factory;  //!< A reference to the sampler factory shared by all threads.
    bool done;                    //!< Has the tile been sampled completely?
    Vec2i currentPixel;           //!< Coordinates of the currently sampled pixel.
    int currentSample;            //!< Index of current sample in current pixel.
    int currentSet;               //!< Index of the precomputed sample set that is used for current pixel.
  };

  /*! The sampler factory precomputes samples for usage by multiple samlper threads. */
  class SamplerFactory : public RefCount 
  {
    friend class Sampler;
  public:
    /*! Construction from parameters. */
    SamplerFactory(const Parms& parms);

    /*! Construction from number of samples per pixel and number of precomputed sets. */
    SamplerFactory(const unsigned samplesPerPixel = 1,
                   const unsigned sampleSets = 64);
    
    /*! Destructor */
    virtual ~SamplerFactory();

    /*! Request additional 1D samples per pixel sample. */
    int request1D(int num = 1);

    /*! Request additional 2D samples per pixel sample. */
    int request2D(int num = 1);

    /*! Request a precomputed light sample. */
    int requestLightSample(int baseSample, const Ref<Light>& light);

    /*! Reset the sampler factory. Delete all precomputed samples. */
    void reset();

    /*! Initialize the factory for a given iteration and precompute all samples. */
    void init(int iteration = 0, const Ref<Filter> filter = NULL);

    /*! Create a sampler thread using this factory. */
    Sampler* create();

  public:
    int numSamples1D;                  //!< Number of additional 1D samples per pixel sample.
    int numSamples2D;                  //!< Number of additional 2D samples per pixel sample.
    int numLightSamples;               //!< Number of precomputed light samples per pixel sample.
    std::vector<Ref<Light> > lights;   //!< References to all light sources.
    std::vector<int> lightBaseSamples; //!< Base samples for light sample precomputation.

    int samplesPerPixel;               //!< Number of samples per pixel.
    int sampleSets;                    //!< Number of precomputed sample sets.
    PrecomputedSample** samples;       //!< All precomputed samples.
    int iteration;                     //!< Current iteration.
  };
}

#endif
