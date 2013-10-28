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

#include "math/permutation.h"
#include "samplers/sampler.h"
#include "samplers/patterns.h"

namespace embree
{
  SamplerFactory::SamplerFactory(const Parms& parms)
    : numSamples1D(0), numSamples2D(0), numLightSamples(0),
      samplesPerPixel(1), sampleSets(64), samples(NULL)
  {
    samplesPerPixel = parms.getInt("sampler.spp",1);
    sampleSets      = parms.getInt("sampler.sets",64);
  }

  SamplerFactory::SamplerFactory(const unsigned samplesPerPixel,
                                 const unsigned sampleSets)
    : numSamples1D(0), numSamples2D(0), numLightSamples(0),
      samplesPerPixel(samplesPerPixel), sampleSets(sampleSets), samples(NULL) 
  {
  }

  SamplerFactory::~SamplerFactory() {
    reset();
  }

  int SamplerFactory::request1D(int num)
  {
    int dim = numSamples1D;
    numSamples1D += num;
    return dim;
  }

  int SamplerFactory::request2D(int num)
  {
    int dim = numSamples2D;
    numSamples2D += num;
    return dim;
  }

  int SamplerFactory::requestLightSample(int baseSample, const Ref<Light>& light)
  {
    numLightSamples++;
    lights.push_back(light);
    lightBaseSamples.push_back(baseSample);
    return numLightSamples-1;
  }

  void SamplerFactory::reset()
  {
    if (samples) {
      for (int set = 0; set < sampleSets; set++) {
        for (int s = 0; s < samplesPerPixel; s++) {
          delete[] samples[set][s].samples1D;
          delete[] samples[set][s].samples2D;
          delete[] samples[set][s].lightSamples;
        }
        delete[] samples[set];
      }
      delete[] samples;
    }
    samples = NULL;
    numSamples1D = 0;
    numSamples2D = 0;
    numLightSamples = 0;
  }

  void SamplerFactory::init(int iteration, const Ref<Filter> filter)
  {
    this->iteration = iteration;
    samples = new PrecomputedSample*[sampleSets];
    if (samplesPerPixel != (1 << __bsf(samplesPerPixel)))
      throw std::runtime_error("Number of samples per pixel have to be a power of two.");

    int chunkSize = max((int)samplesPerPixel,64);
    int currentChunk = int(iteration*samplesPerPixel) / chunkSize;
    int offset = (iteration*samplesPerPixel) % chunkSize;
    Random rng;
    rng.setSeed(currentChunk * 5897);

    Vec2f* pixel = new Vec2f[chunkSize];
    float* time = new float[chunkSize];
    Vec2f* lens = new Vec2f[chunkSize];
    float* samples1D = new float[chunkSize];
    Vec2f* samples2D = new Vec2f[chunkSize];

    for (int set = 0; set < sampleSets; set++)
    {
      samples[set] = new PrecomputedSample[samplesPerPixel];

      /*! Generate pixel and lens samples. */
      multiJittered(pixel, chunkSize, rng);
      jittered(time, chunkSize, rng);
      multiJittered(lens, chunkSize, rng);
      for (int s = 0; s < samplesPerPixel; s++) {
        samples[set][s].pixel = pixel[offset + s];
        samples[set][s].time = time[offset + s];
        samples[set][s].lens = lens[offset + s];
        if (filter) {
          samples[set][s].pixel = filter->sample(samples[set][s].pixel) + Vec2f(0.5f, 0.5f);
        }
        samples[set][s].samples1D = new float[SamplerFactory::numSamples1D];
        samples[set][s].samples2D = new Vec2f[SamplerFactory::numSamples2D];
        samples[set][s].lightSamples = new LightSample[SamplerFactory::numLightSamples];
      }

      /*! Generate requested 1D samples. */
      for (int d = 0; d < SamplerFactory::numSamples1D; d++) {
        jittered(samples1D, chunkSize, rng);
        for (int s = 0; s < samplesPerPixel; s++) {
          samples[set][s].samples1D[d] = samples1D[offset + s];
        }
      }

      /*! Generate 2D samples. */
      for (int d = 0; d < SamplerFactory::numSamples2D; d++) {
        multiJittered(samples2D, chunkSize, rng);
        for (int s = 0; s < samplesPerPixel; s++) {
          samples[set][s].samples2D[d] = samples2D[offset + s];
        }
      }

      /*! Generate light samples. */
      for (int d = 0; d < SamplerFactory::numLightSamples; d++) {
        for (int s = 0; s < samplesPerPixel; s++) {
          LightSample ls;
          DifferentialGeometry dg;
          ls.L = lights[d]->sample(dg, ls.wi, ls.tMax, samples[set][s].samples2D[lightBaseSamples[d]]);
          samples[set][s].lightSamples[d] = ls;
        }
      }
    }
    
    delete[] pixel;
    delete[] time;
    delete[] lens;
    delete[] samples1D;
    delete[] samples2D;
  }

  Sampler* SamplerFactory::create() {
    return new Sampler(this);
  }

  void Sampler::init(const Vec2i& imageSize, const Vec2i& tileBegin,
                     const Vec2i& tileEnd, int iteration)
  {
    this->imageSize = imageSize;
    this->tileBegin = tileBegin;
    this->tileEnd = tileEnd;
    this->iteration = iteration;
    done = false;
    randomNumberGenerator.setSeed(tileBegin.x * 91711 + tileBegin.y * 81551);
    currentPixel = tileBegin;
    currentSample = 0;
    currentSet = randomNumberGenerator.getInt(factory->sampleSets);
  }

  void Sampler::proceed(PrecomputedSample& sample)
  {
    sample = factory->samples[currentSet][currentSample];
    sample.integerRaster = currentPixel;
    sample.raster.x = currentPixel.x + sample.pixel.x;
    sample.raster.y = currentPixel.y + sample.pixel.y;
    sample.imageSize = imageSize;

    ++currentSample;
    if (currentSample == factory->samplesPerPixel) {
      currentSample = 0;
      currentSet = randomNumberGenerator.getInt(factory->sampleSets);
      ++currentPixel.x;
      if (currentPixel.x > tileEnd.x) {
        currentPixel.x = tileBegin.x;
        ++currentPixel.y;
        if (currentPixel.y > tileEnd.y) {
          done = true;
        }
      }
    }
  }
}
