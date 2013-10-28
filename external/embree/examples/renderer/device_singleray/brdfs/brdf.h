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

#ifndef __EMBREE_BRDF_H__
#define __EMBREE_BRDF_H__

#include "../shapes/differentialgeometry.h"

namespace embree {

    /*! BRDF type can be used as a hint to the integrator to help pick the best integration method  */
    enum BRDFType {

        ALL                   = 0xFFFFFFFF,    /*! all BRDF components for the given surface        */
        DIFFUSE               = 0x000F000F,    /*! all diffuse BRDFs for the given surface          */
        DIFFUSE_REFLECTION    = 0x00000001,    /*! fully diffuse reflection BRDF                    */
        DIFFUSE_TRANSMISSION  = 0x00010000,    /*! fully diffuse transmission BRDF                  */
        GLOSSY                = 0x00F000F0,    /*! all glossy BRDFs for the given surface           */
        GLOSSY_REFLECTION     = 0x00000010,    /*! semi-specular reflection BRDF                    */
        GLOSSY_TRANSMISSION   = 0x00100000,    /*! semi-specular transmission BRDF                  */
        JIMENEZ               = 0x02000002,    /*! all Jimenez BRDFs for the given surface          */
        JIMENEZ_REFLECTION    = 0x00000002,    /*! diffuse reflectance to be blurred [Jimenez 2009] */
        JIMENEZ_TRANSMISSION  = 0x02000000,    /*! approximate diffuse transmission  [Jimenez 2010] */
        NONE                  = 0x00000000,    /*! no BRDF components are set for the given surface */
        REFLECTION            = 0x0000FFFF,    /*! all reflection BRDFs for the given surface       */
        SPECULAR              = 0x0F000F00,    /*! all specular BRDFs for the given surface         */
        SPECULAR_REFLECTION   = 0x00000100,    /*! perfect specular reflection BRDF                 */
        SPECULAR_TRANSMISSION = 0x01000000,    /*! perfect specular transmission BRDF               */
        TRANSMISSION          = 0xFFFF0000     /*! all transmission BRDFs for the given surface     */

    };

    /*! BRDF interface definition.  A BRDF can be evaluated and sampled, and the sampling PDF       */
    /*! evaluated.  Unlike the definition in the literature, our BRDFs contain a cosine term.       */
    /*! For example, the diffuse BRDF in our system is "a / pi * cos(wi, N)" versus "a / pi".       */
    /*! This makes the formulae for reflection and refraction BRDFs more natural.  This also        */
    /*! adds consistency since the sampling functionality of the BRDF class would have to           */
    /*! include sampling of the cosine term anyway.  As an optimization, evaluation of the          */
    /*! cosine term can be skipped in case it is handled through sampling.                          */
    class BRDF {
    public:

        /*! constructor */
        __forceinline BRDF(const BRDFType type) : type(type) {}

        /*! virtual destructor */
        virtual ~BRDF() {}

        /*! evaluate the BRDF */
        virtual Color eval(const Vector3f               & wo,    /*! outgoing light direction          */
                           const DifferentialGeometry& dg,    /*! shade location on a surface       */
                           const Vector3f               & wi)    /*! incoming light direction          */ const {

            /*! specular BRDFs cannot be evaluated and should return zero */
            return(zero);

        }

        /*! sample the BRDF */
        virtual Color sample(const Vector3f               & wo,  /*! outgoing light direction          */
                             const DifferentialGeometry& dg,  /*! shade location on a surface       */
                             Sample3f                  & wi,  /*! sampled light direction and PDF   */
                             const Vec2f               & s)   /*! sample location given by caller   */ const {

            /*! by default we perform a cosine weighted hemisphere sampling */
            return(eval(wo, dg, wi = cosineSampleHemisphere(s.x, s.y, dg.Ns)));

        }

        /*! evaluate sampling PDF */
        virtual float pdf(const Vector3f               & wo,     /*! outgoing light direction          */
                          const DifferentialGeometry& dg,     /*! shade location on a surface       */
                          const Vector3f               & wi)     /*! incoming light direction          */ const {

            /*! return the probability density */
            return(cosineSampleHemispherePDF(wi, dg.Ns));

        }

        /*! BRDF type hint to the integrator */
        BRDFType type;

    };

}

#endif // __EMBREE_BRDF_H__

