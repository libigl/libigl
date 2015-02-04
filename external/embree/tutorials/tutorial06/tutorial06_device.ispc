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

#include "../common/tutorial/tutorial_device.isph"
#include "../common/tutorial/scene_device.isph"
#include "shapesampler.isph"
#include "optics.isph"

#define __device__

#define PARALLEL_COMMIT

struct DifferentialGeometry
{
  Vec3f P;
  Vec3f Ng;
  Vec3f Ns;
};

struct BRDF
{
  float Ns;               /*< specular exponent */
  float Ni;               /*< optical density for the surface (index of refraction) */
  Vec3f Ka;              /*< ambient reflectivity */
  Vec3f Kd;              /*< diffuse reflectivity */
  Vec3f Ks;              /*< specular reflectivity */
  Vec3f Kt;              /*< transmission filter */
};

struct Medium
{
  Vec3f transmission; //!< Transmissivity of medium.
  float eta;             //!< Refraction index of medium.
};

inline Medium make_Medium(const varying Vec3f& transmission, const float eta)
{
  Medium m;
  m.transmission = transmission;
  m.eta = eta;
  return m;
}

inline Medium make_Medium_Vacuum() { 
  return make_Medium(make_Vec3f((varying float)1.0f),1.0f); 
}

inline bool eq(const Medium& a, const Medium& b) {
  return (a.eta == b.eta) && eq(a.transmission, b.transmission);
}

inline Vec3f sample_component2(const Vec3f& c0, const Sample3f& wi0, const Medium& medium0,
                               const Vec3f& c1, const Sample3f& wi1, const Medium& medium1,
                               const Vec3f& Lw, Sample3f& wi_o, Medium& medium_o, const float s)
{
  const Vec3f m0 = Lw*c0/wi0.pdf;
  const Vec3f m1 = Lw*c1/wi1.pdf;

  const float C0 = wi0.pdf == 0.0f ? 0.0f : max(max(m0.x,m0.y),m0.z);
  const float C1 = wi1.pdf == 0.0f ? 0.0f : max(max(m1.x,m1.y),m1.z);
  const float C  = C0 + C1;

  if (C == 0.0f) {
    wi_o = make_Sample3f(make_Vec3f(0,0,0),0);
    return make_Vec3f(0,0,0);
  }

  const float CP0 = C0/C;
  const float CP1 = C1/C;
  if (s < CP0) {
    wi_o = make_Sample3f(wi0.v,wi0.pdf*CP0); 
    medium_o = medium0; return c0;
  } 
  else {
    wi_o = make_Sample3f(wi1.v,wi1.pdf*CP1); 
    medium_o = medium1; return c1;
  }
}

////////////////////////////////////////////////////////////////////////////////
//                             Ambient Light                                  //
////////////////////////////////////////////////////////////////////////////////

inline Vec3f AmbientLight__eval(const uniform ISPCAmbientLight& light, const Vec3f& wo) {
  return make_Vec3f(light.L);
}

inline Vec3f AmbientLight__sample(const uniform ISPCAmbientLight& light, const DifferentialGeometry& dg, varying Sample3f& wi, varying float& tMax, varying const Vec2f& s) 
{
  wi = cosineSampleHemisphere(s.x,s.y,dg.Ns);
  tMax = 1e20f;
  return make_Vec3f(light.L);
}

////////////////////////////////////////////////////////////////////////////////
//                             Point Light                                    //
////////////////////////////////////////////////////////////////////////////////

inline varying Vec3f PointLight__sample(const uniform ISPCPointLight& light, 
					varying const DifferentialGeometry& dg, 
					varying Sample3f& wi,
					varying float& tMax,
					varying const Vec2f& s) 
{
  Vec3f d = make_Vec3f(light.P) - dg.P;
  float distance = length(d);
  wi = make_Sample3f(d*rcp(distance), distance*distance);
  tMax = distance;
  return make_Vec3f(light.I);
}

////////////////////////////////////////////////////////////////////////////////
//                        Directional Light                                   //
////////////////////////////////////////////////////////////////////////////////

inline varying Vec3f DirectionalLight__sample(const uniform ISPCDirectionalLight& light, 
					      varying const DifferentialGeometry& dg, 
					      varying Sample3f& wi,
					      varying float& tMax,
					      varying const Vec2f& s) 
{
  wi = make_Sample3f(neg(normalize(make_Vec3f(light.D))),1.0f); 
  tMax = inf; 
  return make_Vec3f(light.E);
}

////////////////////////////////////////////////////////////////////////////////
//                          Distant Light                                     //
////////////////////////////////////////////////////////////////////////////////

inline varying Vec3f DistantLight__eval(const uniform ISPCDistantLight& light, varying const Vec3f& wo) 
{
  if (-dot(wo,make_Vec3f(light.D)) >= light.cosHalfAngle) return make_Vec3f(light.L);
  return make_Vec3f(0.0f);
}

inline varying Vec3f DistantLight__sample(const uniform ISPCDistantLight& light,
                                   varying const DifferentialGeometry& dg, 
                                   varying Sample3f& wi,
                                   varying float& tMax,
                                   varying const Vec2f& s) 
{
  wi = UniformSampleCone(s.x,s.y,light.radHalfAngle,make_Vec3f((varying Vec3fa)neg(light.D)));
  tMax = 1e20f;
  return make_Vec3f(light.L);
}

////////////////////////////////////////////////////////////////////////////////
//                          Minneart BRDF                                     //
////////////////////////////////////////////////////////////////////////////////

struct Minneart
{
  /*! The reflectance parameter. The vale 0 means no reflection,
   *  and 1 means full reflection. */
  Vec3f R;
  
  /*! The amount of backscattering. A value of 0 means lambertian
   *  diffuse, and inf means maximum backscattering. */
  float b;
};

inline Vec3f Minneart__eval(const varying Minneart* uniform This,
                     const Vec3f &wo, const DifferentialGeometry &dg, const Vec3f &wi) 
{
  const float cosThetaI = clamp(dot(wi,dg.Ns));
  const float backScatter = pow(clamp(dot(wo,wi)), This->b);
  return (backScatter * cosThetaI * one_over_pi) * This->R;
}

inline Vec3f Minneart__sample(const varying Minneart* uniform This,
                       const Vec3f &wo, 
                       const DifferentialGeometry &dg, 
                       Sample3f &wi, 
                       const Vec2f &s)  
{
  wi = cosineSampleHemisphere(s.x,s.y,dg.Ns);
  return Minneart__eval(This, wo, dg, wi.v);
}

inline void Minneart__Constructor(varying Minneart* uniform This, const varying Vec3f& R, const varying float b) 
{
  This->R = R;
  This->b = b;
}

inline varying Minneart make_Minneart(const varying Vec3f& R, const varying float f) { 
  varying Minneart m; Minneart__Constructor(&m,R,f); return m; 
}

////////////////////////////////////////////////////////////////////////////////
//                            Velvet BRDF                                     //
////////////////////////////////////////////////////////////////////////////////

struct Velvety
{
  uniform BRDF base;

  /*! The reflectance parameter. The vale 0 means no reflection,
   *  and 1 means full reflection. */
  Vec3f R;
  
  /*! The falloff of horizon scattering. 0 no falloff,
   *  and inf means maximum falloff. */
  float f;
};

inline Vec3f Velvety__eval(const varying Velvety* uniform This,
                    const Vec3f &wo, const DifferentialGeometry &dg, const Vec3f &wi) 
{
  const float cosThetaO = clamp(dot(wo,dg.Ns));
  const float cosThetaI = clamp(dot(wi,dg.Ns));
  const float sinThetaO = sqrt(1.0f - cosThetaO * cosThetaO);
  const float horizonScatter = pow(sinThetaO, This->f);
  return (horizonScatter * cosThetaI * one_over_pi) * This->R;
}

inline Vec3f Velvety__sample(const varying Velvety* uniform This,
                      const Vec3f &wo, 
                      const DifferentialGeometry &dg, 
                      Sample3f &wi, 
                      const Vec2f &s)  
{
  wi = cosineSampleHemisphere(s.x,s.y,dg.Ns);
  return Velvety__eval(This, wo, dg, wi.v);
}

inline void Velvety__Constructor(varying Velvety* uniform This, const varying Vec3f& R, const varying float f) 
{
  This->R = R;
  This->f = f;
}

inline varying Velvety make_Velvety(const varying Vec3f& R, const varying float f) { 
  varying Velvety m; Velvety__Constructor(&m,R,f); return m; 
}

////////////////////////////////////////////////////////////////////////////////
//                  Dielectric Reflection BRDF                                //
////////////////////////////////////////////////////////////////////////////////

struct DielectricReflection
{
  float eta;
};

inline Vec3f DielectricReflection__eval(const varying DielectricReflection* uniform This, const Vec3f &wo, const DifferentialGeometry &dg, const Vec3f &wi) {
  return make_Vec3f(0.f);
}

inline Vec3f DielectricReflection__sample(const varying DielectricReflection* uniform This, const Vec3f &wo, const DifferentialGeometry &dg, Sample3f &wi, const Vec2f &s)
{
  const float cosThetaO = clamp(dot(wo,dg.Ns));
  wi = reflect_(wo,dg.Ns,cosThetaO);
  return make_Vec3f(fresnelDielectric(cosThetaO,This->eta));
}

inline void DielectricReflection__Constructor(varying DielectricReflection* uniform This,
                                              const varying float etai,
                                              const varying float etat)
{
  This->eta = etai*rcp(etat);
}

inline varying DielectricReflection make_DielectricReflection(const varying float etai, const varying float etat) {
  varying DielectricReflection v; DielectricReflection__Constructor(&v,etai,etat); return v;
}

////////////////////////////////////////////////////////////////////////////////
//                                Lambertian BRDF                             //
////////////////////////////////////////////////////////////////////////////////

struct Lambertian
{
  Vec3f R;
};

inline Vec3f Lambertian__eval(const varying Lambertian* uniform This,
                              const Vec3f &wo, const DifferentialGeometry &dg, const Vec3f &wi) 
{
  return This->R * (1.0f/(float)(M_PI)) * clamp(dot(wi,dg.Ns));
}

inline Vec3f Lambertian__sample(const varying Lambertian* uniform This,
                                const Vec3f &wo, 
                                const DifferentialGeometry &dg, 
                                Sample3f &wi, 
                                const Vec2f &s)  
{
  wi = cosineSampleHemisphere(s.x,s.y,dg.Ns);
  return Lambertian__eval(This, wo, dg, wi.v);
}

inline void Lambertian__Constructor(varying Lambertian* uniform This, const varying Vec3f& R)
{
  This->R = R;
}

inline varying Lambertian make_Lambertian(const varying Vec3f& R) {
  varying Lambertian v; Lambertian__Constructor(&v,R); return v;
}


////////////////////////////////////////////////////////////////////////////////
//              Lambertian BRDF with Dielectric Layer on top                  //
////////////////////////////////////////////////////////////////////////////////

struct DielectricLayerLambertian
{
  Vec3f T;             //!< Transmission coefficient of dielectricum
  float etait;         //!< Relative refraction index etai/etat of both media
  float etati;         //!< relative refraction index etat/etai of both media
  Lambertian ground;   //!< the BRDF of the ground layer
};

inline Vec3f DielectricLayerLambertian__eval(const varying DielectricLayerLambertian* uniform This,
                                             const Vec3f &wo, const DifferentialGeometry &dg, const Vec3f &wi) 
{
  const float cosThetaO = dot(wo,dg.Ns);
  const float cosThetaI = dot(wi,dg.Ns);
  if (cosThetaI <= 0.0f | cosThetaO <= 0.0f) return make_Vec3f(0.f);

  float cosThetaO1; 
  const Sample3f wo1 = refract(wo,dg.Ns,This->etait,cosThetaO,cosThetaO1);
  float cosThetaI1; 
  const Sample3f wi1 = refract(wi,dg.Ns,This->etait,cosThetaI,cosThetaI1);
  const float Fi = 1.0f - fresnelDielectric(cosThetaI,cosThetaI1,This->etait);
  const Vec3f Fg = Lambertian__eval(&This->ground,neg(wo1.v),dg,neg(wi1.v));
  const float Fo = 1.0f - fresnelDielectric(cosThetaO,cosThetaO1,This->etait);
  return Fo * This->T * Fg * This->T * Fi;
}

inline Vec3f DielectricLayerLambertian__sample(const varying DielectricLayerLambertian* uniform This,
                                               const Vec3f &wo, 
                                               const DifferentialGeometry &dg, 
                                               Sample3f &wi, 
                                               const Vec2f &s)  
{
  /*! refract ray into medium */
  float cosThetaO = dot(wo,dg.Ns);
  if (cosThetaO <= 0.0f) return make_Vec3f(0.f);
  float cosThetaO1; Sample3f wo1 = refract(wo,dg.Ns,This->etait,cosThetaO,cosThetaO1);
  
  /*! sample ground BRDF */
  Sample3f wi1 = make_Sample3f(make_Vec3f(0.f),1.f); 
  Vec3f Fg = Lambertian__sample(&This->ground,neg(wo1.v),dg,wi1,s);

  /*! refract ray out of medium */
  float cosThetaI1 = dot(wi1.v,dg.Ns);
  if (cosThetaI1 <= 0.0f) return make_Vec3f(0.f);
  
  float cosThetaI; 
  Sample3f wi0 = refract(neg(wi1.v),neg(dg.Ns),This->etati,cosThetaI1,cosThetaI);
  if (wi0.pdf == 0.0f) return make_Vec3f(0.f);
  
  /*! accumulate contribution of path */
  wi = make_Sample3f(wi0.v,wi1.pdf);
  float Fi = 1.0f - fresnelDielectric(cosThetaI,cosThetaI1,This->etait);
  float Fo = 1.0f - fresnelDielectric(cosThetaO,cosThetaO1,This->etait);
  return Fo * This->T * Fg * This->T * Fi;
}

inline void DielectricLayerLambertian__Constructor(varying DielectricLayerLambertian* uniform This,
                                                   const varying Vec3f& T, 
                                                   const varying float etai, 
                                                   const varying float etat, 
                                                   const varying Lambertian& ground)
{
  This->T = T;
  This->etait = etai*rcp(etat);
  This->etati = etat*rcp(etai);
  This->ground = ground;
}

inline varying DielectricLayerLambertian make_DielectricLayerLambertian(const varying Vec3f& T, 
                                                                        const varying float etai, 
                                                                        const varying float etat, 
                                                                        const varying Lambertian& ground)
{
  varying DielectricLayerLambertian m; 
  DielectricLayerLambertian__Constructor(&m,T,etai,etat,ground);
  return m;
}

////////////////////////////////////////////////////////////////////////////////
//                          Matte Material                                    //
////////////////////////////////////////////////////////////////////////////////

__device__ void MatteMaterial__preprocess(uniform MatteMaterial* uniform material, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

__device__ Vec3f MatteMaterial__eval(MatteMaterial* uniform This, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) 
{
  Lambertian lambertian = make_Lambertian(make_Vec3f((varying Vec3fa)This->reflectance));
  return Lambertian__eval(&lambertian,wo,dg,wi);
}

__device__ Vec3f MatteMaterial__sample(MatteMaterial* uniform This, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  Lambertian lambertian = make_Lambertian(make_Vec3f((varying Vec3fa)This->reflectance));
  return Lambertian__sample(&lambertian,wo,dg,wi_o,s);
}

////////////////////////////////////////////////////////////////////////////////
//                          Mirror Material                                    //
////////////////////////////////////////////////////////////////////////////////

__device__ void MirrorMaterial__preprocess(uniform MirrorMaterial* uniform material, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

__device__ Vec3f MirrorMaterial__eval(MirrorMaterial* uniform This, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) {
  return make_Vec3f(0.0f);
}

__device__ Vec3f MirrorMaterial__sample(MirrorMaterial* uniform This, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  wi_o = reflect_(wo,dg.Ns);
  return make_Vec3f(This->reflectance);
}

////////////////////////////////////////////////////////////////////////////////
//                          OBJ Material                                      //
////////////////////////////////////////////////////////////////////////////////

__device__ void OBJMaterial__preprocess(uniform OBJMaterial* uniform material, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
    float d = material->d;
    //if (material->map_d) { d *= material->map_d.get(s,t); }
    brdf.Ka = make_Vec3f(material->Ka);
    //if (material->map_Ka) { brdf.Ka *= material->map_Ka->get(dg.st); }
    brdf.Kd = d * make_Vec3f(material->Kd);  
    //if (material->map_Kd) brdf.Kd *= material->map_Kd->get(dg.st);  
    brdf.Ks = d * make_Vec3f(material->Ks);  
    //if (material->map_Ks) brdf.Ks *= material->map_Ks->get(dg.st); 
    brdf.Ns = material->Ns;  
    //if (material->map_Ns) { brdf.Ns *= material->map_Ns.get(dg.st); }
    brdf.Kt = (1.0f-d)*make_Vec3f(material->Kt);
    brdf.Ni = material->Ni;
}

__device__ Vec3f OBJMaterial__eval(OBJMaterial* uniform material, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) 
{
  Vec3f R = make_Vec3f(0.0f,0.0f,0.0f);
  const float Md = max(max(brdf.Kd.x,brdf.Kd.y),brdf.Kd.z);
  const float Ms = max(max(brdf.Ks.x,brdf.Ks.y),brdf.Ks.z);
  const float Mt = max(max(brdf.Kt.x,brdf.Kt.y),brdf.Kt.z);
  if (Md > 0.0f) {
    R = R + (1.0f/M_PI) * clamp(dot(wi,dg.Ns)) * brdf.Kd; // FIXME: +=
  }
  if (Ms > 0.0f) {
    const Sample3f refl = reflect_(wo,dg.Ns);
    if (dot(refl.v,wi) > 0.0f) 
      R = R + (brdf.Ns+2) * one_over_two_pi * pow(max(1e-10f,dot(refl.v,wi)),brdf.Ns) * clamp(dot(wi,dg.Ns)) * brdf.Ks; // FIXME: +=
  }
  if (Mt > 0.0f) {
  }
  return R;
}

__device__ Vec3f OBJMaterial__sample(OBJMaterial* uniform material, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  Vec3f cd = make_Vec3f(0.0f); 
  Sample3f wid = make_Sample3f(make_Vec3f(0.0f),0.0f);
  if (max(max(brdf.Kd.x,brdf.Kd.y),brdf.Kd.z) > 0.0f) {
    wid = cosineSampleHemisphere(s.x,s.y,dg.Ns);
    cd = one_over_pi * clamp(dot(wid.v,dg.Ns)) * brdf.Kd;
  }

  Vec3f cs = make_Vec3f(0.0f); 
  Sample3f wis = make_Sample3f(make_Vec3f(0.0f),0.0f);
  if (max(max(brdf.Ks.x,brdf.Ks.y),brdf.Ks.z) > 0.0f)
  {
    const Sample3f refl = reflect_(wo,dg.Ns);
    wis = powerCosineSampleHemisphere(s.x,s.y,refl.v,brdf.Ns);
    cs = (brdf.Ns+2) * one_over_two_pi * pow(dot(refl.v,wis.v),brdf.Ns) * clamp(dot(wis.v,dg.Ns)) * brdf.Ks;
  }

  Vec3f ct = make_Vec3f(0.0f); 
  Sample3f wit = make_Sample3f(make_Vec3f(0.0f),0.0f);
  if (max(max(brdf.Kt.x,brdf.Kt.y),brdf.Kt.z) > 0.0f)
  {
    wit = make_Sample3f(neg(wo),1.0f);
    ct = brdf.Kt;
  }

  const Vec3f md = Lw*cd/wid.pdf;
  const Vec3f ms = Lw*cs/wis.pdf;
  const Vec3f mt = Lw*ct/wit.pdf;

  const float Cd = wid.pdf == 0.0f ? 0.0f : max(max(md.x,md.y),md.z);
  const float Cs = wis.pdf == 0.0f ? 0.0f : max(max(ms.x,ms.y),ms.z);
  const float Ct = wit.pdf == 0.0f ? 0.0f : max(max(mt.x,mt.y),mt.z);
  const float C  = Cd + Cs + Ct;

  if (C == 0.0f) {
    wi_o = make_Sample3f(make_Vec3f(0,0,0),0);
    return make_Vec3f(0,0,0);
  }

  const float CPd = Cd/C;
  const float CPs = Cs/C;
  const float CPt = Ct/C;

  if (s.x < CPd) {
    wi_o = make_Sample3f(wid.v,wid.pdf*CPd);
    return cd;
  } 
  else if (s.x < CPd + CPs)
  {
    wi_o = make_Sample3f(wis.v,wis.pdf*CPs);
    return cs;
  }
  else 
  {
    wi_o = make_Sample3f(wit.v,wit.pdf*CPt);
    return ct;
  }
}

////////////////////////////////////////////////////////////////////////////////
//                        Metal Material                                      //
////////////////////////////////////////////////////////////////////////////////

__device__ void MetalMaterial__preprocess(uniform MetalMaterial* uniform material, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

__device__ Vec3f MetalMaterial__eval(MetalMaterial* uniform This, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) 
{
  const FresnelConductor fresnel = make_FresnelConductor(make_Vec3f(This->eta),make_Vec3f(This->k));
  const PowerCosineDistribution distribution = make_PowerCosineDistribution(rcp(This->roughness));

  const float cosThetaO = dot(wo,dg.Ns);
  const float cosThetaI = dot(wi,dg.Ns);
  if (cosThetaI <= 0.0f | cosThetaO <= 0.0f) return make_Vec3f(0.f);
  const Vec3f wh = normalize(wi+wo);
  const float cosThetaH = dot(wh, dg.Ns);
  const float cosTheta = dot(wi, wh); // = dot(wo, wh);
  const Vec3f F = eval(fresnel,cosTheta);
  const float D = eval(distribution,cosThetaH);
  const float G = min(1.0f, min(2.0f * cosThetaH * cosThetaO / cosTheta, 
                                2.0f * cosThetaH * cosThetaI / cosTheta));
  return (make_Vec3f(This->reflectance)*F) * D * G * rcp(4.0f*cosThetaO);
}

__device__ Vec3f MetalMaterial__sample(MetalMaterial* uniform This, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  const PowerCosineDistribution distribution = make_PowerCosineDistribution(rcp(This->roughness));

  if (dot(wo,dg.Ns) <= 0.0f) return make_Vec3f(0.0f);
  sample(distribution,wo,dg.Ns,wi_o,s);
  if (dot(wi_o.v,dg.Ns) <= 0.0f) return make_Vec3f(0.0f);
  return MetalMaterial__eval(This,brdf,wo,dg,wi_o.v);
}

////////////////////////////////////////////////////////////////////////////////
//                        ReflectiveMetal Material                            //
////////////////////////////////////////////////////////////////////////////////

__device__ void ReflectiveMetalMaterial__preprocess(uniform ReflectiveMetalMaterial* uniform material, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  {
}

__device__ Vec3f ReflectiveMetalMaterial__eval(ReflectiveMetalMaterial* uniform This, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) {
  return make_Vec3f(0.0f);
}

__device__ Vec3f ReflectiveMetalMaterial__sample(ReflectiveMetalMaterial* uniform This, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  wi_o = reflect_(wo,dg.Ns);
  return make_Vec3f(This->reflectance) * fresnelConductor(dot(wo,dg.Ns),make_Vec3f((varying Vec3fa)This->eta),make_Vec3f((varying Vec3fa)This->k));
}

////////////////////////////////////////////////////////////////////////////////
//                        Velvet Material                                     //
////////////////////////////////////////////////////////////////////////////////

__device__ void VelvetMaterial__preprocess(uniform VelvetMaterial* uniform material, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

__device__ Vec3f VelvetMaterial__eval(VelvetMaterial* uniform This, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) 
{
  Minneart minneart; Minneart__Constructor(&minneart,(varying Vec3f)make_Vec3f(This->reflectance),This->backScattering);
  Velvety velvety; Velvety__Constructor (&velvety,make_Vec3f((varying Vec3fa)This->horizonScatteringColor),This->horizonScatteringFallOff);
  return Minneart__eval(&minneart,wo,dg,wi) + Velvety__eval(&velvety,wo,dg,wi);
}

__device__ Vec3f VelvetMaterial__sample(VelvetMaterial* uniform This, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  Minneart minneart; Minneart__Constructor(&minneart,make_Vec3f((varying Vec3fa)This->reflectance),This->backScattering);
  Velvety velvety; Velvety__Constructor (&velvety,make_Vec3f((varying Vec3fa)This->horizonScatteringColor),This->horizonScatteringFallOff);

  Sample3f wi0; Vec3f c0 = Minneart__sample(&minneart,wo,dg,wi0,s);
  Sample3f wi1; Vec3f c1 = Velvety__sample(&velvety,wo,dg,wi1,s);
  return sample_component2(c0,wi0,medium,c1,wi1,medium,Lw,wi_o,medium,s.x);
}

////////////////////////////////////////////////////////////////////////////////
//                          Dielectric Material                               //
////////////////////////////////////////////////////////////////////////////////

__device__ void DielectricMaterial__preprocess(uniform DielectricMaterial* uniform material, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

__device__ Vec3f DielectricMaterial__eval(DielectricMaterial* uniform material, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) {
  return make_Vec3f(0.0f);
}

__device__ Vec3f DielectricMaterial__sample(DielectricMaterial* uniform material, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  float eta = 0.0f;
  Medium mediumOutside = make_Medium(make_Vec3f((varying Vec3fa)material->transmissionOutside),material->etaOutside);
  Medium mediumInside  = make_Medium(make_Vec3f((varying Vec3fa)material->transmissionInside ),material->etaInside );
  Medium mediumFront, mediumBack;
  if (eq(medium,mediumInside)) {
    eta = material->etaInside/material->etaOutside;
    mediumFront = mediumInside;
    mediumBack = mediumOutside;
  }
  else {
    eta = material->etaOutside/material->etaInside;
    mediumFront = mediumOutside;
    mediumBack = mediumInside;
  }

  float cosThetaO = clamp(dot(wo,dg.Ns));
  float cosThetaI; Sample3f wit = refract(wo,dg.Ns,eta,cosThetaO,cosThetaI);
  Sample3f wis = reflect_(wo,dg.Ns);
  float R = fresnelDielectric(cosThetaO,cosThetaI,eta);
  Vec3f cs = make_Vec3f(R);
  Vec3f ct = make_Vec3f(1.0f-R);
  return sample_component2(cs,wis,mediumFront,ct,wit,mediumBack,Lw,wi_o,medium,s.x);
}

////////////////////////////////////////////////////////////////////////////////
//                          ThinDielectric Material                               //
////////////////////////////////////////////////////////////////////////////////

__device__ void ThinDielectricMaterial__preprocess(uniform ThinDielectricMaterial* uniform This, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}
 
__device__ Vec3f ThinDielectricMaterial__eval(ThinDielectricMaterial* uniform This, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) {
  return make_Vec3f(0.0f);
}

__device__ Vec3f ThinDielectricMaterial__sample(ThinDielectricMaterial* uniform This, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  float cosThetaO = clamp(dot(wo,dg.Ns));
  if (cosThetaO <= 0.0f) return make_Vec3f(0.0f);
  float R = fresnelDielectric(cosThetaO,rcp(This->eta));
  Sample3f wit = make_Sample3f(neg(wo),1.0f);
  Sample3f wis = reflect_(wo,dg.Ns);
  Vec3f ct = exp(make_Vec3f(This->transmission)*rcp(cosThetaO))*make_Vec3f(1.0f-R);
  Vec3f cs = make_Vec3f(R);
  return sample_component2(cs,wis,medium,ct,wit,medium,Lw,wi_o,medium,s.x);
}

////////////////////////////////////////////////////////////////////////////////
//                     MetallicPaint Material                                 //
////////////////////////////////////////////////////////////////////////////////

__device__ void MetallicPaintMaterial__preprocess(uniform MetallicPaintMaterial* uniform material, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

__device__ Vec3f MetallicPaintMaterial__eval(MetallicPaintMaterial* uniform This, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi) 
{
  DielectricReflection reflection; DielectricReflection__Constructor(&reflection, 1.0f, This->eta);
  DielectricLayerLambertian lambertian; DielectricLayerLambertian__Constructor(&lambertian, make_Vec3f((varying float)1.0f), 1.0f, This->eta, make_Lambertian(make_Vec3f((varying Vec3fa)This->shadeColor)));
  return DielectricReflection__eval(&reflection,wo,dg,wi) + DielectricLayerLambertian__eval(&lambertian,wo,dg,wi);
}

__device__ Vec3f MetallicPaintMaterial__sample(MetallicPaintMaterial* uniform This, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  DielectricReflection reflection; DielectricReflection__Constructor(&reflection, 1.0f, This->eta);
  DielectricLayerLambertian lambertian; DielectricLayerLambertian__Constructor(&lambertian, make_Vec3f((varying float)1.0f), 1.0f, This->eta, make_Lambertian(make_Vec3f((varying Vec3fa)This->shadeColor)));
  Sample3f wi0; Vec3f c0 = DielectricReflection__sample(&reflection,wo,dg,wi0,s);
  Sample3f wi1; Vec3f c1 = DielectricLayerLambertian__sample(&lambertian,wo,dg,wi1,s);
  return sample_component2(c0,wi0,medium,c1,wi1,medium,Lw,wi_o,medium,s.x);
}

////////////////////////////////////////////////////////////////////////////////
//                              Material                                      //
////////////////////////////////////////////////////////////////////////////////

inline void Material__preprocess(ISPCMaterial* uniform materials, int materialID, uniform int numMaterials, BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
  foreach_unique (id in materialID)
  {
    if (id < 0 || id >= numMaterials) continue;
  ISPCMaterial* uniform material = &materials[id];
  switch (material->ty) {
  case MATERIAL_OBJ  : OBJMaterial__preprocess  ((uniform OBJMaterial* uniform)  material,brdf,wo,dg,medium); break;
  case MATERIAL_METAL: MetalMaterial__preprocess((uniform MetalMaterial* uniform)material,brdf,wo,dg,medium); break;
  case MATERIAL_REFLECTIVE_METAL: ReflectiveMetalMaterial__preprocess((uniform ReflectiveMetalMaterial* uniform)material,brdf,wo,dg,medium); break;
  case MATERIAL_VELVET: VelvetMaterial__preprocess((uniform VelvetMaterial* uniform)material,brdf,wo,dg,medium); break;
  case MATERIAL_DIELECTRIC: DielectricMaterial__preprocess((uniform DielectricMaterial* uniform)material,brdf,wo,dg,medium); break;
  case MATERIAL_METALLIC_PAINT: MetallicPaintMaterial__preprocess((uniform MetallicPaintMaterial* uniform)material,brdf,wo,dg,medium); break;
  case MATERIAL_MATTE: MatteMaterial__preprocess((uniform MatteMaterial* uniform)material,brdf,wo,dg,medium); break;
  case MATERIAL_MIRROR: MirrorMaterial__preprocess((uniform MirrorMaterial* uniform)material,brdf,wo,dg,medium); break;
  case MATERIAL_THIN_DIELECTRIC: ThinDielectricMaterial__preprocess((uniform ThinDielectricMaterial* uniform)material,brdf,wo,dg,medium); break;
  default: break;
  }
  }
}

inline Vec3f Material__eval(ISPCMaterial* uniform materials, int materialID, uniform int numMaterials, const BRDF& brdf, const Vec3f& wo, const DifferentialGeometry& dg, const Vec3f& wi)
{
  Vec3f c = make_Vec3f(0.0f);
  foreach_unique (id in materialID)
  {
    if (id < 0 || id >= numMaterials) continue;
  ISPCMaterial* uniform material = &materials[id];
  switch (material->ty) {
  case MATERIAL_OBJ  : c = OBJMaterial__eval  ((uniform OBJMaterial* uniform)  material, brdf, wo, dg, wi); break;
  case MATERIAL_METAL: c = MetalMaterial__eval((uniform MetalMaterial* uniform)material, brdf, wo, dg, wi); break;
  case MATERIAL_REFLECTIVE_METAL: c = ReflectiveMetalMaterial__eval((uniform ReflectiveMetalMaterial* uniform)material, brdf, wo, dg, wi); break;
  case MATERIAL_VELVET: c = VelvetMaterial__eval((uniform VelvetMaterial* uniform)material, brdf, wo, dg, wi); break;
  case MATERIAL_DIELECTRIC: c = DielectricMaterial__eval((uniform DielectricMaterial* uniform)material, brdf, wo, dg, wi); break;
  case MATERIAL_METALLIC_PAINT: c = MetallicPaintMaterial__eval((uniform MetallicPaintMaterial* uniform)material, brdf, wo, dg, wi); break;
  case MATERIAL_MATTE: c = MatteMaterial__eval((uniform MatteMaterial* uniform)material, brdf, wo, dg, wi); break;
  case MATERIAL_MIRROR: c = MirrorMaterial__eval((uniform MirrorMaterial* uniform)material, brdf, wo, dg, wi); break;
  case MATERIAL_THIN_DIELECTRIC: c = ThinDielectricMaterial__eval((uniform ThinDielectricMaterial* uniform)material, brdf, wo, dg, wi); break;
  default: c = make_Vec3f(0.0f); 
  }
  }
  return c;
}

inline Vec3f Material__sample(ISPCMaterial* uniform uniform materials, int materialID, uniform int numMaterials, const BRDF& brdf, const Vec3f& Lw, const Vec3f& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{  
  Vec3f c = make_Vec3f(0.0f);
  foreach_unique (id in materialID)
  {
    if (id < 0 || id >= numMaterials) continue;
  ISPCMaterial* uniform material = &materials[id];
  switch (material->ty) {
  case MATERIAL_OBJ  : c = OBJMaterial__sample  ((uniform OBJMaterial* uniform)  material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_METAL: c = MetalMaterial__sample((uniform MetalMaterial* uniform)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_REFLECTIVE_METAL: c = ReflectiveMetalMaterial__sample((uniform ReflectiveMetalMaterial* uniform)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_VELVET: c = VelvetMaterial__sample((uniform VelvetMaterial* uniform)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_DIELECTRIC: c = DielectricMaterial__sample((uniform DielectricMaterial* uniform)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_METALLIC_PAINT: c = MetallicPaintMaterial__sample((uniform MetallicPaintMaterial* uniform)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_MATTE: c = MatteMaterial__sample((uniform MatteMaterial* uniform)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_MIRROR: c = MirrorMaterial__sample((uniform MirrorMaterial* uniform)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_THIN_DIELECTRIC: c = ThinDielectricMaterial__sample((uniform ThinDielectricMaterial* uniform)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  default: c = make_Vec3f(0.0f); 
  }
  }
  return c;
}

////////////////////////////////////////////////////////////////////////////////
//                               Scene                                        //
////////////////////////////////////////////////////////////////////////////////

/* scene data */
extern uniform ISPCScene* uniform g_ispc_scene;
RTCScene g_scene = NULL;
void* uniform* uniform geomID_to_mesh = NULL;
uniform int* uniform geomID_to_type = NULL;

/* render function to use */
renderPixelFunc renderPixel;

/* occlusion filter function */
void occlusionFilterReject(void* uniform ptr, RTCRay& ray) {
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
}

/* rtcCommitThread called by all ISPC worker threads to enable parallel build */
#if defined(PARALLEL_COMMIT)
task void parallelCommit(RTCScene scene) {
  rtcCommitThread (scene,threadIndex,threadCount); 
}
#endif

/* error reporting function */
void error_handler(const uniform RTCError code, const uniform int8* uniform str)
{
  print("Embree: ");
  switch (code) {
  case RTC_UNKNOWN_ERROR    : print("RTC_UNKNOWN_ERROR"); break;
  case RTC_INVALID_ARGUMENT : print("RTC_INVALID_ARGUMENT"); break;
  case RTC_INVALID_OPERATION: print("RTC_INVALID_OPERATION"); break;
  case RTC_OUT_OF_MEMORY    : print("RTC_OUT_OF_MEMORY"); break;
  case RTC_UNSUPPORTED_CPU  : print("RTC_UNSUPPORTED_CPU"); break;
  default                   : print("invalid error code"); break;
  }
  if (str) { 
    print(" ("); 
    while (*str) putchar(*str++); 
    print(")\n"); 
  }
  abort();
} // error handler

/* accumulation buffer */
uniform Vec3fa* uniform g_accu = NULL;
uniform unsigned int g_accu_width = 0;
uniform unsigned int g_accu_height = 0;
uniform unsigned int g_accu_count = 0;
uniform Vec3f g_accu_vx;
uniform Vec3f g_accu_vy;
uniform Vec3f g_accu_vz;
uniform Vec3f g_accu_p;
extern uniform bool g_changed;

/* called by the C++ code for initialization */
export void device_init (uniform int8* uniform cfg)
{
  /* initialize last seen camera */
  g_accu_vx = make_Vec3f(0.0f);
  g_accu_vy = make_Vec3f(0.0f);
  g_accu_vz = make_Vec3f(0.0f);
  g_accu_p  = make_Vec3f(0.0f);

  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* set start render mode */
  renderPixel = renderPixelStandard;
  //  renderPixel = renderPixelEyeLight;

} // device_init

void convertTriangleMeshes(uniform ISPCScene* uniform scene_in, RTCScene scene_out, uniform size_t numGeometries)
{
  /* add all meshes to the scene */
  for (uniform int i=0; i<scene_in->numMeshes; i++)
  {
    /* get ith mesh */
    uniform ISPCMesh* uniform mesh = scene_in->meshes[i];

    /* create a triangle mesh */
    uniform unsigned int geomID = rtcNewTriangleMesh (scene_out, RTC_GEOMETRY_STATIC, mesh->numTriangles, mesh->numVertices);
    assert(geomID < numGeometries);
    geomID_to_mesh[geomID] = mesh;
    geomID_to_type[geomID] = 0;
    
    /* set vertices */
    uniform Vertex* uniform vertices = (uniform Vertex* uniform) rtcMapBuffer(scene_out,geomID,RTC_VERTEX_BUFFER); 
    for (uniform int j=0; j<mesh->numVertices; j++) {
      vertices[j].x = mesh->positions[j].x;
      vertices[j].y = mesh->positions[j].y;
      vertices[j].z = mesh->positions[j].z;
    }
    rtcUnmapBuffer(scene_out,geomID,RTC_VERTEX_BUFFER); 

    /* set triangles */
    uniform Triangle* uniform triangles = (uniform Triangle* uniform) rtcMapBuffer(scene_out,geomID,RTC_INDEX_BUFFER);
    for (uniform int j=0; j<mesh->numTriangles; j++) {
      triangles[j].v0 = mesh->triangles[j].v0;
      triangles[j].v1 = mesh->triangles[j].v1;
      triangles[j].v2 = mesh->triangles[j].v2;
    }
    rtcUnmapBuffer(scene_out,geomID,RTC_INDEX_BUFFER);

    uniform bool allOpaque = true;
    uniform bool allTransparent = true;
    for (uniform size_t j=0; j<mesh->numTriangles; j++) {
      uniform ISPCTriangle triangle = mesh->triangles[j];
      if (g_ispc_scene->materials[triangle.materialID].ty == MATERIAL_DIELECTRIC ||
	  g_ispc_scene->materials[triangle.materialID].ty == MATERIAL_THIN_DIELECTRIC)
	allOpaque = false;
      else 
	allTransparent = false;
    }
    if (allTransparent)
      rtcSetOcclusionFilterFunction(scene_out,geomID,(RTCFilterFuncVarying)&occlusionFilterReject);
  }
}

void convertSubdivMeshes(uniform ISPCScene* uniform scene_in, RTCScene scene_out, uniform size_t numGeometries)
{
  for (uniform size_t i=0; i<g_ispc_scene->numSubdivMeshes; i++)
  {
    uniform ISPCSubdivMesh* uniform mesh = g_ispc_scene->subdiv[i];
    uniform unsigned int geomID = rtcNewSubdivisionMesh(scene_out, RTC_GEOMETRY_STATIC, mesh->numFaces, mesh->numEdges, mesh->numVertices, 
						mesh->numEdgeCreases, mesh->numVertexCreases, mesh->numHoles);
    assert(geomID < numGeometries);
    geomID_to_mesh[geomID] = mesh;
    geomID_to_type[geomID] = 1;

    for (size_t i=0; i<mesh->numEdges; i++) mesh->subdivlevel[i] = 16;
    rtcSetBuffer(scene_out, geomID, RTC_VERTEX_BUFFER, mesh->positions, 0, sizeof(uniform Vec3fa  ));
    rtcSetBuffer(scene_out, geomID, RTC_LEVEL_BUFFER,  mesh->subdivlevel, 0, sizeof(uniform float));
    rtcSetBuffer(scene_out, geomID, RTC_INDEX_BUFFER,  mesh->position_indices  , 0, sizeof(uniform unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_FACE_BUFFER,   mesh->verticesPerFace, 0, sizeof(uniform unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_HOLE_BUFFER,   mesh->holes, 0, sizeof(uniform unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_EDGE_CREASE_INDEX_BUFFER,    mesh->edge_creases,          0, 2*sizeof(uniform unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_EDGE_CREASE_WEIGHT_BUFFER,   mesh->edge_crease_weights,   0, sizeof(uniform float));
    rtcSetBuffer(scene_out, geomID, RTC_VERTEX_CREASE_INDEX_BUFFER,  mesh->vertex_creases,        0, sizeof(uniform unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_VERTEX_CREASE_WEIGHT_BUFFER, mesh->vertex_crease_weights, 0, sizeof(uniform float));
  }
}      

typedef void* uniform void_ptr;

RTCScene convertScene(uniform ISPCScene* uniform scene_in)
{
  uniform size_t numGeometries = scene_in->numMeshes + scene_in->numSubdivMeshes;
  geomID_to_mesh = uniform new void_ptr[numGeometries];
  geomID_to_type = uniform new uniform int[numGeometries];

  /* create scene */
  RTCScene scene_out = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT, RTC_INTERSECT_VARYING);
  convertTriangleMeshes(scene_in,scene_out,numGeometries);
  convertSubdivMeshes(scene_in,scene_out,numGeometries);

  /* commit changes to scene */
#if !defined(PARALLEL_COMMIT)
  rtcCommit (scene_out);
#else
  launch[ getNumHWThreads() ] parallelCommit(scene_out); sync;
#endif

  return scene_out;
} // convertScene

/* for details about this random number generator see: P. L'Ecuyer,
   "Maximally Equidistributed Combined Tausworthe Generators",
   Mathematics of Computation, 65, 213 (1996), 203--213:
   http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps */

struct rand_state {
  unsigned int s1, s2, s3;
};

inline unsigned int irand(rand_state& state)
{
  state.s1 = ((state.s1 & 4294967294U) << 12U) ^ (((state.s1<<13U)^state.s1)>>19U);
  state.s2 = ((state.s2 & 4294967288U) <<  4U) ^ (((state.s2<< 2U)^state.s2)>>25U);
  state.s3 = ((state.s3 & 4294967280U) << 17U) ^ (((state.s3<< 3U)^state.s3)>>11U);
  return state.s1 ^ state.s2 ^ state.s3;
}

inline void init_rand(rand_state& state, unsigned int x, unsigned int y, unsigned int z)
{
  state.s1 = x >=   2 ? x : x +   2;
  state.s2 = y >=   8 ? y : y +   8;
  state.s3 = z >=  16 ? z : z +  16;
  for (uniform int i=0; i<10; i++) irand(state);
}

inline float frand(rand_state& state) {
  return irand(state)*2.3283064365386963e-10f;
}

inline Vec3f face_forward(const Vec3f& dir, const Vec3f& _Ng) {
  const Vec3f Ng = _Ng;
  return dot(dir,Ng) < 0.0f ? Ng : neg(Ng);
}

#if 0
inline Vec3f interpolate_normal(RTCRay& ray)
{
#if 0 // FIXME: pointer gather not implemented on ISPC for Xeon Phi
  uniform ISPCMesh* varying mesh = g_ispc_scene->meshes[ray.geomID];
  uniform ISPCTriangle* varying tri = &mesh->triangles[ray.primID];

  /* load material ID */
  int materialID = tri->materialID;

  /* interpolate shading normal */
  if (mesh->normals) {
    Vec3f n0 = make_Vec3f(mesh->normals[tri->v0]);
    Vec3f n1 = make_Vec3f(mesh->normals[tri->v1]);
    Vec3f n2 = make_Vec3f(mesh->normals[tri->v2]);
    float u = ray.u, v = ray.v, w = 1.0f-ray.u-ray.v;
    return normalize(w*n0 + u*n1 + v*n2);
  } else {
    return normalize(ray.Ng);
  }

#else

  Vec3f Ns = make_Vec3f(0.0f);
  int materialID = 0;
  foreach_unique (geomID in ray.geomID) 
  {
    if (geomID >= 0 && geomID < g_ispc_scene->numMeshes)  { // FIXME: workaround for ISPC bug

    uniform ISPCMesh* uniform mesh = g_ispc_scene->meshes[geomID];
    
    foreach_unique (primID in ray.primID) 
    {
      uniform ISPCTriangle* uniform tri = &mesh->triangles[primID];
      
      /* load material ID */
      materialID = tri->materialID;

      /* interpolate shading normal */
      if (mesh->normals) {
        Vec3f n0 = make_Vec3f(mesh->normals[tri->v0]);
        Vec3f n1 = make_Vec3f(mesh->normals[tri->v1]);
        Vec3f n2 = make_Vec3f(mesh->normals[tri->v2]);
        float u = ray.u, v = ray.v, w = 1.0f-ray.u-ray.v;
        Ns = w*n0 + u*n1 + v*n2;
      } else {
        Ns = normalize(ray.Ng);
      }
    }
    }
  }
  return normalize(Ns);
#endif
}
#endif

Vec3f renderPixelFunction(float x, float y, rand_state& state, const uniform Vec3f& vx, const uniform Vec3f& vy, const uniform Vec3f& vz, const uniform Vec3f& p)
{
  /* radiance accumulator and weight */
  Vec3f L = make_Vec3f(0.0f);
  Vec3f Lw = make_Vec3f(1.0f);
  Medium medium = make_Medium_Vacuum();

  /* initialize ray */
  RTCRay ray = make_Ray(p,normalize(x*vx + y*vy + vz),0.0f,inf);

  /* iterative path tracer loop */
  for (uniform int i=0; i<8; i++)
  {
    /* terminate if contribution too low */
    if (max(Lw.x,max(Lw.y,Lw.z)) < 0.01f)
      break;

    /* intersect ray with scene */ 
    rtcIntersect(g_scene,ray);
    const Vec3f wo = neg(ray.dir);
    
    /* invoke environment lights if nothing hit */
    if (ray.geomID == RTC_INVALID_GEOMETRY_ID) 
    {
#if 0
      /* iterate over all ambient lights */
      for (uniform size_t i=0; i<g_ispc_scene->numAmbientLights; i++)
        L = L + Lw*AmbientLight__eval(g_ispc_scene->ambientLights[i],ray.dir); // FIXME: +=
#endif

#if 0
      /* iterate over all distant lights */
      for (uniform size_t i=0; i<g_ispc_scene->numDistantLights; i++)
        L = L + Lw*DistantLight__eval(g_ispc_scene->distantLights[i],ray.dir); // FIXME: +=
#endif

      break;
    }

    /* compute differential geometry */
    DifferentialGeometry dg;
    dg.P  = ray.org+ray.tfar*ray.dir;
    dg.Ng = face_forward(ray.dir,normalize(ray.Ng));
    //Vec3f _Ns = interpolate_normal(ray);
    Vec3f _Ns = normalize(ray.Ng);
    dg.Ns = face_forward(ray.dir,_Ns);

    /* shade all rays that hit something */
#if 0 // FIXME: pointer gather not implemented in ISPC for Xeon Phi
    int materialID = 0;
    if (geomID_to_type[ray.geomID] == 0)
      materialID = ((ISPCMesh* uniform) geomID_to_mesh[ray.geomID])->triangles[ray.primID].materialID; 
    else 
      materialID = ((ISPCSubdivMesh* uniform) geomID_to_mesh[ray.geomID])->materialID; 
#else
    int materialID = 0;
    foreach_unique (geomID in ray.geomID) {
      if (geomID >= 0 && geomID < g_ispc_scene->numMeshes) { // FIXME: workaround for ISPC bug
	if (geomID_to_type[geomID] == 0) 
	  materialID = ((ISPCMesh* uniform) geomID_to_mesh[geomID])->triangles[ray.primID].materialID; 
	else                             
	  materialID = ((ISPCSubdivMesh* uniform) geomID_to_mesh[geomID])->materialID; 
      }
    }
#endif

    /*! Compute  simple volumetric effect. */
    Vec3f c = make_Vec3f(1.0f);
    const Vec3f transmission = medium.transmission;
    if (ne(transmission,make_Vec3f(1.0f)))
      c = c * pow(transmission,ray.tfar);
    
    /* calculate BRDF */ // FIXME: avoid gathers
    BRDF brdf;
    uniform int numMaterials = g_ispc_scene->numMaterials;
    //uniform ISPCMaterial* material = &g_ispc_scene->materials[materialID];
    uniform ISPCMaterial* uniform material_array = &g_ispc_scene->materials[0];
    Material__preprocess(material_array,materialID,numMaterials,brdf,wo,dg,medium);

    /* sample BRDF at hit point */
    Sample3f wi1;
    c = c * Material__sample(material_array,materialID,numMaterials,brdf,Lw, wo, dg, wi1, medium, make_Vec2f(frand(state),frand(state)));

    /* iterate over ambient lights */
    for (uniform size_t i=0; i<g_ispc_scene->numAmbientLights; i++)
    {
#if 1
      Vec3f L0 = make_Vec3f(0.0f);
      Sample3f wi0; float tMax0;
      Vec3f Ll0 = AmbientLight__sample(g_ispc_scene->ambientLights[i],dg,wi0,tMax0,make_Vec2f(frand(state),frand(state)));
      if (wi0.pdf > 0.0f) {
        RTCRay shadow = make_Ray(dg.P,wi0.v,0.001f,tMax0);
        rtcOccluded(g_scene,shadow);
        if (shadow.geomID == RTC_INVALID_GEOMETRY_ID) {
          L0 = Ll0/wi0.pdf*Material__eval(material_array,materialID,numMaterials,brdf,wo,dg,wi0.v);
        }
        L = L + Lw*L0;
      }
#endif

#if 0
      Vec3f L1 = make_Vec3f(0.0f);
      Vec3f Ll1 = AmbientLight__eval(g_ispc_scene->ambientLights[i],wi1.v);
      if (wi1.pdf > 0.0f) {
        RTCRay shadow = make_Ray(dg.P,wi1.v,0.001f,inf);
        rtcOccluded(g_scene,shadow);
        if (shadow.geomID == RTC_INVALID_GEOMETRY_ID) {
          L1 = Ll1/wi1.pdf*c;
        }
        L = L + Lw*L1;
      }
#endif

#if 0
      float s = wi0.pdf*wi0.pdf + wi1.pdf*wi1.pdf;
      if (s > 0) {
        float w0 = 0;
        float w1 = 1;
        //float w0 = wi0.pdf*wi0.pdf/s;
        //float w1 = wi1.pdf*wi1.pdf/s;
        L = L + Lw*(w0*L0+w1*L1);
      }
#endif
    }

    Sample3f wi; float tMax;

    /* iterate over point lights */
    for (uniform size_t i=0; i<g_ispc_scene->numPointLights; i++)
    {
      Vec3f Ll = PointLight__sample(g_ispc_scene->pointLights[i],dg,wi,tMax,make_Vec2f(frand(state),frand(state)));
      if (wi.pdf <= 0.0f) continue;
      RTCRay shadow = make_Ray(dg.P,wi.v,0.001f,tMax);
      rtcOccluded(g_scene,shadow);
      if (shadow.geomID != RTC_INVALID_GEOMETRY_ID) continue;
      L = L + Lw*Ll/wi.pdf*Material__eval(material_array,materialID,numMaterials,brdf,wo,dg,wi.v); // FIXME: +=
    }

    /* iterate over directional lights */
    for (uniform size_t i=0; i<g_ispc_scene->numDirectionalLights; i++)
    {
      Vec3f Ll = DirectionalLight__sample(g_ispc_scene->dirLights[i],dg,wi,tMax,make_Vec2f(frand(state),frand(state)));
      if (wi.pdf <= 0.0f) continue;
      RTCRay shadow = make_Ray(dg.P,wi.v,0.001f,tMax);
      rtcOccluded(g_scene,shadow);
      if (shadow.geomID != RTC_INVALID_GEOMETRY_ID) continue;
      L = L + Lw*Ll/wi.pdf*Material__eval(material_array,materialID,numMaterials,brdf,wo,dg,wi.v); // FIXME: +=
    }

    /* iterate over distant lights */
    for (uniform size_t i=0; i<g_ispc_scene->numDistantLights; i++)
    {
      Vec3f Ll = DistantLight__sample(g_ispc_scene->distantLights[i],dg,wi,tMax,make_Vec2f(frand(state),frand(state)));
      if (wi.pdf <= 0.0f) continue;
      RTCRay shadow = make_Ray(dg.P,wi.v,0.001f,tMax);
      rtcOccluded(g_scene,shadow);
      if (shadow.geomID != RTC_INVALID_GEOMETRY_ID) continue;
      L = L + Lw*Ll/wi.pdf*Material__eval(material_array,materialID,numMaterials,brdf,wo,dg,wi.v); // FIXME: +=
    }

    if (wi1.pdf <= 0.0f) break;
    Lw = Lw*c/wi1.pdf; // FIXME: *=

    /* setup secondary ray */
    ray = make_Ray(dg.P,normalize(wi1.v),0.001f,inf);
  }
  return L;
}

/* task that renders a single screen tile */
Vec3f renderPixelStandard(float x, float y, const uniform Vec3f& vx, const uniform Vec3f& vy, const uniform Vec3f& vz, const uniform Vec3f& p)
{
  rand_state state;
  init_rand(state,
            253*x+35*y+152*g_accu_count+54,
            1253*x+345*y+1452*g_accu_count+564,
            10253*x+3435*y+52*g_accu_count+13);

  Vec3f L = make_Vec3f(0.0f,0.0f,0.0f);
  //for (uniform int i=0; i<16; i++) {
  L = L + renderPixelFunction(x,y,state,vx,vy,vz,p); // FIXME: +=
  //}
  //L = L*(1.0f/16.0f);
  return L;
}
  
/* task that renders a single screen tile */
task void renderTile(uniform int* uniform pixels,
                     const uniform int width,
                     const uniform int height, 
                     const uniform float time,
                     const uniform Vec3f& vx, 
                     const uniform Vec3f& vy, 
                     const uniform Vec3f& vz, 
                     const uniform Vec3f& p,
                     const uniform int numTilesX, 
                     const uniform int numTilesY)
{
  const uniform int tileY = taskIndex / numTilesX;
  const uniform int tileX = taskIndex - tileY * numTilesX;
  const uniform int x0 = tileX * TILE_SIZE_X;
  const uniform int x1 = min(x0+TILE_SIZE_X,width);
  const uniform int y0 = tileY * TILE_SIZE_Y;
  const uniform int y1 = min(y0+TILE_SIZE_Y,height);

  foreach (y = y0 ... y1, x = x0 ... x1)
  {
    //if (x != 200 || y != 450) continue;

    /* calculate pixel color */
    Vec3f color = renderPixel(x,y,vx,vy,vz,p);

    /* write color to framebuffer */
    Vec3fa* dst = &g_accu[y*width+x];
    *dst = *dst + make_Vec3fa(color.x,color.y,color.z,1.0f); // FIXME: use += operator
    float f = rcp(max(0.001f,dst->w));
    unsigned int r = (unsigned int) (255.0f * clamp(dst->x*f,0.0f,1.0f));
    unsigned int g = (unsigned int) (255.0f * clamp(dst->y*f,0.0f,1.0f));
    unsigned int b = (unsigned int) (255.0f * clamp(dst->z*f,0.0f,1.0f));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
} // renderTile

/* called by the C++ code to render */
export void device_render (uniform int* uniform pixels,
                           const uniform int width,
                           const uniform int height, 
                           const uniform float time,
                           const uniform Vec3f& vx, 
                           const uniform Vec3f& vy, 
                           const uniform Vec3f& vz, 
                           const uniform Vec3f& p)
{
  /* create scene */
  if (g_scene == NULL)
    g_scene = convertScene(g_ispc_scene);

  /* create accumulator */
  if (g_accu_width != width || g_accu_height != height) {
    delete[] g_accu;
    g_accu = uniform new uniform Vec3fa[width*height];
    g_accu_width = width;
    g_accu_height = height;
    memset(g_accu,0,width*height*sizeof(uniform Vec3fa));
  }

  /* reset accumulator */
  uniform bool camera_changed = g_changed; g_changed = false;
  camera_changed |= ne(g_accu_vx,vx); g_accu_vx = vx; // FIXME: use != operator
  camera_changed |= ne(g_accu_vy,vy); g_accu_vy = vy; // FIXME: use != operator
  camera_changed |= ne(g_accu_vz,vz); g_accu_vz = vz; // FIXME: use != operator
  camera_changed |= ne(g_accu_p,  p); g_accu_p  = p;  // FIXME: use != operator
  g_accu_count++;
  if (camera_changed) {
    g_accu_count=0;
    memset(g_accu,0,width*height*sizeof(uniform Vec3fa));
  }

  /* render image */
  const uniform int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const uniform int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch[numTilesX*numTilesY] renderTile(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); sync;
  rtcDebug();
} // device_render

/* called by the C++ code for cleanup */
export void device_cleanup ()
{
  delete[] g_accu;
  rtcDeleteScene (g_scene);
  rtcExit();
} // device_cleanup
