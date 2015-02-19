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

#include "../common/tutorial/tutorial_device.h"
#include "../common/tutorial/scene_device.h"
#include "shapesampler.h"
#include "optics.h"





struct DifferentialGeometry
{
  Vec3fa P;
  Vec3fa Ng;
  Vec3fa Ns;
};

struct BRDF
{
  float Ns;               /*< specular exponent */
  float Ni;               /*< optical density for the surface (index of refraction) */
  Vec3fa Ka;              /*< ambient reflectivity */
  Vec3fa Kd;              /*< diffuse reflectivity */
  Vec3fa Ks;              /*< specular reflectivity */
  Vec3fa Kt;              /*< transmission filter */
};

struct Medium
{
  Vec3fa transmission; //!< Transmissivity of medium.
  float eta;             //!< Refraction index of medium.
};

inline Medium make_Medium(const Vec3fa& transmission, const float eta)
{
  Medium m;
  m.transmission = transmission;
  m.eta = eta;
  return m;
}

inline Medium make_Medium_Vacuum() { 
  return make_Medium(Vec3fa((float)1.0f),1.0f); 
}

inline bool eq(const Medium& a, const Medium& b) {
  return (a.eta == b.eta) && eq(a.transmission, b.transmission);
}

inline Vec3fa sample_component2(const Vec3fa& c0, const Sample3f& wi0, const Medium& medium0,
                               const Vec3fa& c1, const Sample3f& wi1, const Medium& medium1,
                               const Vec3fa& Lw, Sample3f& wi_o, Medium& medium_o, const float s)
{
  const Vec3fa m0 = Lw*c0/wi0.pdf;
  const Vec3fa m1 = Lw*c1/wi1.pdf;

  const float C0 = wi0.pdf == 0.0f ? 0.0f : max(max(m0.x,m0.y),m0.z);
  const float C1 = wi1.pdf == 0.0f ? 0.0f : max(max(m1.x,m1.y),m1.z);
  const float C  = C0 + C1;

  if (C == 0.0f) {
    wi_o = Sample3f(Vec3fa(0,0,0),0);
    return Vec3fa(0,0,0);
  }

  const float CP0 = C0/C;
  const float CP1 = C1/C;
  if (s < CP0) {
    wi_o = Sample3f(wi0.v,wi0.pdf*CP0); 
    medium_o = medium0; return c0;
  } 
  else {
    wi_o = Sample3f(wi1.v,wi1.pdf*CP1); 
    medium_o = medium1; return c1;
  }
}

////////////////////////////////////////////////////////////////////////////////
//                             Ambient Light                                  //
////////////////////////////////////////////////////////////////////////////////

inline Vec3fa AmbientLight__eval(const ISPCAmbientLight& light, const Vec3fa& wo) {
  return Vec3fa(light.L);
}

inline Vec3fa AmbientLight__sample(const ISPCAmbientLight& light, const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) 
{
  wi = cosineSampleHemisphere(s.x,s.y,dg.Ns);
  tMax = 1e20f;
  return Vec3fa(light.L);
}

////////////////////////////////////////////////////////////////////////////////
//                             Point Light                                    //
////////////////////////////////////////////////////////////////////////////////

inline Vec3fa PointLight__sample(const ISPCPointLight& light, 
					const DifferentialGeometry& dg, 
					Sample3f& wi,
					float& tMax,
					const Vec2f& s) 
{
  Vec3fa d = Vec3fa(light.P) - dg.P;
  float distance = length(d);
  wi = Sample3f(d*rcp(distance), distance*distance);
  tMax = distance;
  return Vec3fa(light.I);
}

////////////////////////////////////////////////////////////////////////////////
//                        Directional Light                                   //
////////////////////////////////////////////////////////////////////////////////

inline Vec3fa DirectionalLight__sample(const ISPCDirectionalLight& light, 
					      const DifferentialGeometry& dg, 
					      Sample3f& wi,
					      float& tMax,
					      const Vec2f& s) 
{
  wi = Sample3f(neg(normalize(Vec3fa(light.D))),1.0f); 
  tMax = inf; 
  return Vec3fa(light.E);
}

////////////////////////////////////////////////////////////////////////////////
//                          Distant Light                                     //
////////////////////////////////////////////////////////////////////////////////

inline Vec3fa DistantLight__eval(const ISPCDistantLight& light, const Vec3fa& wo) 
{
  if (-dot(wo,Vec3fa(light.D)) >= light.cosHalfAngle) return Vec3fa(light.L);
  return Vec3fa(0.0f);
}

inline Vec3fa DistantLight__sample(const ISPCDistantLight& light,
                                   const DifferentialGeometry& dg, 
                                   Sample3f& wi,
                                   float& tMax,
                                   const Vec2f& s) 
{
  wi = UniformSampleCone(s.x,s.y,light.radHalfAngle,Vec3fa((Vec3fa)neg(light.D)));
  tMax = 1e20f;
  return Vec3fa(light.L);
}

////////////////////////////////////////////////////////////////////////////////
//                          Minneart BRDF                                     //
////////////////////////////////////////////////////////////////////////////////

struct Minneart
{
  /*! The reflectance parameter. The vale 0 means no reflection,
   *  and 1 means full reflection. */
  Vec3fa R;
  
  /*! The amount of backscattering. A value of 0 means lambertian
   *  diffuse, and inf means maximum backscattering. */
  float b;
};

inline Vec3fa Minneart__eval(const Minneart* This,
                     const Vec3fa &wo, const DifferentialGeometry &dg, const Vec3fa &wi) 
{
  const float cosThetaI = clamp(dot(wi,dg.Ns));
  const float backScatter = pow(clamp(dot(wo,wi)), This->b);
  return (backScatter * cosThetaI * float(one_over_pi)) * This->R;
}

inline Vec3fa Minneart__sample(const Minneart* This,
                       const Vec3fa &wo, 
                       const DifferentialGeometry &dg, 
                       Sample3f &wi, 
                       const Vec2f &s)  
{
  wi = cosineSampleHemisphere(s.x,s.y,dg.Ns);
  return Minneart__eval(This, wo, dg, wi.v);
}

inline void Minneart__Constructor(Minneart* This, const Vec3fa& R, const float b) 
{
  This->R = R;
  This->b = b;
}

inline Minneart make_Minneart(const Vec3fa& R, const float f) { 
  Minneart m; Minneart__Constructor(&m,R,f); return m; 
}

////////////////////////////////////////////////////////////////////////////////
//                            Velvet BRDF                                     //
////////////////////////////////////////////////////////////////////////////////

struct Velvety
{
  BRDF base;

  /*! The reflectance parameter. The vale 0 means no reflection,
   *  and 1 means full reflection. */
  Vec3fa R;
  
  /*! The falloff of horizon scattering. 0 no falloff,
   *  and inf means maximum falloff. */
  float f;
};

inline Vec3fa Velvety__eval(const Velvety* This,
                    const Vec3fa &wo, const DifferentialGeometry &dg, const Vec3fa &wi) 
{
  const float cosThetaO = clamp(dot(wo,dg.Ns));
  const float cosThetaI = clamp(dot(wi,dg.Ns));
  const float sinThetaO = sqrt(1.0f - cosThetaO * cosThetaO);
  const float horizonScatter = pow(sinThetaO, This->f);
  return (horizonScatter * cosThetaI * float(one_over_pi)) * This->R;
}

inline Vec3fa Velvety__sample(const Velvety* This,
                      const Vec3fa &wo, 
                      const DifferentialGeometry &dg, 
                      Sample3f &wi, 
                      const Vec2f &s)  
{
  wi = cosineSampleHemisphere(s.x,s.y,dg.Ns);
  return Velvety__eval(This, wo, dg, wi.v);
}

inline void Velvety__Constructor(Velvety* This, const Vec3fa& R, const float f) 
{
  This->R = R;
  This->f = f;
}

inline Velvety make_Velvety(const Vec3fa& R, const float f) { 
  Velvety m; Velvety__Constructor(&m,R,f); return m; 
}

////////////////////////////////////////////////////////////////////////////////
//                  Dielectric Reflection BRDF                                //
////////////////////////////////////////////////////////////////////////////////

struct DielectricReflection
{
  float eta;
};

inline Vec3fa DielectricReflection__eval(const DielectricReflection* This, const Vec3fa &wo, const DifferentialGeometry &dg, const Vec3fa &wi) {
  return Vec3fa(0.f);
}

inline Vec3fa DielectricReflection__sample(const DielectricReflection* This, const Vec3fa &wo, const DifferentialGeometry &dg, Sample3f &wi, const Vec2f &s)
{
  const float cosThetaO = clamp(dot(wo,dg.Ns));
  wi = reflect_(wo,dg.Ns,cosThetaO);
  return Vec3fa(fresnelDielectric(cosThetaO,This->eta));
}

inline void DielectricReflection__Constructor(DielectricReflection* This,
                                              const float etai,
                                              const float etat)
{
  This->eta = etai*rcp(etat);
}

inline DielectricReflection make_DielectricReflection(const float etai, const float etat) {
  DielectricReflection v; DielectricReflection__Constructor(&v,etai,etat); return v;
}

////////////////////////////////////////////////////////////////////////////////
//                                Lambertian BRDF                             //
////////////////////////////////////////////////////////////////////////////////

struct Lambertian
{
  Vec3fa R;
};

inline Vec3fa Lambertian__eval(const Lambertian* This,
                              const Vec3fa &wo, const DifferentialGeometry &dg, const Vec3fa &wi) 
{
  return This->R * (1.0f/(float)(float(pi))) * clamp(dot(wi,dg.Ns));
}

inline Vec3fa Lambertian__sample(const Lambertian* This,
                                const Vec3fa &wo, 
                                const DifferentialGeometry &dg, 
                                Sample3f &wi, 
                                const Vec2f &s)  
{
  wi = cosineSampleHemisphere(s.x,s.y,dg.Ns);
  return Lambertian__eval(This, wo, dg, wi.v);
}

inline void Lambertian__Constructor(Lambertian* This, const Vec3fa& R)
{
  This->R = R;
}

inline Lambertian make_Lambertian(const Vec3fa& R) {
  Lambertian v; Lambertian__Constructor(&v,R); return v;
}


////////////////////////////////////////////////////////////////////////////////
//              Lambertian BRDF with Dielectric Layer on top                  //
////////////////////////////////////////////////////////////////////////////////

struct DielectricLayerLambertian
{
  Vec3fa T;             //!< Transmission coefficient of dielectricum
  float etait;         //!< Relative refraction index etai/etat of both media
  float etati;         //!< relative refraction index etat/etai of both media
  Lambertian ground;   //!< the BRDF of the ground layer
};

inline Vec3fa DielectricLayerLambertian__eval(const DielectricLayerLambertian* This,
                                             const Vec3fa &wo, const DifferentialGeometry &dg, const Vec3fa &wi) 
{
  const float cosThetaO = dot(wo,dg.Ns);
  const float cosThetaI = dot(wi,dg.Ns);
  if (cosThetaI <= 0.0f | cosThetaO <= 0.0f) return Vec3fa(0.f);

  float cosThetaO1; 
  const Sample3f wo1 = refract(wo,dg.Ns,This->etait,cosThetaO,cosThetaO1);
  float cosThetaI1; 
  const Sample3f wi1 = refract(wi,dg.Ns,This->etait,cosThetaI,cosThetaI1);
  const float Fi = 1.0f - fresnelDielectric(cosThetaI,cosThetaI1,This->etait);
  const Vec3fa Fg = Lambertian__eval(&This->ground,neg(wo1.v),dg,neg(wi1.v));
  const float Fo = 1.0f - fresnelDielectric(cosThetaO,cosThetaO1,This->etait);
  return Fo * This->T * Fg * This->T * Fi;
}

inline Vec3fa DielectricLayerLambertian__sample(const DielectricLayerLambertian* This,
                                               const Vec3fa &wo, 
                                               const DifferentialGeometry &dg, 
                                               Sample3f &wi, 
                                               const Vec2f &s)  
{
  /*! refract ray into medium */
  float cosThetaO = dot(wo,dg.Ns);
  if (cosThetaO <= 0.0f) return Vec3fa(0.f);
  float cosThetaO1; Sample3f wo1 = refract(wo,dg.Ns,This->etait,cosThetaO,cosThetaO1);
  
  /*! sample ground BRDF */
  Sample3f wi1 = Sample3f(Vec3fa(0.f),1.f); 
  Vec3fa Fg = Lambertian__sample(&This->ground,neg(wo1.v),dg,wi1,s);

  /*! refract ray out of medium */
  float cosThetaI1 = dot(wi1.v,dg.Ns);
  if (cosThetaI1 <= 0.0f) return Vec3fa(0.f);
  
  float cosThetaI; 
  Sample3f wi0 = refract(neg(wi1.v),neg(dg.Ns),This->etati,cosThetaI1,cosThetaI);
  if (wi0.pdf == 0.0f) return Vec3fa(0.f);
  
  /*! accumulate contribution of path */
  wi = Sample3f(wi0.v,wi1.pdf);
  float Fi = 1.0f - fresnelDielectric(cosThetaI,cosThetaI1,This->etait);
  float Fo = 1.0f - fresnelDielectric(cosThetaO,cosThetaO1,This->etait);
  return Fo * This->T * Fg * This->T * Fi;
}

inline void DielectricLayerLambertian__Constructor(DielectricLayerLambertian* This,
                                                   const Vec3fa& T, 
                                                   const float etai, 
                                                   const float etat, 
                                                   const Lambertian& ground)
{
  This->T = T;
  This->etait = etai*rcp(etat);
  This->etati = etat*rcp(etai);
  This->ground = ground;
}

inline DielectricLayerLambertian make_DielectricLayerLambertian(const Vec3fa& T, 
                                                                        const float etai, 
                                                                        const float etat, 
                                                                        const Lambertian& ground)
{
  DielectricLayerLambertian m; 
  DielectricLayerLambertian__Constructor(&m,T,etai,etat,ground);
  return m;
}

////////////////////////////////////////////////////////////////////////////////
//                          Matte Material                                    //
////////////////////////////////////////////////////////////////////////////////

 void MatteMaterial__preprocess(MatteMaterial* material, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

 Vec3fa MatteMaterial__eval(MatteMaterial* This, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) 
{
  Lambertian lambertian = make_Lambertian(Vec3fa((Vec3fa)This->reflectance));
  return Lambertian__eval(&lambertian,wo,dg,wi);
}

 Vec3fa MatteMaterial__sample(MatteMaterial* This, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  Lambertian lambertian = make_Lambertian(Vec3fa((Vec3fa)This->reflectance));
  return Lambertian__sample(&lambertian,wo,dg,wi_o,s);
}

////////////////////////////////////////////////////////////////////////////////
//                          Mirror Material                                    //
////////////////////////////////////////////////////////////////////////////////

 void MirrorMaterial__preprocess(MirrorMaterial* material, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

 Vec3fa MirrorMaterial__eval(MirrorMaterial* This, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) {
  return Vec3fa(0.0f);
}

 Vec3fa MirrorMaterial__sample(MirrorMaterial* This, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  wi_o = reflect_(wo,dg.Ns);
  return Vec3fa(This->reflectance);
}

////////////////////////////////////////////////////////////////////////////////
//                          OBJ Material                                      //
////////////////////////////////////////////////////////////////////////////////

 void OBJMaterial__preprocess(OBJMaterial* material, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
    float d = material->d;
    //if (material->map_d) { d *= material->map_d.get(s,t); }
    brdf.Ka = Vec3fa(material->Ka);
    //if (material->map_Ka) { brdf.Ka *= material->map_Ka->get(dg.st); }
    brdf.Kd = d * Vec3fa(material->Kd);  
    //if (material->map_Kd) brdf.Kd *= material->map_Kd->get(dg.st);  
    brdf.Ks = d * Vec3fa(material->Ks);  
    //if (material->map_Ks) brdf.Ks *= material->map_Ks->get(dg.st); 
    brdf.Ns = material->Ns;  
    //if (material->map_Ns) { brdf.Ns *= material->map_Ns.get(dg.st); }
    brdf.Kt = (1.0f-d)*Vec3fa(material->Kt);
    brdf.Ni = material->Ni;
}

 Vec3fa OBJMaterial__eval(OBJMaterial* material, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) 
{
  Vec3fa R = Vec3fa(0.0f,0.0f,0.0f);
  const float Md = max(max(brdf.Kd.x,brdf.Kd.y),brdf.Kd.z);
  const float Ms = max(max(brdf.Ks.x,brdf.Ks.y),brdf.Ks.z);
  const float Mt = max(max(brdf.Kt.x,brdf.Kt.y),brdf.Kt.z);
  if (Md > 0.0f) {
    R = R + (1.0f/float(pi)) * clamp(dot(wi,dg.Ns)) * brdf.Kd; // FIXME: +=
  }
  if (Ms > 0.0f) {
    const Sample3f refl = reflect_(wo,dg.Ns);
    if (dot(refl.v,wi) > 0.0f) 
      R = R + (brdf.Ns+2) * float(one_over_two_pi) * pow(max(1e-10f,dot(refl.v,wi)),brdf.Ns) * clamp(dot(wi,dg.Ns)) * brdf.Ks; // FIXME: +=
  }
  if (Mt > 0.0f) {
  }
  return R;
}

 Vec3fa OBJMaterial__sample(OBJMaterial* material, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  Vec3fa cd = Vec3fa(0.0f); 
  Sample3f wid = Sample3f(Vec3fa(0.0f),0.0f);
  if (max(max(brdf.Kd.x,brdf.Kd.y),brdf.Kd.z) > 0.0f) {
    wid = cosineSampleHemisphere(s.x,s.y,dg.Ns);
    cd = float(one_over_pi) * clamp(dot(wid.v,dg.Ns)) * brdf.Kd;
  }

  Vec3fa cs = Vec3fa(0.0f); 
  Sample3f wis = Sample3f(Vec3fa(0.0f),0.0f);
  if (max(max(brdf.Ks.x,brdf.Ks.y),brdf.Ks.z) > 0.0f)
  {
    const Sample3f refl = reflect_(wo,dg.Ns);
    wis = powerCosineSampleHemisphere(s.x,s.y,refl.v,brdf.Ns);
    cs = (brdf.Ns+2) * float(one_over_two_pi) * pow(dot(refl.v,wis.v),brdf.Ns) * clamp(dot(wis.v,dg.Ns)) * brdf.Ks;
  }

  Vec3fa ct = Vec3fa(0.0f); 
  Sample3f wit = Sample3f(Vec3fa(0.0f),0.0f);
  if (max(max(brdf.Kt.x,brdf.Kt.y),brdf.Kt.z) > 0.0f)
  {
    wit = Sample3f(neg(wo),1.0f);
    ct = brdf.Kt;
  }

  const Vec3fa md = Lw*cd/wid.pdf;
  const Vec3fa ms = Lw*cs/wis.pdf;
  const Vec3fa mt = Lw*ct/wit.pdf;

  const float Cd = wid.pdf == 0.0f ? 0.0f : max(max(md.x,md.y),md.z);
  const float Cs = wis.pdf == 0.0f ? 0.0f : max(max(ms.x,ms.y),ms.z);
  const float Ct = wit.pdf == 0.0f ? 0.0f : max(max(mt.x,mt.y),mt.z);
  const float C  = Cd + Cs + Ct;

  if (C == 0.0f) {
    wi_o = Sample3f(Vec3fa(0,0,0),0);
    return Vec3fa(0,0,0);
  }

  const float CPd = Cd/C;
  const float CPs = Cs/C;
  const float CPt = Ct/C;

  if (s.x < CPd) {
    wi_o = Sample3f(wid.v,wid.pdf*CPd);
    return cd;
  } 
  else if (s.x < CPd + CPs)
  {
    wi_o = Sample3f(wis.v,wis.pdf*CPs);
    return cs;
  }
  else 
  {
    wi_o = Sample3f(wit.v,wit.pdf*CPt);
    return ct;
  }
}

////////////////////////////////////////////////////////////////////////////////
//                        Metal Material                                      //
////////////////////////////////////////////////////////////////////////////////

 void MetalMaterial__preprocess(MetalMaterial* material, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

 Vec3fa MetalMaterial__eval(MetalMaterial* This, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) 
{
  const FresnelConductor fresnel = make_FresnelConductor(Vec3fa(This->eta),Vec3fa(This->k));
  const PowerCosineDistribution distribution = make_PowerCosineDistribution(rcp(This->roughness));

  const float cosThetaO = dot(wo,dg.Ns);
  const float cosThetaI = dot(wi,dg.Ns);
  if (cosThetaI <= 0.0f | cosThetaO <= 0.0f) return Vec3fa(0.f);
  const Vec3fa wh = normalize(wi+wo);
  const float cosThetaH = dot(wh, dg.Ns);
  const float cosTheta = dot(wi, wh); // = dot(wo, wh);
  const Vec3fa F = eval(fresnel,cosTheta);
  const float D = eval(distribution,cosThetaH);
  const float G = min(1.0f, min(2.0f * cosThetaH * cosThetaO / cosTheta, 
                                2.0f * cosThetaH * cosThetaI / cosTheta));
  return (Vec3fa(This->reflectance)*F) * D * G * rcp(4.0f*cosThetaO);
}

 Vec3fa MetalMaterial__sample(MetalMaterial* This, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  const PowerCosineDistribution distribution = make_PowerCosineDistribution(rcp(This->roughness));

  if (dot(wo,dg.Ns) <= 0.0f) return Vec3fa(0.0f);
  sample(distribution,wo,dg.Ns,wi_o,s);
  if (dot(wi_o.v,dg.Ns) <= 0.0f) return Vec3fa(0.0f);
  return MetalMaterial__eval(This,brdf,wo,dg,wi_o.v);
}

////////////////////////////////////////////////////////////////////////////////
//                        ReflectiveMetal Material                            //
////////////////////////////////////////////////////////////////////////////////

 void ReflectiveMetalMaterial__preprocess(ReflectiveMetalMaterial* material, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  {
}

 Vec3fa ReflectiveMetalMaterial__eval(ReflectiveMetalMaterial* This, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) {
  return Vec3fa(0.0f);
}

 Vec3fa ReflectiveMetalMaterial__sample(ReflectiveMetalMaterial* This, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  wi_o = reflect_(wo,dg.Ns);
  return Vec3fa(This->reflectance) * fresnelConductor(dot(wo,dg.Ns),Vec3fa((Vec3fa)This->eta),Vec3fa((Vec3fa)This->k));
}

////////////////////////////////////////////////////////////////////////////////
//                        Velvet Material                                     //
////////////////////////////////////////////////////////////////////////////////

 void VelvetMaterial__preprocess(VelvetMaterial* material, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

 Vec3fa VelvetMaterial__eval(VelvetMaterial* This, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) 
{
  Minneart minneart; Minneart__Constructor(&minneart,(Vec3fa)Vec3fa(This->reflectance),This->backScattering);
  Velvety velvety; Velvety__Constructor (&velvety,Vec3fa((Vec3fa)This->horizonScatteringColor),This->horizonScatteringFallOff);
  return Minneart__eval(&minneart,wo,dg,wi) + Velvety__eval(&velvety,wo,dg,wi);
}

 Vec3fa VelvetMaterial__sample(VelvetMaterial* This, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  Minneart minneart; Minneart__Constructor(&minneart,Vec3fa((Vec3fa)This->reflectance),This->backScattering);
  Velvety velvety; Velvety__Constructor (&velvety,Vec3fa((Vec3fa)This->horizonScatteringColor),This->horizonScatteringFallOff);

  Sample3f wi0; Vec3fa c0 = Minneart__sample(&minneart,wo,dg,wi0,s);
  Sample3f wi1; Vec3fa c1 = Velvety__sample(&velvety,wo,dg,wi1,s);
  return sample_component2(c0,wi0,medium,c1,wi1,medium,Lw,wi_o,medium,s.x);
}

////////////////////////////////////////////////////////////////////////////////
//                          Dielectric Material                               //
////////////////////////////////////////////////////////////////////////////////

 void DielectricMaterial__preprocess(DielectricMaterial* material, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

 Vec3fa DielectricMaterial__eval(DielectricMaterial* material, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) {
  return Vec3fa(0.0f);
}

 Vec3fa DielectricMaterial__sample(DielectricMaterial* material, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  float eta = 0.0f;
  Medium mediumOutside = make_Medium(Vec3fa((Vec3fa)material->transmissionOutside),material->etaOutside);
  Medium mediumInside  = make_Medium(Vec3fa((Vec3fa)material->transmissionInside ),material->etaInside );
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
  Vec3fa cs = Vec3fa(R);
  Vec3fa ct = Vec3fa(1.0f-R);
  return sample_component2(cs,wis,mediumFront,ct,wit,mediumBack,Lw,wi_o,medium,s.x);
}

////////////////////////////////////////////////////////////////////////////////
//                          ThinDielectric Material                               //
////////////////////////////////////////////////////////////////////////////////

 void ThinDielectricMaterial__preprocess(ThinDielectricMaterial* This, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}
 
 Vec3fa ThinDielectricMaterial__eval(ThinDielectricMaterial* This, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) {
  return Vec3fa(0.0f);
}

 Vec3fa ThinDielectricMaterial__sample(ThinDielectricMaterial* This, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  float cosThetaO = clamp(dot(wo,dg.Ns));
  if (cosThetaO <= 0.0f) return Vec3fa(0.0f);
  float R = fresnelDielectric(cosThetaO,rcp(This->eta));
  Sample3f wit = Sample3f(neg(wo),1.0f);
  Sample3f wis = reflect_(wo,dg.Ns);
  Vec3fa ct = exp(Vec3fa(This->transmission)*rcp(cosThetaO))*Vec3fa(1.0f-R);
  Vec3fa cs = Vec3fa(R);
  return sample_component2(cs,wis,medium,ct,wit,medium,Lw,wi_o,medium,s.x);
}

////////////////////////////////////////////////////////////////////////////////
//                     MetallicPaint Material                                 //
////////////////////////////////////////////////////////////////////////////////

 void MetallicPaintMaterial__preprocess(MetallicPaintMaterial* material, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
}

 Vec3fa MetallicPaintMaterial__eval(MetallicPaintMaterial* This, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi) 
{
  DielectricReflection reflection; DielectricReflection__Constructor(&reflection, 1.0f, This->eta);
  DielectricLayerLambertian lambertian; DielectricLayerLambertian__Constructor(&lambertian, Vec3fa((float)1.0f), 1.0f, This->eta, make_Lambertian(Vec3fa((Vec3fa)This->shadeColor)));
  return DielectricReflection__eval(&reflection,wo,dg,wi) + DielectricLayerLambertian__eval(&lambertian,wo,dg,wi);
}

 Vec3fa MetallicPaintMaterial__sample(MetallicPaintMaterial* This, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{
  DielectricReflection reflection; DielectricReflection__Constructor(&reflection, 1.0f, This->eta);
  DielectricLayerLambertian lambertian; DielectricLayerLambertian__Constructor(&lambertian, Vec3fa((float)1.0f), 1.0f, This->eta, make_Lambertian(Vec3fa((Vec3fa)This->shadeColor)));
  Sample3f wi0; Vec3fa c0 = DielectricReflection__sample(&reflection,wo,dg,wi0,s);
  Sample3f wi1; Vec3fa c1 = DielectricLayerLambertian__sample(&lambertian,wo,dg,wi1,s);
  return sample_component2(c0,wi0,medium,c1,wi1,medium,Lw,wi_o,medium,s.x);
}

////////////////////////////////////////////////////////////////////////////////
//                              Material                                      //
////////////////////////////////////////////////////////////////////////////////

inline void Material__preprocess(ISPCMaterial* materials, int materialID, int numMaterials, BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Medium& medium)  
{
  
  {
    
  ISPCMaterial* material = &materials[materialID];
  switch (material->ty) {
  case MATERIAL_OBJ  : OBJMaterial__preprocess  ((OBJMaterial*)  material,brdf,wo,dg,medium); break;
  case MATERIAL_METAL: MetalMaterial__preprocess((MetalMaterial*)material,brdf,wo,dg,medium); break;
  case MATERIAL_REFLECTIVE_METAL: ReflectiveMetalMaterial__preprocess((ReflectiveMetalMaterial*)material,brdf,wo,dg,medium); break;
  case MATERIAL_VELVET: VelvetMaterial__preprocess((VelvetMaterial*)material,brdf,wo,dg,medium); break;
  case MATERIAL_DIELECTRIC: DielectricMaterial__preprocess((DielectricMaterial*)material,brdf,wo,dg,medium); break;
  case MATERIAL_METALLIC_PAINT: MetallicPaintMaterial__preprocess((MetallicPaintMaterial*)material,brdf,wo,dg,medium); break;
  case MATERIAL_MATTE: MatteMaterial__preprocess((MatteMaterial*)material,brdf,wo,dg,medium); break;
  case MATERIAL_MIRROR: MirrorMaterial__preprocess((MirrorMaterial*)material,brdf,wo,dg,medium); break;
  case MATERIAL_THIN_DIELECTRIC: ThinDielectricMaterial__preprocess((ThinDielectricMaterial*)material,brdf,wo,dg,medium); break;
  default: break;
  }
  }
}

inline Vec3fa Material__eval(ISPCMaterial* materials, int materialID, int numMaterials, const BRDF& brdf, const Vec3fa& wo, const DifferentialGeometry& dg, const Vec3fa& wi)
{
  Vec3fa c = Vec3fa(0.0f);
  
  {
    
  ISPCMaterial* material = &materials[materialID];
  switch (material->ty) {
  case MATERIAL_OBJ  : c = OBJMaterial__eval  ((OBJMaterial*)  material, brdf, wo, dg, wi); break;
  case MATERIAL_METAL: c = MetalMaterial__eval((MetalMaterial*)material, brdf, wo, dg, wi); break;
  case MATERIAL_REFLECTIVE_METAL: c = ReflectiveMetalMaterial__eval((ReflectiveMetalMaterial*)material, brdf, wo, dg, wi); break;
  case MATERIAL_VELVET: c = VelvetMaterial__eval((VelvetMaterial*)material, brdf, wo, dg, wi); break;
  case MATERIAL_DIELECTRIC: c = DielectricMaterial__eval((DielectricMaterial*)material, brdf, wo, dg, wi); break;
  case MATERIAL_METALLIC_PAINT: c = MetallicPaintMaterial__eval((MetallicPaintMaterial*)material, brdf, wo, dg, wi); break;
  case MATERIAL_MATTE: c = MatteMaterial__eval((MatteMaterial*)material, brdf, wo, dg, wi); break;
  case MATERIAL_MIRROR: c = MirrorMaterial__eval((MirrorMaterial*)material, brdf, wo, dg, wi); break;
  case MATERIAL_THIN_DIELECTRIC: c = ThinDielectricMaterial__eval((ThinDielectricMaterial*)material, brdf, wo, dg, wi); break;
  default: c = Vec3fa(0.0f); 
  }
  }
  return c;
}

inline Vec3fa Material__sample(ISPCMaterial* materials, int materialID, int numMaterials, const BRDF& brdf, const Vec3fa& Lw, const Vec3fa& wo, const DifferentialGeometry& dg, Sample3f& wi_o, Medium& medium, const Vec2f& s)  
{  
  Vec3fa c = Vec3fa(0.0f);
  
  {
    
  ISPCMaterial* material = &materials[materialID];
  switch (material->ty) {
  case MATERIAL_OBJ  : c = OBJMaterial__sample  ((OBJMaterial*)  material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_METAL: c = MetalMaterial__sample((MetalMaterial*)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_REFLECTIVE_METAL: c = ReflectiveMetalMaterial__sample((ReflectiveMetalMaterial*)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_VELVET: c = VelvetMaterial__sample((VelvetMaterial*)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_DIELECTRIC: c = DielectricMaterial__sample((DielectricMaterial*)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_METALLIC_PAINT: c = MetallicPaintMaterial__sample((MetallicPaintMaterial*)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_MATTE: c = MatteMaterial__sample((MatteMaterial*)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_MIRROR: c = MirrorMaterial__sample((MirrorMaterial*)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  case MATERIAL_THIN_DIELECTRIC: c = ThinDielectricMaterial__sample((ThinDielectricMaterial*)material, brdf, Lw, wo, dg, wi_o, medium, s); break;
  default: c = Vec3fa(0.0f); 
  }
  }
  return c;
}

////////////////////////////////////////////////////////////////////////////////
//                               Scene                                        //
////////////////////////////////////////////////////////////////////////////////

/* scene data */
extern "C" ISPCScene* g_ispc_scene;
RTCScene g_scene = NULL;
void** geomID_to_mesh = NULL;
int* geomID_to_type = NULL;

/* render function to use */
renderPixelFunc renderPixel;

/* occlusion filter function */
void occlusionFilterReject(void* ptr, RTCRay& ray) {
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
}

/* rtcCommitThread called by all ISPC worker threads to enable parallel build */
#if defined(PARALLEL_COMMIT)
task void parallelCommit(RTCScene scene) {
  rtcCommitThread (scene,threadIndex,threadCount); 
}
#endif

/* error reporting function */
void error_handler(const RTCError code, const int8* str)
{
  printf("Embree: ");
  switch (code) {
  case RTC_UNKNOWN_ERROR    : printf("RTC_UNKNOWN_ERROR"); break;
  case RTC_INVALID_ARGUMENT : printf("RTC_INVALID_ARGUMENT"); break;
  case RTC_INVALID_OPERATION: printf("RTC_INVALID_OPERATION"); break;
  case RTC_OUT_OF_MEMORY    : printf("RTC_OUT_OF_MEMORY"); break;
  case RTC_UNSUPPORTED_CPU  : printf("RTC_UNSUPPORTED_CPU"); break;
  default                   : printf("invalid error code"); break;
  }
  if (str) { 
    printf(" ("); 
    while (*str) putchar(*str++); 
    printf(")\n"); 
  }
  abort();
} // error handler

/* accumulation buffer */
Vec3fa* g_accu = NULL;
unsigned int g_accu_width = 0;
unsigned int g_accu_height = 0;
unsigned int g_accu_count = 0;
Vec3fa g_accu_vx;
Vec3fa g_accu_vy;
Vec3fa g_accu_vz;
Vec3fa g_accu_p;
extern "C" bool g_changed;

/* called by the C++ code for initialization */
extern "C" void device_init (int8* cfg)
{
  /* initialize last seen camera */
  g_accu_vx = Vec3fa(0.0f);
  g_accu_vy = Vec3fa(0.0f);
  g_accu_vz = Vec3fa(0.0f);
  g_accu_p  = Vec3fa(0.0f);

  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* set start render mode */
  renderPixel = renderPixelStandard;
  //  renderPixel = renderPixelEyeLight;

} // device_init

void convertTriangleMeshes(ISPCScene* scene_in, RTCScene scene_out, size_t numGeometries)
{
  /* add all meshes to the scene */
  for (int i=0; i<scene_in->numMeshes; i++)
  {
    /* get ith mesh */
    ISPCMesh* mesh = scene_in->meshes[i];

    /* create a triangle mesh */
    unsigned int geomID = rtcNewTriangleMesh (scene_out, RTC_GEOMETRY_STATIC, mesh->numTriangles, mesh->numVertices);
    assert(geomID < numGeometries);
    geomID_to_mesh[geomID] = mesh;
    geomID_to_type[geomID] = 0;
    
    /* set vertices */
    Vertex* vertices = (Vertex*) rtcMapBuffer(scene_out,geomID,RTC_VERTEX_BUFFER); 
    for (int j=0; j<mesh->numVertices; j++) {
      vertices[j].x = mesh->positions[j].x;
      vertices[j].y = mesh->positions[j].y;
      vertices[j].z = mesh->positions[j].z;
    }
    rtcUnmapBuffer(scene_out,geomID,RTC_VERTEX_BUFFER); 

    /* set triangles */
    Triangle* triangles = (Triangle*) rtcMapBuffer(scene_out,geomID,RTC_INDEX_BUFFER);
    for (int j=0; j<mesh->numTriangles; j++) {
      triangles[j].v0 = mesh->triangles[j].v0;
      triangles[j].v1 = mesh->triangles[j].v1;
      triangles[j].v2 = mesh->triangles[j].v2;
    }
    rtcUnmapBuffer(scene_out,geomID,RTC_INDEX_BUFFER);

    bool allOpaque = true;
    bool allTransparent = true;
    for (size_t j=0; j<mesh->numTriangles; j++) {
      ISPCTriangle triangle = mesh->triangles[j];
      if (g_ispc_scene->materials[triangle.materialID].ty == MATERIAL_DIELECTRIC ||
	  g_ispc_scene->materials[triangle.materialID].ty == MATERIAL_THIN_DIELECTRIC)
	allOpaque = false;
      else 
	allTransparent = false;
    }
    if (allTransparent)
      rtcSetOcclusionFilterFunction(scene_out,geomID,(RTCFilterFunc)&occlusionFilterReject);
  }
}

void convertSubdivMeshes(ISPCScene* scene_in, RTCScene scene_out, size_t numGeometries)
{
  for (size_t i=0; i<g_ispc_scene->numSubdivMeshes; i++)
  {
    ISPCSubdivMesh* mesh = g_ispc_scene->subdiv[i];
    unsigned int geomID = rtcNewSubdivisionMesh(scene_out, RTC_GEOMETRY_STATIC, mesh->numFaces, mesh->numEdges, mesh->numVertices, 
						mesh->numEdgeCreases, mesh->numVertexCreases, mesh->numHoles);
    assert(geomID < numGeometries);
    geomID_to_mesh[geomID] = mesh;
    geomID_to_type[geomID] = 1;

    for (size_t i=0; i<mesh->numEdges; i++) mesh->subdivlevel[i] = 16;
    rtcSetBuffer(scene_out, geomID, RTC_VERTEX_BUFFER, mesh->positions, 0, sizeof(Vec3fa  ));
    rtcSetBuffer(scene_out, geomID, RTC_LEVEL_BUFFER,  mesh->subdivlevel, 0, sizeof(float));
    rtcSetBuffer(scene_out, geomID, RTC_INDEX_BUFFER,  mesh->position_indices  , 0, sizeof(unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_FACE_BUFFER,   mesh->verticesPerFace, 0, sizeof(unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_HOLE_BUFFER,   mesh->holes, 0, sizeof(unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_EDGE_CREASE_INDEX_BUFFER,    mesh->edge_creases,          0, 2*sizeof(unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_EDGE_CREASE_WEIGHT_BUFFER,   mesh->edge_crease_weights,   0, sizeof(float));
    rtcSetBuffer(scene_out, geomID, RTC_VERTEX_CREASE_INDEX_BUFFER,  mesh->vertex_creases,        0, sizeof(unsigned int));
    rtcSetBuffer(scene_out, geomID, RTC_VERTEX_CREASE_WEIGHT_BUFFER, mesh->vertex_crease_weights, 0, sizeof(float));
  }
}      

typedef void* void_ptr;

RTCScene convertScene(ISPCScene* scene_in)
{
  size_t numGeometries = scene_in->numMeshes + scene_in->numSubdivMeshes;
  geomID_to_mesh = new void_ptr[numGeometries];
  geomID_to_type = new int[numGeometries];

  /* create scene */
  RTCScene scene_out = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT, RTC_INTERSECT1);
  convertTriangleMeshes(scene_in,scene_out,numGeometries);
  convertSubdivMeshes(scene_in,scene_out,numGeometries);

  /* commit changes to scene */
#if !defined(PARALLEL_COMMIT)
  rtcCommit (scene_out);
#else
  launch[ getNumHWThreads() ] parallelCommit(scene_out); 
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
  for (int i=0; i<10; i++) irand(state);
}

inline float frand(rand_state& state) {
  return irand(state)*2.3283064365386963e-10f;
}

inline Vec3fa face_forward(const Vec3fa& dir, const Vec3fa& _Ng) {
  const Vec3fa Ng = _Ng;
  return dot(dir,Ng) < 0.0f ? Ng : neg(Ng);
}

#if 0
inline Vec3fa interpolate_normal(RTCRay& ray)
{
#if 1 // FIXME: pointer gather not implemented on ISPC for Xeon Phi
  ISPCMesh* mesh = g_ispc_scene->meshes[ray.geomID];
  ISPCTriangle* tri = &mesh->triangles[ray.primID];

  /* load material ID */
  int materialID = tri->materialID;

  /* interpolate shading normal */
  if (mesh->normals) {
    Vec3fa n0 = Vec3fa(mesh->normals[tri->v0]);
    Vec3fa n1 = Vec3fa(mesh->normals[tri->v1]);
    Vec3fa n2 = Vec3fa(mesh->normals[tri->v2]);
    float u = ray.u, v = ray.v, w = 1.0f-ray.u-ray.v;
    return normalize(w*n0 + u*n1 + v*n2);
  } else {
    return normalize(ray.Ng);
  }

#else

  Vec3fa Ns = Vec3fa(0.0f);
  int materialID = 0;
  foreach_unique (geomID in ray.geomID) 
  {
    if (geomID >= 0 && geomID < g_ispc_scene->numMeshes)  { // FIXME: workaround for ISPC bug

    ISPCMesh* mesh = g_ispc_scene->meshes[geomID];
    
    foreach_unique (primID in ray.primID) 
    {
      ISPCTriangle* tri = &mesh->triangles[primID];
      
      /* load material ID */
      materialID = tri->materialID;

      /* interpolate shading normal */
      if (mesh->normals) {
        Vec3fa n0 = Vec3fa(mesh->normals[tri->v0]);
        Vec3fa n1 = Vec3fa(mesh->normals[tri->v1]);
        Vec3fa n2 = Vec3fa(mesh->normals[tri->v2]);
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

Vec3fa renderPixelFunction(float x, float y, rand_state& state, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* radiance accumulator and weight */
  Vec3fa L = Vec3fa(0.0f);
  Vec3fa Lw = Vec3fa(1.0f);
  Medium medium = make_Medium_Vacuum();

  /* initialize ray */
  RTCRay ray = RTCRay(p,normalize(x*vx + y*vy + vz),0.0f,inf);

  /* iterative path tracer loop */
  for (int i=0; i<8; i++)
  {
    /* terminate if contribution too low */
    if (max(Lw.x,max(Lw.y,Lw.z)) < 0.01f)
      break;

    /* intersect ray with scene */ 
    rtcIntersect(g_scene,ray);
    const Vec3fa wo = neg(ray.dir);
    
    /* invoke environment lights if nothing hit */
    if (ray.geomID == RTC_INVALID_GEOMETRY_ID) 
    {
#if 0
      /* iterate over all ambient lights */
      for (size_t i=0; i<g_ispc_scene->numAmbientLights; i++)
        L = L + Lw*AmbientLight__eval(g_ispc_scene->ambientLights[i],ray.dir); // FIXME: +=
#endif

#if 0
      /* iterate over all distant lights */
      for (size_t i=0; i<g_ispc_scene->numDistantLights; i++)
        L = L + Lw*DistantLight__eval(g_ispc_scene->distantLights[i],ray.dir); // FIXME: +=
#endif

      break;
    }

    /* compute differential geometry */
    DifferentialGeometry dg;
    dg.P  = ray.org+ray.tfar*ray.dir;
    dg.Ng = face_forward(ray.dir,normalize(ray.Ng));
    //Vec3fa _Ns = interpolate_normal(ray);
    Vec3fa _Ns = normalize(ray.Ng);
    dg.Ns = face_forward(ray.dir,_Ns);

    /* shade all rays that hit something */
#if 1 // FIXME: pointer gather not implemented in ISPC for Xeon Phi
    int materialID = 0;
    if (geomID_to_type[ray.geomID] == 0)
      materialID = ((ISPCMesh*) geomID_to_mesh[ray.geomID])->triangles[ray.primID].materialID; 
    else 
      materialID = ((ISPCSubdivMesh*) geomID_to_mesh[ray.geomID])->materialID; 
#else
    int materialID = 0;
    foreach_unique (geomID in ray.geomID) {
      if (geomID >= 0 && geomID < g_ispc_scene->numMeshes) { // FIXME: workaround for ISPC bug
	if (geomID_to_type[geomID] == 0) 
	  materialID = ((ISPCMesh*) geomID_to_mesh[geomID])->triangles[ray.primID].materialID; 
	else                             
	  materialID = ((ISPCSubdivMesh*) geomID_to_mesh[geomID])->materialID; 
      }
    }
#endif

    /*! Compute  simple volumetric effect. */
    Vec3fa c = Vec3fa(1.0f);
    const Vec3fa transmission = medium.transmission;
    if (ne(transmission,Vec3fa(1.0f)))
      c = c * pow(transmission,ray.tfar);
    
    /* calculate BRDF */ // FIXME: avoid gathers
    BRDF brdf;
    int numMaterials = g_ispc_scene->numMaterials;
    //ISPCMaterial* material = &g_ispc_scene->materials[materialID];
    ISPCMaterial* material_array = &g_ispc_scene->materials[0];
    Material__preprocess(material_array,materialID,numMaterials,brdf,wo,dg,medium);

    /* sample BRDF at hit point */
    Sample3f wi1;
    c = c * Material__sample(material_array,materialID,numMaterials,brdf,Lw, wo, dg, wi1, medium, Vec2f(frand(state),frand(state)));

    /* iterate over ambient lights */
    for (size_t i=0; i<g_ispc_scene->numAmbientLights; i++)
    {
#if 1
      Vec3fa L0 = Vec3fa(0.0f);
      Sample3f wi0; float tMax0;
      Vec3fa Ll0 = AmbientLight__sample(g_ispc_scene->ambientLights[i],dg,wi0,tMax0,Vec2f(frand(state),frand(state)));
      if (wi0.pdf > 0.0f) {
        RTCRay shadow = RTCRay(dg.P,wi0.v,0.001f,tMax0);
        rtcOccluded(g_scene,shadow);
        if (shadow.geomID == RTC_INVALID_GEOMETRY_ID) {
          L0 = Ll0/wi0.pdf*Material__eval(material_array,materialID,numMaterials,brdf,wo,dg,wi0.v);
        }
        L = L + Lw*L0;
      }
#endif

#if 0
      Vec3fa L1 = Vec3fa(0.0f);
      Vec3fa Ll1 = AmbientLight__eval(g_ispc_scene->ambientLights[i],wi1.v);
      if (wi1.pdf > 0.0f) {
        RTCRay shadow = RTCRay(dg.P,wi1.v,0.001f,inf);
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
    for (size_t i=0; i<g_ispc_scene->numPointLights; i++)
    {
      Vec3fa Ll = PointLight__sample(g_ispc_scene->pointLights[i],dg,wi,tMax,Vec2f(frand(state),frand(state)));
      if (wi.pdf <= 0.0f) continue;
      RTCRay shadow = RTCRay(dg.P,wi.v,0.001f,tMax);
      rtcOccluded(g_scene,shadow);
      if (shadow.geomID != RTC_INVALID_GEOMETRY_ID) continue;
      L = L + Lw*Ll/wi.pdf*Material__eval(material_array,materialID,numMaterials,brdf,wo,dg,wi.v); // FIXME: +=
    }

    /* iterate over directional lights */
    for (size_t i=0; i<g_ispc_scene->numDirectionalLights; i++)
    {
      Vec3fa Ll = DirectionalLight__sample(g_ispc_scene->dirLights[i],dg,wi,tMax,Vec2f(frand(state),frand(state)));
      if (wi.pdf <= 0.0f) continue;
      RTCRay shadow = RTCRay(dg.P,wi.v,0.001f,tMax);
      rtcOccluded(g_scene,shadow);
      if (shadow.geomID != RTC_INVALID_GEOMETRY_ID) continue;
      L = L + Lw*Ll/wi.pdf*Material__eval(material_array,materialID,numMaterials,brdf,wo,dg,wi.v); // FIXME: +=
    }

    /* iterate over distant lights */
    for (size_t i=0; i<g_ispc_scene->numDistantLights; i++)
    {
      Vec3fa Ll = DistantLight__sample(g_ispc_scene->distantLights[i],dg,wi,tMax,Vec2f(frand(state),frand(state)));
      if (wi.pdf <= 0.0f) continue;
      RTCRay shadow = RTCRay(dg.P,wi.v,0.001f,tMax);
      rtcOccluded(g_scene,shadow);
      if (shadow.geomID != RTC_INVALID_GEOMETRY_ID) continue;
      L = L + Lw*Ll/wi.pdf*Material__eval(material_array,materialID,numMaterials,brdf,wo,dg,wi.v); // FIXME: +=
    }

    if (wi1.pdf <= 0.0f) break;
    Lw = Lw*c/wi1.pdf; // FIXME: *=

    /* setup secondary ray */
    ray = RTCRay(dg.P,normalize(wi1.v),0.001f,inf);
  }
  return L;
}

/* task that renders a single screen tile */
Vec3fa renderPixelStandard(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  rand_state state;
  init_rand(state,
            253*x+35*y+152*g_accu_count+54,
            1253*x+345*y+1452*g_accu_count+564,
            10253*x+3435*y+52*g_accu_count+13);

  Vec3fa L = Vec3fa(0.0f,0.0f,0.0f);
  //for (int i=0; i<16; i++) {
  L = L + renderPixelFunction(x,y,state,vx,vy,vz,p); // FIXME: +=
  //}
  //L = L*(1.0f/16.0f);
  return L;
}
  
/* task that renders a single screen tile */
void renderTile(int taskIndex, int* pixels,
                     const int width,
                     const int height, 
                     const float time,
                     const Vec3fa& vx, 
                     const Vec3fa& vy, 
                     const Vec3fa& vz, 
                     const Vec3fa& p,
                     const int numTilesX, 
                     const int numTilesY)
{
  const int tileY = taskIndex / numTilesX;
  const int tileX = taskIndex - tileY * numTilesX;
  const int x0 = tileX * TILE_SIZE_X;
  const int x1 = min(x0+TILE_SIZE_X,width);
  const int y0 = tileY * TILE_SIZE_Y;
  const int y1 = min(y0+TILE_SIZE_Y,height);

  for (int y = y0; y<y1; y++) for (int x = x0; x<x1; x++)
  {
    //if (x != 200 || y != 450) continue;

    /* calculate pixel color */
    Vec3fa color = renderPixel(x,y,vx,vy,vz,p);

    /* write color to framebuffer */
    Vec3fa* dst = &g_accu[y*width+x];
    *dst = *dst + Vec3fa(color.x,color.y,color.z,1.0f); // FIXME: use += operator
    float f = rcp(max(0.001f,dst->w));
    unsigned int r = (unsigned int) (255.0f * clamp(dst->x*f,0.0f,1.0f));
    unsigned int g = (unsigned int) (255.0f * clamp(dst->y*f,0.0f,1.0f));
    unsigned int b = (unsigned int) (255.0f * clamp(dst->z*f,0.0f,1.0f));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
} // renderTile

/* called by the C++ code to render */
extern "C" void device_render (int* pixels,
                           const int width,
                           const int height, 
                           const float time,
                           const Vec3fa& vx, 
                           const Vec3fa& vy, 
                           const Vec3fa& vz, 
                           const Vec3fa& p)
{
  /* create scene */
  if (g_scene == NULL)
    g_scene = convertScene(g_ispc_scene);

  /* create accumulator */
  if (g_accu_width != width || g_accu_height != height) {
    alignedFree(g_accu);
    g_accu = (Vec3fa*) alignedMalloc(width*height*sizeof(Vec3fa));
    g_accu_width = width;
    g_accu_height = height;
    memset(g_accu,0,width*height*sizeof(Vec3fa));
  }

  /* reset accumulator */
  bool camera_changed = g_changed; g_changed = false;
  camera_changed |= ne(g_accu_vx,vx); g_accu_vx = vx; // FIXME: use != operator
  camera_changed |= ne(g_accu_vy,vy); g_accu_vy = vy; // FIXME: use != operator
  camera_changed |= ne(g_accu_vz,vz); g_accu_vz = vz; // FIXME: use != operator
  camera_changed |= ne(g_accu_p,  p); g_accu_p  = p;  // FIXME: use != operator
  g_accu_count++;
  if (camera_changed) {
    g_accu_count=0;
    memset(g_accu,0,width*height*sizeof(Vec3fa));
  }

  /* render image */
  const int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch_renderTile(numTilesX*numTilesY,pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); 
  rtcDebug();
} // device_render

/* called by the C++ code for cleanup */
extern "C" void device_cleanup ()
{
  alignedFree(g_accu);
  rtcDeleteScene (g_scene);
  rtcExit();
} // device_cleanup
