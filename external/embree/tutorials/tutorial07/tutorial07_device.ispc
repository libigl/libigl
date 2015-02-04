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

#define PARALLEL_COMMIT

#if defined(__XEON_PHI__) // FIXME: gather of pointers not working in ISPC for Xeon Phi
#define renderPixelTestEyeLight renderPixelStandard
#else
#define renderPixelPathTrace renderPixelStandard
#endif

/* accumulation buffer */
uniform Vec3fa* uniform g_accu = NULL;
uniform size_t g_accu_width = 0;
uniform size_t g_accu_height = 0;
uniform size_t g_accu_count = 0;
uniform Vec3f g_accu_vx;
uniform Vec3f g_accu_vy;
uniform Vec3f g_accu_vz;
uniform Vec3f g_accu_p;
extern uniform bool g_changed;

/* light settings */
extern uniform Vec3fa g_dirlight_direction;
extern uniform Vec3fa g_dirlight_intensity;
extern uniform Vec3fa g_ambient_intensity;

/* hair material */
uniform Vec3f hair_K;
uniform Vec3f hair_dK;
uniform Vec3f hair_Kr;    //!< reflectivity of hair
uniform Vec3f hair_Kt;    //!< transparency of hair

void filterDispatch(void* uniform ptr, struct RTCRay2& ray);

/* scene data */
extern uniform ISPCScene* uniform g_ispc_scene;
RTCScene g_scene = NULL;

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
}

/* rtcCommitThread called by all ISPC worker threads to enable parallel build */
#if defined(PARALLEL_COMMIT)
task void parallelCommit(RTCScene scene) {
  rtcCommitThread (scene,threadIndex,threadCount); 
}
#endif

/* render function to use */
renderPixelFunc renderPixel;

Vec3f renderPixelTestEyeLight(float x, float y, const uniform Vec3f& vx, const uniform Vec3f& vy, const uniform Vec3f& vz, const uniform Vec3f& p);

/*! random number generator for floating point numbers in range [0,1] */
inline float frand(int& seed) {
  seed = 7 * seed + 5;
  return (seed & 0xFFFF)/(float)0xFFFF;
}

/*! Uniform hemisphere sampling. Up direction is the z direction. */
Vec3fa sampleSphere(const float u, const float v) 
{
  const float phi = 2.0f*(float)pi * u;
  const float cosTheta = 1.0f - 2.0f * v, sinTheta = 2.0f * sqrt(v * (1.0f - v));
  return make_Vec3fa(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta, one_over_four_pi);
}

RTCScene convertScene(uniform ISPCScene* uniform scene_in)
{
  //scene_in->numHairSets = 0;
  //scene_in->numMeshes = 0;

  /* create scene */
  RTCScene scene_out = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT, RTC_INTERSECT_VARYING);

  /* add all hair sets to the scene */
  for (uniform int i=0; i<scene_in->numHairSets; i++)
  {
    uniform ISPCHairSet* uniform hair = scene_in->hairs[i];
    uniform unsigned int geomID = rtcNewHairGeometry (scene_out, RTC_GEOMETRY_STATIC, hair->numHairs, hair->numVertices, hair->v2 ? 2 : 1);
    rtcSetBuffer(scene_out,geomID,RTC_VERTEX_BUFFER,hair->v,0,sizeof(uniform Vertex));
    if (hair->v2) rtcSetBuffer(scene_out,geomID,RTC_VERTEX_BUFFER1,hair->v2,0,sizeof(uniform Vertex));
    rtcSetBuffer(scene_out,geomID,RTC_INDEX_BUFFER,hair->hairs,0,sizeof(uniform ISPCHair));
    rtcSetOcclusionFilterFunction(scene_out,geomID,(RTCFilterFuncVarying)&filterDispatch);
  }

  /* add all triangle meshes to the scene */
  for (uniform int i=0; i<scene_in->numMeshes; i++)
  {
    uniform ISPCMesh* uniform mesh = scene_in->meshes[i];
    uniform unsigned int geomID = rtcNewTriangleMesh (scene_out, RTC_GEOMETRY_STATIC, mesh->numTriangles, mesh->numVertices, mesh->positions2 ? 2 : 1);
    rtcSetBuffer(scene_out,geomID,RTC_VERTEX_BUFFER,mesh->positions,0,sizeof(uniform Vertex));
    if (mesh->positions2) rtcSetBuffer(scene_out,geomID,RTC_VERTEX_BUFFER1,mesh->positions2,0,sizeof(uniform Vertex));
    rtcSetBuffer(scene_out,geomID,RTC_INDEX_BUFFER,mesh->triangles,0,sizeof(uniform ISPCTriangle));
    rtcSetOcclusionFilterFunction(scene_out,geomID,(RTCFilterFuncVarying)&filterDispatch);
  }

  /* commit changes to scene */
#if !defined(PARALLEL_COMMIT)
  rtcCommit (scene_out);
#else
  launch[ getNumHWThreads() ] parallelCommit(scene_out); sync;
#endif

  return scene_out;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


/* called by the C++ code for initialization */
export void device_init (uniform int8* uniform cfg)
{
  /* initialize last seen camera */
  g_accu_vx = make_Vec3f(0.0f);
  g_accu_vy = make_Vec3f(0.0f);
  g_accu_vz = make_Vec3f(0.0f);
  g_accu_p  = make_Vec3f(0.0f);

  /* initialize hair colors */
  hair_K  = make_Vec3f(0.8f,0.57f,0.32f);
  hair_dK = make_Vec3f(0.1f,0.12f,0.08f);
  hair_Kr = 0.2f*hair_K;    //!< reflectivity of hair
  hair_Kt = 0.8f*hair_K;    //!< transparency of hair

  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* set start render mode */
  renderPixel = renderPixelStandard;

}

#if !defined(__XEON_PHI__)

/*! Anisotropic power cosine microfacet distribution. */
struct AnisotropicBlinn {
  Vec3f dx;       //!< x-direction of the distribution.
  float nx;        //!< Glossiness in x direction with range [0,infinity[ where 0 is a diffuse surface.
  Vec3f dy;       //!< y-direction of the distribution.
  float ny;        //!< Exponent that determines the glossiness in y direction.
  Vec3f dz;       //!< z-direction of the distribution.
  float norm1;     //!< Normalization constant for calculating the pdf for sampling.
  float norm2;     //!< Normalization constant for calculating the distribution.
  Vec3f Kr,Kt; 
  float side;
};

  /*! Anisotropic power cosine distribution constructor. */
inline void AnisotropicBlinn__Constructor(varying AnisotropicBlinn* uniform This, const uniform Vec3f& Kr, const uniform Vec3f& Kt, 
                                          const Vec3f& dx, float nx, const Vec3f& dy, float ny, const Vec3f& dz) 
{
  This->Kr = Kr;
  This->Kt = Kt;
  This->dx = dx;
  This->nx = nx;
  This->dy = dy;
  This->ny = ny;
  This->dz = dz;
  This->norm1 = sqrtf((nx+1)*(ny+1)) * one_over_two_pi;
  This->norm2 = sqrtf((nx+2)*(ny+2)) * one_over_two_pi;
  This->side = reduce_max(Kr)/(reduce_max(Kr)+reduce_max(Kt));
}

/*! Evaluates the power cosine distribution. \param wh is the half
 *  vector */
inline float AnisotropicBlinn__eval(const varying AnisotropicBlinn* uniform This, const Vec3f& wh)  
{
  const float cosPhiH   = dot(wh, This->dx);
  const float sinPhiH   = dot(wh, This->dy);
  const float cosThetaH = dot(wh, This->dz);
  const float R = sqr(cosPhiH)+sqr(sinPhiH);
  if (R == 0.0f) return This->norm2;
  const float n = (This->nx*sqr(cosPhiH)+This->ny*sqr(sinPhiH))*rcp(R);
  return This->norm2 * pow(abs(cosThetaH), n);
}

/*! Samples the distribution. \param s is the sample location
 *  provided by the caller. */
inline Vec3fa AnisotropicBlinn__sample(const varying AnisotropicBlinn* uniform This, const float sx, const float sy)
{
  const float phi = two_pi*sx;
  const float sinPhi0 = sqrtf(This->nx+1)*sinf(phi);
  const float cosPhi0 = sqrtf(This->ny+1)*cosf(phi);
  const float norm = rsqrt(sqr(sinPhi0)+sqr(cosPhi0));
  const float sinPhi = sinPhi0*norm;
  const float cosPhi = cosPhi0*norm;
  const float n = This->nx*sqr(cosPhi)+This->ny*sqr(sinPhi);
  const float cosTheta = powf(sy,rcp(n+1));
  const float sinTheta = cos2sin(cosTheta);
  const float pdf = This->norm1*powf(cosTheta,n);
  const Vec3f h = make_Vec3f(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
  const Vec3f wh = h.x*This->dx + h.y*This->dy + h.z*This->dz;
  return make_Vec3fa(wh,pdf);
}

inline Vec3f AnisotropicBlinn__eval(const varying AnisotropicBlinn* uniform This, const Vec3f& wo, const uniform Vec3f& wi) 
{
  const float cosThetaI = dot(wi,This->dz);
  
  /* reflection */
  if (cosThetaI > 0.0f) {
    const Vec3f wh = normalize(wi + wo);
    return This->Kr * AnisotropicBlinn__eval(This,wh) * abs(cosThetaI);
  } 
  
  /* transmission */
  else {
    const Vec3f wh = normalize(reflect(wi,This->dz) + wo);
    return This->Kt * AnisotropicBlinn__eval(This,wh) * abs(cosThetaI);
  }
}

inline Vec3f AnisotropicBlinn__sample(const varying AnisotropicBlinn* uniform This, const Vec3f& wo, Vec3fa& wi, const float sx, const float sy, const float sz) 
{
  //wi = Vec3f(reflect(normalize(wo),normalize(dz)),1.0f); return Kr;
  //wi = Vec3f(neg(wo),1.0f); return Kt;
  const Vec3fa wh = AnisotropicBlinn__sample(This,sx,sy);
  //if (dot(wo,wh) < 0.0f) return Vec3f(0.0f,0.0f);
  
  /* reflection */
  if (sz < This->side) {
    wi = make_Vec3fa(reflect(wo,make_Vec3f(wh)),wh.w*This->side);
    const float cosThetaI = dot(make_Vec3f(wi),This->dz);
    return This->Kr * AnisotropicBlinn__eval(This,make_Vec3f(wh)) * abs(cosThetaI);
  }
  
  /* transmission */
  else {
    wi = make_Vec3fa(reflect(reflect(wo,make_Vec3f(wh)),This->dz),wh.w*(1-This->side));
    const float cosThetaI = dot(make_Vec3f(wi),This->dz);
    return This->Kt * AnisotropicBlinn__eval(This,make_Vec3f(wh)) * abs(cosThetaI);
  }
}

typedef uniform Vec3fa* uniform uniform_Vec3fa_ptr;

inline Vec3fa evalBezier(const int geomID, const int primID, const float t)
{
  const float t0 = 1.0f - t, t1 = t;
  const ISPCHairSet* hair = g_ispc_scene->hairs[geomID];
  const Vec3fa* vertices = hair->v;
  const ISPCHair* hairs = hair->hairs;
  
  const int i = hairs[primID].vertex;
  const Vec3fa p00 = vertices[i+0];
  const Vec3fa p01 = vertices[i+1];
  const Vec3fa p02 = vertices[i+2];
  const Vec3fa p03 = vertices[i+3];
  
  const Vec3fa p10 = p00 * t0 + p01 * t1;
  const Vec3fa p11 = p01 * t0 + p02 * t1;
  const Vec3fa p12 = p02 * t0 + p03 * t1;
  const Vec3fa p20 = p10 * t0 + p11 * t1;
  const Vec3fa p21 = p11 * t0 + p12 * t1;
  const Vec3fa p30 = p20 * t0 + p21 * t1;
  
  return p30;
  //tangent = p21-p20;
}

#endif 

/* extended ray structure that includes total transparency along the ray */
struct RTCRay2
{
  Vec3f org;     //!< Ray origin
  Vec3f dir;     //!< Ray direction
  float tnear;   //!< Start of ray segment
  float tfar;    //!< End of ray segment
  float time;    //!< Time of this ray for motion blur.
  int mask;      //!< used to mask out objects during traversal
  Vec3f Ng;      //!< Geometric normal.
  float u;       //!< Barycentric u coordinate of hit
  float v;       //!< Barycentric v coordinate of hit
  int geomID;    //!< geometry ID
  int primID;    //!< primitive ID
  int instID;    //!< instance ID

  // ray extensions
  RTCFilterFuncVarying filter;
  Vec3f transparency; //!< accumulated transparency value
};

uniform bool enableFilterDispatch = false;

/* filter dispatch function */
void filterDispatch(void* uniform ptr, RTCRay2& ray) {
  if (!enableFilterDispatch) return;
  if (ray.filter) ray.filter(ptr,*((varying RTCRay*)&ray)); // FIXME: use RTCRay& cast
}

#if !defined(__XEON_PHI__)

/* occlusion filter function */
void occlusionFilter(void* uniform ptr, RTCRay2& ray)
{
  /* make all surfaces opaque */
  if (ray.geomID >= g_ispc_scene->numHairSets) {
    ray.transparency = make_Vec3f(0.0f);
    return;
  }
  Vec3f T = hair_Kt;
  T = T * ray.transparency; // FIXME: use *= operator
  ray.transparency = T;
  if (ne(T,make_Vec3f(0.0f))) ray.geomID = RTC_INVALID_GEOMETRY_ID; // FIXME: use != operator
}

Vec3f occluded(RTCScene scene, RTCRay2& ray)
{
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.filter = (RTCFilterFuncVarying) &occlusionFilter;
  ray.transparency = make_Vec3f(1.0f);
  rtcOccluded(scene,*((varying RTCRay*)&ray)); // FIXME: use (RTCRay&) cast
  return ray.transparency;
}

/* task that renders a single screen tile */
Vec3f renderPixelPathTrace(float x, float y, const uniform Vec3f& vx, const uniform Vec3f& vy, const uniform Vec3f& vz, const uniform Vec3f& p)
{
  int seed = 21344*x+121233*y+234532*g_accu_count;
  float time = frand(seed);

  /* initialize ray */
  RTCRay2 ray;
  ray.org = p;
  ray.dir = normalize(x*vx + y*vy + vz);
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = time;
  ray.filter = NULL; 
  
  Vec3f color = make_Vec3f(0.0f);
  Vec3f weight = make_Vec3f(1.0f);
  size_t depth = 0;

  while (true)
  {
    /* terminate ray path */
    if (reduce_max(weight) < 0.01 || depth > 20) 
      return color;

    /* intersect ray with scene and gather all hits */
    rtcIntersect(g_scene,*((varying RTCRay*)&ray)); // FIXME: use (RTCRay&) cast
    
    /* exit if we hit environment */
    if (ray.geomID == RTC_INVALID_GEOMETRY_ID) 
      return color + weight*make_Vec3f(g_ambient_intensity);

    /* calculate transmissivity of hair */
    AnisotropicBlinn brdf;
    float tnear_eps = 0.0001f;

    if (ray.geomID < g_ispc_scene->numHairSets) 
    {
      /* calculate tangent space */
      const Vec3f dx = normalize(ray.Ng);
      const Vec3f dy = normalize(cross(ray.dir,dx));
      const Vec3f dz = normalize(cross(dy,dx));

      /* generate anisotropic BRDF */
      AnisotropicBlinn__Constructor(&brdf,hair_Kr,hair_Kt,dx,20.0f,dy,2.0f,dz);
      brdf.Kr = hair_Kr;
      Vec3fa p = evalBezier(ray.geomID,ray.primID,ray.u);
      tnear_eps = 1.1f*p.w;
    }
    else 
    {
      int meshID = ray.geomID-g_ispc_scene->numHairSets;
      ISPCMesh* mesh = g_ispc_scene->meshes[meshID];
      ISPCTriangle* triangle = &mesh->triangles[ray.primID];
      OBJMaterial* material = (OBJMaterial*) &g_ispc_scene->materials[triangle->materialID];
      if (material->illum == 1)
      {
        /* calculate tangent space */
        const Vec3f dx = normalize(make_Vec3f(mesh->normals[triangle->v0]));
        const Vec3f dy = normalize(cross(ray.dir,dx));
        const Vec3f dz = normalize(cross(dy,dx));
        
        /* generate anisotropic BRDF */
        AnisotropicBlinn__Constructor(&brdf,hair_Kr,hair_Kt,dx,20.0f,dy,2.0f,dz);
        brdf.Kr = hair_Kr;
        tnear_eps = 1.1f*mesh->texcoords[triangle->v0].x;
      }
      else
      {
        if (dot(ray.dir,ray.Ng) > 0) ray.Ng = neg(ray.Ng);
        
        /* calculate tangent space */
        const Vec3f dz = normalize(ray.Ng);
        const Vec3f dx = normalize(cross(dz,ray.dir));
        const Vec3f dy = normalize(cross(dz,dx));
        
        /* generate isotropic BRDF */
        AnisotropicBlinn__Constructor(&brdf,make_Vec3f(1.0f),make_Vec3f(0.0f),dx,1.0f,dy,1.0f,dz);
      }
    }
    
    /* sample directional light */
    RTCRay2 shadow;
    shadow.org = ray.org + ray.tfar*ray.dir;
    shadow.dir = neg(make_Vec3f(g_dirlight_direction));
    shadow.tnear = tnear_eps;
    shadow.tfar = inf;
    shadow.time = time;
    Vec3f T = occluded(g_scene,shadow);
    Vec3f c = AnisotropicBlinn__eval(&brdf,neg(ray.dir),neg(make_Vec3f(g_dirlight_direction)));
    color = color + weight*c*T*make_Vec3f(g_dirlight_intensity); // FIXME: use += operator

#if 1
    /* sample BRDF */
    Vec3fa wi;
    c = AnisotropicBlinn__sample(&brdf,neg(ray.dir),wi,frand(seed),frand(seed),frand(seed));
    if (wi.w <= 0.0f) return color;

    /* calculate secondary ray and offset it out of the hair */
    float sign = dot(make_Vec3f(wi),brdf.dz) < 0.0f ? -1.0f : 1.0f;
    ray.org = ray.org + ray.tfar*ray.dir + sign*tnear_eps*brdf.dz;
    ray.dir = make_Vec3f(wi);
    ray.tnear = 0.001f;
    ray.tfar = inf;
    ray.geomID = RTC_INVALID_GEOMETRY_ID;
    ray.primID = RTC_INVALID_GEOMETRY_ID;
    ray.mask = -1;
    ray.time = time;
    ray.filter = NULL;
    weight = weight * c/wi.w; // FIXME: use *= operator

#else    

    /* continue with transparency ray */
    ray.geomID = RTC_INVALID_GEOMETRY_ID;
    ray.tnear = 1.001f*ray.tfar;
    ray.tfar = inf;
    weight *= brdf.Kt;
    
#endif

    depth++;
  }
  return color;
}

#endif

Vec3f renderPixelTestEyeLight(float x, float y, const uniform Vec3f& vx, const uniform Vec3f& vy, const uniform Vec3f& vz, const uniform Vec3f& p)
{
  /* initialize ray */
  RTCRay2 ray;
  ray.org = p;
  ray.dir = normalize(x*vx + y*vy + vz);
  Vec3f dir1 = normalize((x+1)*vx + (y+1)*vy + vz);
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  Vec3f color = make_Vec3f(0.0f);
  float weight = 1.0f;

  rtcIntersect(g_scene,*((varying RTCRay*)&ray)); // FIXME: use (RTCRay&) cast
  ray.filter = NULL; 

  if (ray.primID == -1)
    return make_Vec3f(0.0f);
  
  Vec3f Ng;
  if (ray.geomID < g_ispc_scene->numHairSets) 
  {
    const Vec3f dx = normalize(ray.Ng);
    const Vec3f dy = normalize(cross(ray.dir,dx));
    const Vec3f dz = normalize(cross(dy,dx));
    Ng = dz;
  }
  else 
  {
    if (dot(ray.dir,ray.Ng) > 0) ray.Ng = neg(ray.Ng);
    const Vec3f dz = normalize(ray.Ng);
    const Vec3f dx = normalize(cross(dz,ray.dir));
    const Vec3f dy = normalize(cross(dz,dx));
    Ng = dz;
  }

  color = color + make_Vec3f(0.2f + 0.5f * abs(dot(ray.dir,Ng))); // FIXME: use += operator
  return color;
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
  //int seed = tileY*numTilesX+tileX+programIndex + g_accu_count;
  int seed = (tileY*numTilesX+tileX+programIndex) * g_accu_count;

  foreach (y = y0 ... y1, x = x0 ... x1)
  {
   int seed = (y*width+x+programIndex) * g_accu_count;

    /* calculate pixel color */
     float subpixel_x = frand(seed);
     float subpixel_y = frand(seed);
     float fx = x + subpixel_x;
     float fy = y + subpixel_y;
     Vec3f color = renderPixel(fx,fy,vx,vy,vz,p);
    
    //Vec3f color = renderPixelTestEyeLight(fx,fy,vx,vy,vz,p);

    /* write color to framebuffer */
    Vec3fa* dst = &g_accu[y*width+x];
    *dst = *dst + make_Vec3fa(color.x,color.y,color.z,1.0f); // FIXME: use += operator
    float f = rcp(max(0.001f,dst->w));
    unsigned int r = (unsigned int) (255.0f * clamp(dst->x*f,0.0f,1.0f));
    unsigned int g = (unsigned int) (255.0f * clamp(dst->y*f,0.0f,1.0f));
    unsigned int b = (unsigned int) (255.0f * clamp(dst->z*f,0.0f,1.0f));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
}

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

  /* render frame */
  uniform const int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  uniform const int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  enableFilterDispatch = renderPixel == renderPixelStandard; 
  launch[numTilesX*numTilesY] renderTile(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); sync;
  enableFilterDispatch = false;
  rtcDebug();
}

/* called by the C++ code for cleanup */
export void device_cleanup ()
{
  rtcDeleteScene (g_scene);
  rtcExit();
}
