// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "sys/platform.h"
#include "sys/filename.h"
#include "image/image.h"
#include "lexers/streamfilters.h"
#include "lexers/parsestream.h"
#include "loaders/loaders.h"
#include "glutdisplay.h"

namespace embree
{
  /******************************************************************************/
  /*                                  State                                     */
  /******************************************************************************/

  /* camera settings */
  Vec3f g_camPos    = Vec3f(0.0f,0.0f,0.0f);
  Vec3f g_camLookAt = Vec3f(1.0f,0.0f,0.0f);
  Vec3f g_camUp     = Vec3f(0,1,0);
  float g_camFieldOfView = 64.0f;
  float g_camRadius = 0.0f;

  /* rendering device and global handles */
  Device* g_device = NULL;
  Handle<Device::RTRenderer> g_renderer = NULL;
  Handle<Device::RTToneMapper> g_tonemapper = NULL;
  Handle<Device::RTFrameBuffer> g_frameBuffer = NULL;
  Handle<Device::RTImage> g_backplate = NULL;
  Handle<Device::RTScene> g_render_scene;
  std::vector<Handle<Device::RTPrimitive> > g_prims;

  void clearGlobalObjects() {
    g_renderer = null;
    g_tonemapper = null;
    g_frameBuffer = null;
    g_backplate = null;
    g_render_scene = null;
    g_prims.clear();
    delete g_device; 
    g_device = NULL;
  }

  /* rendering settings */
  std::string g_accel = "default";
  std::string g_tri = "default";

  int g_depth = -1;                       //!< recursion depth
  int g_spp = 1;                          //!< samples per pixel for ordinary rendering

  /* output settings */
  int g_numBuffers = 1;                   //!< number of buffers of the framebuffer
  bool g_rendered = false;                //!< set to true after rendering
  int g_refine = 1;                       //!< refinement mode
  float g_gamma = 1.0f;
  bool g_vignetting = true;
  bool g_fullscreen = false;
  size_t g_width = 512, g_height = 512;
  size_t g_numThreads = 0;

  /* regression testing mode */
  bool g_regression = false;

  
  /******************************************************************************/
  /*                            Object Creation                                 */
  /******************************************************************************/

  Handle<Device::RTCamera> createCamera(const AffineSpace3f& space)
  {
    /*! pinhole camera */
    if (g_camRadius == 0.0f) 
    {
      Handle<Device::RTCamera> camera = g_device->rtNewCamera("pinhole");
      g_device->rtSetTransform(camera, "local2world", copyToArray(space));
      g_device->rtSetFloat1(camera, "angle", g_camFieldOfView);
      g_device->rtSetFloat1(camera, "aspectRatio", float(g_width) / float(g_height));
      g_device->rtCommit(camera);
      return camera;
    }
    /*! depth of field camera */
    else
    {
      Handle<Device::RTCamera> camera = g_device->rtNewCamera("depthoffield");
      g_device->rtSetTransform(camera, "local2world", copyToArray(space));
      g_device->rtSetFloat1(camera, "angle", g_camFieldOfView);
      g_device->rtSetFloat1(camera, "aspectRatio", float(g_width) / float(g_height));
      g_device->rtSetFloat1(camera, "lensRadius", g_camRadius);
      g_device->rtSetFloat1(camera, "focalDistance", length(g_camLookAt - g_camPos));
      g_device->rtCommit(camera);
      return camera;
    }
  }

  Handle<Device::RTScene> createScene() {
    return g_device->rtNewScene((g_accel+" "+g_tri).c_str(), (Device::RTPrimitive*)(g_prims.size() == 0 ? NULL : &g_prims[0]), g_prims.size());
  }

  void createGlobalObjects()
  {
    g_renderer = g_device->rtNewRenderer("pathtracer");
    if (g_depth >= 0) g_device->rtSetInt1(g_renderer, "maxDepth", g_depth);
    g_device->rtSetInt1(g_renderer, "sampler.spp", g_spp);
    g_device->rtCommit(g_renderer);
    
    g_tonemapper = g_device->rtNewToneMapper("default");
    g_device->rtSetFloat1(g_tonemapper, "gamma", g_gamma);
    g_device->rtCommit(g_tonemapper);

    g_frameBuffer = g_device->rtNewFrameBuffer("RGB_FLOAT32", g_width, g_height, 1);
    g_backplate = NULL;
  }

  /******************************************************************************/
  /*                      Command line parsing                                  */
  /******************************************************************************/

  static void parseDebugRenderer(Ref<ParseStream> cin, const FileName& path)
  {
    g_renderer = g_device->rtNewRenderer("debug");
    if (g_depth >= 0) g_device->rtSetInt1(g_renderer, "maxDepth", g_depth);
    g_device->rtSetInt1(g_renderer, "sampler.spp", g_spp);

    if (cin->peek() != "{") goto finish;
    cin->drop();

    while (cin->peek() != "}") {
      std::string tag = cin->getString();
      cin->force("=");
      if (tag == "depth") g_device->rtSetInt1(g_renderer, "maxDepth", cin->getInt());
      else std::cout << "unknown tag \"" << tag << "\" in debug renderer parsing" << std::endl;
    }
    cin->drop();

  finish:
    g_device->rtCommit(g_renderer);
  }

  static void parsePathTracer(Ref<ParseStream> cin, const FileName& path)
  {
    g_renderer = g_device->rtNewRenderer("pathtracer");
    if (g_depth >= 0) g_device->rtSetInt1(g_renderer, "maxDepth", g_depth);
    g_device->rtSetInt1(g_renderer, "sampler.spp", g_spp);
    if (g_backplate) g_device->rtSetImage(g_renderer, "backplate", g_backplate);

    if (cin->peek() != "{") goto finish;
    cin->drop();

    while (cin->peek() != "}") {
      std::string tag = cin->getString();
      cin->force("=");
      if      (tag == "depth"          ) g_device->rtSetInt1  (g_renderer, "maxDepth"       , cin->getInt()  );
      else if (tag == "spp"            ) g_device->rtSetInt1  (g_renderer, "sampler.spp"    , cin->getInt()  );
      else if (tag == "minContribution") g_device->rtSetFloat1(g_renderer, "minContribution", cin->getFloat());
      else if (tag == "backplate"      ) g_device->rtSetImage (g_renderer, "backplate", loadImage(path + cin->getFileName(), g_device));
      else std::cout << "unknown tag \"" << tag << "\" in debug renderer parsing" << std::endl;
    }
    cin->drop();

  finish:
    g_device->rtCommit(g_renderer);
  }

  static void displayMode()
  {
    if (!g_renderer) throw std::runtime_error("no renderer set");
    AffineSpace3f camSpace = AffineSpace3f::lookAtPoint(g_camPos, g_camLookAt, g_camUp);
    float speed = 0.02f * length(g_camLookAt - g_camPos);
    GLUTDisplay(camSpace, speed, createScene());
    g_rendered = true;
  }

  static void outputMode(const FileName& fileName)
  {
    if (!g_renderer) throw std::runtime_error("no renderer set");

    /* render image */
    Handle<Device::RTCamera> camera = createCamera(AffineSpace3f::lookAtPoint(g_camPos, g_camLookAt, g_camUp));
    Handle<Device::RTScene> scene = createScene();
    g_device->rtSetInt1(g_renderer, "showprogress", 1);
    g_device->rtCommit(g_renderer);
    g_device->rtRenderFrame(g_renderer, camera, scene, g_tonemapper, g_frameBuffer, 0);
    for (int i=0; i < g_numBuffers; i++) g_device->rtSwapBuffers(g_frameBuffer);

    /* store to disk */
    Vec4f* rgba = (Vec4f*) g_device->rtMapFrameBuffer(g_frameBuffer);
	if (sizeof(Col3f) == sizeof(Vec4f)) {
      Ref<Image3f> image = new Image3f(g_width, g_height, (Col3f*)rgba);
      storeImage(image.cast<Image>(), fileName);
	} else {
      Col3f* ptr = new Col3f[g_width*g_height];
	  for (size_t i=0; i<g_width*g_height; i++)
		ptr[i] = Col3f(rgba[i].x,rgba[i].y,rgba[i].z);
      Ref<Image3f> image = new Image3f(g_width, g_height, (Col3f*)ptr);
      storeImage(image.cast<Image>(), fileName);
      delete[] ptr;
	}
    g_device->rtUnmapFrameBuffer(g_frameBuffer);
    g_rendered = true;
  }

  static void parseCommandLine(Ref<ParseStream> cin, const FileName& path)
  {
    while (true)
    {
      std::string tag = cin->getString();
      if (tag == "") return;

      /* set number of local threads to use */
      else if (tag == "-connect") {
        clearGlobalObjects();
        g_device = Device::rtCreateDevice(cin->getString().c_str(),g_numThreads);
        createGlobalObjects();
      }

      /* set number of local threads to use */
      else if (tag == "-threads") {
        g_numThreads = cin->getInt();
        clearGlobalObjects();
        g_device = Device::rtCreateDevice("default",g_numThreads);
        createGlobalObjects();
      }

      /* parse command line parameters from a file */
      else if (tag == "-c") {
        FileName file = path + cin->getFileName();
        parseCommandLine(new ParseStream(new LineCommentFilter(file, "#")), file.path());
      }

      /* read model from file */
      else if (tag == "-i") {
        std::vector<Device::RTPrimitive> prims = loadScene(path + cin->getFileName(), g_device);
        g_prims.insert(g_prims.end(), prims.begin(), prims.end());
      }
      
      /* triangulated sphere */
      else if (tag == "-trisphere")
      {
        Handle<Device::RTShape> sphere = g_device->rtNewShape("sphere");
        const Vec3f P = cin->getVec3f();
        g_device->rtSetFloat3(sphere, "P", P.x, P.y, P.z);
        g_device->rtSetFloat1(sphere, "r", cin->getFloat());
        g_device->rtSetInt1(sphere, "numTheta", cin->getInt());
        g_device->rtSetInt1(sphere, "numPhi", cin->getInt());
        g_device->rtCommit(sphere);

        Handle<Device::RTMaterial> material = g_device->rtNewMaterial("matte");
        g_device->rtSetFloat3(material, "reflection", 1.0f, 0.0f, 0.0f);
        g_device->rtCommit(material);
        g_prims.push_back(g_device->rtNewShapePrimitive(sphere, material, NULL));
      }

      /* ambient light source */
      else if (tag == "-ambientlight") {
        Handle<Device::RTLight> light = g_device->rtNewLight("ambientlight");
        const Col3f L = cin->getCol3f();
        g_device->rtSetFloat3(light, "L", L.r, L.g, L.b);
        g_device->rtCommit(light);
        g_prims.push_back(g_device->rtNewLightPrimitive(light, NULL));
      }

      /* point light source */
      else if (tag == "-pointlight") {
        Handle<Device::RTLight> light = g_device->rtNewLight("pointlight");
        const Vec3f P = cin->getVec3f();
        const Col3f I = cin->getCol3f();
        g_device->rtSetFloat3(light, "P", P.x, P.y, P.z);
        g_device->rtSetFloat3(light, "I", I.r, I.g, I.b);
        g_device->rtCommit(light);
        g_prims.push_back(g_device->rtNewLightPrimitive(light, NULL));
      }

      /* directional light source */
      else if (tag == "-directionallight" || tag == "-dirlight") {
        Handle<Device::RTLight> light = g_device->rtNewLight("directionallight");
        const Vec3f D = cin->getVec3f();
        const Col3f E = cin->getCol3f();
        g_device->rtSetFloat3(light, "D", D.x, D.y, D.z);
        g_device->rtSetFloat3(light, "E", E.r, E.g, E.b);
        g_device->rtCommit(light);
        g_prims.push_back(g_device->rtNewLightPrimitive(light, NULL));
      }

      /* distant light source */
      else if (tag == "-distantlight") {
        Handle<Device::RTLight> light = g_device->rtNewLight("distantlight");
        const Vec3f D = cin->getVec3f();
        const Col3f L = cin->getCol3f();
        g_device->rtSetFloat3(light, "D", D.x, D.y, D.z);
        g_device->rtSetFloat3(light, "L", L.r, L.g, L.b);
        g_device->rtSetFloat1(light, "halfAngle", cin->getFloat());
        g_device->rtCommit(light);
        g_prims.push_back(g_device->rtNewLightPrimitive(light, NULL));
      }

      /* triangular light source */
      else if (tag == "-trianglelight") {
        Vec3f P = cin->getVec3f();
        Vec3f U = cin->getVec3f();
        Vec3f V = cin->getVec3f();
        Vec3f L = cin->getVec3f();
        
        Handle<Device::RTLight> light = g_device->rtNewLight("trianglelight");
        g_device->rtSetFloat3(light, "v0", P.x, P.y, P.z);
        g_device->rtSetFloat3(light, "v1", P.x + U.x, P.y + U.y, P.z + U.z);
        g_device->rtSetFloat3(light, "v2", P.x + V.x, P.y + V.y, P.z + V.z);
        g_device->rtSetFloat3(light, "L",  L.x, L.y, L.z);
        g_device->rtCommit(light);
        g_prims.push_back(g_device->rtNewLightPrimitive(light, NULL));
      }

      /* quad light source */
      else if (tag == "-quadlight")
      {
        Vec3f P = cin->getVec3f();
        Vec3f U = cin->getVec3f();
        Vec3f V = cin->getVec3f();
        Vec3f L = cin->getVec3f();

        Handle<Device::RTLight> light0 = g_device->rtNewLight("trianglelight");
        g_device->rtSetFloat3(light0, "v0", P.x + U.x + V.x, P.y + U.y + V.y, P.z + U.z + V.z);
        g_device->rtSetFloat3(light0, "v1", P.x + U.x, P.y + U.y, P.z + U.z);
        g_device->rtSetFloat3(light0, "v2", P.x, P.y, P.z);
        g_device->rtSetFloat3(light0, "L",  L.x, L.y, L.z);
        g_device->rtCommit(light0);
        g_prims.push_back(g_device->rtNewLightPrimitive(light0, NULL));

        Handle<Device::RTLight> light1 = g_device->rtNewLight("trianglelight");
        g_device->rtSetFloat3(light1, "v0", P.x + U.x + V.x, P.y + U.y + V.y, P.z + U.z + V.z);
        g_device->rtSetFloat3(light1, "v1", P.x, P.y, P.z);
        g_device->rtSetFloat3(light1, "v2", P.x + V.x, P.y + V.y, P.z + V.z);
        g_device->rtSetFloat3(light1, "L",  L.x, L.y, L.z);
        g_device->rtCommit(light1);
        g_prims.push_back(g_device->rtNewLightPrimitive(light1, NULL));
      }

      /* HDRI light source */
      else if (tag == "-hdrilight")
      {
        Handle<Device::RTLight> light = g_device->rtNewLight("hdrilight");
        const Col3f L = cin->getCol3f();
        g_device->rtSetFloat3(light, "L", L.r, L.g, L.b);
        g_device->rtSetImage(light, "image", loadImage(path + cin->getFileName(), g_device));
        g_device->rtCommit(light);
        g_prims.push_back(g_device->rtNewLightPrimitive(light, NULL));
      }

      /* parse camera parameters */
      else if (tag == "-vp")     g_camPos         = Vec3f(cin->getVec3f());
      else if (tag == "-vi")     g_camLookAt      = Vec3f(cin->getVec3f());
      else if (tag == "-vd")     g_camLookAt      = g_camPos + cin->getVec3f();
      else if (tag == "-vu")     g_camUp          = cin->getVec3f();
      else if (tag == "-angle")  g_camFieldOfView = cin->getFloat();
      else if (tag == "-fov")    g_camFieldOfView = cin->getFloat();
      else if (tag == "-radius") g_camRadius      = cin->getFloat();

      /* frame buffer size */
      else if (tag == "-size") {
        g_width = cin->getInt();
        g_height = cin->getInt();
        g_frameBuffer = g_device->rtNewFrameBuffer("RGB_FLOAT32", g_width, g_height, 1);
      }

      /* full screen mode */
      else if (tag == "-fullscreen") g_fullscreen = true;

      /* refine rendering when not moving */
      else if (tag == "-refine") g_refine = cin->getInt();
      
      /* acceleration structure to use */
      else if (tag == "-accel") g_accel = cin->getString();

      /* triangle representation to use */
      else if (tag == "-tri") g_tri = cin->getString();

      /* set renderer */
      else if (tag == "-renderer")
      {
        std::string renderer = cin->getString();
        if      (renderer == "debug"     ) parseDebugRenderer(cin, path);
        else if (renderer == "pt"        ) parsePathTracer(cin, path);
        else if (renderer == "pathtracer") parsePathTracer(cin, path);
        else throw std::runtime_error("unknown renderer: " + renderer);
      }

      /* set gamma */
      else if (tag == "-gamma") {
        g_device->rtSetFloat1(g_tonemapper, "gamma", g_gamma = cin->getFloat());
        g_device->rtCommit(g_tonemapper);
      }

      /* set gamma */
      else if (tag == "-vignetting") {
        g_device->rtSetBool1(g_tonemapper, "vignetting", g_vignetting = cin->getInt());
        g_device->rtCommit(g_tonemapper);
      }

      /* set recursion depth */
      else if (tag == "-depth") {
        g_device->rtSetInt1(g_renderer, "maxDepth", g_depth = cin->getInt());
        g_device->rtCommit(g_renderer);
      }

      /* set samples per pixel */
      else if (tag == "-spp") {
        g_device->rtSetInt1(g_renderer, "sampler.spp", g_spp = cin->getInt());
        g_device->rtCommit(g_renderer);
      }

      /* set the backplate */
      else if (tag == "-backplate") {
        g_device->rtSetImage(g_renderer, "backplate", g_backplate = loadImage(path + cin->getFileName(), g_device));
        g_device->rtCommit(g_renderer);
      }

      /* render frame */
      else if (tag == "-o") outputMode(path + cin->getFileName());

      /* display image */
      else if (tag == "-display") displayMode();

      /* regression testing */
      else if (tag == "-regression")
      {
        g_refine = false;
        g_regression = true;
        GLUTDisplay(AffineSpace3f::lookAtPoint(g_camPos, g_camLookAt, g_camUp), 0.01f);
      }

      /* print version information */
      else if (tag == "-version") {
        std::cout << "Embree renderer version 1.1" << std::endl;
        exit(1);
      }

      /* print help */
      else if (tag == "-h" || tag == "-?" || tag == "-help" || tag == "--help")
      {
        std::cout << std::endl;
        std::cout << "Embree Version 1.1" << std::endl;
        std::cout << std::endl;
        std::cout << "  usage: embree -i model.obj -renderer debug" << std::endl;
        std::cout << "         embree -i model.obj -renderer pathtracer -o out.tga" << std::endl;
        std::cout << "         embree -c model.ecs" << std::endl;
        std::cout << std::endl;
        std::cout << "-renderer [debug,pathtracer]" << std::endl;
        std::cout << "  Sets the renderer to use." << std::endl;
        std::cout << std::endl;
        std::cout << "-c file" << std::endl;
        std::cout << "  Parses command line parameters from file." << std::endl;
        std::cout << std::endl;
        std::cout << "-i file" << std::endl;
        std::cout << "  Loads a scene from file." << std::endl;
        std::cout << std::endl;
        std::cout << "-o file" << std::endl;
        std::cout << "  Renders and outputs the image to the file." << std::endl;
        std::cout << std::endl;
        std::cout << "-display" << std::endl;
        std::cout << "  Interactively displays the rendering into a window." << std::endl;
        std::cout << std::endl;
        std::cout << "-vp x y z" << std::endl;
        std::cout << "  Sets camera position to the location (x,y,z)." << std::endl;
        std::cout << std::endl;
        std::cout << "-vi x y z" << std::endl;
        std::cout << "  Sets camera lookat point to the location (x,y,z)." << std::endl;
        std::cout << std::endl;
        std::cout << "-vd x y z" << std::endl;
        std::cout << "  Sets camera viewing direction to (x,y,z)." << std::endl;
        std::cout << std::endl;
        std::cout << "-vu x y z" << std::endl;
        std::cout << "  Sets camera up direction to (x,y,z)." << std::endl;
        std::cout << std::endl;
        std::cout << "-fov angle" << std::endl;
        std::cout << "  Sets camera field of view in y direction to angle." << std::endl;
        std::cout << std::endl;
        std::cout << "-size width height" << std::endl;
        std::cout << "  Sets the width and height of image to render." << std::endl;
        std::cout << std::endl;
        std::cout << "-fullscreen" << std::endl;
        std::cout << "  Enables full screen display mode." << std::endl;
        std::cout << std::endl;
        std::cout << "-accel [bvh2,bvh4,bvh4.spatial]" << std::endl;
        std::cout << "  Sets the spatial index structure to use." << std::endl;
        std::cout << std::endl;
        std::cout << "-gamma v" << std::endl;
        std::cout << "  Sets gamma correction to v (only pathtracer)." << std::endl;
        std::cout << std::endl;
        std::cout << "-depth i" << std::endl;
        std::cout << "  Sets the recursion depth to i (default 16)" << std::endl;
        std::cout << std::endl;
        std::cout << "-spp i" << std::endl;
        std::cout << "  Sets the number of samples per pixel to i (default 1) (only pathtracer)." << std::endl;
        std::cout << std::endl;
        std::cout << "-backplate" << std::endl;
        std::cout << "  Sets a high resolution back ground image. (default none) (only pathtracer)." << std::endl;
        std::cout << std::endl;

        std::cout << "-ambientlight r g b" << std::endl;
        std::cout << "  Creates an ambient light with intensity (r,g,b)." << std::endl;
        std::cout << std::endl;
        std::cout << "-pointlight px py pz r g b" << std::endl;
        std::cout << "  Creates a point light with intensity (r,g,b) at position (px,py,pz)." << std::endl;
        std::cout << std::endl;
        std::cout << "-distantlight dx dy dz r g b halfAngle" << std::endl;
        std::cout << "  Creates a distant sun light with intensity (r,g,b) shining into " << std::endl;
        std::cout << "  direction (dx,dy,dz) from the cone spanned by halfAngle." << std::endl;
        std::cout << std::endl;
        std::cout << "-trianglelight px py pz ux uy uz vx vy vz r g b" << std::endl;
        std::cout << "  Creates a triangle-light with intensity (r,g,b) spanned by the point " << std::endl;
        std::cout << "  (px,py,pz) and the vectors (vx,vy,vz) and (ux,uy,uz)." << std::endl;
        std::cout << std::endl;
        std::cout << "-quadlight px py pz ux uy uz vx vy vz r g b" << std::endl;
        std::cout << "  Creates a quad-light with intensity (r,g,b) spanned by the point " << std::endl;
        std::cout << "  (px,py,pz) and the vectors (vx,vy,vz) and (ux,uy,uz)." << std::endl;
        std::cout << std::endl;
        std::cout << "-hdrilight r g b file" << std::endl;
        std::cout << "  Creates a high dynamic range environment light from the image " << std::endl;
        std::cout << "  file. The intensities are multiplies by (r,g,b)." << std::endl;
        std::cout << std::endl;
        std::cout << "-trisphere px py pz r theta phi" << std::endl;
        std::cout << "  Creates a triangulated sphere with radius r at location (px,py,pz) " << std::endl;
        std::cout << "  and triangulation rates theta and phi." << std::endl;
        std::cout << std::endl;
        std::cout << "-[no]refine" << std::endl;
        std::cout << "  Enables (default) or disables the refinement display mode." << std::endl;
        std::cout << std::endl;
        std::cout << "-regression" << std::endl;
        std::cout << "  Runs a stress test of the system." << std::endl;
        std::cout << std::endl;
        std::cout << "-version" << std::endl;
        std::cout << "  Prints version number." << std::endl;
        std::cout << std::endl;
        std::cout << "-h, -?, -help, --help" << std::endl;
        std::cout << "  Prints this help." << std::endl;
        exit(1);
      }

      /* skip unknown command line parameter */
      else {
        std::cerr << "unknown command line parameter: " << tag << " ";
        while (cin->peek() != "" && cin->peek()[0] != '-') std::cerr << cin->getString() << " ";
        std::cerr << std::endl;
      }
    }
  }

  /* main function in embree namespace */
  int main(int argc, char** argv) 
  {
    /*! create embree device */
    g_device = Device::rtCreateDevice("default",g_numThreads);
    createGlobalObjects();
    
    /*! parse command line */  
    parseCommandLine(new ParseStream(new CommandLineStream(argc, argv)), FileName());

    /*! if we did no render yet but have loaded a scene, switch to display mode */
    if (!g_rendered && g_prims.size()) displayMode();

    return(0);

  }
}

/******************************************************************************/
/*                               Main Function                                */
/******************************************************************************/

int main(int argc, char** argv)
{
  try {
    return embree::main(argc, argv);
  }
  catch (const std::exception& e) {
    std::cout << "Error: " << e.what() << std::endl;
    return 1;
  }
}
