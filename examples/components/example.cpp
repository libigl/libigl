#include <igl/C_STR.h>
#include <igl/Camera.h>
#include <igl/REDRUM.h>
#include <igl/components.h>
#include <igl/opengl/create_shader_program.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/get_seconds.h>
#include <igl/hsv_to_rgb.h>
#include <igl/opengl/init_render_to_texture.h>
#include <igl/jet.h>
#include <igl/per_face_normals.h>
#include <igl/randperm.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/rgb_to_hsv.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/snap_to_fixed_up.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/write_triangle_mesh.h>
#include <igl/anttweakbar/ReAntTweakBar.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#ifndef GLUT_WHEEL_UP
#define GLUT_WHEEL_UP    3
#endif
#ifndef GLUT_WHEEL_DOWN
#define GLUT_WHEEL_DOWN  4
#endif
#ifndef GLUT_WHEEL_RIGHT
#define GLUT_WHEEL_RIGHT 5
#endif
#ifndef GLUT_WHEEL_LEFT
#define GLUT_WHEEL_LEFT  6
#endif
#ifndef GLUT_ACTIVE_COMMAND
#define GLUT_ACTIVE_COMMAND 8
#endif

#include <ctime>
#include <string>
#include <vector>
#include <stack>
#include <iostream>

int cc_hover = -1;

Eigen::MatrixXd V;
Eigen::VectorXd Vmid,Vmin,Vmax;
double bbd = 1.0;
Eigen::MatrixXi F;
Eigen::VectorXi CC;
Eigen::MatrixXd N;
struct State
{
  igl::Camera camera;
  Eigen::VectorXf I;
  Eigen::Matrix<GLubyte,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> selected;
  GLuint mask_id;
} s;
std::string out_filename;

GLuint pick_tex = 0;
GLuint pick_fbo = 0;
GLuint pick_dfbo = 0;

// See README for descriptions
enum RotationType
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type;

enum CenterType
{
  CENTER_TYPE_ORBIT = 0,
  CENTER_TYPE_FPS  = 1,
  NUM_CENTER_TYPES = 2,
} center_type = CENTER_TYPE_ORBIT;


std::stack<State> undo_stack;
std::stack<State> redo_stack;

bool wireframe_visible = false;
bool fill_visible = true;

bool is_rotating = false;
int down_x,down_y;
igl::Camera down_camera;

bool is_animating = false;
double animation_start_time = 0;
double ANIMATION_DURATION = 0.5;
Eigen::Quaterniond animation_from_quat;
Eigen::Quaterniond animation_to_quat;

int width,height;
Eigen::Vector4f light_pos(-0.1,-0.1,0.9,0);

#define REBAR_NAME "temp.rbr"
igl::anttweakbar::ReTwBar rebar;

// Forward
void init_components();
void init_relative();

void push_undo()
{
  undo_stack.push(s);
  // Clear
  redo_stack = std::stack<State>();
}

void TW_CALL set_rotation_type(const void * value, void * clientData)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  const RotationType old_rotation_type = rotation_type;
  rotation_type = *(const RotationType *)(value);
  if(rotation_type == ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP &&
    old_rotation_type != ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP)
  {
    animation_from_quat = s.camera.m_rotation_conj;
    snap_to_fixed_up(animation_from_quat,animation_to_quat);
    // start animation
    animation_start_time = get_seconds();
    is_animating = true;
  }
}
void TW_CALL get_rotation_type(void * value, void *clientData)
{
  RotationType * rt = (RotationType *)(value);
  *rt = rotation_type;
}

void reshape(int width, int height)
{
  ::width = width;
  ::height = height;
  glViewport(0,0,width,height);
  // Send the new window size to AntTweakBar
  TwWindowSize(width, height);
  s.camera.m_aspect = (double)width/(double)height;
  igl::opengl::init_render_to_texture(width,height, pick_tex, pick_fbo, pick_dfbo);
  igl::opengl::report_gl_error("init_render_to_texture: ");
  glutPostRedisplay();
}

void push_scene()
{
  using namespace igl;
  using namespace std;
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  auto & camera = s.camera;
  glMultMatrixd(camera.projection().data());
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(
    camera.eye()(0), camera.eye()(1), camera.eye()(2),
    camera.at()(0), camera.at()(1), camera.at()(2),
    camera.up()(0), camera.up()(1), camera.up()(2));
  glScaled(2./bbd,2./bbd,2./bbd);
  glTranslated(-Vmid(0),-Vmid(1),-Vmid(2));
}

void pop_scene()
{
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void draw_mesh(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & N,
  const Eigen::VectorXf & S,
  const GLuint & S_loc)
{
  using namespace Eigen;
  using namespace std;
  static Matrix<float,Dynamic,3,RowMajor> VR,NR;
  static Matrix<int,Dynamic,3,RowMajor> FR;
  static Matrix<float,Dynamic,1,ColMajor> SR;
  static GLuint ibo,vbo,sbo,nbo;
  static bool scene_dirty = true;
  if(scene_dirty)
  {
    VR.resize(F.rows()*3,3);
    NR.resize(F.rows()*3,3);
    SR.resize(F.rows()*3,1);
    FR.resize(F.rows(),3);
    for(int f = 0;f<F.rows();f++)
    {
      for(int c = 0;c<3;c++)
      {
        VR.row(3*f+c) = V.row(F(f,c)).cast<float>();
        SR(3*f+c) = S(F(f,c));
        NR.row(3*f+c) = N.row(f).cast<float>();
        FR(f,c) = 3*f+c;
      }
    }

    glGenBuffers(1,&ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(GLuint)*FR.size(),FR.data(),GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glGenBuffers(1,&vbo);
    glGenBuffers(1,&nbo);
    glGenBuffers(1,&sbo);

    glBindBuffer(GL_ARRAY_BUFFER,vbo);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float)*VR.size(),VR.data(),GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER,nbo);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float)*NR.size(),NR.data(),GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER,sbo);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float)*SR.size(),SR.data(),GL_STATIC_DRAW);
    igl::opengl::report_gl_error("glBindBuffer: ");

    scene_dirty = false;
  }

  glEnableClientState(GL_VERTEX_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER,vbo);
  glVertexPointer(3,GL_FLOAT,0,0);
  glEnableClientState(GL_NORMAL_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER,nbo);   
  glNormalPointer(GL_FLOAT,0,0);

  glBindBuffer(GL_ARRAY_BUFFER,sbo);   
  glVertexAttribPointer(S_loc, 1, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(S_loc);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ibo);
  glDrawElements(GL_TRIANGLES,FR.size(),GL_UNSIGNED_INT,0);
  glBindBuffer(GL_ARRAY_BUFFER,0);
}

// Set up double-sided lights
void lights()
{
  using namespace std;
  using namespace Eigen;
  glEnable(GL_LIGHTING);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glEnable(GL_LIGHT0);
  float WHITE[4] = {1,1,1,1.};
  float BLACK[4] = {0.,0.,0.,1.};
  Vector4f pos = light_pos;
  glLightfv(GL_LIGHT0,GL_AMBIENT,BLACK);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT0,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT0,GL_POSITION,pos.data());
  //glEnable(GL_LIGHT1);
  //pos(0) *= -1;
  //pos(1) *= -1;
  //pos(2) *= -1;
  //glLightfv(GL_LIGHT1,GL_AMBIENT,BLACK);
  //glLightfv(GL_LIGHT1,GL_DIFFUSE,NEAR_BLACK);
  //glLightfv(GL_LIGHT1,GL_SPECULAR,BLACK);
  //glLightfv(GL_LIGHT1,GL_POSITION,pos.data());
}

template <int Rows, int Cols>
GLuint generate_1d_texture(
  const Eigen::Matrix<GLubyte,Rows,Cols,Eigen::RowMajor> & colors)
{
  assert(colors.cols() == 3 && "Seems colors.cols() must be 3");
  GLuint tex_id = 0;
  glGenTextures(1,&tex_id);
  glBindTexture(GL_TEXTURE_1D,tex_id);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage1D(GL_TEXTURE_1D, 0, colors.cols(),colors.rows(),
    0,GL_RGB, GL_UNSIGNED_BYTE,
    colors.data());
  igl::opengl::report_gl_error("glTexImage1D: ");
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  igl::opengl::report_gl_error("texture: ");
  return tex_id;
}
 
GLuint color_shader(const size_t max_ids, GLuint & scalar_loc, GLuint & tex_id)
{
  std::string vertex_shader = R"(
#version 120
attribute float scalar_in;
varying float scalar_out;
void main()
{
  scalar_out = scalar_in;
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
)";
  std::string fragment_shader = R"(
#version 120
varying float scalar_out;
uniform float cmin;
uniform float cmax;
uniform sampler1D color_map;
void main()
{
  float scalar_normalized = max(min((scalar_out-cmin)/(cmax-cmin),1.0),0.0);
  gl_FragColor = texture1D(color_map,scalar_normalized);
}
)";
  Eigen::Matrix<GLubyte,Eigen::Dynamic,3,Eigen::RowMajor> colors(max_ids,3);
  for(size_t id = 0;id<max_ids;id++)
  {
    size_t index = id;
    size_t re = (index)%(256*256);
    colors(id,0) = (index-re)/(256*256);
    index = re;
    re = index%(256);
    colors(id,1) = (index-re)/(256);
    colors(id,2) = re;
  }
  tex_id = generate_1d_texture(colors);
  return igl::opengl::create_shader_program(
    vertex_shader.c_str(), 
    fragment_shader.c_str(),
    {{"scalar_in",scalar_loc}}
    );
}


void display()
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  glClearColor(0.8,0.8,0.8,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if(is_animating)
  {
    double t = (get_seconds() - animation_start_time)/ANIMATION_DURATION;
    if(t > 1)
    {
      t = 1;
      is_animating = false;
    }
    Quaterniond q = animation_from_quat.slerp(t,animation_to_quat).normalized();
    auto & camera = s.camera;
    switch(center_type)
    {
      default:
      case CENTER_TYPE_ORBIT:
        camera.orbit(q.conjugate());
        break;
      case CENTER_TYPE_FPS:
        camera.turn_eye(q.conjugate());
        break;
    }
  }

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  lights();
  push_scene();


  const auto & color_components_shader = [](
     const GLuint scalar_loc,
     GLuint & tex_id)->GLuint
  {
  std::string vertex_shader = R"(
#version 120
attribute float scalar_in;
varying vec3 normal;
varying float scalar_out;
void main()
{
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
  normal = normalize(gl_NormalMatrix * gl_Normal);
  scalar_out = scalar_in;
}
)";
  std::string fragment_shader = R"(
#version 120
varying vec3 normal;
varying float scalar_out;
uniform float cmin;
uniform float cmax;
uniform float cc_hover;
uniform sampler1D color_map;
uniform sampler1D selected_mask;
void main()
{
  float scalar_normalized = max(min((scalar_out-cmin)/(cmax-cmin),1.0),0.0);
  vec4 texture_color = texture1D(color_map,scalar_normalized);
  bool is_selected = texture1D(selected_mask,scalar_normalized).x > 0.5;
  const vec4 selected_color = vec4(1,0.2,0.2,1);
  if(scalar_out==cc_hover)
  {
    texture_color = 0.5*(texture_color + selected_color);
  }
  if(is_selected)
  {
    texture_color = selected_color;
  }
  const float num_lights = 1.0;
  vec4 diffuse = (1.0/num_lights)*(gl_LightSource[0].diffuse);
  vec4 ambient = vec4(0,0,0,0);
  ambient += (1.0/num_lights)*(gl_FrontMaterial.ambient * gl_LightSource[0].ambient);
  ambient += (1.0/num_lights)*(gl_LightModel.ambient * gl_FrontMaterial.ambient);
  vec4 color = ambient;
  // Phong
  vec3 lightDir = normalize(vec3(gl_LightSource[0].position));
  vec3 halfVector = gl_LightSource[0].halfVector.xyz;
  vec3 n = normalize(normal);
  float NdotL = max(abs(dot(n.xyz,lightDir)), 0.0);  
  vec4 specular = vec4(0.0,0.0,0.0,0.0);
  if (NdotL > 0.0) {
      color += diffuse * NdotL;
      vec3 halfV = normalize(halfVector);
      float NdotHV = max(abs(dot(n,halfV)),0.0);
      specular += gl_FrontMaterial.specular * gl_LightSource[0].specular * pow(NdotHV, gl_FrontMaterial.shininess);
  }
  gl_FragColor = color * texture_color + specular;
}

)";

    typedef Matrix<GLubyte,64,3,RowMajor> Matrix64_3_R_ubyte;
    typedef Matrix<float,64,3,RowMajor> Matrix64_3_R_float;
    Matrix64_3_R_ubyte colors;
    {
      Matrix64_3_R_float rgb = (Matrix64_3_R_ubyte()<<
        255,   0,   0,
        255,  24,   0,
        255,  48,   0,
        255,  72,   0,
        255,  96,   0,
        255, 120,   0,
        255, 143,   0,
        255, 167,   0,
        255, 191,   0,
        255, 215,   0,
        255, 239,   0,
        247, 255,   0,
        223, 255,   0,
        199, 255,   0,
        175, 255,   0,
        151, 255,   0,
        128, 255,   0,
        104, 255,   0,
         80, 255,   0,
         56, 255,   0,
         32, 255,   0,
          8, 255,   0,
          0, 255,  16,
          0, 255,  40,
          0, 255,  64,
          0, 255,  88,
          0, 255, 112,
          0, 255, 135,
          0, 255, 159,
          0, 255, 183,
          0, 255, 207,
          0, 255, 231,
          0, 255, 255,
          0, 231, 255,
          0, 207, 255,
          0, 183, 255,
          0, 159, 255,
          0, 135, 255,
          0, 112, 255,
          0,  88, 255,
          0,  64, 255,
          0,  40, 255,
          0,  16, 255,
          8,   0, 255,
         32,   0, 255,
         56,   0, 255,
         80,   0, 255,
        104,   0, 255,
        128,   0, 255,
        151,   0, 255,
        175,   0, 255,
        199,   0, 255,
        223,   0, 255,
        247,   0, 255,
        255,   0, 239,
        255,   0, 215,
        255,   0, 191,
        255,   0, 167,
        255,   0, 143,
        255,   0, 120,
        255,   0,  96,
        255,   0,  72,
        255,   0,  48,
        255,   0,  24).finished().cast<float>()/255.f;

      Matrix64_3_R_float H;
      rgb_to_hsv(rgb,H);
      H.col(1) *= 0.1;
      H.col(2) = (H.col(2).array() + 0.1*(1.-H.col(2).array())).eval();
      hsv_to_rgb(H,rgb);
      colors = (rgb*255.).cast<GLubyte>();
    }

    tex_id = generate_1d_texture(colors);

    GLuint prog_id = igl::opengl::create_shader_program(
      vertex_shader.c_str(), 
      fragment_shader.c_str(),
      {{"scalar_in",scalar_loc}}
      );
    igl::opengl::report_gl_error("create_shader_program: ");
    return prog_id;
  };
  static GLuint scalar_loc = 1;
  static GLuint tex_id = 0;
  static GLuint color_components_prog = 
    color_components_shader(scalar_loc,tex_id);

  // Set material properties
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,(const GLfloat[]){1,1,1,1});
  if(wireframe_visible)
  {
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    if(fill_visible)
    {
      glColor3f(0,0,0);
      glUseProgram(0);
      draw_mesh(V,F,N,s.I,scalar_loc);
    }else
    {
      glUseProgram(color_components_prog);
      igl::opengl::report_gl_error("UseProgram: ");
      draw_mesh(V,F,N,s.I,scalar_loc);
    }
  }

  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glUseProgram(color_components_prog);
    igl::opengl::report_gl_error("use: ");
  glUniform1f(glGetUniformLocation(color_components_prog,"cmin"),s.I.minCoeff());
  glUniform1f(glGetUniformLocation(color_components_prog,"cmax"),s.I.maxCoeff());
  //glUniform1f(glGetUniformLocation(color_components_prog,"cc_selected"),cc_selected);
  glUniform1f(glGetUniformLocation(color_components_prog,"cc_hover"),cc_hover);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_1D, tex_id);
  glUniform1i(glGetUniformLocation(color_components_prog,"color_map"),0);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_1D, s.mask_id);
  glUniform1i(glGetUniformLocation(color_components_prog,"selected_mask"),1);

    igl::opengl::report_gl_error("unif: ");
  if(fill_visible)
  {
    glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
    glPolygonOffset(1.0, 0);
  }
  draw_mesh(V,F,N,s.I,scalar_loc);
  glPopAttrib();
  glUseProgram(0);

  // Draw a nice floor
  glPushMatrix();
  const double floor_offset =
    -2./bbd*(V.col(1).maxCoeff()-Vmid(1));
  glTranslated(0,floor_offset,0);
  const float GREY[4] = {0.5,0.5,0.6,1.0};
  const float DARK_GREY[4] = {0.2,0.2,0.3,1.0};
  igl::opengl2::draw_floor(GREY,DARK_GREY);
  glPopMatrix();

  pop_scene();

  TwDraw();
  glutSwapBuffers();
  if(is_animating)
  {
    glutPostRedisplay();
  }
}

void mouse_wheel(int wheel, int direction, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  if(wheel == 0 && TwMouseMotion(mouse_x, viewport[3] - mouse_y))
  {
    static double mouse_scroll_y = 0;
    const double delta_y = 0.125*direction;
    mouse_scroll_y += delta_y;
    TwMouseWheel(mouse_scroll_y);
    return;
  }

  auto & camera = s.camera;
  switch(center_type)
  {
    case CENTER_TYPE_ORBIT:
      if(wheel==0)
      {
        // factor of zoom change
        double s = (1.-0.01*direction);
        //// FOV zoom: just widen angle. This is hardly ever appropriate.
        //camera.m_angle *= s;
        //camera.m_angle = min(max(camera.m_angle,1),89);
        camera.push_away(s);
      }else
      {
        // Dolly zoom:
        camera.dolly_zoom((double)direction*1.0);
      }
      break;
    default:
    case CENTER_TYPE_FPS:
      // Move `eye` and `at`
      camera.dolly((wheel==0?Vector3d(0,0,1):Vector3d(-1,0,0))*0.1*direction);
      break;
  }
  glutPostRedisplay();
}

bool pick(const int x, const int y, int & cc_selected)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  static GLuint scalar_loc = 1;
  static GLuint tex_id = 0;
  static const size_t max_ids = s.I.maxCoeff()+1;
  static GLuint color_shader_prog = color_shader(max_ids,scalar_loc,tex_id);
  const int pick_s = 0;
  const int pick_w = pick_s;
  GLint old_vp[4];
  glGetIntegerv(GL_VIEWPORT,old_vp);
  const double pick_ratio = double(pick_w)/double(old_vp[2]);
  // ceil, cause might otherwise round down to 0
  const int pick_h = ceil(double(old_vp[3])*pick_ratio);
  glViewport(
    x-pick_w,
    old_vp[3]-y-pick_h,2*pick_w+1,2*pick_h+1);
  glMatrixMode(GL_PROJECTION);
  Matrix4d proj;
  glGetDoublev(GL_PROJECTION_MATRIX,proj.data());
  glPushMatrix();
  glLoadIdentity();
  gluPickMatrix(
    x,
    old_vp[3]-y, 
    pick_w*2+1,
    pick_h*2+1,
    old_vp);
  glMultMatrixd(proj.data());
  glMatrixMode(GL_MODELVIEW);
  // Activate color shader
  glUseProgram(color_shader_prog);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,pick_fbo);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT,pick_dfbo);
  // Clear screen
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glEnable(GL_TEXTURE_1D);
  glBindTexture(GL_TEXTURE_1D, tex_id);
  glUniform1f(glGetUniformLocation(color_shader_prog,"cmin"),s.I.minCoeff());
  glUniform1f(glGetUniformLocation(color_shader_prog,"cmax"),s.I.maxCoeff());
  draw_mesh(V,F,N,s.I,scalar_loc);
  glPopAttrib();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glViewport(old_vp[0],old_vp[1],old_vp[2],old_vp[3]);

  Matrix<GLubyte,1,4> pixel;
  glReadPixels(x,old_vp[3]-y,1,1,GL_RGBA,GL_UNSIGNED_BYTE,pixel.data());

  glUseProgram(0);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT,0);

  if(pixel(3) == 0)
  {
    cc_selected = -1;
    return false;
  }
  cc_selected = pixel(0)*256*256+pixel(1)*256+pixel(2);
  return true;
}

void regenerate_mask()
{
  if(glIsTexture(s.mask_id))
  {
    glDeleteTextures(1,&s.mask_id);
  }
  s.mask_id = generate_1d_texture(s.selected);
}

void mouse(int glutButton, int glutState, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  bool tw_using = TwEventMouseButtonGLUT(glutButton,glutState,mouse_x,mouse_y);
  int mod = glutGetModifiers();
  switch(glutButton)
  {
    case GLUT_RIGHT_BUTTON:
    {
      switch(glutState)
      {
        case 1:
          // up
          glutSetCursor(GLUT_CURSOR_INHERIT);
          is_rotating = false;
          break;
        case 0:
          glutSetCursor(GLUT_CURSOR_CYCLE);
          // collect information for trackball
          is_rotating = true;
          down_camera = s.camera;
          down_x = mouse_x;
          down_y = mouse_y;
        break;
      }
      break;
    }
    case GLUT_LEFT_BUTTON:
    {
      switch(glutState)
      {
        case 1:
          // up
          glutSetCursor(GLUT_CURSOR_INHERIT);
          is_rotating = false;
          break;
        case 0:
          if(!tw_using)
          {
            push_scene();
            int cc_selected=-1;
            if(pick(mouse_x,mouse_y,cc_selected))
            {
              push_undo();
              if(!(mod & GLUT_ACTIVE_SHIFT))
              {
                s.selected.setConstant(0);
              }
              s.selected(cc_selected,0) = 255;
              regenerate_mask();
            }else
            {
              glutSetCursor(GLUT_CURSOR_CYCLE);
              // collect information for trackball
              is_rotating = true;
              down_camera = s.camera;
              down_x = mouse_x;
              down_y = mouse_y;
            }
            pop_scene();
          }
        break;
      }
      break;
    }
    // Scroll down
    case GLUT_WHEEL_DOWN:
    {
      mouse_wheel(0,-1,mouse_x,mouse_y);
      break;
    }
    // Scroll up
    case GLUT_WHEEL_UP:
    {
      mouse_wheel(0,1,mouse_x,mouse_y);
      break;
    }
    // Scroll left
    case GLUT_WHEEL_LEFT:
    {
      mouse_wheel(1,-1,mouse_x,mouse_y);
      break;
    }
    // Scroll right
    case GLUT_WHEEL_RIGHT:
    {
      mouse_wheel(1,1,mouse_x,mouse_y);
      break;
    }
  }
  glutPostRedisplay();
}

void mouse_move(int mouse_x, int mouse_y)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  bool tw_using = TwMouseMotion(mouse_x,mouse_y);
  push_scene();
  pick(mouse_x,mouse_y,cc_hover);
  pop_scene();
  glutPostRedisplay();
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  bool tw_using = TwMouseMotion(mouse_x,mouse_y);

  if(is_rotating)
  {
    glutSetCursor(GLUT_CURSOR_CYCLE);
    Quaterniond q;
    auto & camera = s.camera;
    switch(rotation_type)
    {
      case ROTATION_TYPE_IGL_TRACKBALL:
      {
        // Rotate according to trackball
        igl::trackball<double>(
          width,
          height,
          2.0,
          down_camera.m_rotation_conj.coeffs().data(),
          down_x,
          down_y,
          mouse_x,
          mouse_y,
          q.coeffs().data());
          break;
      }
      case ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP:
      {
        // Rotate according to two axis valuator with fixed up vector
        two_axis_valuator_fixed_up(
          width, height,
          2.0,
          down_camera.m_rotation_conj,
          down_x, down_y, mouse_x, mouse_y,
          q);
        break;
      }
      default:
        break;
    }
    switch(center_type)
    {
      default:
      case CENTER_TYPE_ORBIT:
        camera.orbit(q.conjugate());
        break;
      case CENTER_TYPE_FPS:
        camera.turn_eye(q.conjugate());
        break;
    }
  }
  glutPostRedisplay();
}

void init_relative()
{
  using namespace Eigen;
  using namespace igl;
  per_face_normals(V,F,N);
  Vmax = V.colwise().maxCoeff();
  Vmin = V.colwise().minCoeff();
  Vmid = 0.5*(Vmax + Vmin);
  bbd = (Vmax-Vmin).norm();
}

void init_components()
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  components(F,CC);
  s.I = CC.cast<float>();
  s.selected = Matrix<GLubyte,Dynamic,Dynamic>::Zero(s.I.maxCoeff()+1,3);
  cout<<"s.selected: "<<s.selected.rows()<<endl;
  regenerate_mask();
}

void undo()
{
  using namespace std;
  if(!undo_stack.empty())
  {
    redo_stack.push(s);
    s = undo_stack.top();
    undo_stack.pop();
  }
  regenerate_mask();
}

void redo()
{
  using namespace std;
  if(!redo_stack.empty())
  {
    undo_stack.push(s);
    s = redo_stack.top();
    redo_stack.pop();
  }
  regenerate_mask();
}

bool save(const std::string & out_filename)
{
  using namespace std;
  using namespace igl;
  if(write_triangle_mesh(out_filename,V,F))
  {
    cout<<GREENGIN("Saved mesh to `"<<out_filename<<"` successfully.")<<endl;
    return true;
  }else
  {
    cout<<REDRUM("Failed to save mesh to `"<<out_filename<<"`.")<<endl;
    return false;
  }
}

void TW_CALL saveCB(void * /*clientData*/)
{
  save(out_filename);
}

void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  int mod = glutGetModifiers();
  switch(key)
  {
    // ESC
    case char(27):
      rebar.save(REBAR_NAME);
    // ^C
    case char(3):
      exit(0);
    case 'z':
    case 'Z':
      if(mod & GLUT_ACTIVE_COMMAND)
      {
        if(mod & GLUT_ACTIVE_SHIFT)
        {
          redo();
        }else
        {
          undo();
        }
      }else
      {
        Quaterniond q;
        snap_to_canonical_view_quat(s.camera.m_rotation_conj,1.0,q);
        switch(center_type)
        {
          default:
          case CENTER_TYPE_ORBIT:
            s.camera.orbit(q.conjugate());
            break;
          case CENTER_TYPE_FPS:
            s.camera.turn_eye(q.conjugate());
            break;
        }
      }
      break;
    case 'u':
        mouse_wheel(0, 1,mouse_x,mouse_y);
        break;
    case 'j':
        mouse_wheel(0,-1,mouse_x,mouse_y);
        break;
    default:
      if(!TwEventKeyboardGLUT(key,mouse_x,mouse_y))
      {
        cout<<"Unknown key command: "<<key<<" "<<int(key)<<endl;
      }
  }

  glutPostRedisplay();
}

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  string filename = "../shared/truck.obj";
  switch(argc)
  {
    case 3:
      out_filename = argv[2];
    case 2:
      // Read and prepare mesh
      filename = argv[1];
      break;
    default:
      cerr<<"Usage:"<<endl<<"    ./example input.obj (output.obj)"<<endl;
      cout<<endl<<"Opening default mesh..."<<endl;
      break;
  }

  // print key commands
  cout<<"[Click] and [drag]  Rotate model using trackball."<<endl;
  cout<<"[Z,z]               Snap rotation to canonical view."<<endl;
  cout<<"[Command+Z]         Undo."<<endl;
  cout<<"[Shift+Command+Z]   Redo."<<endl;
  cout<<"[^C,ESC]            Exit."<<endl;

  read_triangle_mesh(filename,V,F);


  // Init glut
  glutInit(&argc,argv);
  if( !TwInit(TW_OPENGL, NULL) )
  {
    // A fatal error occured
    fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
    return 1;
  }
  // Create a tweak bar
  rebar.TwNewBar("bar");
  TwDefine("bar label='Components' size='200 550' text=light alpha='200' color='68 68 68'");
  rebar.TwAddVarRW("camera_rotation", TW_TYPE_QUAT4D,
    s.camera.m_rotation_conj.coeffs().data(), "open readonly=true");
  TwType RotationTypeTW = igl::anttweakbar::ReTwDefineEnumFromString("RotationType",
    "igl_trackball,two-axis-valuator-fixed-up");
  rebar.TwAddVarCB( "rotation_type", RotationTypeTW,
    set_rotation_type,get_rotation_type,NULL,"keyIncr=] keyDecr=[");
  TwType CenterTypeTW = igl::anttweakbar::ReTwDefineEnumFromString("CenterType","orbit,fps");
  rebar.TwAddVarRW("center_type", CenterTypeTW,&center_type,
    "keyIncr={ keyDecr=}");

  rebar.TwAddVarRW("wireframe_visible",TW_TYPE_BOOLCPP,&wireframe_visible,"key=l");
  rebar.TwAddVarRW("fill_visible",TW_TYPE_BOOLCPP,&fill_visible,"key=f");
  if(out_filename != "")
  {
    rebar.TwAddButton("save",
      saveCB,NULL,
      C_STR("label='Save to `"<<out_filename<<"`' "<<
      "key=s"));
  }
  rebar.load(REBAR_NAME);


  animation_from_quat = Quaterniond(1,0,0,0);
  s.camera.m_rotation_conj = animation_from_quat;
  animation_start_time = get_seconds();

  glutInitDisplayString( "rgba depth double samples>=8");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT)/2.0);
  glutCreateWindow("components");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc(mouse_move);

  init_components();
  init_relative();
  regenerate_mask();

  std::cout<<"OpenGL version: "<<glGetString(GL_VERSION)<<std::endl;
  glutMainLoop();

  return 0;
}
