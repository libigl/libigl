#include "lens_flare.h"

#include "C_STR.h"
#include "unproject.h"
#include "project.h"
#include "shine_textures.h"
#include "flare_textures.h"

// http://www.opengl.org/archives/resources/features/KilgardTechniques/LensFlare/glflare.c

static void setup_texture(
  const uint8_t * texture, 
  const int width,
  const int height,
  GLuint texobj,
  GLenum minFilter, GLenum maxFilter)
{
  glBindTexture(GL_TEXTURE_2D, texobj);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, maxFilter);
  glTexImage2D(GL_TEXTURE_2D, 0, 1, width, height, 0,
    GL_LUMINANCE, GL_UNSIGNED_BYTE, texture);
}

void igl::lens_flare_load_textures(
  std::vector<GLuint> & shine_id,
  std::vector<GLuint> & flare_id)
{
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  shine_id.resize(10);
  glGenTextures(10,&shine_id[0]);
  for (int i = 0; i < (int)shine_id.size(); i++) {
    setup_texture(
      SHINE_TEXTURES[i],
      SHINE_TEXTURE_WIDTHS[i],
      SHINE_TEXTURE_HEIGHTS[i],
      shine_id[i], GL_LINEAR, GL_LINEAR);
  }
  flare_id.resize(6);
  glGenTextures(6,&flare_id[0]);
  for (int i = 0; i < (int)flare_id.size(); i++) {
    setup_texture(
      FLARE_TEXTURES[i],
      FLARE_TEXTURE_WIDTHS[i],
      FLARE_TEXTURE_HEIGHTS[i],
      flare_id[i], GL_LINEAR, GL_LINEAR);
  }
}

void igl::lens_flare_create(
  const float * A,
  const float * B,
  const float * C,
  std::vector<igl::Flare> & flares)
{
  using namespace igl;
  flares.resize(12);
  /* Shines */
  flares[0] = Flare(-1, 1.0f, 0.1f, C, 1.0);
  flares[1] = Flare(-1, 1.0f, 0.15f, B, 1.0);
  flares[2] = Flare(-1, 1.0f, 0.35f, A, 1.0);

  /* Flares */
  flares[3] =  Flare(2, 1.3f, 0.04f, A, 0.6);
  flares[4] =  Flare(3, 1.0f, 0.1f, A, 0.4);
  flares[5] =  Flare(1, 0.5f, 0.2f, A, 0.3);
  flares[6] =  Flare(3, 0.2f, 0.05f, A, 0.3);
  flares[7] =  Flare(0, 0.0f, 0.04f, A, 0.3);
  flares[8] =  Flare(5, -0.25f, 0.07f, A, 0.5);
  flares[9] =  Flare(5, -0.4f, 0.02f, A, 0.6);
  flares[10] = Flare(5, -0.6f, 0.04f, A, 0.4);
  flares[11] = Flare(5, -1.0f, 0.03f, A, 0.2);
}

void igl::lens_flare_draw(
  const std::vector<igl::Flare> & flares,
  const std::vector<GLuint> & shine_ids,
  const std::vector<GLuint> & flare_ids,
  const Eigen::Vector3f & light,
  const float near_clip,
  int & shine_tic)
{
  bool ot2,ob,odt;
  int obsa,obda;
  ot2 = glIsEnabled(GL_TEXTURE_2D);
  ob = glIsEnabled(GL_BLEND);
  odt = glIsEnabled(GL_DEPTH_TEST);
  glGetIntegerv(GL_BLEND_SRC_ALPHA,&obsa);
  glGetIntegerv(GL_BLEND_DST_ALPHA,&obda);

  glEnable(GL_TEXTURE_2D);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ONE);

  using namespace Eigen;
  using namespace igl;

  // Assuming scene is pushed
  Vector3f win_o,from,at;
  project(Vector3f(0,0,0),win_o);
  win_o(2) = -1e26;
  unproject(win_o,from);
  win_o(2) = 1;
  unproject(win_o,at);

  Vector3f view_dir, tmp, light0, light_dir, position, dx, dy,
    center, axis, sx, sy;
  GLfloat dot, global_scale = 1.0;
  int i;

  /* view_dir = normalize(at-from) */
  view_dir =  at -  from;
  view_dir.normalize();

  /* center = from + near_clip * view_dir */
  tmp =  near_clip* view_dir;
  center =  from +  tmp;

  /* light_dir = normalize(light-from) */
  light_dir =  light -  from;
  light_dir.normalize();

  /* light0 = from + dot(light,view_dir)*near_clip*light_dir */
  dot = light_dir.dot( view_dir);
  tmp =  near_clip / dot* light_dir;
  light0 =  from + light_dir;

  /* axis = light - center */
  axis =  light0 -  center;
  dx =  axis;

  /* dx = normalize(axis) */
  dx.normalize();

  /* dy = cross(dx,view_dir) */
  dy =  dx.cross( view_dir);


  for (i = 0; i < (int)flares.size(); i++)
  {
    sx =  flares[i].scale * global_scale* dx;
    sy =  flares[i].scale * global_scale* dy;

    glColor3fv(flares[i].color);
    if (flares[i].type < 0) {
      glBindTexture(GL_TEXTURE_2D, shine_ids[shine_tic]);
      shine_tic = (shine_tic + 1) % shine_ids.size();
    } else {
      glBindTexture(GL_TEXTURE_2D, flare_ids[flares[i].type]);
    }

    /* position = center + flare[i].loc * axis */
    tmp =  flares[i].loc* axis;
    position =  center +  tmp;

    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0);
    tmp =  position +  sx;
    tmp =  tmp +  sy;
    glVertex3fv(tmp.data());

    glTexCoord2f(1.0, 0.0);
    tmp =  position -  sx;
    tmp =  tmp +  sy;
    glVertex3fv(tmp.data());

    glTexCoord2f(1.0, 1.0);
    tmp =  position -  sx;
    tmp =  tmp -  sy;
    glVertex3fv(tmp.data());

    glTexCoord2f(0.0, 1.0);
    tmp =  position +  sx;
    tmp =  tmp -  sy;
    glVertex3fv(tmp.data());
    glEnd();
  }
  ot2?glEnable(GL_TEXTURE_2D):glDisable(GL_TEXTURE_2D);
  ob?glEnable(GL_BLEND):glDisable(GL_BLEND);
  odt?glEnable(GL_DEPTH_TEST):glDisable(GL_DEPTH_TEST);
  glBlendFunc(obsa,obda);
}
