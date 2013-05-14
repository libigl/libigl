#ifndef IGL_MATERIAL_COLORS_H
#define IGL_MATERIAL_COLORS_H
// Define constant material colors for use with opengl glMaterialfv
// Most of these colors come from IGL publications
namespace igl
{
  // Gold/Silver used in BBW/MONO/STBS/FAST
  const float GOLD_AMBIENT[4] =   {  51.0/255.0, 43.0/255.0,33.3/255.0,1.0f };
  const float GOLD_DIFFUSE[4] =   { 255.0/255.0,228.0/255.0,58.0/255.0,1.0f };
  const float GOLD_SPECULAR[4] =  { 255.0/255.0,235.0/255.0,80.0/255.0,1.0f };
  const float SILVER_AMBIENT[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
  const float SILVER_DIFFUSE[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
  const float SILVER_SPECULAR[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
  // Blue/Cyan more similar to Jovan Popovic's blue than to Mario Botsch's blue
  const float CYAN_AMBIENT[4] =   {  59.0/255.0, 68.0/255.0,255.0/255.0,1.0f };
  const float CYAN_DIFFUSE[4] =   {  94.0/255.0,185.0/255.0,238.0/255.0,1.0f };
  const float CYAN_SPECULAR[4] =   { 163.0/255.0,221.0/255.0,255.0/255.0,1.0f };
  const float DENIS_PURPLE_DIFFUSE[4] =   { 80.0/255.0,64.0/255.0,255.0/255.0,1.0f };
  const float LADISLAV_ORANGE_DIFFUSE[4] = {1.0f, 125.0f / 255.0f, 19.0f / 255.0f, 0.0f};
  // FAST armadillos colors
  const float FAST_GREEN_DIFFUSE[4] = { 113.0f/255.0f, 239.0f/255.0f,  46.0f/255.0f, 1.0f};
  const float FAST_RED_DIFFUSE[4]   = { 255.0f/255.0f,  65.0f/255.0f,  46.0f/255.0f, 1.0f};
  const float FAST_BLUE_DIFFUSE[4]  = { 106.0f/255.0f, 106.0f/255.0f, 255.0f/255.0f, 1.0f};
  const float FAST_GRAY_DIFFUSE[4]  = { 150.0f/255.0f, 150.0f/255.0f, 150.0f/255.0f, 1.0f};
  // Basic colors
  const float WHITE_AMBIENT[4] =   { 255.0/255.0,255.0/255.0,255.0/255.0,1.0f };
  const float WHITE_DIFFUSE[4] =   { 255.0/255.0,255.0/255.0,255.0/255.0,1.0f };
  const float WHITE_SPECULAR[4] =  { 255.0/255.0,255.0/255.0,255.0/255.0,1.0f };
  const float BBW_POINT_COLOR[4] = {239./255.,213./255.,46./255.,255.0/255.0};
  const float BBW_LINE_COLOR[4] = {106./255.,106./255.,255./255.,255./255.};
  const float MIDNIGHT_BLUE_DIFFUSE[4]  = { 21.0f/255.0f, 27.0f/255.0f, 84.0f/255.0f, 1.0f};
  // Winding number colors
  const float EASTER_RED_DIFFUSE[4] = {0.603922,0.494118f,0.603922f,1.0f};
  const float WN_OPEN_BOUNDARY_COLOR[4] = {154./255.,0./255.,0./255.,1.0f};
  const float WN_NON_MANIFOLD_EDGE_COLOR[4] = {201./255., 51./255.,255./255.,1.0f};
}
#endif
