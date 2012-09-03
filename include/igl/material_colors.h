#ifndef IGL_MATERIAL_COLORS_H
#define IGL_MATERIAL_COLORS_H
// Define som constant material colors for use with opengl glMaterialfv
namespace igl
{
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
  const float WHITE_AMBIENT[4] =   { 255.0/255.0,255.0/255.0,255.0/255.0,1.0f };
  const float WHITE_DIFFUSE[4] =   { 255.0/255.0,255.0/255.0,255.0/255.0,1.0f };
  const float WHITE_SPECULAR[4] =  { 255.0/255.0,255.0/255.0,255.0/255.0,1.0f };
}
#endif
