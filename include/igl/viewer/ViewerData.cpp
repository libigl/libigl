#include "ViewerData.h"

void igl::ViewerData::InitSerialization()
{
  #ifdef ENABLE_XML_SERIALIZATION
  xmlSerializer->Add(V,"V");
  xmlSerializer->Add(F,"F");
  xmlSerializer->Add(F_normals,"F_normals");

  xmlSerializer->Add(F_material_ambient,"F_material_ambient");
  xmlSerializer->Add(F_material_diffuse,"F_material_diffuse");
  xmlSerializer->Add(F_material_specular,"F_material_specular");

  xmlSerializer->Add(V_normals,"V_normals");
  xmlSerializer->Add(V_material_ambient,"V_material_ambient");
  xmlSerializer->Add(V_material_diffuse,"V_material_diffuse");
  xmlSerializer->Add(V_material_specular,"V_material_specular");

  xmlSerializer->Add(V_uv,"V_uv");
  xmlSerializer->Add(F_uv,"F_uv");
  xmlSerializer->Add(texture_R,"texture_R");
  xmlSerializer->Add(texture_G,"texture_G");
  xmlSerializer->Add(texture_B,"texture_B");
  xmlSerializer->Add(lines,"lines");
  xmlSerializer->Add(points,"points");

  xmlSerializer->Add(labels_positions,"labels_positions");
  xmlSerializer->Add(labels_strings,"labels_strings");
  #endif
}
