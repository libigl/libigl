#ifndef IGL_VIEWER_SERIALIZATION_H
#define IGL_VIEWER_SERIALIZATION_H

#ifdef IGL_VIEWER_WITH_NANOGUI_SERIALIZATION

#include <igl/serialize.h>

namespace igl {
  namespace serialization {

    // Viewer members
    IGL_INLINE void serialization(bool s,igl::viewer::Viewer& obj,std::vector<char>& buffer)
    {
      obj.data_buffer[obj.active_data_id] = obj.data;

      SERIALIZE_MEMBER(core);
      SERIALIZE_MEMBER(data_buffer);
      SERIALIZE_MEMBER(data_ids);
      SERIALIZE_MEMBER(active_data_id);
    }

    template<>
    IGL_INLINE void serialize(const igl::viewer::Viewer& obj,std::vector<char>& buffer)
    {
      serialization(true,const_cast<igl::viewer::Viewer&>(obj),buffer);
    }

    template<>
    IGL_INLINE void deserialize(igl::viewer::Viewer& obj,const std::vector<char>& buffer)
    {
      serialization(false,obj,const_cast<std::vector<char>&>(buffer));

      obj.data = obj.data_buffer[obj.active_data_id];
    }

    // ViewerData members
    IGL_INLINE void serialization(bool s,igl::viewer::ViewerData& obj,std::vector<char>& buffer)
    {
      SERIALIZE_MEMBER(model_translation);
      SERIALIZE_MEMBER(model);

      SERIALIZE_MEMBER(object_scale);

      SERIALIZE_MEMBER(V);
      SERIALIZE_MEMBER(F);

      SERIALIZE_MEMBER(F_normals);
      SERIALIZE_MEMBER(F_material_ambient);
      SERIALIZE_MEMBER(F_material_diffuse);
      SERIALIZE_MEMBER(F_material_specular);

      SERIALIZE_MEMBER(V_normals);
      SERIALIZE_MEMBER(V_material_ambient);
      SERIALIZE_MEMBER(V_material_diffuse);
      SERIALIZE_MEMBER(V_material_specular);

      SERIALIZE_MEMBER(V_uv);
      SERIALIZE_MEMBER(F_uv);
      SERIALIZE_MEMBER(model_translation_uv);
      SERIALIZE_MEMBER(model_zoom_uv);

      SERIALIZE_MEMBER(texture_R);
      SERIALIZE_MEMBER(texture_G);
      SERIALIZE_MEMBER(texture_B);

      SERIALIZE_MEMBER(lines);
      SERIALIZE_MEMBER(points);

      SERIALIZE_MEMBER(labels_positions);
      SERIALIZE_MEMBER(labels_strings);
      SERIALIZE_MEMBER(labels_colors);

      SERIALIZE_MEMBER(dirty);

      SERIALIZE_MEMBER(face_based);

      SERIALIZE_MEMBER(show_overlay);
      SERIALIZE_MEMBER(show_overlay_depth);
      SERIALIZE_MEMBER(show_texture);
      SERIALIZE_MEMBER(show_faces);
      SERIALIZE_MEMBER(show_lines);
      SERIALIZE_MEMBER(show_vertid);
      SERIALIZE_MEMBER(show_faceid);
      SERIALIZE_MEMBER(invert_normals);
      SERIALIZE_MEMBER(depth_test);
      SERIALIZE_MEMBER(visible);

      SERIALIZE_MEMBER(point_size);
      SERIALIZE_MEMBER(line_width);
    }

    template<>
    IGL_INLINE void serialize(const igl::viewer::ViewerData& obj,std::vector<char>& buffer)
    {
      serialization(true,const_cast<igl::viewer::ViewerData&>(obj),buffer);
    }

    template<>
    IGL_INLINE void deserialize(igl::viewer::ViewerData& obj,const std::vector<char>& buffer)
    {
      serialization(false,obj,const_cast<std::vector<char>&>(buffer));
      obj.dirty = igl::viewer::ViewerData::DIRTY_ALL;
    }

    // ViewerCore members
    IGL_INLINE void serialization(bool s,igl::viewer::ViewerCore& obj,std::vector<char>& buffer)
    {
      SERIALIZE_MEMBER(background_color);
      SERIALIZE_MEMBER(line_color);

      SERIALIZE_MEMBER(shininess);
      SERIALIZE_MEMBER(light_position);
      SERIALIZE_MEMBER(lighting_factor);

      SERIALIZE_MEMBER(trackball_angle);
      SERIALIZE_MEMBER(rotation_type);
      SERIALIZE_MEMBER(global_translation);

      SERIALIZE_MEMBER(camera_zoom);
      SERIALIZE_MEMBER(orthographic);
      SERIALIZE_MEMBER(camera_eye);
      SERIALIZE_MEMBER(camera_up);
      SERIALIZE_MEMBER(camera_center);
      SERIALIZE_MEMBER(camera_view_angle);
      SERIALIZE_MEMBER(camera_dnear);
      SERIALIZE_MEMBER(camera_dfar);

      SERIALIZE_MEMBER(is_animating);
      SERIALIZE_MEMBER(animation_max_fps);

      SERIALIZE_MEMBER(viewport);

      SERIALIZE_MEMBER(view);
      SERIALIZE_MEMBER(proj);
    }

    template<>
    IGL_INLINE void serialize(const igl::viewer::ViewerCore& obj,std::vector<char>& buffer)
    {
      serialization(true,const_cast<igl::viewer::ViewerCore&>(obj),buffer);
    }

    template<>
    IGL_INLINE void deserialize(igl::viewer::ViewerCore& obj,const std::vector<char>& buffer)
    {
      serialization(false,obj,const_cast<std::vector<char>&>(buffer));
    }
  }
}

#endif

#endif