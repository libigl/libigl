// Main class of the Viewer

#ifndef IGL_VIEWER_H
#define IGL_VIEWER_H

#include <AntTweakBar.h>

#include <vector>
#include <string>
#include <cstdint>

#define IGL_MOD_SHIFT           0x0001
#define IGL_MOD_CONTROL         0x0002
#define IGL_MOD_ALT             0x0004
#define IGL_MOD_SUPER           0x0008

#ifdef ENABLE_XML_SERIALIZATION
  #include <igl/xml/XMLSerializer.h>
  #include <igl/xml/XMLSerialization.h>
#endif

#include <Eigen/Core>
#include <igl/viewer/OpenGL_shader.h>
#include <igl/viewer/ViewerData.h>
#include <igl/viewer/OpenGL_state.h>
#include <igl/viewer/ViewerPlugin.h>

namespace igl
{

  class Viewer
  {
  public:

    int launch(std::string filename = "");
    void init();

    class Options
    #ifdef ENABLE_XML_SERIALIZATION
    : public ::igl::XMLSerialization
    #endif
    {
    public:
      Options()
      #ifdef ENABLE_XML_SERIALIZATION
      : XMLSerialization("Options")
      #endif
      {};
      void InitSerialization();

      // Shape material
      float shininess;

      // Colors
      Eigen::Vector3f background_color;
      Eigen::Vector3f line_color;

      // Lighting
      Eigen::Vector3f light_position;
      float lighting_factor;

      // Trackball angle (quaternion)
      Eigen::Vector4f trackball_angle;

      // Model viewing parameters
      float model_zoom;
      Eigen::Vector3f model_translation;

      // Model viewing paramters (uv coordinates)
      float model_zoom_uv;
      Eigen::Vector3f model_translation_uv;

      // Camera parameters
      float camera_zoom;
      bool orthographic;
      Eigen::Vector3f camera_eye;
      Eigen::Vector3f camera_up;
      Eigen::Vector3f camera_center;
      float camera_view_angle;
      float camera_dnear;
      float camera_dfar;

      // Visualization options
      bool show_overlay;
      bool show_overlay_depth;
      bool show_texture;
      bool show_faces;
      bool show_lines;
      bool show_vertid;
      bool show_faceid;
      bool invert_normals;

      // Point size / line width
      float point_size;
      float line_width;

      // Enable per-face colors and normals
      bool face_based;

      // Animation
      bool is_animating;
      double animation_max_fps;
    };

    // Stores all the viewing options
    Options options;

    // Stores all the data that should be visualized
    igl::ViewerData data;

    // Stores the vbos indices and opengl related settings
    igl::OpenGL_state opengl;

    // List of registered plugins
    std::vector<ViewerPlugin*> plugins;
    void init_plugins();
    void shutdown_plugins();

    // Temporary data stored when the mouse button is pressed
    Eigen::Vector4f down_rotation;
    int current_mouse_x;
    int current_mouse_y;
    int down_mouse_x;
    int down_mouse_y;
    float down_mouse_z;
    Eigen::Vector3f down_translation;
    bool down;
    bool hack_never_moved;

    // Anttweak bar
    TwBar* bar;

    // Window size
    int width;
    int height;

    // Keep track of the global position of the scrollwheel
    float scroll_position;

    // Useful functions
    void compute_normals(); // Computes the normals of the mesh
    void uniform_colors(Eigen::Vector3d ambient, Eigen::Vector3d diffuse, Eigen::Vector3d specular); // assign uniform colors to all faces/vertices
    void grid_texture(); // Generate a default grid texture

    void clear_mesh();      // Clear the mesh data
    void align_camera_center(); // Adjust the view to see the entire model

    // Change the visualization mode, invalidating the cache if necessary
    void set_face_based(bool newvalue);

    // Helpers that can draw the most common meshes
    void set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
    void set_vertices(const Eigen::MatrixXd& V);
    void set_normals(const Eigen::MatrixXd& N);
    // Set the color of the mesh
    //
    // Inputs:
    //   C  #V|#F|1 by 3 list of colors
    void set_colors(const Eigen::MatrixXd &C);
    void set_uv(const Eigen::MatrixXd& UV);
    void set_uv(const Eigen::MatrixXd& UV_V, const Eigen::MatrixXi& UV_F);
    void set_texture(
                      const Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& R,
                      const Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& G,
                      const Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& B);

    void add_points(const Eigen::MatrixXd& P,  const Eigen::MatrixXd& C);
    // Sets edges given a list of edge vertices and edge indices. In constrast
    // to `add_edges` this will (purposefully) clober existing edges.
    //
    // Inputs:
    //   P  #P by 3 list of vertex positions
    //   E  #E by 2 list of edge indices into P
    //   C  #E|1 by 3 color(s)
    void set_edges (const Eigen::MatrixXd& P, const Eigen::MatrixXi& E, const Eigen::MatrixXd& C);
    void add_edges (const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2, const Eigen::MatrixXd& C);
    void add_label (const Eigen::VectorXd& P,  const std::string& str);

    // Save the OpenGL transformation matrices used for the previous rendering pass
    Eigen::Matrix4f view;
    Eigen::Matrix4f model;
    Eigen::Matrix4f proj;

    Eigen::Vector4f viewport;

    // UI Enumerations
    enum MouseButton {IGL_LEFT, IGL_MIDDLE, IGL_RIGHT};
    enum MouseMode { NOTHING, ROTATION, ZOOM, PAN, TRANSLATE} mouse_mode;
    enum KeyModifier { NO_KEY = TW_KMOD_NONE, SHIFT = TW_KMOD_SHIFT, CTRL =TW_KMOD_CTRL, ALT = TW_KMOD_ALT } key_modifier;

    Viewer();
    ~Viewer();

    // Mesh IO
    bool load_mesh_from_file(const char* mesh_file_name);
    bool save_mesh_to_file(const char* mesh_file_name);

    // Callbacks
    bool key_down(unsigned char key, int modifier);
    bool key_up(unsigned char key, int modifier);

    bool mouse_down(MouseButton button, int modifier);
    bool mouse_up(MouseButton button, int modifier);

    bool mouse_move(int mouse_x, int mouse_y);
    bool mouse_scroll(float delta_y);

    // Scene IO
    bool load_scene();
    bool save_scene();

    // Determines how much to zoom and shift such that the mesh fills the unit
    // box (centered at the origin)
    static void get_scale_and_shift_to_fit_mesh(
      const Eigen::MatrixXd& V,
      const Eigen::MatrixXi& F,
      float & zoom,
      Eigen::Vector3f& shift);

    // Draw everything
    void draw();

    // OpenGL context resize
    void resize(int w, int h);


    // C-style callbacks
    bool (*callback_pre_draw)(Viewer& viewer);
    bool (*callback_post_draw)(Viewer& viewer);
    bool (*callback_mouse_down)(Viewer& viewer, int button, int modifier);
    bool (*callback_mouse_up)(Viewer& viewer, int button, int modifier);
    bool (*callback_mouse_move)(Viewer& viewer, int mouse_x, int mouse_y);
    bool (*callback_mouse_scroll)(Viewer& viewer, float delta_y);
    bool (*callback_key_down)(Viewer& viewer, unsigned char key, int modifiers);
    bool (*callback_key_up)(Viewer& viewer, unsigned char key, int modifiers);

    // Pointers to per-callback data
    void* callback_pre_draw_data;
    void* callback_post_draw_data;
    void* callback_mouse_down_data;
    void* callback_mouse_up_data;
    void* callback_mouse_move_data;
    void* callback_mouse_scroll_data;
    void* callback_key_down_data;
    void* callback_key_up_data;


    /********* AntTweakBar callbacks *********/
    static void TW_CALL snap_to_canonical_quaternion_cb(void *clientData);
    static void TW_CALL save_scene_cb(void *clientData);
    static void TW_CALL load_scene_cb(void *clientData);
    static void TW_CALL open_dialog_mesh(void *clientData);
    static void TW_CALL align_camera_center_cb(void *clientData);
    static void TW_CALL set_face_based_cb(const void *param, void *clientData);
    static void TW_CALL get_face_based_cb(void *param, void *clientData);
    static void TW_CALL set_invert_normals_cb(const void *param, void *clientData);
    static void TW_CALL get_invert_normals_cb(void *param, void *clientData);
  public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

} // end namespace

#ifndef IGL_STATIC_LIBRARY
#  include "Viewer.cpp"
#endif

#endif
