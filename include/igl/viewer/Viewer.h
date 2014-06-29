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

namespace igl
{
  // Forward declaration of the viewer_plugin class
  class Viewer_plugin;

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

    enum DirtyFlags
    {
      DIRTY_NONE           = 0x0000,
      DIRTY_POSITION       = 0x0001,
      DIRTY_UV             = 0x0002,
      DIRTY_NORMAL         = 0x0004,
      DIRTY_AMBIENT        = 0x0008,
      DIRTY_DIFFUSE        = 0x0010,
      DIRTY_SPECULAR       = 0x0020,
      DIRTY_TEXTURE        = 0x0040,
      DIRTY_FACE           = 0x0080,
      DIRTY_MESH           = 0x00FF,
      DIRTY_OVERLAY_LINES  = 0x0100,
      DIRTY_OVERLAY_POINTS = 0x0200,
      DIRTY_ALL            = 0x03FF
    };

    class OpenGL_state
    {
    public:
      typedef unsigned int GLuint;

      GLuint vao_mesh;
      GLuint vao_overlay_lines;
      GLuint vao_overlay_points;
      OpenGL_shader shader_mesh;
      OpenGL_shader shader_overlay_lines;
      OpenGL_shader shader_overlay_points;

      GLuint vbo_V; // Vertices of the current mesh (#V x 3)
      GLuint vbo_V_uv; // UV coordinates for the current mesh (#V x 2)
      GLuint vbo_V_normals; // Vertices of the current mesh (#V x 3)
      GLuint vbo_V_ambient; // Ambient material  (#V x 3)
      GLuint vbo_V_diffuse; // Diffuse material  (#V x 3)
      GLuint vbo_V_specular; // Specular material  (#V x 3)

      GLuint vbo_F; // Faces of the mesh (#F x 3)
      GLuint vbo_tex; // Texture

      GLuint vbo_lines_F;         // Indices of the line overlay
      GLuint vbo_lines_V;         // Vertices of the line overlay
      GLuint vbo_lines_V_colors;  // Color values of the line overlay
      GLuint vbo_points_F;        // Indices of the point overlay
      GLuint vbo_points_V;        // Vertices of the point overlay
      GLuint vbo_points_V_colors; // Color values of the point overlay

      // Temporary copy of the content of each VBO
      Eigen::MatrixXf V_vbo;
      Eigen::MatrixXf V_normals_vbo;
      Eigen::MatrixXf V_ambient_vbo;
      Eigen::MatrixXf V_diffuse_vbo;
      Eigen::MatrixXf V_specular_vbo;
      Eigen::MatrixXf V_uv_vbo;
      Eigen::MatrixXf lines_V_vbo;
      Eigen::MatrixXf lines_V_colors_vbo;
      Eigen::MatrixXf points_V_vbo;
      Eigen::MatrixXf points_V_colors_vbo;

      int tex_u;
      int tex_v;
      Eigen::Matrix<char,Eigen::Dynamic,1> tex;

      Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> F_vbo;
      Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> lines_F_vbo;
      Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> points_F_vbo;

      // Marks dirty buffers that need to be uploaded to OpenGL
      uint32_t dirty;

      // Create a new set of OpenGL buffer objects
      void init();

      // Update contents from a 'Data' instance
      void set_data(const igl::ViewerData &data, bool face_based, bool invert_normals);

      // Bind the underlying OpenGL buffer objects for subsequent mesh draw calls
      void bind_mesh();

      /// Draw the currently buffered mesh (either solid or wireframe)
      void draw_mesh(bool solid);

      // Bind the underlying OpenGL buffer objects for subsequent line overlay draw calls
      void bind_overlay_lines();

      /// Draw the currently buffered line overlay
      void draw_overlay_lines();

      // Bind the underlying OpenGL buffer objects for subsequent point overlay draw calls
      void bind_overlay_points();

      /// Draw the currently buffered point overlay
      void draw_overlay_points();

      // Release the OpenGL buffer objects
      void free();
    };

    // Stores all the viewing options
    Options options;

    // Stores all the data that should be visualized
    igl::ViewerData data;

    // Stores the vbos indices and opengl related settings
    OpenGL_state opengl;

    // List of registered plugins
    std::vector<Viewer_plugin*> plugins;
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


    // Init opengl shaders and VBOs
    void init_opengl();
    void free_opengl();

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


  // Abstract class for plugins
  // All plugins MUST have this class as their parent and implement all the callbacks
  // For an example of a basic plugins see plugins/skeleton.h
  //
  // Return value of callbacks: returning true to any of the callbacks tells Preview3D that the event has been
  // handled and that it should not be passed to other plugins or to other internal functions of Preview3D

  class Viewer_plugin
  #ifdef ENABLE_XML_SERIALIZATION
  : public ::igl::XMLSerialization
  #endif
  {
  public:
    Viewer_plugin()
    #ifdef ENABLE_XML_SERIALIZATION
    : XMLSerialization("dummy")
    #endif
    {plugin_name = "dummy";};

    ~Viewer_plugin(){};

    // This function is called when the viewer is initialized (no mesh will be loaded at this stage)
    virtual void init(igl::Viewer *_viewer)
    {
      viewer = _viewer;
    }

    // This function is called before shutdown
    virtual void shutdown()
    {
    }

    // This function is called before a mesh is loaded
    virtual bool load(std::string filename)
    {
      return false;
    }

    // This function is called before a mesh is saved
    virtual bool save(std::string filename)
    {
      return false;
    }

    // Runs immediately after a new mesh had been loaded.
    virtual bool post_load()
    {
      return false;
    }

    // This function is called before the draw procedure of Preview3D
    virtual bool pre_draw()
    {
      return false;
    }

    // This function is called after the draw procedure of Preview3D
    virtual bool post_draw()
    {
      return false;
    }

    // This function is called when the mouse button is pressed
    // - button can be GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON or GLUT_RIGHT_BUTTON
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool mouse_down(int button, int modifier)
    {
      return false;
    }

    // This function is called when the mouse button is released
    // - button can be GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON or GLUT_RIGHT_BUTTON
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool mouse_up(int button, int modifier)
    {
      return false;
    }

    // This function is called every time the mouse is moved
    // - mouse_x and mouse_y are the new coordinates of the mouse pointer in screen coordinates
    virtual bool mouse_move(int mouse_x, int mouse_y)
    {
      return false;
    }

    // This function is called every time the scroll wheel is moved
    // Note: this callback is not working with every glut implementation
    virtual bool mouse_scroll(float delta_y)
    {
      return false;
    }

    // This function is called when a keyboard key is pressed
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool key_down(unsigned char key, int modifiers)
    {
      return false;
    }

    // This function is called when a keyboard key is release
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool key_up(unsigned char key, int modifiers)
    {
      return false;
    }

    // Priority of the plugin (use only positive numbers, negative are reserved for internal use)
    // The plugins will be initialized in increasing priority
    virtual int priority()
    {
      return 0;
    }

    std::string plugin_name;
  protected:
    // Pointer to the main Preview3D class
    Viewer *viewer;
  public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

} // end namespace

#ifndef IGL_STATIC_LIBRARY
#  include "Viewer.cpp"
#endif

#endif
