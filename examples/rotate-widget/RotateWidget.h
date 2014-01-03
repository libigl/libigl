#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
// 3D Rotate tool widget 
class RotateWidget
{
  public:
    static Eigen::Quaterniond axis_q[3];
    static Eigen::Vector3d view_direction(const int x, const int y);
    static Eigen::Vector3d view_direction(const Eigen::Vector3d & pos);
    Eigen::Vector3d pos;
    Eigen::Quaterniond rot,down_rot;
    Eigen::Vector2d down_xy,drag_xy,down_dir;
    Eigen::Vector3d udown,udrag;
    double outer_radius_on_screen;
    double outer_over_inner;
    enum DownType
    {
      DOWN_TYPE_X = 0,
      DOWN_TYPE_Y = 1,
      DOWN_TYPE_Z = 2,
      DOWN_TYPE_OUTLINE = 3,
      DOWN_TYPE_TRACKBALL = 4,
      DOWN_TYPE_NONE = 5,
      NUM_DOWN_TYPES = 6
    } down_type, selected_type;
    std::vector<Eigen::Vector3d> debug_points;
    RotateWidget();
    // Vector from origin to mouse click "Unprojected" onto plane with depth of
    // origin and scale to so that outer radius is 1
    // 
    // Inputs:
    //   x  mouse x position
    //   y  mouse y position
    // Returns vector
    Eigen::Vector3d unproject_onto(const int x, const int y);
    // Shoot ray from mouse click to sphere
    //
    // Inputs:
    //   x  mouse x position
    //   y  mouse y position
    // Outputs:
    //   hit  position of hit
    // Returns true only if there was a hit
    bool intersect(const int x, const int y, Eigen::Vector3d & hit);
    double unprojected_inner_radius();
    bool down(const int x, const int y);
    bool drag(const int x, const int y);
    bool up(const int x, const int y);
    bool is_down();
    void draw();
    void draw_guide();
};
