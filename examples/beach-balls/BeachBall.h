#ifndef BEACHBALL_H
#define BEACHBALL_H
class BeachBall
{
  public:
    static double * scene_rot;
  public:
    // Rigid transformation: rotation as quaternion
    double r[4];
    // Rigid transformation: translation as quaternion
    double t[4];
    double radius;
  protected:
    bool is_hover;
    bool is_down;
    // trackball
    bool trackball_on;
    // Rigid transformation at down
    double down_r[4];
    double down_t[4];
    // query position at down
    int down_x;
    int down_y;
    // modelview matrix at draw 
    double mv[16];
  public:
    BeachBall();
    // Push the latest modelview matrix seen at draw time onto the stack
    void pushmv() const;
    // Pop the stack
    void popmv() const;
    // Push rigid transformation onto the stack
    void push() const;
    // Pop the stack
    void pop() const;
    // Draw and cache the modelview matrix
    void draw();
    // Return whether (x,y) in screenspace is "inside" (hitting)
    bool in(const int x,const int y) const;
    // Called on passive mouse move to (x,y)
    //
    // Return whether (x,y) is hovering over, set is_hover accordingly
    bool hover(const int x,const int y);
    // Called on mouse down at (x,y)
    //
    // Return whether (x,y) is down, set is_down accordingly
    bool down(const int x,const int y);
    // Called on mouse drag to (x,y)
    //
    // Returns whether still is_down
    bool drag(const int x,const int y);
    // Called on right mouse down at (x,y)
    //
    // Return whether (x,y) is down, set is_down accordingly
    bool right_down(const int x,const int y);
    // Called on mouse up at (x,y)
    //
    // Returns whether still is_down
    bool up(const int x,const int y);
};
#endif
