#ifndef UI_HPP
#define UI_HPP

#include "Utilities.hpp"
#include "command_line.hpp"

#include <cmath>
#define _USE_MATH_DEFINES

#include <GL/glut.h>

#include <Eigen/Eigen>

using namespace Eigen;

struct Camera {
    // Camera transformation information
    Vec3f position;
    Vec3f axis;
    float angle;

    // Frustum information
    float near;
    float far;
    float fov;
    float aspect;

    // Constructors
    Camera() = default;
    Camera(float *position, float *axis, float angle, float near, float far,
        float fov, float aspect);

    // Accessors
    Vector3f getPosition();
    Vector3f getAxis();
    float getAngle();
    float getNear();
    float getFar();
    float getFov();
    float getAspect();
};

// Screen resolution and camera object
static const int default_xres = 1000;
static const int default_yres = 1000;

// Per-pixel or per-vertex shading?
static const float default_shader_mode = 1.0;
// Scale of the scene
static const float default_scene_scale = 1.0;
// Should we rebuild the scene?
static const bool default_rebuild_scene = false;
// Is wireframe mode on?
static const bool default_wireframe_mode = false;
// Are we drawing the normals in wireframe mode?
static const bool default_normal_mode = false;
// Are we drawing the IO test?
static const bool default_io_mode = true;
// Are we drawing the intersection of the camera look vector with a primitive?
static const bool default_intersect_mode = true;

// Does the arcball rotate the lights as well as the objects?
static const bool default_arcball_scene = true;
// Is the left mouse button down?
static const bool default_mouse_down = false;

// Arcball rotation matrix for the current update
static const Matrix4f default_arcball_rotate_mat = Matrix4f::Identity();
// Overall arcball object rotation matrix
static const Matrix4f default_arcball_object_mat = Matrix4f::Identity();
// Overall arcball light rotation matrix
static const Matrix4f default_arcball_light_mat = Matrix4f::Identity();

class UI {
    public:
        static UI *singleton;
        static UI *getSingleton();
        static UI *getSingleton(int xres, int yres);

        static void handleMouseButton(int button, int state, int x, int y);
        static void handleMouseMotion(int x, int y);
        static void handleKeyPress(unsigned char key, int x, int y);

        int xres;
        int yres;
        Camera camera;

        float shader_mode;
        float scene_scale;
        bool rebuild_scene;
        bool raytrace_scene;
        bool wireframe_mode;
        bool normal_mode;
        bool io_mode;
        bool intersect_mode;

        Matrix4f arcball_object_mat;
        Matrix4f arcball_light_mat;

        UI(int xres, int yres);

        void reshape(int xres, int yres);

    private:
        int mouse_x;
        int mouse_y;

        bool arcball_scene;
        bool mouse_down;

        Matrix4f arcball_rotate_mat;

        void createCamera();
        Vector3f getArcballVector(int x, int y);
};

#endif