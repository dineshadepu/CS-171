#include "UI.hpp"

UI *UI::singleton = NULL;

/* Constructs the program's camera. */
Camera::Camera(float *position, float *axis, float angle, float near, float far,
    float fov, float aspect)
{
    for (int i = 0; i < 3; i++) {
        this->position = Vec3f(position);
        this->axis = Vec3f(axis);
    }
    this->angle = angle;
    this->near = near;
    this->far = far;
    this->fov = fov;
    this->aspect = aspect;
}

Vector3f Camera::getPosition() {
    return Vector3f(this->position.x, this->position.y, this->position.z);
}

Vector3f Camera::getAxis() {
    return Vector3f(this->axis.x, this->axis.y, this->axis.z);
}

float Camera::getAngle() {
    return this->angle;
}

float Camera::getNear() {
    return this->near;
}

float Camera::getFar() {
    return this->far;
}

float Camera::getFov() {
    return this->fov;
}

float Camera::getAspect() {
    return this->aspect;
}

/* Initializes the UI with the given resolution, and creates the camera. */
UI::UI(int xres, int yres) : 
    xres(xres),
    yres(yres),
    shader_mode(default_shader_mode),
    scene_scale(default_scene_scale),
    rebuild_scene(default_rebuild_scene),
    wireframe_mode(default_wireframe_mode),
    normal_mode(default_normal_mode),
    io_mode(default_io_mode),
    intersect_mode(default_intersect_mode),
    arcball_object_mat(default_arcball_object_mat),
    arcball_light_mat(default_arcball_light_mat),

    arcball_scene(default_arcball_scene),
    mouse_down(default_mouse_down),
    arcball_rotate_mat(default_arcball_rotate_mat)
{
    this->createCamera();
}

/* Sets up the scene's camera. */
void UI::createCamera() {
    float position[3] = {0.0, 0.0, 10.0};
    float axis[3] = {0.0, 0.0, 1.0};
    float angle = 0.0;
    float near = 0.1;
    float far = 500.0;
    float fov = 60.0;
    float aspect = (float) xres / yres;

    // Create the Camera struct
    this->camera = Camera(position, axis, angle, near, far, fov, aspect);
}

/* Returns/sets up the singleton instance of the class. */
UI *UI::getSingleton() {
    if (!UI::singleton) {
        UI::singleton = new UI(default_xres, default_yres);
    }

    return UI::singleton;
}

/* Returns/sets up the singleton instance of the class. */
UI *UI::getSingleton(int xres, int yres) {
    if (!UI::singleton) {
        UI::singleton = new UI(xres, yres);
    }
    else {
        UI::singleton->reshape(xres, yres);
    }

    return UI::singleton;
}

/* Handles mouse click events. */
void UI::handleMouseButton(int button, int state, int x, int y) {
    UI *ui = UI::getSingleton();

    // If the action pertains to the left button...
    if (button == GLUT_LEFT_BUTTON) {
        // If the button is being held down, update the start of the arcball
        // rotation, and store that the button is currently down
        ui->mouse_down = state == GLUT_DOWN;
        if (state == GLUT_DOWN) {
            ui->mouse_x = x;
            ui->mouse_y = y;
            ui->mouse_down = true;
        }
        // Otherwise, store that the button is up
        else
            ui->mouse_down = false;
    }
}

/* Handles mouse motion events. */
void UI::handleMouseMotion(int x, int y) {
    UI *ui = UI::getSingleton();

    // If the left button is being clicked, and the mouse has moved, and the
    // mouse is in the window, then update the arcball UI
    if (ui->mouse_down && (x != ui->mouse_x || y != ui->mouse_y) &&
        (x >= 0 && x < ui->xres && y >= 0 && y < ui->yres))
    {
        // Set up some matrices we need
        Matrix4f camera_to_ndc = Matrix4f();
        Matrix4f world_to_camera = Matrix4f();
        glGetFloatv(GL_PROJECTION_MATRIX, camera_to_ndc.data());
        glGetFloatv(GL_MODELVIEW_MATRIX, world_to_camera.data());
        Matrix3f ndc_to_world = camera_to_ndc.topLeftCorner(3, 3).inverse();

        // Get the two arcball vectors by transforming from NDC to camera
        // coordinates, ignoring translation components
        Vector3f va =
            (ndc_to_world * ui->getArcballVector(ui->mouse_x, ui->mouse_y)).normalized();
        Vector3f vb = (ndc_to_world * ui->getArcballVector(x, y)).normalized();
        
        // Compute the angle between them and the axis to rotate around
        // (this time rotated into world space, where the matrix is applied)
        Vector3f arcball_axis =
            (world_to_camera.topLeftCorner(3, 3).transpose() * va.cross(vb)).normalized();
        float arcball_angle = acos(fmax(fmin(va.dot(vb), 1.0), -1.0));

        // Update current arcball rotation and overall object rotation matrices
        makeRotateMat(ui->arcball_rotate_mat.data(), arcball_axis[0],
            arcball_axis[1], arcball_axis[2], arcball_angle);
        ui->arcball_object_mat =
            ui->arcball_rotate_mat * ui->arcball_object_mat;
        
        // If the arcball should rotate the entire scene, update the light
        // rotation matrix too 
        if (ui->arcball_scene) {
            ui->arcball_light_mat =
                ui->arcball_rotate_mat * ui->arcball_light_mat;
        }
        
        // Update the arcball start position
        ui->mouse_x = x;
        ui->mouse_y = y;
        
        // Update the image
        glutPostRedisplay();
    }
}

/* Handles keyboard events. */
void UI::handleKeyPress(unsigned char key, int x, int y) {
    UI *ui = UI::getSingleton();

    // Q quits the program
    if (key == 'q' || key == 27) {
        CommandLine::run();
    }
    // B toggles whether the arcball rotates the lights along with the objects
    else if (key == 'b') {
        ui->arcball_scene = !ui->arcball_scene;
    }
    // M toggles between per-pixel and per-vertex shading
    else if (key == 'm') {
        ui->shader_mode = !((int) ui->shader_mode);
    }
    // W blows up the scene
    else if (key == 'w') {
        ui->scene_scale *= 1.1;
    }
    // S shrinks the scene
    else if (key == 's') {
        ui->scene_scale /= 1.1;
    }
    // R forces a refresh of the scene's entities
    else if (key == 'r') {
        ui->rebuild_scene = true;
    }
    // T toggles wireframe mode
    else if (key == 't') {
        ui->wireframe_mode = !ui->wireframe_mode;
    }
    // N toggles normal mode
    else if (key == 'n') {
        ui->normal_mode = !ui->normal_mode;
    }
    // O toggles IO mode
    else if (key == 'o') {
        ui->io_mode = !ui->io_mode;
    }
    // I toggles intersect mode
    else if (key == 'i') {
        ui->intersect_mode = !ui->intersect_mode;
    }
    // P raytraces the scene
    else if (key == 'p') {
        ui->raytrace_scene = true;
    }
    glutPostRedisplay();
}

/* Gets an arcball vector from a position in the program's window. */
Vector3f UI::getArcballVector(int x, int y) {
    // Convert from screen space to NDC
    Vector3f p(2.0 * x / this->xres - 1.0, -(2.0 * y / this->yres - 1.0), 0.0);
    // Compute the appropriate z-coordinate, and make sure it's normalized
    float squared_length = p.squaredNorm();
    if (squared_length < 1.0)
        p[2] = -sqrtf(1.0 - squared_length);
    else
        p /= sqrtf(squared_length);

    return p;
}

/* Handles window resize events. */
void UI::reshape(int xres, int yres) {
    // Update the internal resolution
    this->xres = xres;
    this->yres = yres;
}
