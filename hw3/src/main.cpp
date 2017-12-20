/* The following 2 headers contain all the main functions, data structures, and
 * variables that allow for OpenGL development.
 */
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

/* '_USE_MATH_DEFINES' line allows you to use the syntax 'M_PI'
 * to represent pi to double precision in C++. 
 * Besides the use of 'M_PI', the trigometric functions also show up a lot in
 * graphics computations.
 */
#include <math.h>
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>
#include "parsing.h"
#include "structs.h"
#include <cassert>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////

/* MAJOR FUNCTIONS IN ORDER */

void init(void);
void reshape(int width, int height);
void display(void);

void init_lights();
void set_lights();
void draw_objects();

PointNDC ScreenToNDC(Point screen);
Eigen::Matrix4d Compute_Rotation_Matrix(Point p_current, Point p_start);
Quaternion Compute_Rotation_Quaternion(Point p_current, Point p_start);
Eigen::Matrix4d Get_Current_Rotation();
Quaternion Get_Current_Rotation_Quat();

void mouse_pressed(int button, int state, int x, int y);
void mouse_moved(int x, int y);
void key_pressed(unsigned char key, int x, int y);

///////////////////////////////////////////////////////////////////////////////////////////////////

/* STRUCTS DEFINED IN SEPARATE HEADER FILE */

///////////////////////////////////////////////////////////////////////////////////////////////////

/* CAMERA STORED IN SEPARATE STRUCTURE BECAUSE OF PARSING METHOD */

///////////////////////////////////////////////////////////////////////////////////////////////////

/* GLOBAL VARIABLES */
 
vector<Light> lights; // vector of lights
vector<Object> objects; //  vector of objects
vector<Object> objectsCopies; // vector of all copies
Eigen::Matrix4d last_rotation, current_rotation;
Quaternion last_rotation_quat, current_rotation_quat;
Point p_start, p_current; // points
int xres, yres;
Camera camera;

int mouse_x, mouse_y;
float mouse_scale_x, mouse_scale_y;

const float step_size = 0.2;
const float x_view_step = 90.0, y_view_step = 90.0;
float x_view_angle = 0, y_view_angle = 0;

bool is_pressed = false;
bool wireframe_mode = false;

bool quatMode = true; // Toggle for using quaternions or axis-angle rotation matrices

///////////////////////////////////////////////////////////////////////////////////////////////////

void init(void)
{
    /* Initializing OpenGL features */
    glShadeModel(GL_SMOOTH);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    
    /* Edit projection matrix */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(camera.left, camera.right,
              camera.bottom, camera.top,
              camera.near, camera.far);
    
    /* Edit modelview matrix for transformations */
    glMatrixMode(GL_MODELVIEW);
    
    init_lights();

    /* Setting rotations to identity */
    if (quatMode == false)
    {
        last_rotation.setIdentity(4, 4);
        current_rotation.setIdentity(4, 4);
    }
    else
    {
        last_rotation_quat.qs = 1;
        last_rotation_quat.qx = 0;
        last_rotation_quat.qy = 0;
        last_rotation_quat.qz = 0;

        current_rotation_quat.qs = 1;
        current_rotation_quat.qx = 0;
        current_rotation_quat.qy = 0;
        current_rotation_quat.qz = 0;
    }
}

/*The 'reshape' function is supposed to tell your program how to react
 * whenever the program window is resized. */
void reshape(int width, int height)
{
    height = (height == 0) ? 1 : height;
    width = (width == 0) ? 1 : width;
    
    /* The 'glViewport' function tells OpenGL to determine how to convert from
     * NDC to screen coordinates given the dimensions of the window. */
    glViewport(0, 0, width, height);
    
    /* The following two lines are specific to updating our mouse interface
     * parameters. Details will be given in the 'mouse_moved' function */
    mouse_scale_x = (float) (camera.right - camera.left) / (float) width;
    mouse_scale_y = (float) (camera.top - camera.bottom) / (float) height;
    
    glutPostRedisplay();
}

/* The 'display' function is supposed to handle all the processing of points
 * in world and camera space. */
void display(void)
{
    /* Resetting the "color buffer" is equivalent to clearing the program
     * window so that it only displays a black background. Resetting the 
     * depth buffer sets all the values in the depth buffer back to a very high number. 
     */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    /* To do so, our first step is to "reset" the Modelview Matrix by setting it
     * to the identity matrix: */
    glLoadIdentity();

    /* Our next step is to specify the inverse rotation of the camera by its
     * orientation angle about its orientation axis: */
    glRotatef(-camera.orientation_angle * 180.0 / M_PI,
              camera.orientation_axis[0], 
              camera.orientation_axis[1], 
              camera.orientation_axis[2]);

    /* We then specify the inverse translation of the camera by its position */
    glTranslatef(-camera.position[0], -camera.position[1],  -camera.position[2]);

    /* Get the current rotation matrix and rotate the scene */
    Eigen::Matrix4d rotation;
    Quaternion q = Get_Current_Rotation_Quat();
    if (quatMode == false)
        rotation = Get_Current_Rotation();
    else
        rotation = q.Get_Rotation_Matrix_Quat();

    glMultMatrixd(rotation.data());

    set_lights();
    draw_objects();

    glutSwapBuffers();
}

/* This function has OpenGL enable its built-in lights to represent our point
 * lights. */
void init_lights()
{
    glEnable(GL_LIGHTING);
    
    int num_lights = lights.size();
    
    for(int i = 0; i < num_lights; ++i)
    {
        int light_id = GL_LIGHT0 + i;
        
        glEnable(light_id);
        
        glLightfv(light_id, GL_AMBIENT, lights[i].color);
        glLightfv(light_id, GL_DIFFUSE, lights[i].color);
        glLightfv(light_id, GL_SPECULAR, lights[i].color);
        
        glLightf(light_id, GL_QUADRATIC_ATTENUATION, lights[i].k);
    }
}

/* While the 'init_lights' function enables and sets the colors of the lights,
 * the 'set_lights' function is supposed to position the lights. */
void set_lights()
{
    int num_lights = lights.size();
    
    for(int i = 0; i < num_lights; ++i)
    {
        int light_id = GL_LIGHT0 + i;
        
        glLightfv(light_id, GL_POSITION, lights[i].position);
    }
}

/* This function has OpenGL render our objects to the display screen. */
void draw_objects()
{
    int num_objects = objectsCopies.size();
    Transform t;
    
    for(int i = 0; i < num_objects; ++i)
    {
        glPushMatrix();
        {
            int num_transform_sets = objectsCopies[i].transformations.size();
            
            /* The loop tells OpenGL to modify our modelview matrix
             * in the REVERSE order of how the transformations are specified
             * because OpenGL edits our modelview using post-multiplication */
            for (int j = num_transform_sets - 1; j >= 0; --j)
            {
                t = objectsCopies[i].transformations[j];
                switch (t.label)
                {
                    case "rot":
                        glRotatef(t.rotation_angle * 180.0 / M_PI,
                          t.transformation[0],
                          t.transformation[1],
                          t.transformation[2]);
                        break;
                    case "scal":
                        glScalef(t.transformation[0],
                          t.transformation[1],
                          t.transformation[2]);
                        break;
                    case "trans":
                        glTranslatef(t.transformation[0],
                          t.transformation[1],
                          t.transformation[2]);
                        break;
                }
            }
            
            /* The 'glMaterialfv' and 'glMaterialf' functions tell OpenGL
             * the material properties of the surface we want to render. */
            glMaterialfv(GL_FRONT, GL_AMBIENT, 
                objectsCopies[i].ambient);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, 
                objectsCopies[i].diffuse);
            glMaterialfv(GL_FRONT, GL_SPECULAR, 
                objectsCopies[i].specular);
            glMaterialf(GL_FRONT, GL_SHININESS, 
                objectsCopies[i].shininess);
            
            /* How we tell OpenGL to render geometry for us. 
             * Use points to buffer arrays */
            glVertexPointer(3, GL_FLOAT, 0, &objectsCopies[i].vertex_buffer[0]);
            glNormalPointer(GL_FLOAT, 0, &objectsCopies[i].normal_buffer[0]);
            
            int buffer_size = objectsCopies[i].vertex_buffer.size();
            
            if(!wireframe_mode)
                /* Finally, we tell OpenGL to render everything with the
                 * 'glDrawArrays' function. */
                glDrawArrays(GL_TRIANGLES, 0, buffer_size);
            else
                /* If we are in "wireframe mode" (see the 'key_pressed'
                 * function for more information), then we want to render
                 * lines instead of triangle surfaces. */
                for(int j = 0; j < buffer_size; j += 3)
                    glDrawArrays(GL_LINE_LOOP, j, 3);
        }
        glPopMatrix();
    }
}

/* Reverse maps a screen point to NDC */
PointNDC ScreenToNDC(Point screen)
{
    PointNDC ndc;
    ndc.z = 0.0;
    ndc.x = (float(screen.x + 1) - (xres/2))*2/(xres);
    ndc.y = (float(screen.y + 1) - (yres/2))*2/(-yres);
    return ndc;
}

/* Compute rotation matrix from two points
 * for axis angle rotations */
Eigen::Matrix4d Compute_Rotation_Matrix(Point p_current, Point p_start)
{
    // declare vars
    PointNDC p1, p2;
    Eigen::Vector3d p1vec, p2vec, u;
    float theta, temp;
    float x, y, z, r;

    // convert from Screen x and y to NDC x and y
    p1 = ScreenToNDC(p_start);
    p2 = ScreenToNDC(p_current);

    // calculate z, value is already set to 0 in case x2 + y2 > 1
    if (p1.x*p1.x + p1.y*p1.y <= 1.0)
        p1.z = +sqrt(1 - p1.x*p1.x - p1.y*p1.y);
    if (p2.x*p2.x + p2.y*p2.y <= 1.0)
        p2.z = +sqrt(1 - p2.x*p2.x - p2.y*p2.y);

    // make vector of points
    p1vec << p1.x, p1.y, p1.z;
    p2vec << p2.x, p2.y, p2.z;

    // get theta and u from two points
    temp = min(1.0, p1vec.dot(p2vec)/(p1vec.norm()*p2vec.norm()));
    theta = acos(temp);
    u = p1vec.cross(p2vec);

    // Normalize
    temp = sqrt(u(0)*u(0) + u(1)*u(1) + u(2)*u(2));
    x = u(0)/temp;
    y = u(1)/temp;
    z = u(2)/temp;
    r = theta;

    // Make axis angle rotation matrix and return it
    Eigen::Matrix4d m;
    float cosr = cos(r);
    float sinr = sin(r);
    m << x*x+(1-x*x)*cosr, x*y*(1-cosr)-z*sinr, 
    x*z*(1-cosr)+y*sinr, 0,
    y*x*(1-cosr)+z*sinr, y*y+(1-y*y)*cosr, 
    y*z*(1-cosr)-x*sinr, 0,
    z*x*(1-cosr)-y*sinr, z*y*(1-cosr)+x*sinr, 
    z*z+(1-z*z)*cosr, 0,
    0, 0, 0, 1;
    return m;
}

Quaternion Compute_Rotation_Quaternion(Point p_current, Point p_start)
{
    PointNDC p1, p2;
    Quaternion q1, q2, q;
    Eigen::Vector3d p1vec, p2vec, u;
    float theta, temp;
    float x, y, z, r;

    // convert from Screen x and y to NDC x and y
    p1 = ScreenToNDC(p_start);
    p2 = ScreenToNDC(p_current);

    // calculate z, value is already set to 0 in case x2 + y2 > 1
    if (p1.x*p1.x + p1.y*p1.y <= 1.0)
        p1.z = +sqrt(1 - p1.x*p1.x - p1.y*p1.y);
    if (p2.x*p2.x + p2.y*p2.y <= 1.0)
        p2.z = +sqrt(1 - p2.x*p2.x - p2.y*p2.y);

    // represent points as quaternions
    q1.qs = 0;
    q1.qx = p1.x;
    q1.qy = p1.y;
    q1.qz = p1.z;

    q2.qs = 0;
    q2.qx = p2.x;
    q2.qy = p2.y;
    q2.qz = p2.z;

    p1vec << p1.x, p1.y, p1.z;
    p2vec << p2.x, p2.y, p2.z;

    temp = min(1.0, p1vec.dot(p2vec)/(p1vec.norm()*p2vec.norm()));
    theta = acos(temp);
    u = p1vec.cross(p2vec);

    // Normalize
    temp = sqrt(u(0)*u(0) + u(1)*u(1) + u(2)*u(2));
    x = u(0)/temp;
    y = u(1)/temp;
    z = u(2)/temp;
    r = theta;

    // create unit rotation quaternion
    q.qs = cos(r/2);
    q.qx = x*sin(r/2);
    q.qy = y*sin(r/2);
    q.qz = z*sin(r/2);

    return q;
}

/* Get rotation matrix by multiplying the current rotation matrix
 * with the last one */
Eigen::Matrix4d Get_Current_Rotation()
{
    return current_rotation*last_rotation;
}

/* Using current and last rotation quaternions, get rotation
 * quaternion by multiplying them together */
Quaternion Get_Current_Rotation_Quat()
{
    // declare variables
    Quaternion q;
    Eigen::Vector3d va, vb, v;
    float sa, sb, s;

    // make new quaternion out of the two quaternions
    sa = current_rotation_quat.qs;
    sb = last_rotation_quat.qs;

    va << current_rotation_quat.qx,
        current_rotation_quat.qy,
        current_rotation_quat.qz;
    vb << last_rotation_quat.qx,
        last_rotation_quat.qy,
        last_rotation_quat.qz;

    s = sa*sb - va.dot(vb);
    v = sa*vb + sb*va + va.cross(vb);

    q.qs = s;
    q.qx = v(0);
    q.qy = v(1);
    q.qz = v(2);

    return q;
}

/* This function is meant to respond to mouse clicks and releases. */
void mouse_pressed(int button, int state, int x, int y)
{
    /* If the left-mouse button was clicked down, then...
     */
    if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        p_start.x = x;
        p_start.y = y;
        /* Store the mouse position in our global variables.
         */
        mouse_x = x;
        mouse_y = y;
        
        /* Since the mouse is being pressed down, we set our 'is_pressed"
         * boolean indicator to true.
         */
        is_pressed = true;
    }
    /* If the left-mouse button was released up, then...
     */
    else if(button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    {
        /* Mouse is no longer being pressed, so set our indicator to false.
         */
        is_pressed = false;

        /* Since mouse is unpressed, set the last quaternion to current
         * one and set the current one to identity */
        if (quatMode == false)
        {
            last_rotation = Get_Current_Rotation();
            current_rotation.setIdentity(4, 4);
        }
        else
        {
            last_rotation_quat = Get_Current_Rotation_Quat();
            current_rotation_quat.qs = 1;
            current_rotation_quat.qx = 0;
            current_rotation_quat.qy = 0;
            current_rotation_quat.qz = 0;
        }
    }
}

/* This function is meant to respond to when the mouse is being moved.*/
void mouse_moved(int x, int y)
{
    /* If the left-mouse button is being clicked down...
     */
    if(is_pressed)
    {
        p_current.x = x;
        p_current.y = y;
        if (quatMode == false)
            current_rotation = Compute_Rotation_Matrix(p_current, p_start);
        else
            current_rotation_quat = Compute_Rotation_Quaternion(p_current, p_start);
        
        mouse_x = x;
        mouse_y = y;
        
        glutPostRedisplay();
    }
}

/* 'deg2rad' function:
 * 
 * Converts given angle in degrees to radians.
 */
float deg2rad(float angle)
{
    return angle * M_PI / 180.0;
}

/* This function is meant to respond to key pressed on the keyboard. */
void key_pressed(unsigned char key, int x, int y)
{
    /* If 'q' is pressed, quit the program.
     */
    if(key == 'q')
    {
        exit(0);
    }
    /* If 't' is pressed, toggle our 'wireframe_mode' boolean to make OpenGL
     * render our cubes as surfaces of wireframes.
     */
    else if(key == 't')
    {
        wireframe_mode = !wireframe_mode;
        glutPostRedisplay();
    }
    else
    {        
        float x_view_rad = deg2rad(x_view_angle);
        
        /* 'w' for step forward
         */
        if(key == 'w')
        {
            camera.position[0] += step_size * sin(x_view_rad);
            camera.position[2] -= step_size * cos(x_view_rad);
            glutPostRedisplay();
        }
        /* 'a' for step left
         */
        else if(key == 'a')
        {
            camera.position[0] -= step_size * cos(x_view_rad);
            camera.position[2] -= step_size * sin(x_view_rad);
            glutPostRedisplay();
        }
        /* 's' for step backward
         */
        else if(key == 's')
        {
            camera.position[0] -= step_size * sin(x_view_rad);
            camera.position[2] += step_size * cos(x_view_rad);
            glutPostRedisplay();
        }
        /* 'd' for step right
         */
        else if(key == 'd')
        {
            camera.position[0] += step_size * cos(x_view_rad);
            camera.position[2] += step_size * sin(x_view_rad);
            glutPostRedisplay();
        }
    }
}

int main(int argc, char* argv[])
{
    assert(argc == 4);

    /*** DECLARE LOCAL VARIABLES ***/

    string line;
    bool checker = false;
    bool checker2;

    fstream myfile(argv[1]);
    istringstream(argv[2]) >> xres;
    istringstream(argv[3]) >> yres;

    /*** PARSE FILE AND PUT INTO DATA STRUCTURES ***/

    if (myfile.is_open())
    {
        while (!myfile.eof())
        {
            getline(myfile, line);
            parse_all(line, objects, objectsCopies, checker, checker2,
                camera, lights);
        }
    }
    myfile.close();

    /* 'glutInit' intializes the GLUT (Graphics Library Utility Toolkit) library. */
    glutInit(&argc, argv);

    /* The following line of code tells OpenGL that we need a double buffer,
     * a RGB pixel buffer, and a depth buffer. */
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    /* The following line tells OpenGL to create a program window of size
     * 'xres' by 'yres'. */
    glutInitWindowSize(xres, yres);

    /* The following line tells OpenGL to set the program window in the top-left
     * corner of the computer screen (0, 0). */
    glutInitWindowPosition(0, 0);

    /* The following line tells OpenGL to name the program window */
    glutCreateWindow("Render");
    
    /* Call our 'init' function */
    init();

    /* Specify to OpenGL our display function.
     */
    glutDisplayFunc(display);
    /* Specify to OpenGL our reshape function.
     */
    glutReshapeFunc(reshape);
    /* Specify to OpenGL our function for handling mouse presses.
     */
    glutMouseFunc(mouse_pressed);
    /* Specify to OpenGL our function for handling mouse movement.
     */
    glutMotionFunc(mouse_moved);
    /* Specify to OpenGL our function for handling key presses.
     */
    glutKeyboardFunc(key_pressed);
    /* The following line tells OpenGL to start the "event processing loop". This
     * is an infinite loop where OpenGL will continuously use our display, reshape,
     * mouse, and keyboard functions to essentially run our program.
     */
    glutMainLoop();
}
