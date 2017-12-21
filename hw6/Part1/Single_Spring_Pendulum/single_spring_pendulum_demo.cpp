/* CS/CNS 171
 * Fall 2015
 * Written by Kevin (Kevli) Li (Class of 2016)
 *
 * This program is meant to simulate a 2D spring pendulum using the discrete
 * Lagrangian. The main part of the program that you need to worry about is
 * marked with a TODO, though you may need to skim through the file to find
 * the variables that you need to know in order to do the TODO. You should not
 * need to edit any other part of this file other than the TODO.
 *
 * The file is sparsely commented since you are not required to understand
 * everything in this program. If you have any questions, feel free to send
 * me an email.
 */

/* Usual includes and namespaces */
#include <GL/glew.h>
#include <GL/glut.h>

#include <cfloat>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

/* Represents a point for our path tracing of the pendulum */
struct Point
{
    float x;
    float y;
};

/* Represents a single color; used for the color pointer later */
struct Color
{
    float r;
    float g;
    float b;
};

/* Represents a spring pendulum for our simulation */
struct Spring_Pendulum
{
    float m; // mass
    
    float x; // x position
    float y; // y position
    
    float px; // x-momentum
    float py; // y-momentum
    
    float k; // spring constant
    float rl; // rest length of spring
};

/* Camera and perspective parameters */
const float cam_position[] = {0, 0, 2};
const float near_param = 1, far_param = 4,
            left_param = -8, right_param = 8,
            top_param = 6, bottom_param = -10;

/* Arbitrary light */
const float light_color[3] = {1, 1, 1};
const float light_position[3] = {0, 0, -2};

/* Arbitrary material properties for our pendulum */
const float ambient_reflect[3] = {0.3, 0.2, 0.4};
const float diffuse_reflect[3] = {0.7, 0.2, 0.8};
const float specular_reflect[3] = {1, 1, 1};
const float shininess = 0.1;

/* Time step and gravity constant for simulation */
const float dt = 0.01; // our time step
const float g = -9.8; // acceleration of gravity

/* Keeps track of total time of simulation */
float t = 0;

/* Max number of points to keep in our traced path */
const int max_num_points = 1000;

/* Kinetic energy, potential energy, and min/max total energy variables */
float ke, pe;
float min_total = FLT_MAX, max_total = -FLT_MAX;

/* Initial value of Lagrangian for color value calculations */
float lagrangian_0;

/* Our pendulum */
Spring_Pendulum m1;

/* For drawing the path of the pendulum over time */
vector<Point> path1;
vector<Color> path_colors1;

/* Some key prototype functions */
void trace_path();
void draw_spring_pendulum();
void draw_text();

void init(void)
{
    glEnableClientState(GL_VERTEX_ARRAY);

    // Similar to how we have vertex arrays and normal arrays,
    // we also have color arrays for setting the colors of each
    // vertex
    glEnableClientState(GL_COLOR_ARRAY);

    // Set the size of our points for the paths we trace
    glPointSize(2);

    // Set lighting
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_color);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_color);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    // Set material properties of our pendulum
    glMaterialfv(GL_FRONT, GL_AMBIENT, ambient_reflect);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_reflect);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular_reflect);
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    // We use glOrtho instead of glFrustum, since we're only doing 2D displays
    glOrtho(left_param, right_param,
            bottom_param, top_param,
            near_param, far_param);
    
    glMatrixMode(GL_MODELVIEW);
}

void reshape(int width, int height)
{
    height = (height == 0) ? 1 : height;
    width = (width == 0) ? 1 : width;
    
    glViewport(0, 0, width, height);
    
    glutPostRedisplay();
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glTranslatef(-cam_position[0], -cam_position[1], -cam_position[2]);
    
    // Draw the path that the pendulum has taken (without lighting effects)
    trace_path();
    
    // Draw the pendulum system
    glEnable(GL_LIGHTING);
    draw_spring_pendulum();
    glDisable(GL_LIGHTING);

    // Draw the numbers displayed on the screen
    draw_text();
    
    glutSwapBuffers();
}

/* Computes the continuous Lagrangian, kinetic energy, potential energy, and total
 * values for coloring the traced path and for displaying on the screen */
float compute_lagrangian()
{
    ke = 1.0 / 2.0 * m1.m * ((m1.px / m1.m) * 
        (m1.px / m1.m) + (m1.py / m1.m) * (m1.py / m1.m));
    pe = 1.0 / 2.0 * m1.k * 
    (sqrt(m1.x * m1.x + m1.y * m1.y) - m1.rl) *
    (sqrt(m1.x * m1.x + m1.y * m1.y) - m1.rl)
         - m1.m * g * m1.y;

    float total = ke + pe;
    min_total = (total < min_total) ? total : min_total;
    max_total = (total > max_total) ? total : max_total;
    
    return ke - pe;
}

/* Updates the path of the pendulum when it moves to a new point */
void update_path()
{    
    // Delete the oldest point in the path when the paths get too big to
    // prevent slowdown
    if(path1.size() == max_num_points)
    {
        path1.erase(path1.begin());
        path_colors1.erase(path_colors1.begin());
    }
    
    // Make new point
    Point point1;
    point1.x = m1.x;
    point1.y = m1.y;
    
    // Compute a "normalized" Lagrangian value by dividing the current Lagrangian
    // by the initial value of the Lagrangian; this isn't exactly a normalization,
    // but it suffices for generating nice colors
    float lagrangian_norm = abs(compute_lagrangian() / lagrangian_0);
    lagrangian_norm = (lagrangian_norm > 1.0) ? 1.0 : lagrangian_norm;
    
    Color color1;
    color1.r = lagrangian_norm;
    color1.g = lagrangian_norm;
    color1.b = 1.0 - lagrangian_norm;
    
    // Add new point and corresponding color to our lists
    path1.push_back(point1);
    path_colors1.push_back(color1);
}

void update_pendulum()
{
    /******************************* TODO *******************************/

    /* Your task is to write some lines of code to update:
     * 
     *     m1.x
     *     m1.y
     *     m1.px
     *     m1.py
     *
     * (not necessarily in that order)
     * To start, you should write the continuous Lagrangian for the spring
     * pendulum system. Then, you should write the discrete analog. From
     * there, use the discrete Euler-Lagrangian equations to solve for
     * the correct update rules.
     *
     * The variables that you'll need (in addition to the above-listed
     * ones) are:
     *
     *     dt
     *     m1.m
     *     m1.k
     *     m1.rl
     *     g
     * 
     */

    /* The continuous Lagrangian for the spring pendulum system,
     * with a pendulum of mass m, stiffness k, rest length l:
     * L = T - U = 0.5 * m * (v_x^2 + v_y^2) 
                   - 0.5 * k * (sqrt(x^2 + y^2) - l)^2
                   + m * g * y
     * the discrete analog: FINISH WRITING THIS OUT
     * L_d (q_k, q_k+1, dt) = dt * (0.5 * m)...
     * Using discrete analog and DEL, we can compute px, py, x, and y
     */

    // Calculate x_k+1 and y_k+1

    float tempx = m1.x + dt / m1.m * 
        (m1.px - (dt * m1.k * m1.x * (sqrt(m1.x * m1.x + m1.y * m1.y) - m1.rl))
        / sqrt(m1.x * m1.x + m1.y * m1.y));

    float tempy = m1.y + dt / m1.m * 
        (m1.py - (dt * m1.k * m1.y * (sqrt(m1.x * m1.x + m1.y * m1.y) - m1.rl))
        / sqrt(m1.x * m1.x + m1.y * m1.y) + g * m1.m * dt);

    // Using x_k and x_k+1, get p_x_k+1. Same for y

    m1.px = m1.m * (tempx - m1.x) / dt;
    m1.py = m1.m * (tempy - m1.y) / dt;

    // Update x_k and y_k with x_k+1 and y_k+1
    m1.x = tempx;
    m1.y = tempy;

    /****************************** END TODO ****************************/

    t += dt;
}

/* Updates our pendulum's position and its traced path */
void update()
{
    update_path();
    update_pendulum();
    
    glutPostRedisplay();
}

/* Draws the path of the pendulum using OpenGL's vertex and color pointers */
void trace_path()
{
    glVertexPointer(2, GL_FLOAT, 0, &path1[0]);
    glColorPointer(3, GL_FLOAT, 0, &path_colors1[0]);
    glDrawArrays(GL_POINTS, 0, path1.size());
}

/* Draws our spring pendulum system */
void draw_spring_pendulum()
{
    // Draws line connecting anchor point and pendulum bob
    glBegin(GL_LINES);
    glVertex2f(0, 0);
    glVertex2f(m1.x, m1.y);
    glEnd();

    float pendulum_radius = 0.3;

    // Draw anchor point
    glutSolidSphere(pendulum_radius * 0.5, 20, 20);
    
    // Draw pendulum bob in correct position
    glPushMatrix();
    {
        glTranslatef(m1.x, m1.y, 0);
        glutSolidSphere(pendulum_radius, 20, 20);
    }
    glPopMatrix();
}

/* Converts given input to a string */
template<typename T>
string tostr(const T& t)
{
    ostringstream os;
    os << t;
    return os.str();
}

/* Draws the text on the screen */
void draw_text()
{
    // Set drawing color to green
    glColor3f(0, 1, 0);

    // Initialize strings to draw
    string ke_str = "KE: " + tostr(ke);
    string pe_str = "PE: " + tostr(pe);
    string total_str = "Total: " + tostr(ke + pe);
    string min_str = "Min Total: " + tostr(min_total);
    string max_str = "Max Total: " + tostr(max_total);
    string t_str = "Time: " + tostr(t);

    // The following lines of code show one way to draw text in OpenGL.
    // glRasterPos sets the x and y coordinates of where you want to
    // begin the text (and if 3D, you set z as well). Then, the function,
    // glutBitmapCharacter, draws a single character; the function takes
    // in a font specification and the character you want to draw.
    glRasterPos2f(-7.4,5);
    for(int i = 0; i < ke_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, ke_str[i]);

    glRasterPos2f(-7.4,4);
    for(int i = 0; i < pe_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, pe_str[i]);

    glRasterPos2f(3.64,5);
    for(int i = 0; i < total_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, total_str[i]);

    glRasterPos2f(-7.4,-8.3);
    for(int i = 0; i < min_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, min_str[i]);

    glRasterPos2f(-7.4,-9.3);
    for(int i = 0; i < max_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, max_str[i]);

    glRasterPos2f(3.64,-9.3);
    for(int i = 0; i < t_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, t_str[i]);
}

void key_pressed(unsigned char key, int x, int y)
{
    if(key == 'q')
    {
        exit(0);
    }
}

int main(int argc, char* argv[])
{
    if(argc != 5)
    {
        cerr << "\nERROR: Incorrect number of arguments." << endl;
        exit(1);
    }
    
    // xres and yres are given command line arguments
    int xres = atoi(argv[1]);
    int yres = atoi(argv[2]);
    
    // Set mass to 1 for simplicity
    m1.m = 1;

    // Initial position of pendulum bob is given by command line arguments
    m1.x = atof(argv[3]);
    m1.y = atof(argv[4]);

    // No velocity in the beginning
    m1.px = 0;
    m1.py = 0;

    // Arbitrary spring constant
    m1.k = 20;

    // Set rest length to 1 for simplicity
    m1.rl = 1;
    
    // Set initial Lagrangian value for generating the color effect of the
    // pendulum's path
    lagrangian_0 = compute_lagrangian();
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(xres, yres);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Spring Pendulum Demo");
    
    init();
    glutDisplayFunc(display);
    glutIdleFunc(update);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(key_pressed);
    glutMainLoop();
}
