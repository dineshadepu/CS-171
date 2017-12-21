/* CS/CNS 171
 * Fall 2015
 * Written by Kevin (Kevli) Li (Class of 2016)
 *
 * This program is meant to simulate elasticity in a 2D mesh. The only
 * part of the program that you need to worry about is marked with a
 * TODO. You should not need to edit any other code in this file.
 *
 * The file is sparsely commented since you are not required to understand
 * everything in this program. If you have any questions, feel free to send
 * me an email.
 */

/* Usual includes and namespaces */
#include <GL/glew.h>
#include <GL/glut.h>
#include <Eigen/Dense>

#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

/* Separate header file */
#include "utils.h"

using namespace std;
using namespace Eigen;

/* Represents a vertex in our mesh for the simulation */
struct Vertex
{
    float x, y; // position of vertex
    float restx, resty; // positions of vertex when mesh is at rest
    float vx, vy; // velocity of vertex
    float fx, fy; // force acting on vertex
    float mass; // mass, which we're going to set to 1 for simplicity
};

/* Represents a triangular face in our mesh */
struct Triangle
{
    int p1, p2, p3; // indices of the vertices that make up the face
    float restarea; // the area of the triangle when the mesh is at rest
};

/* Perspective parameters */
GLfloat nearParam = -1.0, farParam = 1.0,
        leftParam = -2, rightParam = 2,
        bottomParam = -2, topParam = 2;

/* Time step and total time variables for the simulation */
const float tstep = 0.2;
float total_time = 0.0;

/* Kinetic energy, potential energy, and min/max total energy variables */
float ke, pe;
float min_total = FLT_MAX, max_total = -FLT_MAX;

/* Screen resolutions */
const int xres = 800, yres = 800;

/* For mouse UI */
const float clickRad = 0.05;
int mouseX, mouseY;
float mouseScaleX = (rightParam - leftParam) / (float) xres,
      mouseScaleY = (topParam - bottomParam) / (float) yres;
int draggedPoint;
bool dragging = false;

/* Our list of vertices */
vector<Vertex> vertices;
/* Our list of faces */
vector<Triangle> triangles;
/* List of "rest matrices"; needed for simulation */
vector<MatrixXd> restmats;

/* Parses in an OBJ file and fills vertices, triangles, and restmats above */
static bool parse_OBJ(std::string file_name)
{
    std::string read_line;
    std::ifstream file_in(file_name.c_str());

    if(file_in.is_open())
    {
        while(getline(file_in, read_line))
        {
            std::vector<std::string> tokens = split(read_line, ' ');
            
            if(tokens.size() == 0)
                continue;

            // Read in vertex, initialize fields, and push into list
            if(tokens[0].compare("v") == 0)
            {
                Vertex v;
                v.x = stof(tokens[1]);
                v.y = stof(tokens[2]);

                v.restx = v.x;
                v.resty = v.y;
    
                v.vx = 0;
                v.vy = 0;
    
                v.fx = 0;
                v.fy = 0;
    
                v.mass = 1.0;
    
                vertices.push_back(v);
            }
            // Read in face, initialize fields and corresponding rest matrix,
            // and push them onto the appropriate lists
            else if(tokens[0].compare("f") == 0)
            {
                Triangle t;
    
                t.p1 = stoi(split(tokens[1], '/').at(0)) - 1;
                t.p2 = stoi(split(tokens[2], '/').at(0)) - 1;
                t.p3 = stoi(split(tokens[3], '/').at(0)) - 1;
    
                // Don't worry about this matrix; it's needed for the simulation
                MatrixXd restmat(2,2);
                restmat << vertices.at(t.p1).restx - vertices.at(t.p3).restx,
                           vertices.at(t.p2).restx - vertices.at(t.p3).restx,
                           vertices.at(t.p1).resty - vertices.at(t.p3).resty,
                           vertices.at(t.p2).resty - vertices.at(t.p3).resty;
                restmats.push_back(restmat);
    
                // Also don't worry about the rest area
                t.restarea = 1.0 / 2.0 * abs(restmat(0,0) * restmat(1,1) -
                                             restmat(1,0) * restmat(0,1));
    
                triangles.push_back(t);
            }
        }
        
        file_in.close();
        return 1;
    }
    else
    {
        std::cerr << "ERROR: Cannot find specified input file: " << file_name << "!\n";
        return 0;
    }
}

void init(void)
{   
    // Tells OpenGL to draw our points with their fastest method
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    
    // Tells OpenGL to draw our lines with their fastest method
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // We use glOrtho instead of glFrustum, since we're only doing 2D displays
    glOrtho(leftParam, rightParam,
            bottomParam, topParam,
            nearParam, farParam);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void reshape(int width, int height)
{
    height = (height == 0) ? 1 : height;
    width = (width == 0) ? 1 : width;

    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    glutPostRedisplay();
}

/* This function updates the forces that act on each vertex due to the elasticity
 * of the mesh. You DO NOT need to understand what is going on in this function
 * (though you are free to explore it for extra credit). Once this function is
 * called, each vertex in vertices will have their fx and fy fields with the
 * correct forces for the next time step. The potential energy variable, pe,
 * is also updated.
 */
void updateforces(int index)
{
    Triangle t = triangles[index];
    
    MatrixXd ds(2, 2);
    ds << vertices[t.p1].x - vertices[t.p3].x, vertices[t.p2].x - vertices[t.p3].x,
          vertices[t.p1].y - vertices[t.p3].y, vertices[t.p2].y - vertices[t.p3].y;
    
    MatrixXd dminverse = (restmats[index]).inverse();
    MatrixXd f = ds * dminverse;
    MatrixXd finvt = (f.inverse()).transpose();
    
    float i1 = (f.transpose() * f).trace();
    float j = f(0,0) * f(1,1) - f(1,0) * f(0,1);
    
    MatrixXd p = f - finvt + 1.0 / 2.0 * log(j * j) * finvt;
    MatrixXd h = -t.restarea * p * dminverse.transpose();
    
    vertices[t.p1].fx += h(0,0);
    vertices[t.p2].fx += h(0,1);
    vertices[t.p3].fx += -h(0,0) - h(0,1);
    
    vertices[t.p1].fy += h(1,0);
    vertices[t.p2].fy += h(1,1);
    vertices[t.p3].fy += -h(1,0) - h(1,1);
    
    pe += t.restarea * ( (i1 - 2) - log(j) + 1.0 / 2.0 * log(j) * log(j) );
}

void drawTriangle(int index)
{
    Triangle t = triangles[index];

    // Set drawing color to cyan
    glColor3f(0, 1, 1);
 
    // In practice, it's always better to use vertex pointers for drawing; but
    // since this is meant to be a quick demo, we're being lazy here and using
    // glVertex instead. glVertex has been deprecated in a lot of new versions
    // of OpenGL, and using it is slower than using vertex pointers.
    glBegin(GL_TRIANGLES);
    {
        glVertex2f(vertices[t.p1].x, vertices[t.p1].y);
        glVertex2f(vertices[t.p2].x, vertices[t.p2].y);
        glVertex2f(vertices[t.p3].x, vertices[t.p3].y);
    }
    glEnd();
    
    // Set drawing color to black
    glColor3f(0, 0, 0);
    
    glBegin(GL_LINE_LOOP);
    {
        glVertex2f(vertices[t.p1].x, vertices[t.p1].y);
        glVertex2f(vertices[t.p2].x, vertices[t.p2].y);
        glVertex2f(vertices[t.p3].x, vertices[t.p3].y);
    }
    glEnd();
}

/* Converts given input to a string */
template<typename T>
string tostr(const T& t)
{
    ostringstream os; os<<t;
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
    string t_str = "Time: " + tostr(total_time);

    // The following lines of code show one way to draw text in OpenGL.
    // glRasterPos sets the x and y coordinates of where you want to
    // begin the text (and if 3D, you set z as well). Then, the function,
    // glutBitmapCharacter, draws a single character; the function takes
    // in a font specification and the character you want to draw.

    glRasterPos2f(-1.93, 1.87);
    for(int i = 0; i < ke_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, ke_str[i]);

    glRasterPos2f(-1.93, 1.73);
    for(int i = 0; i < pe_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, pe_str[i]);

    glRasterPos2f(1.1, 1.87);
    for(int i = 0; i < total_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, total_str[i]);

    glRasterPos2f(-1.93, -1.93);
    for(int i = 0; i < min_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, min_str[i]);

    glRasterPos2f(-1.93, -1.79);
    for(int i = 0; i < max_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, max_str[i]);

    glRasterPos2f(1.32, -1.93);
    for(int i = 0; i < t_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, t_str[i]);
}

void display(void)
{ 
    glClear(GL_COLOR_BUFFER_BIT);
    
    // Reset our kinetic and potential energies, since we're going to recompute
    // them for the new iteration.
    ke = 0.0;
    pe = 0.0;
    
    // Reset our forces to 0, since we're going to recompute them for the new
    // iteration.
    for(int i = 0; i < vertices.size(); i++)
    {
        vertices[i].fx = 0;
        vertices[i].fy = 0;
    }
    
    // Draw the last iterations's mesh
    for(int i = 0; i < triangles.size(); i++)
        drawTriangle(i);
    
    // Update the forces for the next iteration
    for(int i = 0; i < triangles.size(); i++)
        updateforces(i);
    
    // Update the positions and velocities of each vertex in our mesh
    for(int i = 0; i < vertices.size(); i++)
    {
        // Ignore the vertex that we're dragging; that should stay where
        // our mouse is
        if(dragging && draggedPoint == i)
            continue;
        
        /******************************* TODO *******************************/

        /* Your task is to write some lines of code to correctly update:
         *
         *     vertices[i].x
         *     vertices[i].y
         *     vertices[i].vx
         *     vertices[i].vy
         *
         * (not necessarily in that order)
         * You are to update these variables knowing that this system follows
         * the standard Lagrangian:
         *
         *     L = 1/2 * m * ( v_x^2 + v_y^2 ) - U( x, y )
         *
         * Remember that:
         *
         *     d/dx U( x, y ) = -f_x = -vertices[i].fx
         *     d/dy U( x, y ) = -f_y = -vertices[i].fy
         *
         * where vertices[i].fx and vertices[i].fy are already computed for you via
         * the updateforces() function.
         *
         * Besides the above mentioned variables, you'll find the variable, tstep,
         * necessary.
         *
         * As always, write the discrete analog to the Lagrangian, L, above, and
         * use the discrete Euler-Lagrangian equations to solve for the update rules.
         * HINT: this is not that hard; you might even see this as a trick question...
         */

        Vertex v = vertices[i];

        float tempx = v.x + tstep / v.mass * 
            (v.mass * v.vx + tstep * v.fx);

        float tempy = v.y + tstep / v.mass *
            (v.mass * v.vy + tstep * v.fy);

        vertices[i].vx = (tempx - v.x) / tstep;
        vertices[i].vy = (tempy - v.y) / tstep;

        vertices[i].x = tempx;
        vertices[i].y = tempy;

         

        /****************************** END TODO ****************************/
        
        ke += 1.0 / 2.0 * ( vertices[i].vx * vertices[i].vx
                            + vertices[i].vy * vertices[i].vy );
    }

    // Set our total variables, which are to be drawn as text soon
    float total = ke + pe;
    min_total = (total < min_total) ? total : min_total;
    max_total = (total > max_total) ? total : max_total;

    total_time += tstep;

    draw_text();

    glutSwapBuffers();
}

void keyPressed(unsigned char key, int x, int y)
{
    if(key == 'q')
    {
        exit(0);
    }
}

/* Checks if mouse click is near a "clickrad" radius of a vertex, and if
 * it is, return that vertex's index in the list */
int getClickedPoint(float x, float y)
{    
    for(int i = 0; i < vertices.size(); i++)
    {
        if(vertices[i].x > x - clickRad && vertices[i].x < x + clickRad &&
           vertices[i].y > y - clickRad && vertices[i].y < y + clickRad)
            return i;
    }
    
    return -1;
}

void mousePressed(int button, int state, int x, int y)
{
    // Gets the vertex that we click and designate it as being dragged
    if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        int index = getClickedPoint(((float) x - (float) xres / 2.0) * mouseScaleX,
                                    -((float) y - (float) yres / 2.0) * mouseScaleY);
        
        if(index == -1)
            return;
        
        mouseX = x;
        mouseY = y;
        
        draggedPoint = index;
        dragging = true;
    }
    else if(button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    {
        dragging = false;
    }
}

void mouseMoved(int x, int y)
{
    // Move the dragged vertex appropriate
    if(dragging)
    {
        float dx = (float) (x - mouseX) * mouseScaleX;
        float dy = (float) (y - mouseY) * mouseScaleY;
        
        vertices[draggedPoint].x += dx;
        vertices[draggedPoint].y -= dy;
        
        mouseX = x;
        mouseY = y;
    }
}

int main(int argc, char* argv[])
{
    std::cout << "Parsing in obj file: " << argv[1] << "\n";
    bool success = parse_OBJ(argv[1]);
    
    if(!success)
    {
        exit(1);
    }
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
    glutInitWindowSize(xres, yres);
    glutInitWindowPosition(0, 0);
    
    glutCreateWindow("Elastic Demo");
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyPressed);
    glutMouseFunc(mousePressed);
    glutMotionFunc(mouseMoved);

    glutIdleFunc(display);
    glutMainLoop();
}
