/* The following 2 headers contain all the main functions, data structures, and
 * variables that allow for OpenGL development.
 */
#include <GL/glew.h>
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
#include <Eigen/Sparse>
#include "halfedge.h"

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
Quaternion Compute_Rotation_Quaternion(Point p_current, Point p_start);
Quaternion Get_Current_Rotation_Quat();

void mouse_pressed(int button, int state, int x, int y);
void mouse_moved(int x, int y);
void key_pressed(unsigned char key, int x, int y);

void calc_vertex_normal(HEV* vertex);
void halfedgeCalcs();
Eigen::SparseMatrix<double> build_F_operator(vector<HEV*> *vertices);
void solve(vector<HEV*> *vertices);
double calcArea(HE *he);
double calcAlpha(HE *he);
double calcBeta(HE *he);
void storeBuffers();
void deleteMeshes();

///////////////////////////////////////////////////////////////////////////////////////////////////

/* STRUCTS DEFINED IN SEPARATE HEADER FILE */

///////////////////////////////////////////////////////////////////////////////////////////////////

/* CAMERA STORED IN SEPARATE STRUCTURE */

///////////////////////////////////////////////////////////////////////////////////////////////////

/* GLOBAL VARIABLES */

Camera camera;
vector<Light> lights; // vector of lights
vector<Mesh_Data*> *meshes; //  vector of objects
vector<Object> objectsCopies; // vector of all copies
Eigen::Matrix4d last_rotation, current_rotation;
Quaternion last_rotation_quat, current_rotation_quat;
Point p_start, p_current; // points
int xres, yres;
float h; // time step

int mouse_x, mouse_y;
float mouse_scale_x, mouse_scale_y;

const float step_size = 0.2;
const float x_view_step = 90.0, y_view_step = 90.0;
float x_view_angle = 0, y_view_angle = 0;

bool is_pressed = false;
bool wireframe_mode = false;
bool initial_toggle = false;
int reserve = 7;
double minArea = 1e-9;

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

    /* Setting quaternion rotations to identity */
    last_rotation_quat.qs = 1;
    last_rotation_quat.qx = 0;
    last_rotation_quat.qy = 0;
    last_rotation_quat.qz = 0;

    current_rotation_quat.qs = 1;
    current_rotation_quat.qx = 0;
    current_rotation_quat.qy = 0;
    current_rotation_quat.qz = 0;
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
    rotation = Get_Current_Rotation_Quat().Get_Rotation_Matrix_Quat();

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
        last_rotation_quat = Get_Current_Rotation_Quat();
        current_rotation_quat.qs = 1;
        current_rotation_quat.qx = 0;
        current_rotation_quat.qy = 0;
        current_rotation_quat.qz = 0;
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
        deleteMeshes();
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
    else if(key == 'i') // FOR IMPLICIT FAIRING
    {
        // Smooth down based on h
        // builds F operater, solves given current vertices
        // then recomputes normals
        halfedgeCalcs();
        // clears current buffers and stores new vertices
        // and normals
        storeBuffers();
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

/* Calculate the normal of a vertex using neighboring face normals */
void calc_vertex_normal(HEV* v)
{
    // initialize a normal vector structure to the zero vector
    Vec3f normal;
    normal.x = 0.0;
    normal.y = 0.0;
    normal.z = 0.0;

    // go through all faces and calculate face normals
    // and add to normal vector for vertex normal
    HE *he = v->out; 
    do
    {
        HEF *f = he->face;

        HEV *v1 = f->edge->vertex; 
        HEV *v2 = f->edge->next->vertex;
        HEV *v3 = f->edge->next->next->vertex; 

        // compute normal of plane of face using cross product
        // of (v2 - v1) x (v3 - v1)
        // with a vector from Eigen to ease computation
        Eigen::Vector3d a, b, face_normal, temp;
        float n;
        a << v2->x - v1->x, v2->y - v1->y, v2->z - v1->z;
        b << v3->x - v1->x, v3->y - v1->y, v3->z - v1->z;
        temp = a.cross(b);
        n = temp.norm();
        face_normal = temp / n; // normalize face normal

        // compute area of the triangle face
        // which is half the length of the norm
        float face_area = 0.5 * n;

        // accumulate onto our normal vector of the vertex
        normal.x += face_normal(0) * face_area;
        normal.y += face_normal(1) * face_area;
        normal.z += face_normal(2) * face_area;

        he = he->flip->next;
    }
    while (he != v->out);

    // normalize n
    float n = sqrt(normal.x * normal.x + normal.y * normal.y
        + normal.z * normal.z);
    normal.x = normal.x / n;
    normal.y = normal.y / n;
    normal.z = normal.z / n;

    v->normal = normal;
}

void halfedgeCalcs()
{
    // do these calculations for each mesh
    for (int i = 0; i < meshes->size(); i++)
    {
        // once mesh->vertices and mesh->vertices are already filled
        vector<HEV*> *hevs = new vector<HEV*>();
        vector<HEF*> *hefs = new vector<HEF*>();

        cerr << "Building Halfedges..." << endl;

        build_HE(meshes->at(i), hevs, hefs);

        cerr << "Finished Halfedge Creation..." << endl;

        // if toggle is true, that means we have done the initial mesh display
        // and now need to solve the operator and smooth the mesh with new points
        if (initial_toggle)
        {
            cerr << "Building and Solving the Operator..." << endl;

            solve(hevs);

            cerr << "Built Operator...Now Recomputing Normals and Updating Meshes..." << endl;
        }
        else
        {
            cerr << "Now Doing Initial Normal Calculations and Storing in Buffers..." << endl;
            // Mesh normals vector is being set to size of vertices
            meshes->at(i)->normals->resize(hevs->size());
        }
        
        // For each vertex, do calculations
        for (int j = 1; j < hevs->size(); j++)
        {
            Vec3f *normal;
            // Calculate the normal of the vertex
            calc_vertex_normal(hevs->at(j));

            // if the toggle is on, that means initial normals have already
            // been calculated, so now we just need to update the existing
            // normals and update vertices with new hev vertices
            if (initial_toggle)
            {
                // pull from stored normal
                normal = meshes->at(i)->normals->at(hevs->at(j)->index);

                normal->x = hevs->at(j)->normal.x;
                normal->y = hevs->at(j)->normal.y;
                normal->z = hevs->at(j)->normal.z;

                // Since operator and new vertices were already solved in hev vector
                // replace old coordinates with new ones. mesh will be redrawn later
                meshes->at(i)->vertices->at(j)->x = hevs->at(j)->x;
                meshes->at(i)->vertices->at(j)->y = hevs->at(j)->y;
                meshes->at(i)->vertices->at(j)->z = hevs->at(j)->z;
            }
            else // before we have implicit faired, just initializing scene
            {
                // create a normal for each vertex
                normal = new Vec3f;

                normal->x = hevs->at(j)->normal.x;
                normal->y = hevs->at(j)->normal.y;
                normal->z = hevs->at(j)->normal.z;

                // store normal from halfedge into mesh->normals
                meshes->at(i)->normals->at(hevs->at(j)->index) = normal;
            }
        }

        // done with normal and operator calculations as well as 
        // vertex and normal updating, so done with these halfedges
        delete_HE(hevs, hefs);

        cerr << "Deleted Halfedges..." << endl;
    }
}

// function to construct our F operator in matrix form 
Eigen::SparseMatrix<double> build_F_operator(vector<HEV*> *vertices) 
{
    // recall that due to 1-indexing of obj files
    int num_vertices = vertices->size()-1;

    // initialize a sparse matrix to represent our B operator 
    Eigen::SparseMatrix<double> F(num_vertices, num_vertices);

    // reserve room for 7 non-zeros per row of B 
    F.reserve(Eigen::VectorXi::Constant(num_vertices, reserve)); 

    // Build the operator by each row
    for (int i = 1; i < vertices->size(); ++i)
    { 
        HE *he = vertices->at(i)->out;
        double A = 0; // neighbor area sum

        // first do area calculations
        // iterating over all neighboring vertices
        do
        { 
            A += calcArea(he);
            he = he->flip->next;
        }
        while (he != vertices->at(i)->out);

        // if area is close to 0, don't do any calcs
        // and leave as 0 because we could divide by 0
        // otherwise do build row
        // HERE IS HOW THE OPERATOR IS BUILT:
        // Distributing the equation for the discrete Laplacian
        // and F = I-hdx, we get that each x_j gets the values
        // gets -h * 1/2A * cota+cotB for that x_j/x_i pairing
        // while x_i gets 1 - (h * 1/2A * (cota+cotB for all x_j))
        // where the 1 comes from identity matrix I
        if (A > minArea)
        {
            double total = 0;
            double val;
            // now create matrix operator
            do // iterate over all vertices adjacent to v_i
            { 
                int j = he->next->vertex->index; // get adjacent vertex to v_i 

                // do cotangent calculations for inner angles
                // using dot product of vectors coming out of angle
                // over the normal of the cross product
                double cotAlpha = calcAlpha(he);
                double cotBeta = calcBeta(he);

                // insert value at in operator row then add to total for -x_i
                val = -(h/(2 * A))*(cotAlpha + cotBeta);
                F.insert(i-1, j-1) = val;

                total += val;
 
                he = he->flip->next; 
            }
            while (he != vertices->at(i)->out);

            // insert total value of laplacian operator at x_i for (I-hdeltax)

            F.insert(i-1, i-1) = 1.0-total;
        }
    }

    F.makeCompressed();

    return F;
}

/* Calculate area of a face given it's halfedge */
double calcArea(HE *he)
{
    HEF *f = he->face;

    HEV *v1 = f->edge->vertex; 
    HEV *v2 = f->edge->next->vertex;
    HEV *v3 = f->edge->next->next->vertex;
    Eigen::Vector3d a, b, face_normal;
    a << v2->x - v1->x, v2->y - v1->y, v2->z - v1->z;
    b << v3->x - v1->x, v3->y - v1->y, v3->z - v1->z;
    face_normal = a.cross(b);

    // return area, which is 1/2 * magnitude of normal
    return 0.5*face_normal.norm();
}

/* Return cota of inner angle for a given half edge */
double calcAlpha(HE *he)
{
    Eigen::Vector3d a, b;

    HEV *v1 = he->vertex; 
    HEV *v2 = he->next->vertex;
    HEV *v3 = he->next->next->vertex;
    a << v2->x - v3->x, v2->y - v3->y, v2->z - v3->z;
    b << v1->x - v3->x, v1->y - v3->y, v1->z - v3->z;

    return a.dot(b) / a.cross(b).norm();
}

/* Return cotB of inner angle for a given half edge */
double calcBeta(HE *he)
{
    Eigen::Vector3d a, b;

    HEV *v1 = he->flip->vertex; 
    HEV *v2 = he->flip->next->vertex;
    HEV *v3 = he->flip->next->next->vertex;
    a << v2->x - v3->x, v2->y - v3->y, v2->z - v3->z;
    b << v1->x - v3->x, v1->y - v3->y, v1->z - v3->z;

    return a.dot(b) / a.cross(b).norm();
}

/* Use operator to solve new vertices for smoothed mesh */
void solve(vector<HEV*> *vertices)
{
    // get our matrix representation of F
    Eigen::SparseMatrix<double> F = build_F_operator(vertices); 

    // initialize Eigen's sparse solver 
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver; 

    // the following two lines essentially tailor our solver to our operator F
    solver.analyzePattern(F);
    solver.factorize(F); 

    int num_vertices = vertices->size() - 1;

    // initialize our column vectors for each axis of the vertices
    Eigen::VectorXd x0_vector(num_vertices),
        y0_vector(num_vertices), z0_vector(num_vertices);

    // populate vectors from vertex coordinates 0 indexed
    for (int i = 1; i < vertices->size(); i++)
    {
        x0_vector(i - 1) = vertices->at(i)->x;
        y0_vector(i - 1) = vertices->at(i)->y;
        z0_vector(i - 1) = vertices->at(i)->z;
    }

    // have Eigen solve for our final vertex coordinate vectors
    Eigen::VectorXd xh_vector(num_vertices),
        yh_vector(num_vertices), zh_vector(num_vertices); 

    xh_vector = solver.solve(x0_vector);
    yh_vector = solver.solve(y0_vector);
    zh_vector = solver.solve(z0_vector);

    // replaces old HEV vertex coordinates
    // with new coordinates, to be updated in meshes later
    for (int i = 1; i < vertices->size(); i++)
    {
        vertices->at(i)->x = xh_vector(i - 1);
        vertices->at(i)->y = yh_vector(i - 1);
        vertices->at(i)->z = zh_vector(i - 1);
    }
}

/* Take vertex/normal arrays of current object, 
 * and based on faces, store into large
 * buffer arrays to be fed to OpenGL */
void storeBuffers()
{
    for (int j = 0; j < objectsCopies.size(); j++)
    {
        Mesh_Data *mesh = objectsCopies[j].mesh;

        // clear buffers so that we don't keep building
        // on original mesh
        objectsCopies[j].vertex_buffer.clear();
        objectsCopies[j].normal_buffer.clear();

        // store buffer for each face
        for (int i = 0; i < mesh->faces->size(); i++)
        {
            Face *face = mesh->faces->at(i);

            // since meshes are pointers, store values at the pointers and statically
            // allocate the vertices/normals into buffer to be passed to OpenGL
            objectsCopies[j].vertex_buffer.push_back(*mesh->vertices->at(face->idx1));
            objectsCopies[j].vertex_buffer.push_back(*mesh->vertices->at(face->idx2));
            objectsCopies[j].vertex_buffer.push_back(*mesh->vertices->at(face->idx3));

            objectsCopies[j].normal_buffer.push_back(*mesh->normals->at(face->idx1));
            objectsCopies[j].normal_buffer.push_back(*mesh->normals->at(face->idx2));
            objectsCopies[j].normal_buffer.push_back(*mesh->normals->at(face->idx3));
        }
    }
}

void deleteMeshes()
{
    for (int i = 0; i < meshes->size(); i++)
    {
        // delete each vertex in vertex array of that mesh
        for (int j = 1; j < meshes->at(i)->vertices->size(); j++)
        {
            delete meshes->at(i)->vertices->at(j);
            delete meshes->at(i)->normals->at(j); // and normal of same size array
        }
        // delete each face in vertex array of that mesh
        for (int j = 0; j < meshes->at(i)->faces->size(); j++)
            delete meshes->at(i)->faces->at(j);

        // delete the vectors
        delete meshes->at(i)->vertices;
        delete meshes->at(i)->faces;
        delete meshes->at(i)->normals;
        delete meshes->at(i);
    }

    // outermost: delete the meshes array itself
    delete meshes;
}

int main(int argc, char* argv[])
{
    assert(argc == 5);

    /*** DECLARE LOCAL VARIABLES ***/

    string line;
    bool checker = false;
    bool checker2;

    fstream myfile(argv[1]);
    istringstream(argv[2]) >> xres;
    istringstream(argv[3]) >> yres;
    istringstream(argv[4]) >> h;

    meshes = new vector<Mesh_Data*>(); //  vector of objects
    //objectsCopies = new vector<Object*>(); // vector of all copies

    cerr << "Parsing Files..." << endl;

    /*** PARSE FILE AND PUT INTO DATA STRUCTURES ***/

    if (myfile.is_open())
    {
        while (!myfile.eof())
        {
            getline(myfile, line);
            parse_all(line, meshes, objectsCopies, checker, checker2, 
                camera, lights);
        }
    }
    myfile.close();

    cerr << "Done Parsing...Doing Initial Halfedge Calculations..." << endl;

    halfedgeCalcs();
    initial_toggle = true;

    cerr << "Storing Initial Vertex and Normal Buffers..." << endl;

    // populate vertex buffer and normal buffer of each object
    // in objectCopies using the now complete meshes
    storeBuffers();

    cerr << "Rendering and Displaying..." << endl;

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
