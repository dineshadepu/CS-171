/* CS/CNS 171
 * Fall 2017
 * Written by Eshan Govil
 *
 * The following program animates the bunny smoothing
 */

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <math.h>
#define _USE_MATH_DEFINES

#include <cmath>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "structs.h"
#include "halfedge.h"

using namespace std;

/*** Hardcode camera, frustum, lights, and material from scene_bunny.txt ***/

float cam_position[] = {0, 0, 3};
float cam_orientation[] = {0, 1, 0};
float cam_orientation_angle = 0;

const float near_param = 1.0, far_param = 10.0,
            left_param = -0.5, right_param = 0.5,
            top_param = 0.5, bottom_param = -0.5;

const float light1_color[3] = {0.0, 0.6, 0.8};
const float light1_position[4] = {5, 5, 5, 1};
const float light1_attenuation = 0.0;

const float light2_color[3] = {0.0, 0.1, 0.2};
const float light2_position[4] = {-5, -5, -5, 1};
const float light2_attenuation = 0.0;

const float ambient_reflect[3] = {0.5, 0.5, 0.5};
const float diffuse_reflect[3] = {0.5, 0.5, 0.5};
const float specular_reflect[3] = {0.5, 0.5, 0.5};
const float shininess = 0.1;

/* Global variables */

const float fstep = 1; // frame step
const float t_val = 0; // for Catmull-Rom splines
float currFrame = 0; // current frame that will be updated and used to render
float totFrames = 20; // total bunny frames

vector<int> frameVals; // indices of keyframes
vector<Mesh_Data*> *meshes; //  vectors of all keyframe objects and interpolated objects
vector<Face*> *faces;

// for storing the object that will be rendered at that frame
vector<Vertex> vertex_buffer;
vector<Vec3f> normal_buffer;

Eigen::Matrix4d last_rotation, current_rotation;
Quaternion last_rotation_quat, current_rotation_quat;
Point p_start, p_current; // points
int xres, yres;

int mouse_x, mouse_y;
float mouse_scale_x, mouse_scale_y;

const float step_size = 0.2;
const float x_view_step = 90.0, y_view_step = 90.0;
float x_view_angle = 0, y_view_angle = 0;

bool is_pressed = false;
bool wireframe_mode = false;

/* Declare relevant functions */
void draw_text();

void storeBuffers();
void deleteMeshes();
void interpolateVals();
void parse_obj(string &line, int &meshindex);
void initMeshes();
void storeBunnies();
void storeObjFiles();

void init_lights();
void set_lights();
void draw_objects();
Quaternion Compute_Rotation_Quaternion(Point p_current, Point p_start);
Quaternion Get_Current_Rotation_Quat();
void calc_vertex_normal(HEV* v);
void halfedgeCalcs();

Quaternion Get_Current_Rotation_Quat();
Quaternion Compute_Rotation_Quaternion(Point p_current, Point p_start);

void init(void)
{
    glShadeModel(GL_SMOOTH);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    glFrustum(left_param, right_param,
            bottom_param, top_param,
            near_param, far_param);

    init_lights();

    // Set material properties of bunny, hard coded
    glMaterialfv(GL_FRONT, GL_AMBIENT, ambient_reflect);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_reflect);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular_reflect);
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);
    
    glMatrixMode(GL_MODELVIEW);

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

void reshape(int width, int height)
{
    height = (height == 0) ? 1 : height;
    width = (width == 0) ? 1 : width;
    
    glViewport(0, 0, width, height);

    /* The following two lines are specific to updating our mouse interface
     * parameters. Details will be given in the 'mouse_moved' function */
    mouse_scale_x = (float) (right_param - left_param) / (float) width;
    mouse_scale_y = (float) (top_param - bottom_param) / (float) height;
    
    glutPostRedisplay();
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glRotatef(-cam_orientation_angle * 180.0 / M_PI,
              cam_orientation[0], 
              cam_orientation[1], 
              cam_orientation[2]);
    glTranslatef(-cam_position[0], -cam_position[1], -cam_position[2]);

    // set the current quaternion rotation of the arcball
    // then set the lights and draw the bunny based on current frame
    glPushMatrix();
    glMultMatrixd(Get_Current_Rotation_Quat().Get_Rotation_Matrix_Quat().data());
    set_lights();
    draw_objects();
    glPopMatrix();

    draw_text();
    
    glutSwapBuffers();
}

/* Hard code the lights in from bunny scene */
void init_lights()
{
    glEnable(GL_LIGHTING);
    
    int light_id = GL_LIGHT0 + 0;
    glEnable(light_id);

    glLightfv(light_id, GL_AMBIENT, light1_color);
    glLightfv(light_id, GL_DIFFUSE, light1_color);
    glLightfv(light_id, GL_SPECULAR, light1_color);
    glLightf(light_id, GL_QUADRATIC_ATTENUATION, light1_attenuation);

    light_id = GL_LIGHT0 + 1;
    glEnable(light_id);
    
    glLightfv(light_id, GL_AMBIENT, light2_color);
    glLightfv(light_id, GL_DIFFUSE, light2_color);
    glLightfv(light_id, GL_SPECULAR, light2_color);
    glLightf(light_id, GL_QUADRATIC_ATTENUATION, light2_attenuation);
}

/* Hard code the light positions */
void set_lights()
{
    int light_id = GL_LIGHT0 + 0;
    glLightfv(light_id, GL_POSITION, light1_position);

    light_id = GL_LIGHT0 + 1;
    glLightfv(light_id, GL_POSITION, light2_position);
}

/* This function has OpenGL render our bunny to the display screen. */
void draw_objects()
{
    glPushMatrix();
    
    /* Hard code the initial bunny transformations from the bunny scene:
       translate, rotate, scale */
    glTranslatef(0.0, 0.0, 0.0);
    glRotatef(0.0, 1.0, 0.0, 0.0);
    glScalef(1.0, 1.0, 1.0);
    
    /* How we tell OpenGL to render geometry for us. 
     * Use points to buffer arrays */
    glVertexPointer(3, GL_FLOAT, 0, &vertex_buffer[0]);
    glNormalPointer(GL_FLOAT, 0, &normal_buffer[0]);
    
    int buffer_size = vertex_buffer.size();
    
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

    glPopMatrix();
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

// Forward incrementing the frame with no wrap around
void updateF()
{
    /* FRAME UPDATE RULES */
    if (currFrame < totFrames)
    {
        currFrame += fstep;
        storeBuffers();
        glutPostRedisplay();
    }
}

// Backward incrementing the frame with no wrap around
void updateB()
{
    /* FRAME UPDATE RULES */
    if (currFrame > 0)
    {
        currFrame -= fstep;
        // given the current frame, store the appropriate buffer
        storeBuffers();
        glutPostRedisplay();
    }
}

template<typename T>
string tostr(const T& t)
{
    ostringstream os;
    os << t;
    return os.str();
}

void draw_text()
{
    glColor3f(0, 1, 0);

    string frame_str = "Frame: " + tostr(currFrame);

    glRasterPos2f(1,1.3);
    for(int i = 0; i < frame_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, frame_str[i]);
}

void key_pressed(unsigned char key, int x, int y)
{
    if(key == 'q')
    {
        deleteMeshes();
        exit(0);
    }
    else if (key == 'i')
    {
        updateF();
    }
    else if (key == 'o')
    {
        updateB();
    }
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
            cam_position[0] += step_size * sin(x_view_rad);
            cam_position[2] -= step_size * cos(x_view_rad);
            glutPostRedisplay();
        }
        /* 'a' for step left
         */
        else if(key == 'a')
        {
            cam_position[0] -= step_size * cos(x_view_rad);
            cam_position[2] -= step_size * sin(x_view_rad);
            glutPostRedisplay();
        }
        /* 's' for step backward
         */
        else if(key == 's')
        {
            cam_position[0] -= step_size * sin(x_view_rad);
            cam_position[2] += step_size * cos(x_view_rad);
            glutPostRedisplay();
        }
        /* 'd' for step right
         */
        else if(key == 'd')
        {
            cam_position[0] += step_size * cos(x_view_rad);
            cam_position[2] += step_size * sin(x_view_rad);
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

// do this for all meshes
// assuming all vertices for all meshes have been calculated and stored
void halfedgeCalcs()
{
    cerr << "doing half edge calculations for all meshes" << endl;

    for (int i = 0; i < meshes->size(); i++)
    {
        cerr << ".";
        // once mesh->vertices and mesh->vertices are already filled
        vector<HEV*> *hevs = new vector<HEV*>();
        vector<HEF*> *hefs = new vector<HEF*>();

        build_HE(meshes->at(i), faces, hevs, hefs);

        // Mesh normals vector is being set to size of vertices
        meshes->at(i)->normals->resize(hevs->size());
        
        // For each vertex, do calculations
        for (int j = 1; j < hevs->size(); j++)
        {
            Vec3f *normal;
            // Calculate the normal of the vertex
            calc_vertex_normal(hevs->at(j));

            // create a normal for each vertex
            normal = new Vec3f;

            normal->x = hevs->at(j)->normal.x;
            normal->y = hevs->at(j)->normal.y;
            normal->z = hevs->at(j)->normal.z;

            // store normal from halfedge into mesh->normals
            meshes->at(i)->normals->at(hevs->at(j)->index) = normal;
        }

        // done with normal and operator calculations as well as 
        // vertex and normal updating, so done with these halfedges
        delete_HE(hevs, hefs);
    }
    cerr << endl;

    cerr << "finished half edge calculations" << endl;
}


/* Take vertex/normal arrays of current object, 
 * and based on faces, store into large
 * buffer arrays to be fed to OpenGL */
void storeBuffers()
{
    // clear buffers so that we don't keep building
    // on original mesh
    vertex_buffer.clear();
    normal_buffer.clear();

    // store buffer for each face
    for (int i = 0; i < faces->size(); i++)
    {
        Face *face = faces->at(i);

        // since meshes are pointers, store values at the pointers and statically
        // allocate the vertices/normals into buffer to be passed to OpenGL
        vertex_buffer.push_back(*meshes->at(currFrame)->vertices->at(face->idx1));
        vertex_buffer.push_back(*meshes->at(currFrame)->vertices->at(face->idx2));
        vertex_buffer.push_back(*meshes->at(currFrame)->vertices->at(face->idx3));

        normal_buffer.push_back(*meshes->at(currFrame)->normals->at(face->idx1));
        normal_buffer.push_back(*meshes->at(currFrame)->normals->at(face->idx2));
        normal_buffer.push_back(*meshes->at(currFrame)->normals->at(face->idx3));
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

        // delete the vectors
        delete meshes->at(i)->vertices;
        delete meshes->at(i)->normals;
        delete meshes->at(i);
    }

    // outermost: delete the meshes array itself
    delete meshes;

    // delete each face
    for (int i = 0; i < faces->size(); i++)
        delete faces->at(i);
    delete faces;
}

void interpolateVals()
{
    // make B matrix for t = 0
    float s = 0.5*(1 - t_val);
    Eigen::Matrix4d B;
    B << 0, 1, 0, 0,
         -s, 0, s, 0,
         2*s, s-3, 3-2*s, -s,
         -s, 2-s, s-2, s;    

    // do this for each vertex in array
    // since the order is consistent across meshes
    for (int j = 1; j < meshes->at(0)->vertices->size(); j++)
    {
        // since we don't wrap around, don't need to do calculations from frame 20
        for (int i = 0; i < frameVals.size()-1; i++)
        {
            // vectors of point vectors
            Eigen::Vector4d transvx, transvy, transvz;

            // indices for primary keyframe positions
            int iminus1, icurrent, iplus1, iplus2;
            icurrent = frameVals[i];

            // handle index assignment for cases before general case
            // by just having the edge case repeat itself because
            // we don't wrap around
            // that is, before frame 5, first and second control points is 0
            // and for frame 15, third and fourth keyframe is 20
            if (i == 0)
            {
                iminus1 = frameVals[i];
                iplus1 = frameVals[i+1];
                iplus2 = frameVals[i+2];
            }
            else if (i == frameVals.size()-2)
            {
                iminus1 = frameVals[i-1];
                iplus1 = frameVals[i+1];
                iplus2 = frameVals[i+1];
            }
            else
            {
                iminus1 = frameVals[i-1];
                iplus1 = frameVals[i+1];
                iplus2 = frameVals[i+2];
            }

            // now store values into point vectors of four control points
            // for all x, y, and z of that column of vertices
            transvx << meshes->at(iminus1)->vertices->at(j)->x, 
                meshes->at(icurrent)->vertices->at(j)->x,
                meshes->at(iplus1)->vertices->at(j)->x,
                meshes->at(iplus2)->vertices->at(j)->x;

            transvy << meshes->at(iminus1)->vertices->at(j)->y, 
                meshes->at(icurrent)->vertices->at(j)->y,
                meshes->at(iplus1)->vertices->at(j)->y,
                meshes->at(iplus2)->vertices->at(j)->y;

            transvz << meshes->at(iminus1)->vertices->at(j)->z, 
                meshes->at(icurrent)->vertices->at(j)->z,
                meshes->at(iplus1)->vertices->at(j)->z,
                meshes->at(iplus2)->vertices->at(j)->z;

            // get current frame to start at and end frame
            float du, u, normer, qs, qx, qy, qz;

            int start = icurrent;
            int end = iplus1;

            // to go from 0 to 1 across all frames between start and end
            // move in increments of du = 1/(end-start)
            du = 1.0 / (end-start);
            u = du;

            // interpolate all values in between start and end
            for (int k = start+1; k < end; k++)
            {
                Eigen::Vector4d uvec, umultB;
                
                // make u vector and start curve function calculation
                uvec << 1, u, u*u, u*u*u;
                umultB = uvec.transpose() * B;

                // finish the cardinal curve function using dot product and store f in vertex vectors
                Vertex *vert = new Vertex;
                vert->x = umultB.dot(transvx);
                vert->y = umultB.dot(transvy);
                vert->z = umultB.dot(transvz);

                meshes->at(k)->vertices->push_back(vert);

                u += du;
            }
        }
    }
}

void parse_obj(string &line, int &meshindex)
{
    string vf, vn1, vn2, vn3;
    istringstream iss(line);

    if (line[0] == 'v') // if it's a vertex line
    {
        Vertex *vert = new Vertex;
        iss >> vf >> vert->x >> vert->y >> vert->z;
        meshes->at(meshindex)->vertices->push_back(vert);
    }
    // since faces are the same for all frames
    // store faces just once for the first mesh
    else if (line[0] == 'f' && meshindex == 0) // if it's a face line
    {
        Face *face = new Face;
        iss >> vf >> face->idx1 >> face->idx2 >> face->idx3;
        faces->push_back(face); // add face to vector in object
    }
}

void initMeshes()
{
    meshes = new vector<Mesh_Data*>();
    faces = new vector<Face*>();
    // do 21 meshes because we don't wrap around so keyframe 20 is separate
    // and we have 5 keyframes + 16 missing frames
    for (int i = 0; i < 21; i++)
    {
        // dynamically allocate the mesh structure
        Mesh_Data *mesh = new Mesh_Data;
        mesh->vertices = new vector<Vertex*>();
        mesh->normals = new vector<Vec3f*>();

        mesh->vertices->push_back(NULL); // 1 index
        mesh->normals->push_back(NULL); // 1 index

        meshes->push_back(mesh);
    }    
}

void storeBunnies()
{
    for (int i = 0; i < 21; i += 5)
    {
        string filename, objline;

        frameVals.push_back(i);

        if (i < 10)
            filename = "keyframes/bunny0" + tostr(i) + ".obj";
        else
            filename = "keyframes/bunny" + tostr(i) + ".obj";

        fstream myobjfile(filename.c_str());

        if(myobjfile.is_open())
        {
            while(getline(myobjfile, objline))
            {
                parse_obj(objline, i);
            }
        }
        myobjfile.close();
    }
}

void storeObjFiles()
{
    for (int i = 0; i < meshes->size(); i++)
    {
        ofstream myfile;
        string filename;
        filename = "myframes/bunny" + tostr(i) + ".obj";
        myfile.open(filename.c_str());

        // write that meshes vertices and faces to the files
        for (int j = 1; j < meshes->at(i)->vertices->size(); j++)
        {
            Vertex *vert = meshes->at(i)->vertices->at(j);
            myfile << "v " << vert->x << " " << vert->y
                << " " << vert->z << "\n";
        }
        for (int j = 0; j < faces->size(); j++)
        {
            Face *f = faces->at(j);
            myfile << "f " << f->idx1 << " " << f->idx2
                << " " << f->idx3 << "\n";
        }

        myfile.close();
    }
}

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        cerr << "\nERROR: Incorrect number of arguments." << endl;
        exit(1);
    }

    //fstream myobjfile(argv[1]);
    istringstream(argv[1]) >> xres;
    istringstream(argv[2]) >> yres;

    initMeshes();

    string objline;

    // Parse object and store in vectors
    cerr << "Loading all bunny keyframes" << endl;

    storeBunnies();

    cerr << "Done object parsing and storing...interpolating values..." << endl;

    interpolateVals();

    cerr << "done interpolation...calculate normals and store..." << endl;

    halfedgeCalcs();

    storeObjFiles();

    cerr << "Storing Initial Vertex and Normal Buffers..." << endl;

    storeBuffers();

    cerr << "Rendering and Displaying..." << endl;
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(xres, yres);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Bunny Smoother Display");
    
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse_pressed);
    glutMotionFunc(mouse_moved);
    glutKeyboardFunc(key_pressed);
    glutMainLoop();
}
