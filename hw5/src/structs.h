#ifndef STRUCTS_H
#define STRUCTS_H

#include <vector>
#include <cstdlib>
#include <string>
#include <Eigen/Dense>

using namespace std;

/* Primary structs for normals, vertices and faces */
struct Vec3f
{
	float x, y, z;
};

struct Vertex
{
    float x, y, z;
};

struct Face
{
    int idx1, idx2, idx3;
};

struct Mesh_Data
{
    vector<Vertex*> *vertices;
    vector<Vec3f*> *normals;

    string label;
};

/* The following struct is used for representing a point light. */
struct Light
{
    /* Index 0 has the x-coordinate
     * Index 1 has the y-coordinate
     * Index 2 has the z-coordinate
     * Index 3 has the w-coordinate
     */
    float position[4];
    
    /* Index 0 has the r-component
     * Index 1 has the g-component
     * Index 2 has the b-component
     */
    float color[3];

    float k;
};

struct Point
{
    int x, y;
};

struct PointNDC
{
    float x, y, z;
};

struct Transform
{
    /* if label is "trans"
     * we have translation[3]
     * if label is "rot"
     * we have rotation[3] and rot_angle
     * if label is "scal"
     * we have scaling[3]
     */
    string label;
    float transformation[3];
    float rotation_angle; // only used for rot
};

struct Object
{
    /* See the note above and the comments in the 'draw_objects' and
     * 'create_cubes' functions for details about these buffer vectors.
     */
    Mesh_Data *mesh;

    /* Vectors with all the faces represented by vertices
     * and normals in order to be fed to OpenGL to render */
    vector<Vertex> vertex_buffer;
    vector<Vec3f> normal_buffer;
    vector<Transform> transformations;
    
    /* Index 0 has the r-component
     * Index 1 has the g-component
     * Index 2 has the b-component
     */
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float shininess;  
};

struct Quaternion
{
    float qs, qx, qy, qz;
    /* Using current rotation quaternion, make rotation matrix */
    Eigen::Matrix4d Get_Rotation_Matrix_Quat()
    {
        Eigen::Matrix4d m;
        m << 1-2*qy*qy-2*qz*qz, 2*(qx*qy-qz*qs),
            2*(qx*qz+qy*qs), 0,
            2*(qx*qy+qz*qs), 1-2*qx*qx-2*qz*qz,
            2*(qy*qz-qx*qs), 0,
            2*(qx*qz-qy*qs), 2*(qy*qz+qx*qs),
            1-2*qx*qx-2*qy*qy, 0, 
            0, 0, 0, 1;
        return m;
    }
};

/* The following are the camera specifications and parameters. 
 */
struct Camera
{
    /* Index 0 has the x-coordinate
     * Index 1 has the y-coordinate
     * Index 2 has the z-coordinate
     */
    float position[3];
    float orientation_axis[3];

    /* Angle in degrees. */ 
    float orientation_angle;

    float near, far, left, right, top, bottom;
};

#endif
