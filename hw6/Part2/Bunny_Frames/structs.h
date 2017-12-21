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
    vector<Face*> *faces;
    vector<Vec3f*> *normals;
};

struct Point
{
    int x, y;
};

struct PointNDC
{
    float x, y, z;
};

struct Transformation
{
    float x, y, z;
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


#endif
