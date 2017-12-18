#include <Eigen/Dense>
#include "matrixFunctions.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <Eigen/Dense>
#include <cmath>
#include "structs.h"

using namespace std;

/* Translation Matrix */
Eigen::Matrix4d createTransMat(float x, float y, float z)
{
    Eigen::Matrix4d m;
    m << 1, 0, 0, x,
    0, 1, 0, y,
    0, 0, 1, z,
    0, 0, 0, 1;
    return m;
}

/* Rotation matrix */
Eigen::Matrix4d createRotMat(float x, float y, float z, float r)
{
	// Normalize
    float temp = sqrt(x*x + y*y + z*z);
    x = x/temp;
    y = y/temp;
    z = z/temp;

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

/* Scaling matrix */
Eigen::Matrix4d createScaleMat(float x, float y, float z)
{
    Eigen::Matrix4d m;
    m << x, 0, 0, 0,
    0, y, 0, 0,
    0, 0, z, 0,
    0, 0, 0, 1;
    return m;
}

/* Perspective Projection Matrix */
Eigen::Matrix4d makePerspProj(float n, float f, float l, float r, 
    float t, float b)
{
    Eigen::Matrix4d m;
    m << 2*n/(r-l), 0, (r+l)/(r-l), 0,
    0, 2*n/(t-b), (t+b)/(t-b), 0,
    0, 0, -(f+n)/(f-n), -2*f*n/(f-n),
    0, 0, -1, 0;
    return m;
}

/* Apply a geometric transformation on a point */
Vertex geomTransform(Vertex vert, Eigen::Matrix4d tm)
{
    Vertex v;
    Eigen::Vector4d vec;
    Eigen::Vector4d transfPoint;

    vec << vert.x, vert.y, vert.z, 1;
    // transform point vector by geometric transformation passed in
    transfPoint = tm*vec;

    v.x = transfPoint(0);
    v.y = transfPoint(1);
    v.z = transfPoint(2);

    return v;
}

/* Apply a geometric transformation on a point */
Vertex normalTransform(Vertex sn, Eigen::Matrix4d tm)
{
    Vertex normal;
    Eigen::Vector4d vec, transfPoint;
    Eigen::Vector3d vec2;

    vec << sn.x, sn.y, sn.z, 1;

    // transform point vector by geometric transformation passed in
    transfPoint = tm*vec;

    // calculate norm
    float normal_norm = transfPoint.head(3).norm();

    // normalize each point and store
    normal.x = transfPoint(0)/normal_norm;
    normal.y = transfPoint(1)/normal_norm;
    normal.z = transfPoint(2)/normal_norm;

    return normal;
}

/* Convert a given point from world space to camera to cartesian NDC */
Vertex worldToNDC(Vertex vert, Eigen::Matrix4d perspProj, Eigen::Matrix4d w2c)
{
    Vertex v;
    Eigen::Vector4d vec;

    vec << vert.x, vert.y, vert.z, 1;

    // Transformation from world to camera space
    // to have [x_c  y_c  z_c  1]
    // then transform to homogeneous NDC using perspective projection 
    // which is [x_ndc  y_ndc  z_ndc  w_ndc] where w_ndc = -z_c
    Eigen::Vector4d transfPoint = perspProj*w2c*vec;

    // Store Cartesian NDC by scaling by -z_c
    v.x = transfPoint(0)/transfPoint(3);
    v.y = transfPoint(1)/transfPoint(3);
    v.z = transfPoint(2)/transfPoint(3);
    
    return v;
}

Point NDCToScreen(Vertex vert, int xres, int yres)
{
	// Map NDC points to screen coordinates
    // Since the box is in the range [-1,1] in both x and y,
    // points will be in this range. We can use this range to normalize
    // points against xres and yres.
    // We scale to the pixel gride dimensions while keeping 0,0 center
    // then shifting so that 0,0 is corner and using int to round
    // Also flips y axis for screen.
	Point result;
	result.x = int(vert.x*(xres/2)+(xres/2))-1;
	result.y = int(vert.y*(-yres/2)+(yres/2))-1;
	return result;
}

