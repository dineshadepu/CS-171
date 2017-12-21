#include <Eigen/Dense>
#include "matrixFunctions.hpp"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <Eigen/Dense>
#include <cmath>
#include "Utilities.hpp"

using namespace std;

/* Translation Matrix */
Eigen::Matrix4f createTransMat(float x, float y, float z)
{
    Eigen::Matrix4f m;
    m << 1, 0, 0, x,
    0, 1, 0, y,
    0, 0, 1, z,
    0, 0, 0, 1;
    return m;
}

// /* Rotation matrix */
Eigen::Matrix4f createRotMat(float x, float y, float z, float r)
{
	// Normalize
	float temp = sqrtf(x*x + y*y + z*z);
	x = x/temp;
	y = y/temp;
	z = z/temp;

    Eigen::Matrix4f m;
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
Eigen::Matrix4f createScaleMat(float x, float y, float z)
{
    Eigen::Matrix4f m;
    m << x, 0, 0, 0,
    0, y, 0, 0,
    0, 0, z, 0,
    0, 0, 0, 1;
    return m;
}