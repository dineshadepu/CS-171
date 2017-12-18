#include <Eigen/Dense>
#include "matrixFunctions.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <Eigen/Dense>
#include <cmath>

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
	0, 0, -(f+n)/(f-n), 2*f*n/(f-n),
	0, 0, -1, 0;
	return m;
}