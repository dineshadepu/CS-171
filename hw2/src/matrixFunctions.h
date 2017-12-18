#include <Eigen/Dense>
#include "structs.h"

/*** Functions for creating matrices ***/

/* Translation Matrix */
Eigen::Matrix4d createTransMat(float x, float y, float z);

/* Rotation matrix */
Eigen::Matrix4d createRotMat(float x, float y, float z, float r);

/* Scaling matrix */
Eigen::Matrix4d createScaleMat(float x, float y, float z);

/* Perspective Projection Matrix */
Eigen::Matrix4d makePerspProj(float n, float f, float l, float r, 
							  float t, float b);

/* Applying a Geometric Transformation to a point */
Vertex geomTransform(Vertex vert, Eigen::Matrix4d tm);

/* Applying the inverse transpose geometric transformation to a surface normal */
Vertex normalTransform(Vertex sn, Eigen::Matrix4d tm);

/* Converting a point from World Space to Cartesian NDC */
Vertex worldToNDC(Vertex vert, Eigen::Matrix4d perspProj, Eigen::Matrix4d w2c);

/* Mapping a vertex to screen coordinates */
Point NDCToScreen(Vertex vert, int xres, int yres);