#include <Eigen/Dense>

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