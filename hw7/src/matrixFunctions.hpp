#include <Eigen/Dense>

/*** Functions for creating matrices ***/

/* Translation Matrix */
Eigen::Matrix4f createTransMat(float x, float y, float z);

/* Rotation matrix with quaternion*/
// Eigen::Matrix4f createRotMat(float qx, float qy, float qz, float qs);
Eigen::Matrix4f createRotMat(float x, float y, float z, float r);

/* Scaling matrix */
Eigen::Matrix4f createScaleMat(float x, float y, float z);