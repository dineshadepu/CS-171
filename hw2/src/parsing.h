#include <iostream>
#include <cstdio>
#include <string>
#include "structs.h"
#include <Eigen/Dense>

/* Parse object file for object parameters, store them, 
and return a corresponding object */
Object parse_obj(std::string line);

/* Take line with vector parameters, store as a 
transformation matrix, return it */
Eigen::Matrix4d parse_mat(std::string line);

/* Parse each object copy for its parameters and store */
void parse_objCopy(Object &object, std::string line);

/* parse camera section and store specific parameters */
void parse_camera(Camera &cam, std::vector<Light> &lights, 
    std::string line, bool &checker);