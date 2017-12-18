#include <iostream>
#include <cstdio>
#include <string>
#include "structs.h"
#include <Eigen/Dense>

Object parse_obj(std::string line);
Eigen::Matrix4d parse_mat(std::string line);