#include <iostream>
#include <cstdio>
#include <string>
#include "structs.h"

using namespace std;

/* Parse object file for object parameters, store them, 
and return a corresponding object */
Mesh_Data* parse_obj(string &line);

/* Take line with vector parameters, store as a 
transformation matrix, return it */
Transform parse_mat(string &line);

/* Parse each object copy for its parameters and store */
void parse_objCopy(Object &object, string &line);

/* parse camera section and store specific parameters */
void parse_camera(Camera &cam, vector<Light> &lights, 
    string &line, bool &checker);

void parse_all(string &line, vector<Mesh_Data*> *meshes, vector<Object> &objectsCopies,
    bool &checker, bool &checker2, Camera &camera, vector<Light> &lights);