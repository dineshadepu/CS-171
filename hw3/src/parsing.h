#include <iostream>
#include <cstdio>
#include <string>
#include "structs.h"

using namespace std;

/* Take vertex/normal arrays of current object, 
 * and based on faces, store into large
 * buffer arrays to be fed to OpenGL */
void storeBuffers(Object &object);

/* Parse object file for object parameters, store them, 
and return a corresponding object */
Object parse_obj(string line);

/* Take line with vector parameters, store as a 
transformation matrix, return it */
Transform parse_mat(string line);

/* Parse each object copy for its parameters and store */
void parse_objCopy(Object &object, string line);

/* parse camera section and store specific parameters */
void parse_camera(Camera &cam, vector<Light> &lights, 
    string line, bool &checker);

void parse_all(string &line, vector<Object> &objects, vector<Object> &objectsCopies,
    bool &checker, bool &checker2, Camera &camera, vector<Light> &lights);