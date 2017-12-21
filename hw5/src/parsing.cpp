#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "structs.h"
#include "parsing.h"
#include <string>

using namespace std;

/* Parse object file for object parameters, store them, 
and return a corresponding object */
Mesh_Data* parse_obj(string &line)
{

    // dynamically allocate a mesh structure
    Mesh_Data *mesh = new Mesh_Data;

    string vf, objline, buf;
    string vn1, vn2, vn3;
    vector<string> tokens;
    int pos;
    string delim = "/"; // for the center // delimiting of faces

    istringstream iss(line);

    // using string stream, split line by whitespaces and store
    while(iss >> buf) 
        tokens.push_back(buf);

    // dynamically allocate the mesh vectors
    mesh->vertices = new vector<Vertex*>();
    mesh->faces = new vector<Face*>();
    mesh->normals = new vector<Vec3f*>();

    mesh->vertices->push_back(NULL); // 1 index
    mesh->normals->push_back(NULL); // 1 index

    string filename = tokens[1];
    mesh->label = tokens[0]; // identify that object with the label
    fstream myobjfile(filename.c_str()); // open the obj file

    /* for each line in the obj file, put into vert/face objects */
    if(myobjfile.is_open())
    {
        while(getline(myobjfile,objline))
        {
            istringstream iss2(objline);
            if (objline[0] == 'v') // if it's a vertex line
            {
                Vertex *vert = new Vertex;
                iss2 >> vf >> vert->x >> vert->y >> vert->z;
                
                mesh->vertices->push_back(vert);
            }
            else if (objline[0] == 'f') // if it's a face line
            {
                Face *face = new Face;
                iss2 >> vf >> face->idx1 >> face->idx2 >> face->idx3;

                mesh->faces->push_back(face); // add face to vector in object
            }
        }
    }
    myobjfile.close();

    return mesh;
}

/* Take line with vector parameters, store as a 
transformation matrix, return it */
Transform parse_mat(string &line)
{
    Transform t;
    istringstream iss(line);
    string trs;
    float x, y, z, theta;

    // for translation or scaling
    switch (line[0])
    {
        case 't':
            t.label = "trans";
            iss >> trs >> x >> y >> z;
            t.transformation[0] = x;
            t.transformation[1] = y;
            t.transformation[2] = z;
            break;
        case 's':
            t.label = "scal";
            iss >> trs >> x >> y >> z;
            t.transformation[0] = x;
            t.transformation[1] = y;
            t.transformation[2] = z;
            break;
        case 'r':
            t.label = "rot";
            iss >> trs >> x >> y >> z >> theta;
            t.transformation[0] = x;
            t.transformation[1] = y;
            t.transformation[2] = z;
            t.rotation_angle = theta;
            break;
    }

    return t;
}

/* Parse each object copy for its parameters and store */
void parse_objCopy(Object &object, string &line)
{
    istringstream iss(line);
    string temp;
    float r, g, b;

    if (line.find("ambient") != string::npos)
    {
        iss >> temp >> r >> g >> b;
        object.ambient[0] = r;
        object.ambient[1] = g;
        object.ambient[2] = b;
    }
    else if (line.find("diffuse") != string::npos)
    {
        iss >> temp >> r >> g >> b;
        object.diffuse[0] = r;
        object.diffuse[1] = g;
        object.diffuse[2] = b;
    }
    else if (line.find("specular") != string::npos)
    {
        iss >> temp >> r >> g >> b;
        object.specular[0] = r;
        object.specular[1] = g;
        object.specular[2] = b;
    }
    else if (line.find("shininess") != string::npos)
    {
        iss >> temp >> r;
        object.shininess = r;
    }
    else
    {
        // Take line and make into transformation 
        // and push that into object transformations
        object.transformations.push_back(parse_mat(line));
    }
}

/* parse camera section and store specific parameters */
void parse_camera(Camera &cam, vector<Light> &lights, string &line, bool &checker)
{
    istringstream iss(line);
    string temp;
    float x, y, z, theta;

    if (line.find("position") != string::npos)
    {
        iss >> temp >> x >> y >> z;
        cam.position[0] = x;
        cam.position[1] = y;
        cam.position[2] = z;
    }
    else if (line.find("orientation") != string::npos)
    {
        iss >> temp >> x >> y >> z >> theta;
        cam.orientation_angle = theta;
        cam.orientation_axis[0] = x;
        cam.orientation_axis[1] = y;
        cam.orientation_axis[2] = z;
    }
    else if (line.find("near") != string::npos)
    {
        iss >> temp >> cam.near;
    }
    else if (line.find("far") != string::npos)
    {
        iss >> temp >> cam.far;
    }
    else if (line.find("left") != string::npos)
    {
        iss >> temp >> cam.left;
    }
    else if (line.find("right") != string::npos)
    {
        iss >> temp >> cam.right;
    }
    else if (line.find("top") != string::npos)
    {
        iss >> temp >> cam.top;
    }
    else if (line.find("bottom") != string::npos)
    {
        iss >> temp >> cam.bottom;
    }
    else if (line.find("light") != string::npos)
    {
        Light light;
        iss >> temp >> light.position[0] >> light.position[1] >> light.position[2]
            >> temp >> light.color[0] >> light.color[1] >> light.color[2]
            >> temp >> light.k;
        light.position[3] = 1;
        lights.push_back(light);
    }
    else if (line.find("objects:") != string::npos)
    {
        checker = true; // incrementing counter because objects: 
        // was found so there's no longer a need to loop for c 
        // and p parameters
    }
}

void parse_all(string &line, vector<Mesh_Data*> *meshes, vector<Object> &objectsCopies,
    bool &checker, bool &checker2, Camera &camera, vector<Light> &lights)
{
    // Load variables for camera and perspective parameters 
    if (checker == false)
    {
        parse_camera(camera, lights, line, checker);
    }
    else if (line.find(".obj") != string::npos)
    {
        // parse object files and push it back into vector
        Mesh_Data *mesh = parse_obj(line);
        meshes->push_back(mesh); // add that object to object vector
    }
    else if(line.find_first_not_of(" ") != string::npos)
    {
        // parse object copies for label/parameters and push into copies vector
        checker2 = false;
        for(int i = 0; i < meshes->size(); i++)
        {
            if (line == meshes->at(i)->label)
            {
                checker2 = true;
                Object object;
                object.mesh = meshes->at(i);
                objectsCopies.push_back(object);
            }
        }
        if (checker2 == false)
        {
            // cerr << "Parsing Object Copies..." << endl;
            parse_objCopy(objectsCopies.back(), line);
        }
    }
}