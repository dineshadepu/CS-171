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
#include <vector>

using namespace std;

/* Take vertex/normal arrays of current object, 
 * and based on faces, store into large
 * buffer arrays to be fed to OpenGL */
void storeBuffers(Object &object)
{
    Face face;
    for (int i = 0; i < object.faces.size(); i++)
    {
        face = object.faces[i];

        object.vertex_buffer.push_back(object.vertices[face.v1]);
        object.vertex_buffer.push_back(object.vertices[face.v2]);
        object.vertex_buffer.push_back(object.vertices[face.v3]);

        object.normal_buffer.push_back(object.normals[face.sn1]);
        object.normal_buffer.push_back(object.normals[face.sn2]);
        object.normal_buffer.push_back(object.normals[face.sn3]);
    }
}

/* Parse object file for object parameters, store them, 
and return a corresponding object */
Object parse_obj(string line)
{
    Object object; // create object for obj file
    string vf, objline, buf;
    string vn1, vn2, vn3;
    vector<string> tokens;
    int pos;
    string delim = "/"; // for the center // delimiting of faces

    istringstream iss(line);

    // using string stream, split line by whitespaces and store
    while(iss >> buf) 
        tokens.push_back(buf);

    object.vertices.push_back(Triple()); // 1 index
    object.normals.push_back(Triple()); // 1 index

    string filename = "data/" + tokens[1];
    object.label = tokens[0]; // identify that object with the label
    fstream myobjfile(filename.c_str()); // open the obj file

    /* for each line in the obj file, put into vert/face objects */
    if(myobjfile.is_open())
    {
        while(getline(myobjfile,objline))
        {
            istringstream iss2(objline);

            if (objline[0] == 'v') // if it's a vertex line
            {
                if (objline[1] != 'n') // just a vertex
                {
                    Triple vert;
                    iss2 >> vf >> vert.x >> vert.y >> vert.z; 
                    // add vertex to vertex vector
                    object.vertices.push_back(vert);
                }
                else // it's a surface normal
                {
                    Triple sn;
                    iss2 >> vf >> sn.x >> sn.y >> sn.z; 
                    // add vertex to original vertex vector
                    object.normals.push_back(sn);
                    // add vertex to vertices to be transformed later
                }
            }
            else if (objline[0] == 'f') // if it's a face line
            {
                Face face;
                vector<int> vec;
                
                while (iss2 >> vf)
                    if (vf.find("f") == string::npos)
                    {
                        // store number before and after the //
                        pos = vf.find(delim);
                        vec.push_back(atoi(vf.substr(0,pos).c_str()));
                        vec.push_back(atoi(vf.substr(pos+2, vf.length()).c_str()));
                    }

                // based on numbers in vectors, store in face fields
                face.v1 = vec[0];
                face.v2 = vec[2];
                face.v3 = vec[4];
                face.sn1 = vec[1];
                face.sn2 = vec[3];
                face.sn3 = vec[5];

                object.faces.push_back(face); // add face to vector in object
            }
        }
    }
    myobjfile.close();

    // all the faces, vertices, normals have been store
    // so go ahead and construct buffer arrays
    storeBuffers(object);

    return object;
}

/* Take line with vector parameters, store as a 
transformation matrix, return it */
Transform parse_mat(std::string line)
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
void parse_objCopy(Object &object, string line)
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
        if (line[0] != 't')
            object.Ntransformations.push_back(parse_mat(line));
    }
}

/* parse camera section and store specific parameters */
void parse_camera(Camera &cam, vector<Light> &lights, string line, bool &checker)
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

void parse_all(string &line, vector<Object> &objects, vector<Object> &objectsCopies,
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
        Object object = parse_obj(line);
        objects.push_back(object); // add that object to object vector
    }
    else if(line.find_first_not_of(" ") != string::npos)
    {
        // parse object copies for label/parameters and push into copies vector
        checker2 = false;
        for(int i = 0; i < objects.size(); i++)
        {
            if (line == objects[i].label)
            {
                checker2 = true;
                objectsCopies.push_back(objects[i]);
            }
        }

        if (checker2 == false)
            parse_objCopy(objectsCopies.back(), line);
    }
}