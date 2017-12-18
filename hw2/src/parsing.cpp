#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "structs.h"
#include "parsing.h"
#include "matrixFunctions.h"
#include <string>
#include <Eigen/Dense>

using namespace std;

/* Parse object file for object parameters, store them, 
and return a corresponding object */
Object parse_obj(string line)
{
    Object object; // create object for obj file
    string vf, objline, buf;
    string vn1, vn2, vn3;
    vector<string> tokens;
    Eigen::Matrix4d m;
    int pos;
    string delim = "/"; // for the center // delimiting of faces

    istringstream iss(line);

    // using string stream, split line by whitespaces and store
    while(iss >> buf) 
        tokens.push_back(buf);

    object.vertices.push_back(Vertex()); // 1 index
    object.v_transf.push_back(Vertex());
    object.normals.push_back(Vertex()); // 1 index
    object.n_transf.push_back(Vertex());
    m.setIdentity(4,4); // for initializing transformation matrices
    object.transfMatPoints = m; // for geom transformations
    object.transfMatNormals = m;

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
                    Vertex vert;
                    iss2 >> vf >> vert.x >> vert.y >> vert.z; 
                    // add vertex to original vertex vector
                    object.vertices.push_back(vert);
                    // add vertex to vertices to be transformed later
                    object.v_transf.push_back(vert);
                }
                else // it's a surface normal
                {
                    Vertex sn;
                    iss2 >> vf >> sn.x >> sn.y >> sn.z; 
                    // add vertex to original vertex vector
                    object.normals.push_back(sn);
                    // add vertex to vertices to be transformed later
                    object.n_transf.push_back(sn);
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
                        vec.push_back(stoi(vf.substr(0,pos)));
                        vec.push_back(stoi(vf.substr(pos+2, vf.length())));
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

    return object;
}

/* Take line with vector parameters, store as a 
transformation matrix, return it */
Eigen::Matrix4d parse_mat(std::string line)
{
    istringstream iss(line);
    Eigen::Matrix4d m;
    string trs;
    float x, y, z, theta;

    if(line[0] == 't') // translation
    {
        iss >> trs >> x >> y >> z;
        m = createTransMat(x, y, z);
    }
    else if(line[0] == 's') // scaling
    {
        iss >> trs >> x >> y >> z;
        m = createScaleMat(x, y, z);
    }
    else // rotation
    {
        iss >> trs >> x >> y >> z >> theta;
        m = createRotMat(x, y, z, theta);
    }

    return m;
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
        object.material.ambient.r = r;
        object.material.ambient.g = g;
        object.material.ambient.b = b;
    }
    else if (line.find("diffuse") != string::npos)
    {
        iss >> temp >> r >> g >> b;
        object.material.diffuse.r = r;
        object.material.diffuse.g = g;
        object.material.diffuse.b = b;
    }
    else if (line.find("specular") != string::npos)
    {
        iss >> temp >> r >> g >> b;
        object.material.specular.r = r;
        object.material.specular.g = g;
        object.material.specular.b = b;
    }
    else if (line.find("shininess") != string::npos)
    {
        iss >> temp >> r;
        object.material.shininess = r;
    }
    else
    {
        // Take vector on each line and make into transformation 
        // and push that into the matrix vector of that object 
        Eigen::Matrix4d m = parse_mat(line);
        Eigen::Matrix4d m2 = m;

        /* Combine new transformation matrix into existing one
        for that object */
        m *= object.transfMatPoints;
        object.transfMatPoints = m;
        if (line[0] != 't')
        {
            //cout << line[0] << endl;
            m2 *= object.transfMatNormals;
            object.transfMatNormals = m2;
        }
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
        cam.position = createTransMat(x, y, z);
    }
    else if (line.find("orientation") != string::npos)
    {
        iss >> temp >> x >> y >> z >> theta;
        cam.orientation = createRotMat(x, y, z, theta);
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
        iss >> temp >> light.position.x >> light.position.y >> light.position.z
            >> temp >> light.color.r >> light.color.g >> light.color.b
            >> temp >> light.k;
        lights.push_back(light);
    }
    else if (line.find("objects:") != string::npos)
    {
        checker = true; // incrementing counter because objects: 
        // was found so there's no longer a need to loop for c 
        // and p parameters
    }
}



