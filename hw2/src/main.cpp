#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>
#include "matrixFunctions.h"
#include "structs.h"
#include "parsing.h"
#include "LSalgorithms.h"
#include <limits>

using namespace std;

int main(int argc, char* argv[])
{
    assert(argc == 5);

    /*** DECLARE MAJOR VARIABLES ***/
    vector<Object> objects; // vector of the loaded objects
    vector<Object> objectsCopies; // vector of objects including copies
    vector<Light> lights; // vector of lights
    vector<vector<Color>> gridCol; // screen grid of colors
    vector<vector<int>> gridBool; // screen grid of bools
    vector<vector<float>> buffer; // buffer grid

    /* Open file and store args from command line */
    int xres, yres, mode;
    fstream myfile(argv[1]);
    istringstream(argv[2]) >> xres;
    istringstream(argv[3]) >> yres;
    istringstream(argv[4]) >> mode;
    assert(xres > 0 && yres > 0 && (mode == 0 || mode == 1));

    /* Use xres and yres to make 2d vectors:
            screen grid and fill with the colors
            identical screen grid filled with zeros to serve as an indicator
            identical buffer grid */
    Color emptyColor;
    gridCol.resize(yres, vector<Color>(xres, emptyColor));
    gridBool.resize(yres, vector<int>(xres, 0));
    buffer.resize(yres, vector<float>(xres, 
        numeric_limits<float>::max()));

    bool checker = false;
    bool checker2;
    string line;
    /*** PARSING AND PUTTING DATA INTO DATA STRUCTURES ***/
    if (myfile.is_open())
    {
        while (!myfile.eof())
        {
            getline(myfile, line);
            istringstream iss(line);

            /* Load variables for camera and perspective parameters
             till we see "objects:" since that comes at the beginning
             of the file after "camera:" . Store in respective vars

             The logic is that a counter will remain 0 till we see
             "objects:", upon which the counter increments so we no longer
             check for parameters from c and p */
            if (checker == false)
            {
                /* Check each line for specification and
                load to appropriate variables */
                Camera camera;
                parse_camera(camera, lights, line, checker);
            }
            else if (line.find(".obj") != string::npos)
            {
                /* Importing lines with .obj in it, meaning importing
                objects from corresponding .obj files and storing it
                in the vector. Store objects in objects vector*/
                Object object = parse_obj(line);
                objects.push_back(object); // add that object to object vector
            }
            else if(line.find_first_not_of(" ") != string::npos)
            {
                /* Handle remaining lines that don't have .obj in it and aren't
                white space, meaning it's either the label of the objects or the
                transformation vector. Store objects and transformation
                matrices in objectsCopies vector */
                checker2 = false; // counter to see if the line 
                // is object or the matrix vectors
                for(int i = 0; i < objects.size(); i++)
                {
                    if (line == objects[i].label) // for finding which 
                        // of the objects in 
                        // our object vector corresponds to this line
                    {
                        checker2 = true; // incremented because we found an 
                        // object label line
                        objectsCopies.push_back(objects[i]); // for each 
                        // copy push objects in
                    }
                }

                if (checker2 == false) // meaning that the line isn't an 
                    // object label
                {
                    // function parses the line of that object copy
                    // for abient, diffuse, specular, or shininess,
                    // and if it doesn't find those, it instead constructs
                    // the transformation matrices using provided specs
                    // of s, t, and r
                    parse_objCopy(objectsCopies.back(), line);
                }
            }
        }
    }
    myfile.close();


    /*** CREATE TRANSFORMATION AND PROJECTION MATRICES ***/

    Eigen::Matrix4d C = camera.position*camera.orientation;
    Eigen::Matrix4d world2cam = C.inverse();

    Eigen::Matrix4d P = makePerspProj(camera.near, camera.far, camera.left, 
        camera.right, camera.top, camera.bottom); // Perspective projection matrix


    /*** TRANSFORM EACH POINT BY GEOMETRIC/NORMAL TRANSFORMATIONS ***/
    /*** AND TRANSFORM SURFACE NORMALS TO NORMAL TRANSFORMATIONS ***/

    for (int i = 0; i < objectsCopies.size(); i++)
    {
        /* First take inverse and transpose of normal transformation matrix */
        m = objectsCopies[i].transfMatNormals.inverse().transpose();
        objectsCopies[i].transfMatNormals = m;

        for (int j = 1; j < objectsCopies[i].vertices.size(); j++)
        {
            /* Transform each vertex in the object copy with the geometric transformation
            and store it in the transformed vertex vector */
            objectsCopies[i].v_transf[j] = geomTransform(objectsCopies[i].vertices[j],
                objectsCopies[i].transfMatPoints);
        }
        for (int j = 1; j < objectsCopies[i].normals.size(); j++)
        {
            /* Transform each surface normal by surface normal transf mat,
            which is inverse transpose of vertex transformation matrices
            without translations */
            objectsCopies[i].n_transf[j] = normalTransform(objectsCopies[i].normals[j],
                objectsCopies[i].transfMatNormals);
        }
    }


    /*** RENDER OBJECTS AND THEIR LIGHTS/SHADING BASED ONE MODE ***/
    Vertex v_1, v_2, v_3, sn_1, sn_2, sn_3;
    for (int i = 0; i < objectsCopies.size(); i++)
    {
        for (int j = 0; j < objectsCopies[i].faces.size(); j++)
        {
            /* Going through each object and using face indices, 
            find corresponding vertices and normals */
            v_1 = objectsCopies[i].v_transf[objectsCopies[i].faces[j].v1];
            v_2 = objectsCopies[i].v_transf[objectsCopies[i].faces[j].v2];
            v_3 = objectsCopies[i].v_transf[objectsCopies[i].faces[j].v3];
            sn_1 = objectsCopies[i].n_transf[objectsCopies[i].faces[j].sn1];
            sn_2 = objectsCopies[i].n_transf[objectsCopies[i].faces[j].sn2];
            sn_3 = objectsCopies[i].n_transf[objectsCopies[i].faces[j].sn3];

            /* Based on mode, call appropriate Shading */
            if (mode == 0)
                Gouraud_Shading(v_1, v_2, v_3, sn_1, sn_2, sn_3,
                    objectsCopies[i].material, lights, camera.position,
                    gridCol, gridBool, xres, yres, buffer, P, world2cam);
            else if (mode == 1)
                Phong_Shading(v_1, v_2, v_3, sn_1, sn_2, sn_3,
                    objectsCopies[i].material, lights, camera.position,
                    gridCol, gridBool, xres, yres, buffer, P, world2cam);
        }
    }


    /*** OUTPUT PPM ***/

    cout << "P3" << endl
    << xres << " " << yres << endl
    << "255" << endl;
    for (int i = 0; i < yres; i++)
    {
        for (int j = 0; j < xres; j++)
        {
            if (gridBool[i][j] == 1) // point is filled according to bool grid
            {
                /* Displaying from color grid */
                cout << (int)(gridCol[i][j].r*255.0) << " "
                    << (int)(gridCol[i][j].g*255.0) << " "
                    << (int)(gridCol[i][j].b*255.0) << endl;
            }
            else
            {
                cout << "0 0 0" << endl;
            }
        }
    }

    return 0;
}

