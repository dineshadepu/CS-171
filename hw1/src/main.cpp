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
#include "bresenham.h"
#include "parsing.h"

using namespace std;

int main(int argc, char* argv[])
{
    assert(argc == 4);

    vector<Object> objects; // vector of the loaded objects
    vector<Object> objectsCopies; // vector of objects including copies
    vector<Point> points;

    /* DECLARE IMPORTANT VARIABLES for parsing */
    float n, f, l, r, t, b, x, y, z, theta;
    int checker = false; int checker2;
    Eigen::Matrix4d position, orientation, m, m2;
    string line, trs, temp, objline;
    
    /* Open file and store from command line */
    fstream myfile(argv[1]);
    int xres, yres;
    istringstream(argv[2]) >> xres;
    istringstream(argv[3]) >> yres;
    assert(xres > 0 && yres > 0);

    /* Use xres and yres to make 2d vector for screen grid and fill with 0's*/
    vector<vector<int>> grid;
    grid.resize(yres, vector<int>(xres, 0)); // 0's mean screen is blank


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
                load to appropriate variables 

                Need to do this in main function and not in separate
                function because of the hard coded nature of the 
                variables we need to extract from the file */
                if (line.find("position") != string::npos)
                {
                    iss >> temp >> x >> y >> z;
                    position = createTransMat(x, y, z);
                }
                else if (line.find("orientation") != string::npos)
                {
                    iss >> temp >> x >> y >> z >> theta;
                    orientation = createRotMat(x, y, z, theta);
                }
                else if (line.find("near") != string::npos)
                {
                    iss >> temp >> n;
                }
                else if (line.find("far") != string::npos)
                {
                    iss >> temp >> f;
                }
                else if (line.find("left") != string::npos)
                {
                    iss >> temp >> l;
                }
                else if (line.find("right") != string::npos)
                {
                    iss >> temp >> r;
                }
                else if (line.find("top") != string::npos)
                {
                    iss >> temp >> t;
                }
                else if (line.find("bottom") != string::npos)
                {
                    iss >> temp >> b;
                }
                else if (line.find("objects:") != string::npos)
                {
                    checker = true; // incrementing counter because objects: 
                    // was found so there's no longer a need to loop for c 
                    // and p parameters
                }
            }
            /* Importing lines with .obj in it, meaning importing
            objects from corresponding .obj files and storing it
            in the vector. Store objects in objects vector*/
            else if (line.find(".obj") != string::npos)
            {
                Object object = parse_obj(line);
                objects.push_back(object); // add that object to object vector
                
            }
            /* Handle remaining lines that don't have .obj in it and aren't
            white space, meaning it's either the label of the objects or the
            transformation vector. Store objects and transformation
            matrices in objectsCopies vector */
            else if(line.find_first_not_of(" ") != string::npos)
            {
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
                    /* Take vector on each line and make into transformation 
                    and push that into the matrix vector of that object */
                    m = parse_mat(line);
                    /* Combine new transformation matrix into existing one
                    for that object */
                    m *= objectsCopies.back().transfMat;
                    objectsCopies.back().transfMat = m;
                }
            }
        }
    }
    myfile.close();


    /*** CREATE TRANSFORMATION AND PROJECTION MATRICES ***/
    Eigen::Matrix4d C = position*orientation;
    Eigen::Matrix4d world2cam = C.inverse();
    Eigen::Matrix4d P = makePerspProj(n, f, l, r, t, b); 


    /*** TRANSFORM EACH POINT AND MAP TO SCREEN COORDS***/
    Eigen::Vector4d vec, transfPoint;
    for (int i = 0; i < objectsCopies.size(); i++)
    {
        for (int j = 1; j < objectsCopies[i].vertices.size(); j++)
        {   
            // make vector of point
            vec << objectsCopies[i].vertices[j].x, 
            objectsCopies[i].vertices[j].y, 
            objectsCopies[i].vertices[j].z, 1;

            // First, geometric transformation in world space
            // then transformation from world to camera space
            // to have [x_c  y_c  z_c  1]
            // then transform to homogeneous NDC using perspective projection 
            // which is [x_ndc  y_ndc  z_ndc  w_ndc] where w_ndc = -z_c
            transfPoint = P*world2cam*objectsCopies[i].transfMat*vec;

            // Store Cartesian NDC by scaling by -z_c
            // in separate vertex vector of transformed points
            // that can always be overwritten later.
            // DOES NOT CHANGE ORIGINAL UNTRANSFORMED VERTICES
            objectsCopies[i].v_transf[j].x = transfPoint(0)/transfPoint(3);
            objectsCopies[i].v_transf[j].y = transfPoint(1)/transfPoint(3);
            objectsCopies[i].v_transf[j].z = transfPoint(2)/transfPoint(3);

            // Map NDC points to screen coordinates
            // Since the box is in the range [-1,1] in both x and y,
            // points will be in this range. We can use this range to normalize
            // points against xres and yres.
            // We scale to the pixel gride dimensions while keeping 0,0 center
            // then shifting so that 0,0 is corner and using int to round
            objectsCopies[i].v_transf[j].screen_x = 
            int(objectsCopies[i].v_transf[j].x*(xres/2)+(xres/2));
            objectsCopies[i].v_transf[j].screen_y = 
            int(objectsCopies[i].v_transf[j].y*(-yres/2)+(yres/2));
        }
    }


    /*** FILL GRID USING FACES FROM OBJECTS ***/
    Vertex v_1, v_2, v_3;
    for (int i = 0; i < objectsCopies.size(); i++)
    {
        for (int j = 0; j < objectsCopies[i].faces.size(); j++)
        {
            // from each face, extract the vertex indices
            // to find the corresponding vertices in that object
            v_1 = objectsCopies[i].v_transf[objectsCopies[i].faces[j].v1];
            v_2 = objectsCopies[i].v_transf[objectsCopies[i].faces[j].v2];
            v_3 = objectsCopies[i].v_transf[objectsCopies[i].faces[j].v3];

            // Check to make sure vertex isn't out of cube bounds
            // or it won't render upon making the line
            if (v_1.x < 1 && v_1.x > -1
                && v_1.y < 1 && v_1.y > -1
                && v_2.x < 1 && v_2.x > -1
                && v_2.y < 1 && v_2.y > -1
                && v_3.x < 1 && v_3.x > -1
                && v_3.y < 1 && v_3.y > -1)
            {
                fillGrid(v_1.screen_x, v_1.screen_y, 
                    v_2.screen_x, v_2.screen_y, 
                    grid);
                fillGrid(v_2.screen_x, v_2.screen_y, 
                    v_3.screen_x, v_3.screen_y, 
                    grid);
                fillGrid(v_3.screen_x, v_3.screen_y, 
                    v_1.screen_x, v_1.screen_y, 
                    grid);
            }
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
            if (grid[i][j] == 1) // point is filled
            {
                cout << "255 255 255" << endl;
            }
            else
            {
                cout << "0 0 0" << endl;
            }
        }
    }

    return 0;
}