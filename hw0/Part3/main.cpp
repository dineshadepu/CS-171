#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
#include <cmath>

using namespace std;

/* Important structs */
struct Vertex
{
    float x, y, z;
};

struct Face
{
    int v1, v2, v3;
};

struct Object
{
    vector<Vertex> vertices;
    vector<Face> faces;
    Eigen::Matrix4d transfMat;
    string label;
};

/* Functions for creating matrices */
Eigen::Matrix4d createTransMat(float x, float y, float z)
{
    Eigen::Matrix4d m;
    m << 1, 0, 0, x,
         0, 1, 0, y,
         0, 0, 1, z,
         0, 0, 0, 1;
    return m;
}

Eigen::Matrix4d createRotMat(float x, float y, float z, float r)
{
    Eigen::Matrix4d m;
    m << x*x+(1-x*x)*cos(r), x*y*(1-cos(r))-z*sin(r), x*z*(1-cos(r))+y*sin(r), 0,
         y*x*(1-cos(r))+z*sin(r), y*y+(1-y*y)*cos(r), y*z*(1-cos(r))-x*sin(r), 0,
         z*x*(1-cos(r))-y*sin(r), z*y*(1-cos(r))+x*sin(r), z*z+(1-z*z)*cos(r), 0,
         0, 0, 0, 1;
    return m;
}

Eigen::Matrix4d createScaleMat(float x, float y, float z)
{
    Eigen::Matrix4d m;
    m << x, 0, 0, 0,
         0, y, 0, 0,
         0, 0, z, 0,
         0, 0, 0, 1;
    return m;
}

/* Main loop inputs file, outputs transformed objects */
int main(int argc, char* argv[])
{
    cout << endl;
	fstream myfile(argv[1]); // open the text file

    vector<Object> objects; // vector of the loaded objects
    vector<Object> objectsCopies; // vector of objects including copies

    string trs;
    float x, y, z, r;
    
    // go through each line in the file
    if(myfile.is_open())
    {
    	while(!myfile.eof())
    	{
            string line;
    		getline(myfile,line);

            // check if the line contains .obj, because if so, this is a line
            // where we need to load data from
            if (line.find(".obj") != string::npos)
            {
                Object object; // create object for obj file

                // using string stream, split line by whitespaces and store
                string buf;
                istringstream iss(line);
                vector<string> tokens;
                while(iss >> buf) 
                    tokens.push_back(buf);
                
                object.label = tokens[0]; // identify that object with the given label
                fstream myobjfile(tokens[1].c_str()); // open the obj file itself

                // initialize some variables to use later
                string objline;
                string v, f;

                object.vertices.push_back(Vertex()); // push back a null to 1 index
                m.setIdentity(4,4);
                object.transfMat = m; // initially have identity matrix
                // for each line in the obj file, put vertices/faces into object for it
                if(myobjfile.is_open())
                {
                    while(!myobjfile.eof())
                    {
                        getline(myobjfile,objline);
                        istringstream iss(objline);
                        if(objline[0] == 'v') // if it's a vertex line
                        {
                            Vertex vert;
                            iss >> v >> vert.x >> vert.y >> vert.z; // put the data into variables
                            object.vertices.push_back(vert); // add that vertex to the vertex vector in the object
                        }
                        else if(objline[0] == 'f') // if it's a face line
                        {
                            Face face;
                            iss >> f >> face.v1 >> face.v2 >> face.v3;
                            object.faces.push_back(face); // add that face to the face vector in object
                        }
                    }
                }
                objects.push_back(object); // add that object to the base object vector

                myobjfile.close();
            }
            // handle the lines that don't have .obj in it and aren't white space,
            // meaning it's either the label of the object or the transformation vector
            else if(line.find_first_not_of(" ") != string::npos)
            {
                int n = 0; // counter to see if the line is object or the matrix vectors
                for(int i = 0; i < objects.size(); i++)
                {
                    if (line == objects[i].label) // for finding which of the objects in 
                        // our object vector corresponds to this line
                    {
                        n++; // incremented because we found an object label line
                        objectsCopies.push_back(objects[i]); // for each copy push objects in
                    }
                }

                if (n == 0) // meaning that the line isn't an object label
                {
                    // take vector on each line and make into transformation matrix
                    // and push that into the matrix vector of that object
                    Eigen::Matrix4d m;
                    istringstream iss(line);
                    switch (line[0])
                    {
                        case 't':
                            // translation
                            iss >> trs >> x >> y >> z;
                            m = createTransMat(x,y,z);
                            break;
                        case 's':
                            // scaling
                            iss >> trs >> x >> y >> z;
                            m = createScaleMat(x,y,z);
                            break;
                        case 'r':
                            // rotation
                            iss >> trs >> x >> y >> z >> r;
                            m = createRotMat(x,y,z,r);
                            break;
                    }
                    Eigen::Matrix4d m2 = objectsCopies.back().transfMat;
                    m *= m2;
                    objectsCopies.back().transfMat = m;
                }
            }
    	}
    }
    myfile.close();

    /* First, relabel all the copies, then transform the points of each object */
    int n=0;
    string lbl;
    for (int i=0; i<objectsCopies.size(); i++)
    {
        if (n == 0) // so first copy
        {   
            lbl = objectsCopies[i].label;
            n++;  
        }
        else if(objectsCopies[i].label == lbl) // another copy of the previous
        {
            n++;
        }
        else // a new label that isnt a copy of the previous
        {
            n=1;
            lbl = objectsCopies[i].label;
        }
        string s;
        ostringstream convert;
        convert << n;
        s = convert.str();
        objectsCopies[i].label = lbl+"_copy"+s; // the new name that will be printed later

        // // make final transformation matrix for that object using its other matrices
        // Eigen::Matrix4d m;
        // for(int j = 0; j < objectsCopies[i].transfMats.size(); j++)
        // {
        //     if(j == 0)
        //         m = objectsCopies[i].transfMats[j];
        //     else
        //     {
        //         Eigen::Matrix4d m2 = objectsCopies[i].transfMats[j];
        //         m2 *= m;
        //         m = m2;
        //     }
        // }

        // transform the vertices of that object and re-store it in object
        for(int j = 1; j < objectsCopies[i].vertices.size(); j++)
        {
            Eigen::Vector4d vec;
            vec << objectsCopies[i].vertices[j].x, objectsCopies[i].vertices[j].y, 
                   objectsCopies[i].vertices[j].z, 1;
            Eigen::Vector4d finvec = objectsCopies[i].transfMat*vec;
            objectsCopies[i].vertices[j].x = finvec(0);
            objectsCopies[i].vertices[j].y = finvec(1);
            objectsCopies[i].vertices[j].z = finvec(2);
        }
    }

    // print out each object copy and its new points
    for (int i=0; i<objectsCopies.size(); i++)
    {
        cout << objectsCopies[i].label << endl;
        for(int j=1; j<objectsCopies[i].vertices.size(); j++)
            cout << objectsCopies[i].vertices[j].x << " " 
                 << objectsCopies[i].vertices[j].y << " "
                 << objectsCopies[i].vertices[j].z << endl;
        cout << endl;
    }
    
	return 0;
}