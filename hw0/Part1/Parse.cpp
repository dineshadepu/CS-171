/*************************************************************************
Parse.cpp

Parses OBJ files and creates a shape object containing a vector of its
vertices and faces.

Author: Eshan Govil
Assignment: Hw 0, part 1
*************************************************************************/

#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "Vertex.hpp" // classes made to represent vertex/face objects
#include "Face.hpp"

using namespace std;

/* Initially made classes for both. Structs will work more easily,
so I'll use structs going ahead in Part 3 but keep classes for now.

struct Vertex
{
	float x, y, z;
};

struct Face
{
	int v1, v2, v3;
};
*/

/* Struct for object of vectors */
struct Object
{
	vector<Vertex> vertices;
	vector<Face> faces;
};

int main(int argc, char* argv[]) {
	// creating vector of objects for all command line arg .obj files
	vector<Object> objects; 

	/* Parse all object files and store in data structures */
	for (int i = 1; i < argc; i++)
	{
	    /* Open obj file and save lines to structs */
		
	    fstream myfile(argv[i]); // opens an object file from the command line args

	    // declare variables for lines being read, data on each line
	    // and the vectors holding the items
	    string line;
	    string v, f;
	    float x, y, z;
	    int v1, v2, v3;

	    Object object; // Initializing an object
	    object.vertices.push_back(Vertex()); // push back a null to 1 index vertices

	    // for each line in the file, put variables into appropriate struct
	    if(myfile.is_open())
	    {
	    	while(!myfile.eof())
	    	{
	    		getline(myfile,line);
	    		istringstream iss(line);
	    		switch (line[0])
	    		{
	    			case 'v':
	    				// vertex info, store floats
		    			iss >> v >> x >> y >> z; // stream the data into variables
		    			Vertex vert = Vertex(x, y, z);
		    			object.vertices.push_back(vert);
		    			break;
		    		case 'f':
		    			// face info, store ints
			    		iss >> f >> v1 >> v2 >> v3;
		    			Face face = Face(v1, v2, v3);
		    			object.faces.push_back(face);
		    			break;

	    		}
	    	}
	    }
	    myfile.close();
	    objects.push_back(object);
	}

	/* Print the vertex/face contents of stored data structures */
	for(int i = 0; i < objects.size(); i++)
	{
		// First, print file name of file being read
		string obj_file = argv[i+1];
		istringstream iss(obj_file);
		string first;
		getline(iss, first, '.'); // splitting tokens by period
		cout << first << ":\n" << endl;

		// print from the vectors within the object struct, first vector then face
	    for(int j = 1; j < objects[i].vertices.size(); j++)
	    {
	    	cout << "v " << objects[i].vertices[j].getX() << " " 
	    		 << objects[i].vertices[j].getY() << " " 
	    		 << objects[i].vertices[j].getZ() << endl;
	    }
	    for(int j = 0; j < objects[i].faces.size(); j++)
	    {
	    	cout << "f " << objects[i].faces[j].getvOne() << " " 
	    		 << objects[i].faces[j].getvTwo() << " " 
	    		 << objects[i].faces[j].getvThree() << endl;
	    }

	    cout << endl;
	}

	
    return 0;
}