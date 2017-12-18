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

Object parse_obj(string line)
{
	Object object; // create object for obj file
	string vf, objline, buf;
	vector<string> tokens;
	Eigen::Matrix4d m;

	istringstream iss(line);

	// using string stream, split line by whitespaces and store
	while(iss >> buf) 
		tokens.push_back(buf);

	object.vertices.push_back(Vertex()); // 1 index
	object.v_transf.push_back(Vertex());
	m.setIdentity(4,4);
	object.transfMat = m; // initially have identity matrix
	// for combining transformation matrices to later

	object.label = tokens[0]; // identify that object with the label
	fstream myobjfile(tokens[1].c_str()); // open the obj file 

	/* for each line in the obj file, put into vert/face objects */
	if(myobjfile.is_open())
	{
		while(getline(myobjfile,objline))
		{
			istringstream iss2(objline);

			if(objline[0] == 'v') // if it's a vertex line
			{
				Vertex vert;
				// put the data into variables
				iss2 >> vf >> vert.x >> vert.y >> vert.z; 
				// add vertex to original vertex vector
				object.vertices.push_back(vert);
				// add vertex to vertices to be transformed later
				object.v_transf.push_back(vert);
			}
			else if(objline[0] == 'f') // if it's a face line
			{
				Face face;
				iss2 >> vf >> face.v1 >> face.v2 >> face.v3;
				object.faces.push_back(face); // add face to vector in object
			}
		}
	}
	myobjfile.close();

	return object;
}

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