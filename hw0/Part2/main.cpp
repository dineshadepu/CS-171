/*************************************************************************
main.cpp

Inputs a single txt file representing transformations and outputting
the overall transformation matrix product inverted using Eigen.

Author: Eshan Govil
Assignment: Hw 0, part 2
*************************************************************************/

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

int main(int argc, char* argv[])
{
	fstream myfile(argv[1]); // open the single inputted txt file

    // declare variables to store data later
	string line;
    string trs;
    float x, y, z, r;
   	vector<Eigen::Matrix4d> mats; // vector of all matrices

    // for each line in the file, create transformation matrix and store it
    if(myfile.is_open())
    {
    	while(!myfile.eof())
    	{
    		Eigen::Matrix4d m;
    		getline(myfile,line);
    		istringstream iss(line);
    		if(line[0] == 't') // translation vector
    		{
    			iss >> trs >> x >> y >> z;
    			m << 1, 0, 0, x,
	    			 0, 1, 0, y,
	    			 0, 0, 1, z,
	    			 0, 0, 0, 1;
    			mats.push_back(m); // push back the vector-made-matrix
    		}
    		else if(line[0] == 's') // scaling vector
    		{
    			iss >> trs >> x >> y >> z;
    			m << x, 0, 0, 0,
	    			 0, y, 0, 0,
	    			 0, 0, z, 0,
	    			 0, 0, 0, 1;
    			mats.push_back(m);
    		}
    		else if(line[0] == 'r') // rotation vector
    		{
    			iss >> trs >> x >> y >> z >> r;
    			m << x*x+(1-x*x)*cos(r), x*y*(1-cos(r))-z*sin(r), 
                     x*z*(1-cos(r))+y*sin(r), 0,
	    		     y*x*(1-cos(r))+z*sin(r), y*y+(1-y*y)*cos(r), 
                     y*z*(1-cos(r))-x*sin(r), 0,
	    			 z*x*(1-cos(r))-y*sin(r), z*y*(1-cos(r))+x*sin(r), 
                     z*z+(1-z*z)*cos(r), 0,
	    			 0, 0, 0, 1;
    			mats.push_back(m);
    		}
    	}
    }
    myfile.close();

    // go through vector of matrices and create a product in correct order
    Eigen::Matrix4d m = mats[0];
    for(int i = 1; i < mats.size(); i++)
    {
		Eigen::Matrix4d m2 = mats[i];
		m2 *= m;
		m = m2;
    }

    Eigen::Matrix4d m_inv = m.inverse(); // invert the matrix in the final step

    cout << m_inv << endl; // print the inversed transformation matrix

	return 0;
}
