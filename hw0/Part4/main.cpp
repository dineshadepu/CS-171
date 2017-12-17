/*************************************************************************
main.cpp

Outputs an white ricle image in ppm format through the std output
given input image resolution.

Author: Eshan Govil
Assignment: Hw 0, part 4
*************************************************************************/

#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{
	int xres;
	int yres;

	// store the command line xres/yres into int variables
	istringstream(argv[1]) >> xres; 
	istringstream(argv[2]) >> yres;

	double r = min(xres,yres)/4; // the radius

	cout << "P3" << endl
		 << xres << " " << yres << endl
		 << "255" << endl;

	for (int i=0; i<yres; i++)
	{
		for (int j=0; j<xres; j++)
		{
			// calculating relative to the center of the circle
			if((j-xres/2)*(j-xres/2)+(i-yres/2)*(i-yres/2) <= r*r)
				cout << "255 0 255" << endl;
			else
				cout << "0 255 0" << endl;
		}
	}

	return 0;
}