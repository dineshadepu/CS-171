#include "Face.hpp"

Face::Face()
{
	v1 = 0;
	v2 = 0;
	v3 = 0;
}

Face::Face(int vert1, int vert2, int vert3)
{
	v1 = vert1;
	v2 = vert2;
	v3 = vert3;
}

Face::~Face(){}

float Face::getvOne()
{
	return v1;
}

float Face::getvTwo()
{
	return v2;
}

float Face::getvThree()
{
	return v3;
}