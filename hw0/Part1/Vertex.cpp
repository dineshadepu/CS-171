#include "Vertex.hpp"

Vertex::Vertex()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

Vertex::Vertex(float x1, float x2, float x3)
{
	x = x1;
	y = x2;
	z = x3;
}

Vertex::~Vertex(){}

float Vertex::getX()
{
	return x;
}

float Vertex::getY()
{
	return y;
}

float Vertex::getZ()
{
	return z;
}