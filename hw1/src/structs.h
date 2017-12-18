#ifndef STRUCTS_H
#define STRUCTS_H

#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>

struct Vertex
{
	float x, y, z;
	int screen_x, screen_y;
};

struct Face
{
	int v1, v2, v3;
};

struct Point
{
	int x, y;
};

struct Object
{
	std::vector<Vertex> vertices;
	std::vector<Vertex> v_transf;
	std::vector<Face> faces;
	Eigen::Matrix4d transfMat;
	std::string label;
};

#endif