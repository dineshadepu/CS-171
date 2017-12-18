#ifndef STRUCTS_H
#define STRUCTS_H

#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>

struct Color
{
    float r, g, b;
};

struct Point
{
    int x, y;
};

struct Vertex
{
    float x, y, z;
};

struct Face
{
    int v1, v2, v3;
    int sn1, sn2, sn3;
};

struct Camera
{
    Eigen::Matrix4d position, orientation;
    float near, far, left, right, top, bottom;
};

struct Light
{
    Vertex position;
    Color color;
    float k;
};

struct Material
{
    Color ambient, diffuse, specular;
    float shininess;
};

struct Object
{
    std::string label;

    std::vector<Vertex> vertices;
    std::vector<Vertex> v_transf;
    std::vector<Vertex> normals;
    std::vector<Vertex> n_transf;
    std::vector<Face> faces;

    Material material;

    Eigen::Matrix4d transfMatPoints;
    Eigen::Matrix4d transfMatNormals;
};

#endif