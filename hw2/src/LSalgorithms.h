#include <iostream>
#include <cstdio>
#include <string>
#include "structs.h"
#include <Eigen/Dense>

using namespace std;

/* Normalizing a vector */
Eigen::Vector3d normalize(Eigen::Vector3d vec);

/* Color a point based on the material, the lights, and its position */
Color Lighting(Vertex P, Vertex sn, Material material, 
    vector<Light> lights, Eigen::Matrix4d position);

/* function to assist in calculating barycentric coordinates */
int funcij(Point i, Point j, Point p);

float ComputeAlpha(Point a, Point b, Point c, Point p);
float ComputeBeta(Point a, Point b, Point c, Point p);
float ComputeGamma(Point a, Point b, Point c, Point p);

/* Check if NDC coordinates are in NDC projection cube */
bool inNDCCube(float alpha, float beta, float gamma, 
    Vertex a, Vertex b, Vertex c);

/* Fill the grids with the color passed in at the given point */
void Fill(Point p, vector<vector<Color>> &gridCol, 
    vector<vector<int>> &gridBool, Color color);

/* Take cross product of vertices to backface cull */
float BackfaceCull(Vertex a, Vertex b, Vertex c);

/* Implemented gourard shading, which calls Lighting and RCT */
void Gouraud_Shading(Vertex v_a, Vertex v_b, Vertex v_c,
    Vertex sn_a, Vertex sn_b, Vertex sn_c,
    Material material, vector<Light> lights, Eigen::Matrix4d position,
    vector<vector<Color>> &gridCol, vector<vector<int>> &gridBool,
    int xres, int yres, vector<vector<float>> &buffer,
    Eigen::Matrix4d perspProj, Eigen::Matrix4d w2c);

/* Renders color/shading in a triangle using barycentric coordinates */
void Raster_Colored_Triangle(Vertex NDC_a, Vertex NDC_b, Vertex NDC_c, 
    Color color_a, Color color_b, Color color_c, 
    vector<vector<Color>> &gridCol, vector<vector<int>> &gridBool,
    int xres, int yres, vector<vector<float>> &buffer);

void Phong_Shading(Vertex v_a, Vertex v_b, Vertex v_c,
    Vertex sn_a, Vertex sn_b, Vertex sn_c,
    Material material, vector<Light> lights, Eigen::Matrix4d position,
    vector<vector<Color>> &gridCol, vector<vector<int>> &gridBool,
    int xres, int yres, vector<vector<float>> &buffer,
    Eigen::Matrix4d perspProj, Eigen::Matrix4d w2c);
