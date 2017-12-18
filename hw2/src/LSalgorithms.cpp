#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <cassert>
#include "matrixFunctions.h"
#include "structs.h"
#include "LSalgorithms.h"
#include <algorithm>

using namespace std;

/* Normalizing a vector */
Eigen::Vector3d normalize(Eigen::Vector3d vec)
{
    Eigen::Vector3d result;
    float normm = vec.norm();

    result = vec/normm;
    
    return result;
}

/* Color a point based on the material, the lights, and its position */
Color Lighting(Vertex point, Vertex sn, Material material, 
    vector<Light> lights, Eigen::Matrix4d position)
{
    Color c;

    Eigen::Vector4d vec;
    vec.setZero(4, 1);
    vec(3) = 1;
    Eigen::Vector3d e = (position*vec).head(3); // making a vector of the point

    Eigen::Vector3d c_d, c_a, c_s, n;
    n << sn.x, sn.y, sn.z; // representing surface normal as a vector
    // representing diffuse, ambient, and specular colors as vectors
    c_d << material.diffuse.r, material.diffuse.g, material.diffuse.b;
    c_a << material.ambient.r, material.ambient.g, material.ambient.b;
    c_s << material.specular.r, material.specular.g, material.specular.b;
    float p = material.shininess;

    Eigen::Vector3d diffuse_sum, specular_sum;
    diffuse_sum.setZero(3,1);
    specular_sum.setZero(3,1);

    Eigen::Vector3d P;
    P << point.x, point.y, point.z; // point represented as a vector

    // directional vector from camera (the eye) to point
    Eigen::Vector3d e_direction = normalize(e - P);

    /* Use lights to calculate diffuse and specular sums */
    for (int i = 0; i < lights.size(); i++)
    {
        // store the light position and color
        Eigen::Vector3d l_p l_c;
        l_p << lights[i].position.x, lights[i].position.y, 
            lights[i].position.z;
        l_c << lights[i].color.r, lights[i].color.g, lights[i].color.b;

        // make the directional vector from light to point
        Eigen::Vector3d temp = l_p - P;
        Eigen::Vector3d l_direction = normalize(temp);

        // attenuate light color
        float d = temp.norm();
        l_c *= 1/(1+lights[i].k*pow(d,2));

        // calculate diffuse light and add to sum
        Eigen::Vector3d l_diffuse = l_c*(max(0., n.dot(l_direction)));
        diffuse_sum += l_diffuse;

        // calculate specular light and add to sum
        Eigen::Vector3d l_specular = l_c*pow((max(0., 
            n.dot(normalize(e_direction + l_direction)))), p);
        specular_sum += l_specular;
    }

    Eigen::Vector3d onesvec;
    onesvec.setOnes(3,1);
    temp = onesvec.cwiseMin(c_a + diffuse_sum.cwiseProduct(c_d) 
        + specular_sum.cwiseProduct(c_s));

    // finish color assignment
    c.r = temp(0); c.g = temp(1); c.b = temp(2);

    return c;
}

/* function to assist in calculating barycentric coordinates */
int funcij(Point i, Point j, Point p)
{
    return (i.y - j.y)*p.x + (j.x - i.x)*p.y + i.x*j.y - j.x*i.y;
}

/* Three functions to compute barycentric coordinates */

float ComputeAlpha(Point a, Point b, Point c, Point p)
{
    return (float)(funcij(b, c, p))/(float)(funcij(b, c, a));
}

float ComputeBeta(Point a, Point b, Point c, Point p)
{
    return (float)(funcij(a, c, p))/(float)(funcij(a, c, b));
}

float ComputeGamma(Point a, Point b, Point c, Point p)
{
    return (float)(funcij(a, b, p))/(float)(funcij(a, b, c));
}

/* Check if NDC coordinates are in NDC projection cube */
bool inNDCCube(float alpha, float beta, float gamma, 
    Vertex a, Vertex b, Vertex c)
{
    float xdim = alpha*a.x + beta*b.x + gamma*c.x;
    float ydim = alpha*a.y + beta*b.y + gamma*c.y;
    float zdim = alpha*a.z + beta*b.z + gamma*c.z;

    if (xdim > -1.0 && xdim < 1.0
        && ydim > -1.0 && ydim < 1.0
        && zdim > -1.0 && zdim < 1.0)
        return true;
    return false;
}

/* Fill the grids with the color passed in at the given point */
void Fill(Point p, vector<vector<Color>> &gridCol, 
    vector<vector<int>> &gridBool, Color color)
{
    gridCol[p.y][p.x] = color;
    gridBool[p.y][p.x] = 1;
}

/* Take cross product of vertices to backface cull */
float BackfaceCull(Vertex a, Vertex b, Vertex c)
{
    Eigen::Vector3d va, vb, vc, cross, temp1, temp2;

    va << a.x, a.y, a.z;
    vb << b.x, b.y, b.z;
    vc << c.x, c.y, c.z;

    temp1 = vc - vb;
    temp2 = va - vb;
    cross = temp1.cross(temp2);

    return cross(2); // returning cross.z
}

/* Implemented gourard shading, which calls Lighting and RCT */
void Gouraud_Shading(Vertex v_a, Vertex v_b, Vertex v_c,
    Vertex sn_a, Vertex sn_b, Vertex sn_c,
    Material material, std::vector<Light> lights, Eigen::Matrix4d position,
    vector<vector<Color>> &gridCol, vector<vector<int>> &gridBool,
    int xres, int yres, vector<vector<float>> &buffer,
    Eigen::Matrix4d perspProj, Eigen::Matrix4d w2c)
{
    Vertex NDCa, NDCb, NDCc;
    Color colora, colorb, colorc;

    // get colors of vertices
    colora = Lighting(v_a, sn_a, material, lights, position);
    colorb = Lighting(v_b, sn_b, material, lights, position);
    colorc = Lighting(v_c, sn_c, material, lights, position);

    // convert from world space to cartesian NDC
    NDCa = worldToNDC(v_a, perspProj, w2c);
    NDCb = worldToNDC(v_b, perspProj, w2c);
    NDCc = worldToNDC(v_c, perspProj, w2c);

    // rasterize that face from those vertices with the colors
    Raster_Colored_Triangle(NDCa, NDCb, NDCc, colora, colorb, colorc,
        gridCol, gridBool, xres, yres, buffer);
}

/* Renders color/shading in a triangle using barycentric coordinates */
void Raster_Colored_Triangle(Vertex NDC_a, Vertex NDC_b, Vertex NDC_c, 
    Color color_a, Color color_b, Color color_c, 
    vector<vector<Color>> &gridCol, vector<vector<int>> &gridBool,
    int xres, int yres, vector<vector<float>> &buffer)
{
    /* Do Backface culling */
    float crossz = BackfaceCull(NDC_a, NDC_b, NDC_c);
    if (crossz >= 0)
    {
        Color color;

        // convert from ndc coordinates to screen coords
        Point a = NDCToScreen(NDC_a, xres, yres);
        Point b = NDCToScreen(NDC_b, xres, yres);
        Point c = NDCToScreen(NDC_c, xres, yres);

        int xmin = min(min(a.x, b.x), c.x);
        int ymin = min(min(a.y, b.y), c.y);
        int xmax = max(max(a.x, b.x), c.x);
        int ymax = max(max(a.y, b.y), c.y);

        // pass over all the x and y coordinates composing the triangle
        // made by the points, and fill colors
        for (int x = xmin; x <= xmax; x++)
        {
            for (int y = ymin; y <= ymax; y++)
            {
                // make point vertex from current x and y
                p.x = x; p.y = y;

                float alpha = ComputeAlpha(a, b, c, p);
                float beta = ComputeBeta(a, b, c, p);
                float gamma = ComputeGamma(a, b, c, p);

                // check if alphas are in range and points are in NDC Cube
                if (alpha >= 0.0 && alpha <= 1.0
                    && beta >= 0.0 && beta <= 1.0
                    && gamma >= 0.0 && gamma <= 1.0
                    && inNDCCube(alpha, beta, gamma, NDC_a, NDC_b, NDC_c))
                {
                    // Check z coordinate against buffer to make sure
                    // it should be rendered
                    float NDCz = alpha*NDC_a.z + beta*NDC_b.z + gamma*NDC_c.z;
                    if (!(NDCz > buffer[p.y][p.x]))
                    {
                        buffer[p.y][p.x] = NDCz;
                        // calculate color using barycentric coords
                        color.r = alpha*color_a.r + 
                            beta*color_b.r + gamma*color_c.r;
                        color.g = alpha*color_a.g + 
                            beta*color_b.g + gamma*color_c.g;
                        color.b = alpha*color_a.b + 
                            beta*color_b.b + gamma*color_c.b;
                        Fill(p, gridCol, gridBool, color); // fill grid with color
                    }
                }
            }
        }
    }   
}

/* Implemented phong shading, which colors by pixels instead of gradient */
void Phong_Shading(Vertex v_a, Vertex v_b, Vertex v_c,
    Vertex sn_a, Vertex sn_b, Vertex sn_c,
    Material material, std::vector<Light> lights, Eigen::Matrix4d position,
    vector<vector<Color>> &gridCol, vector<vector<int>> &gridBool,
    int xres, int yres, vector<vector<float>> &buffer,
    Eigen::Matrix4d perspProj, Eigen::Matrix4d w2c)
{
    // convert vertices from world to NDC coords to backface cull
    Vertex NDCa = worldToNDC(v_a, perspProj, w2c);
    Vertex NDCb = worldToNDC(v_b, perspProj, w2c);
    Vertex NDCc = worldToNDC(v_c, perspProj, w2c);

    float crossz = BackfaceCull(NDCa, NDCb, NDCc);
    if (crossz >= 0)
    {
        Color color;

        // convert to screen coordinates
        Point a = NDCToScreen(NDCa, xres, yres);
        Point b = NDCToScreen(NDCb, xres, yres);
        Point c = NDCToScreen(NDCc, xres, yres);

        int xmin = min(min(a.x, b.x), c.x);
        int ymin = min(min(a.y, b.y), c.y);
        int xmax = max(max(a.x, b.x), c.x);
        int ymax = max(max(a.y, b.y), c.y);

        float alpha, beta, gamma, NDCz;
        Vertex sn, v;
        for (int x = xmin; x <= xmax; x++)
        {
            for (int y = ymin; y <= ymax; y++)
            {
                // create point vertex of current x and y
                p.x = x; p.y = y;

                // calculate barycentric coordinates
                float alpha = ComputeAlpha(a, b, c, p);
                float beta = ComputeBeta(a, b, c, p);
                float gamma = ComputeGamma(a, b, c, p);

                // check if alphas are in range and points are in NDC Cube
                if (alpha >= 0.0 && alpha <= 1.0
                    && beta >= 0.0 && beta <= 1.0
                    && gamma >= 0.0 && gamma <= 1.0
                    && inNDCCube(alpha, beta, gamma, NDCa, NDCb, NDCc))
                {
                    // Check z coordinate against buffer to make sure
                    // it should be rendered
                    float NDCz = alpha*NDCa.z + beta*NDCb.z + gamma*NDCc.z;
                    if (!(NDCz > buffer[p.y][p.x]))
                    {
                        buffer[p.y][p.x] = NDCz;

                        // calculate current point/normal
                        // from barycentric coordinates
                        v.x = alpha*v_a.x + beta*v_b.x + gamma*v_c.x;
                        v.y = alpha*v_a.y + beta*v_b.y + gamma*v_c.y;
                        v.z = alpha*v_a.z + beta*v_b.z + gamma*v_c.z;

                        sn.x = alpha*sn_a.x + beta*sn_b.x + gamma*sn_c.x;
                        sn.y = alpha*sn_a.y + beta*sn_b.y + gamma*sn_c.y;
                        sn.z = alpha*sn_a.z + beta*sn_b.z + gamma*sn_c.z;
                        
                        // calculate color at that point given the lighting
                        color = Lighting(v, sn, material, lights, position);

                        // fill the grid
                        Fill(p, gridCol, gridBool, color);
                    }
                }
            }
        }
    }
}


