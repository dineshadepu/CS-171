#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#define _USE_MATH_DEFINES

#include <Eigen/Eigen>

using namespace Eigen;

struct Vec3i {
    float i, j, k;

    Vec3i() = default;
    Vec3i(int i, int j, int k);
    Vec3i(int *v);
};

struct Vec3f {
    float x, y, z;

    Vec3f() = default;
    Vec3f(float x, float y, float z);
    Vec3f(float *v);
};

struct Vec4f {
    float x, y, z, w;

    Vec4f() = default;
    Vec4f(float x, float y, float z, float w);
    Vec4f(float *v);
};

float degToRad(float angle);
float radToDeg(float angle);
int sign(float x);
float pSin(float u, float p);
float pCos(float u, float p);

void makeRotateMat(float *matrix, float x, float y, float z, float angle);

char *readFile(char *file_name);

#endif