#ifndef SCENE_HPP
#define SCENE_HPP

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <Eigen/Eigen>

#include "command_line.hpp"
#include "model.hpp"
#include "Utilities.hpp"

using namespace std;
using namespace Eigen;

struct PointLight {
    float position[4];
    float color[3];
    float k; // Attenuation coefficient

    PointLight() = default;
    PointLight(float *position, float *color, float k);
};

// struct Primitive {
//     float e;
//     float n;
//     Matrix4f rotate;
//     Vec3f scale;
//     Vec3f translate;

//     int ures = 100;
//     int vres = 50;

//     Primitive() = default;
//     Primitive(float e, float n, float *scale, float *rotate, float theta,
//         float *translate);
//     Primitive(float e, float n, float *scale, float *rotate, float theta,
//         float *translate, int ures, int vres);

//     Vec3f getVertex(float u, float v);
//     Vec3f getNormal(Vec3f *vertex);
// };

// struct Object {
//     vector<Vec3f> vertices;
//     vector<Vec3f> normals;
//     vector<Vec3i> faces;

//     Object();
// };

static const bool default_needs_update = false;

class Scene {
    public:
        static Scene *singleton;
        static Scene *getSingleton();

        vector<Object *> root_objs;

        unordered_map<Primitive*, unsigned int> prm_tessellation_start;
        // static vector<Object> objects;
        vector<PointLight> lights;

        vector<Vector3f> vertices;
        vector<Vector3f> normals;

        Scene();

        void createLights();
        int getLightCount();

        void update();

    private:
        bool needs_update;

        void generateVertex(Primitive *prm, float u, float v);
        void tessellatePrimitive(Primitive *prm);
        void tessellateObject(Object *obj);
        // static void setupObject(Object *object);
};

#endif