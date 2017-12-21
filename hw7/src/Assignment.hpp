#ifndef ASSIGNMENT_HPP
#define ASSIGNMENT_HPP

#include <vector>

class Camera;
class Scene;

namespace Assignment {
    void drawIOTest();
    void drawIntersectTest(Camera *camera);
    void raytrace(Camera camera, Scene scene);
};

#endif