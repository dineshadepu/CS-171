#ifndef RENDERER_HPP
#define RENDERER_HPP

#include <thread>
#include <vector>

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "Scene.hpp"

class Camera;
class Shader;
class UI;

class Renderer {
    public:
        static Renderer *singleton;
        static Renderer *getSingleton();
        static Renderer *getSingleton(int xres, int yres);

        Renderer(int xres, int yres);

        int getXRes();
        int getYRes();
        Camera *getCamera();
        std::vector<PointLight> *getLights();

        void init();
        void start();
        
        // static void addPrimitive(float e, float n, float *scale, float *rotate,
        //     float theta, float *translate);
        // static void addObject(char *file_name);
        void updateScene();

    private:
        static void display();
        static void reshape(int xres, int yres);

        GLuint display_list;
        GLuint vb_array;
        GLuint vb_objects[2];

        Scene *scene;
        Shader *shader;
        UI *ui;

        void initLights();
        void setupLights();

        void checkUIState();
        // static uint drawObjects(uint start);
        void draw(Renderable* ren, int depth);
        void drawPrimitive(Primitive* prm);
        void drawObject(Object* obj, int depth);
        void drawAxes();
        void transform(const Transformation& trans);
};

#endif