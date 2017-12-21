#ifndef SHADER_HPP
#define SHADER_HPP

#include "Utilities.hpp"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

static const float default_mode = 1.0;

class Shader {
    public:
        static Shader *singleton;
        static Shader *getSingleton();
        static Shader *getSingleton(float mode);

        float mode;
        GLenum program;

        Shader(float mode);
        
        void compileShaders();
        void linkf(float f, char *name);

    private:
        static const char *vert_prog_file_name;
        static const char *frag_prog_file_name;
};

#endif