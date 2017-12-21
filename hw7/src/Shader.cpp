#include "Shader.hpp"

#define GL_GLEXT_PROTOTYPES 1

using namespace std;

Shader *Shader::singleton;

const char *Shader::vert_prog_file_name = "vertex.glsl";
const char *Shader::frag_prog_file_name = "fragment.glsl";

/* Creates a shader set to a rendering (per-vertex or per-pixel) mode. */
Shader::Shader(float mode) {
    this->mode = mode;
}

/* Returns/sets up the singleton instance of the class. */
Shader *Shader::getSingleton() {
    if (!Shader::singleton) {
        Shader::singleton = new Shader(default_mode);
    }

    return Shader::singleton;
}

/* Returns/sets up the singleton instance of the class. */
Shader *Shader::getSingleton(float mode) {
    if (!Shader::singleton) {
        Shader::singleton = new Shader(mode);
    }
    else {
        Shader::singleton->mode = mode;
    }

    return Shader::singleton;
}

/* Compiles the shader source files into an OpenGL program. */
void Shader::compileShaders() {
    // Read the files
    char *vert_prog = readFile((char *) vert_prog_file_name);
    if (!vert_prog)
        cerr << "Error opening vertex shader program\n";
    char *frag_prog = readFile((char *) frag_prog_file_name);
    if (!frag_prog)
        cerr << "Error opening fragment shader program\n";

    // Shader handles
    GLenum vert_shader, frag_shader;

    // Bind the first one to a vertex shader, compiled from source
    vert_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vert_shader, 1, (const GLchar **) &vert_prog, NULL);
    glCompileShader(vert_shader);

    // Check that it compiled, and if not, print the error log and exit
    GLint compiled = 0;
    glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &compiled);
    if (!compiled) {
        GLint max_length = 0;
        glGetShaderiv(vert_shader, GL_INFO_LOG_LENGTH, &max_length);
        
        // The maxLength includes the NULL character
        vector<GLchar> error_log(max_length);
        glGetShaderInfoLog(vert_shader, max_length, &max_length, &error_log[0]);
        
        // Provide the infolog in whatever manor you deem best.
        // Exit with failure.
        for (uint i = 0; i < error_log.size(); i++)
            cout << error_log[i];
        
        glDeleteShader(vert_shader); // Don't leak the shader.

        return;
    }

    // Bind the second one to a fragment shader, compiled from source
    frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(frag_shader, 1, (const GLchar **) &frag_prog, NULL);
    glCompileShader(frag_shader);

    // Repeat the error check
    compiled = 0;
    glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &compiled);
    if (!compiled) {
        GLint max_length = 0;
        glGetShaderiv(frag_shader, GL_INFO_LOG_LENGTH, &max_length);
        
        // The maxLength includes the NULL character
        vector<GLchar> error_log(max_length);
        glGetShaderInfoLog(frag_shader, max_length, &max_length, &error_log[0]);
        
        // Provide the infolog in whatever manor you deem best.
        // Exit with failure.
        for (uint i = 0; i < error_log.size(); i++)
            cout << error_log[i];
        
        glDeleteShader(frag_shader); // Don't leak the shader.

        return;
    }

    // Create an OpenGL program and attach the two shaders
    this->program = glCreateProgram();
    glAttachShader(this->program, vert_shader);
    glAttachShader(this->program, frag_shader);

    // Link the program and set it as active
    glLinkProgram(this->program);
    glUseProgram(this->program);

    // Link the current mode setting into the program 
    char mode_name[5] = "mode";
    linkf(this->mode, mode_name);
}

/*
 * Links a float into a GLSL program, assigning its value to the variable with
 * the supplied name.
 */
void Shader::linkf(float f, char *name) {
    glUniform1f(glGetUniformLocation(this->program, name), f);
}
