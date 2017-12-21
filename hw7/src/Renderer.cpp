#include "Renderer.hpp"

#include "Utilities.hpp"

#include "Shader.hpp"
#include "UI.hpp"
#include "Assignment.hpp"

#define MAX_RECURSION_DEPTH 25

Renderer *Renderer::singleton;

/*
 * Creates the data structure used to store the scene, the shader program
 * manager, and the ui engine.
 */
Renderer::Renderer(int xres, int yres) {
    this->scene = Scene::getSingleton();
    this->shader = Shader::getSingleton();
    this->ui = UI::getSingleton(xres, yres);
}

/* Returns/sets up the singleton instance of the class. */
Renderer *Renderer::getSingleton() {
    if (!Renderer::singleton) {
        Renderer::singleton = new Renderer(default_xres, default_yres);
    }

    return Renderer::singleton;
}

/* Returns/sets up the singleton instance of the class. */
Renderer *Renderer::getSingleton(int xres, int yres) {
    if (!Renderer::singleton) {
        Renderer::singleton = new Renderer(xres, yres);
    }
    else {
        Renderer::reshape(xres, yres);
    }

    return Renderer::singleton;
}

int Renderer::getXRes() {
    return this->ui->xres;
}

int Renderer::getYRes() {
    return this->ui->yres;
}

Camera *Renderer::getCamera() {
    return &this->ui->camera;
}

vector<PointLight> *Renderer::getLights() {
    return &this->scene->lights;
}

void Renderer::updateScene() {
    this->scene->update();
}

/* Initializes the renderer. */
void Renderer::init() {
    // I (Nailen) needed this to get shaders to work
    glewExperimental = GL_TRUE; 
    glewInit();

    // Set the background color and enable the depth buffer, backface culling
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Set up two vertex attribute buffers. This is important because it allows
    // us to copy all of the vertices and normals in the scene to faster memory
    // in the graphics hardware (GPU) before rendering, as opposed to reading
    // each value sequentially from system memory (much slower).
    glGenVertexArrays(1, &this->vb_array);
    glBindVertexArray(this->vb_array);
    glGenBuffers(2, this->vb_objects);

    // Bind the v_v and n_v variables in vertex.glsl to the two buffers
    glBindAttribLocation(this->shader->program, 0, "v_v");
    glBindAttribLocation(this->shader->program, 1, "n_v");

    // Set up the perspective projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(this->ui->camera.fov, this->ui->camera.aspect,
        this->ui->camera.near, this->ui->camera.far);

    // Set up the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-this->ui->camera.position.x, -this->ui->camera.position.y,
        -this->ui->camera.position.z);
    glRotatef(radToDeg(-this->ui->camera.angle), this->ui->camera.axis.x,
        this->ui->camera.axis.y, this->ui->camera.axis.z);

    // Initialize the lights, set the material properties for every polygon in
    // the scene, and then compile and activate the shader program
    initLights();
    this->shader->compileShaders();
}

/* Initializes the OpenGL light parameters from the scene-> */
void Renderer::initLights() {
    // Enable built-in lights
    glEnable(GL_LIGHTING);
    
    // Loop through all the lights and bind their properties to OpenGL's lights
    for (uint i = 0; i < this->scene->lights.size(); i++) {
        GLint light_id = GL_LIGHT0 + i;
        glEnable(light_id);
        
        glLightfv(light_id, GL_AMBIENT, this->scene->lights[i].color);
        glLightfv(light_id, GL_DIFFUSE, this->scene->lights[i].color);
        glLightfv(light_id, GL_SPECULAR, this->scene->lights[i].color);

        glLightf(light_id, GL_QUADRATIC_ATTENUATION, this->scene->lights[i].k);
    }
}

/* Places the OpenGL lights in position for the scene-> */
void Renderer::setupLights() {
    // Loop through all the lights and bind their positions to OpenGL's lights
    for (uint i = 0; i < this->scene->lights.size(); i++)
        glLightfv(GL_LIGHT0 + i, GL_POSITION, this->scene->lights[i].position);
}

/* Starts the main OpenGL loop. */
void Renderer::start() {
    // Functions are pretty self-explanatory, but don't be afraid to to use
    // www.opengl.org/sdk/docs/ and www.opengl.org/wiki/
    glutReshapeFunc(Renderer::reshape);
    glutDisplayFunc(Renderer::display);
    glutMouseFunc(UI::handleMouseButton);
    glutMotionFunc(UI::handleMouseMotion);
    glutKeyboardFunc(UI::handleKeyPress);
    glutMainLoop();
}

/* Checks for and applies pending this->ui updates. */
void Renderer::checkUIState() {
    // If the shader mode has changed, re-link the mode variable
    if (this->shader->mode != this->ui->shader_mode) {
        this->shader->mode = this->ui->shader_mode;
        char mode_name[5] = "mode";
        this->shader->linkf(this->shader->mode, mode_name);
    }
    // If we need to rebuild the vertex buffers, do so
    if (this->ui->rebuild_scene) {
        this->ui->rebuild_scene = false;
        this->scene->update();
    }
    // If we need to raytrace the scene, do so
    if (this->ui->raytrace_scene) {
        this->ui->raytrace_scene = false;
        // Spawn the raytracer thread and detach it to keep doing OpenGL stuff
        printf("Raytracing operation started on a new thread...\n");
        thread(Assignment::raytrace, this->ui->camera, Scene(*this->scene))
            .detach();
    }
}


void Renderer::draw(Renderable* ren, int depth) {
    assert(ren);

    if (ren->getType() == PRM) {
        drawPrimitive(dynamic_cast<Primitive*>(ren));
    } else if (ren->getType() == OBJ) {
        if (depth <= MAX_RECURSION_DEPTH)
            drawObject(dynamic_cast<Object*>(ren), depth);
    } else {
        fprintf(stderr, "Renderer::draw ERROR invalid RenderableType %d\n",
            ren->getType());
        exit(1);
    }
}

/* FIX HERE
 * Draws the primitives present in the scene, assuming their vertices and
 * normals start at the given offset in the respective vertex buffers.
 */
void Renderer::drawPrimitive(Primitive *prm) {
    assert(this->scene->prm_tessellation_start.find(prm) !=
        this->scene->prm_tessellation_start.end());

    const RGBf& color = prm->getColor();
    const float ambient = prm->getAmbient();
    const float diffuse = prm->getDiffuse();
    const float specular = prm->getSpecular();
    float ambientColor[3] = {color.r * ambient, 
                        color.g * ambient, 
                        color.b * ambient};
    float diffuseColor[3] = {color.r * diffuse, 
                        color.g * diffuse, 
                        color.b * diffuse};
    float specularColor[3] = {color.r * specular, 
                        color.g * specular, 
                        color.b * specular};

    // Set the built-in material properties
    glMaterialfv(GL_FRONT, GL_AMBIENT, ambientColor);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuseColor);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specularColor);
    glMaterialf(GL_FRONT, GL_SHININESS, prm->getReflected());

    uint start = this->scene->prm_tessellation_start[prm];

    int ures = prm->getPatchX();
    int vres = prm->getPatchY();

    uint offset;
    Vector3f *vertex, *normal, endpoint;
    // Draw all vres - 2 of its circumferential triangle strips
    offset = 2 * (ures + 1);
    for (int j = 1; j < vres - 1; j++) {
        if (!this->ui->wireframe_mode)
            glDrawArrays(GL_TRIANGLE_STRIP, start, offset);
        // If we're in wireframe mode, we have to manually specify the loops
        else {
            for (uint i = start; i < start + offset - 3; i += 2) {
                glBegin(GL_LINE_LOOP);
                vertex = &this->scene->vertices.at(i);
                normal = &this->scene->normals.at(i);
                glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);

                vertex = &this->scene->vertices.at(i + 1);
                normal = &this->scene->normals.at(i + 1);
                glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);

                vertex = &this->scene->vertices.at(i + 2);
                normal = &this->scene->normals.at(i + 2);
                glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);
                glEnd();

                glBegin(GL_LINE_LOOP);
                vertex = &this->scene->vertices.at(i + 1);
                normal = &this->scene->normals.at(i + 1);
                glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);

                vertex = &this->scene->vertices.at(i + 3);
                normal = &this->scene->normals.at(i + 3);
                glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);

                vertex = &this->scene->vertices.at(i + 2);
                normal = &this->scene->normals.at(i + 2);
                glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);
                glEnd();

                // If we're in normal mode, draw them sticking out of the surface
                if (this->ui->normal_mode) {
                    glBegin(GL_LINES);
                    vertex = &this->scene->vertices.at(i);
                    normal = &this->scene->normals.at(i);
                    glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                    endpoint = *vertex + *normal * 0.25;
                    glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
                    glEnd();

                    glBegin(GL_LINES);
                    vertex = &this->scene->vertices.at(i + 1);
                    normal = &this->scene->normals.at(i + 1);
                    glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                    endpoint = *vertex + *normal * 0.25;
                    glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
                    glEnd();
                }
            }
        }
        start += offset;
    }

    // Draw the triangle fans that meet at the superquadric's poles
    offset = ures + 2;
    if (!this->ui->wireframe_mode)
        glDrawArrays(GL_TRIANGLE_FAN, start, offset);
    // If we're in wireframe mode, we have to manually specify the loops
    else {
        for (uint i = start; i < start + offset - 1; i++) {
            glBegin(GL_LINE_LOOP);
            vertex = &this->scene->vertices.at(start);
            normal = &this->scene->normals.at(start);
            glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
            glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);

            vertex = &this->scene->vertices.at(i);
            normal = &this->scene->normals.at(i);
            glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
            glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);

            vertex = &this->scene->vertices.at(i + 1);
            normal = &this->scene->normals.at(i + 1);
            glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
            glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);
            glEnd();

            // If we're in normal mode, draw it sticking out of the surface
            if (this->ui->normal_mode) {
                glBegin(GL_LINES);
                vertex = &this->scene->vertices.at(i);
                normal = &this->scene->normals.at(i);
                glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                endpoint = *vertex + *normal * 0.25;
                glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
                glEnd();
            }
        }
    }
    start += offset;
    if (!this->ui->wireframe_mode)
        glDrawArrays(GL_TRIANGLE_FAN, start, offset);
    // If we're in wireframe mode, we have to manually specify the loops
    else {
        for (uint i = start + 1; i < start + offset - 1; i++) {
            glBegin(GL_LINE_LOOP);
            vertex = &this->scene->vertices.at(start);
            normal = &this->scene->normals.at(start);
            glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
            glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);

            vertex = &this->scene->vertices.at(i);
            normal = &this->scene->normals.at(i);
            glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
            glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);

            vertex = &this->scene->vertices.at(i + 1);
            normal = &this->scene->normals.at(i + 1);
            glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
            glNormal3f((*normal)[0], (*normal)[1], (*normal)[2]);
            glEnd();

            // If we're in normal mode, draw it sticking out of the surface
            if (this->ui->normal_mode) {
                glBegin(GL_LINES);
                vertex = &this->scene->vertices.at(i);
                normal = &this->scene->normals.at(i);
                glVertex3f((*vertex)[0], (*vertex)[1], (*vertex)[2]);
                endpoint = *vertex + *normal * 0.25;
                glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
                glEnd();
            }
        }
    }
}

void Renderer::drawObject(Object* obj, int depth) {
    glPushMatrix();

    const vector<Transformation>& overall_trans =
        obj->getOverallTransformation();
    for (int i = overall_trans.size() - 1; i >= 0; i--) {
        transform(overall_trans.at(i));
    }

    for (auto& child_it : obj->getChildren()) {
        glPushMatrix();

        const vector<Transformation>& child_trans = 
            child_it.second.transformations;
        for (int i = child_trans.size() - 1; i >= 0; i--) {
            transform(child_trans.at(i));
        }
        draw(Renderable::get(child_it.second.name), depth + 1);

        glPopMatrix();
    }

    glPopMatrix();
}

void Renderer::drawAxes() {
    const float red_bright[3] = {1.0, 0.0, 0.0};
    const float green_bright[3] = {0.0, 1.0, 0.0};
    const float blue_bright[3] = {0.0, 0.0, 1.0};
    const float red_light[3] = {0.2, 0.0, 0.0};
    const float green_light[3] = {0.0, 0.2, 0.0};
    const float blue_light[3] = {0.0, 0.0, 0.2};
    const float black[3] = {0.0, 0.0, 0.0};

    glMaterialfv(GL_FRONT, GL_DIFFUSE, black);
    glMaterialfv(GL_FRONT, GL_SPECULAR, black);
    glMaterialf(GL_FRONT, GL_SHININESS, 1.0);


    // positive x-axis
    glMaterialfv(GL_FRONT, GL_AMBIENT, red_bright);
    glLineWidth(2.5);
    glBegin(GL_LINES);
    {
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(100.0, 0.0, 0.0);
    }
    glEnd();
    // negative x-axis
    glMaterialfv(GL_FRONT, GL_AMBIENT, red_light);
    glLineWidth(0.5);
    glBegin(GL_LINES);
    {
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(-100.0, 0.0, 0.0);
    }
    glEnd();

    // positive y-axis
    glMaterialfv(GL_FRONT, GL_AMBIENT, green_bright);
    glLineWidth(2.5);
    glBegin(GL_LINES);
    {
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 100.0, 0.0);
    }
    glEnd();
    // negative y-axis
    glMaterialfv(GL_FRONT, GL_AMBIENT, green_light);
    glLineWidth(0.5);
    glBegin(GL_LINES);
    {
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, -100.0, 0.0);
    }
    glEnd();

    // positive z-axis
    glMaterialfv(GL_FRONT, GL_AMBIENT, blue_bright);
    glLineWidth(2.5);
    glBegin(GL_LINES);
    {
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 100.0);
    }
    glEnd();
    // negative z-axis
    glMaterialfv(GL_FRONT, GL_AMBIENT, blue_light);
    glLineWidth(0.5);
    glBegin(GL_LINES);
    {
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, -100.0);
    }
    glEnd();
}

void Renderer::transform(const Transformation& trans) {
    switch (trans.type) {
        case TRANS:
            glTranslatef(trans.trans[0], trans.trans[1], trans.trans[2]);
            break;
        case ROTATE:
            glRotatef(
                radToDeg(trans.trans[3]),
                trans.trans[0],
                trans.trans[1],
                trans.trans[2]);
            break;
        case SCALE:
            glScalef(trans.trans[0], trans.trans[1], trans.trans[2]);
            break;
        default:
            fprintf(stderr, "Renderer::transform ERROR invalid TransformationType %d\n",
                trans.type);
            exit(1);
    }
}

/* Renders the scene in its current state. */
void Renderer::display() {
    Renderer *renderer = Renderer::getSingleton();
    UI *ui = renderer->ui;
    Scene *scene = renderer->scene;

    // Clear the color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Make sure there aren't any pending ui changes
    renderer->checkUIState();

    // Preserve the current modelview matrix and apply our scene scaling
    glPushMatrix();
    glScalef(ui->scene_scale, ui->scene_scale, ui->scene_scale);

    // Preserve the current modelview matrix, then rotate the lights using their
    // arcball matrix, popping back to the old modelview matrix
    glPushMatrix();
    glMultMatrixf(ui->arcball_light_mat.data());
    renderer->setupLights();

    const float black[3] = {0.0, 0.0, 0.0};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, black);
    glMaterialfv(GL_FRONT, GL_SPECULAR, black);
    glMaterialf(GL_FRONT, GL_SHININESS, 1.0);

    glPointSize(5.0);
    for (uint i = 0; i < scene->lights.size(); i++) {
        glMaterialfv(GL_FRONT, GL_AMBIENT,
            (const float *) &scene->lights[i].color);
        glBegin(GL_POINTS);
        glVertex3f(scene->lights[i].position[0], scene->lights[i].position[1],
            scene->lights[i].position[2]);
        glEnd();
    }

    renderer->drawAxes();
    
    glPopMatrix();

    // Do the same for the objects and their matrix
    glPushMatrix();
    glMultMatrixf(ui->arcball_object_mat.data());

    // Bind the first vertex buffer as the active one
    glBindBuffer(GL_ARRAY_BUFFER, renderer->vb_objects[0]);
    // Copy the array of vertex values to the buffer
    glBufferData(GL_ARRAY_BUFFER, scene->vertices.size() * sizeof(Vec3f),
        scene->vertices.data(), GL_STATIC_DRAW);
    // Set the buffer as a list of groups of 3 floats, and enable it
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    // Do the same for the second buffer, filling it with the normal vectors
    glBindBuffer(GL_ARRAY_BUFFER, renderer->vb_objects[1]);
    glBufferData(GL_ARRAY_BUFFER, scene->normals.size() * sizeof(Vec3f),
        scene->normals.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);
    
    // Draw the .obj entities in the scene, then draw the primitives starting
    // from the end of the .obj data in the vertex buffers

    // Traverse the tree of objects and draw them all
    if (scene->root_objs.size() != 0) {
        for (Object* obj : scene->root_objs) {
            renderer->drawObject(obj, 1);
        }
    } else {
        for (auto& prm_it : scene->prm_tessellation_start) {
            renderer->drawPrimitive(prm_it.first);
        }
    }

    if (ui->wireframe_mode) {
        if (ui->io_mode) {
            glPushMatrix();
            Assignment::drawIOTest();
            glPopMatrix();
        }
    }
    else if (ui->intersect_mode) {
        glPushMatrix();
        Assignment::drawIntersectTest(&ui->camera);
        glPopMatrix();
    }

    // Pop the arcball and scaling matrices
    glPopMatrix();
    glPopMatrix();

    // Display the current scene
    glutSwapBuffers();

    if (CommandLine::active()) {
        printf("> ");
        CommandLine::readLine(cin);
        scene->update();
        glutPostRedisplay();
    }
}

/* Reshapes the window. */
void Renderer::reshape(int xres, int yres) {
    UI *ui = Renderer::getSingleton()->ui;

    ui->reshape(xres, yres);
    glViewport(0, 0, xres, yres);

    // Reset the perspective projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    ui->camera.aspect = (float) xres / yres;
    gluPerspective(ui->camera.fov, ui->camera.aspect, ui->camera.near,
        ui->camera.far);

    glMatrixMode(GL_MODELVIEW);
    glutPostRedisplay();
}
