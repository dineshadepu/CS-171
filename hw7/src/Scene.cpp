#include "Scene.hpp"

using namespace std;

Scene *Scene::singleton;

/* Constructs a light in our scene. */
PointLight::PointLight(float *position, float *color, float k) {
    for (int i = 0; i < 3; i++) {
        this->position[i] = position[i];
        this->color[i] = color[i];
    }
    this->position[3] = position[3];
    this->k = k;
}

/* Initializes the scene's data structures. */
Scene::Scene() : needs_update(default_needs_update) {
    this->prm_tessellation_start = unordered_map<Primitive*, unsigned int>();
    // Scene::objects = vector<Object>();
    this->vertices = vector<Vector3f>();
    this->normals = vector<Vector3f>();
    this->lights = vector<PointLight>();
    this->root_objs = vector<Object *>();

    this->createLights();
}

/* Returns/sets up the singleton instance of the class. */
Scene *Scene::getSingleton() {
    if (!Scene::singleton) {
        Scene::singleton = new Scene();
    }

    return Scene::singleton;
}

/* Creates the scene's lights. */
void Scene::createLights() {
    float position[4] = {3.0, 4.0, 5.0, 1.0};
    float color[3] = {1.0, 1.0, 1.0};
    float k = 0.2;

    // Set up a single light
    this->lights.emplace_back(position, color, k);
}

/*
 * Adds a vertex and normal to the respective buffers, calculated from a 
 * parametric (u, v) point on a superquadric's surface.
 */
void Scene::generateVertex(Primitive *prm, float u, float v) {
    this->vertices.push_back(prm->getVertex(u, v));
    this->normals.push_back(prm->getNormal(this->vertices.back()));
}

/* Tesselates a primitive, adding the vertices and normals to their buffers. */
void Scene::tessellatePrimitive(Primitive *prm) {
    if (this->prm_tessellation_start.find(prm) != this->prm_tessellation_start.end()) {
        return;
    }

    this->prm_tessellation_start.insert({prm, vertices.size()});

    // Less typing, computation at runtime
    int ures = prm->getPatchX();
    int vres = prm->getPatchY();
    float u, v;
    float half_pi = M_PI / 2;
    float du = 2 * M_PI / ures, dv = M_PI / vres;

    // Create GL_TRIANGLE_STRIPs by moving circumferentially around each of its
    // vres - 2 non-polar latitude ranges
    v = dv - half_pi;
    for (int j = 1; j < vres - 1; j++) {
        // U sweeps counterclockwise from -pi to pi, so the first edge should
        // point down in order for the right-hand rule to make the normal
        // point out of the primitive
        u = -M_PI;
        for (int i = 0; i < ures; i++) {
            this->generateVertex(prm, u, v + dv);
            this->generateVertex(prm, u, v);
            u += du;
        }
        // Connect back to the beginning
        this->generateVertex(prm, -M_PI, v + dv);
        this->generateVertex(prm, -M_PI, v);

        v += dv;
    }

    // Draw the primitive's bottom by filling in its southernmost latitude range
    // with a GL_TRIANGLE_FAN centered on the south pole
    u = M_PI;
    v = dv - half_pi;
    this->generateVertex(prm, u, -half_pi);
    for (int i = 0; i < ures; i++) {
        // U sweeps clockwise to make the normals point out
        this->generateVertex(prm, u, v);
        u -= du;
    }
    // Connect back to the beginning
    u = -M_PI;
    this->generateVertex(prm, u, v);

    // Now we tessellate its top by doing the same at the north pole
    v *= -1;
    this->generateVertex(prm, u, half_pi);
    for (int i = 0; i < ures; i++) {
        // U sweeps counterclockwise to make the normals point out
        this->generateVertex(prm, u, v);
        u += du;
    }
    // Connect back to the beginning
    u = M_PI;
    this->generateVertex(prm, u, v);
}

void Scene::tessellateObject(Object *obj) {
    for (auto& child_it : obj->getChildren()) {
        Renderable* ren = Renderable::get(child_it.second.name);
        switch (ren->getType()) {
            case OBJ: {
                Object* obj = dynamic_cast<Object*>(ren);
                this->tessellateObject(obj);
                break;
            }
            case PRM: {
                Primitive* prm = dynamic_cast<Primitive*>(ren);
                this->tessellatePrimitive(prm);
                break;
            }
            default:
                fprintf(stderr, "Scene::tessellateObject ERROR invalid Renderable type %s\n",
                    toCstr(ren->getType()));
                exit(1);
        }
    }
}

/*
 * Regenerates the scene's vertex and normal buffers based on currently selected
 * Renderable
 */
void Scene::update() {
    this->root_objs.clear();
    this->prm_tessellation_start.clear();
    this->vertices.clear();
    this->normals.clear();

    // Add each of the objects, then each of the primitives
    // for (uint i = 0; i < objects.size(); i++)
    //     setupObject(objects.data() + i);

    const Line* cur_state = CommandLine::getState();

    if (cur_state) {
        switch (cur_state->toCommandID()) {
            case Commands::primitive_get_cmd_id: {
                Renderable* ren = Renderable::get(cur_state->tokens[1]);
                assert(ren->getType() == PRM);
                this->tessellatePrimitive(dynamic_cast<Primitive*>(ren));
                break;
            }
            case Commands::object_get_cmd_id: {
                Renderable* ren = Renderable::get(cur_state->tokens[1]);
                assert(ren->getType() == OBJ);
                this->root_objs.push_back(dynamic_cast<Object*>(ren));
                break;
            }
            default:
                fprintf(stderr, "ERROR Commands:info invalid state CommandID %d from current state\n",
                    cur_state->toCommandID());
                exit(1);
        }
    }

    for (Object* obj : root_objs) {
        this->tessellateObject(obj);
    }
}




// /* Initializes an object's data structures. */
// Object::Object() {
//     Object::vertices = vector<Vec3f>();
//     Object::normals = vector<Vec3f>();
//     Object::faces = vector<Vec3i>();
// }


// /* Adds an object to the current scene's vertex buffers. */
// void Scene::setupObject(Object *object) {
//     // Object faces are stored as alternating triples of vertex and normal
//     // indices, so we loop through them all and put each pair into their buffers
//     Vec3i *f = object->faces.data();
//     for (uint i = 0; i < object->faces.size(); i += 2) {
//         vertices.push_back(object->vertices[f->i]);
//         vertices.push_back(object->vertices[f->j]);
//         vertices.push_back(object->vertices[f->k]);
//         f++;

//         normals.push_back(object->normals[f->i]);
//         normals.push_back(object->normals[f->j]);
//         normals.push_back(object->normals[f->k]);
//         f++;
//     }
// }

// /* Reads a .obj file and adds its contents to the scene. */
// void Scene::addObject(char *file_name) {
//     // Initialize a new Object at the end of the vector
//     objects.emplace_back(Object());
//     Object *object = &objects.back();

//     // Open the file
//     ifstream obj(file_name);
//     char type;
//     float x, y, z;
//     int v0, v1, v2, n0, n1, n2;
//     // While we can read the next character...
//     while (obj.get(type)) {
//         // If the first character of a line is 'v', this is either a vertex or
//         // a normal
//         if (type == 'v') {
//             // When it's followed by 'n', it's a normal
//             obj.get(type);
//             if (type != 'n') {
//                 obj >> x >> y >> z;
//                 object->vertices.emplace_back(x, y, z);
//             }
//             else {
//                 obj >> x >> y >> z;
//                 object->normals.emplace_back(x, y, z);
//             }
//         }
//         // Otherwise, check if this is a face ('f')
//         else if (type == 'f') {
//             // Read each pair of vertex//normal indices
//             obj >> v0;
//             obj.ignore(1, '/');
//             obj.ignore(1, '/');
//             obj >> n0;

//             obj >> v1;
//             obj.ignore(1, '/');
//             obj.ignore(1, '/');
//             obj >> n1;

//             obj >> v2;
//             obj.ignore(1, '/');
//             obj.ignore(1, '/');
//             obj >> n2;

//             // Convert from 1- to 0-indexing
//             v0--;
//             v1--;
//             v2--;
//             n0--;
//             n1--;
//             n2--;

//             // Add each triple to the face vector, so that it represents
//             // alternating sets of vertices and normals
//             object->faces.emplace_back(v0, v1, v2);
//             object->faces.emplace_back(n0, n1, n2);
//         }
//     }

//     // Objects are expected to precede primitives, so we have to rebuild the
//     // vertex buffers to put this one at the beginning
//     rebuildScene();
// }