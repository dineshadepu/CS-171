#include "model.hpp"

/********************************* RGBf Struct ********************************/

RGBf::RGBf(const float r, const float g, const float b) {
    this->r = r;
    this->g = g;
    this->b = b;
}

/**************************** Transformation Struct ***************************/

Transformation::Transformation(TransformationType type, const Vector4f& trans) :
    type(type), trans(trans)
{
    
}
Transformation::Transformation(
    TransformationType type,
    const float x,
    const float y,
    const float z,
    const float w) : type(type), trans()
{
    trans << x, y, z, w;
}

/********************************* Name Struct ********************************/

Name::Name(const char* name) : name() {
    strcpy(this->name, name);
}

const bool Name::operator==(const Name& rhs) const {
    for (unsigned int i = 0; i < name_buffer_size; i++) {
        if (this->name[i] != rhs.name[i]) {
            return false;
        }
    }
    return true;
}
const bool Name::operator!=(const Name& rhs) const {
    for (unsigned int i = 0; i < name_buffer_size; i++) {
        if (this->name[i] != rhs.name[i]) {
            return true;
        }
    }
    return false;
}

/******************************* Renderable Class *****************************/

unordered_map<Name, Renderable*, NameHasher> Renderable::renderables;

Renderable::Renderable() {
    
}
Renderable::~Renderable() {
    
}

// static instance controller functions
Renderable* Renderable::create(RenderableType type, const Name& name) {
    if (exists(name)) {
        fprintf(stderr, "Renderable::create ERROR name %s already exists with type %s\n",
            name.name, toCstr(type));
        return get(name);
    }

    Renderable* new_renderable = NULL;
    switch (type) {
        case PRM:
            printf("creating new Primitive with name %s\n", name.name);
            new_renderable = new Primitive();
            break;
        case OBJ:
            printf("creating new Object with name %s\n", name.name);
            new_renderable = new Object();
            break;
        // case OBJ:
        //     printf("creating new Object with name %s\n", name.name);
        //     new_renderable = new Object();
        //     break;
        default:
            fprintf(stderr, "ERROR Renderable::create invalid RenderableType %d\n",
                type);
    }
    renderables[name] = new_renderable;

    return new_renderable;
}
Renderable* Renderable::get(const Name& name) {
    if (exists(name)) {
        return renderables[name];
    }
    return NULL;
}
bool Renderable::exists(const Name& name) {
    if (renderables.find(name) != renderables.end()) {
        return true;
    }
    return false;
}
void Renderable::clear() {
    renderables.clear();
}
const unordered_map<Name, Renderable*, NameHasher>&
    Renderable::getActiveRenderables()
{
    return renderables;
}

/******************************** Primitive Class *****************************/

Primitive::Primitive() :
    Renderable(),
    coeff(),
    exp0(default_exp0),
    exp1(default_exp1),
    patch_x(default_patch_x),
    patch_y(default_patch_y),
    color(default_color),
    ambient(default_ambient),
    reflected(default_reflected),
    refracted(default_refracted),
    gloss(default_gloss),
    diffuse(default_diffuse),
    specular(default_specular)
{
    coeff << default_coeff_x, default_coeff_y, default_coeff_z;
}

Primitive::~Primitive() {

}

// modifiers for private variables
void Primitive::setCoeff(const Vector3f& coeff) {
    this->coeff = coeff;
}
void Primitive::setCoeff(const float x, const float y, const float z) {
    this->coeff << x, y, z;
}
void Primitive::setExponents(const float exp0, const float exp1) {
    this->exp0 = exp0;
    this->exp1 = exp1;
}
void Primitive::setPatch(
    const unsigned int patch_x,
    const unsigned int patch_y)
{
    this->patch_x = patch_x;
    this->patch_y = patch_y;
}
void Primitive::setColor(const RGBf& color) {
    this->color = color;
}
void Primitive::setColor(const float r, const float g, const float b) {
    this->color.r = r;
    this->color.g = g;
    this->color.b = b;
}
void Primitive::setAmbient(const float ambient) {
    this->ambient = ambient;
}
void Primitive::setReflected(const float reflected) {
    this->reflected = reflected;
}
void Primitive::setRefracted(const float refracted) {
    this->refracted = refracted;
}
void Primitive::setGloss(const float gloss) {
    this->gloss = gloss;
}
void Primitive::setDiffuse(const float diffuse) {
    this->diffuse = diffuse;
}
void Primitive::setSpecular(const float specular) {
    this->specular = specular;
}

// accessors for private variables
const Vector3f& Primitive::getCoeff() const {
    return this->coeff;
}
const float Primitive::getExp0() const {
    return this->exp0;
}
const float Primitive::getExp1() const {
    return this->exp1;
}
const unsigned int Primitive::getPatchX() const {
    return this->patch_x;
}
const unsigned int Primitive::getPatchY() const {
    return this->patch_y;
}
const RGBf& Primitive::getColor() const {
    return this->color;
}
float Primitive::getAmbient() const {
    return this->ambient;
}
float Primitive::getReflected() const {
    return this->reflected;
}
float Primitive::getRefracted() const {
    return this->refracted;
}
float Primitive::getGloss() const {
    return this->gloss;
}
float Primitive::getDiffuse() const {
    return this->diffuse;
}
float Primitive::getSpecular() const {
    return this->specular;
}

/*
 * Converts from a parametric point (u, v) on a primitive's surface to world
 * space.
 */ 
const Vector3f Primitive::getVertex(float u, float v) {
    float cos_v = pCos(v, this->exp1);
    float sin_v = pSin(v, this->exp1);
    Vector3f p;
    p << this->coeff[0] * cos_v * pCos(u, this->exp0),
        this->coeff[1] * cos_v * pSin(u, this->exp0),
        this->coeff[2] * sin_v;
    return p; 
}

/* Gets a primitive's surface normal at a location (on its surface). */
const Vector3f Primitive::getNormal(const Vector3f& vertex) {
    float x = vertex(0) / this->coeff[0];
    float y = vertex(1) / this->coeff[1];
    float z = vertex(2) / this->coeff[2];
    Vector3f normal;

    normal(0) = (x == 0.0) ? 0.0 :
        2.0 * x * powf(x * x, 1.0 / this->exp0 - 1.0) *
        powf(powf(x * x, 1.0 / this->exp0) + powf(y * y, 1.0 / this->exp0),
            this->exp0 / this->exp1 - 1.0) / this->exp1;
    normal(1) = (y == 0.0) ? 0.0 :
        2.0 * y * powf(y * y, 1.0 / this->exp0 - 1.0) *
        powf(powf(x * x, 1.0 / this->exp0) + powf(y * y, 1.0 / this->exp0),
            this->exp0 / this->exp1 - 1.0) / this->exp1;
    normal(2) = (z == 0.0) ? 0.0 :
        2.0 * z * powf(z * z, 1.0 / this->exp1 - 1.0) / this->exp1;

    normal[0] /= this->coeff[0];
    normal[1] /= this->coeff[1];
    normal[2] /= this->coeff[2];
    normal.normalize();

    return normal;
}

/********************************* Object Class *******************************/

Child::Child() : name("default child"), transformations() {
    fprintf(stderr, "Child ERROR Child default constructor called\n");
    exit(1);
}
Child::Child(const Name& name) : name(name), transformations() {

}

Object::Object() :
    Renderable(),
    // display(false),
    transformations(),
    children(),
    cursor(default_cursor)
{

}

Object::~Object() {
    // this->children_obj.clear();
    // this->children_prm.clear();
    // ^ can't do this. this will cause infinite loop. we need a static class
    // wide destruction behavior
}

// vector<Object*> Object::getRootObjs() {
//     vector<Object*> root_objs;

//     for (const auto& ren_it : Renderable::getActiveRenderables()) {
//         if (ren_it.second->getType() == OBJ) {
//             Object* obj = dynamic_cast<Object*>(ren_it.second);
//             if (obj->getDisplayStatus() == true) {
//                 root_objs.push_back(obj);
//             }
//         }
//     }

//     return root_objs;
// }

// modifiers for private variables

// void Object::show() {
//     this->display = true;
// }
// void Object::hide() {
//     this->display = false;
// }
// bool Object::getDisplayStatus() {
//     return this->display;
// }

bool Object::aliasExists(const Name& name) {
    return this->children.find(name) != this->children.end();
}

// for overall transformation
void Object::overallTranslate(const float x, const float y, const float z) {
    this->transformations.push_back(Transformation(TRANS, x, y, z, 1));
}
void Object::overallRotate(
    const float x, 
    const float y, 
    const float z,
    const float theta)
{
    Vector3f rotate;
    rotate << x, y, z;
    rotate.normalize();
    this->transformations.push_back(
        Transformation(ROTATE, rotate[0], rotate[1], rotate[2], theta));
}
void Object::overallScale(const float x, const float y, const float z) {
    this->transformations.push_back(Transformation(SCALE, x, y, z, 1));
}

const vector<Transformation>& Object::getOverallTransformation() const {
    return this->transformations;
}

void Object::addChild(const Name& name, const Name& alias) {
    if (Renderable::exists(name)) {
        this->children.insert({alias, Child(name)});
        this->cursor = alias;
    } else {
        fprintf(stderr, "Object::addChild ERROR Renderable with name %s does not exist\n",
            name.name);
    }
}

const unordered_map<Name, Child, NameHasher>& Object::getChildren() const
{
    return this->children;
}

// set cursor to given alias
void Object::setCursor(const Name& alias) {
    if (this->children.find(alias) == this->children.end()) {
        fprintf(stderr, "Object::validateCursor ERROR child alias %s does not exist\n",
            alias.name);
    } else {
        this->cursor = alias;
    }
}

const Name& Object::getCursor() const {
    return this->cursor;
}

bool Object::validateCursor() {
    if (this->children.find(this->cursor) == this->children.end()) {
        fprintf(stderr, "Object::validateCursor ERROR child with alias %s does not exist\n",
            this->cursor.name);
        return false;
    }
    return true;
}
void Object::cursorTranslate(const float x, const float y, const float z) {
    if (this->validateCursor()) {
        this->children[this->cursor].transformations.push_back(
            Transformation(TRANS, x, y, z, 1));
    }
}
void Object::cursorRotate(
    const float x, 
    const float y, 
    const float z,
    const float theta)
{
    if (this->validateCursor()) {
        Vector3f rotate;
        rotate << x, y, z;
        rotate.normalize();
        this->children[this->cursor].transformations.push_back(
            Transformation(ROTATE, rotate[0], rotate[1], rotate[2], theta));
    }
}
void Object::cursorScale(const float x, const float y, const float z) {
    if (this->validateCursor()) {
        this->children[this->cursor].transformations.push_back(
            Transformation(SCALE, x, y, z, 1));
    }
}

/*************************** PrintInfo Helper Functions ***********************/

void printSceneInfo(int indent) {
    printf("printing scene info\n");

    // print primitives
    unordered_map<Name, Renderable*, NameHasher>::const_iterator ren_it;
    const unordered_map<Name, Renderable*, NameHasher>& renderables =
        Renderable::getActiveRenderables();
    if (renderables.size() > 0) {
        printIndent(indent + 1);
        printf("currently active Primitive(s):\n");
        for (ren_it = renderables.begin();
            ren_it != renderables.end();
            ren_it++)
        {
            if (ren_it->second->getType() == PRM) {
                printIndent(indent + 2);
                printf("%s\n", ren_it->first.name);
            }
        }
    }

    // print objects
    if (renderables.size() > 0) {
        printIndent(indent + 1);
        printf("currently active Object(s):\n");
        for (ren_it = renderables.begin();
            ren_it != renderables.end();
            ren_it++)
        {
            if (ren_it->second->getType() == OBJ) {
                printIndent(indent + 2);
                printf("%s\n", ren_it->first.name);
            }
        }
    }

    printf("DONE\n");
}

void printInfo(const Renderable* ren, int indent) {
    assert(ren);
    switch (ren->getType()) {
        case PRM:
        {
            const Primitive* prm = dynamic_cast<const Primitive*>(ren);
            printIndent(indent);
            const Vector3f& coeff = prm->getCoeff();
            printf("coeff: %f %f %f\n", coeff[0], coeff[1], coeff[2]);
            printIndent(indent);
            printf("exponents: %f %f\n", prm->getExp0(), prm->getExp1());
            printIndent(indent);
            printf("patch: %d %d\n", prm->getPatchX(), prm->getPatchY());
            printInfo(prm->getColor(), indent);
            printIndent(indent);
            printf("ambient: %f\n", prm->getAmbient());
            printIndent(indent);
            printf("reflected: %f\n", prm->getReflected());
            printIndent(indent);
            printf("refracted: %f\n", prm->getRefracted());
            printIndent(indent);
            printf("gloss: %f\n", prm->getGloss());
            printIndent(indent);
            printf("diffuse: %f\n", prm->getDiffuse());
            printIndent(indent);
            printf("specular: %f\n", prm->getSpecular());
            break;
        }
        case OBJ:
        {
            const Object* obj = dynamic_cast<const Object*>(ren);
            // print overall transformation
            printIndent(indent);
            printf("overall transformation(s):\n");
            printInfo(obj->getOverallTransformation(), indent + 1);

            // print children objects
            printIndent(indent);
            printf("child object(s):\n");
            unordered_map<Name, Child, NameHasher>::const_iterator child_it;    // FIX HERE use const auto& for each loop
            for (child_it = obj->getChildren().begin();
                child_it != obj->getChildren().end();
                child_it++)
            {
                assert(Renderable::get(child_it->second.name));
                if (Renderable::get(child_it->second.name)->getType() == OBJ) {
                    printIndent(indent + 1);
                    printf("OBJ %s", child_it->second.name.name);
                    if (child_it->second.name != child_it->first.name) {
                        printf(" aliased as %s", child_it->first.name);
                    }
                    printf("\n");
                }
            }

            // print children primitives
            printIndent(indent);
            printf("child primitive(s):\n");
            for (child_it = obj->getChildren().begin();
                child_it != obj->getChildren().end();
                child_it++)
            {
                assert(Renderable::get(child_it->second.name));
                if (Renderable::get(child_it->second.name)->getType() == PRM) {
                    printIndent(indent + 1);
                    printf("PRM %s", child_it->second.name.name);
                    if (child_it->second.name != child_it->first.name) {
                        printf(" aliased as %s", child_it->first.name);
                    }
                    printf("\n");
                }
            }

            // print cursor info
            printIndent(indent);
            printf("current cursor set at: %s\n", obj->getCursor().name);
            if (obj->getCursor() != default_cursor) {
                printIndent(indent);
                printf("transformation(s) on cursor renderable:\n");
                printInfo(
                    obj->getChildren().at(obj->getCursor()).transformations,
                    indent + 1);
            }
            break;
        }
        default:
            fprintf(stderr, "printInfo ERROR invalid Renderable type %d\n",
                ren->getType());
            exit(1);
    }
}

void printInfo(const RGBf& color, int indent) {
    printIndent(indent);
    printf("color: %f %f %f\n", color.r, color.g, color.b);
}

void printInfo(const vector<Transformation>& trans, int indent) {
    vector<Transformation>::const_iterator it;
    for (it = trans.begin(); it != trans.end(); it++) {
        printIndent(indent);
        switch (it->type) {
            case TRANS:
                assert(it->trans[3] == 1);
                printf("TRANS %f %f %f\n",
                    it->trans[0], it->trans[1], it->trans[2]);
                break;
            case ROTATE:
            {
                printf("ROTATE %f %f %f %f\n",
                    it->trans[0],
                    it->trans[1],
                    it->trans[2],
                    it->trans[3] * 180.0 / M_PI);
                break;
            }
            case SCALE:
                assert(it->trans[3] == 1);
                printf("SCALE %f %f %f\n",
                    it->trans[0], it->trans[1], it->trans[2]);
                break;
            default:
                fprintf(stderr, "ERROR printInfo invalid TransformationType %d\n",
                    it->type);
                exit(1);
        }
    }
}
