#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <utility>
#include <stdexcept>

#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;


inline int sign(float x) {
    return (x > 0) ? 1 : -1;
}

/* Calculates the parametric superquadric sine analog */
inline float pSin(float u, float p) {
    float sin_u = sinf(u);
    return sign(sin_u) * powf(fabs(sin_u), p);
}

/* Calculates the parametric superquadric cosine analog */
inline float pCos(float u, float p) {
    float cos_u = cosf(u);
    return sign(cos_u) * powf(fabs(cos_u), p);
}

// colors
struct RGBf {
    float r;
    float g;
    float b;

    RGBf(const float r, const float g, const float b);
};
static const RGBf default_color(1.0, 1.0, 1.0);

// transformations stored as 3 3f vectors
enum TransformationType {TRANS, SCALE, ROTATE};
struct Transformation {
    TransformationType type;
    Vector4f trans;
    // translation defined by a 3f displacement vector with trailing 1
    // scaling defined by a 3f scaling vector with trailing 1
    // rotation defined by a quaternion (x, y, z, theta)
    Transformation(TransformationType type, const Vector4f& trans);
    Transformation(
        TransformationType type,
        const float x,
        const float y,
        const float z,
        const float w);
    // Fix Eigen packing issue with fixed-size Eigen objects:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

static const unsigned int name_buffer_size = 64;
struct Name {
private:
    // disable default constructor
    Name();

public:
    char name[name_buffer_size];

    Name(const char* name);

    const bool operator==(const Name& rhs) const;
    const bool operator!=(const Name& rhs) const;
};
struct NameHasher {
    unsigned int operator()(const Name& name) const {
        unsigned int hash = 5381;
        for (unsigned int i = 0; i < name_buffer_size; i++) {
            hash = ((hash << 5) + hash) + (unsigned int) name.name[i];
        }

        return hash;
    }
};


enum RenderableType {PRM, OBJ, MSH};
inline const char* toCstr(RenderableType type) {
    switch (type) {
        case PRM:
            return "PRM";
        case OBJ:
            return "OBJ";
        case MSH:
            return "MSH";
        default:
            fprintf(stderr, "ERROR toCstr RenderableType invalid arg %d\n", type);
            exit(1);
    }
}

class Renderable {
private:
    static unordered_map<Name, Renderable*, NameHasher> renderables;

protected:
    explicit Renderable();
    virtual ~Renderable();

public:
    // static instance controller functions
    static Renderable* create(RenderableType type, const Name& name);
    static Renderable* get(const Name& name);
    static bool exists(const Name& name);
    static void clear();
    static const unordered_map<Name, Renderable*, NameHasher>&
        getActiveRenderables();

    virtual const RenderableType getType() const = 0;
};

// class for Primitive information from <modeling language>
// default values
static const float default_coeff_x = 1.0;
static const float default_coeff_y = 1.0;
static const float default_coeff_z = 1.0;
static const float default_exp0 = 1.0;
static const float default_exp1 = 1.0;
static const unsigned int default_patch_x = 15;
static const unsigned int default_patch_y = 15;
static const float default_ambient = 0.1;
static const float default_reflected = 0.3;
static const float default_refracted = 0.5;
static const float default_gloss = 0.3;
static const float default_diffuse = 0.8;
static const float default_specular = 0.1;
class Primitive : public Renderable {
private:
    // coefficients
    Vector3f coeff;
    // exponents
    float exp0;
    float exp1;
    // patch count
    unsigned int patch_x;
    unsigned int patch_y;
    // colors
    RGBf color;
    float ambient;
    float reflected;
    float refracted;
    float gloss;
    float diffuse;
    float specular;

public:
    explicit Primitive();
    ~Primitive();

    const RenderableType getType() const {
        return PRM;
    }

    // modifiers for private variables
    void setCoeff(const Vector3f& coeff);
    void setCoeff(const float x, const float y, const float z);
    void setExponents(const float exp0, const float exp1);
    void setPatch(const unsigned int patch_x, const unsigned int patch_y);
    void setColor(const RGBf& color);
    void setColor(const float r, const float g, const float b);
    void setAmbient(const float ambient);
    void setReflected(const float reflected);
    void setRefracted(const float refracted);
    void setGloss(const float gloss);
    void setDiffuse(const float diffuse);
    void setSpecular(const float specular);

    // accessors for private variables
    const Vector3f& getCoeff() const;
    const float getExp0() const;
    const float getExp1() const;
    const unsigned int getPatchX() const;
    const unsigned int getPatchY() const;
    const RGBf& getColor() const;
    float getAmbient() const;
    float getReflected() const;
    float getRefracted() const;
    float getGloss() const;
    float getDiffuse() const;
    float getSpecular() const;

    const Vector3f getVertex(float u, float v);
    const Vector3f getNormal(const Vector3f& vertex);
};

struct Child {
    Name name;
    vector<Transformation> transformations;

    Child();
    Child(const Name& name);
};

static const Name default_cursor("[NONE]");
// class for Object information from <modeling language>
class Object : public Renderable {
private:
    // overall transformation
    vector<Transformation> transformations;

    // all child objects and primitives
    unordered_map<Name, Child, NameHasher> children;

    // cursor for modifying children
    Name cursor;

public:
    explicit Object();
    ~Object();

    const RenderableType getType() const {
        return OBJ;
    }

    bool aliasExists(const Name& name);

    // overall transformation
    void overallTranslate(const float x, const float y, const float z);
    void overallRotate(
        const float x,
        const float y,
        const float z,
        const float theta);
    void overallScale(const float x, const float y, const float z);
    const vector<Transformation>& getOverallTransformation() const;

    // children objects and primitives
    void addChild(const Name& name, const Name& alias);
    const unordered_map<Name, Child, NameHasher>& getChildren() const;

    // cursor for modifying children
    void setCursor(const Name& alias);
    const Name& getCursor() const;
    bool validateCursor();
    // for object or primitive indicated by cursor
    void cursorTranslate(const float x, const float y, const float z);
    void cursorRotate(
        const float x,
        const float y,
        const float z,
        const float theta);
    void cursorScale(const float x, const float y, const float z);
};

inline void printIndent(int indent) {
    for (int i = 0; i < indent; i++) {
        printf("\t");
    }
}

void printSceneInfo(int indent);
void printInfo(const Renderable* ren, int indent);
void printInfo(const RGBf& color, int indent);
void printInfo(const vector<Transformation>& trans, int indent);

#endif