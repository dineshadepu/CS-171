#include "Utilities.hpp"

/* Constructs a triple of integers. */
Vec3i::Vec3i(int i, int j, int k) {
    this->i = i;
    this->j = j;
    this->k = k;
}

/* Constructs a triple of integers. */
Vec3i::Vec3i(int *t) {
    *this = *((Vec3i *) t);
}

/* Constructs a point in 3D space. */
Vec3f::Vec3f(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

/* Constructs a point in 3D space. */
Vec3f::Vec3f(float *v) {
    *this = *((Vec3f *) v);
}

/* Constructs a point in homogeneous space. */
Vec4f::Vec4f(float x, float y, float z, float w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

/* Constructs a point in homogeneous space. */
Vec4f::Vec4f(float *v) {
    *this = *((Vec4f *) v);
}

/* Converts a value in degrees to one in radians. */
float degToRad(float angle) {
    return angle * M_PI / 180.0;
}

/* Converts a value in radians to one in degrees. */
float radToDeg(float angle) {
    return angle * 180.0 / M_PI;
}

/*
 * Creates a rotation matrix in place, given the rotation's orientation axis and
 * angle.
 */
void makeRotateMat(float *matrix, float x, float y, float z, float angle) {
    // Normalize the rotation axis
    float length = sqrtf(x * x + y * y + z * z);
    x = x / length;
    y = y / length;
    z = z / length;
    
    // Set the matrix's values
    float sin_a = sinf(angle);
    float cos_a = cosf(angle);
    matrix[0] = x * x + (1 - x * x) * cos_a;
    matrix[1] = y * x * (1 - cos_a) + z * sin_a;
    matrix[2] = z * x * (1 - cos_a) - y * sin_a;
    matrix[3] = 0.0;
    matrix[4] = x * y * (1 - cos_a) - z * sin_a;
    matrix[5] = y * y + (1 - y * y) * cos_a;
    matrix[6] = z * y * (1 - cos_a) + x * sin_a;
    matrix[7] = 0.0;
    matrix[8] = x * z * (1 - cos_a) + y * sin_a;
    matrix[9] = y * z * (1 - cos_a) - x * sin_a;
    matrix[10] = z * z + (1 - z * z) * cos_a;
    matrix[11] = 0.0;
    matrix[12] = 0.0;
    matrix[13] = 0.0;
    matrix[14] = 0.0;
    matrix[15] = 1.0;
}

/* Reads the contents of a file into a char array. */
char *readFile(char *file_name) {
    // Open the file as read-only
    FILE *f = fopen(file_name, "rb");
    if (!f)
        return NULL;

    // Seek to the end of the file and figure out how long it is, then go back
    // to the beginning
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    // Allocate memory to store the file's contents, then copy them, checking
    // that the file wasn't empty
    char *contents = (char *) malloc(fsize + 1);
    assert(fread(contents, fsize, 1, f) > 0);
    // Set the last character to 0 (EOF), and close the file
    contents[fsize] = 0;
    fclose(f);

    // Return a pointer to the array
    return contents;
}
