#include "PNGMaker.hpp"

#define MAX_INTENSITY 255   // Maximum color intensity

/* Creates a PNGMaker with the specified resoltuion. */
PNGMaker::PNGMaker(int xres, int yres) {
    this->screen = (Vec3f *) malloc(xres * yres * sizeof(Vec3f));
    this->xres = xres;
    this->yres = yres;
}

/* Frees the PNGMaker's pixel grid. */
PNGMaker::~PNGMaker() {
    free(this->screen);
}

/* Sets the color of a pixel in the PNGMaker's screen. */
void PNGMaker::setPixel(int x, int y, float r, float g, float b) {
    Vec3f &pixel = this->screen[y * this->xres + x];
    pixel.x = r;
    pixel.y = g;
    pixel.z = b;
}

/* Saves the PNGMaker's current stored screen to a PNG image file. */
int PNGMaker::saveImage() {
    FILE *fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_byte **row_pointers = NULL;
    
    // Open a handle for the file we're going to write to
    fp = fopen("./rt.png", "wb");
    if (!fp) {
        return 1;   
    }

    // Create the data and info structures used by libpng
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info_ptr = png_create_info_struct(png_ptr);
    if (!png_ptr || !info_ptr)
        return 1;
    
    // Set up libpng to write a XRES x YRES image with 8 bits of RGB color depth
    png_set_IHDR(png_ptr, info_ptr, this->xres, this->yres, 8,
        PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT);

    // Uncomment this if you want to save hard drive space; it takes a lot
    // longer for the CPU to write the PNGs, but they're significantly smaller
    // png_set_compression_level(png_ptr, Z_BEST_COMPRESSION);

    // Allocate an array of pointers to each row of pixels in the PNG file
    row_pointers = (png_byte **) png_malloc(png_ptr,
        this->yres * sizeof(png_byte *));
    for (int y = 0; y < this->yres; y++) {
        // Allocate a row of pixels in the PNG file
        png_byte *row = (png_byte *) png_malloc(png_ptr,
            3 * this->xres * sizeof(uint8_t));
        // The screen matrix has (x, y) = (0, 0) in the lower left so that it's
        // easy to translate between it and the "camera sensor" grid, so we have
        // to go through the rows in reverse because of convention
        row_pointers[this->yres - 1 - y] = row;
        for (int x = 0; x < this->xres; x++) {
            // Get the pixel's location and fill in its red, green, and blue
            // values, bounding them between 0 and MAX_INTENSITY = 255
            Vec3f &pixel = this->screen[y * this->xres + x];
            *row++ = (uint8_t) (fmin(fmax(0.0, pixel.x), 1.0) * MAX_INTENSITY);
            *row++ = (uint8_t) (fmin(fmax(0.0, pixel.y), 1.0) * MAX_INTENSITY);
            *row++ = (uint8_t) (fmin(fmax(0.0, pixel.z), 1.0) * MAX_INTENSITY);
        }
    }

    // Write the PNG to file
    png_init_io(png_ptr, fp);
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    
    // Free each row of the PNG that we allocated
    for (int y = 0; y < this->yres; y++) {
        png_free(png_ptr, row_pointers[y]);
    }
    // Free the array of row pointers
    png_free(png_ptr, row_pointers);
    fclose(fp);

    return 0;
}