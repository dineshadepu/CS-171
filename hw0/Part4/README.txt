### Part 4 ###

A program that outputs a PPM image of a colored circle centered on a different colored background. Use ImageMagick to the view the image that the program produces.

# Usage:
make
./ppm_test xres yres
./ppm_test xres yres | display -
./ppm_test xres yres | convert - my_image_name.png

# Output:
Of the format:
P3  
3 4  
255  
255 255 255  
128 128 128  
0 0 0  
255 0 0  
0 255 0  
0 0 255  
255 255 0  
255 0 255  
0 255 255  
128 0 0  
0 128 0  
0 0 128

