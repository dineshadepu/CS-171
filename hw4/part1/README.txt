### Part 1 ###

Implementing Phong Lighting/Shading in OpenGL

---------------
Usage

make
./opengl_renderer [scene_description_file.txt] [xres] [yres] [mode]

---------------
Input

xres, yres = integers for opengl display window resolution
scene file is in format seen before
mode = 0 (Gouraud) or 1 (Phong)

---------------
Info

This program implements per-pixel Phong shading in OpenGL by overwriting the vertex and fragment shaders (look at readShaders function and the separate GLSL shader program files).

Included files:
	parsing.cpp, parsing.h - for parsing the txt file
	structs.h - to hold structs, both given from demo and created
	data - holds a bunch of txt files and objects

I changed the main function so that we just use one parse all function for cleanliness 
and pass all variables in by reference. Arcball (quaternion) works fine.