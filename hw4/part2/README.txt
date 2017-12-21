### Part 2 ###

Rendering shading with texturing

---------------
Usage

make
./opengl_renderer [color_texture.png] [normal_map.png]

---------------
Input

xres, yres = integers for opengl display window resolution
color map is for color, normal map is for texture

---------------
Info

This program takes a color texture and normal map and renders the given textures onto a flat polygon.

Run makefile to create executable /opengl_renderer
usage: ./opengl_renderer [color.png] [normals.png]

Included files:
	structs.h - to hold structs, both given from demo and created
	readpng.cpp - from demo code
	vertexProgramEshan.glsl, fragmentProgramEshan.glsl - shaders

I hard coded the camera (which I ended up not using it, I just use default camera 
parameters and use gluPerspective to set things), the light(s), the vertices and texture
 coordinates, and the tangents. The normals were derived from normals.png, and the 
 bitangent was taken from the cross product of the tangent and the normals. Tangents 
 were in a tangent buffer and passed into the shader as an attribute. I also hard code
  what I set the material as, and I hardcode the square that the image is rendered onto
   using two triangles.

In the vertex shader, I pass in the tangent as an attribute, and set varying normal, 
tangent, and bitangent vectors, as well as eye direction, texture coordinate, and gl 
Position. In the fragment shader, I make the TBN matrix, use the normal map to create 
the normal, remap it, and normalize it, , then calculate the color using the lights and 
phong shading. Inside the phong shading loop, TBN matrix is used to transform the light 
and eye to surface coordinates. Because the material diffuse messes with the calculations,
 I hardcode the material diffuse by using a multiplying var instead, and use that in the 
 color calculation.

Comments and potential issues: Because of all the hard coding, the light and material 
properties can seem a bit weird, but in the end the normal mapping and color texture 
mapping comes out correct (with somewhat weird lighting). Also, the arcball doesnt work. 
It could be that the sensitivity of the arcball is too high, or the square that the 
texture is rendering onto isnt a cube so the movement properties are off, or some 
other issue, but that arcball implementation is the same from previous homeworks, and 
my quaternion arcball worked correctly in the past, so the arcball not working in this 
assignment shouldnt be an issue. If you still want to manipulate the image, move slowly, 
and you can see the two sides and light reflection.