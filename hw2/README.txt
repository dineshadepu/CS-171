### Hw 2 ###

3D Shaded Surface Renderer

---------------
Usage

make
./shaded [data/scene_description_file.txt] [xres] [yres] [mode]

---------------
Input

xres, yres = ints for resolution (optimal 800x800)
mode = 0 or 1 boolean for choosing shading


scene txt file format example:
camera:  
position -2 0 5  
orientation 0 1 0 0  
near 1  
far 10  
left -0.5  
right 0.5  
top 0.5  
bottom -0.5  
 
light -0.8 0 1 , 1 1 0 , 0  
light 0.15 0.85 0.7 , 1 0 1 , 0  
light 0.5 -0.5 0.85 , 0 1 1 , 0.5  
 
objects:  
cube cube.obj  
 
cube  
ambient 0.2 0 0  
diffuse 0.8 0 0  
specular 0 0 0  
shininess 0.2  
s 1 1 1  
r 0 1 0 0.6  
t 0 0 0  
s 1 1 1  
r 1 0 0 0.75  
t -2.7 0 -2 
...


Obj file format example:
v -1 -1 1  
v 1 -1 1  
v 1 1 1  
v -1 1 1  
v -1 -1 -1  
v 1 -1 -1  
v 1 1 -1  
v -1 1 -1
... 
vn 0 0 1  
vn 0 0 -1  
vn 0 1 0  
vn 0 -1 0  
vn 1 0 0  
vn -1 0 0
... 
f 1//1 2//1 3//1  
f 1//1 3//1 4//1  
f 6//2 5//2 7//2  
f 7//2 5//2 8//2  
f 2//5 6//5 3//5  
f 3//5 6//5 7//5  
f 5//6 4//6 8//6  
f 4//6 5//6 1//6  
f 4//3 3//3 8//3  
f 7//3 8//3 3//3  
f 1//4 5//4 2//4  
f 2//4 5//4 6//4
...

---------------
Info

My code:
- Parses scene description file and stores data into structs for camera, lights, objects and their copies
- Creates transformation and perspective projection matrices
- Transforms each point by all geometric transformations to object copy
- Transforms surface normals by normal transformations
- Lights and shades the meshes by rasterizing colored triangles using Gouraud or Phong shading
- Outputs pixel grid as ppm image

For more information about implementation, input scene/file structure, and how I'm shading, go to http://courses.cms.caltech.edu/cs171/assignments/hw2/hw2-html/cs171hw2.html





