### Hw 7 ###

Introduction to Ray Tracing

The code in src (run make in src) ray traces a scene. There is a lot of details to this assignment, and to how the command line and viewer works. For more details, check out the assignment page:

http://courses.cms.caltech.edu/cs171/assignments/hw7/hw7-html/cs171hw7.html

In all the given files from src for this assignment, I only added two new funcs:
matrixFunctions.cpp / hpp
Which contain the functions to create transformation, rotation (non quaternion), and
scaling Eigen matrices. These were from hw2, nothing changed in the functions themselves,
just additional functions deleted.

FOLLOW ALONG COMMENTS IN CODE FOR MORE DETAIL ON IMPLEMENTATION.

---------------------

Part 1.1:

My IOTest is simple. I use two additional functions created, one that takes all
the transformations (and scaling from coeff of the primitive) for that object and
makes one combined transformation matrix, and another that calculates the inside out
function S(x,y,z,e,n) using the inverse transformation matrix * the xyz vector and the
prm.

New functions (written above IOtest):

float calcInsideOut - inputs xyz location vector and primitive, returns S(x,y,z,e,n)

Matrix4f getTransfMat - inputs primitive, transformation stack, and a boolean flag,
then does a switch case to combine all transformations in transformation stack
and coeff, returns the combined transformation matrix.

---------------------

Part 1.2:

I make quite a few changes to findIntersection (which is now a void function). This
successfully draws a normal from object to camera on the correct object in any
given scene.

New arguments to findIntersection:
	Camera &camera_ray - ray from camera in world space
	Renderable *ren - current selected ren using block of code adapted from drawIOtest,
		starter for object/primitive tree traversal
	transformation stack vector - holds transformations applied to prm in obj
	float &current_dist - a variable passed by reference that keeps track of current
		shortest distance from the camera position to all found intersections in
		world space
	Ray &current_ray - variable passed by reference that keeps track of the origin
		and normal of the intersection that is currently shortest distance from camera
		in world space
	Primitive ** current_prm - pointer to the pointer for the current primitive at which
		the shortest distance intersection occurs

findIntersection:
I have the skeleton of the code from IOtest that checks the
current renderable to see if its a primitive, and if its an object, recursively runs
on all children. If the renderable is a primitive, I get transformations matrices with
and without translations and transform the incoming camera ray from world space to 
sq space. In the function, that camera ray could be any ray incoming, ie. a light ray
as used in Phong Shading for shadowing later. I then calculate determinant, and if it's
not imaginary, I calculate t_final using a newtonsMethod function I made earlier in the
code if t+ and t- are both positive. I plug in t- to newtons method. If the t final is
positive, I find the intersection coordinates using ray a(t_final)+b in sq space,
transform it to world, normalize it, then take its distance from the world camera
position. If this distance found is shorter than the last shortest distance, reset
it to the new shortest distance and set this intersection point as the intersection ray
and the prm object that is being intersected. Also, after checking distance, I calculate
the intersection normal in sq, change to world, normalize, and make the ray of it. I
also do this all for the case where t+ is positive but t- is negative, in which case
I try newtons method on both t+ and t-, compare distances, and set the distance,
ray, and primitives to the shorter one. If the renderable is an object, I run 
findIntersection recursively on all children.

drawIntersectTest:
I find current renderable for tree, get camera position and direction in world
space, and then set the running variables: transformation stack, distance, prm,
and the intersection ray. These are fed into findIntersection. If a prm is intersected,
I render a thick normal from there.

New functions:

float newtonsMethod - inputs t initial, the primitive, and the direction and positional
vectors from the incoming ray and runs the iterative process. Error is 0.0001. I calculate
initial g(t0) and g'(t0), do a check, then run the iteration process, updating t and
recalculating g and g' using calcInsideOut and calcSQnormal. If it doesn't return t
during the checks, it returns -1.0 which means a miss.

Vector3f calcSQNormal - inputs location and primitive, outputs dS(x,y,z,e,n).

float calcDistance - inputs camera position and intersection ray position in world space,
outputs the distance.

void makeRay - takes a position vector, direction/normal vector, and a ray, and sets
that rays position and normal.

---------------------

Part 2.1:

I create a Phong_Shading function that is called in raytrace.
ISSUES: it doesn't seem to be coloring correctly, overall everything seems right except
that the surfaces facing away from the light source are shaded similar to surfaces
facing towards the light even though they should be diminished. I don't think it's
an issue with the normals, but I think it's either an issue with how I calculate 
material properties and use those to make color, or it might be an issue with how I
handle distance from the light to the point it's color, or might be an issue with
coloring even though another surface of that same primitive is in the way. It also might
be the normals.

raytrace:
I get the current renderable tree, find height and width of frustum using camera,
get the camera position in world space, then run the for loops across the x and y
pixels. Using w, h, x, and y, I find the points at which the camera ray enters the scene,
get the camera directional vector at that point using equations in lecture notes, and 
feed that into findintersection along with running variables. If a primitive is found
intersecting, I run PhongShading on that intersection ray and set the pixel to the
output color from PhongShading. Otherwise, i set the pixel to black. The cameras
directional vector is made from the basis vectors that are converted to world space
after being altered by x and y and near.

New functions:

float Phong_Shading - inputs the intersection ray, the camera ray, the primitive,
the scene, and the current renderable, and outputs the color at that point. In this,
I get position and normal vectors of the intersection ray, position and normal vectors
of the camera ray, set vectors for diffuse, ambient, and specular using prm methods,
then running lighting of Phong Shading on each light. I return the color at the end,
which is a combination of the materials ambient and the sum of diffuse and specular
from the light c wise producted with material diffuse and specular. In the lighting,
I check for shadowing by taking the intersection point minus the light point for the 
vector pointing from light to point, then feed a light ray of the lights position
and lights directional vector into findIntersection. If the resulting prm from this
is different from the original intersecting prm, then don't do lighting calculations,
so it's only colored based on ambient. If it's the same prm, then we do full lighting.

---------------------

Part 2.2:

I successfully implement shadowing in Phong Shading using findIntersection on the light
ray instead of the camera ray to that intersection point.