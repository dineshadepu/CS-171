### Hw 3 ###

Introduction to OpenGL and Arcball implementation

---------------
Usage

make
./opengl_renderer [data/scene.txt] [xres] [yres]

---------------
Input

xres, yres = integers for opengl display window resolution
scene file is in format seen before

---------------
Info

In this program I implement the shaded surface renderer in OpenGL and implemented the ability to rotate scenes using arcball (with quaternions for rotations).

My code:
- Parses scene file and puts data into structures
- Creates shaded surface renderer within an opengl context
- Allows arcball scene manipulation based on quaternion method

Included files:
	src folder:
	parsing.cpp, parsing.h - for parsing the txt file
	structs.h - to hold structs, both given from demo and created
	data folder - holds a bunch of txt files and objects

I normalize the rotation matrices by taking the norm of the u vector and dividing all the elements by that norm, and use this normalized u vector when create the unit quaternion instead of dividing the quaternion by its norm (original * conjugate).

I also decided to keep all of the major functions in the main file because of how the global variables are used. In addition to the functions that came with the demo (some which I modified based on the lecture notes / algorithms or based on how I parsed and stored data structures differently), I added the following functions:
	Screen to NDC
	Computing Rotation Matrix (axis angle and quaternion)
	Get current rotation ( current * last)
	Get rotation matrix from current quaternion
	Get current rotation quaternion (current quat * last quat)

Comments and potential issues: not sure how big the sphere is, because when dragging and rotating from the corners, it doesn’t seem to properly be rotating around the z axis (at least it doesn’t seem like it is). Also unsure if how I calculated unit quaternions will affect the computations negatively.
	