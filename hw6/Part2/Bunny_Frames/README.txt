### Part 2.2 ###

Interpolating a Smoothing Bunny

---------------
Usage

make
./keyframe [xres] [yres]

---------------
Input

xres, yres = integers for opengl display window resolution

---------------
Info

This program interpolates across the keyframes given for the bunny meshes (from unsmoothed to very smoothed) to generate 16 missing frames and create a smooth transition from frame 0 to frame 20 and back.

Additional files to bunnysmoother.cpp (which contains the main loop): halfedge.h, 
structs.h. In bunnysmoother.cpp, I used the opengl program for the double pendulum
as a template and added onto it from hw 5 (implicit fairing) code.

Overall, my program both dumps all my interpolated values into obj files to compare with 
the correct values, and also displays the bunny with the ability to go up/down frames
and display it smoothly. This also has quaternion, so the scene can be rotated and still
render accordingly.

Run at 800 x 800. The reason I don't include the object file names in the command line is
that the executable program from this code should be run from the folder that contains 
subfolders "keyframes" and "interpolatedframes" for the bunny frames. There is also the 
"myframes" folder, in which I write the object files of all 20 frames, to verify that
it is storing the keyframe ones correctly and also storing the interpolated ones that
I calculated. Comparing myframes object files with the interpolatedframes ones after,
running, I found that all my frames are the same to floating point precision, so I
interpolated correctly.

My code has the keyframes destination hard coded in, so I go through all object files
in that folder first and store them in meshes (in the function storeBunnies()). The 
global meshes vector will hold all meshes for all frames. After storing, I interpolate
across all vertices and all dimensions, and store all the calc'd vertices in the
appropriate vectors/frames in the meshes (in the function interpolateVals()). After that,
I go through each frame mesh and calculate all the normals from the vertices/faces (I only
have one main face vector because each object has the same face structure in the same
order), and store those normals for each frame (in the function halfedgeCalcs()). After,
I go through all frame meshes and store their vertices and faces in object files which
I dump into a folder called "myframes". After that, I display.

Pressing key "i" will increment the current frame and redisplay the new object
by storing that frames' vertices and normals into the buffer and then using gl
to display it. It will cap at frame 20 (does not wrap around). Pressing key "o" will
decrement the current frame, and has a lower limit at 0.

Since this is fairly hard coded for the given bunny keyframes and task, I allocate
everything and calculate everything initially before rendering and cycling through
the frames because the memory load isn't prohibitive. If we had a much greater number
of frames to work with, and many more objects, then I would compute on the fly
and not store as such.

A lot of my logic and structure is in the comments as well.