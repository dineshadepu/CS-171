### Part 2.1 ###

Flying I-Bar animation

---------------
Usage

make
./keyframe [script file] [xres] [yres]

---------------
Input

xres, yres = integers for opengl display window resolution
script file = test script file that gives I bar frames

---------------
Info


Main loop in keyframes.cpp. There are no other files to worry about. My keyframes.cpp
is based off of the double pendulum opengl program as a template, with additional code
added from hw 3 (rendering with arcball even though I removed arcball).

I add two structures: Transformation (translation and scaling) and Quaternion
(rotation) which I declare at the top of the file. Got these from Hw 3.

Important global vars: vectors to hold all interpolated and keyframe values.

Before rendering and running the main loop, I initialize and do all my interpolation
calculations and store them in vectors for translations, scalings and quaternions that
are the length of the number of frames in the scene (specified in the first line).
Before that, I parsed the script file and stored the keyframe transformations.
I also correctly normalize and set my quaternions as well as translations and 
scalings during my interpolation function.

I use key 'i' to increment one frame, and use key 'o' to decrement one frame
since we have all the frame data stored in vectors and aren't calculating on the fly,
so it's easy to access values in the vectors and use them to transform the scene
based on the frame number. I know that storing all the frame data in vectors to start
instead of calculating interpolated values on the fly can get memory expensive,
but it's less computationally expensive when running the loops several times and it
was intuitive for me to interpolate and store everything at first, but I can see how to
compute it on the fly. Every time we increment, I redisplay, in which I call
applyFrameTransf() to the scene that pulls that transformations from the stored vectors
using current Frame as the index, before rendering and redisplaying the cylinders
that make up the I bar.

Press 'q' to quit. Nothing is dynamically allocated so there is no need to have a 
delete function. A lot of my logic can be found in the comments in my code.