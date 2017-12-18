### Hw 1 ###

Wireframe Renderer

---------------
Usage

make
./wireframe [scene_description_file.txt] [xres] [yres]

---------------
Input

xres, yres = ints

file format example:
camera:  
position 0 0 5  
orientation 0 1 0 0  
near 1  
far 10  
left -1  
right 1  
top 1  
bottom -1  
 
objects:  
cube cube.obj  
cube2 cube2.obj  
 
cube  
s 1 1 1  
r 0 0 1 0  
t -3 -3 -1  
 
cube2  
s 2 -2 2  
r 0 0 1 0  
t -3 3 -1
...

---------------
Info

I have my main file and separate functions in separate files. I do a lot of functionalities in my main loop as well, which I change and simplify in later assignments by moving more blocks of codes into separate functions and storing those functions in separate files. Ie. I do things like outputting ppm images and certain parsing measures in the main loop, but those can be moved into separate files.

My code:
- Parses scene description file, stores in data structures
- Creates matrices for world to camera and perspective projection transformations
- applies geometric transformations to all object copies, transforms points to NDC
- Maps NDC to screen coordinates in a pixel grid
- Generalizes Bresenham's line algorithm to rasterize lines in pixel rid
- Outputs pixel grid as ppm image

To generalize and implement Bresenham's, I split the algorithm into two functions. 

My first function, which I call Bresenham's, is pretty much a direct implementation on the non-naive algorithm in the lecture notes (algorithm 2). The only difference is that instead of Fill function in the pseudocode, I store that x and y in a point struct and push that struct onto a points list of all the points that will eventually be filled into the grid. The point is in the coordinate dimensions of the screen.

In the second function, called fillGrid, I pass the screen grid by reference and do a bunch of if else statements / switches. Notice that in the first octant, the slope is >= 0 and <= 1, and that the line is calculated by incrementing x and y upwards. From there, I had two options: to generalize to other octants, I could make different functions that decrement instead of increment, or I could use axes to my advantage to flip points such that they have equivalent representations in the first octant. I went the second route; depending on the octant, I identify which axes would need to be flipped to get the mirror of the points from the pointslist in the first octant, then input the start and end points into Bresenham's algorithm having done this flips. For example, first octant is (x,y), and the flip equivalent (so that it has the same slope) in the second octant is (y,x). Or, fourth octant is (-x,y), and so on. Then, Bresenham's gives the pointslist of the rasterized line, which I plug in to the grid to fill the line which will eventually be rendered. I fill the grid by reversing the axes flips on each point so that we get the original line again instead of the mirror.






