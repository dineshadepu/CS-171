### Hw 5 ###

Introduction to Geometry Processing (Implicit fairing)

---------------
Usage

make
./smooth [scene_file.txt] [xres] [yres] [h]

---------------
Input

xres, yres = integers for opengl display window resolution
h = float time step
scene description file is like we've seen in the past

---------------
Info

This program takes a scene file, resolutions, and the time step, and displays the mesh with a smoothing operation at a key command (if you press "i"). This builds upon code from Hw 3. It first renders the scene in OpenGL but without given normals, so I compute the normal vectors for each point using the halfedge data structure. I think implement a sparse matrix operator that uses Discrete Laplacian and the time step to smooth the mesh (input all vertices and output the new vertex locations after smoothing in each xyz direction).

The program does not render smoothings of various time steps, but rather reapplies the h passed in the command line with each iteration. That is, once my mesh updates, I further smooth on that mesh, and not on the original mesh. This may be the cause of a "bug" in OpenGL, where as the area becomes extremely small at various points, we end up with a column of zeroes in the sparse matrix and the resultant solver becomes singular and unable to solve so OpenGL/Eigen throws an exception and dumps the core. You'll notice this after running enough times on either bunny or armadillo at any large enough h. I could fix this if I were to reimplement such that I smoothed on the original mesh each time and added h or multiplied the time step or whatever, but chose not to do that because A) my program works as is B) Brian said don't sweat it :D and C) that would take a decent amount of modification from my current code which isn't entirely necessary. Either way, I can smooth correctly given pretty much any time step at least 1 time if not several times.

Also, my mesh is rendered on a black background and looks like the normals might be off or that the rendering is facing issues in comparison to how the mesh looks on the course site (grey background), but Brian's looks and behaves the same so I'm not sure if this is an issue.

On top of the dynamically allocated HEVs and HEFs, which I make sure to delete, I also dynamically allocate the meshes vector and all the meshes of the potentially numerous copies. I have created a function that deletes the meshes upon pressing "q" and quitting the program, so I should not have any memory leaks to my knowledge.

New functions and what they do:
	calc_vertex_normal - calculates normal for vertex given adjacent vertices and faces.
	halfedgeCalcs - do normal and operator computations. More details in comments in function and function itself.
	build_F_operator - builds operator F = I-hdx and returns it. The logic for it can be found in the comments within the function.
	solve - uses operator to solve for axis vectors of each point
	calcArea/Alpha/Beta
	storeBuffers - takes meshes with updated vertices/normals and stores into buffer to be fed to OpenGL
	deleteMeshes - memory leaks are a bad thing indeed

I only make one modification to halfedge.h, which is to assign the index of the vertex to hev->index.

I made edits to my parsing files/functions to handle the new meshes.

The boolean initial_toggle is such that I can change the usage of the function halfedgeCalcs() so that it handles smoothing instead of initial normal calculation and dispaly.


	
