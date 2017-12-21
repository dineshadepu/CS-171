### Part 1 ###

Single Spring Pendulum, Double Spring Pendulum, and Elasticity Demo

---------------
Usage

make
1.1: ./single_pendulum [xres] [yres] [x_start] [y_start]
1.2: ./double_pendulum [xres] [yres] [x_start_1] [y_start_1] [x_start_2] [y_start_2]
1.3: ./simulate man.obj

---------------
Input

xres, yres = integers for opengl display window resolution
x_start, y_start = starting positions for pendulum mass(es)

---------------
Info

These programs use discrete Lagrangians and DEL equations and ODEs to calculate the overall system energy at every time step and render pendulum/elasticity animations that follow realistic physical principles. As a note, the vertices are not meant to be dragged a large distance in the elastic demo.

For each of Double_Spring_Pendulum, Elasticity, and Single_Spring_Pendulum, I do majority
of my Lagrangian calculations (deriving the discrete analog and isolating for unknown
variable after setting up the logic such as finding the discrete analog or putting
together the entire analog with and without individual components of q) using 
Mathematica and on paper.

In terms of logic for each, I first calculate k+1 for each of x and y and store in 
temporary variables, then calculate k+1 for all px and py and update originals, then
update original x and y with the temporary variables.