### Part 1
A program that takes as input (from the command line) one or more .obj files; stores the vertex and face data from each file into workable C++ data structures; and outputs the data from each file to standard output from the data structures.

#### Structure of .obj file
v x1 y1 z1  
v x2 y2 z2  
v x3 y3 z3  
...  
v xm ym zm  
f face0v1 face0v2 face0v3  
f face1v1 face1v2 face1v3  
f face2v1 face2v2 face2v3  
...  
f facenv1 facenv2 facenv3

Vertex data is xyz 3D float coordinates, face data is 3 integers specifying three vertices that make up a triangular face.

#### Usage:
make
./my_prog obj_file1.obj obj_file2.obj ... obj_fileN.obj

#### Output:
(printing data as stored in data structures)

obj_file1:  
 
v x1 y1 z1  
...  
v xm ym zm  
f face0v1 face0v2 face0v3  
...  
f facenv1 facenv2 facenv3  
 
obj_file2:  
 
v x1 y1 z1  
...  
v xm ym zm  
f face0v1 face0v2 face0v3  
...  
f facenv1 facenv2 facenv3  
 
...  
 
obj_fileN:  
 
v x1 y1 z1  
...  
v xm ym zm  
f face0v1 face0v2 face0v3  
...  
f facenv1 facenv2 facenv3
