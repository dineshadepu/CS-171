### Part 3 ###
A program in C++ that takes as input (from the command line) a single text file; loads all the data into appropriate data structures; applies all the specified transformations to their respective objects; and outputs, for each object, the transformed vertices to standard output.

# Structure of input txt file
object1 object1_filename.obj  
object2 object2_filename.obj  
object3 object3_filename.obj  
 
object1  
t tx ty tz  
r rx ry rz angle_in_radians  
 
object1  
t tx ty tz  
s sx sy sz  
 
object2  
r rx ry rz angle_in_radians  
s sx sy sz  
 
object3  
t tx ty tz  
r rx ry rz angle_in_radians  
s sx sy sz

# Usage:
make
./my_prog file.txt

# Output:
Vertex coordinates are transformed coordinates.

object1_copy1  
x1 y1 z1  
...  
xm ym zm  
 
object1_copy2  
x1 y1 z1  
...  
xm ym zm  
 
object2_copy1  
x1 y1 z1  
...  
xm ym zm  
 
object3_copy1  
x1 y1 z1  
...  
xm ym zm

