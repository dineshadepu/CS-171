### Part 2 ###

A program that takes as input (from the command line) a single text file that contains a list of translation, rotation, and scaling vectors; creates the corresponding translation, rotation, and scaling matrices for the specified vectors; and outputs to standard output the inverse of the product of all the matrices.

# Structure of input file
t tx ty tz
r rx ry rz angle_in_radians
s sx sy sz

where a line beginning with t indicates a translation vector, a line beginning with r indicates a rotation vector, and a line beginning with s indicates a scaling vector.

# Usage:
make
./my_prog transform_data.txt

# Output:
Eigen matrix printed into std output.