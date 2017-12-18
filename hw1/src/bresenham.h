#include <vector>
#include "structs.h"

/* Compute Bresenham's algorithm and line rasterize by returning all the points
it fills along the line instead of actually filling in the grid */
void Bresenham(int x0, int y0, int x1, int y1, std::vector<Point> &pointsList);

/* Use Bresenham's algorithm to get the points to fill in the grid. The if else 
switches determine which octant the given points are in based on the slope
of the line and whether or not dy is positive. Grid is filled because it is
passed by reference */
void fillGrid(int x0, int y0, int x1, int y1, std::vector<std::vector<int>> &grid);