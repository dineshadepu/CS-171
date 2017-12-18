#include "bresenham.h"
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
#include <cmath>

using namespace std;

/* Compute Bresenham's algorithm and line rasterize by returning all the points
it fills along the line instead of actually filling in the grid */
void Bresenham(int x0, int y0, int x1, int y1, vector<Point> &pointsList)
{
	int err = 0;
	int y = y0;
	int dx = x1 - x0;
	int dy = y1 - y0;
	for (int x = x0; x <= x1; x++)
	{
		Point p; // Use point struct to store coordinate
		p.x = x;
		p.y = y;
		pointsList.push_back(p);
		if (2*(err + dy) < dx)
			err += dy;
		else
		{
			err += dy - dx;
			y++;
		}
	}
}

/* Use Bresenham's algorithm to get the points to fill in the grid. The if else 
switches determine which octant the given points are in based on the slope
of the line and whether or not dy is positive. Grid is filled because it is
passed by reference */
void fillGrid(int x0, int y0, int x1, int y1, vector<vector<int>> &grid)
{
	vector<Point> points;
	int dy = y1-y0;
	int dx = x1-x0;

	if (dx != 0) // not a vertical line so slope is not indefinite
	{
		float m = float(dy)/float(dx);
		{
			if (m >=0 && m <= 1) // first octant
			{
				if (dy >= 0) // first octant
					Bresenham(x0, y0, x1, y1, points);
				else // dy is neg so actually fifth octant, flip points to positive
					// because order you render points in doesn't matter
					Bresenham(x1, y1, x0, y0, points);
				for (int i = 0; i<points.size(); i++)
					grid[points[i].y][points[i].x] = 1;
			}
			else if (m > 1) // second octant
			{
				if (dy >= 0)
					// flip over xy axis so x and y switch places
					Bresenham(y0, x0, y1, x1, points);
				else // flip points then flip xy axis, sixth octant
					Bresenham(y1, x1, y0, x0, points);
				for (int i = 0; i<points.size(); i++)
					grid[points[i].x][points[i].y] = 1;
			}
			else if (m < -1) // third octant
			{
				if (dy >= 0)
					// flip y axis then xy axis
					Bresenham(y0, -x0, y1, -x1, points);
				else // flip points, flip y axis, flip xy axis
					// seventh octant
					Bresenham(y1, -x1, y0, -x0, points);
				for (int i = 0; i<points.size(); i++)
					grid[points[i].x][-points[i].y] = 1;
			}
			else if (m >= -1 && m < 0) // fourth octant
			{
				if (dy >= 0)
					// flip y axis
					Bresenham(-x0, y0, -x1, y1, points);
				else // flip points then flip y axis
					// eighth octant
					Bresenham(-x1, y1, -x0, y0, points);
				for (int i = 0; i<points.size(); i++)
					grid[points[i].y][-points[i].x] = 1;
			}
		}
	}
	else // dx == 0 so vertical line
	{
		if (dy >= 0) // treat like second octant
			Bresenham(y0, x0, y1, x1, points);
		else // flip then treat like second octant
			Bresenham(y1, x1, y0, x0, points);
		for (int i = 0; i<points.size(); i++)
			grid[points[i].x][points[i].y] = 1;
	}
}