###############################################################################
# CS/CNS 171 Fall 2015
#
# This is a template Makefile for OpenGL programs. Edit it however you find
# convenient.
# 
# The current version of this file should compile OpenGL programs just fine on
# Debian-based Linux operating systems.
#
# If you run Mac OS or other distributions of Linux, then you may have to
# fool around with the 'INCLUDE' and 'LIBDIR' lines below before the Makefile
# can compile the OpenGL parts successfully.
#
# The C++ matrix library, Eigen, is assumed to be in the directory right
# above the directory that contains this Makefile.
###############################################################################
CC = g++
FLAGS = -g -std=c++11

INCLUDE = -I/usr/X11R6/include -I/usr/include/GL -I/usr/include -I../
LIBDIR = -L/usr/X11R6/lib -L/usr/local/lib
SOURCES = *.h *.cpp
LIBS = -lGLEW -lGL -lGLU -lglut -lm

EXENAME = simulate

all: $(SOURCES)
	$(CC) $(FLAGS) -o $(EXENAME) $(INCLUDE) $(LIBDIR) $(SOURCES) $(LIBS)

clean:
	rm -f *.o $(EXENAME)

.PHONY: all clean

