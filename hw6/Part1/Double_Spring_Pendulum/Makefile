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
###############################################################################
CC = g++
FLAGS = -g

INCLUDE = -I/usr/X11R6/include -I/usr/include/GL -I/usr/include
LIBDIR = -L/usr/X11R6/lib -L/usr/local/lib
SOURCES = *.cpp
LIBS = -lGLEW -lGL -lGLU -lglut -lm

EXENAME = double_pendulum

all: $(SOURCES)
	$(CC) $(FLAGS) -o $(EXENAME) $(INCLUDE) $(LIBDIR) $(SOURCES) $(LIBS)

clean:
	rm -f *.o $(EXENAME)

.PHONY: all clean

