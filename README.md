# CS-171
## Caltech Computer Graphics Lab with Professor Alan Barr

Final Grade received: A

You will need to use a unix-based operating system to run this code, such as Mac OS or Linux. I have a Mac OS but I set up and ran a Ubuntu 14.04 Linux VM since Linux allows easy installation of required versions of OpenGL and ImageMagick.

### Installing g++ on Ubuntu Linux:
sudo apt-get update
sudo apt-get install build-essential

### Installing the older, stabler version of OpenGL used for this class on Ubuntu Linux:
sudo apt-get install freeglut3
sudo apt-get install freeglut3-dev
sudo apt-get install libglew1.5
sudo apt-get install libglew1.5-dev
sudo apt-get install libglu1-mesa
sudo apt-get install libglu1-mesa-dev
sudo apt-get install libgl1-mesa-dev

For Mac and other distributions of Linux, you will need to find the installation processes online.

### You can test your OpenGL installation using the OpenGL demo.
Inside the OpenGL_demo folder, run Make and execute ./opengl_demo. You can then test the following:

- Clicking in the window with the left-mouse button and dragging the mouse around rotates your view of the scene in the direction you drag. Note that you cannot rotate your view up or down vertically more than 90 degrees. This mouse interaction scheme is supposed to mimic a real-life first-person-view of a scene, and of course, in real life, our heads cannot rotate vertically more than about 90 degrees.
- Pressing the w, s, a, and d keys on your keyboard moves your view of the scene forward, backward, to the right, and to the left respectively. As you move away from the cubes, parts of the cubes will disappear into the black distance. If you move far enough, you may go off the blue ground and find yourself floating in black space. You can also move through the cubes.
- Pressing the q key quits the program.
- Pressing the t key on your keyboard toggles the cubes back and forth from being solid to “wireframes” as shown below:
