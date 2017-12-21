/* CS/CNS 171
 * Fall 2017
 * Written by Eshan Govil
 *
 * The following program animates the I-bar
 */

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <math.h>
#define _USE_MATH_DEFINES

#include <cmath>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

using namespace std;

struct Transformation
{
    float x, y, z;
};

struct Quaternion
{
    float qs, qx, qy, qz;
    /* Using current rotation quaternion, make rotation matrix */
    Eigen::Matrix4d Get_Rotation_Matrix_Quat()
    {
        Eigen::Matrix4d m;
        m << 1-2*qy*qy-2*qz*qz, 2*(qx*qy-qz*qs),
            2*(qx*qz+qy*qs), 0,
            2*(qx*qy+qz*qs), 1-2*qx*qx-2*qz*qz,
            2*(qy*qz-qx*qs), 0,
            2*(qx*qz-qy*qs), 2*(qy*qz+qx*qs),
            1-2*qx*qx-2*qy*qy, 0, 
            0, 0, 0, 1;
        return m;
    }
};

const float cam_position[] = {0, 0, 40};

const float near_param = 1.0, far_param = 60.0,
            left_param = -1.0, right_param = 1.0,
            top_param = 1.0, bottom_param = -1.0;

const float fstep = 1; // frame step
const float t_val = 0; // for Catmull-Rom splines
float currFrame = 0; // current frame that will be updated and used to render
float totFrames; // total frames in script

vector<int> frameVals; // indices of keyframes
vector<Transformation> translations; // vector of all interpolated translations
vector<Transformation> scalings; // vector of all interpolated scalings
vector<Quaternion> rotations; // vector of all interpolated rotations

GLUquadricObj *quadratic;

void applyFrameTransf();
void drawIBar();
void draw_text();
void parse_script(string &line, bool &b, int &frameNum);

void init(void)
{
    // don't initialize lights or material
    // because Ibar renders colored material for you

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    glFrustum(left_param, right_param,
            bottom_param, top_param,
            near_param, far_param);
    
    glMatrixMode(GL_MODELVIEW);

    quadratic = gluNewQuadric();
}

void reshape(int width, int height)
{
    height = (height == 0) ? 1 : height;
    width = (width == 0) ? 1 : width;
    
    glViewport(0, 0, width, height);
    
    glutPostRedisplay();
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glTranslatef(-cam_position[0], -cam_position[1], -cam_position[2]);

    // first doing keyframe scene transformations to whole scene
    // based on current frame
    // then doing drawing and transforming each cylinder
    glPushMatrix();
    applyFrameTransf();
    drawIBar();
    glPopMatrix();

    draw_text();
    
    glutSwapBuffers();
}

// Forward incrementing the frame with wrap around
void updateF()
{
    // /* FRAME UPDATE RULES */
    currFrame += fstep;
    if (currFrame == totFrames)
        currFrame = 0;
    glutPostRedisplay();
}

// Backward incrementing the frame with wrap around
void updateB()
{
    // /* FRAME UPDATE RULES */
    currFrame -= fstep;
    if (currFrame == -1)
        currFrame = totFrames-1;
    glutPostRedisplay();
}

void applyFrameTransf()
{
    // these transformations are in order as given in script
    glTranslatef(translations[currFrame].x, translations[currFrame].y, 
        translations[currFrame].z);
    glScalef(scalings[currFrame].x, scalings[currFrame].y, 
        scalings[currFrame].z);
    glMultMatrixd(rotations[currFrame].Get_Rotation_Matrix_Quat().data());
}

void drawIBar()
{
    /* Parameters for drawing the cylinders */
    float cyRad = 0.2, cyHeight = 1.0;
    int quadStacks = 4, quadSlices = 4;

    glPushMatrix();
    glColor3f(0, 0, 1);
    glTranslatef(0, cyHeight, 0);
    glRotatef(90, 1, 0, 0);
    gluCylinder(quadratic, cyRad, cyRad, 2.0 * cyHeight, quadSlices, quadStacks);
    glPopMatrix();
    
    glPushMatrix();
    glColor3f(0, 1, 1);
    glTranslatef(0, cyHeight, 0);
    glRotatef(90, 0, 1, 0);
    gluCylinder(quadratic, cyRad, cyRad, cyHeight, quadSlices, quadStacks);
    glPopMatrix();
    
    glPushMatrix();
    glColor3f(1, 0, 1);
    glTranslatef(0, cyHeight, 0);
    glRotatef(-90, 0, 1, 0);
    gluCylinder(quadratic, cyRad, cyRad, cyHeight, quadSlices, quadStacks);
    glPopMatrix();
    
    glPushMatrix();
    glColor3f(1, 1, 0);
    glTranslatef(0, -cyHeight, 0);
    glRotatef(-90, 0, 1, 0);
    gluCylinder(quadratic, cyRad, cyRad, cyHeight, quadSlices, quadStacks);
    glPopMatrix();
    
    glPushMatrix();
    glColor3f(0, 1, 0);
    glTranslatef(0, -cyHeight, 0);
    glRotatef(90, 0, 1, 0);
    gluCylinder(quadratic, cyRad, cyRad, cyHeight, quadSlices, quadStacks);
    glPopMatrix();

}

template<typename T>
string tostr(const T& t)
{
    ostringstream os;
    os << t;
    return os.str();
}

void draw_text()
{
    glColor3f(0, 1, 0);

    string frame_str = "Frame: " + tostr(currFrame);

    glRasterPos2f(25,35);
    for(int i = 0; i < frame_str.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, frame_str[i]);
}

void key_pressed(unsigned char key, int x, int y)
{
    if(key == 'q')
    {
        exit(0);
    }
    else if (key == 'i')
    {
        updateF();
    }
    else if (key == 'o')
    {
        updateB();
    }
}

void parse_script(string &line, bool &b, int &frameNum)
{
    istringstream iss(line);
    string temp;

    if (b == false)
    {
        iss >> totFrames;

        rotations.resize(totFrames);
        translations.resize(totFrames);
        scalings.resize(totFrames);

        b = true;
    }
    else if (line.find("Frame") != string::npos)
    {
        iss >> temp >> frameNum;
        frameVals.push_back(frameNum);
    }
    else if (line.find("translation") != string::npos)
    {
        iss >> temp >> translations[frameNum].x
            >> translations[frameNum].y
            >> translations[frameNum].z;
    }
    else if (line.find("scale") != string::npos)
    {
        iss >> temp >> scalings[frameNum].x
            >> scalings[frameNum].y
            >> scalings[frameNum].z;
    }
    else if (line.find("rotation") != string::npos)
    {
        float x, y, z, qs, qx, qy, qz, theta, n;
        iss >> temp >> x >> y >> z >> theta;

        // for the rotation matrix
        // convert input into a quaternion and store

        // first, normalize the rotation axis
        n = sqrt(x*x + y*y + z*z);
        x = x/n;
        y = y/n;
        z = z/n;

        // set quaternion values using x, y, z, theta
        qs = cos((theta * M_PI / 180.0) / 2.0);
        qx = x*sin((theta * M_PI / 180.0) / 2.0);
        qy = y*sin((theta * M_PI / 180.0) / 2.0);
        qz = z*sin((theta * M_PI / 180.0) / 2.0);

        // normalize the quaternions again, just in case. Not necessary
        float normer = sqrt(qs*qs + qx*qx + qy*qy + qz*qz);
        rotations[frameNum].qs = qs/normer;
        rotations[frameNum].qx = qx/normer;
        rotations[frameNum].qy = qy/normer;
        rotations[frameNum].qz = qz/normer;
    }
}

void interpolateVals()
{
    // make B matrix for t = 0
    float s = 0.5*(1 - t_val);
    Eigen::Matrix4d B;
    B << 0, 1, 0, 0,
         -s, 0, s, 0,
         2*s, s-3, 3-2*s, -s,
         -s, 2-s, s-2, s;    

    for (int i = 0; i < frameVals.size(); i++)
    {
        // vectors of point vectors
        Eigen::Vector4d transx, transy, transz,
            scalx, scaly, scalz, rotqs,
            rotqx, rotqy, rotqz;

        // indices for primary keyframe positions
        int iminus1, icurrent, iplus1, iplus2;
        icurrent = frameVals[i];

        // handle index assignment for wrap around cases before general case
        if (i == 0)
        {
            iminus1 = frameVals[frameVals.size()-1];
            iplus1 = frameVals[i+1];
            iplus2 = frameVals[i+2];
        }
        else if (i == frameVals.size()-2)
        {
            iminus1 = frameVals[i-1];
            iplus1 = frameVals[i+1];
            iplus2 = frameVals[0];
        }
        else if (i == frameVals.size()-1)
        {
            iminus1 = frameVals[i-1];
            iplus1 = frameVals[0];
            iplus2 = frameVals[1];
        }
        else
        {
            iminus1 = frameVals[i-1];
            iplus1 = frameVals[i+1];
            iplus2 = frameVals[i+2];
        }

        // now store values into point vectors of four control points
        // for all translation, rotation, scaling components
        transx << translations[iminus1].x, translations[icurrent].x,
            translations[iplus1].x, translations[iplus2].x;
        transy << translations[iminus1].y, translations[icurrent].y,
            translations[iplus1].y, translations[iplus2].y;
        transz << translations[iminus1].z, translations[icurrent].z,
            translations[iplus1].z, translations[iplus2].z;

        scalx << scalings[iminus1].x, scalings[icurrent].x,
            scalings[iplus1].x, scalings[iplus2].x;
        scaly << scalings[iminus1].y, scalings[icurrent].y,
            scalings[iplus1].y, scalings[iplus2].y;
        scalz << scalings[iminus1].z, scalings[icurrent].z,
            scalings[iplus1].z, scalings[iplus2].z;

        rotqs << rotations[iminus1].qs, rotations[icurrent].qs,
            rotations[iplus1].qs, rotations[iplus2].qs;
        rotqx << rotations[iminus1].qx, rotations[icurrent].qx,
            rotations[iplus1].qx, rotations[iplus2].qx;
        rotqy << rotations[iminus1].qy, rotations[icurrent].qy,
            rotations[iplus1].qy, rotations[iplus2].qy;
        rotqz << rotations[iminus1].qz, rotations[icurrent].qz,
            rotations[iplus1].qz, rotations[iplus2].qz;

        // get current frame to start at and end frame
        int start, end;
        float du, u, normer, qs, qx, qy, qz;

        start = icurrent;
        if (i == frameVals.size()-1)
            end = totFrames;
        else
            end = iplus1;

        // to go from 0 to 1 across all frames between start and end
        // move in increments of du = 1/(end-start)
        du = 1.0 / (end-start);
        u = du;

        // interpolate all values in between start and end
        for (int j = start+1; j < end; j++)
        {
            Eigen::Vector4d uvec, umultB;
            
            // make u vector and start curve function calculation
            uvec << 1, u, u*u, u*u*u;
            umultB = uvec.transpose() * B;

            // finish the cardinal curve function using dot product and store f in transf vectors

            translations[j].x = umultB.dot(transx);
            translations[j].y = umultB.dot(transy);
            translations[j].z = umultB.dot(transz);

            scalings[j].x = umultB.dot(scalx);
            scalings[j].y = umultB.dot(scaly);
            scalings[j].z = umultB.dot(scalz);

            qs = umultB.dot(rotqs);
            qx = umultB.dot(rotqx);
            qy = umultB.dot(rotqy);
            qz = umultB.dot(rotqz);

            //after calculating the interpolated quaternion, normalize it
            normer = sqrt(qs*qs + qx*qx + qy*qy + qz*qz);
            rotations[j].qs = qs/normer;
            rotations[j].qx = qx/normer;
            rotations[j].qy = qy/normer;
            rotations[j].qz = qz/normer;

            u += du;
        }
    }
}

int main(int argc, char* argv[])
{
    if(argc != 4)
    {
        cerr << "\nERROR: Incorrect number of arguments." << endl;
        exit(1);
    }
    
    fstream myfile(argv[1]);
    int xres = atoi(argv[2]);
    int yres = atoi(argv[3]);

    string line;
    bool foundTotFrames = false;
    int frameNum;

    // Parse script and store in vectors
    cerr << "Loading script" << endl;
    if (myfile.is_open())
    {
        while (!myfile.eof())
        {
            getline(myfile, line);
            parse_script(line, foundTotFrames, frameNum);   
        }
    }
    myfile.close();

    cerr << "Now interpolating values" << endl;
    interpolateVals();
    cerr << "Done interpolating values. Rendering" << endl;
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(xres, yres);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("I Bar Display");
    
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(key_pressed);
    glutMainLoop();
}
