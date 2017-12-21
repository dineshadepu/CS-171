#include "Assignment.hpp"

#include "model.hpp"
#include "UI.hpp"
#include "Scene.hpp"
#include "PNGMaker.hpp"
#include "matrixFunctions.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

#include <math.h>
#define _USE_MATH_DEFINES

#define XRES 800
#define YRES 800

using namespace std;

namespace Assignment {
    /* Assignment Part A, Part 1 */
    static const float k_min_io_test_x = -7.0;
    static const float k_max_io_test_x = 7.0;
    static const float k_min_io_test_y = -7.0;
    static const float k_max_io_test_y = 7.0;
    static const float k_min_io_test_z = -7.0;
    static const float k_max_io_test_z = 7.0;
    static const int k_io_test_x_count = 30;
    static const int k_io_test_y_count = 30;
    static const int k_io_test_z_count = 30;

    // return S(x, y, z, e, n) given point and superquadric
    float calcInsideOut(Eigen::Vector4f loc, Primitive *prm)
    {
        float e = prm->getExp0();
        float n = prm->getExp1();
        float x = loc(0);
        float y = loc(1);
        float z = loc(2);

        return powf(powf(x * x, 1.0 / e) + powf(y * y, 1.0 / e), e / n) 
            + powf(z * z, 1.0 / n) - 1.0;
    }

    // returing dS(x, y, z, e, n) as an alternative for prm->getNormal()
    Eigen::Vector3f calcSQNormal(Eigen::Vector4f loc, Primitive *prm)
    {
        float e = prm->getExp0();
        float n = prm->getExp1();
        float x = loc(0);
        float y = loc(1);
        float z = loc(2);
        Eigen::Vector3f dS;

        dS(0) = (x == 0.0) ? 0.0 :
            (2.0 * x * powf(x * x, 1.0 / e - 1.0) * 
            powf(powf(x * x, 1.0 / e) + powf(y * y, 1.0 / e),
                e / n - 1.0)) / n;
        dS(1) = (y == 0.0) ? 0.0 :
            (2.0 * y * powf(y * y, 1.0 / e - 1.0) * 
            powf(powf(x * x, 1.0 / e) + powf(y * y, 1.0 / e),
                e / n - 1.0)) / n;
        dS(2) = (z == 0.0) ? 0.0 :
            (2.0 * z * powf(z * z, 1.0 / n - 1.0)) / n;

        return dS;
    }

    // calculate distance from camera position to a point of intersection
    float calcDistance(Eigen::Vector4f &camera_pos, Eigen::Vector4f &intersection_pos)
    {
        // ALL IN WORLD SPACE
        float x1 = camera_pos(0);
        float y1 = camera_pos(1);
        float z1 = camera_pos(2);
        float x2 = intersection_pos(0);
        float y2 = intersection_pos(1);
        float z2 = intersection_pos(2);

        // return the distance between the two points
        return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
    }

    // get transformation matrix from all transformations applied to that primitive
    Eigen::Matrix4f getTransfMat(Primitive *prm, vector<Transformation> &transformation_stack, bool normal)
    {
        Eigen::Matrix4f m;
        m.setIdentity(4,4);

        for (unsigned int i = 0; i < transformation_stack.size(); i++)
        {
            switch (transformation_stack[i].type)
            {
                case TRANS:
                    if (!normal) // if false, this isnt for normals so include translations
                    {
                        m *= createTransMat(transformation_stack[i].trans(0),
                        transformation_stack[i].trans(1),
                        transformation_stack[i].trans(2));
                    }
                    break;
                case SCALE:
                    m *= createScaleMat(transformation_stack[i].trans(0),
                    transformation_stack[i].trans(1),
                    transformation_stack[i].trans(2));
                    break;
                case ROTATE:
                    m *= createRotMat(transformation_stack[i].trans(0),
                    transformation_stack[i].trans(1),
                    transformation_stack[i].trans(2),
                    transformation_stack[i].trans(3));
                    break;
            }
        }

        // add scaling for coeff
        Eigen::Vector3f coeffVec = prm->getCoeff();
        m *= createScaleMat(coeffVec(0), coeffVec(1), coeffVec(2));

        return m;
    }

    bool IOTest(
        Renderable *ren,
        vector<Transformation> &transformation_stack,
        float x,
        float y,
        float z)
    {
        if (ren->getType() == PRM) {
            Primitive *prm = dynamic_cast<Primitive*>(ren);

            Eigen::Vector4f vec;
            vec << x, y, z, 1.0;

            Eigen::Matrix4f m = getTransfMat(prm, transformation_stack, false);

            // calculate general inside outside fuction of superquadrics in sq space
            return (calcInsideOut(m.inverse() * vec, prm) <= 0) ? true : false;

        } else if (ren->getType() == OBJ) {
            Object *obj = dynamic_cast<Object*>(ren);
            const vector<Transformation>& overall_trans =
                obj->getOverallTransformation();
            for (int i = overall_trans.size() - 1; i >= 0; i--) {
                transformation_stack.push_back(overall_trans.at(i));
            }

            bool IO_result = false;
            for (auto& child_it : obj->getChildren()) {
                const vector<Transformation>& child_trans = 
                    child_it.second.transformations;
                for (int i = child_trans.size() - 1; i >= 0; i--) {
                    transformation_stack.push_back(child_trans.at(i));
                }
                IO_result |= IOTest(
                    Renderable::get(child_it.second.name),
                    transformation_stack,
                    x, y, z);
                transformation_stack.erase(
                    transformation_stack.end() - child_trans.size(),
                    transformation_stack.end());
            }

            transformation_stack.erase(
                transformation_stack.end() - overall_trans.size(),
                transformation_stack.end());
            return IO_result;
        } else {
            fprintf(stderr, "Renderer::draw ERROR invalid RenderableType %d\n",
                ren->getType());
            exit(1);
        }

        return true;
    }

    void drawIOTest() {
        const Line* cur_state = CommandLine::getState();
        Renderable* current_selected_ren = NULL;

        if (cur_state) {
            current_selected_ren = Renderable::get(cur_state->tokens[1]);
        }

        if (current_selected_ren == NULL) {
            return;
        }

        const float IO_test_color[3] = {0.5, 0.0, 0.0};
        glMaterialfv(GL_FRONT, GL_AMBIENT, IO_test_color);
        for (int x = 0; x < k_io_test_x_count; x++) {
            for (int y = 0; y < k_io_test_y_count; y++) {
                for (int z = 0; z < k_io_test_z_count; z++) {
                    float test_x = k_min_io_test_x
                        + x * (k_max_io_test_x - k_min_io_test_x)
                            / (float) k_io_test_x_count;
                    float test_y = k_min_io_test_y
                        + y * (k_max_io_test_y - k_min_io_test_y)
                            / (float) k_io_test_y_count;
                    float test_z = k_min_io_test_z
                        + z * (k_max_io_test_z - k_min_io_test_z)
                            / (float) k_io_test_z_count;

                    vector<Transformation> transformation_stack;
                    if (IOTest(
                            current_selected_ren,
                            transformation_stack,
                            test_x,
                            test_y,
                            test_z))
                    {
                        glPushMatrix();
                        glTranslatef(test_x, test_y, test_z);
                        glutWireSphere(0.05, 4, 4);
                        glPopMatrix();
                    }
                }
            }
        }
    }

    /* Assignment Part A, Part 2 */
    struct Ray {
        float origin_x;
        float origin_y;
        float origin_z;

        float direction_x;
        float direction_y;
        float direction_z;

        Vector3f getLocation(float t) {
            Vector3f loc;
            loc << origin_x + t * direction_x,
                origin_y + t * direction_y,
                origin_z + t * direction_z;
            return loc;
        }
    };

    // make a ray of a point's origin and normal
    void makeRay(Eigen::Vector4f &pos, Eigen::Vector3f &dir, Ray &ray)
    {
        ray.origin_x = pos(0);
        ray.origin_y = pos(1);
        ray.origin_z = pos(2);
        ray.direction_x = dir(0);
        ray.direction_y = dir(1);
        ray.direction_z = dir(2);
    }

    float newtonsMethod(float t_initial, Primitive *prm, 
        Eigen::Vector4f &dir, Eigen::Vector4f &pos)
    {
        float epsilon = 0.0001;
        float t = t_initial;

        // the iterative equation is t_new = t - g(t) / g'(t)
        // for g(t), using S(at+b, e, n) where at+b is the ray
        // for g'(t), using dir dot dS(at+b, e, n)

        // initial calculations for initial check
        float g = calcInsideOut(dir * t + pos, prm);
        float gprime = dir.head(3).dot(calcSQNormal(dir * t + pos, prm));
        if (gprime > epsilon) // initial check to make sure we aren't already past sq
            return -1.0; // bad t

        // cerr << "initial ray point" << dir * t + pos << endl;
        // cerr << "g(t0): " << g << endl;
        // cerr << "g'(t0): " << gprime << endl;

        float newgprime;
        // run the iteration
        while (gprime <= epsilon)
        {
            if (abs(g) < epsilon) // |S(at+b)| ~ 0 so return t
                return t;
            if (abs(gprime) < epsilon) // S is not close to 0 but dS is, sign change case
                return -1.0;
            t = t - g / gprime;
            g = calcInsideOut(dir * t + pos, prm);
            newgprime = dir.head(3).dot(calcSQNormal(dir * t + pos, prm));
            if (newgprime < gprime)
                break;
            gprime = newgprime;
            // cerr << "updated t: " << t << endl;
            // cerr << "g(t): " << g << endl;
            // cerr << "g'(t): " << gprime << endl;
        }

        // cerr << endl;

        // if we exit out the while loop and still havent returned anything
        // that means the sign changed of the derivative, so return bad t
        return -1.0;
    }

    void findIntersection(const Ray &camera_ray, 
        Renderable *ren, 
        vector<Transformation> &transformation_stack,
        float &current_dist,
        Ray &current_ray,
        Primitive ** current_prm) 
    {
        if (ren->getType() == PRM)
        {
            Primitive *prm = dynamic_cast<Primitive*>(ren);

            // get transformation matrices, both with and without translations
            Eigen::Matrix4f m = getTransfMat(prm, transformation_stack, false);
            Eigen::Matrix4f m2 = getTransfMat(prm, transformation_stack, true);

            // world -> sq space inverse transformation of incoming camera ray
            Eigen::Vector4f dir_world, pos_world;
            dir_world << camera_ray.direction_x, camera_ray.direction_y,
                camera_ray.direction_z, 1.0;
            pos_world << camera_ray.origin_x, camera_ray.origin_y,
                camera_ray.origin_z, 1.0;

            Eigen::Vector4f pos = m.inverse() * pos_world;
            Eigen::Vector4f dir = m2.inverse() * dir_world;
            pos /= pos(3); // normalize
            dir /= dir(3);

            // calculate a, b, c for t initial
            float a = dir.head(3).dot(dir.head(3));
            float b = 2 * dir.head(3).dot(pos.head(3));
            float c = pos.head(3).dot(pos.head(3)) - 3.0;

            float determinant = b * b - 4.0 * a * c;

            // check if ray misses bounding sphere altogether using determinant
            if (determinant >= 0.0)
            {
                // calculate initial t's using rays in superquadric space
                float t_plus = (-b + sqrt(determinant)) / (2.0 * a);
                float t_minus = (-b - sqrt(determinant)) / (2.0 * a);

                if (t_plus >= 0.0 && t_minus >= 0.0)
                {
                    // use t minus
                    float t_final = newtonsMethod(t_minus, prm, dir, pos);

                    // t_final positive means we hit sq face towards us, 
                    // so check distance and update prm/ray
                    // t_final negative means we miss, superquadric cant be seen
                    if (t_final >= 0.0)
                    {
                        // find intersection in sq space using t_final, 
                        // transform it to world space
                        Eigen::Vector4f intersection_sq = dir * t_final + pos;
                        intersection_sq(3) = 1.0;
                        Eigen::Vector4f intersection_world = m * intersection_sq;
                        intersection_world /= intersection_world(3); // normalize

                        float dist = calcDistance(pos_world, intersection_world);

                        // check if intersection makes a new shortest distance to 
                        // trace to from camera to object
                        if (dist < current_dist)
                        {
                            // calculate normal in sq space then transform 
                            // to world space and normalize
                            Eigen::Vector3f intersection_normal_sq = 
                                prm->getNormal(intersection_sq.head(3));
                            Eigen::Vector4f temp;
                            temp << intersection_normal_sq(0), 
                                intersection_normal_sq(1), 
                                intersection_normal_sq(2), 
                                1.0;
                            Eigen::Vector4f intersection_normal_world = 
                                m2.inverse().transpose() * temp;
                            intersection_normal_world /= intersection_normal_world(3);
                            Eigen::Vector3f normal = intersection_normal_world.head(3);
                            normal.normalize();

                            // reset the running variables
                            makeRay(intersection_world, normal, current_ray);
                            *current_prm = prm;
                            current_dist = dist;
                        }
                    }
                }
                else if (t_plus >= 0.0 && t_minus < 0.0)
                {   
                    // try both
                    float t_final_plus = newtonsMethod(t_plus, prm, dir, pos);
                    float t_final_minus = newtonsMethod(t_minus, prm, dir, pos);

                    /* COULD JUST USE T_MINUS, DOING BOTH TO BE THOROUGH */

                    // if both values are positive, the superquadric is also 
                    // in bounding sphere
                    if (t_final_plus >= 0.0 && t_final_minus >= 0.0)
                    {
                        // calculate intersecting points, transform sq space -> world space
                        Eigen::Vector4f intersection_sq1 = dir * t_final_plus + pos;
                        Eigen::Vector4f intersection_sq2 = dir * t_final_minus + pos;

                        intersection_sq1(3) = 1.0;
                        intersection_sq2(3) = 1.0;

                        Eigen::Vector4f intersection_world1 = m * intersection_sq1;
                        Eigen::Vector4f intersection_world2 = m * intersection_sq2;

                        intersection_world1 /= intersection_world1(3);
                        intersection_world2 /= intersection_world2(3);

                        // calculate distances of both intersecting points 
                        // to camera and compare
                        float dist1 = calcDistance(pos_world, intersection_world1);
                        float dist2 = calcDistance(pos_world, intersection_world2);

                        // closer to camera point is side of SQ facing the camera
                        float dist = (dist1 <= dist2) ? dist1 : dist2;

                        if (dist < current_dist) // we have a new closer sq
                        {
                            Eigen::Vector3f intersection_normal_sq;
                            Eigen::Vector4f temp, intersection_normal_world;
                            if (dist == dist1)
                            {
                                // calculate normal in sq space, transform it to world space
                                Eigen::Vector3f intersection_normal_sq = 
                                    prm->getNormal(intersection_sq1.head(3));
                                temp << intersection_normal_sq(0), 
                                    intersection_normal_sq(1), 
                                    intersection_normal_sq(2), 
                                    1.0;
                                Eigen::Vector4f intersection_normal_world = 
                                    m2 * temp;
                                intersection_normal_world /= intersection_normal_world(3);
                                Eigen::Vector3f normal = 
                                    intersection_normal_world.head(3);
                                normal.normalize();
                                makeRay(intersection_world1, 
                                    normal, current_ray);
                            }
                            else // dist == dist2
                            {
                                // calculate normal in sq space, transform it to world space
                                Eigen::Vector3f intersection_normal_sq = 
                                    prm->getNormal(intersection_sq2.head(3));
                                temp << intersection_normal_sq(0), 
                                    intersection_normal_sq(1), 
                                    intersection_normal_sq(2), 
                                    1.0;
                                Eigen::Vector4f intersection_normal_world = 
                                    m2 * temp;
                                intersection_normal_world /= intersection_normal_world(3);
                                Eigen::Vector3f normal = 
                                    intersection_normal_world.head(3);
                                normal.normalize();
                                makeRay(intersection_world2, 
                                    normal, current_ray);
                            }
                            // reset the running variables
                            *current_prm = prm;
                            current_dist = dist;
                        }
                    }
                }
            }
        }
        else if (ren->getType() == OBJ)
        {
            Object *obj = dynamic_cast<Object*>(ren);
            const vector<Transformation>& overall_trans =
                obj->getOverallTransformation();
            for (int i = overall_trans.size() - 1; i >= 0; i--) {
                transformation_stack.push_back(overall_trans.at(i));
            }

            for (auto& child_it : obj->getChildren())
            {
                const vector<Transformation>& child_trans = 
                    child_it.second.transformations;
                for (int i = child_trans.size() - 1; i >= 0; i--)
                    transformation_stack.push_back(child_trans.at(i));

                // run the void function on each child of the object
                findIntersection(camera_ray,
                    Renderable::get(child_it.second.name),
                    transformation_stack,
                    current_dist,
                    current_ray,
                    current_prm);

                transformation_stack.erase(
                    transformation_stack.end() - child_trans.size(),
                    transformation_stack.end());
            }

            transformation_stack.erase(
                transformation_stack.end() - overall_trans.size(),
                transformation_stack.end());
        }
        else
        {
            fprintf(stderr, "Renderer::draw ERROR invalid RenderableType %d\n",
                ren->getType());
            exit(1);
        }
    }

    void drawIntersectTest(Camera *camera) {
        // get current renderable in scene
        const Line* cur_state = CommandLine::getState();
        Renderable* current_selected_ren = NULL;

        if (cur_state) {
            current_selected_ren = Renderable::get(cur_state->tokens[1]);
        }

        if (current_selected_ren == NULL) {
            return;
        }

        Ray camera_ray;

        // direction is initial direction (0,0,-1)
        Eigen::Vector3f pos = camera->getPosition(); // b position
        Eigen::Vector4f dir;
        dir << 0.0, 0.0, -1.0, 1.0; // a direction

        // create rotation matrix from rotation axis and angle of camera
        // to transform direction from camera space to world space
        Eigen::Vector3f axis = camera->getAxis();
        Eigen::Matrix4f m = createRotMat(axis(0), axis(1), axis(2), camera->getAngle());
        Eigen::Vector4f newdir = m * dir;
        newdir /= newdir(3); // normalize

        // camera ray set with world space location and direction
        camera_ray.origin_x = pos(0);
        camera_ray.origin_y = pos(1);
        camera_ray.origin_z = pos(2);
        camera_ray.direction_x = newdir(0);
        camera_ray.direction_y = newdir(1);
        camera_ray.direction_z = newdir(2);

        // variables to be fed to findIntersection
        vector<Transformation> transformation_stack;
        float distance = numeric_limits<float>::infinity();
        Primitive *prm = NULL;
        Ray intersection_ray;

        // this function will set the primitive pointer and intersection ray
        // if an intersection is found based on the closest distance
        // otherwise primitive pointer is null and nothing will render
        findIntersection(camera_ray, current_selected_ren, transformation_stack, 
            distance, intersection_ray, &prm);

        if (prm != NULL)
        {
            const float IO_test_color[3] = {1.0, 1.0, 1.0};
            glMaterialfv(GL_FRONT, GL_AMBIENT, IO_test_color);
            glLineWidth(10.0);
            glBegin(GL_LINES);
            glVertex3f(
                intersection_ray.origin_x,
                intersection_ray.origin_y,
                intersection_ray.origin_z);
            Vector3f endpoint = intersection_ray.getLocation(1.0);
            glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
            glEnd();
        }
        else
            cerr << "nothing found" << endl;
        
    }

    /* Assignment Part B */

    /* Implemented phong shading, which colors that pixel */
    float* Phong_Shading(Ray &intersection_ray, Ray &camera_ray, Primitive *prm, Scene &scene,
        Renderable *ren, float color[])
    {
        // make vectors of camera origin and direction, already normalized
        Eigen::Vector3f e, e_direction;
        e << camera_ray.origin_x, camera_ray.origin_y, 
            camera_ray.origin_z; 
        e_direction << camera_ray.direction_x, camera_ray.direction_y, 
            camera_ray.direction_z;

        // make vectors of intersection point position and normal, already normalized
        Eigen::Vector3f P, n;
        P << intersection_ray.origin_x, intersection_ray.origin_y, 
            intersection_ray.origin_z;
        n << intersection_ray.direction_x, intersection_ray.direction_y, 
            intersection_ray.direction_z;

        // set diffuse, ambient, specular and shininess using prm's color * property
        // because each primitive object has a native color and properties 
        // that affect that coloration
        Eigen::Vector3f c_d, c_a, c_s;
        RGBf prmcolor = prm->getColor();
        float diffuse = prm->getDiffuse();
        float ambient = prm->getAmbient();
        float specular = prm->getSpecular();
        c_d << prmcolor.r * diffuse, prmcolor.g * diffuse, prmcolor.b * diffuse;
        c_a << prmcolor.r * ambient, prmcolor.g * ambient, prmcolor.b * ambient;
        c_s << prmcolor.r * specular, prmcolor.g * specular, prmcolor.b * specular;
        float p = prm->getGloss();

        // set diffuse sum and specular sum
        Eigen::Vector3f diffuse_sum, specular_sum;
        diffuse_sum.setZero(3,1);
        specular_sum.setZero(3,1);

        PointLight l;
        Ray light_ray;
        Ray inter_ray;
        vector<Transformation> transformation_stack;
        Primitive *prm2;
        Eigen::Vector3f l_p, l_c, l_direction, l_diffuse, l_specular, temp;
        float distance;
        double ndotl, ndotdir;
        /* Use lights to calculate diffuse and specular sums */
        for (unsigned int i = 0; i < scene.lights.size(); i++)
        {
            l = scene.lights[i];
            l_p << l.position[0], l.position[1], l.position[2];
            l_c << l.color[0], l.color[1], l.color[2];

            // vector pointing from light at the point
            l_direction = P - l_p;
            // now normalize the direction vector from light to point
            l_direction.normalize();

            // attenuate using distance (norm of vector from light to point)
            // using 0.0002 instead of l.k because it's too dark
            // and scaling up a bit
            l_c *= 1.0 / (1.0 + 0.0002 * pow(l_direction.norm(), 2.0));

            // make light ray using l direction and l origin
            // check if light intersects with other object before intersection point
            // by comparing original intersection point with found intersection from this call
            // if different, it hits something before intersection point so shadow it
            light_ray.origin_x = l_p(0);
            light_ray.origin_y = l_p(1);
            light_ray.origin_z = l_p(2);
            light_ray.direction_x = l_direction(0);
            light_ray.direction_y = l_direction(1);
            light_ray.direction_z = l_direction(2);

            transformation_stack.clear();
            distance = numeric_limits<float>::infinity();
            prm2 = NULL;

            // cerr << " Checking in light " << endl;
            findIntersection(light_ray, ren, transformation_stack,
               distance, inter_ray, &prm2);

            if (prm2 == prm)
            {
                ndotl = n.dot(l_direction);
                l_diffuse = l_c * (max(0.0, ndotl));
                diffuse_sum += l_diffuse;

                temp = e_direction - l_direction;
                temp.normalize();
                ndotdir = n.dot(temp);
                l_specular = l_c * pow(max(0.0, ndotdir), p);
                specular_sum += l_specular;
            } 
        }

        Eigen::Vector3f onesvec;
        onesvec.setOnes(3,1);

        // if it's being shadowed, diffuse sum and specular sum are 0, so only ambient affects lighting
        temp = onesvec.cwiseMin(c_a + diffuse_sum.cwiseProduct(c_d) + specular_sum.cwiseProduct(c_s));

        // finish color assignment
        color[0] = temp(0); color[1] = temp(1); color[2] = temp(2);

        return color;
    }

    /* Ray traces the scene. */
    void raytrace(Camera camera, Scene scene)
    {
        // LEAVE THIS UNLESS YOU WANT TO WRITE YOUR OWN OUTPUT FUNCTION
        PNGMaker png = PNGMaker(XRES, YRES);

        const Line* cur_state = CommandLine::getState();
        Renderable* current_selected_ren = NULL;
        if (cur_state)
            current_selected_ren = Renderable::get(cur_state->tokens[1]);
        if (current_selected_ren == NULL)
            return;

        Ray camera_ray;

        // calculate width and height
        // fov is in degrees, convert to radians
        float h = 2.0 * camera.getNear() * tan((camera.getFov() * M_PI / 180.0) / 2.0);
        float w = camera.getAspect() * h;

        // get e1 for camera basis vector that is direction camera is looking
        // and e2 and e3
        // then calculate a at each position we want

        Eigen::Vector4f temp, dir, e1, e2, e3;
        Eigen::Vector3f axis = camera.getAxis();
        Eigen::Matrix4f m = createRotMat(axis(0), axis(1), axis(2), camera.getAngle());

        Eigen::Vector3f pos = camera.getPosition(); // position (b), already in world space
        camera_ray.origin_x = pos(0);
        camera_ray.origin_y = pos(1);
        camera_ray.origin_z = pos(2);

        e1 << 0.0, 0.0, -1.0, 1.0; // direction camera is looking in camera space for basis
        e2 << 1.0, 0.0, 0.0, 1.0; // pointing right relative to camera in camera space for basis
        e3 << 0.0, 1.0, 0.0, 1.0; // pointing up relative to camera in camera space for basis

        // iterate through x by y grid and convert to plane using w and h
        float xi, yi;
        vector<Transformation> transformation_stack;
        float distance;
        Primitive *prm;
        Ray intersection_ray;
        float color[3];
        float *finalcolor;
        for (int i = 0; i < XRES; i ++) {
            for (int j = 0; j < YRES; j++) {
                xi = (i - XRES / 2.0) * w / XRES;
                yi = (j - YRES / 2.0) * h / YRES;

                // create direction vector (a) for at+b
                temp = camera.getNear() * e1 + xi * e2 + yi * e3;
                temp(3) = 1.0;
                dir = m * temp; // transform to world space
                dir /= dir(3); // normalize

                camera_ray.direction_x = dir(0);
                camera_ray.direction_y = dir(1);
                camera_ray.direction_z = dir(2);

                // reset variables and feed camera ray + vars into find intersection function
                transformation_stack.clear();
                distance = numeric_limits<float>::infinity();
                prm = NULL;

                findIntersection(camera_ray, current_selected_ren, transformation_stack,
                    distance, intersection_ray, &prm);

                if (prm != NULL)
                {                        
                    // do phong shading at the found intersection
                    // using the lights and the normal of the intersection ray and prm
                    finalcolor = Phong_Shading(intersection_ray, camera_ray, prm, scene,
                        current_selected_ren, color);
                    png.setPixel(i, j, finalcolor[0], finalcolor[1], finalcolor[2]);
                }
                else // hit no object so color black
                    png.setPixel(i, j, 0.0, 0.0, 0.0);

                cerr << '\r' << "x: " << i << "  y: " << j << flush;
            }
        }
        cerr << endl;

        // LEAVE THIS UNLESS YOU WANT TO WRITE YOUR OWN OUTPUT FUNCTION
        if (png.saveImage())
            fprintf(stderr, "Error: couldn't save PNG image\n");
        else
            printf("DONE!\n");
    }
};
