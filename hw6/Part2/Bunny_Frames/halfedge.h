/* By Kevin (Kevli) Li (Class of 2016)
 * 
 * This header file contains a quick and dirty implementation of the Halfedge
 * data structure for the purpose of processing meshes WITHOUT boundary.
 *
 * Every part of the implementation is defined loosely as structs or standalone
 * functions, which makes this header file easy to use and portable. Of course,
 * this also means that this implementation is not as robust or has as many
 * safe-checks as professional software.
 *
 * This halfedge implementation was meant to be USED, not edited. Hence, there
 * are only comments here detailing how to use the halfedge, but none on the
 * actual code that implements the data structure.
 *
 * The main function of interest in this header file is:
 *
 *     build_HE(Mesh_Data *mesh,
 *              std::vector<HEV*> *hevs,
 *              std::vector<HEF*> *hefs);
 * 
 * This function constructs the halfedge data structure and takes in three
 * arguments:
 *
 *     Mesh_Data *mesh is a reference to a Mesh_Data struct that contains
 *      std::vector lists of our vertices and faces. See the struct.h file
 *      to see the code for the Mesh_Data struct. It should be fairly straight-
 *      forward to figure out how to use it based on the code, since the structs
 *      are all very simple. One thing to note is that the first element pushed
 *      onto the vertices vector should be a NULL; this is just filler since
 *      .obj files 1-index vertices.
 *
 *     std::vector<HEV*> *hevs is just a new, empty vector of HEV structs.
 *      This vector is meant to be initialized and populated by our build_HE
 *      function. The idea is for hevs to be a list of our "halfedge vertices",
 *      which are vertex structs that also each contain a pointer to a halfedge.
 *      Note that this list is also 1-indexed, so the first element is NULL.
 *
 *     std::vector<HEF*> *hefs is similar to hevs in that it should just be
 *      a new, empty vector of HEF structs. The vector will be initialized
 *      and populated by our build_HE function with "halfedge faces", which are
 *      face structs that each contain a pointer to a halfedge. These face structs
 *      are only meant for accessing halfedges; they don't actually contain
 *      the indices of the vertices that make up the face as we can see in
 *      the code below. They don't need to anyway, since we can find all
 *      the vertices in the face by traversing through the halfedges.
 *
 * The following is a code snippet of how we might build a halfedge using
 * the build_HE function:
 *
 *     Mesh_Data *mesh_data = new Mesh_Data;
 *     parse_OBJ(mesh_data, filename); // your parsing function
 *
 *     // at this point, mesh_data contains our vertices and faces for the mesh
 *
 *     std::vector<HEV*> *hevs = new std::vector<HEV*>();
 *     std::vector<HEF*> *hefs = new std::vector<HEF*>();
 *
 *     build_HE(mesh_data, hevs, hefs);
 *     
 *     calculate normals for hevs
 *     use hevs indices to build normal list in mesh data in same order as vertices
 *     build normal/vertex buffers, and so on...
 *     delete_HE
 *
 * If the build is successful, then hevs and hefs should contain lists of
 * vertex and face structs that also have pointer data to halfedges.
 *
 * Now, to use our hevs and hefs lists for applications like computing
 * the area-weighted vertex normal, we can write code like the following:
 *
 * Vec3f calc_vertex_normal(HEV *vertex)
   {
       Vec3f normal;
       normal.x = 0;
       normal.y = 0;
       normal.z = 0;

       HE* he = vertex->out; // get outgoing halfedge from given vertex

       do
       {
           // compute the normal of the plane of the face
           Vec3f face_normal = calc_normal(he->face);
           // compute the area of the triangular face
           double face_area = calc_area(he->face);

           // accummulate onto our normal vector
           normal.x += face_normal.x * face_area;
           normal.y += face_normal.y * face_area;
           normal.z += face_normal.z * face_area;

           // gives us the halfedge to the next adjacent vertex
           he = he->flip->next;
    }
    while(he != vertex->out);
 *
 * The code above follows the loop structure that we explain in the lecture notes.
 * Refer there for details of why this loop traverses the neighborhood region
 * of the given vertex correctly. We can imagine writing a similar loop to help
 * us compute parts of the discrete Laplacian, since it is defined to be a sum
 * over vertices adjacent to the current vertex of interest.
 *
 * When we're done using our halfedge, we should call the delete_HE() function,
 * which takes in the vector of hevs and hefs built by build_HE. This delete function
 * frees all the memory that we set aside for our halfedge.
 *
 * Realize that the hevs and hefs vectors are meant to exist IN ADDITION to your regular
 * list of vertices and faces (i.e. the ones you passed into the build_HE function). The
 * idea is to mainly access the hevs and hefs vectors for halfedge-related computations
 * and still use the regular list of vertices and faces for the usual tasks like drawing
 * the mesh in OpenGL.
 *
 * If you have any questions or difficulties working with this code, don't
 * be afraid to send me an email at kevli@caltech.edu.
 */

#ifndef HALFEDGE_H
#define HALFEDGE_H

#include <algorithm>
#include <cassert>
#include <map>
#include <iostream>
#include <utility>
#include <vector>

#include "structs.h"

/* After this point, the comments stop. You shouldn't really need to know the
 * details of the following functions to know how to use this halfedge implementation.
 */

/* Halfedge structs */

struct HE // HE for halfedge
{
    // the vertex that this halfedge comes out off
    struct HEV *vertex;
    // the face adjacent to this halfedge
    struct HEF *face;
    // the flip and next halfedge as described in the lecture notes
    struct HE *flip, *next;

    // we omit the pointer to the adjacent edge (as well as a "halfedge edge"
    // struct) because it is not necessary for the assignment
};

struct HEF // HEF for halfedge face
{
    // the halfedge associated with this face
    struct HE *edge;
    // this variable is used to help orientate the halfedge when building it;
    // you don't have to worry about this
    bool oriented;
};

struct HEV // HEV for halfedge vertex
{
    // the coordinates of the vertex in the mesh
    double x, y, z;
    // the halfedge going out off this vertex
    struct HE *out;
    // use this to store your index for this vertex when you index the vertices
    // before building the operator for implicit fairing
    int index;
    // you can use this to store the normal vector for the vertex
    Vec3f normal;
};

/* Function prototypes */

static std::pair<int, int> get_edge_key(int x, int y);
static void hash_edge(std::map<std::pair<int, int>, HE*> &edge_hash,
                      std::pair<int, int> edge_key,
                      HE *edge);

static bool check_flip(HE *edge);
static bool check_edge(HE *edge);
static bool check_face(HEF *face);

static bool orient_flip_face(HE *edge);
static bool orient_face(HEF *face);

static bool build_HE(Mesh_Data *mesh,
                     std::vector<Face*> *f,
                     std::vector<HEV*> *hevs,
                     std::vector<HEF*> *hefs);

static void delete_HE(std::vector<HEV*> *hevs, std::vector<HEF*> *hefs);



/* Function implementations */

static std::pair<int, int> get_edge_key(int x, int y)
{
    assert(x != y);
    return std::pair<int, int>(std::min(x, y), std::max(x, y));
}

static void hash_edge(std::map<std::pair<int, int>, HE*> &edge_hash,
                     std::pair<int, int> edge_key,
                     HE *edge)
{
    if(edge_hash.count(edge_key) != 0)
    {
        HE *flip = edge_hash[edge_key];
        flip->flip = edge;
        edge->flip = flip;
    }
    else
        edge_hash[edge_key] = edge;
}

static bool check_flip(HE *edge)
{
    return edge->flip == NULL || edge->flip->vertex != edge->vertex;
}

static bool check_edge(HE *edge)
{
    return check_flip(edge) || !edge->face->oriented || !edge->flip->face->oriented;
}

static bool check_face(HEF *face)
{
    bool b1 = check_edge(face->edge);
    bool b2 = (face->edge->next != NULL) ? check_edge(face->edge->next) : 1;
    bool b3;

    if(face->edge->next != NULL)
        b3 = (face->edge->next->next != NULL) ? check_edge(face->edge->next->next) : 1;
    else
        b3 = 1;
    
    return b1 && b2 && b3;
}

static bool orient_flip_face(HE *edge)
{
    if(edge->flip == NULL)
        return 1;

    HE *flip = edge->flip;
    HEF *face = flip->face;

    if(face->oriented)
        return check_face(face);

    if (!check_flip(edge))
    {
        HEV *v1 = face->edge->vertex;
        HEV *v2 = face->edge->next->vertex;
        HEV *v3 = face->edge->next->next->vertex;
        
        assert(v1 != v2 && v1 != v3 && v2 != v3);

        HE *e3 = face->edge;
        HE *e1 = face->edge->next;
        HE *e2 = face->edge->next->next;

        assert(e3->vertex == v1);
        assert(e1->vertex == v2);
        assert(e2->vertex == v3);

        e3->vertex = v3;
        e3->next = e2;

        e1->vertex = v1;
        e1->next = e3;

        e2->vertex = v2;
        e2->next = e1;

        v1->out = e3;
        v2->out = e1;
        v3->out = e2;

        assert(face->edge->next->next->next == face->edge);
    }
    
    face->oriented = 1;

    assert(check_flip(edge));
    assert(check_face(face));

    return check_face(face) && orient_face(face);
}

static bool orient_face(HEF *face)
{
    assert(face->oriented);
    return orient_flip_face(face->edge)
           && orient_flip_face(face->edge->next)
           && orient_flip_face(face->edge->next->next)
           && check_face(face);
}

static bool build_HE(Mesh_Data *mesh,
                     std::vector<Face*> *f,
                     std::vector<HEV*> *hevs,
                     std::vector<HEF*> *hefs)
{
    std::vector<Vertex*> *vertices = mesh->vertices;
    std::vector<Face*> *faces = f;

    hevs->push_back(NULL);
    std::map<std::pair<int, int>, HE*> edge_hash;

    int size_vertices = vertices->size();

    for(int i = 1; i < size_vertices; ++i)
    {
        HEV *hev = new HEV;
        hev->x = vertices->at(i)->x;
        hev->y = vertices->at(i)->y;
        hev->z = vertices->at(i)->z;
        hev->out = NULL;
        hev->index = i;

        hevs->push_back(hev);
    }

    HEF *first_face = NULL;
    int num_faces = faces->size();
    
    for (int i = 0; i < num_faces; ++i)
    {
        Face *f = faces->at(i);

        HE *e1 = new HE;
        HE *e2 = new HE;
        HE *e3 = new HE;

        e1->flip = NULL;
        e2->flip = NULL;
        e3->flip = NULL;

        HEF *hef = new HEF;

        hef->oriented = 0;
        hef->edge = e1;

        e1->face = hef;
        e2->face = hef;
        e3->face = hef;
        
        e1->next = e2;
        e2->next = e3;
        e3->next = e1;

        e1->vertex = hevs->at(f->idx1);
        e2->vertex = hevs->at(f->idx2);
        e3->vertex = hevs->at(f->idx3);

        hevs->at(f->idx1)->out = e1;
        hevs->at(f->idx2)->out = e2;
        hevs->at(f->idx3)->out = e3;

        hash_edge(edge_hash, get_edge_key(f->idx1, f->idx2), e1);
        hash_edge(edge_hash, get_edge_key(f->idx2, f->idx3), e2);
        hash_edge(edge_hash, get_edge_key(f->idx3, f->idx1), e3);

        hefs->push_back(hef);

        if(first_face == NULL)
        {
            first_face = hef;
            first_face->oriented = 1;
        }
    }

    return orient_face(first_face);
}

static void delete_HE(std::vector<HEV*> *hevs, std::vector<HEF*> *hefs)
{
    int hev_size = hevs->size();
    int num_hefs = hefs->size();

    for(int i = 1; i < hev_size; ++i)
        delete hevs->at(i);

    for(int i = 0; i < num_hefs; ++i)
    {
        delete hefs->at(i)->edge->next->next;
        delete hefs->at(i)->edge->next;
        delete hefs->at(i)->edge;
        delete hefs->at(i);
    }

    delete hevs;
    delete hefs;
}

#endif
