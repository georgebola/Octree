#include <iostream>
// The Octree library resides in a single header file.
#include "Octree.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <stdlib.h>

#include "obj_reader/obj_reader.cpp"
#include "Stopwatch.h"

#include <pangolin/pangolin.h>
using namespace std;

objReader pointcloud;

// Any 3D vector implementation can be used, as long as it can be converted
// to float[3]. Here we use MyPoint.
struct MyPoint
{
    float x; 
    float y; 
    float z;
    MyPoint(const MyPoint& p2): x(p2.x), y(p2.y), z(p2.z) {}
    MyPoint& operator=(const MyPoint& p2) { x = p2.x; y = p2.y; z = p2.z; return *this;}
    MyPoint(float in_x, float in_y, float in_z): x(in_x), y(in_y), z(in_z) {}
    MyPoint(const float p2[3]): x(p2[0]), y(p2[1]), z(p2[2]) {}
    operator float*() { return &x; }
    operator const float*() const { return &x; }
    MyPoint operator+(const MyPoint& p2) const { return MyPoint(x+p2.x, y+p2.y, z+p2.z); }
    MyPoint operator-(const MyPoint& p2) const { return MyPoint(x-p2.x, y-p2.y, z-p2.z); }
    MyPoint operator*(float f) const { return MyPoint(x*f, y*f, z*f); }
    float lengthsq() { return x*x + y*y + z*z; }
};

// Data that we're manipulating.
struct Particle
{
    MyPoint pos;
    Particle(float x, float y, float z): pos(x, y, z) {}
};

// A node in the quadtree.
struct Node
{
    vector<Particle*> particles;
};

#define NUM_PARTICLES 10000000
#define R 1.0f
#define EPSILON 0.0001f

// This class does the work of counting particles while traversing the octree.
class CallbackTraverse: public Octree<Node>::Callback
{
public:
    int count;
    Particle* p0;
    virtual bool operator()(const float min[3], const float max[3], Node& n)
    {
        MyPoint pmin(min), pmax(max);
        float cellSizeSq = (pmax - pmin).lengthsq();
        float maxDist = (sqrtf(cellSizeSq) * 0.5f) + R + EPSILON;
        
        MyPoint center = (pmin + pmax) * 0.5f;
        MyPoint vectCenter = center - p0->pos;
        float distCenterSq = vectCenter.lengthsq();
        if (distCenterSq > maxDist * maxDist)
            return false; // Too far; don't subdivide cell.
        
        // Iterate through particles in this cell.
        vector<Particle*>::const_iterator it = n.particles.begin();
        for (; it != n.particles.end(); it++)
        {
            Particle* p = *it;
            if (p == p0)
                continue;
            float dsq = (p->pos - p0->pos).lengthsq();
            // If particle is within the radius, increment counter.
            if (dsq <= R * R)
                count++;
        }
        // Subdivide cell if needed.
        return true;
    }
};


objReader readObj(){
    objReader testobj;

    char* filename = "Z6.obj";
    testobj.objLoadFile(filename);
    testobj.objLoadModel();
    cout<<"No. of vertices: " << testobj.nVertex << endl;

    return testobj;
}


int main()
{
    
    pointcloud = readObj();
    const int nPoints = pointcloud.nVertex;
    cout << "Number of particles: " << nPoints << endl;
    cout << "Initializing particles, uniform distribution in cube ([0,1], [0,1], [0,1])." << endl;


    vector<Particle> myParticles;


    for (int i = 0; i < nPoints; i++){
        float x = (float)pointcloud.vertexArray[i].x;
        float y = (float)pointcloud.vertexArray[i].y;
        float z = (float)pointcloud.vertexArray[i].z;

        myParticles.push_back(Particle(x , y ,z ));

    }


    cout << endl;
    cout << "Count number of particles within a radius of " << R << "." << endl;
    cout << "Output format: 'p: c' where c is the number of particles" << endl;
    cout << "around particle p." << endl << endl;
    
    cout << "Computing result with brute force:" << endl;

    double start = stopwatch();
    for (int i = 0; i < 120; i++)
    {
        cout << setw(3) << i << ": ";
        int count = 0;
        MyPoint pi(myParticles[i].pos);
        for (int j = 0; j < nPoints; j++)
        {
            if (i == j)
                continue;
            MyPoint pj(myParticles[j].pos);
            float dsq = (pj - pi).lengthsq();
            if (dsq <= R * R)
                count++;
        }
        cout << count << "  ";
        if (i % 8 == 7)
            cout << endl;
    }

    double T = stopwatch() - start;
    printf("Brute Force in %.5f sec.\n", T);
    
    cout << endl << "Building octree." << endl;
    // Initialize octree.
    // Minimum coordinates.
    float min[3] = {0.0f, 0.0f, 0.0f};
    // Maximum coordinates.
    float max[3] = {211130.0, 211130.0, 211130.0};
    // Minimum size to use when subdividing cells.
    float cellSize[3] = {0.1, 0.1, 0.1};
    Octree<Node> octree(min, max, cellSize);
    
    // Add all particles to the octree.
    for (int i = 0; i < nPoints; i++)
    {
        Node& n = octree.getCell(myParticles[i].pos);
        n.particles.push_back(&myParticles[i]);
    }
    cout << "Building octree done." << endl << endl;
    
    cout << "Computing result with octree:" << endl;

    double start_octree = stopwatch();

    for (int i = 0; i < 120; i++)
    {
        cout << setw(3) << i << ": ";
        // Prepare the callback class.
        CallbackTraverse ct;
        ct.count = 0;
        ct.p0 = &myParticles[i];
        // Invoke the traversal.
        octree.traverse(&ct);
        cout << ct.count << "  ";
        if (i % 8 == 7)
            cout << endl;
    }
    double T_octree = stopwatch() - start_octree;
    printf("Octree in %.5f sec.\n", T_octree);
    return 0;
}
