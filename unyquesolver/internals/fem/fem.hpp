#ifndef __FEM_H
#define __FEM_H

#include <vector>
#include <fstream>
#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <pyublas/numpy.hpp>
#include "ublas.hpp"
#include "function.hpp"
namespace bp = boost::python;

// Define point type for FEM_Point2D
#define INTERIOR 0
#define BOUNDARY 1

namespace fem {

  class FEM_Element {
  public:
    int id, node1, node2, node3, node4, node5, node6;
    int reg; // Region Marker for each element
    FEM_Element(int eid, int region,
		int n1, int n2, int n3, int n4, int n5, int n6);
  };
  //----------------------------------------------------------------------------
  class FEM_Edge {
  public:
    // Edge parameters - node3 only to store mid node for quadratic FEM basis
    int id, node1, node2, node3, eno, eid;
    int bmarker; // Boundary Marker for each edge
    double normal;
    FEM_Edge(int edid, int elementnum, int edgetype, int bm, double norm,
	     int n1, int n2, int n3);
  };
  //----------------------------------------------------------------------------
  class FEM_Point2D {
  public:
    double x, y; // deformed coordinates, not used in electrostatics
    int id;
    int bmarker; // Boundary Marker for each point
    FEM_Point2D(double ix, double iy);
    FEM_Point2D(int pid, int bm, double ix, double iy);
  };
  //----------------------------------------------------------------------------
  class FEM_Surface2D {
  public:
    vector<FEM_Point2D *> Nodes, RefNodes;
    vector<FEM_Point2D *> BNodes;
    vector<FEM_Element *> Elements;
    vector<FEM_Edge *> Edges;
    vector<FEM_Edge *> BEdges;

    //nnode: number of nodes, nbnode: number of boundary nodes
    //nelem: number of elements
    //nedge: number of edges, nbedge: number of edges on boundary
    int nnode,nbnode,nelem,nedge,nbedge;

    int id;

    // Mechanical displacements
    pyublas::numpy_vector<double> U, V;

    // Temperature
    pyublas::numpy_vector<double> T, Told;

    // Electric potential
    pyublas::numpy_vector<double> Phi;
    double potential; // Potential of equipotential surfaces

    // BEM electrostatics
    pyublas::numpy_vector<double> BPhi, BdPhidn, SCharge;
    unyque::DMatrix Ent;

  public:
    FEM_Surface2D();
    ~FEM_Surface2D() {};

    // Initialization
    void SetID(int a) {id = a;};
    void SetMesh(bp::list nodes, bp::list edges, bp::list elements);
    void SetNodes(bp::list nodes);
    void SetEdges(bp::list edges);
    void SetElements(bp::list elements);
    void InitDOFs();

    // Move mesh to include parametric variation of geometry
    void MoveMesh(unyque::DMatrix &displacement);

  };

}
//------------------------------------------------------------------------------
class FEM_Common {

public:
  // Switch controlling whether debug output is printed to stdout
  int DEBUG;

  // Multiplier to convert potential to non-dimensional form
  double Phi_mult;

  // Potential at infinity for electrostatics
  double Phi_inf;

  // Functions used to describe nonlinear properties
  Function *functions;

public:
  FEM_Common();
  ~FEM_Common() {};

};
#endif
