#ifndef __FEM_H
#define __FEM_H

#include <boost/python.hpp> // Python requirement: this has to be included first
#include <boost/python/list.hpp>
#include <vector>
#include <fstream>
#include <pyublas/numpy.hpp>
#include "ublas.hpp"
#include "function.hpp"
namespace bp = boost::python;

// Define point type for FEM_Point
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

  struct FEM_Element_pickle_suite : bp::pickle_suite {
    static bp::tuple getinitargs(const FEM_Element& e) {
      return bp::make_tuple(e.id, e.reg, e.node1, e.node2, e.node3, e.node4,
			    e.node5, e.node6);
    }
    static bp::tuple getstate(bp::object o);
    static void setstate(bp::object o, bp::tuple state);
    static bool getstate_manages_dict() { return true; }
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

  struct FEM_Edge_pickle_suite : bp::pickle_suite {
    static bp::tuple getinitargs(const FEM_Edge& e) {
      return bp::make_tuple(e.id, e.eno, e.eid, e.bmarker, e.normal,
			    e.node1, e.node2, e.node3);
    }
    static bp::tuple getstate(bp::object o);
    static void setstate(bp::object o, bp::tuple state);
    static bool getstate_manages_dict() { return true; }
  };
  //----------------------------------------------------------------------------
  class FEM_Point {
  public:
    double x, y; // deformed coordinates, not used in electrostatics
    int id;
    int bmarker; // Boundary Marker for each point
    FEM_Point(double ix, double iy);
    FEM_Point(int pid, int bm, double ix, double iy);
  };

  struct FEM_Point_pickle_suite : bp::pickle_suite {
    static bp::tuple getinitargs(const FEM_Point& p) {
      return bp::make_tuple(p.id, p.bmarker, p.x, p.y);
    }
    static bp::tuple getstate(bp::object o);
    static void setstate(bp::object o, bp::tuple state);
    static bool getstate_manages_dict() { return true; }
  };
  //----------------------------------------------------------------------------
  class FEM_Domain {
  public:
    vector<FEM_Point *> Nodes, RefNodes;
    vector<FEM_Point *> BNodes;
    vector<FEM_Element *> Elements;
    vector<FEM_Edge *> Edges;
    vector<FEM_Edge *> BEdges;

    //nnode: number of nodes, nbnode: number of boundary nodes
    //nelem: number of elements
    //nedge: number of edges, nbedge: number of edges on boundary
    int nnode,nbnode,nelem,nedge,nbedge;

    int id;

  public:
    FEM_Domain();
    virtual ~FEM_Domain() {};

    // Initialization
    void SetID(int a) {id = a;};
    void SetMesh(bp::list nodes, bp::list edges, bp::list elements);
    void SetNodes(bp::list nodes);
    void SetEdges(bp::list edges);
    void SetElements(bp::list elements);
    virtual void InitDOFs() = 0;

    // Move mesh to include parametric variation of geometry
    void MoveMesh(unyque::DMatrix &displacement);

  };
  //------------------------------------------------------------------------------
  class FEM_PhysicalDomain: public FEM_Domain {
  public:
    // Mechanical displacements
    pyublas::numpy_vector<double> U, V, Uold, Vold;

    // Mechanical velocity and acceleration
    pyublas::numpy_vector<double> Ud, Vd, Udd, Vdd, Udold, Vdold, Uddold, Vddold;

    // Traction due to fluid pressure
    pyublas::numpy_vector<double> Pf;

    // Temperature
    pyublas::numpy_vector<double> T, Told;

    // Electric potential
    pyublas::numpy_vector<double> Phi;
    double potential; // Potential of equipotential surfaces

    // BEM electrostatics
    pyublas::numpy_vector<double> BPhi, BdPhidn, SCharge;
    pyublas::numpy_matrix<double> Ent;

  public:
    FEM_PhysicalDomain();
    ~FEM_PhysicalDomain() {};

    // Initialization
    void InitDOFs();

  };
  //------------------------------------------------------------------------------
  class FEM_FluidDomain: public FEM_Domain {
  public:
    // Fluid pressure fluctuations
    pyublas::numpy_vector<double> P, Pold;

    // Projection of mechanical displacement on fluid domain
    pyublas::numpy_vector<double> U;

    // Gap height
    pyublas::numpy_vector<double> H, Hold;

  public:
    FEM_FluidDomain();
    ~FEM_FluidDomain() {};

    // Initialization
    void InitDOFs();

  };
  //------------------------------------------------------------------------------
  class FEM_Common {

  public:
    // Switch controlling whether debug output is printed to stdout
    int DEBUG;

    // Multiplier to convert potential to non-dimensional form
    double Phi_mult;

    // Potential at infinity for electrostatics
    double Phi_inf;

    // Original and new values of electrode gap
    double original_gap, new_gap;

    // Current time, stop time and time-step
    double t, t_stop, dt;

    // Functions used to describe nonlinear properties
    Function *functions;

  public:
    FEM_Common();
    ~FEM_Common() {};

  };

  struct FEM_Common_pickle_suite : bp::pickle_suite {
    static bp::tuple getstate(bp::object o);
    static void setstate(bp::object o, bp::tuple state);
    static bool getstate_manages_dict() { return true; }
  };

}
#endif
