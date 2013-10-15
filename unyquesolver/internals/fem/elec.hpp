#ifndef __ELEC_HPP
#define __ELEC_HPP

#include "fem.hpp"
#include "ublas.hpp"

class Elec {

public:
  boost::shared_ptr<fem::FEM_PhysicalDomain> s;
  boost::shared_ptr<fem::FEM_Common> c;
  unyque::SparseMatrix K;
  unyque::DVector RHS, dU;
  unyque::IMatrix ENC;

  // Number of nodes, elements, boundary nodes, edges, boundary edges and b.c's
  int nnode, nelem, nbnode, nedge, nbedge, enode, nbc;

  // Mapping nodal variables to global degrees of freedom and vice versa
  unyque::IVector L2G;
  int ndof;

  // Temporary datastructures
  unyque::IMatrix BCtype; // Boundary id and type
  unyque::DMatrix BCvals; // BCvals - boundary condition values
  double detF, detJ, slen, tlen, vlen;

  // Material properties
  double EC, Qj;
  void ReadElec(char *filename);

  // Numerical integration
  int Gnquad, Genquad;
  unyque::DVector Gs, Gt, Gw;  //Quad pt (s,t) & weight w for 2D integration
  unyque::DVector Ges, Get, Gew; //Quad pt (s,t) & weight w for 1D integration
  void Init_GIntegration();

  // Constructors
  Elec();
  Elec(boost::shared_ptr<fem::FEM_PhysicalDomain> is,
       boost::shared_ptr<fem::FEM_Common> ic);

  // Global initialization
  void GenerateENC();
  void Init();
  void MapDOFs();

  // Preparation of global matrices and solving
  void SolveStatic();
  void ApplyInhomogeneousDBC();
  void CompDomIntegrals();
  void CompN(double s, double t, unyque::DVector & N, unyque::DMatrix & dN);
  unyque::DMatrix CompJandB(unyque::DMatrix& dN, unyque::DMatrix& Ecoor);
  unyque::DMatrix CompF(int eid, unyque::DMatrix &B);
  void ApplyBC();
  void ApplyNBC(int bcno, fem::FEM_Edge *ed);
  void ApplyDBC(int nid);
  void UpdateGlobalPotential();
  void PrintResults();
};
#endif
