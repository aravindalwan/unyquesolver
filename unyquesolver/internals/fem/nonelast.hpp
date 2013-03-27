#ifndef __NONELAST_H
#define __NONELAST_H

#include "fem.hpp"
#include "ublas.hpp"

#define PLANE_STRESS 0
#define PLANE_STRAIN 1

class NonElast {
public:
  fem::FEM_Domain *s;
  FEM_Common *c;
  unyque::SparseMatrix K;
  unyque::DVector RHS, dU;
  unyque::IMatrix ENC;

  // Number of nodes, elements, boundary nodes, edges, boundary edges,
  // nodes in an element and boundary conditions
  int nnode, nelem, nbnode, nedge, nbedge, enode, nbc;

  // Mapping nodal variables to global degrees of freedom and vice versa
  unyque::IVector L2G;
  int ndof;

  // Temp datastructures
  unyque::IMatrix BCtype; //Boundary id and type of boundary condition
  unyque::DMatrix BCvals, BL, BNL, BR0, BR1; //BCvals - b.c. values
  unyque::DMatrix dp, Finv, S, Smat;
  unyque::DVector hatSvec;
  double detF, detJ, slen, tlen, vlen, TEcoeff;

  // Material properties
  double EM, MU, ALPHA, T_REF;
  int iconstitutive, ithermoelasticity, ielecforce;
  unyque::DMatrix D;
  void ReadNonelast(char *filename);
  void SetD();

  // Numerical integration
  int Gnquad, Genquad;
  unyque::DVector Gs, Gt, Gw;  // Quad pt (s,t) and weight w for 2D integration
  unyque::DVector Ges, Get, Gew; // Quad pt (s) and weight w for 1D integration
  void Init_GIntegration();

  // Constructors
  NonElast();
  NonElast(fem::FEM_Domain *is, FEM_Common *ic);

  // Global initialization
  void GenerateENC();
  void Init();
  void MapDOFs();

  // Preperation of global matrices and solving
  void SolveStatic();
  void ApplyInhomogeneousDBC();
  void CompK();
  void CompN(double s, double t, unyque::DVector &N, unyque::DMatrix &dN);
  void CompJacobian(unyque::DMatrix &dN, unyque::DMatrix &Ecoor);
  void CompBandS(unyque::DMatrix &Ue);
  void CompTEcoeff(int eid, unyque::DVector &N);
  void ApplyBC();
  void ApplyNBC(int bcno, fem::FEM_Edge *ed);
  void CompF(int eid);
  void ApplyElecBC(fem::FEM_Edge *ed);
  void CompBR(unyque::DVector &N, unyque::DVector nX);
  void ApplyDBC(int nid);
  void ConstructGlobalU();

  // Post-processing
  void PrintResults();
  double MaxAbsDisp(int direction);
  int MaxAbsDispPoint(int direction);
  bp::list DispBoundaryEdge(int bmarker, int direction);
};
#endif
