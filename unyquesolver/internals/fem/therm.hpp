#ifndef __THERM_HPP
#define __THERM_HPP

#include "fem.hpp"
#include "ublas.hpp"

class Therm {

public:
  fem::FEM_PhysicalDomain *s;
  FEM_Common *c;
  unyque::SparseMatrix K;
  unyque::DVector RHS, dU;
  unyque::IMatrix ENC;
  double t, dt;

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
  int JouleHeating;
  double TC, dTCdT, Qext, Q, Ht, Hb, Tinf, Thickness, RhoCp;
  void ReadTherm(char *filename);

  // Numerical integration
  int Gnquad, Genquad;
  unyque::DVector Gs, Gt, Gw;  // Quad pt (s,t) & weight w for 2D integration
  unyque::DVector Ges, Get, Gew; // Quad pt (s,t) & weight w for 1D integration
  void Init_GIntegration();

  // Constructors
  Therm();
  Therm(fem::FEM_PhysicalDomain *is, FEM_Common *ic);

  // Global initialization
  void GenerateENC();
  void Init();
  void MapDOFs();

  // Preparation of global matrices and solving
  void SolveStatic();
  void SolveDynamic(double tn, double dtn);
  void ApplyInhomogeneousDBC();
  void CompDomIntegrals();
  void CompN(double s, double t, unyque::DVector &N, unyque::DMatrix &dN);
  unyque::DMatrix CompJandB(unyque::DMatrix &dN, unyque::DMatrix &Ecoor);
  unyque::DMatrix CompF(int eid, unyque::DMatrix &B);
  double CompQjoule(int eid, unyque::DMatrix &B, double T);
  void ApplyBC();
  void ApplyNBC(int bcno, fem::FEM_Edge *ed);
  void ApplyRBC(int bcno, fem::FEM_Edge *ed);
  void ApplyDBC(int nid);
  void UpdateGlobalTemp();
  void PrintResults();
  double MaxAbsTemp();
};
#endif
