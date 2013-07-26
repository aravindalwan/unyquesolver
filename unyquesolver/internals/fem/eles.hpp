#ifndef __ELES_HPP
#define __ELES_HPP

#include "fem.hpp"
#include "ublas.hpp"

class ElEs {

public:
  fem::FEM_PhysicalDomain *s;
  fem::FEM_Common *c;
  unyque::SparseMatrix K;
  unyque::DVector RHS, dU;
  unyque::IMatrix ENC;

  // Number of nodes, elements, boundary nodes, edges, boundary edges and b.c's
  int nnode, nelem, nbnode, nedge, nbedge, enode, nbc;

  // Mapping nodal variables to global degrees of freedom and vice versa
  unyque::IVector L2G_Phi, L2G_BPhi, L2G_BdPhidn;
  int nPhi, nBPhi, nBdPhidn;

  // Temporary datastructures
  unyque::IMatrix BCtype; //Boundary id and type
  // RegAttr is for reg. attributes, BCvals is for b.c values
  unyque::DMatrix BCvals, RegAttr;
  unyque::DVector Lengths;
  double detF, detJ, slen, tlen, vlen;

  // Material properties
  double EPS0, EC, EC0;
  void ReadElEs(char *filename);

  // Numerical integration
  int Gnquad, Genquad;
  unyque::DVector Gs, Gt, Gw;  // Quad pt (s,t) & weight w for 2D integration
  unyque::DVector Ges, Get, Gew; // Quad pt (s,t) & weight w for 1D integration
  void Init_GIntegration();

  // Constructors
  ElEs();
  ElEs(fem::FEM_PhysicalDomain *is, fem::FEM_Common *ic);

  // Global initialization
  void GenerateENC();
  void Init();
  void MapDOFs();

  // Preparation of global matrices and solving
  void SolveStatic();
  void ApplyInhomogeneousDBC();
  void CompKandRHS();
  void CompN(double s, double t, unyque::DVector &N, unyque::DMatrix &dN);
  unyque::DMatrix CompJandB(unyque::DMatrix &dN, unyque::DMatrix &Ecoor);
  unyque::DMatrix CompF(int eid, unyque::DMatrix &B);
  unyque::DVector CompIGGstar(int p, int q);
  unyque::DVector CompIGGstar2(int p, int q);
  void ApplyBC();
  void ApplyIBC(fem::FEM_Edge *ed);
  void ApplyDBC(int nid);
  void UpdateGlobalPotentials();
  void PrintResults();

private:
  void GaussLegendre(double x1, double x2, double* x, double* w, int n);

};
#endif
