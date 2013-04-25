#ifndef __FLUID_HPP
#define __FLUID_HPP

#include "fem.hpp"
#include "ublas.hpp"
#include "quad.hpp"

class Fluid {

public:
  fem::FEM_PhysicalDomain *s;
  fem::FEM_FluidDomain *sf;
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
  double ETA, LAMBDA, P_ATM, KN;
  int MOVING_EDGE, FIXED_EDGE; // IDs of moving and fixed edges in physical domain
  void ReadFluid(char *filename);

  // Numerical integration
  int Gnquad, Genquad, Gbnquad;
  unyque::DVector Gs, Gt, Gw;  // Quad pt (s,t) & weight w for 2D integration
  unyque::DVector Ges, Get, Gew; // Quad pt (s,t) & weight w for 1D integration
  unyque::DVector Gbs, Gbw; // 1D Quad pt and weight along breadth of domain
  double ZMIN, ZMAX;
  void Init_GIntegration();

  // Constructors
  Fluid();
  Fluid(fem::FEM_PhysicalDomain *is, fem::FEM_FluidDomain *isf, FEM_Common *ic);

  // Global initialization
  void GenerateENC();
  void Init();
  void MapDOFs();

  // Preparation of global matrices and solving
  void PreProcess();
  void SolveDynamic(double tn, double dtn);
  void MapPhysicalToFluid();
  void ApplyInhomogeneousDBC();
  void CompDomIntegrals();
  void CompN(double s, double t, unyque::DVector &N, unyque::DMatrix &dN);
  unyque::DMatrix CompJandB(unyque::DMatrix &dN, unyque::DMatrix &Ecoor);
  unyque::DMatrix CompF(int eid, unyque::DMatrix &B);
  void UpdateGlobalPressure();
  void PrintResults();
  void MapFluidToPhysical();
  pyublas::numpy_vector<double> Pressure();
};
#endif
