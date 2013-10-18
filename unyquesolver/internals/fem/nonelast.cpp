#include "nonelast.hpp"

//------------------------------------------------------------------------------
NonElast::NonElast() {
}
//------------------------------------------------------------------------------
NonElast::NonElast(fem::FEM_PhysicalDomain *is, FEM_Common *ic) {
  s = is;
  c = ic;
  nelem = is->nelem;
  nnode = is->nnode; nbnode = is->nbnode;
  nedge = is->nedge; nbedge = is->nbedge;
  enode = 6;
}
//------------------------------------------------------------------------------
void NonElast::Init() {
  char filename[100];

  nelem = s->nelem;
  nnode = s->nnode; nbnode = s->nbnode;
  nedge = s->nedge; nbedge = s->nbedge;

  ENC = unyque::IMatrix_zero(nelem,7);

  // Initialize and set the integration points and weights
  Init_GIntegration();

  t = 0; dt = 0;
  D = unyque::DMatrix_zero(3,3);
  BL = unyque::DMatrix_zero(3,2*enode); BNL = unyque::DMatrix_zero(4,2*enode);
  BR0 = unyque::DMatrix_zero(2,2*enode); BR1 = unyque::DMatrix_zero(2,2*enode);
  Finv = unyque::DMatrix_zero(2,2); S = unyque::DMatrix_zero(2,2);
  Smat = unyque::DMatrix_zero(4,4);
  dp = unyque::DMatrix_zero(2, enode);
  hatSvec = unyque::DVector_zero(3);

  sprintf(filename,"./conf/mechanical.conf");
  ReadNonelast(filename);
  SetD();
  GenerateENC();
  MapDOFs();
}
//------------------------------------------------------------------------------
void NonElast::Init_GIntegration() {

  Genquad = 2; // Number of integration points for 1D edge integration
  Gnquad = 4;  // Number of integration points for 2D element integration

  Ges = unyque::DVector_zero(3*Genquad);
  Get = unyque::DVector_zero(3*Genquad);
  Gew = unyque::DVector_scalar(3*Genquad,0.5);

  // Setting the integration point for 1D edge integration
  Ges(0) = 0.211324865405187;
  Ges(1) = 0.788675134594813;
  Ges(2) = 0.788675134594813;
  Ges(3) = 0.211324865405187;
  Ges(4) = 0.0;
  Ges(5) = 0.0;

  Get(0) = 0.0;
  Get(1) = 0.0;
  Get(2) = 0.211324865405187;
  Get(3) = 0.788675134594813;
  Get(4) = 0.788675134594813;
  Get(5) = 0.211324865405187;

  Gs = unyque::DVector_zero(Gnquad);
  Gt = unyque::DVector_zero(Gnquad);
  Gw = unyque::DVector_zero(Gnquad);

  switch (Gnquad) {
  case 1:
    Gs(0) = 0.33333333333333; Gt(0) = 0.33333333333333;
    Gw(0) = 0.50000000000000;
    break;

  case 3:
    Gs(0) = 0.50000000000000; Gt(0) = 0.50000000000000;
    Gs(1) = 0.50000000000000; Gt(1) = 0.0;
    Gs(2) = 0.0;              Gt(2) = 0.50000000000000;

    Gw(0) = 0.16666666666667;
    Gw(1) = 0.16666666666667;
    Gw(2) = 0.16666666666667;
    break;

  case 4:
    Gs(0) = 0.33333; Gt(0) = 0.33333;
    Gs(1) = 0.6;     Gt(1) = 0.2;
    Gs(2) = 0.2;     Gt(2) = 0.6;
    Gs(3) = 0.2;     Gt(3) = 0.2;

    Gw(0) = -0.28125;
    Gw(1) = 0.26042;
    Gw(2) = 0.26042;
    Gw(3) = 0.26042;
    break;

  case 7:
    Gs(0) = 0.33333;  Gt(0) = 0.33333;
    Gs(1) = 0.79743;  Gt(1) = 0.10129;
    Gs(2) = 0.10129;  Gt(2) = 0.79743;
    Gs(3) = 0.10129;  Gt(3) = 0.10129;
    Gs(4) = 0.47014;  Gt(4) = 0.059716;
    Gs(5) = 0.059716; Gt(5) = 0.47014;
    Gs(6) = 0.47014;  Gt(6) = 0.47014;

    Gw(0) = 0.1125;
    Gw(1) = 0.06297;
    Gw(2) = 0.06297;
    Gw(3) = 0.06297;
    Gw(4) = 0.066197;
    Gw(5) = 0.066197;
    Gw(6) = 0.066197;
    break;
  }
}
//------------------------------------------------------------------------------
void NonElast::ReadNonelast(char *filename) {
  string key, line, commentChars = "%#;";
  char first;
  ifstream fp(filename);

  if (fp.is_open()) {
    fp >> noskipws;
    while ( !fp.eof() ) {
      fp >> ws;
      first = fp.peek();
      if (commentChars.find(first) != string::npos) {
	getline(fp, line);
	continue;
      }
      fp >> key;

      if (key.compare("CONSTITUTIVE") == 0)
	fp >> ws >> iconstitutive;
      if (key.compare("E") == 0)
	fp >> ws >> EM;
      if (key.compare("POISSON") == 0)
	fp >> ws >> MU;
      if (key.compare("DENSITY") == 0)
	fp >> ws >> RHO;
      if (key.compare("THERMOELASTICITY") == 0)
	fp >> ws >> ithermoelasticity;
      if (key.compare("ELECTROSTATIC_FORCE") == 0)
	fp >> ws >> ielecforce;
      if (key.compare("T_REF") == 0)
	fp >> ws >> T_REF;
      if (key.compare("NBC") == 0) {
	fp >> ws >> nbc;
	BCvals = unyque::DMatrix_zero(nbc,2);
	BCtype = unyque::IMatrix_zero(nbc,2);
	for (int i = 0; i < nbc; i++) {
	  fp >> ws >> BCtype(i,0) >> ws >> BCtype(i,1)	\
	     >> ws >> BCvals(i,0) >> ws >> BCvals(i,1);
	  if (BCtype(i,1) == 0)
	    BCvals(i,0) = BCvals(i,1) = 0.0;
	}
      }

      fp >> ws;
    }
    fp.close();
  } else
    cout << "Unable to open file" << endl;

}
//------------------------------------------------------------------------------
void NonElast::SetD() {
  double Y = 1.0, v = MU, coef;
  if (iconstitutive == 1) { // Plane stress
    coef = Y/(1.0-v*v);
    D(0,0) = coef;
    D(0,1) = v*coef;
    D(1,0) = D(0,1);
    D(1,1) = D(0,0);
    D(2,2) = 0.5*(1.0-v)*coef;
    D(0,2) = 0.0; D(1,2) = 0.0; D(2,0) = 0.0; D(2,1) = 0.0;
  }
  else if (iconstitutive == 0) {
    coef = Y/(1.0+v)/(1.0-2.0*v);
    D(0,0) = (1.0-v)*coef;
    D(0,1) = v*coef;
    D(1,0) = D(0,1);
    D(1,1) = D(0,0);
    D(2,2) = 0.5*(1.0-2.0*v)*coef;
    D(0,2) = 0.0; D(1,2) = 0.0; D(2,0) = 0.0; D(2,1) = 0.0;
  }
  else{
    if (c->DEBUG)
      cout<<"Invalid constitutive relation specified!\n"<<endl;
    exit(1);
  }
}
//------------------------------------------------------------------------------
void NonElast::GenerateENC() {
  fem::FEM_Element *ee;
  for (int i = 1; i <= nelem; i++) {
    ee = s->Elements[i];
    ENC(i-1,0) = i;
    ENC(i-1,1) = ee->node1;
    ENC(i-1,2) = ee->node2;
    ENC(i-1,3) = ee->node3;
    ENC(i-1,4) = ee->node4;
    ENC(i-1,5) = ee->node5;
    ENC(i-1,6) = ee->node6;
  }
}
//------------------------------------------------------------------------------
void NonElast::MapDOFs() {

  fem::FEM_Edge *ed;
  int bcno;

  L2G = unyque::IVector_zero(nnode);

  // First add all nodes to the global set
  for (int i = 0; i < nelem; i++) {
    for (int j = 1; j <= enode; j++) // Loop over nodes in that element
      L2G(ENC(i,j)-1) = 1;
  }

  // Loop over boundary edges
  for (int eid = 0; eid < nbedge; eid++) {

    ed = s->BEdges[eid+1];
    bcno = -1;
    // Find the b.c. no: corresponding to the edge's marker (Default -1)
    for (int i = 0; i < nbc; i++) {
      if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
    }

    // If it is a Dirichlet b.c. then remove nodes from global set
    if (BCtype(bcno,1) == 1) {
      L2G(ed->node1 - 1) = 0;
      L2G(ed->node2 - 1) = 0;
      L2G(ed->node3 - 1) = 0;
    }

  }  // End of loop over boundary edges

  ndof = 0;

  // Count number of dofs that are included
  for (int i = 0; i < nnode; i++) {
    ndof += L2G(i);
    L2G(i) *= 2*ndof-1; // Store position in global matrix
  }

}
//------------------------------------------------------------------------------
void NonElast::SolveStatic() {
  int iter = 0;
  double err, eps = 1e-10;

  ApplyInhomogeneousDBC();
  do {
    K.resize(2*ndof,2*ndof);
    RHS = unyque::DVector_zero(2*ndof);
    dU = unyque::DVector_zero(2*ndof);
    CompK();
    ApplyBC();
    dU = unyque::umfpackSolve(K, RHS);
    err = ublas::norm_inf(dU);
    ConstructGlobalU();
    iter++;
    if (c->DEBUG) cout<<"iteration: "<<iter<<"  Error: "<<err<<endl;
  } while (err > eps);
  if (c->DEBUG) PrintResults();
}
//------------------------------------------------------------------------------
void NonElast::PreProcess() {

  // Store values from previous time step
  for (int i = 0; i < nnode; i++) {
    (s->Uold)(i) = (s->U)(i);
    (s->Vold)(i) = (s->V)(i);
    (s->Udold)(i) = (s->Ud)(i);
    (s->Vdold)(i) = (s->Vd)(i);
    (s->Uddold)(i) = (s->Udd)(i);
    (s->Vddold)(i) = (s->Vdd)(i);
  }

}
//------------------------------------------------------------------------------
void NonElast::SolveDynamic(double tn, double dtn) {
  int iter = 0;
  double err, eps = 1e-10;
  t = tn, dt = dtn;

  ApplyInhomogeneousDBC();
  do {
    K.resize(2*ndof,2*ndof);
    RHS = unyque::DVector_zero(2*ndof);
    dU = unyque::DVector_zero(2*ndof);
    CompK();
    ApplyBC();
    dU = unyque::umfpackSolve(K, RHS);
    err = ublas::norm_inf(dU);
    ConstructGlobalU();
    iter++;
    if (c->DEBUG) cout<<"iteration: "<<iter<<"  Error: "<<err<<endl;
  } while (err > eps);
  ConstructGlobalUDyn();
  if (c->DEBUG) PrintResults();
}
//------------------------------------------------------------------------------
void NonElast::ApplyInhomogeneousDBC() {

  int bcno;
  fem::FEM_Edge *ed;

  // Loop over boundary edges
  for (int eid = 0; eid < nbedge; eid++) {

    ed = s->BEdges[eid+1];
    bcno = -1;
    // Find the b.c. no: corresponding to the edge's marker (Default -1)
    for (int i = 0; i < nbc; i++) {
      if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
    }

    // If it is a Dirichlet b.c. then initialize the corresponding nodes in the
    // solution vector with the value of the prescribed displacement
    if (BCtype(bcno,1) == 1) {
      (s->U)(ed->node1 - 1) = BCvals(bcno,0);
      (s->V)(ed->node1 - 1) = BCvals(bcno,1);
      (s->U)(ed->node2 - 1) = BCvals(bcno,0);
      (s->V)(ed->node2 - 1) = BCvals(bcno,1);
      (s->U)(ed->node3 - 1) = BCvals(bcno,0);
      (s->V)(ed->node3 - 1) = BCvals(bcno,1);
    }

  }  // End of loop over boundary edges

}
//------------------------------------------------------------------------------
void NonElast::CompK() {

  int Gnode, iglobal, jglobal;
  double si, ti;
  unyque::DMatrix Ke, Ktemp, Ue, Ecoor, dN;
  unyque::DVector Udde, RHSe, N;

  Ecoor = unyque::DMatrix_zero(enode,2);
  Ue = unyque::DMatrix_zero(enode, 2);
  Udde = unyque::DVector_zero(2*enode);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);

  // Loop over elements

  for (int eid = 0; eid < nelem; eid++) {
    for (int i = 0; i < enode; i++) {
      Gnode = ENC(eid,i+1); // Global node corresponding to (i+1)th local node
      Ue(i, 0) = (s->U)(Gnode-1);
      Ue(i, 1) = (s->V)(Gnode-1);
      Ecoor(i,0) = s->Nodes[Gnode]->x;
      Ecoor(i,1) = s->Nodes[Gnode]->y;
      if (dt > 0) {
	Udde(2*i) = -(s->Uddold)(Gnode-1) - 4./dt*(s->Udold)(Gnode-1) -
	  4./(dt*dt)*((s->Uold)(Gnode-1) - (s->U)(Gnode-1));
	Udde(2*i+1) = -(s->Vddold)(Gnode-1) - 4./dt*(s->Vdold)(Gnode-1) -
	  4./(dt*dt)*((s->Vold)(Gnode-1) - (s->V)(Gnode-1));
      }
    }

    // Loop over the Gauss points
    Ke = unyque::DMatrix_zero(2*enode,2*enode);
    RHSe = unyque::DVector_zero(2*enode);
    for (int ip = 0; ip < Gnquad; ip++) {
      si = Gs(ip); ti = Gt(ip);

      CompN(si, ti, N, dN);
      CompJacobian(dN, Ecoor);

      if (ithermoelasticity == 1)
	CompTEcoeff(eid, N);
      else
	TEcoeff = 0.0;

      CompBandS(Ue);

      // Compute the value of [B]'[D][B] at integration point ip
      // -- [BL]'[D][BL]: material part
      Ktemp = ublas::prod(D, BL);
      Ke += ublas::prod(ublas::trans(BL), Ktemp)*Gw(ip)*detJ;
      // -- [B_NL]'[Smat][B_NL]: geometric part
      Ktemp = ublas::prod(Smat, BNL);
      Ke += ublas::prod(ublas::trans(BNL), Ktemp)*Gw(ip)*detJ;

      // Compute [BL]'hatSvec part of the RHS
      RHSe += -ublas::prod(hatSvec, BL)*Gw(ip)*detJ;

      if (dt > 0) {
	CompF(eid);
	Ktemp = unyque::DMatrix_zero(2*enode, 2*enode);
	for (int i = 0; i < enode; i++) {
	  for (int j = 0; j < enode; j++) {
	    Ktemp(2*i, 2*j) = N(i)*N(j);
	    Ktemp(2*i + 1, 2*j + 1) = N(i)*N(j);
	  }
	}
	// Compute 4*rho*[N]'[N]*(1/dt^2)
	Ke += 1e-12*Ktemp*4*RHO/(dt*dt)*detF*Gw(ip)*detJ/EM;
	// Compute rho*[N]'[N]*Udd
	RHSe -= 1e-12*ublas::prod(Ktemp,Udde)*RHO*detF*Gw(ip)*detJ/EM;
      }

    } // End of loop over integration points

    // Assemble the elemental matrices Ke & RHSe to K and RHS respectively
    for (int i = 0; i < enode; i++) {
      iglobal = L2G(ENC(eid, i+1) - 1);
      if (iglobal > 0) {
	for (int j = 0; j < enode; j++) {
	  jglobal = L2G(ENC(eid,j+1) - 1);
	  if (jglobal > 0) {
	    K.append_element(iglobal-1, jglobal-1, Ke(2*i,2*j));
	    K.append_element(iglobal,   jglobal-1, Ke(2*i+1,2*j));
	    K.append_element(iglobal-1, jglobal,   Ke(2*i,2*j+1));
	    K.append_element(iglobal,   jglobal,   Ke(2*i+1,2*j+1));
	  }
	}
	RHS(iglobal-1) += RHSe(2*i);
	RHS(iglobal)   += RHSe(2*i+1);
      }
    }

  } // End of loop over elements

}
//------------------------------------------------------------------------------
void NonElast::CompN(double s, double t, unyque::DVector &N,
		     unyque::DMatrix &dN) {

  // Compute shape functions and their s- & t- derivatives, for triangle element
  if (enode == 6) {
    N(3) = 4.0*s*(1.0-s-t);
    N(4) = 4.0*s*t;
    N(5) = 4.0*t*(1.0-s-t);
    N(0) = (1.0-s-t)-N(3)/2.0-N(5)/2.0;
    N(1) = s-N(4)/2.0-N(3)/2.0;
    N(2) = t-N(5)/2.0-N(4)/2.0;

    dN(0, 3) = 4.0*(1.0-2.0*s-t);
    dN(0, 4) = 4.0*t;
    dN(0, 5) = -4.0*t;
    dN(0, 0) = -1.0-dN(0, 3)/2.0-dN(0, 5)/2.0;
    dN(0, 1) = 1.0-dN(0, 4)/2.0-dN(0, 3)/2.0;
    dN(0, 2) = -dN(0, 5)/2.0-dN(0, 4)/2.0;

    dN(1, 3) = -4.0*s;
    dN(1, 4) = 4.0*s;
    dN(1, 5) = 4.0*(1.0-s-2.0*t);
    dN(1, 0) = -1.0-dN(1, 3)/2.0-dN(1, 5)/2.0;
    dN(1, 1) = -dN(1, 4)/2.0-dN(1, 3)/2.0;
    dN(1, 2) = 1.0-dN(1, 5)/2.0-dN(1, 4)/2.0;
  }
  else if (enode == 3) {
    N(0) = (1.0-s-t);
    N(1) = s;
    N(2) = t;

    dN(0, 0) =-1.0;
    dN(0, 1) = 1.0;
    dN(0, 2) = 0.0;

    dN(1, 0) = -1.0;
    dN(1, 1) = 0.0;
    dN(1, 2) = 1.0;
  }
  else{
    if (c->DEBUG)
      cout<<"The number of nodes on each element is not 3 or 6!\n"<<endl;
    exit(1);
  }
}
//------------------------------------------------------------------------------
void NonElast::CompJacobian(unyque::DMatrix &dN, unyque::DMatrix &Ecoor) {

  unyque::DMatrix J, InvJ(2,2);

  // Compute Jacobian matrix
  J = ublas::prod(dN, Ecoor);

  // Compute DET,SLEN,TLEN,VLEN
  detJ = J(0,0)*J(1,1)-J(1,0)*J(0,1);
  slen = sqrt(J(0,0)*J(0,0)+J(0,1)*J(0,1));
  tlen = sqrt(J(1,0)*J(1,0)+J(1,1)*J(1,1));
  vlen = sqrt((J(1,0)-J(0,0))*(J(1,0)-J(0,0))+(J(1,1)-J(0,1))*(J(1,1)-J(0,1)));

  // Compute inverse of Jacobian matrix
  double idet;
  idet = 1.0/detJ;
  InvJ(0,0) = J(1,1)*idet;
  InvJ(0,1) = -J(0,1)*idet;
  InvJ(1,0) = -J(1,0)*idet;
  InvJ(1,1) = J(0,0)*idet;

  // Compute x- and y-derivatives of shape functions
  for (int i = 0; i < enode ; i++) {
    dp(0, i) = dN(0, i)*InvJ(0,0)+dN(1, i)*InvJ(0,1);
    dp(1, i) = dN(0, i)*InvJ(1,0)+dN(1, i)*InvJ(1,1);
  }

}
//------------------------------------------------------------------------------
void NonElast::CompTEcoeff(int eid, unyque::DVector &N) {
  int Gnode;
  double T;
  unyque::DVector Te;

  // Get the temperature at every node in the element
  Te = unyque::DVector_zero(enode);
  for (int i = 0; i < enode; i++) {
    Gnode = ENC(eid,i+1) - 1; // Global node corresponding to (i+1)th local node
    Te(i) = (s->T)(Gnode);
  }

  T = ublas::inner_prod(N, Te);

  ALPHA = c->functions->Eval_Alpha(T);

  TEcoeff = 0.0;
  if (iconstitutive == 1) { // Plane stress
    TEcoeff = ALPHA*(T - T_REF)*(1+MU*MU)/(1.0-MU*MU)/(1-2*MU);
  }
  else if (iconstitutive == 0) {
    TEcoeff = ALPHA*(T - T_REF)/(1-2*MU);
  }
}
//------------------------------------------------------------------------------
void NonElast::CompBandS(unyque::DMatrix &Ue) {

  double u11 = 0.0, u12 = 0.0, u21 = 0.0, u22 = 0.0;
  unyque::DMatrix BL0, BL1, F, E;

  BL0 = unyque::DMatrix_zero(3,2*enode); BL1 = unyque::DMatrix_zero(3,2*enode);
  F = unyque::DMatrix_zero(2,2); E = unyque::DMatrix_zero(2,2);

  // BL0
  for (int ib = 0; ib < enode; ib++) {
    BL0(0,2*ib) = dp(0, ib);
    BL0(0,2*ib+1) = 0.0;
    BL0(1,2*ib) = 0.0;
    BL0(1,2*ib+1) = dp(1, ib);
    BL0(2,2*ib) = dp(1, ib);
    BL0(2,2*ib+1) = dp(0, ib);
  }

  // BL1
  for (int ib = 0; ib < enode; ib++) {
    u11 += dp(0, ib)*Ue(ib, 0);
    u22 += dp(1, ib)*Ue(ib, 1);
    u12 += dp(1, ib)*Ue(ib, 0);
    u21 += dp(0, ib)*Ue(ib, 1);
  }
  for (int ib = 0; ib < enode; ib++) {
    BL1(0,2*ib) = u11*dp(0, ib);
    BL1(0,2*ib+1) = u21*dp(0, ib);
    BL1(1,2*ib) = u12*dp(1, ib);
    BL1(1,2*ib+1) = u22*dp(1, ib);
    BL1(2,2*ib) = u11*dp(1, ib)+u12*dp(0, ib);
    BL1(2,2*ib+1) = u21*dp(1, ib)+u22*dp(0, ib);
  }

  // BL = BL0 + BL1
  BL = BL0 + BL1;

  // BNL
  for (int ib = 0; ib < enode; ib++) {
    BNL(0,2*ib) = dp(0, ib);
    BNL(0,2*ib+1) = 0.0;
    BNL(1,2*ib) = dp(1, ib);
    BNL(1,2*ib+1) = 0.0;
    BNL(2,2*ib) = 0.0;
    BNL(2,2*ib+1) = dp(0, ib);
    BNL(3,2*ib) = 0.0;
    BNL(3,2*ib+1) = dp(1, ib);
  }

  // Compute F
  F(0,0) = 1+u11; F(0,1) = u12;
  F(1,0) = u21; F(1,1) = 1+u22;

  // Compute E
  E = 0.5*(ublas::prod(ublas::trans(F), F) - unyque::DMatrix_identity(2));

  // Compute S
  S(0,0) = D(0,0)*E(0,0)+D(0,1)*E(1,1)+2.0*D(0,2)*E(0,1)-TEcoeff;
  S(1,1) = D(1,0)*E(0,0)+D(1,1)*E(1,1)+2.0*D(1,2)*E(0,1)-TEcoeff;
  S(0,1) = D(2,0)*E(0,0)+D(2,1)*E(1,1)+2.0*D(2,2)*E(0,1);
  S(1,0) = S(0,1);

  // Compute Smat
  Smat(0,0) = S(0,0); Smat(0,1) = S(0,1); Smat(0,2) = 0.0; Smat(0,3) = 0.0;
  Smat(1,0) = S(1,0); Smat(1,1) = S(1,1); Smat(1,2) = 0.0; Smat(1,3) = 0.0;
  Smat(2,0) = 0.0; Smat(2,1) = 0.0; Smat(2,2) = S(0,0); Smat(2,3) = S(0,1);
  Smat(3,0) = 0.0; Smat(3,1) = 0.0; Smat(3,2) = S(1,0); Smat(3,3) = S(1,1);

  // Compute hatSvec
  hatSvec(0) = S(0,0);
  hatSvec(1) = S(1,1);
  hatSvec(2) = S(0,1);

}
//------------------------------------------------------------------------------
void NonElast::ApplyBC() {
  int bcno;
  fem::FEM_Edge *ed;

  // Loop over boundary edges - NBCs first
  for (int eid = 0; eid < nbedge; eid++) {

    ed = s->BEdges[eid+1];
    bcno = -1;
    // Find the b.c. no: corresponding to the edge's marker (Default -1)
    for (int i = 0; i < nbc; i++) {
      if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
    }

    // Check if b.c. is of type 2 and apply NBC if it is
    if (BCtype(bcno,1) == 2)
      ApplyNBC(bcno, ed);

    // If b.c. is of type 0 or 2, apply elec & fluid traction BCs
    if ((BCtype(bcno,1) == 0)||(BCtype(bcno,1) == 2)) {
      if (ielecforce == 1)
	ApplyElecBC(ed);
      ApplyFluidPressureBC(ed);
    }

  } //End of loop over boundary edges

  // // Loop over boundary edges - DBCs finally
  // for (int eid = 0; eid < nbedge; eid++) {

  //   ed = s->BEdges[eid+1];
  //   bcno = -1;
  //   // Find the b.c. no: corresponding to the edge's marker (Default -1)
  //   for (int i = 0; i < nbc; i++) {
  //     if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
  //   }

  //   // Check if b.c. is of type 1 and apply DBC if it is
  //   if (BCtype(bcno,1) == 1) {
  //     ApplyDBC((ed->node1 - 1));
  //     ApplyDBC((ed->node2 - 1));
  //     ApplyDBC((ed->node3 - 1));
  //   }

  // }  //End of loop over boundary edges
}
//------------------------------------------------------------------------------
void NonElast::ApplyNBC(int bcno, fem::FEM_Edge *ed) {

  int Gnode, startip, iglobal;
  double si, ti, tx, ty, JFinvTN, Hx, Hy, detJ1D = 0.0;
  unyque::DMatrix Ecoor, dN;
  unyque::DVector n(2), N;

  Ecoor = unyque::DMatrix_zero(enode,2);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);

  // Find the x & y coordinates of each node
  for (int i = 0; i < enode; i++) {
    Gnode = ENC(ed->eno-1, i+1); //global node corr. to (i+1)th local node
    Ecoor(i,0) = s->Nodes[Gnode]->x;
    Ecoor(i,1) = s->Nodes[Gnode]->y;
  }

  // Find the value of the traction vector
  tx = BCvals(bcno,0);
  ty = BCvals(bcno,1);

  // Compute components of normal vector
  n(0) = cos(ed->normal);
  n(1) = sin(ed->normal);

  startip = Genquad*(ed->eid);
  // Loop over the Gauss points for this edge
  for (int ip = startip; ip < (startip + Genquad); ip++) {

    // Compute shape functions and lengths of edges
    si = Ges(ip); ti = Get(ip);
    CompN(si, ti, N, dN);
    CompJacobian(dN, Ecoor);
    CompF(ed->eno-1);

    // Compute [H] = detF * |[Finv]'*[N]| * [t]
    JFinvTN = detF*abs(ublas::norm_2(ublas::prod(ublas::trans(Finv), n)));
    Hx = JFinvTN*tx/EM; Hy = JFinvTN*ty/EM;

    // Find the length of this edge
    switch (ed->eid) {
    case 0:
      detJ1D = slen;
      break;
    case 1:
      detJ1D = vlen;
      break;
    case 2:
      detJ1D = tlen;
      break;
    default:
      if (c->DEBUG) cout<<"Invalid edge ID"<<endl;
    }

    // Assemble the global vector RHS
    for (int i = 0; i < enode; i++) {
      iglobal = L2G(ENC(ed->eno-1, i+1) - 1);
      if (iglobal > 0) {
	RHS(iglobal-1) += N(i)*Hx*Gew(ip)*detJ1D;
	RHS(iglobal)   += N(i)*Hy*Gew(ip)*detJ1D;
      }
    }

  }

}
//------------------------------------------------------------------------------
void NonElast::CompF(int eid) {
  int Gnode;
  double idetF;
  unyque::DMatrix Ue, F;

  Ue = unyque::DMatrix_zero(enode, 2);
  F = unyque::DMatrix_zero(2,2);
  Finv = unyque::DMatrix_zero(2,2);

  for (int i = 0; i < enode; i++) {
    Gnode = ENC(eid, i+1); // Global node corr. to (i+1)th local node
    Ue(i, 0) = (s->U)(Gnode-1);
    Ue(i, 1) = (s->V)(Gnode-1);
  }

  F = unyque::DMatrix_identity(2) + ublas::trans(ublas::prod(dp, Ue));

  detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);
  idetF = 1.0/detF;
  Finv(0,0) = idetF*F(1,1); Finv(0,1) = -idetF*F(0,1);
  Finv(1,0) = -idetF*F(1,0); Finv(1,1) = idetF*F(0,0);
}
//------------------------------------------------------------------------------
void NonElast::ApplyElecBC(fem::FEM_Edge *ed) {

  int Gnode, startip, iglobal, jglobal;
  double scharge, bdphidn, si, ti, FinvTN;
  double En, Et, detJ1D = 0.0;
  unyque::DMatrix Ecoor, Ke, dN;
  unyque::DVector N, nX(2), nx(2), tx(2), f(2), H(2), RHSe;

  Ecoor = unyque::DMatrix_zero(enode,2);
  Ke = unyque::DMatrix_zero(2*enode,2*enode);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);
  RHSe = unyque::DVector_zero(2*enode);

  // Find the x & y coordinates of each node
  for (int i = 0; i < enode; i++) {
    Gnode = ENC(ed->eno-1, i+1); // Global node corr. to (i+1)th local node
    Ecoor(i,0) = s->Nodes[Gnode]->x;
    Ecoor(i,1) = s->Nodes[Gnode]->y;
  }
  bdphidn = c->Phi_mult*(s->BdPhidn)(ed->id-1);
  scharge = c->Phi_mult*(s->SCharge)(ed->id-1);
  En = c->Phi_mult*(s->Ent)(ed->id-1,0);
  Et = c->Phi_mult*(s->Ent)(ed->id-1,1);

  startip = Genquad*(ed->eid);
  // Loop over the Gauss points for this edge
  for (int ip = startip; ip < (startip + Genquad); ip++) {

    // Compute shape functions and lengths of edges
    si = Ges(ip); ti = Get(ip);
    CompN(si, ti, N, dN);
    CompJacobian(dN, Ecoor);
    CompF(ed->eno-1);

    // Find the components of normal & tangential vectors
    nX(0) = cos(ed->normal);
    nX(1) = sin(ed->normal);
    FinvTN = abs(ublas::norm_2(ublas::prod(ublas::trans(Finv), nX)));
    nx = ublas::prod(ublas::trans(Finv), nX)/FinvTN;
    tx(0) = -nx(1); tx(1) = nx(0);

    // Compute electrostatic force
    f = 1e12*0.5*scharge*((En - bdphidn)*nx + 2*Et*tx)/EM;

    // Compute [H] = detF * |[Finv]'*[N]| * [f]
    H = detF*FinvTN*f;

    // Find the length of this edge
    switch (ed->eid) {
    case 0:
      detJ1D = slen;
      break;
    case 1:
      detJ1D = vlen;
      break;
    case 2:
      detJ1D = tlen;
      break;
    default:
      if (c->DEBUG) cout<<"Invalid edge ID"<<endl;
    }

    // Create the elemental matrix Ke
    CompBR(N, nX);
    for (int i = 0; i < 2*enode; i++) {
      for (int j = 0; j < 2*enode; j++) {
 	Ke(i,j) = 1e12*0.5*scharge*(
		    BR0(0,i)*(En - bdphidn)*BR1(0,j) +
		    BR0(0,i)*(2*Et)*BR1(1,j) +
		    BR0(1,i)*(2*Et)*BR1(0,j) +
		    BR0(1,i)*(-En + bdphidn)*BR1(1,j))*Gew(ip)*detJ1D/EM;
      }
    }

    // Create the elemental vector RHSe
    for (int i = 0; i < enode; i++) {
      RHSe(2*i) += N(i)*H(0)*Gew(ip)*detJ1D;
      RHSe(2*i+1) += N(i)*H(1)*Gew(ip)*detJ1D;
    }

  }

  // Assemble the elemental matrices Ke and RHSe to K and RHS respectively
  for (int i = 0; i < enode; i++) {
    iglobal = L2G(ENC(ed->eno-1, i+1)-1);
    if (iglobal > 0) {
      for (int j = 0; j < enode; j++) {
	jglobal = L2G(ENC(ed->eno-1,j+1)-1);
	if (jglobal > 0) {
	  K.append_element(iglobal-1, jglobal-1, Ke(2*i,2*j));
	  K.append_element(iglobal,   jglobal-1, Ke(2*i+1,2*j));
	  K.append_element(iglobal-1, jglobal,   Ke(2*i,2*j+1));
	  K.append_element(iglobal,   jglobal,   Ke(2*i+1,2*j+1));
	}
      }
      RHS(iglobal-1) += RHSe(2*i);
      RHS(iglobal)   += RHSe(2*i+1);
    }
  }

}
//------------------------------------------------------------------------------
void NonElast::CompBR(unyque::DVector &N, unyque::DVector nX) {

  for (int ib = 0; ib < enode; ib++) {

    BR0(0,2*ib)   = N(ib);
    BR0(0,2*ib+1) = 0.0;
    BR0(1,2*ib)   = 0.0;
    BR0(1,2*ib+1) = N(ib);

    BR1(0,2*ib)   = 0.0;
    BR1(0,2*ib+1) = dp(1, ib)*nX(0) - dp(0, ib)*nX(1);
    BR1(1,2*ib)   = dp(1, ib)*nX(0) - dp(0, ib)*nX(1);
    BR1(1,2*ib+1) = 0.0;

  }

}
//------------------------------------------------------------------------------
void NonElast::ApplyFluidPressureBC(fem::FEM_Edge *ed) {

  int Gnode, startip, iglobal, jglobal;
  double si, ti, FinvTN, Pf, detJ1D = 0.0;
  unyque::DMatrix Ecoor, Ke, dN;
  unyque::DVector N, Pfe, nX(2), nx(2), tx(2), f(2), H(2), RHSe;

  Ecoor = unyque::DMatrix_zero(enode,2);
  Ke = unyque::DMatrix_zero(2*enode,2*enode);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);
  Pfe = unyque::DVector_zero(enode);
  RHSe = unyque::DVector_zero(2*enode);

  // Find the x & y coordinates of each node
  for (int i = 0; i < enode; i++) {
    Gnode = ENC(ed->eno-1, i+1); // Global node corr. to (i+1)th local node
    Ecoor(i,0) = s->Nodes[Gnode]->x;
    Ecoor(i,1) = s->Nodes[Gnode]->y;
    Pfe(i) = (s->Pf)(Gnode - 1);
  }

  startip = Genquad*(ed->eid);
  // Loop over the Gauss points for this edge
  for (int ip = startip; ip < (startip + Genquad); ip++) {

    // Compute shape functions and lengths of edges
    si = Ges(ip); ti = Get(ip);
    CompN(si, ti, N, dN);
    CompJacobian(dN, Ecoor);
    CompF(ed->eno-1);

    // Compute the integrated fluid pressure at this point
    Pf = ublas::inner_prod(N, Pfe);

    // Find the components of normal & tangential vectors
    nX(0) = cos(ed->normal);
    nX(1) = sin(ed->normal);
    FinvTN = abs(ublas::norm_2(ublas::prod(ublas::trans(Finv), nX)));
    nx = ublas::prod(ublas::trans(Finv), nX)/FinvTN;
    tx(0) = -nx(1); tx(1) = nx(0);

    // Compute force per unit length due to fluid pressure
    f = (-Pf*nx + 0.*tx)/EM;

    // Compute [H] = detF * |[Finv]'*[N]| * [f]
    H = detF*FinvTN*f;

    // Find the length of this edge
    switch (ed->eid) {
    case 0:
      detJ1D = slen;
      break;
    case 1:
      detJ1D = vlen;
      break;
    case 2:
      detJ1D = tlen;
      break;
    default:
      if (c->DEBUG) cout<<"Invalid edge ID"<<endl;
    }

    // Create the elemental matrix Ke
    CompBR(N, nX);
    for (int i = 0; i < 2*enode; i++) {
      for (int j = 0; j < 2*enode; j++) {
 	Ke(i,j) = (BR0(0,i)*(-Pf)*BR1(0,j) +
		   BR0(0,i)*(0.)*BR1(1,j) +
		   BR0(1,i)*(0.)*BR1(0,j) +
		   BR0(1,i)*(Pf)*BR1(1,j))*Gew(ip)*detJ1D/EM;
      }
    }

    // Create the elemental vector RHSe
    for (int i = 0; i < enode; i++) {
      RHSe(2*i) += N(i)*H(0)*Gew(ip)*detJ1D;
      RHSe(2*i+1) += N(i)*H(1)*Gew(ip)*detJ1D;
    }

  }

  // Assemble the elemental matrices Ke and RHSe to K and RHS respectively
  for (int i = 0; i < enode; i++) {
    iglobal = L2G(ENC(ed->eno-1, i+1)-1);
    if (iglobal > 0) {
      for (int j = 0; j < enode; j++) {
	jglobal = L2G(ENC(ed->eno-1,j+1)-1);
	if (jglobal > 0) {
	  K.append_element(iglobal-1, jglobal-1, Ke(2*i,2*j));
	  K.append_element(iglobal,   jglobal-1, Ke(2*i+1,2*j));
	  K.append_element(iglobal-1, jglobal,   Ke(2*i,2*j+1));
	  K.append_element(iglobal,   jglobal,   Ke(2*i+1,2*j+1));
	}
      }
      RHS(iglobal-1) += RHSe(2*i);
      RHS(iglobal)   += RHSe(2*i+1);
    }
  }

}
//------------------------------------------------------------------------------
void NonElast::ApplyDBC(int nid) {

  // Zero out the rows '2*nid' & '2*nid+1' and columns '2*nid' & '2*nid+1' in K
  for (int i = 0; i<2*nnode; i++) {
    if (K(2*nid, i) != 0.0)
      K.erase_element(2*nid, i);
    if (K(2*nid+1, i) != 0.0)
      K.erase_element(2*nid+1, i);
    if (K(i, 2*nid) != 0.0)
      K.erase_element(i, 2*nid);
    if (K(i, 2*nid+1) != 0.0)
      K.erase_element(i, 2*nid+1);
  }

  // Set K(2*nid,2*nid) & K(2*nid+1,2*nid+1) to one
  K(2*nid, 2*nid) = 1.0;
  K(2*nid+1, 2*nid+1) = 1.0;

  // Set RHS(2*nid) & RHS(2*nid+1) equal to the values of
  // the u & v Dirichlet boundary conditions
  RHS(2*nid) = 0.0;
  RHS(2*nid+1) = 0.0;

}
//------------------------------------------------------------------------------
void NonElast::ConstructGlobalU() {
  for (int i = 0; i < nnode; i++) {
    if (L2G(i) > 0) {
      (s->U)(i) += dU(L2G(i) - 1);
      (s->V)(i) += dU(L2G(i));
    }
  }
}
//------------------------------------------------------------------------------
void NonElast::ConstructGlobalUDyn() {
  for (int i = 0; i < nnode; i++) {

      // Compute accelerations at the new time-step
      (s->Udd)(i) = -(s->Uddold)(i) -
	4.0/dt*((s->Udold)(i) + ((s->Uold)(i)-(s->U)(i))/dt);
      (s->Vdd)(i) = -(s->Vddold)(i) -
	4.0/dt*((s->Vdold)(i) + ((s->Vold)(i)-(s->V)(i))/dt);

      // Compute velocities at the new time-step
      (s->Ud)(i) = (s->Udold)(i) + dt/2*((s->Udd)(i) + (s->Uddold)(i));
      (s->Vd)(i) = (s->Vdold)(i) + dt/2*((s->Vdd)(i) + (s->Vddold)(i));

  }
}
//------------------------------------------------------------------------------
void NonElast::PrintResults() {
  FILE *fp;
  fem::FEM_Point *pp;
  fp = fopen("postproc/UV.dat","w");
  for (int i = 0; i < nnode; i++) {
    pp = s->Nodes[i+1];
    fprintf(fp,"%.14lf  %.14lf %.14lf %.14lf \n",pp->x, pp->y,
	    (s->U)(i), (s->V)(i));
  }
  fclose(fp);
}
//------------------------------------------------------------------------------
double NonElast::MaxAbsDisp(int direction) {
  double rvalue = 0.0, disp;
  if (direction > 0)
    rvalue = ublas::norm_inf(s->U);
  else if (direction < 0)
    rvalue = ublas::norm_inf(s->V);
  else {
    rvalue = 0;
    for (int i = 0; i < nnode; i++) {
      disp = sqrt(pow((s->U)(i),2) + pow((s->V)(i),2));
      if (disp > rvalue)
	rvalue = disp;
    }
  }
  return rvalue;
}
//------------------------------------------------------------------------------
int NonElast::MaxAbsDispPoint(int direction) {
  int maxPoint = 1;
  double disp, maxDisp;
  if (direction > 0)
    maxPoint = ublas::index_norm_inf(s->U) + 1;
  else if (direction < 0)
    maxPoint = ublas::index_norm_inf(s->V) + 1;
  else {
    maxDisp = sqrt(pow((s->U)(maxPoint - 1),2) + pow((s->V)(maxPoint - 1),2));
    for (int j = 2; j <= s->nbnode; j++) {
      disp = sqrt(pow((s->U)(j - 1),2) + pow((s->V)(j - 1),2));
      if (disp > maxDisp) {
	maxPoint = j;
	maxDisp = disp;
      }
    }
  }
  return maxPoint;
}
//------------------------------------------------------------------------------
bp::list NonElast::DispBoundaryEdge(int bmarker, int direction) {

  bp::dict displacements;
  fem::FEM_Edge *ed;
  double disp, x, y;
  bp::list dispvalues, rvalue;
  bp::tuple item;
  int index;

  // Loop over boundary edges
  for (int eid = 0; eid < nbedge; eid++) {

    ed = s->BEdges[eid+1];

    if (ed->bmarker == bmarker) {
      if (direction > 0) { // Displacement along X
	displacements.setdefault(ed->node1, (s->U)(ed->node1 - 1));
	displacements.setdefault(ed->node2, (s->U)(ed->node2 - 1));
	displacements.setdefault(ed->node3, (s->U)(ed->node3 - 1));
      } else if (direction < 0) { // Displacement along Y
	displacements.setdefault(ed->node1, (s->V)(ed->node1 - 1));
	displacements.setdefault(ed->node2, (s->V)(ed->node2 - 1));
	displacements.setdefault(ed->node3, (s->V)(ed->node3 - 1));
      } else { // Magnitude of displacement
	disp = sqrt(pow((s->U)(ed->node1-1),2) + pow((s->V)(ed->node1-1),2));
	displacements.setdefault(ed->node1, disp);
	disp = sqrt(pow((s->U)(ed->node2-1),2) + pow((s->V)(ed->node2-1),2));
	displacements.setdefault(ed->node2, disp);
	disp = sqrt(pow((s->U)(ed->node3-1),2) + pow((s->V)(ed->node3-1),2));
	displacements.setdefault(ed->node3, disp);
      }
    }

  } // End of loop over boundary edges

  dispvalues = (bp::list)displacements.iteritems();
  for (bp::ssize_t i = 0; i < bp::len(dispvalues); i++) {
    item = bp::extract<bp::tuple>(dispvalues[i]);
    index = bp::extract<int>(item[0]);
    x = s->Nodes[index]->x + (s->U)(index - 1);
    y = s->Nodes[index]->y + (s->V)(index - 1);
    disp = bp::extract<double>(item[1]);
    rvalue.append(bp::make_tuple(x, y, disp));
  }

  return rvalue;

}
//------------------------------------------------------------------------------
pyublas::numpy_vector<double> NonElast::Displacement(int direction) {
  pyublas::numpy_vector<double> rval((s->U).size());
  if (direction > 0)
    rval.assign(s->U);
  else if (direction < 0)
    rval.assign(s->V);
  else {
    for (int j = 0; j < s->nbnode; j++)
      rval(j) = sqrt(pow((s->U)(j),2) + pow((s->V)(j),2));
  }
  return rval;
}
