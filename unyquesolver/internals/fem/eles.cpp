#include "eles.hpp"

//------------------------------------------------------------------------------
ElEs::ElEs() {
}
//------------------------------------------------------------------------------
ElEs::ElEs(boost::shared_ptr<fem::FEM_PhysicalDomain> is,
	   boost::shared_ptr<fem::FEM_Common> ic) {
  s = is;
  c = ic;
  nelem = is->nelem;
  nnode = is->nnode; nbnode = is->nbnode;
  nedge = is->nedge; nbedge = is->nbedge;
  enode = 6;
}
//------------------------------------------------------------------------------
void ElEs::Init() {
  char filename[100];

  nelem = s->nelem;
  nnode = s->nnode; nbnode = s->nbnode;
  nedge = s->nedge; nbedge = s->nbedge;

  ENC = unyque::IMatrix_zero(nelem,7);

  // Initialize and set the integration points and weights
  Init_GIntegration();

  EPS0 = 8.85434e-12; EC = 0.0; EC0 = 0.0;
  Lengths = unyque::DVector_zero(nbedge);
  sprintf(filename,"./conf/electricalelectrostatic.conf");
  ReadElEs(filename);
  GenerateENC();
  MapDOFs();
}
//------------------------------------------------------------------------------
void ElEs::Init_GIntegration() {

  Genquad = 10; // Number of points for 1D edge integration - must be even
  Gnquad = 4;  // Number of points for 2D element integration

  Ges = unyque::DVector_zero(3*Genquad);
  Get = unyque::DVector_zero(3*Genquad);
  Gew = unyque::DVector_zero(3*Genquad);

  // Setting the integration point for 1D edge integration
  double gx[Genquad], gw[Genquad];
  GaussLegendre(0.0, 1.0, gx, gw, Genquad);
  for (int i = 0; i < Genquad; i++) {
    Ges(i) = gx[i]; Gew(i) = gw[i];
    Ges(Genquad+i) = 1-gx[i]; Get(Genquad+i) = gx[i]; Gew(Genquad+i) = gw[i];
    Get(2*Genquad+i) = 1-gx[i]; Gew(2*Genquad+i) = gw[i];
  }

  Gs = unyque::DVector_zero(Gnquad);
  Gt = unyque::DVector_zero(Gnquad);
  Gw = unyque::DVector_zero(Gnquad);

  // Setting the integration points for 2D element integration
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
    Gs(1) = 0.6; Gt(1) = 0.2;
    Gs(2) = 0.2; Gt(2) = 0.6;
    Gs(3) = 0.2; Gt(3) = 0.2;

    Gw(0) = -0.28125;
    Gw(1) =  0.26042;
    Gw(2) =  0.26042;
    Gw(3) =  0.26042;
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
void ElEs::ReadElEs(char *filename) {
  string key, line, commentChars = "%#;";
  char first;
  ifstream fp(filename);
  int num, id;

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

      if (key.compare("EC0") == 0)
	fp >> ws >> EC0;
      if (key.compare("PHI_MULT") == 0)
	fp >> ws >> c->Phi_mult;
      if (key.compare("NR") == 0) {
	fp >> ws >> num;
	RegAttr = unyque::DMatrix_zero(num,2);
	for (int i = 0; i < num; i++)
	  fp >> ws >> id >> ws >> RegAttr(i,0) >> ws >> RegAttr(i,1);
      }
      if (key.compare("NBC") == 0) {
	fp >> ws >> nbc;
	BCvals = unyque::DMatrix_zero(nbc,2);
	BCtype = unyque::IMatrix_zero(nbc,2);
	for (int i = 0; i < nbc; i++) {
	  fp >> ws >> BCtype(i,0) >> ws >> BCtype(i,1)	\
	     >> ws >> BCvals(i,0) >> ws >> BCvals(i,1);
	  if (BCtype(i,1) == 0) {
	    BCtype(i,1) = 2; // Treat no B.C. case as type 2 B.C.
	    BCvals(i,0) = BCvals(i,1) = 0.0;
	  }
	}
      }

      fp >> ws;
    }
    fp.close();
  } else
    cout << "Unable to open file" << endl;

}
//------------------------------------------------------------------------------
void ElEs::GenerateENC() {
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
void ElEs::MapDOFs() {

  fem::FEM_Element *el;
  fem::FEM_Edge *ed;
  int bcno;

  L2G_Phi = unyque::IVector_zero(nnode);
  L2G_BPhi = unyque::IVector_zero(nbedge);
  L2G_BdPhidn = unyque::IVector_zero(nbedge);

  // First identify nodes that are not in conducting regions
  for (int i = 0; i < nelem; i++) {
    el = s->Elements[i+1];
    if (RegAttr(el->reg - 1,0) > 0.0) { // If region is not a conductor
      for (int j = 1; j <= enode; j++) { // Loop over nodes in that element
	L2G_Phi(ENC(i,j)-1) = 1;
      }
    }
  }

  // Check B.C.'s on boundary edges
  for (int i = 1; i <= nbedge; i++) {

    ed = s->BEdges[i];

    bcno = -1;
    // Find the b.c. no: corresponding to the edge's marker (Default -1)
    for (int j = 0; j < nbc; j++) {
      if (BCtype(j,0) == ed->bmarker) {bcno = j; break;}
    }

    // Perform action based on b.c. type
    // Discard Phi dofs for boundary types 1 & 3, add BPhi dofs for boundary
    // type 2 and BdPhidn dofs for all boundary types
    switch (BCtype(bcno,1)) {
    case 1: // Dirichlet boundary - Remove nodes from global set
      L2G_Phi(ed->node1 - 1) = 0;
      L2G_Phi(ed->node2 - 1) = 0;
      L2G_Phi(ed->node3 - 1) = 0;
      L2G_BdPhidn(i-1) = 1;
      break;
    case 2: // Semiconductor-air interface
      L2G_BPhi(i-1) = 1;
      L2G_BdPhidn(i-1) = 1;
      break;
    case 3: // Metal-air interface
      L2G_Phi(ed->node1 - 1) = 0;
      L2G_Phi(ed->node2 - 1) = 0;
      L2G_Phi(ed->node3 - 1) = 0;
      L2G_BdPhidn(i-1) = 1;
      break;
    }

  }

  nPhi = 0, nBPhi = 0, nBdPhidn = 0;

  // Count number of Phi dofs that are included
  for (int i = 1; i <= nnode; i++) {
    nPhi += L2G_Phi(i-1);
    L2G_Phi(i-1) *= nPhi; // Store position in global matrix
  }

  // Count number of BPhi dofs that are included
  for (int i = 1; i <= nbedge; i++) {
    nBPhi += L2G_BPhi(i-1);
    L2G_BPhi(i-1) *= (nPhi+nBPhi);
  }

  // Count number of BdPhidn dofs that are included
  for (int i = 1; i <= nbedge; i++) {
    nBdPhidn += L2G_BdPhidn(i-1);
    L2G_BdPhidn(i-1) *= (nPhi+nBPhi+nBdPhidn);
  }

}
//------------------------------------------------------------------------------
void ElEs::SolveStatic() {
  int iter = 0;
  double err;
  ApplyInhomogeneousDBC();
  K.resize(nPhi+nBPhi+nBdPhidn+1,nPhi+nBPhi+nBdPhidn+1);
  RHS = unyque::DVector_zero(nPhi+nBPhi+nBdPhidn+1);
  dU = unyque::DVector_zero(nPhi+nBPhi+nBdPhidn+1);
  CompKandRHS();
  ApplyBC();
  dU = unyque::umfpackSolve(K,RHS);
  err = ublas::norm_inf(dU);
  UpdateGlobalPotentials();
  iter++;
  if (c->DEBUG) cout<<"iteration: "<<iter<<"  Error: "<<err<<endl;
  if (c->DEBUG) cout<<"Potential at infinity: "<<c->Phi_mult*c->Phi_inf<<endl;
  PrintResults();
}
//------------------------------------------------------------------------------
void ElEs::ApplyInhomogeneousDBC() {

  int bcno;
  fem::FEM_Edge *ed;
  fem::FEM_Element *el;

  // Loop over boundary edges
  for (int eid = 0; eid < nbedge; eid++) {

    ed = s->BEdges[eid+1];
    bcno = -1;
    // Find the b.c. no: corresponding to the edge's marker (Default -1)
    for (int i = 0; i < nbc; i++) {
      if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
    }

    // If it is a Dirichlet b.c. then initialize the corresponding nodes in the
    // solution vector with the value of the prescribed potential
    if (BCtype(bcno,1) == 1) {
      (s->Phi)(ed->node1 - 1) = BCvals(bcno,0);
      (s->Phi)(ed->node2 - 1) = BCvals(bcno,0);
      (s->Phi)(ed->node3 - 1) = BCvals(bcno,0);
      (s->BPhi)(eid) = BCvals(bcno,0);
    }
    if (BCtype(bcno,1) == 3) { // Conductor
      el = s->Elements[ed->eno];
      (s->BPhi)(eid) = RegAttr(el->reg - 1,1);
    }

  } // End of loop over boundary edges

  // Loop over elements
  for (int eid = 0; eid < nelem; eid++) {
    el = s->Elements[eid+1];
    if (RegAttr(el->reg - 1,0) == 0.0) { // If region is a conductor
      (s->Phi)(el->node1 - 1) = RegAttr(el->reg - 1,1);
      (s->Phi)(el->node2 - 1) = RegAttr(el->reg - 1,1);
      (s->Phi)(el->node3 - 1) = RegAttr(el->reg - 1,1);
      (s->Phi)(el->node4 - 1) = RegAttr(el->reg - 1,1);
      (s->Phi)(el->node5 - 1) = RegAttr(el->reg - 1,1);
      (s->Phi)(el->node6 - 1) = RegAttr(el->reg - 1,1);
    }
  } // End of loop over elements

}
//------------------------------------------------------------------------------
void ElEs::CompKandRHS() {

  int Gnode, iglobal, jglobal;
  double si, ti, T, bphi, phi_inf, bdphidn;
  unyque::DMatrix Ecoor, dN, B, Finv, FB, BFFB, Ke;
  unyque::DVector Te, Phie, RHSe, N, integrals;
  fem::FEM_Element *el;
  fem::FEM_Edge *ed;

  Ecoor = unyque::DMatrix_zero(enode,2);
  Te = unyque::DVector_zero(enode);
  Phie = unyque::DVector_zero(enode);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);

  // Loop over elements
  for (int eid = 0; eid < nelem; eid++) {

    el = s->Elements[eid+1];
    if (RegAttr(el->reg - 1,0) != 0.0) { // If region is not a conductor
      // Loop over nodes in that element
      for (int i = 0; i < enode; i++) {
	Gnode = ENC(eid,i+1); // Global node corresponding to (i+1)th local node
	Ecoor(i,0) = s->Nodes[Gnode]->x; // X-coordinate at that node
	Ecoor(i,1) = s->Nodes[Gnode]->y; // Y-coordinate at that node
	Te(i) = (s->T)(Gnode - 1); // Get the temperature at that node
	Phie(i) = (s->Phi)(Gnode - 1); // Get the potential at that node
      }

      // Loop over the Gauss points
      Ke = unyque::DMatrix_zero(enode,enode);
      RHSe = unyque::DVector_zero(enode);
      for (int ip = 0; ip < Gnquad; ip++) {
	si = Gs(ip); ti = Gt(ip);
	CompN(si, ti, N, dN);
	B = CompJandB(dN, Ecoor);
	Finv = CompF(eid, B);

	// Compute FB = [Finv]'[B]
	FB = ublas::prod(ublas::trans(Finv), B);

	// Compute the value of [B]'[Finv][Finv]'[B] at integration point ip
	BFFB = ublas::prod(ublas::trans(FB), FB);

	// Compute the temperature at the current integration point
	T = ublas::inner_prod(N, Te);

	// Evaluate the electrical conductivity at this temperature
	EC = c->functions->Eval_EC(T);

	// Compute elemental stiffness matrix by integrating  [BFFB]*TC*detF
	Ke += EC*BFFB*detF*Gw(ip)*detJ;

	// Create the elemental vector RHSe
	RHSe -= EC*ublas::prod(BFFB, Phie)*detF*Gw(ip)*detJ;

      } // End of loop over integration points

      // Assemble the elemental matrices Ke & RHSe to K and RHS respectively
      for (int i = 0; i < enode; i++) {
	iglobal = L2G_Phi(ENC(eid, i+1) - 1);
	if (iglobal > 0) {
	  for (int j = 0; j < enode; j++) {
	    jglobal = L2G_Phi(ENC(eid, j+1) - 1);
	    if (jglobal > 0)
	      K.append_element(iglobal-1,jglobal-1,Ke(i,j));
	  }
	  RHS(iglobal-1) += RHSe(i);
	}
      }

    } // End of loop over elements
  }

  phi_inf = c->Phi_inf; // Get the BEM potential at the bounding surface

  // Loop over boundary edges
  for (int i = 0; i < nbedge; i++) {

    // Get the edge corresponding to source point
    ed = s->BEdges[i+1];

    // Compute the temperature at the source point
    T = (s->T)(ed->node3 - 1);

    // Evaluate the electrical conductivity at this temperature
    EC = c->functions->Eval_EC(T);

    // Assemble first and last terms of BEM equation
    iglobal = L2G_BdPhidn(i); jglobal = L2G_BPhi(i);
    if (jglobal > 0)
      K.append_element(iglobal-1, jglobal-1, 0.5);
    K.append_element(iglobal-1, nPhi+nBPhi+nBdPhidn, -1.0);
    bphi = (s->BPhi)(i); // Get the BEM potential at ith edge
    RHS(iglobal-1) = -0.5*bphi+phi_inf;

    for (int j = 0; j < nbedge; j++) {

      // Compute boundary integrals
      integrals = CompIGGstar2(i, j);

      // Assemble terms of BEM equation
      jglobal = L2G_BPhi(j);
      if (jglobal > 0)
	K.append_element(iglobal-1, jglobal-1, -integrals(1));
      jglobal = L2G_BdPhidn(j);
      K.append_element(iglobal-1, jglobal-1, integrals(0));
      bphi = (s->BPhi)(j); // Get the BEM potential at jth edge
      bdphidn = (s->BdPhidn)(j); // Derivative of BEM potential at jth edge
      RHS(iglobal-1) -= -bphi*integrals(1) + bdphidn*integrals(0);

      // Assemble terms of charge conservation equation
      if (j == i) {
	K.append_element(nPhi+nBPhi+nBdPhidn, L2G_BdPhidn(j)-1,
			 (1-EC0/EC)*integrals(2));
	RHS(nPhi+nBPhi+nBdPhidn) -= bdphidn*(1-EC0/EC)*integrals(2);
	Lengths(i) = integrals(2);
      }

    }
  }

}
//------------------------------------------------------------------------------
void ElEs::CompN(double s, double t, unyque::DVector &N, unyque::DMatrix &dN) {
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
unyque::DMatrix ElEs::CompJandB(unyque::DMatrix &dN, unyque::DMatrix &Ecoor) {

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
  return ublas::prod(InvJ, dN);

}
//------------------------------------------------------------------------------
unyque::DMatrix ElEs::CompF(int eid, unyque::DMatrix &B) {
  int Gnode;
  double detFinv;
  unyque::DMatrix Ue, F, Finv;

  Ue = unyque::DMatrix_zero(enode, 2);
  F = unyque::DMatrix_zero(2,2);
  Finv = unyque::DMatrix_zero(2,2);

  for (int i = 0; i < enode; i++) {
    Gnode = ENC(eid,i+1); // Global node corresponding to (i+1)th local node
    Ue(i, 0) = (s->U)(Gnode-1);
    Ue(i, 1) = (s->V)(Gnode-1);
  }

  F = unyque::DMatrix_identity(2) + ublas::trans(ublas::prod(B, Ue));

  // To disable total Lagrangian formulation, uncomment the line below
  // F = unyque::DMatrix_identity(2);

  detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

  detFinv = 1.0/detF;

  Finv(0,0) = detFinv*F(1,1); Finv(0,1) = -detFinv*F(0,1);
  Finv(1,0) = -detFinv*F(1,0); Finv(1,1) = detFinv*F(0,0);

  return Finv;

}
//------------------------------------------------------------------------------
unyque::DVector ElEs::CompIGGstar(int p, int q) {

  fem::FEM_Edge *edp, *edq;
  double x1, y1, x2, y2, L, d1, d2, nx, ny;
  unyque::DMatrix M;
  unyque::DVector integrals;

  // Get the edges corresponding to source & field point
  edp = s->BEdges[p+1];
  edq = s->BEdges[q+1];

  // Find the length of the edge corresponding to the field point
  x1 = s->Nodes[edq->node1]->x;
  y1 = s->Nodes[edq->node1]->y;
  x2 = s->Nodes[edq->node2]->x;
  y2 = s->Nodes[edq->node2]->y;
  L = abs(sqrt(pow((x1-x2),2)+pow((y1-y2),2)))/2;

  // Get coordinates of the midpoints of both edges
  x1 = s->Nodes[edp->node3]->x;
  y1 = s->Nodes[edp->node3]->y;
  x2 = s->Nodes[edq->node3]->x;
  y2 = s->Nodes[edq->node3]->y;

  // Compute the transformation matrix corr. to the rotation
  nx = -cos(edq->normal); // cos(Theta)
  ny = -sin(edq->normal); // sin(Theta)
  M = unyque::DMatrix_zero(2,2);
  M(0,0) = -ny; // cos(-PI/2-Theta) = -sin(PI/2)sin(Theta) = -sin(Theta)
  M(0,1) = nx;  // -sin(-PI/2-Theta) = sin(PI/2)cos(Theta) = cos(Theta)
  M(1,0) = -nx; // sin(-PI/2-Theta) = -sin(PI/2)cos(Theta) = -cos(Theta)
  M(1,1) = -ny; // cos(-PI/2-Theta) = -sin(PI/2)sin(Theta) = -sin(Theta)
  d1 = M(0,0)*(x2-x1)+M(0,1)*(y2-y1);
  d2 = M(1,0)*(x2-x1)+M(1,1)*(y2-y1);

  if (p != q) {
    integrals(0) = -(L/(4*4*atan2(1.0,1.0)))*			   \
      ( 4*(log(L)-1) + (1-d1/L)*log(pow((1-d1/L),2)+pow(d2/L,2)) + \
	(1+d1/L)*log(pow((1+d1/L),2)+pow(d2/L,2)) + \
	(2*d2/L)*atan2(2*L*d2,pow(d1,2)+pow(d2,2)-pow(L,2)) );
    integrals(1) = -(1/(2*4*atan2(1.0,1.0)))*
      atan2(2*L*d2,pow(d1,2)+pow(d2,2)-pow(L,2));
  }
  else {
    integrals(0) = -(L/(4*atan2(1.0,1.0)))*(log(L)-1);
    integrals(1) = 0.0;
  }

  integrals(2) = 2*L;

  return integrals;

}
//------------------------------------------------------------------------------
unyque::DVector ElEs::CompIGGstar2(int p, int q) {

  fem::FEM_Edge *edp, *edq;
  int Gnode, startip;
  double Rsq, si, ti, detJ1D = 0.0, FinvTN;
  unyque::DMatrix Ecoor, FF, dN, B, Finv;
  unyque::DVector N, integrals, n(2), xp(2), up(2), xq(2), uq(2);

  // Get the edges corresponding to field & source points respectively
  edp = s->BEdges[p+1];
  edq = s->BEdges[q+1];

  Ecoor = unyque::DMatrix_zero(enode,2);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);
  integrals = unyque::DVector_zero(3);

  // Compute components of normal and position coordinates of source point
  n(0) = -cos(edq->normal);
  n(1) = -sin(edq->normal);

  for (int i = 0; i < enode; i++) {
    Gnode = ENC(edq->eno-1, i+1); // Global node corr. to (i+1)th local node
    Ecoor(i,0) = s->Nodes[Gnode]->x; // X-coordinate at that node
    Ecoor(i,1) = s->Nodes[Gnode]->y; // Y-coordinate at that node
  }

  // Compute coordinates of position and displacement vectors of field point
  xp(0) = s->Nodes[edp->node3]->x;
  xp(1) = s->Nodes[edp->node3]->y;
  up(0) = (s->U)(edp->node3 - 1);
  up(1) = (s->V)(edp->node3 - 1);

  startip = Genquad*(edq->eid);
  // Loop over the  Gauss points for this edge
  for (int ip = startip; ip < (startip + Genquad); ip++) {

    // Compute shape functions and lengths of edges
    si = Ges(ip); ti = Get(ip);
    CompN(si, ti, N, dN);
    B = CompJandB(dN, Ecoor);
    Finv = CompF(edq->eno-1, B);

    // Compute FF = [Finv][Finv]'
    FF = ublas::prod(Finv, ublas::trans(Finv));

    // Compute |[Finv]'*N| where N is outward normal in undeformed configuration
    FinvTN = abs(ublas::norm_2(ublas::prod(ublas::trans(Finv), n)));

    // Compute coordinates of source point and square of distance to field point
    xq = ublas::prod(N, Ecoor);
    uq(0) = 0; uq(1) = 0;
    for (int i = 0; i < enode; i++) {
      Gnode = ENC(edq->eno-1, i+1); // Global node corr. to (i+1)th local node
      uq(0) += N(i)*(s->U)(Gnode - 1);
      uq(1) += N(i)*(s->V)(Gnode - 1);
    }
    Rsq = pow(ublas::norm_2((xq+uq)-(xp+up)), 2);

    // Find the length of the edge
    switch (edq->eid) {
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

    integrals(0) += (-1/(4*4*atan(1.0)))*log(Rsq)*detF*FinvTN*Gew(ip)*detJ1D;
    integrals(1) += (1/(2*4*atan(1.0)*Rsq))*
      ublas::inner_prod(xq-xp, ublas::prod(FF, n))*detF*Gew(ip)*detJ1D;
    integrals(2) += detF*FinvTN*Gew(ip)*detJ1D;

  }

  return integrals;

}
//------------------------------------------------------------------------------
void ElEs::ApplyBC() {
  int bcno;
  fem::FEM_Edge *ed;
  // fem::FEM_Element *el;

  // Loop over boundary edges - Apply interface BCs first
  for (int eid = 0; eid < nbedge; eid++) {

    ed = s->BEdges[eid+1];
    bcno = -1;
    // Find the b.c. no: corresponding to the edge's marker (Default -1)
    for (int i = 0; i < nbc; i++) {
      if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
    }

    // Check if b.c. type is 2 and apply interface BC if it is
    if (BCtype(bcno,1) == 2)
      ApplyIBC(ed);

  }  // End of loop over boundary edges

  // // Loop over boundary edges - Apply DBCs finally
  // for (int eid = 0; eid < nbedge; eid++) {

  //   ed = s->BEdges[eid+1];
  //   bcno = -1;
  //   // Find the b.c. no: corresponding to the edge's marker (Default -1)
  //   for (int i = 0; i < nbc; i++) {
  //     if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
  //   }

  //   // Check if b.c. type is 1 or 3 and apply DBC if it is
  //   if (BCtype(bcno,1) == 1) {
  //     ApplyDBC(ed->node1 - 1);
  //     ApplyDBC(ed->node2 - 1);
  //     ApplyDBC(ed->node3 - 1);
  //     ApplyDBC(nnode + 2*(ed->id - 1));
  //   } else if (BCtype(bcno,1) == 3)
  //     ApplyDBC(nnode + 2*(ed->id - 1));

  // }  // End of loop over boundary edges

  // // Loop over elements
  // for (int eid = 0; eid < nelem; eid++) {
  //   el = s->Elements[eid+1];
  //   if (RegAttr(el->reg - 1,0) == 0.0) { // If region is a conductor
  //     ApplyDBC(el->node1 - 1);
  //     ApplyDBC(el->node2 - 1);
  //     ApplyDBC(el->node3 - 1);
  //     ApplyDBC(el->node4 - 1);
  //     ApplyDBC(el->node5 - 1);
  //     ApplyDBC(el->node6 - 1);
  //   }
  // } // End of loop over elements

}
//------------------------------------------------------------------------------
void ElEs::ApplyIBC(fem::FEM_Edge *ed) {

  int Gnode, startip, iglobal, jglobal;
  double si, ti, phi, bPhi, bdPhidn, detJ1D = 0, FinvTN;
  unyque::DMatrix Ecoor, dN, B, Finv;
  unyque::DVector Phie, N, n(2);

  Ecoor = unyque::DMatrix_zero(enode,2);
  Phie = unyque::DVector_zero(enode);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);

  for (int i = 0; i < enode; i++) {
    Gnode = ENC(ed->eno-1, i+1); // Global node corr. to (i+1)th local node
    Ecoor(i,0) = s->Nodes[Gnode]->x; // X-coordinate at that node
    Ecoor(i,1) = s->Nodes[Gnode]->y; // Y-coordinate at that node
    Phie(i) = (s->Phi)(Gnode - 1); // Get the potential at that node
  }
  bPhi = (s->BPhi)(ed->id - 1); // Get BEM potential at this edge
  bdPhidn = (s->BdPhidn)(ed->id - 1); // Get derivative of BEM potential

  // Compute components of normal vector
  n(0) = -cos(ed->normal);
  n(1) = -sin(ed->normal);

  startip = Genquad*(ed->eid);
  // Loop over the  Gauss points for this edge
  for (int ip = startip; ip < (startip + Genquad); ip++) {

    // Compute shape functions and lengths of edges
    si = Ges(ip); ti = Get(ip);
    CompN(si, ti, N, dN);
    B = CompJandB(dN, Ecoor);
    Finv = CompF(ed->eno-1, B);

    // Compute the FEM potential at the current integration point
    phi = ublas::inner_prod(N, Phie);

    // Compute |[Finv]'*N| where N is outward normal in undeformed configuration
    FinvTN = abs(ublas::norm_2(ublas::prod(ublas::trans(Finv), n)));

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

    for (int i = 0; i < enode; i++) {
      iglobal = L2G_Phi(ENC(ed->eno-1, i+1) - 1);
      if (iglobal > 0) {
	jglobal = L2G_BdPhidn(ed->id-1);
	K.append_element(iglobal-1, jglobal-1,
		-N(i)*EC0*detF*FinvTN*Gew(ip)*detJ1D); // Continuity of Jn
	RHS(iglobal-1) +=
	  N(i)*EC0*bdPhidn*detF*FinvTN*Gew(ip)*detJ1D; // Continuity of Jn
	jglobal = L2G_BPhi(ed->id-1);
	if (jglobal > 0)
	  K.append_element(jglobal-1, iglobal-1,
		  -N(i)*detF*FinvTN*Gew(ip)*detJ1D); // Continuity of Phi
      }
    }
    iglobal = L2G_BPhi(ed->id-1);
    if (iglobal > 0) {
      K.append_element(iglobal-1, iglobal-1,
	      detF*FinvTN*Gew(ip)*detJ1D); // Continuity of Phi
      RHS(iglobal-1) -=
	(bPhi - phi)*detF*FinvTN*Gew(ip)*detJ1D; // Continuity of Phi
    }
  }

}
//------------------------------------------------------------------------------
void ElEs::ApplyDBC(int nid) {

  // Zero out the row 'nid' and column 'nid' in K corresponding to the node
  for (int i = 0; i < (nnode+2*nbedge+1); i++) {
    if (K(nid, i) != 0.0)
      K.erase_element(nid, i);
    if (K(i, nid) != 0.0)
      K.erase_element(i, nid);
  }

  // Set K(nid,nid) as one
  K(nid,nid) = 1.0;
  // Set RHS(nid) to the value of the Dirichlet boundary condition
  RHS(nid) = 0.0;

}
//------------------------------------------------------------------------------
void ElEs::UpdateGlobalPotentials() {

  for (int i = 0; i < nnode; i++)
    if (L2G_Phi(i) > 0)
      (s->Phi)(i) += dU(L2G_Phi(i)-1);

  for (int i = 0; i < nbedge; i++) {
    if (L2G_BPhi(i) > 0)
      (s->BPhi)(i) += dU(L2G_BPhi(i)-1);
    if (L2G_BdPhidn(i) > 0)
      (s->BdPhidn)(i) += dU(L2G_BdPhidn(i)-1);
  }

  c->Phi_inf += dU(nPhi+nBPhi+nBdPhidn);

}
//------------------------------------------------------------------------------
void ElEs::PrintResults() {
  FILE *fp = NULL, *fp2 = NULL;
  fem::FEM_Point *pp;
  fem::FEM_Edge *ed;
  unyque::DMatrix Ecoor, dN, B, Finv;
  unyque::DVector Phie, N, nX(2), nx(2), tx(2), EX(2), Ex(2);
  int Gnode;
  double si = 0, ti = 0;
  double FinvTN, T, xval, yval;

  if (c->DEBUG) {
    fp = fopen("postproc/phi.dat","w");
    for (int i = 0; i < nnode; i++) {
      pp = s->Nodes[i+1];
      xval = pp->x + (s->U)(i);
      yval = pp->y + (s->V)(i);
      fprintf(fp,"%.14lf  %.14lf %.14lf \n",xval,yval,
	      c->Phi_mult*(s->Phi)(i));
    }
    fclose(fp);

    fp = fopen("postproc/bphi.dat","w");
    fp2 = fopen("postproc/scharge.dat","w");
  }
  for (int i = 0; i < nbedge; i++) {

    // Find the edge we're currently on
    ed = s->BEdges[i+1];

    // Store the surface charge density and components internal electric field
    if (RegAttr((s->Elements[ed->eno])->reg - 1,0) == 0.0) { // Metal-air
      (s->SCharge)(i) = -(s->BdPhidn)(i)*EPS0;
      (s->Ent)(i,0) = 0.0; (s->Ent)(i,1) = 0.0; // Field zero inside metal
    }
    else { // Semiconductor-air interface

      T = (s->T)(ed->node3 - 1); // Temperature at this point
      EC = c->functions->Eval_EC(T); // Electrical conductivity at this temp
      (s->SCharge)(i) = (s->BdPhidn)(i)*(EC0/EC-1)*EPS0;

      Ecoor = unyque::DMatrix_zero(enode,2);
      Phie = unyque::DVector_zero(enode);
      N = unyque::DVector_zero(enode);
      dN = unyque::DMatrix_zero(2, enode);

      // Find x & y coordinates & value of electrostatic potential at each node
      for (int j = 0; j < enode; j++) {
	Gnode = ENC(ed->eno-1, j+1); // Global node corr. to (i+1)th local node
	Ecoor(j,0) = s->Nodes[Gnode]->x;
	Ecoor(j,1) = s->Nodes[Gnode]->y;
	Phie(j) = (s->Phi)(Gnode-1);
      }

      // Compute shape functions, Jacobian and deformation gradient
      switch (ed->eid) {
      case 0:
	si = 0.5; ti = 0.0;
	break;
      case 1:
	si = 0.5; ti = 0.5;
	break;
      case 2:
	si = 0.0; ti = 0.5;
	break;
      default:
	if (c->DEBUG) cout<<"Invalid edge ID"<<endl;
      }
      CompN(si, ti, N, dN);
      B = CompJandB(dN, Ecoor);
      Finv = CompF(ed->eno-1, B);

      // Find the components of normal & tangential vectors at mid-point of edge
      nX(0) = cos(ed->normal);
      nX(1) = sin(ed->normal);
      FinvTN = abs(ublas::norm_2(ublas::prod(ublas::trans(Finv), nX)));
      nx = ublas::prod(ublas::trans(Finv), nX)/FinvTN;
      tx(0) = -nx(1); tx(1) = nx(0);

      // Compute components of electric field at edge
      EX = -ublas::prod(B, Phie);
      Ex = ublas::prod(ublas::trans(Finv), EX);
      (s->Ent)(i,0) = ublas::inner_prod(Ex, nx);
      (s->Ent)(i,1) = ublas::inner_prod(Ex, tx);

    }

    // Compute coordinates of midpoint of edge for plotting purposes
    pp = s->Nodes[ed->node3];
    xval = pp->x + (s->U)(ed->node3 - 1);
    yval = pp->y + (s->V)(ed->node3 - 1);

    if (c->DEBUG) {
      fprintf(fp,"%d\t%e\t%e\t%e\n",ed->id,xval,yval,
	      c->Phi_mult*(s->BPhi)(i));
      fprintf(fp2,"%d\t%e\t%e\t%e\t%e\t%e\t%e\n", ed->id, xval, yval, nX(0),
	      nX(1), Lengths(i), c->Phi_mult*1e6*(s->SCharge)(i));
    }
  }
  if (c->DEBUG) {
    fclose(fp);
    fclose(fp2);
  }
}
//------------------------------------------------------------------------------
void ElEs::GaussLegendre(double x1, double x2, double* x, double* w, int n) {
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;

  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++) {
    z=cos(3.141592653589798*(i-0.25)/(n+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > 1e-8);
    x[i-1]=xm-xl*z;
    x[n-i]=xm+xl*z;
    w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n-i]=w[i-1];
  }
}
