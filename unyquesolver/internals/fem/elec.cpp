#include "elec.hpp"

//------------------------------------------------------------------------------
Elec::Elec() {
}
//------------------------------------------------------------------------------
Elec::Elec(boost::shared_ptr<fem::FEM_PhysicalDomain> is,
	   boost::shared_ptr<fem::FEM_Common> ic) {
  s = is;
  c = ic;
  nelem = is->nelem;
  nnode = is->nnode; nbnode = is->nbnode;
  nedge = is->nedge; nbedge = is->nbedge;
  enode = 6;
}
//------------------------------------------------------------------------------
void Elec::Init() {
  char filename[100];

  nelem = s->nelem;
  nnode = s->nnode; nbnode = s->nbnode;
  nedge = s->nedge; nbedge = s->nbedge;

  ENC = unyque::IMatrix_zero(nelem,7);

  // Initialize and set the integration points and weights
  Init_GIntegration();

  EC = 0.0; Qj = 0.0;
  sprintf(filename,"./conf/electrical.conf");
  ReadElec(filename);
  GenerateENC();
  MapDOFs();
}
//------------------------------------------------------------------------------
void Elec::Init_GIntegration() {

  Genquad = 2; // Number of integration points for 1D edge integration
  Gnquad = 4;  // Number of integration points for 2D element integration

  Ges = unyque::DVector_zero(3*Genquad);
  Get = unyque::DVector_zero(3*Genquad);
  Gew = unyque::DVector_zero(3*Genquad);

  // Setting the integration point for 1D edge integration
  switch (Genquad) {
  case 2:
    Ges(0) = 0.211324865405187; Get(0) = 0.0;
    Ges(1) = 0.788675134594813; Get(1) = 0.0;
    Ges(2) = 0.788675134594813; Get(2) = 0.211324865405187;
    Ges(3) = 0.211324865405187; Get(3) = 0.788675134594813;
    Ges(4) = 0.0;               Get(4) = 0.788675134594813;
    Ges(5) = 0.0;               Get(5) = 0.211324865405187;
    Gew = unyque::DVector_scalar(3*Genquad, 0.5);
    break;

  case 3:
    Ges(0) = 0.112701665; Get(0) = 0.0;
    Ges(1) = 0.500000000; Get(1) = 0.0;
    Ges(2) = 0.887298335; Get(2) = 0.0;
    Ges(3) = 0.887298335; Get(3) = 0.112701665;
    Ges(4) = 0.500000000; Get(4) = 0.500000000;
    Ges(5) = 0.112701665; Get(5) = 0.887298335;
    Ges(6) = 0.0;         Get(6) = 0.887298335;
    Ges(7) = 0.0;         Get(7) = 0.500000000;
    Ges(8) = 0.0;         Get(8) = 0.112701665;

    Gew(0) = 0.277777778;
    Gew(1) = 0.444444444;
    Gew(2) = 0.277777778;
    Gew(3) = 0.277777778;
    Gew(4) = 0.444444444;
    Gew(5) = 0.277777778;
    Gew(6) = 0.277777778;
    Gew(7) = 0.444444444;
    Gew(8) = 0.277777778;
    break;

  case 4:
    Ges(0)  = 0.069431845; Get(0)  = 0.0;
    Ges(1)  = 0.330009480; Get(1)  = 0.0;
    Ges(2)  = 0.669990520; Get(2)  = 0.0;
    Ges(3)  = 0.930568155; Get(3)  = 0.0;

    Ges(4)  = 0.930568155; Get(4)  = 0.069431845;
    Ges(5)  = 0.669990520; Get(5)  = 0.330009480;
    Ges(6)  = 0.330009480; Get(6)  = 0.669990520;
    Ges(7)  = 0.069431845; Get(7)  = 0.930568155;

    Ges(8)  = 0.0;         Get(8)  = 0.930568155;
    Ges(9)  = 0.0;         Get(9)  = 0.669990520;
    Ges(10) = 0.0;         Get(10) = 0.330009480;
    Ges(11) = 0.0;         Get(11) = 0.069431845;

    Gew(0)  = 0.173927425;
    Gew(1)  = 0.326072575;
    Gew(2)  = 0.326072575;
    Gew(3)  = 0.173927425;
    Gew(4)  = 0.173927425;
    Gew(5)  = 0.326072575;
    Gew(6)  = 0.326072575;
    Gew(7)  = 0.173927425;
    Gew(8)  = 0.173927425;
    Gew(9)  = 0.326072575;
    Gew(10) = 0.326072575;
    Gew(11) = 0.173927425;
    break;
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
void Elec::ReadElec(char *filename) {

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

      if (key.compare("Qj") == 0)
	fp >> ws >> Qj;
      if (key.compare("PHI_MULT") == 0)
	fp >> ws >> c->Phi_mult;
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
void Elec::GenerateENC() {
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
void Elec::MapDOFs() {

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
    L2G(i) *= ndof; // Store position in global matrix
  }

}
//------------------------------------------------------------------------------
void Elec::SolveStatic() {
  int iter = 0;
  double err = 1.0, eps = 1e-10;

  ApplyInhomogeneousDBC();
  do {
    K.resize(ndof,ndof);
    RHS = unyque::DVector_zero(ndof);
    dU = unyque::DVector_zero(ndof);
    CompDomIntegrals();
    ApplyBC();
    dU = unyque::umfpackSolve(K,RHS);
    err = ublas::norm_inf(dU);
    UpdateGlobalPotential();
    iter++;
    if (c->DEBUG) cout<<"iteration: "<<iter<<"  Error: "<<err<<endl;
  } while (err>eps);
  if (c->DEBUG) PrintResults();
}
//------------------------------------------------------------------------------
void Elec::ApplyInhomogeneousDBC() {

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
    // solution vector with the value of the prescribed potential
    if (BCtype(bcno,1) == 1) {
      (s->Phi)(ed->node1 - 1) = BCvals(bcno,0);
      (s->Phi)(ed->node2 - 1) = BCvals(bcno,0);
      (s->Phi)(ed->node3 - 1) = BCvals(bcno,0);
    }

  } // End of loop over boundary edges

}
//------------------------------------------------------------------------------
void Elec::CompDomIntegrals() {

  int Gnode, iglobal, jglobal;
  double si, ti, T;
  unyque::DMatrix Ecoor, dN, B, Finv, FB, BFFB, Ke;
  unyque::DVector Te, Phie, RHSe, N;

  Ecoor = unyque::DMatrix_zero(enode,2);
  Te = unyque::DVector_zero(enode);
  Phie = unyque::DVector_zero(enode);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);

  // Loop over elements
  for (int eid = 0; eid < nelem; eid++) {

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

      // Compute the value of the integrand: [BFFB]*TC*detF
      Ke += EC*BFFB*detF*Gw(ip)*detJ;

      // Create the elemental vector RHSe
      RHSe += (Qj*N - EC*ublas::prod(BFFB, Phie))*detF*Gw(ip)*detJ;

    } // End of loop over integration points

    // Assemble the elemental matrices Ke and RHSe to K and RHS respectively
    for (int i = 0; i < enode; i++) {
      iglobal = L2G(ENC(eid, i+1) - 1);
      if (iglobal > 0) {
	for (int j = 0; j < enode; j++) {
	  jglobal = L2G(ENC(eid, j+1) - 1);
	  if (jglobal > 0)
	    K.append_element(iglobal-1,jglobal-1,Ke(i,j));
	}
	RHS(iglobal-1) += RHSe(i);
      }
    }

  } // End of loop over elements

}
//------------------------------------------------------------------------------
void Elec::CompN(double s, double t, unyque::DVector &N, unyque::DMatrix &dN) {
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
unyque::DMatrix Elec::CompJandB(unyque::DMatrix &dN, unyque::DMatrix &Ecoor) {

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
unyque::DMatrix Elec::CompF(int eid, unyque::DMatrix &B) {

  int Gnode;
  double detFinv;
  unyque::DMatrix Ue, F, Finv;

  Ue = unyque::DMatrix_zero(enode, 2);
  F = unyque::DMatrix_zero(2,2);
  Finv = unyque::DMatrix(2,2);

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

  Finv(0,0) = detFinv*F(1,1); Finv(0,1) =-detFinv*F(0,1);
  Finv(1,0) =-detFinv*F(1,0); Finv(1,1) = detFinv*F(0,0);

  return Finv;

}
//------------------------------------------------------------------------------
void Elec::ApplyBC() {
  int bcno;
  fem::FEM_Edge *ed;

  // Loop over boundary edges - Apply NBCs first
  for (int eid = 0; eid < nbedge; eid++) {

    ed = s->BEdges[eid+1];
    bcno = -1;
    // Find the b.c. no: corresponding to the edge's marker (Default -1)
    for (int i = 0; i < nbc; i++) {
      if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
    }

    // Check if b.c. type is 2 and apply NBC if it is
    if (BCtype(bcno,1) == 2)
      ApplyNBC(bcno, ed);

  } // End of loop over boundary edges

  // // Loop over boundary edges - Apply DBCs finally
  // for (int eid = 0; eid < nbedge; eid++) {

  //   ed = s->BEdges[eid+1];
  //   bcno = -1;
  //   // Find the b.c. no: corresponding to the edge's marker (Default -1)
  //   for (int i = 0; i < nbc; i++) {
  //     if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
  //   }

  //   // Check if b.c. type is 1 and apply DBC if it is
  //   if (BCtype(bcno,1) == 1) {
  //     ApplyDBC(ed->node1 - 1);
  //     ApplyDBC(ed->node2 - 1);
  //     ApplyDBC(ed->node3 - 1);
  //   }

  // } // End of loop over boundary edges
}
//------------------------------------------------------------------------------
void Elec::ApplyNBC(int bcno, fem::FEM_Edge *ed) {

  int Gnode, startip, iglobal;
  double si, ti, detJ1D = 0, FinvTN;
  unyque::DMatrix Ecoor, dN, B, Finv;
  unyque::DVector N, n(2);

  Ecoor = unyque::DMatrix_zero(enode,2);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);

  // Find the x & y coordinates of each node
  for (int i = 0; i < enode; i++) {
    Gnode = ENC(ed->eno-1, i+1); //global node corr. to (i+1)th local node
    Ecoor(i,0) = s->Nodes[Gnode]->x;
    Ecoor(i,1) = s->Nodes[Gnode]->y;
  }

  // Compute components of normal vector
  n(0) = cos(ed->normal);
  n(1) = sin(ed->normal);

  startip = Genquad*(ed->eid);
  // Loop over the Gauss points for this edge
  for (int ip = startip; ip < (startip + Genquad); ip++) {

    // Compute shape functions and lengths of edges
    si = Ges(ip); ti = Get(ip);
    CompN(si, ti, N, dN);
    B = CompJandB(dN, Ecoor);
    Finv = CompF(ed->eno-1, B);

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

    // Create the elemental vector RHSe
    for (int i = 0;i<enode;i++) {
      iglobal = L2G(ENC(ed->eno-1, i+1) - 1);
      if (iglobal > 0)
	RHS(iglobal-1) -= N(i)*BCvals(bcno,0)*FinvTN*detF*Gew(ip)*detJ1D;
    }
  }

}
//------------------------------------------------------------------------------
void Elec::ApplyDBC(int nid) {

  // Zero out the row 'nid' and column 'nid' in K corresponding to the node
  for (int i = 0; i < nnode; i++) {
    if (K(nid, i) != 0.0)
      K.erase_element(nid, i);
    if (K(i, nid) != 0.0)
      K.erase_element(i, nid);
  }

  // Set K(nid,nid) as one
  K(nid,nid) = 1.0;
  // Set RHS(nid) as the value of the Dirichlet boundary condition
  RHS(nid) = 0.0;

}
//------------------------------------------------------------------------------
void Elec::UpdateGlobalPotential() {
  for (int i = 0; i < nnode; i++)
    if (L2G(i) > 0)
      (s->Phi)(i) += dU(L2G(i) - 1);
}
//------------------------------------------------------------------------------
void Elec::PrintResults() {
  FILE *fp;
  fem::FEM_Point *pp;
  double xval, yval;

  fp = fopen("postproc/phi.dat","w");
  for (int i = 0; i < nnode; i++) {
    pp = s->Nodes[i+1];
    xval = pp->x + (s->U)(i);
    yval = pp->y + (s->V)(i);
    fprintf(fp,"%.14lf  %.14lf %.14lf \n", xval, yval,
	    c->Phi_mult*(s->Phi)(i));
  }
  fclose(fp);
}
//------------------------------------------------------------------------------
