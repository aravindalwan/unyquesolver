#include "fluid.hpp"

//------------------------------------------------------------------------------
Fluid::Fluid() {
}
//------------------------------------------------------------------------------
Fluid::Fluid(boost::shared_ptr<fem::FEM_PhysicalDomain> is,
	     boost::shared_ptr<fem::FEM_FluidDomain> isf,
	     boost::shared_ptr<fem::FEM_Common> ic) {
  s = is;
  sf = isf;
  c = ic;
  nelem = isf->nelem;
  nnode = isf->nnode; nbnode = isf->nbnode;
  nedge = isf->nedge; nbedge = isf->nbedge;
  enode = 6;
}
//------------------------------------------------------------------------------
void Fluid::Init() {
  char filename[100];

  nelem = sf->nelem;
  nnode = sf->nnode; nbnode = sf->nbnode;
  nedge = sf->nedge; nbedge = sf->nbedge;

  ENC = unyque::IMatrix_zero(nelem,7);

  // Initialize and set the integration points and weights
  Init_GIntegration();

  t = 0; dt = 0;
  sprintf(filename,"./conf/fluid.conf");
  ReadFluid(filename);
  GenerateENC();
  MapDOFs();
}
//------------------------------------------------------------------------------
void Fluid::Init_GIntegration() {
  Genquad = 2; // Number of integration points for 1D edge integration
  Gnquad = 4;  // Number of integration points for 2D element integration
  Gbnquad = 21; // Number of integration points for 1D integration along breadth

  Ges = unyque::DVector_zero(3*Genquad);
  Get = unyque::DVector_zero(3*Genquad);
  Gew = unyque::DVector_zero(3*Genquad);

  // Setting the integration points for 1D edge integration
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

  Gbs = unyque::DVector_zero(Gbnquad);
  Gbw = unyque::DVector_zero(Gbnquad);
  double gx[Gbnquad], gw[Gbnquad];

  ZMIN = numeric_limits<double>::max(); ZMAX = numeric_limits<double>::min();
  for (int i = 0; i < nnode; i++) {
    if (sf->Nodes[i + 1]->y < ZMIN) { ZMIN = sf->Nodes[i + 1]->y; }
    if (sf->Nodes[i + 1]->y > ZMAX) { ZMAX = sf->Nodes[i + 1]->y; }
  }
  unyque::GaussLegendre(ZMIN, ZMAX, gx, gw, Gbnquad);
  for (int i = 0; i < Gbnquad; i++) {
    Gbs(i) = gx[i]; Gbw(i) = gw[i];
  }

}
//------------------------------------------------------------------------------
void Fluid::ReadFluid(char *filename) {

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

      if (key.compare("ETA") == 0)
	fp >> ws >> ETA;
      if (key.compare("LAMBDA") == 0)
	fp >> ws >> LAMBDA;
      if (key.compare("P_ATM") == 0)
	fp >> ws >> P_ATM;
      if (key.compare("FIXED_EDGE") == 0)
	fp >> ws >> FIXED_EDGE;
      if (key.compare("MOVING_EDGE") == 0)
	fp >> ws >> MOVING_EDGE;
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
void Fluid::GenerateENC() {
  fem::FEM_Element *ee;
  for (int i = 1; i <= nelem; i++) {
    ee = sf->Elements[i];
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
void Fluid::MapDOFs() {

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

    ed = sf->BEdges[eid+1];
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
void Fluid::MapPhysicalToFluid() {

  fem::FEM_Edge *ed;
  unyque::DVector N, v0, v1, v2;
  double x, x1, x2, x3, z, si, ti;
  int node;

  N = unyque::DVector_zero(3);
  Node2MovingEdge = unyque::IVector_zero(nnode);
  Node2FixedEdge = unyque::IVector_zero(nnode);

  for (int i = 0; i < nnode; i++) {
    for (int eid = 0; eid < s->nbedge; eid++) {
      ed = s->BEdges[eid+1];

      if ((ed->bmarker == MOVING_EDGE) || (ed->bmarker == FIXED_EDGE)) {

	// Store x-coordinates of current node and edge nodes
	x = sf->Nodes[i+1]->x;
	x1 = s->Nodes[ed->node1]->x;
	x2 = s->Nodes[ed->node2]->x;
	x3 = s->Nodes[ed->node3]->x;

	// Compute basis functions
	N(0) = (x - x2)*(x - x3)/(x1 - x2)/(x1 - x3);
	N(1) = (x - x1)*(x - x3)/(x2 - x1)/(x2 - x3);
	N(2) = (x - x1)*(x - x2)/(x3 - x1)/(x3 - x2);

      }

      // Check if edge has the same ID as MOVING_EDGE
      if (ed->bmarker == MOVING_EDGE) {
	if ((x - x1)*(x - x2) <= 0) // If node belongs to this edge
	  Node2MovingEdge(i) = eid;
      }

      // Check if edge has the same ID as FIXED_EDGE
      if (ed->bmarker == FIXED_EDGE) {
	if ((x - x1)*(x - x2) <= 0) // If node belongs to this edge
	  Node2FixedEdge(i) = eid;
      }

    } // End loop over edges

  } // End loop over nodes

  v0 = unyque::DVector(2);
  v1 = unyque::DVector(2);
  v2 = unyque::DVector(2);
  IntegrationPoint2Element = unyque::IMatrix(s->nnode, Gbnquad);
  std::fill(IntegrationPoint2Element.data().begin(),
	    IntegrationPoint2Element.data().end(), -1);

  // Loop over edges in physical domain
  for (int eid = 0; eid < s->nbedge; eid++) {
    ed = s->BEdges[eid+1];

    if (ed->bmarker == MOVING_EDGE) {
      for (int i = 0; i < 3; i++) { // Loop over the nodes in the edge

	switch (i) {
	case 0:
	  node = ed->node1 - 1;
	  break;
	case 1:
	  node = ed->node2 - 1;
	  break;
	case 2:
	  node = ed->node3 - 1;
	}

	// If we have already mapped elements for this node, then ignore it
	if (IntegrationPoint2Element(node, 0) >= 0)
	  continue;

	// X-coordinate of this node
	x = s->Nodes[node + 1]->x;

	for (int j = 0; j < Gbnquad; j++) {

	  z = Gbs(j); // Z-coordinate of the corresponding quadrature point

	  // Loop over elements in fluid domain
	  for (int elem = 0; elem < nelem; elem++) {

	    // Compute barycentric coordinates of quad point w.r.t this element
	    // Reference: http://www.blackpawn.com/texts/pointinpoly/
	    v0(0) = sf->Nodes[ENC(elem, 3)]->x - sf->Nodes[ENC(elem, 1)]->x;
	    v0(1) = sf->Nodes[ENC(elem, 3)]->y - sf->Nodes[ENC(elem, 1)]->y;
	    v1(0) = sf->Nodes[ENC(elem, 2)]->x - sf->Nodes[ENC(elem, 1)]->x;
	    v1(1) = sf->Nodes[ENC(elem, 2)]->y - sf->Nodes[ENC(elem, 1)]->y;
	    v2(0) = x - sf->Nodes[ENC(elem, 1)]->x;
	    v2(1) = z - sf->Nodes[ENC(elem, 1)]->y;

	    si = (ublas::inner_prod(v1, v1) * ublas::inner_prod(v0, v2) -
		  ublas::inner_prod(v0, v1) * ublas::inner_prod(v1, v2)) /
	      (ublas::inner_prod(v0, v0) * ublas::inner_prod(v1, v1) -
	       ublas::inner_prod(v0, v1) * ublas::inner_prod(v0, v1));

	    ti = (ublas::inner_prod(v0, v0) * ublas::inner_prod(v1, v2) -
		  ublas::inner_prod(v0, v1) * ublas::inner_prod(v0, v2)) /
	      (ublas::inner_prod(v0, v0) * ublas::inner_prod(v1, v1) -
	       ublas::inner_prod(v0, v1) * ublas::inner_prod(v0, v1));

	    // If point lies inside element, store the mapping
	    // The 1e-12 is a hack for end cases where the Gauss point lies on
	    // the boundary of the domain
	    if ((si >= -1e-12) && (ti >= -1e-12) && (si + ti <= 1 + 1e-12)) {
	      IntegrationPoint2Element(node, j) = elem;
	      break;
	    }

	  } // End of loop over elements in fluid domain

	} // End of loop over Gauss points

      } // End of loop over the nodes in the edge
    } // End if

  } // End of loop over edges

}
//------------------------------------------------------------------------------
void Fluid::PreProcess() {

  // Store current pressure as the pressure at the prev timestep, Pold
  copy((sf->P).begin(), (sf->P).end(), (sf->Pold).begin());

  // Store current gap height as the gap height at the prev timestep, Hold
  copy((sf->H).begin(), (sf->H).end(), (sf->Hold).begin());

}
//------------------------------------------------------------------------------
void Fluid::SolveDynamic(double tn, double dtn) {

  int iter = 0;
  double err = 1.0, eps = 1e-10;
  t = tn; dt = dtn;

  // Compute gap height and X-displacement of projected domain
  CompGapHeight();

  // Solve for the pressure at the end of new timestep
  ApplyInhomogeneousDBC();
  do {
    K.resize(ndof,ndof);
    RHS = unyque::DVector_zero(ndof);
    dU = unyque::DVector_zero(ndof);
    CompDomIntegrals();
    dU = unyque::umfpackSolve(K, RHS);
    err = ublas::norm_inf(dU);
    UpdateGlobalPressure();
    iter++;
    if (c->DEBUG) cout<<"iteration: "<<iter<<"  Error: "<<err<<endl;
  } while (err > eps && iter <= 25);
  if (iter > 25)
    cout << "Warning: Too many Newton-Raphson iterations" << endl;
  if (c->DEBUG) PrintResults();
  if (c->DEBUG)
    cout<<"No. of steps: "<<iter<<"   Max. change in pressure: "<<
      ublas::norm_inf((sf->P) - (sf->Pold))<<endl;

  // Integrate fluid pressure to physical domain as traction
  CompPressureTraction();

}
//------------------------------------------------------------------------------
void Fluid::CompGapHeight() {

  fem::FEM_Edge *ed;
  unyque::DVector N, Quantity;
  double x, x1, x2, x3;

  N = unyque::DVector_zero(3);
  Quantity = unyque::DVector_zero(3);

  for (int i = 0; i < nnode; i++) {

    // Get edge on moving boundary
    ed = s->BEdges[Node2MovingEdge(i)+1];

    // Store x-coordinates of current node and moving edge nodes
    x = sf->Nodes[i+1]->x;
    x1 = s->Nodes[ed->node1]->x;
    x2 = s->Nodes[ed->node2]->x;
    x3 = s->Nodes[ed->node3]->x;

    // Compute basis functions
    N(0) = (x - x2)*(x - x3)/(x1 - x2)/(x1 - x3);
    N(1) = (x - x1)*(x - x3)/(x2 - x1)/(x2 - x3);
    N(2) = (x - x1)*(x - x2)/(x3 - x1)/(x3 - x2);

    // Get y-position of edge nodes
    Quantity(0) = s->Nodes[ed->node1]->y + (s->V)(ed->node1-1);
    Quantity(1) = s->Nodes[ed->node2]->y + (s->V)(ed->node2-1);
    Quantity(2) = s->Nodes[ed->node3]->y + (s->V)(ed->node3-1);

    // Store y-position of edge at x-location corresponding to node
    (sf->H)(i) = ublas::inner_prod(N, Quantity);

    // Get x-displacement of edge nodes
    Quantity(0) = (s->U)(ed->node1-1);
    Quantity(1) = (s->U)(ed->node2-1);
    Quantity(2) = (s->U)(ed->node3-1);

    // Store x-diplacement of edge at x-location corresponding to node
    (sf->U)(i) = ublas::inner_prod(N, Quantity);

    // Get edge on fixed boundary
    ed = s->BEdges[Node2FixedEdge(i)+1];

    // Get y-position of edge nodes
    Quantity(0) = s->Nodes[ed->node1]->y + (s->V)(ed->node1-1);
    Quantity(1) = s->Nodes[ed->node2]->y + (s->V)(ed->node2-1);
    Quantity(2) = s->Nodes[ed->node3]->y + (s->V)(ed->node3-1);

    // Compute gap height at this node
    (sf->H)(i) = abs((sf->H)(i) - ublas::inner_prod(N, Quantity));

  } // End loop over nodes

}
//------------------------------------------------------------------------------
void Fluid::ApplyInhomogeneousDBC() {

  int bcno;
  fem::FEM_Edge *ed;

  // Loop over boundary edges
  for (int eid = 0; eid < nbedge; eid++) {

    ed = sf->BEdges[eid+1];
    bcno = -1;
    // Find the b.c. no: corresponding to the edge's marker (Default -1)
    for (int i = 0; i < nbc; i++) {
      if (BCtype(i,0) == ed->bmarker) {bcno = i; break;}
    }

    // If it is a Dirichlet b.c. then initialize the corresponding nodes in the
    // solution vector with the value of the prescribed pressure
    if (BCtype(bcno,1) == 1) {
      (sf->P)(ed->node1 - 1) = BCvals(bcno,0);
      (sf->P)(ed->node2 - 1) = BCvals(bcno,0);
      (sf->P)(ed->node3 - 1) = BCvals(bcno,0);
    }

  } // End of loop over boundary edges

}
//------------------------------------------------------------------------------
void Fluid::CompDomIntegrals() {

  int Gnode, iglobal, jglobal;
  double si, ti, P, Pold, H, Hold;
  unyque::DMatrix Ecoor, dN, B, Finv, FB, BFFB, Ke;
  unyque::DVector Pe, Pe_old, He, He_old, RHSe, N;

  Ecoor = unyque::DMatrix_zero(enode,2);
  Pe = unyque::DVector_zero(enode);
  Pe_old = unyque::DVector_zero(enode);
  He = unyque::DVector_zero(enode);
  He_old = unyque::DVector_zero(enode);
  N = unyque::DVector_zero(enode);
  dN = unyque::DMatrix_zero(2, enode);

  // Loop over elements
  for (int eid = 0; eid < nelem; eid++) {

    // Loop over nodes in that element
    for (int i = 0; i < enode; i++) {
      Gnode = ENC(eid,i+1); // Global node corresponding to (i+1)th local node
      Ecoor(i,0) = sf->Nodes[Gnode]->x; // X-coordinate at that node
      Ecoor(i,1) = sf->Nodes[Gnode]->y; // Z-coordinate at that node
      Pe(i) = (sf->P)(Gnode - 1); // Pressure at that node at timestep n
      Pe_old(i) = (sf->Pold)(Gnode - 1); // Pressure at timestep (n-1)
      He(i) = (sf->H)(Gnode - 1); // Gap height at that node
      He_old(i) = (sf->Hold)(Gnode - 1); // Gap height at timestep (n-1)
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

      // Compute the pressure and gap height at the current integration point
      P = ublas::inner_prod(N, Pe);
      Pold = ublas::inner_prod(N, Pe_old);
      H = ublas::inner_prod(N, He);
      Hold = ublas::inner_prod(N, He_old);

      // Compute the Knudsen number
      KN = LAMBDA/H/P_ATM;

      // Compute the element matrix Ke
      // {12*ETA*(H/dt + Hd)*[NN] + (1+6*KN)*H^3*[BFFB]*([Pe][N]+P*Inxn)}*detF
      // Ke += (ublas::outer_prod(12*ETA*N*(H/dt + Hd)/101325 +
      // 			       (1+6*KN)*pow(H,3)*ublas::prod(BFFB, Pe), N) +
      // 	     (1+6*KN)*pow(H,3)*P*BFFB)*detF*Gw(ip)*detJ;
      Ke += (ublas::outer_prod(12*ETA*N*(H/dt)/101325 +
      			       (1+6*KN)*pow(H,3)*ublas::prod(BFFB, Pe), N) +
      	     (1+6*KN)*pow(H,3)*(P+P_ATM)*BFFB)*detF*Gw(ip)*detJ;

      // Create the elemental vector RHSe
      // {12*ETA*(H*(P-Pold)/dt + Hd*P)*[N] + (1+6*KN)*H^3*P*[BFFB]*[Pe]}*detF
      // RHSe -= (12*ETA*N*(H*(P - Pold)/dt + Hd*P)/101325 +
      // 	       (1+6*KN)*pow(H,3)*(P*ublas::prod(BFFB, Pe))
      // 	       )*detF*Gw(ip)*detJ;
      RHSe -= (12*ETA*N*(H*(P+P_ATM) - Hold*(Pold+P_ATM))/dt/101325 +
	       (1+6*KN)*pow(H,3)*((P+P_ATM)*ublas::prod(BFFB, Pe))
	       )*detF*Gw(ip)*detJ;

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
void Fluid::CompN(double s, double t, unyque::DVector &N, unyque::DMatrix &dN) {
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
unyque::DMatrix Fluid::CompJandB(unyque::DMatrix &dN, unyque::DMatrix &Ecoor) {

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
unyque::DMatrix Fluid::CompF(int eid, unyque::DMatrix &B) {

  int Gnode;
  double detFinv;
  unyque::DMatrix Ue, F, Finv;

  Ue = unyque::DMatrix_zero(enode, 2);
  F = unyque::DMatrix_zero(2,2);
  Finv = unyque::DMatrix(2,2);

  // Get displacement at nodes. Domain only deforms along x, not z.
  for (int i = 0; i < enode; i++) {
    Gnode = ENC(eid,i+1); // Global node corresponding to (i+1)th local node
    Ue(i, 0) = (sf->U)(Gnode-1);
    Ue(i, 1) = 0;
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
void Fluid::UpdateGlobalPressure() {
  for (int i = 0; i < nnode; i++)
    if (L2G(i) > 0)
      (sf->P)(i) += dU(L2G(i) - 1);
}
//------------------------------------------------------------------------------
void Fluid::PrintResults() {
  FILE *fp;
  fem::FEM_Point *pp;
  double xval, yval;

  fp = fopen("postproc/pressure.dat","w");
  for (int i = 0; i<nnode; i++) {
    pp = sf->Nodes[i+1];
    xval = pp->x + (sf->U)(i);
    yval = pp->y;
    fprintf(fp,"%.14lf  %.14lf %.14lf \n", xval, yval, (sf->P)(i));
  }
  fclose(fp);
}
//------------------------------------------------------------------------------
void Fluid::CompPressureTraction() {

  int elem;
  double x, z, si, ti;
  unyque::DVector v0, v1, v2, N, Pe;

  v0 = unyque::DVector(2);
  v1 = unyque::DVector(2);
  v2 = unyque::DVector(2);
  N  = unyque::DVector(enode);
  Pe = unyque::DVector(enode);

  fill((s->Pf).begin(), (s->Pf).end(), 0.0);

  // Loop over nodes in physical domain
  for (int node = 0; node < s->nnode; node++) {

    // Confirm that this node lies on the moving edge by checking if it has a
    // valid element mapping for the first integration point. If it doesn't,
    // then don't do anything for this node
    if (IntegrationPoint2Element(node, 0) < 0)
      continue;

    // X-coordinate of this node
    x = s->Nodes[node + 1]->x;

    for (int j = 0; j < Gbnquad; j++) {

      z = Gbs(j); // Z-coordinate of the corresponding quadrature point

      elem = IntegrationPoint2Element(node, j);

      // Compute barycentric coordinates of quad point w.r.t this element
      // Reference: http://www.blackpawn.com/texts/pointinpoly/
      v0(0) = sf->Nodes[ENC(elem, 3)]->x - sf->Nodes[ENC(elem, 1)]->x;
      v0(1) = sf->Nodes[ENC(elem, 3)]->y - sf->Nodes[ENC(elem, 1)]->y;
      v1(0) = sf->Nodes[ENC(elem, 2)]->x - sf->Nodes[ENC(elem, 1)]->x;
      v1(1) = sf->Nodes[ENC(elem, 2)]->y - sf->Nodes[ENC(elem, 1)]->y;
      v2(0) = x - sf->Nodes[ENC(elem, 1)]->x;
      v2(1) = z - sf->Nodes[ENC(elem, 1)]->y;

      si = (ublas::inner_prod(v1, v1) * ublas::inner_prod(v0, v2) -
	    ublas::inner_prod(v0, v1) * ublas::inner_prod(v1, v2)) /
	(ublas::inner_prod(v0, v0) * ublas::inner_prod(v1, v1) -
	 ublas::inner_prod(v0, v1) * ublas::inner_prod(v0, v1));

      ti = (ublas::inner_prod(v0, v0) * ublas::inner_prod(v1, v2) -
	    ublas::inner_prod(v0, v1) * ublas::inner_prod(v0, v2)) /
	(ublas::inner_prod(v0, v0) * ublas::inner_prod(v1, v1) -
	 ublas::inner_prod(v0, v1) * ublas::inner_prod(v0, v1));

      // Compute basis functions at the quadrature point
      N(3) = 4.0*si*(1.0-si-ti);
      N(4) = 4.0*si*ti;
      N(5) = 4.0*ti*(1.0-si-ti);
      N(0) = (1.0-si-ti)-N(3)/2.0-N(5)/2.0;
      N(1) = si-N(4)/2.0-N(3)/2.0;
      N(2) = ti-N(5)/2.0-N(4)/2.0;

      // Store nodal pressure values for this element
      for (int k = 0; k < enode; k++)
	Pe(k) = (sf->P)(ENC(elem, k+1) - 1);

      // Add contribution to the traction integral after converting to Pa
      (s->Pf)(node) += 101325*ublas::inner_prod(N, Pe)*Gbw(j);

    } // End of loop over Gauss points

  } // End of loop over boundary nodes

}
//------------------------------------------------------------------------------
pyublas::numpy_vector<double> Fluid::Pressure() {
  pyublas::numpy_vector<double> rval((sf->P).size());
  rval.assign(sf->P);
  return rval;
}
//------------------------------------------------------------------------------
pyublas::numpy_vector<double> Fluid::GapHeight() {
  pyublas::numpy_vector<double> rval((sf->H).size());
  rval.assign(sf->H);
  return rval;
}
//------------------------------------------------------------------------------
