#include "fem.hpp"

//------------------------------------------------------------------------------
fem::FEM_Element::FEM_Element(int eid, int region,
			      int n1, int n2, int n3, int n4, int n5, int n6) {
  id = eid; reg = region;
  node1 = n1; node2 = n2; node3 = n3; node4 = n4; node5 = n5; node6 = n6;
}
//------------------------------------------------------------------------------
fem::FEM_Edge::FEM_Edge(int edid, int elementnum, int edgetype, int bm,
			double norm, int n1, int n2, int n3) {
  id = edid;
  eno = elementnum; eid = edgetype; bmarker = bm; normal = norm;
  node1 = n1; node2 = n2; node3 = n3;
}
//------------------------------------------------------------------------------
fem::FEM_Point::FEM_Point(double ix, double iy) {
  id = 0; x = ix; y = iy; bmarker = 0;
}
//------------------------------------------------------------------------------
fem::FEM_Point::FEM_Point(int pid, int bm, double ix, double iy) {
  id = pid; bmarker = bm; x = ix; y = iy;
}
//------------------------------------------------------------------------------
fem::FEM_Domain::FEM_Domain() {
  Nodes.push_back(NULL);
  RefNodes.push_back(NULL);
  Elements.push_back(NULL);
  Edges.push_back(NULL);
  BEdges.push_back(NULL);
  BNodes.push_back(NULL);
  nnode = nbnode = nelem = nedge = nbedge = 0;
  id = -1;
}
//------------------------------------------------------------------------------
void fem::FEM_Domain::SetMesh(bp::list nodes, bp::list edges,
				 bp::list elements) {

  SetNodes(nodes);
  SetEdges(edges);
  SetElements(elements);
  InitDOFs();

}
//------------------------------------------------------------------------------
void fem::FEM_Domain::SetNodes(bp::list nodes) {

  int bpid = 0;
  fem::FEM_Point *pp, *ppref;

  nnode = bp::len(nodes);

  for (int i = 0; i < nnode; i++) {
    pp = bp::extract<fem::FEM_Point*>(nodes[i])();
    Nodes.push_back(pp);
    ppref = new fem::FEM_Point(pp->x, pp->y);
    ppref->id = pp->id;
    ppref->bmarker = pp->bmarker;
    RefNodes.push_back(ppref);
    if (pp->bmarker == 0) {
      bpid++;
      BNodes.push_back(pp);
    }
  }
  nbnode = bpid;

}
//------------------------------------------------------------------------------
void fem::FEM_Domain::SetEdges(bp::list edges) {

  int bpid = 0;
  fem::FEM_Edge *ed;

  nedge = bp::len(edges);

  for (int i = 0; i < nedge; i++) {
    ed = bp::extract<fem::FEM_Edge*>(edges[i])();
    Edges.push_back(ed);
    if (ed->bmarker != 0) {
      bpid++;
      ed->id = bpid;
      BEdges.push_back(ed);
    }
  }

  nbedge = bpid;

}
//------------------------------------------------------------------------------
void fem::FEM_Domain::SetElements(bp::list elements) {

  fem::FEM_Element *ee;

  nelem = bp::len(elements);

  for (int i = 0; i < nelem; i++) {
    ee = bp::extract<fem::FEM_Element*>(elements[i])();
    Elements.push_back(ee);
  }

}
//------------------------------------------------------------------------------
void fem::FEM_Domain::MoveMesh(unyque::DMatrix &displacement) {

  fem::FEM_Point *pp, *ppref;
  for (int i = 0; i < nnode; i++) {
    pp = Nodes[i+1];
    ppref = RefNodes[i+1];
    pp->x = ppref->x + displacement(i,0);
    pp->y = ppref->y + displacement(i,1);
  }

}
//------------------------------------------------------------------------------
fem::FEM_PhysicalDomain::FEM_PhysicalDomain() : fem::FEM_Domain() {
}
//------------------------------------------------------------------------------
void fem::FEM_PhysicalDomain::InitDOFs() {

  U = unyque::DVector_zero(nnode); V = unyque::DVector_zero(nnode);
  Uold = unyque::DVector_zero(nnode); Vold = unyque::DVector_zero(nnode);
  Ud = unyque::DVector_zero(nnode); Vd = unyque::DVector_zero(nnode);
  Udold = unyque::DVector_zero(nnode); Vdold = unyque::DVector_zero(nnode);
  Udd = unyque::DVector_zero(nnode); Vdd = unyque::DVector_zero(nnode);
  Uddold = unyque::DVector_zero(nnode); Vddold = unyque::DVector_zero(nnode);

  T = unyque::DVector_scalar(nnode, 300.0);
  Told = unyque::DVector_scalar(nnode, 300.0);

  Phi = unyque::DVector_zero(nnode);

  BPhi = unyque::DVector_zero(nbedge);
  BdPhidn = unyque::DVector_zero(nbedge);
  SCharge = unyque::DVector_zero(nbedge);
  Ent = unyque::DMatrix_zero(nbedge, 2);

}
//------------------------------------------------------------------------------
fem::FEM_FluidDomain::FEM_FluidDomain() : fem::FEM_Domain() {
}
//------------------------------------------------------------------------------
void fem::FEM_FluidDomain::InitDOFs() {

  P = unyque::DVector_zero(nnode);
  Pold = unyque::DVector_zero(nnode);

  H = unyque::DVector_scalar(nnode, 2.0);

}
//------------------------------------------------------------------------------
FEM_Common::FEM_Common() {
  DEBUG = 0;
  Phi_mult = 1; Phi_inf = 0;
  functions = new Function();
}
//------------------------------------------------------------------------------
