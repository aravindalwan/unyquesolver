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
fem::FEM_Point2D::FEM_Point2D(double ix, double iy) {
  id = 0; x = ix; y = iy; bmarker = 0;
}
//------------------------------------------------------------------------------
fem::FEM_Point2D::FEM_Point2D(int pid, int bm, double ix, double iy) {
  id = pid; bmarker = bm; x = ix; y = iy;
}
//------------------------------------------------------------------------------
fem::FEM_Surface2D::FEM_Surface2D() {
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
void fem::FEM_Surface2D::SetMesh(bp::list nodes, bp::list edges,
				 bp::list elements) {

  SetNodes(nodes);
  SetEdges(edges);
  SetElements(elements);
  InitDOFs();

}
//------------------------------------------------------------------------------
void fem::FEM_Surface2D::SetNodes(bp::list nodes) {

  int bpid = 0;
  fem::FEM_Point2D *pp, *ppref;

  nnode = bp::len(nodes);

  for (int i = 0; i < nnode; i++) {
    pp = bp::extract<fem::FEM_Point2D*>(nodes[i])();
    Nodes.push_back(pp);
    ppref = new fem::FEM_Point2D(pp->x, pp->y);
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
void fem::FEM_Surface2D::SetEdges(bp::list edges) {

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
void fem::FEM_Surface2D::SetElements(bp::list elements) {

  fem::FEM_Element *ee;

  nelem = bp::len(elements);

  for (int i = 0; i < nelem; i++) {
    ee = bp::extract<fem::FEM_Element*>(elements[i])();
    Elements.push_back(ee);
  }

}
//------------------------------------------------------------------------------
void fem::FEM_Surface2D::InitDOFs() {

  U = unyque::DVector_zero(nnode); V = unyque::DVector_zero(nnode);

  T = unyque::DVector_scalar(nnode, 300.0);
  Told = unyque::DVector_scalar(nnode, 300.0);

  Phi = unyque::DVector_zero(nnode);

  BPhi = unyque::DVector_zero(nbedge);
  BdPhidn = unyque::DVector_zero(nbedge);
  SCharge = unyque::DVector_zero(nbedge);
  Ent = unyque::DMatrix_zero(nbedge, 2);

}
//------------------------------------------------------------------------------
void fem::FEM_Surface2D::MoveMesh(unyque::DMatrix &displacement) {

  fem::FEM_Point2D *pp, *ppref;
  for (int i = 0; i < nnode; i++) {
    pp = Nodes[i+1];
    ppref = RefNodes[i+1];
    pp->x = ppref->x + displacement(i,0);
    pp->y = ppref->y + displacement(i,1);
  }

}
//------------------------------------------------------------------------------
FEM_Common::FEM_Common() {
  DEBUG = 0;
  Phi_mult = 1; Phi_inf = 0;
  functions = new Function();
}
//------------------------------------------------------------------------------
