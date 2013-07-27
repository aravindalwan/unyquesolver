#include "fem.hpp"

//------------------------------------------------------------------------------
fem::FEM_Element::FEM_Element(int eid, int region,
			      int n1, int n2, int n3, int n4, int n5, int n6) {
  id = eid; reg = region;
  node1 = n1; node2 = n2; node3 = n3; node4 = n4; node5 = n5; node6 = n6;
}
//------------------------------------------------------------------------------
bp::tuple fem::FEM_Element_pickle_suite::getstate(bp::object o) {
  FEM_Element const& e = bp::extract<FEM_Element const&>(o)();
  return bp::make_tuple(o.attr("__dict__"), e.id, e.reg, e.node1, e.node2,
			e.node3, e.node4, e.node5, e.node6);
}
//------------------------------------------------------------------------------
void fem::FEM_Element_pickle_suite::setstate(bp::object o, bp::tuple state) {

  FEM_Element& e = bp::extract<FEM_Element&>(o)();

  if (len(state) != 9) {
    PyErr_SetObject(PyExc_ValueError,
		    ("expected 9-item tuple in call to __setstate__; got %s"
		     % state).ptr());
    bp::throw_error_already_set();
  }

  // restore the object's __dict__
  bp::dict d = bp::extract<bp::dict>(o.attr("__dict__"))();
  d.update(state[0]);

  // restore the internal state of the FEM_Point object
  e.id = bp::extract<int>(state[1]);
  e.reg = bp::extract<int>(state[2]);
  e.node1 = bp::extract<int>(state[3]);
  e.node2 = bp::extract<int>(state[4]);
  e.node3 = bp::extract<int>(state[5]);
  e.node4 = bp::extract<int>(state[6]);
  e.node5 = bp::extract<int>(state[7]);
  e.node6 = bp::extract<int>(state[8]);

}
//------------------------------------------------------------------------------
fem::FEM_Edge::FEM_Edge(int edid, int elementnum, int edgetype, int bm,
			double norm, int n1, int n2, int n3) {
  id = edid;
  eno = elementnum; eid = edgetype; bmarker = bm; normal = norm;
  node1 = n1; node2 = n2; node3 = n3;
}
//------------------------------------------------------------------------------
bp::tuple fem::FEM_Edge_pickle_suite::getstate(bp::object o) {
  FEM_Edge const& e = bp::extract<FEM_Edge const&>(o)();
  return bp::make_tuple(o.attr("__dict__"), e.id, e.eno, e.eid, e.bmarker,
			e.normal, e.node1, e.node2, e.node3);
}
//------------------------------------------------------------------------------
void fem::FEM_Edge_pickle_suite::setstate(bp::object o, bp::tuple state) {

  FEM_Edge& e = bp::extract<FEM_Edge&>(o)();

  if (len(state) != 9) {
    PyErr_SetObject(PyExc_ValueError,
		    ("expected 9-item tuple in call to __setstate__; got %s"
		     % state).ptr());
    bp::throw_error_already_set();
  }

  // restore the object's __dict__
  bp::dict d = bp::extract<bp::dict>(o.attr("__dict__"))();
  d.update(state[0]);

  // restore the internal state of the FEM_Point object
  e.id = bp::extract<int>(state[1]);
  e.eno = bp::extract<int>(state[2]);
  e.eid = bp::extract<int>(state[3]);
  e.bmarker = bp::extract<int>(state[4]);
  e.normal = bp::extract<double>(state[5]);
  e.node1 = bp::extract<int>(state[6]);
  e.node2 = bp::extract<int>(state[7]);
  e.node3 = bp::extract<int>(state[8]);

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
bp::tuple fem::FEM_Point_pickle_suite::getstate(bp::object o) {
  FEM_Point const& p = bp::extract<FEM_Point const&>(o)();
  return bp::make_tuple(o.attr("__dict__"), p.id, p.bmarker, p.x, p.y);
}
//------------------------------------------------------------------------------
void fem::FEM_Point_pickle_suite::setstate(bp::object o, bp::tuple state) {

  FEM_Point& p = bp::extract<FEM_Point&>(o)();

  if (len(state) != 5) {
    PyErr_SetObject(PyExc_ValueError,
		    ("expected 5-item tuple in call to __setstate__; got %s"
		     % state).ptr());
    bp::throw_error_already_set();
  }

  // restore the object's __dict__
  bp::dict d = bp::extract<bp::dict>(o.attr("__dict__"))();
  d.update(state[0]);

  // restore the internal state of the FEM_Point object
  p.id = bp::extract<int>(state[1]);
  p.bmarker = bp::extract<int>(state[2]);
  p.x = bp::extract<double>(state[3]);
  p.y = bp::extract<double>(state[4]);

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
  Pf = unyque::DVector_zero(nnode);

  T = unyque::DVector_scalar(nnode, 300.0);
  Told = unyque::DVector_scalar(nnode, 300.0);

  Phi = unyque::DVector_zero(nnode);

  BPhi = unyque::DVector_zero(nbedge);
  BdPhidn = unyque::DVector_zero(nbedge);
  SCharge = unyque::DVector_zero(nbedge);
  Ent = unyque::DMatrix_zero(nbedge, 2);

}
//------------------------------------------------------------------------------
bp::tuple fem::FEM_PhysicalDomain_pickle_suite::getstate(bp::object o) {

  bp::list nodes, edges, elements;
  bp::tuple u, v;

  FEM_PhysicalDomain const& p = bp::extract<FEM_PhysicalDomain const&>(o)();
  for (int i = 1; i <= p.nnode; i++)
    nodes.append(p.Nodes[i]);
  for (int i = 1; i <= p.nedge; i++)
    edges.append(p.Edges[i]);
  for (int i = 1; i <= p.nelem; i++)
    elements.append(p.Elements[i]);

  u = bp::make_tuple(p.U, p.Uold, p.Ud, p.Udold, p.Udd, p.Uddold);
  v = bp::make_tuple(p.V, p.Vold, p.Vd, p.Vdold, p.Vdd, p.Vddold);

  return bp::make_tuple(o.attr("__dict__"), nodes, edges, elements, p.id,
			u, v, p.T, p.Told, p.Phi, p.BPhi, p.BdPhidn, p.SCharge,
			p.Ent);
}
//------------------------------------------------------------------------------
void fem::FEM_PhysicalDomain_pickle_suite::setstate(bp::object o, bp::tuple state) {

  bp::tuple u, v;

  FEM_PhysicalDomain& p = bp::extract<FEM_PhysicalDomain&>(o)();

  if (len(state) != 14) {
    PyErr_SetObject(PyExc_ValueError,
		    ("expected 14-item tuple in call to __setstate__; got %s"
		     % state).ptr());
    bp::throw_error_already_set();
  }

  // restore the object's __dict__
  bp::dict d = bp::extract<bp::dict>(o.attr("__dict__"))();
  d.update(state[0]);

  // restore the internal state of the FEM_FluidDomain object
  p.SetMesh(bp::extract<bp::list>(state[1]), bp::extract<bp::list>(state[2]),
  	    bp::extract<bp::list>(state[3]));
  p.id = bp::extract<int>(state[4]);

  u = bp::extract< bp::tuple >(state[5]);
  p.U = bp::extract< pyublas::numpy_vector<double> >(u[0]);
  p.Uold = bp::extract< pyublas::numpy_vector<double> >(u[1]);
  p.Ud = bp::extract< pyublas::numpy_vector<double> >(u[2]);
  p.Udold = bp::extract< pyublas::numpy_vector<double> >(u[3]);
  p.Udd = bp::extract< pyublas::numpy_vector<double> >(u[4]);
  p.Uddold = bp::extract< pyublas::numpy_vector<double> >(u[5]);

  v = bp::extract< bp::tuple >(state[6]);
  p.V = bp::extract< pyublas::numpy_vector<double> >(v[0]);
  p.Vold = bp::extract< pyublas::numpy_vector<double> >(v[1]);
  p.Vd = bp::extract< pyublas::numpy_vector<double> >(v[2]);
  p.Vdold = bp::extract< pyublas::numpy_vector<double> >(v[3]);
  p.Vdd = bp::extract< pyublas::numpy_vector<double> >(v[4]);
  p.Vddold = bp::extract< pyublas::numpy_vector<double> >(v[5]);

  p.T = bp::extract< pyublas::numpy_vector<double> >(state[7]);
  p.Told = bp::extract< pyublas::numpy_vector<double> >(state[8]);
  p.Phi = bp::extract< pyublas::numpy_vector<double> >(state[9]);
  p.BPhi = bp::extract< pyublas::numpy_vector<double> >(state[10]);
  p.BdPhidn = bp::extract< pyublas::numpy_vector<double> >(state[11]);
  p.SCharge = bp::extract< pyublas::numpy_vector<double> >(state[12]);
  p.Ent = bp::extract< pyublas::numpy_matrix<double> >(state[13]);

}
//------------------------------------------------------------------------------
fem::FEM_FluidDomain::FEM_FluidDomain() : fem::FEM_Domain() {
}
//------------------------------------------------------------------------------
void fem::FEM_FluidDomain::InitDOFs() {

  P = unyque::DVector_zero(nnode);
  Pold = unyque::DVector_zero(nnode);

  U = unyque::DVector_zero(nnode);

  H = unyque::DVector_scalar(nnode, 2.0);
  Hold = unyque::DVector_scalar(nnode, 2.0);

}
//------------------------------------------------------------------------------
bp::tuple fem::FEM_FluidDomain_pickle_suite::getstate(bp::object o) {

  bp::list nodes, edges, elements;

  FEM_FluidDomain const& p = bp::extract<FEM_FluidDomain const&>(o)();
  for (int i = 1; i <= p.nnode; i++)
    nodes.append(p.Nodes[i]);
  for (int i = 1; i <= p.nedge; i++)
    edges.append(p.Edges[i]);
  for (int i = 1; i <= p.nelem; i++)
    elements.append(p.Elements[i]);

  return bp::make_tuple(o.attr("__dict__"), nodes, edges, elements, p.id,
			p.P, p.Pold, p.U, p.H, p.Hold);
}
//------------------------------------------------------------------------------
void fem::FEM_FluidDomain_pickle_suite::setstate(bp::object o, bp::tuple state) {

  FEM_FluidDomain& p = bp::extract<FEM_FluidDomain&>(o)();

  if (len(state) != 10) {
    PyErr_SetObject(PyExc_ValueError,
		    ("expected 10-item tuple in call to __setstate__; got %s"
		     % state).ptr());
    bp::throw_error_already_set();
  }

  // restore the object's __dict__
  bp::dict d = bp::extract<bp::dict>(o.attr("__dict__"))();
  d.update(state[0]);

  // restore the internal state of the FEM_FluidDomain object
  p.SetMesh(bp::extract<bp::list>(state[1]), bp::extract<bp::list>(state[2]),
  	    bp::extract<bp::list>(state[3]));
  p.id = bp::extract<int>(state[4]);
  p.P = bp::extract< pyublas::numpy_vector<double> >(state[5]);
  p.Pold = bp::extract< pyublas::numpy_vector<double> >(state[6]);
  p.U = bp::extract< pyublas::numpy_vector<double> >(state[7]);
  p.H = bp::extract< pyublas::numpy_vector<double> >(state[8]);
  p.Hold = bp::extract< pyublas::numpy_vector<double> >(state[9]);

}
//------------------------------------------------------------------------------
fem::FEM_Common::FEM_Common() {
  DEBUG = 0;
  Phi_mult = 1; Phi_inf = 0;
  original_gap = new_gap = 2.0;
  t = t_stop = dt = 0;
  functions = new Function();
}
//------------------------------------------------------------------------------
bp::tuple fem::FEM_Common_pickle_suite::getstate(bp::object o) {
  FEM_Common const& c = bp::extract<FEM_Common const&>(o)();
  return bp::make_tuple(o.attr("__dict__"), c.DEBUG, c.Phi_mult, c.Phi_inf,
			c.original_gap, c.new_gap, c.t, c.t_stop, c.dt);
}
//------------------------------------------------------------------------------
void fem::FEM_Common_pickle_suite::setstate(bp::object o, bp::tuple state) {

  FEM_Common& c = bp::extract<FEM_Common&>(o)();

  if (len(state) != 9) {
    PyErr_SetObject(PyExc_ValueError,
		    ("expected 9-item tuple in call to __setstate__; got %s"
		     % state).ptr());
    bp::throw_error_already_set();
  }

  // restore the object's __dict__
  bp::dict d = bp::extract<bp::dict>(o.attr("__dict__"))();
  d.update(state[0]);

  // restore the internal state of the FEM_Common object
  c.DEBUG = bp::extract<int>(state[1]);
  c.Phi_mult = bp::extract<double>(state[2]);
  c.Phi_inf = bp::extract<double>(state[3]);
  c.original_gap = bp::extract<double>(state[4]);
  c.new_gap = bp::extract<double>(state[5]);
  c.t = bp::extract<double>(state[6]);
  c.t_stop = bp::extract<double>(state[7]);
  c.dt = bp::extract<double>(state[8]);

}
//------------------------------------------------------------------------------
