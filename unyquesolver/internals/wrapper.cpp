#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <pyublas/numpy.hpp>

#include "solver.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(_internals)
{

  class_<fem::Solver>("Solver", init<int, list>())
    // Reference the PhysicalDomain object pointer as an attribute
    // stackoverflow.com/questions/2541446/exposing-a-pointer-in-boost-python
    .add_property("physical_domain",
		  make_getter(&fem::Solver::s,
			      return_value_policy<reference_existing_object>()))
    .add_property("fluid_domain",
		  make_getter(&fem::Solver::sf,
			      return_value_policy<reference_existing_object>()))
    .add_property("common_parameters",
		  make_getter(&fem::Solver::c,
			      return_value_policy<reference_existing_object>()))
    .def("InitPhysical", &fem::Solver::InitPhysical)
    .def("InitFluid", &fem::Solver::InitFluid)
    .def("Init", &fem::Solver::Init)
    .def("Restore", &fem::Solver::Restore)
    .def("Solve", &fem::Solver::Solve)
    ;

  class_<fem::FEM_Element, fem::FEM_Element*>("Element",
			   init<int, int, int, int, int, int, int, int>())
    .def_readwrite("id", &fem::FEM_Element::id)
    .def_readwrite("region", &fem::FEM_Element::reg)
    .def_readwrite("node1", &fem::FEM_Element::node1)
    .def_readwrite("node2", &fem::FEM_Element::node2)
    .def_readwrite("node3", &fem::FEM_Element::node3)
    .def_readwrite("node4", &fem::FEM_Element::node4)
    .def_readwrite("node5", &fem::FEM_Element::node5)
    .def_readwrite("node6", &fem::FEM_Element::node6)
    .def_pickle(fem::FEM_Element_pickle_suite())
    ;

  class_<std::vector<fem::FEM_Element*> >("ElementVec")
    .def(vector_indexing_suite<std::vector<fem::FEM_Element*> >())
    ;

  class_<fem::FEM_Edge, fem::FEM_Edge*>("Edge",
			init<int, int, int, int, double, int, int, int>())
    .def_readwrite("id", &fem::FEM_Edge::id)
    .def_readwrite("element", &fem::FEM_Edge::eno)
    .def_readwrite("edge_type", &fem::FEM_Edge::eid)
    .def_readwrite("boundary_marker", &fem::FEM_Edge::bmarker)
    .def_readwrite("normal", &fem::FEM_Edge::normal)
    .def_readwrite("node1", &fem::FEM_Edge::node1)
    .def_readwrite("node2", &fem::FEM_Edge::node2)
    .def_readwrite("node3", &fem::FEM_Edge::node3)
    .def_pickle(fem::FEM_Edge_pickle_suite())
    ;

  class_<std::vector<fem::FEM_Edge*> >("EdgeVec")
    .def(vector_indexing_suite<std::vector<fem::FEM_Edge*> >())
    ;

  class_<fem::FEM_Point, fem::FEM_Point*>("Point",
					  init<int, int, double, double>())
    .def_readwrite("id", &fem::FEM_Point::id)
    .def_readwrite("boundary_marker", &fem::FEM_Point::bmarker)
    .def_readwrite("x", &fem::FEM_Point::x)
    .def_readwrite("y", &fem::FEM_Point::y)
    .def_pickle(fem::FEM_Point_pickle_suite())
    ;

  class_<std::vector<fem::FEM_Point*> >("PointVec")
    .def(vector_indexing_suite<std::vector<fem::FEM_Point*> >())
    ;

  class_<fem::FEM_Domain, boost::noncopyable>("Domain", no_init)
    .def_readwrite("Nodes", &fem::FEM_Domain::Nodes)
    .def_readwrite("RefNodes", &fem::FEM_Domain::RefNodes)
    .def_readwrite("BNodes", &fem::FEM_Domain::BNodes)
    .def_readwrite("Elements", &fem::FEM_Domain::Elements)
    .def_readwrite("Edges", &fem::FEM_Domain::Edges)
    .def_readwrite("BEdges", &fem::FEM_Domain::BEdges)
    .def_readwrite("id", &fem::FEM_Domain::id)
    ;

  class_<fem::FEM_PhysicalDomain, bases<fem::FEM_Domain> >("PhysicalDomain")
    .def(pyublas::by_value_ro_member("U", &fem::FEM_PhysicalDomain::U))
    .def(pyublas::by_value_ro_member("Uold", &fem::FEM_PhysicalDomain::Uold))
    .def(pyublas::by_value_ro_member("Ud", &fem::FEM_PhysicalDomain::Ud))
    .def(pyublas::by_value_ro_member("Udold", &fem::FEM_PhysicalDomain::Udold))
    .def(pyublas::by_value_ro_member("Udd", &fem::FEM_PhysicalDomain::Udd))
    .def(pyublas::by_value_ro_member("Uddold", &fem::FEM_PhysicalDomain::Uddold))
    .def(pyublas::by_value_ro_member("V", &fem::FEM_PhysicalDomain::V))
    .def(pyublas::by_value_ro_member("Vold", &fem::FEM_PhysicalDomain::Vold))
    .def(pyublas::by_value_ro_member("Vd", &fem::FEM_PhysicalDomain::Vd))
    .def(pyublas::by_value_ro_member("Vdold", &fem::FEM_PhysicalDomain::Vdold))
    .def(pyublas::by_value_ro_member("Vdd", &fem::FEM_PhysicalDomain::Vdd))
    .def(pyublas::by_value_ro_member("Vddold", &fem::FEM_PhysicalDomain::Vddold))
    .def(pyublas::by_value_ro_member("T", &fem::FEM_PhysicalDomain::T))
    .def(pyublas::by_value_ro_member("Told", &fem::FEM_PhysicalDomain::Told))
    .def(pyublas::by_value_ro_member("Phi", &fem::FEM_PhysicalDomain::Phi))
    .def(pyublas::by_value_ro_member("BPhi", &fem::FEM_PhysicalDomain::BPhi))
    .def(pyublas::by_value_ro_member("BdPhidn", &fem::FEM_PhysicalDomain::BdPhidn))
    .def(pyublas::by_value_ro_member("SCharge", &fem::FEM_PhysicalDomain::SCharge))
    .def(pyublas::by_value_ro_member("Ent", &fem::FEM_PhysicalDomain::Ent))
    .def_pickle(fem::FEM_PhysicalDomain_pickle_suite())
    ;

  class_<fem::FEM_FluidDomain, bases<fem::FEM_Domain> >("FluidDomain")
    .def(pyublas::by_value_ro_member("P", &fem::FEM_FluidDomain::P))
    .def(pyublas::by_value_ro_member("Pold", &fem::FEM_FluidDomain::Pold))
    .def(pyublas::by_value_ro_member("U", &fem::FEM_FluidDomain::U))
    .def(pyublas::by_value_ro_member("H", &fem::FEM_FluidDomain::H))
    .def(pyublas::by_value_ro_member("Hold", &fem::FEM_FluidDomain::Hold))
    .def_pickle(fem::FEM_FluidDomain_pickle_suite())
    ;

  class_<fem::FEM_Common>("CommonParameters")
    .def_readwrite("Phi_mult", &fem::FEM_Common::Phi_mult)
    .def_readwrite("Phi_inf", &fem::FEM_Common::Phi_inf)
    .def_readwrite("original_gap", &fem::FEM_Common::original_gap)
    .def_readwrite("new_gap", &fem::FEM_Common::new_gap)
    .def_readwrite("t", &fem::FEM_Common::t)
    .def_readwrite("t_stop", &fem::FEM_Common::t_stop)
    .def_readwrite("dt", &fem::FEM_Common::dt)
    .def_pickle(fem::FEM_Common_pickle_suite())
    ;
}
