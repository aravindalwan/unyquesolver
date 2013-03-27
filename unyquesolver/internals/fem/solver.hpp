#ifndef __SOLVER_H
#define __SOLVER_H

#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include "fem.hpp"
#include "nonelast.hpp"
#include "therm.hpp"
#include "elec.hpp"
#include "eles.hpp"
namespace bp = boost::python;

// Define parameters
#define YOUNGS_MODULUS        1
#define THERMAL_CONDUCTIVITY  2
#define INTER_ELECTRODE_GAP   3
#define VOLTAGE               4

// Define analysis modes
#define MECHANICS                 1
#define THERMAL                   2
#define ELECTRICAL                3
#define ELECTROTHERMAL            4
#define ELECTROTHERMAL_ACTUATION  5
#define HYBRID_ELECTROSTATICS     6
#define HYBRID_ETM_ACTUATION      7
#define HYBRID_ETM_PULLIN         8

namespace fem {

  class Solver {

  public:

    Solver(int analysisType, bp::list pType);
    ~Solver();

    void Init(bp::list nodes, bp::list edges, bp::list elements);
    bp::object Solve(bp::list params);

    FEM_Common *c;
    FEM_PhysicalDomain *s;

  private:

    bool useNonElast, useTherm, useElec, useElEs;
    int nParams, *paramType;

    NonElast *nelast;
    Therm *therm;
    Elec *elec;
    ElEs *eles;
    double newGap, originalGap;

    bp::object (Solver::*analysis)();

    bp::object NLMechanics();
    bp::object Thermal();
    bp::object Electrical();
    bp::object Electrothermal();
    bp::object ElectrothermalActuation();
    bp::object HybridElectrostatics();
    bp::object HybridETMActuation();
    bp::object HybridETMPullin();

  };

}
#endif
