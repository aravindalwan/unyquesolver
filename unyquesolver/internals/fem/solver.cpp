#include "solver.hpp"

fem::Solver::Solver(int analysisType, bp::list pType) {

  c = new FEM_Common();
  useNonElast = false;
  useTherm = false;
  useElec = false;
  useElEs = false;

  switch(analysisType) {
  case MECHANICS:
    useNonElast = true;
    analysis = &fem::Solver::NLMechanics;
    break;
  case THERMAL:
    useTherm = true;
    analysis = &fem::Solver::Thermal;
    break;
  case ELECTRICAL:
    useElec = true;
    analysis = &fem::Solver::Electrical;
    break;
  case ELECTROTHERMAL:
    useElec = true;
    useTherm = true;
    analysis = &fem::Solver::Electrothermal;
    break;
  case ELECTROTHERMAL_ACTUATION:
    useElec = true;
    useTherm = true;
    useNonElast = true;
    analysis = &fem::Solver::ElectrothermalActuation;
    break;
  case HYBRID_ELECTROSTATICS:
    useElEs = true;
    analysis = &fem::Solver::HybridElectrostatics;
    break;
  case HYBRID_ETM_ACTUATION:
    useElEs = true;
    useTherm = true;
    useNonElast = true;
    analysis = &fem::Solver::HybridETMActuation;
    break;
  case HYBRID_ETM_PULLIN:
    useElEs = true;
    useTherm = true;
    useNonElast = true;
    analysis = &fem::Solver::HybridETMPullin;
    break;
  default:
    break;
  }

  // Store parameter types
  nParams = len(pType);
  paramType = new int[nParams];
  for (int i = 0; i < nParams; i++)
    paramType[i] = bp::extract<int>(pType[i])();

}
//------------------------------------------------------------------------------
fem::Solver::~Solver() {

  if (useNonElast)
    delete nelast;
  if (useTherm)
    delete therm;
  if (useElec)
    delete elec;
  if (useElEs)
    delete eles;

  delete c;
  delete s;
  delete paramType;

}
//------------------------------------------------------------------------------
void fem::Solver::Init(bp::list nodes, bp::list edges, bp::list elements) {

  // Initialize FEM surface object using nodes, edges and elements data
  s = new FEM_PhysicalDomain();
  s->SetID(1);
  s->SetMesh(nodes, edges, elements);

  if (useNonElast) {
    // Create and initialize Nelast
    if (c->DEBUG) cout<<"Initializing nonlinear mechanics module ... ";
    nelast = new NonElast(s, c);
    nelast->Init();
    if (c->DEBUG) cout<<"done"<<endl;
  }

  if (useTherm) {
    // Create and initialize Thermal
    if (c->DEBUG) cout<<"Initializing thermal module ... ";
    therm = new Therm(s, c);
    therm->Init();
    if (c->DEBUG) cout<<"done"<<endl;
  }

  if (useElec) {
    // Create and initialize Electrical
    if (c->DEBUG) cout<<"Initializing electrical module ... ";
    elec = new Elec(s, c);
    elec->Init();
    if (c->DEBUG) cout<<"done"<<endl;
  }

  if (useElEs) {
    //----------create and initialize ElEs
    if (c->DEBUG) cout<<"Initializing hybrid electrostatics module ... ";
    eles = new ElEs(s, c);
    eles->Init();
    if (c->DEBUG) cout<<"done"<<endl;
  }

}
//------------------------------------------------------------------------------
bp::object fem::Solver::Solve(bp::list params) {

  double param;
  int regionToMove;
  unyque::DMatrix displacement;
  FEM_Element *el;

  assert (len(params) == nParams);

  for (int i = 0; i < len(params); i++) {

    param = bp::extract<double>(params[i])();

    switch (paramType[i]) { // Check which parameter it corresponds to

    case YOUNGS_MODULUS: // Set mechanical elastic modulus
      nelast->EM = param*1e9;
      break;
    case THERMAL_CONDUCTIVITY: // Set thermal conductivity
      c->functions->CONST_FUNC = 1;
      c->functions->TC = param;
      break;
    case INTER_ELECTRODE_GAP: // Set electrostatic gap
      regionToMove = 2;
      originalGap = 2.0;
      newGap = param;
      displacement.resize(s->nnode,2);
      displacement = unyque::DMatrix_zero(s->nnode,2);
      for (int eid = 0; eid < s->nelem; eid++) {
        el = s->Elements[eid+1];
        if (el->reg == regionToMove) { // If element lies within regionToMove
          displacement(el->node1 - 1,1) = -(newGap - originalGap);
          displacement(el->node2 - 1,1) = -(newGap - originalGap);
          displacement(el->node3 - 1,1) = -(newGap - originalGap);
          displacement(el->node4 - 1,1) = -(newGap - originalGap);
          displacement(el->node5 - 1,1) = -(newGap - originalGap);
          displacement(el->node6 - 1,1) = -(newGap - originalGap);
        }
      }
      s->MoveMesh(displacement);
      break;
    case VOLTAGE: // Set voltage
      c->Phi_mult = param;
      break;
    default:
      cout << "Unknown random parameter type" << endl;
      break;
    }

  }

  return (this->*analysis)();

}
//------------------------------------------------------------------------------
bp::object fem::Solver::NLMechanics() {

  // Solve the static nonlinear elasticity problem
  nelast->SolveStatic();

  return bp::object(nelast->MaxAbsDisp(-1));

}
//------------------------------------------------------------------------------
bp::object fem::Solver::Thermal() {

  bp::list rvalue;

  //  // Solve the static heat conduction problem
  //  therm->SolveStatic();
  //  return bp::object(therm->MaxAbsTemp());

  // Solve the dynamic heat conduction problem
  double t_start = 0.0, t_end = 1.0, dt = 0.01, t;

  t = t_start;
  while (t <= t_end) {

    therm->SolveDynamic(t, dt);
    rvalue.append(bp::make_tuple(t, therm->MaxAbsTemp()));
    t = t + dt;

  }

  return rvalue;

}
//------------------------------------------------------------------------------
bp::object fem::Solver::Electrical() {

  // Solve the static current conduction problem
  elec->SolveStatic();

  return bp::object();

}
//------------------------------------------------------------------------------
bp::object fem::Solver::Electrothermal() {

  bp::list rvalue;

  // Solve the static current conduction problem
  elec->SolveStatic();

  //  // Solve the static heat conduction problem
  //  therm->SolveStatic();
  //  return bp::object(therm->MaxAbsTemp());

  // Solve the dynamic heat conduction problem
  double t_start = 0.0, t_end = 1.0e-6, dt = 0.002e-6, t;

  t = t_start;
  while (t <= t_end) {

    c->Phi_mult = 5*sin(10e6*4*atan2(1.0,1.0)*t);
    therm->SolveDynamic(t, dt);
    rvalue.append(bp::make_tuple(t, therm->MaxAbsTemp()));
    t = t + dt;

  }

  return rvalue;

}
//------------------------------------------------------------------------------
bp::object fem::Solver::ElectrothermalActuation() {

  // Solve the coupled electrothermal actuation problem
  unyque::DVector oldU, oldV;
  double err, eps = 1e-6, disp, maxValue = c->Phi_mult;
  FEM_Point *maxPoint = NULL;
  int i, counter = 1, numSteps = 1;

  oldU.resize(s->nnode); oldV.resize(s->nnode);
  oldU = unyque::DVector_zero(s->nnode); oldV = unyque::DVector_zero(s->nnode);

  do {
    c->Phi_mult = (((double)counter)/numSteps)*maxValue;
    i = 0;
    do {

      elec->SolveStatic();
      therm->SolveStatic();
      nelast->SolveStatic();

      err = max(ublas::norm_inf((s->U)-oldU), ublas::norm_inf((s->V)-oldV));
      oldU = (s->U); oldV = (s->V);

      i++;

    } while (err > eps);

    maxPoint = s->Nodes[nelast->MaxAbsDispPoint(-1)];
    disp = (s->V)(maxPoint->id - 1);

    counter++;

  } while (counter <= numSteps);

  return bp::object(disp);

}
//------------------------------------------------------------------------------
bp::object fem::Solver::HybridElectrostatics() {

  // Solve the coupled electrical electrostatics problem
  eles->SolveStatic();

  return bp::object();

}
//------------------------------------------------------------------------------
bp::object fem::Solver::HybridETMActuation() {

  // Solve the hybrid electrothermomechanical actuation problem
  unyque::DVector oldU, oldV;
  double err, prevErr = 1.0, eps = 1e-6, disp, maxValue = c->Phi_mult;
  FEM_Point *maxPoint = NULL;
  int i, counter = 1, numSteps = 1;
  bool pulledIn = false;

  oldU.resize(s->nnode); oldV.resize(s->nnode);
  oldU = unyque::DVector_zero(s->nnode); oldV = unyque::DVector_zero(s->nnode);

  do {

    c->Phi_mult = (((double)counter)/numSteps)*maxValue;
    i = 0; pulledIn = 0;
    do {

      eles->SolveStatic();
      therm->SolveStatic();
      nelast->SolveStatic();

      err = max(ublas::norm_inf((s->U)-oldU), ublas::norm_inf((s->V)-oldV));
      oldU = (s->U); oldV = (s->V);

      if ((i > 0) && ((err > prevErr)||(ublas::norm_inf(s->V) > newGap))) {
	pulledIn = true;
	break;
      }

      prevErr = err;

      i++;

    } while (err > eps);

    if (!pulledIn) {
      maxPoint = s->Nodes[nelast->MaxAbsDispPoint(-1)];
      disp = (s->V)(maxPoint->id - 1);
    } else {
      disp = -newGap;
      s->InitDOFs();
    }

    counter++;

  } while (counter <= numSteps);

  return bp::object(nelast->DispBoundaryEdge(1, -1));

}
//------------------------------------------------------------------------------
bp::object fem::Solver::HybridETMPullin() {

  // Solve the hybrid electrothermomechanical actuation problem
  unyque::DVector oldU, oldV;
  double err, prevErr, eps = 1e-6, epsV = 0.01, start = 0, stop = c->Phi_mult;
  bool pulledIn = false;

  oldU.resize(s->nnode); oldV.resize(s->nnode);

  while (stop - start > epsV) {

    oldU = unyque::DVector_zero(s->nnode);
    oldV = unyque::DVector_zero(s->nnode);
    prevErr = 1; pulledIn = false;
    c->Phi_mult = (start + stop)/2;

    do {

      eles->SolveStatic();
      //therm->SolveStatic();
      nelast->SolveStatic();

      err = max(ublas::norm_inf((s->U)-oldU), ublas::norm_inf((s->V)-oldV));
      oldU = (s->U); oldV = (s->V);

      if ((err > prevErr) || (ublas::norm_inf(s->V) > newGap)) {
	pulledIn = true;
	break;
      }

      prevErr = err;

    } while (err > eps);

    if (pulledIn) {
      stop = (start + stop)/2;
      s->InitDOFs();
    } else {
      start = (start + stop)/2;
    }

  }

  return bp::object((start + stop)/2);

}
