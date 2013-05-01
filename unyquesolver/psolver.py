'''Parametric solver that solves a PDE for different sets of parameters'''

__copyright__ = 'Copyright (C) 2011 Aravind Alwan'

__license__ = '''
This file is part of UnyQuE.

UnyQuE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

UnyQuE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For a copy of the GNU General Public License, please see
<http://www.gnu.org/licenses/>.
'''

import mpihelper
import logmanager
import unyquesolver._internals as internals
from mesh import PhysicalDomain

# MPI Tags
PARAMETER_SET = 1 # Tag for messages containing parameter sets
SOLUTION = 2 # Tag for messages containing solution to the PDE

# FEM Solver analysis modes
MECHANICS = 1
THERMAL = 2
ELECTRICAL = 3
ELECTROTHERMAL = 4
ELECTROTHERMAL_ACTUATION = 5
HYBRID_ELECTROSTATICS = 6
HYBRID_ETM_ACTUATION = 7
HYBRID_ETM_PULLIN = 8
HYBRID_ETM_DYNAMIC = 9

# Possible input parameters for FEM Solver
YOUNGS_MODULUS = 1
THERMAL_CONDUCTIVITY = 2
INTER_ELECTRODE_GAP = 3
VOLTAGE = 4
SPATIALLY_VARYING_BOUNDARY = 5

class ParametricSolver(object):
    '''Base class to solve a PDE for different sets of parameters.
    '''

    def __init__(self):
        self.comm = mpihelper.comm
        self.rank = mpihelper.rank
        self.number_of_processors = mpihelper.size

    def _init_solver(self):
        '''Initialize FEM solver assuming that the domain parameters and the UQ
        parameters have been defined.
        '''

        # Initialize physical domain on which PDE is solved
        self.pdomain = PhysicalDomain(*self.domain_parameters)

        # Remove SPATIALLY_VARYING_BOUNDARY from list of parameters before
        # initializing solver. This parameter will be handled separately
        parameters = [param for param in self.parameters
                      if param != SPATIALLY_VARYING_BOUNDARY]
        self._solver = internals.Solver(self.analysis, parameters)
        self._solver.Init(self.pdomain.nodes, self.pdomain.edges,
                          self.pdomain.elements)

    def _get_solution(self, pset):
        '''Compute the solution for the given parameter set.
        '''

        parameter_set = pset + self.fixed_parameter_values
        rvalue = None

        try:

            # Check if SPATIALLY_VARYING_BOUNDARY is one of the parameters.
            # Throws ValueError if it is not present.
            svb_index = self.parameters.index(SPATIALLY_VARYING_BOUNDARY)

            # Displace boundary of physical domain according to the function
            # given in the list of parameters
            self.pdomain.move_boundary(parameter_set[svb_index])

            # Re-initialize solver for the new physical domain
            self._solver.Init(self.pdomain.nodes, self.pdomain.edges,
                          self.pdomain.elements)

            # Remove the spatial variation function from parameter set before
            # calling the solver to compute the solution
            rvalue = self._solver.Solve(
                parameter_set[:svb_index] + parameter_set[(svb_index+1):])

        except ValueError:

            # SPATIALLY_VARYING_BOUNDARY is not one of the parameters, so call
            # the solver on the entire parameter set

            rvalue = self._solver.Solve(parameter_set)

        return rvalue

class ParametricSolverMaster(ParametricSolver):
    '''Class representing the root or master processor, that handles the task of
    distributing parameter sets among the processors and recieving the solutions
    that are returned. Distribution of parameter sets is implemented as a task
    queue.
    '''

    # Initialize loggers used in this class
    _log = logmanager.getLogger('unyquesolver.psolver')
    _results_log = logmanager.getLogger('unyquesolver.results')

    def __init__(self, domain_parameters, analysis, pvariable, pfixed):
        '''Initialize the solver on the master node.

        Arguments:
        domain_parameters -- Parameters used to initialize PhysicalDomain object
        analysis -- Number corresponding to one of the FEM Solver analysis modes
        pvariable -- List of FEM Solver input parameter types
        pfixed -- List of fixed parameters, specified as (type, value) tuples
        that give the type and the value respectively, of the input parameter
        '''

        super( ParametricSolverMaster, self ).__init__()
        if self.rank != mpihelper.MASTER:
            raise ParametricSolverError(
                'ParametricSolverMaster objects should be instantiated only ' +
                'on the root processor, given by mpihelper.MASTER.')

        self.domain_parameters = domain_parameters
        self.analysis = analysis
        self.parameters = pvariable + [param[0] for param in pfixed]
        self.fixed_parameter_values = [param[1] for param in pfixed]

        # Broadcast FEM solver configuration data to workers
        self._log.info('Broadcasting FEM solver parameters to workers')
        self.comm.bcast(self.domain_parameters, root = mpihelper.MASTER)
        self.comm.bcast(self.analysis, root = mpihelper.MASTER)
        self.comm.bcast(self.parameters, root = mpihelper.MASTER)
        self.comm.bcast(self.fixed_parameter_values, root = mpihelper.MASTER)

        # Initialize FEM solver
        self._log.info('Initializing FEM solver')
        self._init_solver()

        self.parameter_sets = None

    def __call__(self, parameter_sets, postprocessor = None,  tag = ''):
        '''Evaluate the solution for each of the parameters sets in the given
        list. The list of parameter sets is treated as a task queue, from which
        tasks are assigned sequentially to workers on a first-come-first-served
        basis. The results returned by the workers are indexed by the position
        variable and the final list of results is returned as a generator.

        Arguments:
        parameter_sets -- list of parameter sets
        postprocessor -- function used to post-process result before logging
        tag -- string used to identify the result at a later stage
        '''

        self._log.info('Working on %d parameter sets' % len(parameter_sets))

        # Notify workers that we have a new set of tasks
        self.comm.bcast(True, root = mpihelper.MASTER)

        self.parameter_sets = list(parameter_sets) # Create a copy

        if self.number_of_processors > 1: # Enable parallel processing

            # Denotes number of active workers
            active_workers = self.number_of_processors - 1

            # First assign a task to each worker
            for i in range(self.number_of_processors):
                if i == mpihelper.MASTER:
                    continue # Master doesn't do any work!
                else:
                    active_workers -= self._send_next_set(i)

            status = mpihelper.mpi.Status()
            while active_workers > 0:

                # Wait until one of the workers completes its task
                self.comm.Probe(
                    mpihelper.mpi.ANY_SOURCE, tag = SOLUTION, status = status)

                # Determine the rank of the worker and get the result
                waiting_worker = status.Get_source()
                position, result = self.comm.recv(source = waiting_worker,
                                                  tag = SOLUTION)

                # Assign the waiting worker a new task
                active_workers -= self._send_next_set(waiting_worker)

                # Log this result
                meta_data = {'replicate': position, 'raw_result': result,
                             'tag': tag}
                logmessage = 'Replicate: %(replicate)3d     Tag: %(tag)s'
                if postprocessor:
                    meta_data['result'] = postprocessor(result)
                    logmessage = logmessage + '     Result: %(result)1.6f'
                self._results_log.log(logmanager.REPLICATE,
                                      logmessage, meta_data, extra = meta_data)

                # Yield the result along with its position in the results list
                yield (position, result)

        else: # Run parametric solver on a single processor

            while self.parameter_sets:

                result = self._get_solution(self.parameter_sets.pop())
                position = len(self.parameter_sets)

                # Log this result
                meta_data = {'replicate': position, 'raw_result': result,
                             'tag': tag}
                logmessage = 'Replicate: %(replicate)3d     Tag: %(tag)s'
                if postprocessor:
                    meta_data['result'] = postprocessor(result)
                    logmessage = logmessage + '     Result: %(result)1.6f'
                self._results_log.log(logmanager.REPLICATE,
                                      logmessage, meta_data, extra = meta_data)

                # Yield the result along with its position in the results list
                yield (position, result)

    def _send_next_set(self, destination):
        '''Send the next task to the worker at the given destination. If all the
        tasks are complete, then send a position flag of -1 to indicate that
        there are no more tasks for that worker and that it may stop asking for
        tasks. Return False if the worker if the worker was assigned a task and
        is still active and True if it was asked to stop.
        '''

        try:
            pset = self.parameter_sets.pop() # Get task at the end of the list
            position = len(self.parameter_sets)
        except IndexError: # No more tasks
            pset = []
            position = -1

        # Send task and its position to destination
        self.comm.send((position,pset), dest = destination,
                       tag = PARAMETER_SET)

        # Return True if the task list is empty
        return position < 0

    def stop(self):
        '''Signal workers to stop working and tell them to shut down.
        '''
        self.comm.bcast(False, root = mpihelper.MASTER)

class ParametricSolverWorker(ParametricSolver):
    '''Class representing the worker processor, that does the actual work of
    solving the PDE.
    '''

    def __init__(self):
        super( ParametricSolverWorker, self ).__init__()
        if self.rank == mpihelper.MASTER:
            raise ParametricSolverError(
                'ParametricSolverWorker objects should be instantiated only ' +
                'on processors other than the master, designated by ' +
                'mpihelper.MASTER.')

        # Receive solver configuration data from master
        self.domain_parameters = self.comm.bcast(root = mpihelper.MASTER)
        self.analysis = self.comm.bcast(root = mpihelper.MASTER)
        self.parameters = self.comm.bcast(root = mpihelper.MASTER)
        self.fixed_parameter_values = self.comm.bcast(root = mpihelper.MASTER)

        # Initialize FEM solver
        self._init_solver()
        self.continue_working = True

    def run(self):
        '''Primary method that is called on the worker processors. Each worker
        continues to listen for incoming parameter set data as long as the
        value of position is non-negative. The continue_working attribute is
        used to establish an outer loop that continues until a False value is
        broadcast by the master processor, when its stop() method is called.
        '''

        # First check whether we need to listen for parameter sets
        self.continue_working = self.comm.bcast(root = mpihelper.MASTER)

        while self.continue_working:

            # Get the first parameter set. Treat position as a flag that is
            # non-negative for valid parameter sets and -1 to signal the
            # worker to stop asking for parameter sets
            position, current_set = self._receive_next_set()

            while not position < 0:

                # Compute the solution for the current set
                result = self._get_solution(current_set)

                # Send the result to master along with the position
                self.comm.send((position, result), dest = mpihelper.MASTER,
                               tag = SOLUTION)

                # Get the next set
                position, current_set = self._receive_next_set()

            # Current sequence of tasks is complete. Check whether to continue
            self.continue_working = self.comm.bcast(root = mpihelper.MASTER)

    def _receive_next_set(self):
        '''Receive a new tuple containing the parameter set and its position.
        '''
        return self.comm.recv(source = mpihelper.MASTER, tag = PARAMETER_SET)

class ParametricSolverError(Exception):
    '''Class of exceptions raised by the parametric solver.
    '''

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message
