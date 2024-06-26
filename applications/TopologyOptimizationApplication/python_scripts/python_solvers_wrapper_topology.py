# Importing Kratos
import KratosMultiphysics

# Other imports
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type = solver_settings["solver_type"].GetString()

    if solver_settings.Has("time_integration_method"):
        time_integration_method = solver_settings["time_integration_method"].GetString()
    else:
        time_integration_method = "implicit" # defaulting to implicit time-integration

    try_import_custom_solver = False

    # Solvers for OpenMP parallelism
    if parallelism == "OpenMP":
        if solver_type.lower == "StaticSIMP" or solver_type == "staticSIMP":
            solver_module_name = "topology_optimization_simp_static_solver"
        else:
            available_solver_types = ["staticsimp"]
            try_import_custom_solver = True

    # Solvers for MPI parallelism
    elif parallelism == "MPI":
        err_msg =  "The requested parallel type MPI is not yet available!\n"
        raise Exception(err_msg)

    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    if try_import_custom_solver:
        KratosMultiphysics.Logger.PrintInfo("MechanicalSolversWrapper", 'Selected "solver_type" "{0}" not available in the python solvers wrapper, attempting to import custom solver from module "{0}"'.format(solver_type))
        try:
            solver = import_module(solver_type).CreateSolver(model, solver_settings)
            KratosMultiphysics.Logger.PrintInfo("MechanicalSolversWrapper", 'Using custom solver "{}", defined in module "{}"'.format(solver.__class__.__name__, solver.__class__.__module__))
            return solver
        except:
            err_msg =  'Importing custom solver from module "{}" failed.\n'.format(solver_type)
            err_msg += 'The requested solver type "{}" is not in the python solvers wrapper\n'.format(solver_type)
            err_msg += "Available options are: {}".format(', '.join(available_solver_types))
            raise Exception(err_msg)

    kratos_module = "KratosMultiphysics.StructuralMechanicsApplication"
    solver = import_module(kratos_module + "." + solver_module_name).CreateSolver(model, solver_settings)
    return solver


def CreateSolver(model, custom_settings):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)
