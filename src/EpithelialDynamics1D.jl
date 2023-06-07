module EpithelialDynamics1D

using SparseArrays
using SciMLBase
using CommonSolve
using DiffEqCallbacks
using Random
using StatsBase
using ForwardDiff
using FiniteVolumeMethod1D
using MovingBoundaryProblems1D
using DataInterpolations

export CellProblem, SteadyCellProblem
export solve
export cell_densities, cell_midpoints, node_densities
export continuum_limit 
export get_knots
export cell_numbers, leading_edges
export integrate_pde

include("problem.jl")
include("equations.jl")
include("proliferation.jl")
include("solve.jl")
include("statistics.jl")
include("continuum.jl")

end # module