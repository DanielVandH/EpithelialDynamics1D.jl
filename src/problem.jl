"""
    CellProblem{F,DF,DP,G,GP}

Struct for representing a cell simulation problem. 

# Fields 
- `force_law::F`: The force law, e.g. `(δ, p) -> p.k * (p.s - δ)`.
- `force_law_parameters::Fp`: The parameters for the force law.
- `proliferation_law::G = (δ, p) -> zero(δ)`: The proliferation law, e.g. `(δ, p) -> p.β`.
- `proliferation_law_parameters::Gp = nothing`: The parameters for the proliferation law.
- `proliferation_period::Float64 = 1e-2`: How often to attempt a proliferation event.
- `initial_time::Float64 = 0.0`: The initial time.
- `final_time::Float64`: The final time.
- `fix_left::Bool = true`: Whether to fix the left endpoint.
- `fix_right::Bool = true`: Whether to fix the right endpoint.
- `damping_constant::Float64`: The damping constant.
- `initial_condition::Vector{Float64}`: The initial cell positions. Must be sorted.

# Solving 

You can solve a `CellProblem` like you would solve a problem from DifferentialEquations.jl, e.g. with 

    solve(prob, Tsit5(), saveat = 0.1)
"""
Base.@kwdef struct CellProblem{F,Fp,G,Gp}
    force_law::F
    force_law_parameters::Fp = nothing
    proliferation_law::G = (δ, p) -> zero(δ)
    proliferation_law_parameters::Gp = nothing
    proliferation_period::Float64 = 1e-2
    initial_time::Float64 = 0.0
    final_time::Float64
    fix_left::Bool = true
    fix_right::Bool = true
    damping_constant::Float64
    initial_condition::Vector{Float64}
end

"""
    SteadyCellProblem{M}

Defines a steady-state problem for a [`CellProblem`](@ref). Only has a 
single field, `prob::CellProblem`.

You can `solve` this problem as you would a `NonlinearProblem`, e.g. with 

    solve(prob, alg)

where `alg` is e.g. `DynamicSS(TRBDF2()))` from SteadyStateDiffEq.jl.
"""
struct SteadyCellProblem{M}
    prob::M
end